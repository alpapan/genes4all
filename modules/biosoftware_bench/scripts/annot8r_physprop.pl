#!/usr/bin/perl -w
#$Id$

=head1 NAME

 annot8er_physprop.pl version 0.2

=head1 USAGE

 Calculate physical properties for peptide sequences

=head1 DEPENDENCIES

 DBD::Pg 
 BioPerl

=head1 AUTHOR & LICENSE

 Last updated 10/04/2005 by Ralf Schmid
 Hacked by Alexie Papanicolaou for automation.

 Copyright (C) 2004 Ralf Schmid
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

=cut

use strict;
use File::stat;
use Term::ANSIColor;
use Bio::SeqIO;
use Bio::Tools::pICalculator;
use Bio::Tools::SeqStats;
use DBI;
use Data::Dumper;
use DBD::Pg;
use Time::localtime;
my $version_number = "0.2";
my @PATH = split( ":", "$ENV{'PATH'}" );
#############################################################################
############################ start main menu ################################
#############################################################################
#### think about configuration file for following generally used variables
my $seq2phys;         # sequences
my $seq2phys_phys;    #tab delimited output file
my (@filelist)=@ARGV;
if   (@filelist) { &calculate_properties();exit; }
else         { &options(); }
##################################################################################
#################  1-Calculate physical properties  ##############################
##################################################################################
sub calculate_properties() {
####
	print colored( "\n\t##### Calculate physical properties #####\n",
				   "white bold", "on_black" );
	print
"\nThis facility uses Bioperl modules to calculate pI, charge at three different";
	print
"\npH-values and Mw for a set peptides and stores the results in a PartiGene database\n";
	opendir( DIR, "." );
	if (!@filelist) {
		@filelist =
		  grep { ( $_ =~ /_pro.fsa/ ) || ( $_ =~ /translations/ ) }
		  readdir(DIR);
	}
	my $files = @filelist;
	print "Found $files files.\n";
	if ( $files eq 0 ) { exit; }
	foreach $seq2phys (@filelist) {
		$seq2phys_phys = $seq2phys;
		$seq2phys_phys =~ s/.fsa$//;
		$seq2phys_phys .= ".phys";
		if (-s $seq2phys_phys){warn "Output $seq2phys_phys already exists\n";next;}
		open( FILEPHYS, ">$seq2phys_phys" )  || die "Can't open $seq2phys_phys\n";
		print FILEPHYS "#### id ####\tlen\tundef\tminmw\tmaxmw\tchpos\tchneg\tiep\t c5\t c7\t c9\n";
		print "\nCalculating properties now. Please wait ...\n";
		my $in = Bio::SeqIO->new(
								  -file     => "$seq2phys",
								  -format   => 'Fasta',
								  -alphabet => 'protein'
		);
		my $calc =
		  Bio::Tools::pICalculator->new( -places => 2,
										 -pKset  => 'EMBOSS' );

		while ( my $seq = $in->next_seq ) {
			my $string = $seq->seq();
			my $length = length $string;
			$calc->seq($seq);
			my $seq_id = $seq->id;

			#print "Processing $seq_id\n";
			my $iep  = $calc->iep;
			my $ch_5 = $calc->charge_at_pH("5");
			my $ch_7 = $calc->charge_at_pH("7");
			my $ch_9 = $calc->charge_at_pH("9");
			my $seqobj =
			  Bio::PrimarySeq->new( -seq      => "$string",
									-alphabet => 'protein' );
			my $weight   = Bio::Tools::SeqStats->get_mol_wt($seqobj);
			my $minmw    = $$weight[0];
			my $maxmw    = $$weight[1];
			my $hash_ref = Bio::Tools::SeqStats->count_monomers($seqobj);
			my %hash = %$hash_ref;    # dereference hash to avoid warning
			my $ch_pos = ( $hash{"K"} + $hash{"R"} );
			my $ch_neg = ( $hash{"D"} + $hash{"E"} );
			my $undef  = '0';

			if ( $hash{"X"} ) {
				$undef = $hash{"X"};
			}
			print FILEPHYS sprintf(
						   "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.1f\t%.1f\t%.1f\n",
						   $seq_id, $length, $undef, $minmw, $maxmw, $ch_pos,
						   $ch_neg, $iep,    $ch_5,  $ch_7,  $ch_9
			);
		}
		close FILEPHYS;
		print "\nPhysical properties have been saved in $seq2phys_phys\n";
	}
}
##################################################################################
#################  2-Databasing  #################################################
##################################################################################
sub databasing() {
#### parse files, create and populate db
	system("clear");
	print colored( "\n\t##### databasing #####\n", "white bold", "on_black" );
	print "\nThis menu populates a PartiGene like database with the\n";
	print "physical properties calculated in option 1\n\n";
#### First some checks...pg running, partigene db, step1 output
	&postmaster_check;
#### check for partigene.conf file and partigene db
	my $filename = "~/.partigene.conf";
	my $pg_database;
	my $pg_flag;
	$filename =~ s{ ^ ~ ( [^/]* ) }
              { $1
                    ? (getpwnam($1))[7]
                    : ( $ENV{HOME} || $ENV{LOGDIR}
                         || (getpwuid($>))[7]
                       )
              }ex;

	if ( -e "$filename" ) {
		open( CONFILE, "$filename" ) || die "Can't open configuration file\n";
		while ( my $line = <CONFILE> ) {
			if ( $line =~ /^DATABASE\=(.+)/i ) { $pg_database = $1; }
		}
		close(CONFILE);
#### user has a Partigene db defined in the config file
		print "\nDo you want to use PartiGene database: $pg_database ";
		print "to store physical property information?";
		$pg_flag = &yes_no();
	}
#### don't use Pg database
	if ( $pg_flag != 1 ) {
		print "What database do you want to create or use?\n";
		$pg_database = <STDIN>;
		chomp $pg_database;
	}
	my $conn =
	  DBI->connect( "dbi:Pg:dbname=$pg_database", "", "", { PrintError => 0 } )
	  ;    #Last two values would be user/pass.
	my $update_flag = 2;
	if ( !$conn ) {    ### Couldn't connect to the database
		print "\nCouldn't connect to the database $pg_database";
		print "\nDo you want to create it?";
		my $answer = &query_for_continue;
		&create_physpropdb( $pg_database, "1", "1" )
		  ;            ### arguments: <db, does db exist, does table exist>
	} else {
		my @table = $conn->tables( '', '', undef, 'TABLE' );
		my $physprop_table_flag = 0;
		for ( my $n = 0 ; $n < @table ; $n++ ) {
			$table[$n] =~ s/public\.//
			  ; #get rid of "public." which is present in some versions of DBD.Pg
			if ( $table[$n] eq "physprop" ) { $physprop_table_flag = 1; }
		}
		if ( $physprop_table_flag == 0 ) {
			print colored(
"\n$pg_database does not contain physprop table, do you want to create it?",
				"red bold"
			);
			my $answer = &yes_no();
			if ( $answer == 1 ) {
				&create_physpropdb( $pg_database, "0", "1" );
			} else {
				&query_for_exit;
			}
			print "\nphysprop table created.\n";
		} else {
			print "\nphysprop entries are already existing for $pg_database.\n";
			$update_flag = 1;
		}

		#}
		$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";
	}
#### file check and get physprop results
	my @list = glob("*.phys");
	if (@list) {
		foreach $seq2phys_phys (@list) {
			open( FILEPHYS, "<$seq2phys_phys" )
			  || die "Can't open $seq2phys_phys\n";
			open( FILEPHYS, "<$seq2phys_phys" )
			  || die "Can't open $seq2phys_phys\n";
			print
"\nData file $seq2phys_phys found. \nNow parsing results. Please wait ...\n\n";
			my $conn = DBI->connect( "dbi:Pg:dbname=$pg_database", "", "", )
			  || die;    #Last two values would be user/pass.
			my $insert = $conn->prepare(
"INSERT INTO physprop (pept_id, num_res, undef_res, minmolw, maxmolw, pos_res, neg_res, pisoel, char_ph5, char_ph7, char_ph9, cluster_code) values (?,?,?,?.?,?,?,?,?,?,?,?)"
			);
			my @data;
			my $n = 0;
			while ( my $line = <FILEPHYS> ) {

				unless ( $line =~ /^###/ ) {
					@data = split( "\t", $line );
					my $cluster_code = substr( $data[0], 0, 3 );
					$cluster_code =~ s/(\w\w)P/$1C/;
					if ( $update_flag == 1 ) {    ### update
						my $check = $conn->do(
							  "DELETE from physprop WHERE pept_id='$data[0]';");
					}

					#debug
					my $print =
"INSERT INTO physprop (pept_id, num_res, undef_res, minmolw, maxmolw, pos_res, neg_res, pisoel, char_ph5, char_ph7, char_ph9, cluster_code) values ('$data[0]','$data[1]','$data[2]','$data[3]','$data[4]','$data[5]','$data[6]','$data[7]','$data[8]','$data[9]','$data[10]','$cluster_code');\n";
					my $success = $conn->do(
"INSERT INTO physprop (pept_id, num_res, undef_res, minmolw, maxmolw, pos_res, neg_res, pisoel, char_ph5, char_ph7, char_ph9, cluster_code) values ('$data[0]','$data[1]','$data[2]','$data[3]','$data[4]','$data[5]','$data[6]','$data[7]','$data[8]','$data[9]','$data[10]','$cluster_code');"
					);

#     my  $success=$insert->execute($data[0],$data[1],$data[2],$data[3],$data[4],$data[5],$data[6],$data[7],$data[8],$data[9],$data[10],$cluster_code);
					if ( !$success ) { print $print; die; }
					$n++;
					printf( "\r%9d sequences processed so far", $n );
				}
			}
			$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";
		}
	} else {
		print "No files found. exit;\n";
	}
	&query_for_exit();
}
###################################################################################
###################################################################################
###                                                                             ###
###                               all the subs                                  ###
###                                                                             ###
###################################################################################
###################################################################################
#########################################################################################################################
sub options () {
#### option selection
	my $answer = 0;
	while ( $answer != 3 ) {
		$answer = &title_page();
		if ( $answer == 1 ) { &calculate_properties(); }
		if ( $answer == 2 ) { &databasing(); }
	}
	system("clear");
	exit();
}
############################################################################################################################
###################################################################################
sub title_page() {
#### intro & sub-menu selection
	print_title();
	print "\n\tBefore proceeding please make sure you are in an appropriate\n";
	print "\tdirectory.\n\n";
	print "\t\t1. Calculate physical properties\n";
	print "\t\t2. Databasing\n";
	print "\t\t3. Quit.\n";
	my $flag = 0;
	my $answer;

	while ( $flag == 0 ) {
		$answer = <>;
		if ( $answer =~ /^[1|2|3]$/ ) { $flag = 1; next; }
		else {
			print
" You have entered $answer This is not an option. Please try again\n";
		}
	}
	return $answer;
}
####################################################################################
############################################################################################################################
sub print_title() {
##### displays title
	print colored( "\n\n\n", "white bold", "on_black" );
	print colored(
			  "\t###########################################################\n",
			  "white bold", "on_black" );
	print colored(
			  "\t###                                                     ###\n",
			  "white bold", "on_black" );
	print colored(
			  "\t###             annot8r_physprop.pl                     ###\n",
			  "white bold", "on_black" );
	print colored(
			  "\t###     a script to calculate physical properties       ###\n",
			  "white bold", "on_black" );
	print colored(
			  "\t###     for peptide sequences Vs $version_number        ###\n",
			  "white bold", "on_black" );
	print colored(
			  "\t###                                                     ###\n",
			  "white bold", "on_black" );
	print colored(
			  "\t###     EGTDC 2004                                      ###\n",
			  "white bold", "on_black" );
	print colored(
			  "\t###                                                     ###\n",
			  "white bold", "on_black" );
	print colored(
			  "\t###     For news and upgrades and help:                 ###\n",
			  "white bold", "on_black" );
	print colored(
			 "\t###     nematode.bioinf\@ed.ac.uk                        ###\n",
			 "white bold", "on_black" );
	print colored(
			  "\t###                                                     ###\n",
			  "white bold", "on_black" );
	print colored(
			  "\t###     Help for EG-Awardees:                           ###\n",
			  "white bold", "on_black" );
	print colored(
			 "\t###     helpdesk\@envgen.nox.ac.uk                       ###\n",
			 "white bold", "on_black" );
	print colored(
			  "\t###                                                     ###\n",
			  "white bold", "on_black" );
	print colored(
		  "\t###########################################################\n\n\n",
		  "white bold", "on_black" );
}
###########################################################################################################################
###################################################################################
sub get_file() {
#### get file from user input, check existence
	my $flag = 0;
	my $file_name;
	while ( $flag == 0 ) {
		$file_name = <STDIN>;
		chomp $file_name;
		$file_name =~ s/\s+//g;
		if ( -e $file_name ) { $flag = 1; }
		else {
			print "\nCan't find $file_name, do you want to try again?";
			my $answer = &yes_no();
			if ( $answer != 1 ) { exit(); }
		}
	}
	return "$file_name";
}
###################################################################################
####################################################################################
sub yes_no() {
#### returns 1 for y
	my $yflag = 0;
	print colored( " [y/n] : ", "green bold" );
	my $input = '';
	while ( $input !~ /y|n/i ) {
		print "\b";
		$input = <STDIN>;
		chomp $input;
	}
	if ( $input =~ /^y/i ) { $yflag = 1; }
	return $yflag;
}
####################################################################################
#######################################################################################
sub find_program() {
#### search for argument in path, exit if not found
	my $prog     = $_[0];
	my $pathflag = 0;
	my $path;
	my $finalpath;
	foreach $path (@PATH) {
		if ( -f "$path/$prog" ) {
			$pathflag  = 1;
			$finalpath = $path;
			last;
		}
	}
	if ( $pathflag == 0 ) {
		print colored( "\nCan't find the $prog utility on your system\n",
					   "red bold" );
		print colored( "Please ensure it is installed and in your path\n\n",
					   "red bold" );
		exit();
	} else {
		return "$finalpath/$prog";
	}
}
###############################################################################################################
#########################################################################################################################
sub postmaster_check() {
#### check for postmaster/postgresql process
	my $postmaster = `ps -e|grep postmaster`;  ### See if the process is running
	if ( !$postmaster ) {
		print colored( "\n#### Postmaster is not running ####\n", "red bold" );
		print colored(
"Please ensure that postgreSQL is correctly installed and running\n",
			"red bold"
		);
		exit();
	} else {
		print "CHECK POSTGRESQL   OK => Postmaster running\n\n";
	}
#### check whether user does exist
	my $user_status = system("psql -l > /dev/null")
	  ;    ### command will fail unless (postgresql) user exists
	unless ( $user_status == 0 ) {
		my $username = `whoami`;
		chomp $username;
		print colored( "\n\t#### CONNECTION TO POSTGRESQL FAILED ####\n",
					   "red bold" );
		print colored(
			 "Most likely you have forgotten to run \"createuser $username\"\n",
			 "red bold"
		);
		print colored( "during the postgreSQL setup\n", "red bold" );
		exit();
	}
}
######################################################################################################################
###########################################################################################################################
sub query_for_exit() {

	#exits program if 'n' entered back to main for y
	print "\nWould you like to continue? ";
	print colored( " [y/n] : ", "green bold" );
	my $input = '';
	while ( $input !~ /y|n/i ) {
		print "\b";
		$input = <>;
		chomp $input;
	}
	if ( $input =~ /^y/i ) { print "Back to main menu\n";   &options; }
	if ( $input =~ /^n/i ) { print "Exiting the program\n"; exit(); }
}
####################################################################################
###########################################################################################################################
sub query_for_continue() {

	#exits program if 'n' entered carries on for 'y'
	#   print "\nWould you like to continue? ";
	print colored( " [y/n] : ", "green bold" );
	my $input = '';
	while ( $input !~ /y|n/i ) {
		print "\b";
		$input = <>;
		chomp $input;
	}
	if ( $input =~ /^n/i ) { print "Exiting the program\n"; exit(); }
}
####################################################################################
##############################################################################################
sub create_physpropdb() {
#### one basic table
	my $database     = shift;
	my $new_db       = shift;
	my $table_flag   = shift;
	my $createdb_exe = &find_program("createdb");
	if ( $new_db == 1 ) { system("$createdb_exe $database >& /dev/null"); }
	my $conn = DBI->connect( "dbi:Pg:dbname=$database", "", "" );
	if ( $table_flag == 1 ) {
		my $result = $conn->do(
"create table physprop (pept_id varchar(20) not null primary key, num_res int, undef_res int, minmolw double precision, 
    maxmolw double precision, pos_res int, neg_res int, pisoel double precision,  char_ph5 double precision, char_ph7 double precision, char_ph9 double precision,cluster_code char(3));"
		);
	}
	my $errorMessage = $conn->errstr;
	if ($errorMessage) { print "$errorMessage\n"; }
	$conn->disconnect or warn "Disconnection failed: $DBI::errstr\n";
}
###############################################################################################################
###   fin de siecle
###############################################################################################################
