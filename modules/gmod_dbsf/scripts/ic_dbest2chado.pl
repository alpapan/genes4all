#!/usr/bin/perl -w

=head1 NAME 

ic_dbest2chado.pl

=head1 USAGE

 * 'type:s'	=> What kind of input to expect: dbEST, SRA or tracexml (case insensitive). Only one type per run
 * 'db|dbestfile:s'	=> DBEST file with ESTs and libraries. Required unless an XML is given
 * 'xml:s{,}'		=> SRA or NCBI TRACEINFO XML files (one or more). Required unless a dbEST file is given above
 * 'dsn|database_connection:s' => Chado database connection.
   'lineage:s'		=> Only process those that match a taxonomy group. Defaults to Arthropoda, set to 0 to deactivate
   'verbose'			=> Debug output
   'phylogeny'	=> Optionally build phylogeny in Chado. Can also be run standalone without XML/dbEST
   'seqtech:s'	=> Sequencing technology for this run: sanger or 454. Assumes sanger for dbEST and 454 for SRA but can be overruled here.
   'libname:s'	=> Library name if parsing SRA or NCBI Traceinfo XML
   'orgname:s'	=> Organism name if not present in input files or in chado
    'flat:ncbi_flatfiledir:s'   => Directory containing flatfiles of NCBI Taxonomy (def. /db/ncbi_taxonomy/) 


=head1 AUTHORS

 Alexie Papanicolaou 1 2

	1 Max Planck Institute for Chemical Ecology, Germany
	2 Centre for Ecology and Conservation, University of Exeter, UK
	alexie@butterflybase.org

=head1 TIPS

 An example of a DSN file:

 a dsn file is simple text file which has one line like this:
 dbi:Pg:dbname=InsectaCentral;host=localhost;port=15432;user=postgres;password=123321
 or instead of a file, pass the above string as the option. 

 The NCBI Tax flatfile can be acquired from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

 We will look for specific fields and therefore feel free to grep a subfile of the larger dbestfile (which is more than 200 Gb). 
 This will speed up things (by creating a ca 5Gb file) but it might take ca an hour to create it.
 Note that you will need the original file if you wish to include the description tag (which could just be fetched from dbEST rather than cluttering your database.
 
 $ egrep "^dbEST Id:|^GenBank Acc:|^dbEST lib id:" "file"  > "newfile"
  
=head1 DISCLAIMER & LICENSE

This software is released under the GNU General Public License version 3 (GPLv3).
It is provided "as is" without warranty of any kind. You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html. 
Please note that incorporating the whole software or parts of its code in proprietary software is prohibited under the current license.

=head1 BUGS & LIMITATIONS

 None known. Please report them.

 You will need a preloaded chado database in order to get get the correct CV terms for GO

=cut

use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Time::Progress;
use Bio::Taxon;
use Bio::DB::Taxonomy;
use DBI;
use XML::Simple qw(:strict);
use IO::File;
$| = 1;
my ( $dbestfile, $dsn, @xmls, $request_phylogeny, $user_libname, $lib454id,
	 $user_orgname );
my $flatfile_dir='/db/ncbi_taxonomy/';
my $type         = "dbEST";
my $lineage_user = "Arthropoda";
my $debug        = 0;
my $seqtech;
GetOptions(
	   'db|dbestfile:s'              => \$dbestfile,
	   'd|dsn|database_connection:s' => \$dsn,
	   'lineage:s'                   => \$lineage_user,
	   'verbose:i'                   => \$debug,               # 0 to switch off
	   'xml:s{,}'                    => \@xmls,
	   'type:s'                      => \$type,
	   'phylogeny'                   => \$request_phylogeny,
	   'seqtech:s'                   => \$seqtech,
	   'libname:s'                   => \$user_libname,
	   'orgname:s'                   => \$user_orgname,
	   'flat|ncbi_flatfiledir:s'     => \$flatfile_dir,
	   
);
unless ( (@xmls) || ( $dbestfile && -s $dbestfile ) || $request_phylogeny ) {
	pod2usage;
	die;
}
if ( $type =~ /^dbEST/i ) {
	$type    = 'dbEST';
	$seqtech = 'Sanger';
} elsif ( $type =~ /^SRA/i ) {
	$type    = 'SRA';
	$seqtech = '454';
} elsif ( $type =~ /^tracexml/i ) {
	$type = 'TRACEXML';
	if ( !$seqtech ) {
		die "Please specify sequencing technology seqtech (one tech per XML)\n";
	}
	if ( !$user_libname ) { die "Please specify the name of the library\n"; }
#	if ( !$user_orgname ) { die "Please specify organism name\n"; }
	if ( !$dsn ) { die "Please specify DSN file for chado connectivity.\n"; }
} else {
	pod2usage;
}
print "Processing $type";
if ($seqtech)      { print " with $seqtech technology"; }
if ($user_libname) { print " with library name as '$user_libname'"; }
if ($user_orgname) { print " with organism $user_orgname"; }
print "\n";
if ( length($lineage_user) < 2 ) {
	warn "Deactivating lineage search\n" if $debug;
	undef($lineage_user);
}

my $taxondb = Bio::DB::Taxonomy->new(-source => 'entrez');
if ($flatfile_dir && -d $flatfile_dir && -s $flatfile_dir.'/nodes.dmp' && -s $flatfile_dir.'/names.dmp'){
  $taxondb = Bio::DB::Taxonomy->new(-source => 'flatfile',-nodesfile => $flatfile_dir.'/nodes.dmp',-namesfile => $flatfile_dir.'/names.dmp',-directory => $flatfile_dir);
}


my $timer = new Time::Progress;
my $counter;
my %hash_est;    # the one to be processed at the end

#retired:
my %hash_est_temp
  ;    # the one to be transfered to hash_est if and only if passes criteria.
my %hash_lib;    # build slowly as we process above
my %skipped_libs
  ;    # libraries which didn't make the lineage cutoff are skipped next time
my %rank_cache;    # cached cvterms for taxonomy ranks
my $xmlurl =
"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucest&retmode=xml&tool=est2assembly&id=";
mkdir('ncbi_ests');
my $date = `date -R`;
chomp($date);
my %get_libs;
my $dbh = &connecttodb($dsn) if $dsn;

if ($dbestfile) {
	my @data;
	my $est_id;    # global
	my $genbank;
	my $total_lines = `wc -l $dbestfile`;
	$total_lines =~ /^\s*(\d+)/;
	$total_lines = $1;
	$timer->attr( min => 0, max => $total_lines );
	$timer->restart;
	open( DBEST, $dbestfile );

# the idea of parsing the dbest file is retired because of lack of sufficient information regarding the organism.
# instead we will download each est's dbest record as xml (i.e. hammer dbest) using the genbank record
# we will also store the lib id so that we don't have to fetch the xml for ESTs coming from the same library.
	while ( my $line = <DBEST> ) {
		$counter++;
		if ( $counter =~ /0000$/ ) {
			print $timer->report( "eta: %E min, %40b %p\r", $counter );
		}
		chomp($line);
		@data = split( ":", $line );
		if ( $data[1] ) {
			$data[1] =~ s/^\s+|\s+$//;
			if ( $data[0] eq 'dbEST Id' ) {
				$est_id = $data[1];
				undef(@data);
				next;
			}
			if ($est_id) {
				if ( $data[0] eq 'GenBank Acc' ) {
					unless ( $data[1] ) {
						warn "No genbank accession for $est_id. Skipping..."
						  if $debug;
						undef($est_id);
						undef(@data);
						next;
					}
					$genbank = $data[1];
					undef(@data);
					next;
				} elsif ( $data[0] eq 'dbEST lib id' ) {
					unless ( $data[1] ) {
						warn "No db est library id for $est_id. Skipping..."
						  if $debug;
						undef($est_id);
						undef(@data);
						next;
					}
					$get_libs{ $data[1] } = $genbank if !$get_libs{ $data[1] };
					undef(@data);
					next;
				} else {
					undef(@data);
					next;
				}
			}
		}
	}
	close(DBEST);
	print "\n";
	&get_ncbi_dbest();
## do it again to get the ESTs that actually belong to the libraries we asked for.
	if ($dbestfile) {
		my @data;
		my $est_id;    # global
		my $genbank;
		$counter = 0;
		my $total_lines = `wc -l $dbestfile`;
		$total_lines =~ /^\s*(\d+)/;
		$total_lines = $1;
		$timer->attr( min => 0, max => $total_lines );
		$timer->restart;
		open( DBEST, $dbestfile );

# the idea of parsing the dbest file is retired because of lack of sufficient information regarding the organism.
# instead we will download each est's dbest record as xml (i.e. hammer dbest) using the genbank record
# we will also store the lib id so that we don't have to fetch the xml for ESTs coming from the same library.
		while ( my $line = <DBEST> ) {
			$counter++;
			if ( $counter =~ /0000$/ ) {
				print $timer->report( "eta: %E min, %40b %p\r", $counter );
			}
			chomp($line);
			@data = split( ":", $line );
			if ( $data[1] ) {
				$data[1] =~ s/^\s+|\s+$//;
				if ( $data[0] eq 'dbEST Id' ) {
					$est_id = $data[1];
					undef(@data);
					next;
				}
				if ($est_id) {

					# first is genbank and then is library info
					if ( $data[0] eq 'GenBank Acc' ) {
						unless ( $data[1] ) {
							warn "No genbank accession for $est_id. Skipping..."
							  if $debug;
							delete $hash_est{$est_id};
							undef($est_id);
							undef(@data);
							next;
						}
						$hash_est{$est_id}{'genbank'} = $data[1];
						undef(@data);
						next;
					} elsif ( $data[0] eq 'dbEST lib id' ) {
						unless ( $data[1] && $hash_lib{ $data[1] } ) {
							delete $hash_est{$est_id};
							undef($est_id);
							undef(@data);
							next;
						}
						$hash_est{$est_id}{'dbEST_libid'} = $data[1];
						$hash_est{$est_id}{'dbEST_id'}    = $est_id;
						undef(@data);
						next;
					} else {
						undef(@data);
						next;
					}
				}
			}
		}
		close(DBEST);
		print "\n";
	}
	if (%hash_est) {
		print "Acquiring data from NCBI...\n";
		print "Storing EST information...\n";
		print "\t in XML\n";
		&store_in_xml("est");
		undef(%hash_est);
	}

	# process library
	if (%hash_lib) {
		warn Dumper %hash_lib if $debug;
		print "Storing library information...\n";
		print "\t in XML\n";
		&store_in_xml("lib");
	}
} elsif (@xmls) {
	if ( $type eq 'dbEST' ) {
		my $xml_type;
		foreach my $xmlfile (@xmls) {
			unless ( -s $xmlfile ) { warn "Cannot find $xmlfile.\n"; next; }
			if    ( $xmlfile =~ /est\./ ) { $xml_type = "est"; }
			elsif ( $xmlfile =~ /lib\./ ) { $xml_type = "lib"; }
			else {
				warn "I don't know what type this file is : $xmlfile\n";
				next;
			}
			print "Processing $xmlfile as $xml_type\n";
			if ( $xml_type eq "est" ) {
				my $xs =
				  XML::Simple->new(
									NormaliseSpace => 2,
									KeyAttr        => ['dbEST_id'],
									ForceArray     => 0,
				  );
				my $ref = $xs->XMLin($xmlfile)->{'EST'};
				foreach my $dbest_id ( keys %{$ref} ) {
					$hash_est{$dbest_id}{'libid'} =
					  $ref->{$dbest_id}->{'dbEST_libid'};
					$hash_est{$dbest_id}{'genbank'} =
					  $ref->{$dbest_id}->{'genbank'};
				}
			} elsif ( $xml_type eq "lib" ) {
				my $xs =
				  XML::Simple->new(
									NormaliseSpace => 2,
									KeyAttr        => ['dbEST_libid'],
									ForceArray     => 0,
				  );
				my $ref = $xs->XMLin($xmlfile)->{'library'};
				foreach my $dbEST_libid ( keys %{$ref} ) {
					$hash_lib{$dbEST_libid}{'ncbi_taxid'} =
					  $ref->{$dbEST_libid}->{'ncbi_taxid'};
					$hash_lib{$dbEST_libid}{'libname'} =
					  $ref->{$dbEST_libid}->{'clone_lib'}
					  if $ref->{$dbEST_libid}->{'clone_lib'};
					$hash_lib{$dbEST_libid}{'orgname'} =
					  $ref->{$dbEST_libid}->{'organism'}
					  if $ref->{$dbEST_libid}->{'organism'};
					foreach my $key ( keys %{ $ref->{$dbEST_libid} } ) {
						if (    $key eq 'note'
							 && $ref->{$dbEST_libid}->{$key} =~ /\S+:.+;/ )
						{
							my $note = $ref->{$dbEST_libid}->{$key};
							while ( $note =~ s/(\S+):\s([^;]+);\s//g ) {
								my $newtag   = $1;
								my $newvalue = $2;
								if ( $newtag && $newvalue ) {
									$hash_lib{$dbEST_libid}{'tags'}{$newtag} =
									  $newvalue;
								}
							}
							$hash_lib{$dbEST_libid}{'tags'}{$key} = $note;
						}
						unless (    $key eq 'note'
								 || $key eq 'clone'
								 || $key eq 'clone_lib'
								 || $key eq 'ncbi_taxid'
								 || $key eq 'dbEST_libid'
								 || $key eq 'organism'
								 || $key eq 'lineage'
								 || $key eq 'mol_type' )
						{
							$hash_lib{$dbEST_libid}{'tags'}{$key} =
							  $ref->{$dbEST_libid}->{$key}
							  if $ref->{$dbEST_libid}->{$key};
						}
					}
				}
			}
		}
	}
################
	elsif ( $type eq 'TRACEXML' ) {

		#we have: $seqtech,$user_libname,$user_orgname
		print "Checking for library in chado:";
		$lib454id = &getlibid($user_libname);
		unless ($lib454id) {
			die "\tdid not find library id for $user_libname in database.\n";
		} else {
			print "\tfound!\n";
		}
		foreach my $xmlfile (@xmls) {
			print "Processing $xmlfile\n";
			my $xs = XML::Simple->new(
									   NormaliseSpace => 2,
									   KeyAttr        => ['trace'],
									   ForceArray     => [],
			);
			my $ref = $xs->XMLin($xmlfile)->{'trace'};
			print "Parsing...\n";
            $hash_est{$ref->{'trace_name'}}{'libid'}= $lib454id;
#			foreach my $id ( %{$ref} ) {
#				$hash_est{$id}{'libid'} = $lib454id;
#			}
		}
	}
######################
	elsif ( $type eq 'SRA' ) {
		my $xml_type;
		foreach my $xmlfile (@xmls) {
			if ( $xmlfile =~ /est\./ || $xmlfile =~ /454\./ ) {
				$xml_type = "est";
			} elsif ( $xmlfile =~ /lib\./ ) {
				$xml_type = "lib";
			} else {
				warn "I don't know what type this file is : $xmlfile\n";
				next;
			}
			print "Processing $xmlfile as $xml_type\n";
			if ( $xml_type eq "lib" ) {
				my $xs =
				  XML::Simple->new(
									NormaliseSpace => 2,
									KeyAttr        => [],
									ForceArray     => [],
				  );
				my $ref = $xs->XMLin($xmlfile);    #array of hashes
				print "Parsing...\n";
				foreach my $experiment ( @{ $ref->{'EXPERIMENT_PACKAGE'} } ) {
					my $id = $experiment->{'EXPERIMENT'}->{'accession'};
					my $note =
					  $experiment->{'EXPERIMENT'}->{'DESIGN'}
					  ->{'LIBRARY_DESCRIPTOR'}
					  ->{'LIBRARY_CONSTRUCTION_PROTOCOL'};
					my $name =
					  $experiment->{'EXPERIMENT'}->{'DESIGN'}
					  ->{'LIBRARY_DESCRIPTOR'}->{'LIBRARY_NAME'};
					my $ncbi_taxid =
					  $experiment->{'SAMPLE'}->{'SAMPLE_NAME'}->{'TAXON_ID'};
					$hash_lib{$id}{'tags'}{'note'} = $note;
					$hash_lib{$id}{'libname'}      = $name;
					$hash_lib{$id}{'ncbi_taxid'}   = $ncbi_taxid;
					$hash_lib{$id}{'orgname'} =
					  $experiment->{'SAMPLE'}->{'SAMPLE_NAME'}->{'COMMON_NAME'
					  }; # only needed if organism's ncbi taxid is not stored in chado
				}
			} elsif ( $xml_type eq "est" ) {
				unless ($user_libname) {
					die
"I'm sorry, but you need to specify the library name as it is/will be in the SRA\n";
				}
				unless ($dsn) {
					die
"You need to give a database connection -dsn in order to properly parse the SRA EST XMLs\n";
				}
				$lib454id = &getlibid($user_libname);
				unless ($lib454id) {
					die
"Did not find library id for $user_libname in database.\n";
				}
				my $xs =
				  XML::Simple->new(
									NormaliseSpace => 2,
									KeyAttr        => [],
									ForceArray     => 0,
				  );
				my $ref = $xs->XMLin($xmlfile)->{'trace'};
				print "Parsing...\n";
				foreach my $key ( @{$ref} ) {
					my $id = $key->{'trace_name'};
					$hash_est{$id}{'libid'}   = $lib454id;
					$hash_est{$id}{'genbank'} = $id
					  ; # actually used by the database internally. consider renaming
				}
			}
		}
	}
}
if ($dbh) {
	&create_taxonomy_cv();
	print "\t in chado\n";
	&store_in_chado();
	&set_phylogeny() if $request_phylogeny;
}
&disconnectdb($dbh) if $dbh;
##############################################################
sub getlibid ($) {
	my $libname                = shift;
	my $select_libname_sql     = "SELECT uniquename from library where name=?";
	if ($libname=~/^\d+$/){
		print " (integer: assuming library name $libname is an internal database id) ";
		$select_libname_sql    = "SELECT uniquename from library where library_id=?";
	}
	my $select_libname_prepare = $dbh->prepare($select_libname_sql);
	$select_libname_prepare->execute($libname);
	my ($id) = $select_libname_prepare->fetchrow_array;
	if (!$id){die "Library $libname was not found in the database\n";}
	if (!$user_orgname){
		my $select_org_sql="SELECT genus||' '||species from organism join library on organism.organism_id=library.organism_id where library.uniquename=?";
		my $select_org_prepare=$dbh->prepare($select_org_sql);
		$select_org_prepare->execute($id);
		($user_orgname) = $select_org_prepare->fetchrow_array;
	}
	# Now make sure that it exists as a NCBI TAXID
	#my ($ncbi_taxid, $taxgrp, $xmlfilename)= &get_ncbi_taxid_lineage($user_orgname);
	my ( $ncbi_taxid, $species, $ranks_ref) = &get_ncbi_taxid($user_orgname);
	
	&insert_organism_chado( $ncbi_taxid, $user_orgname );
	return ($id);
}

sub insert_organism_chado($$) {
	my $ncbi_taxid             = shift;
	my $org_name_from_external = shift;
	my $species_latin;
	my $dbxref_id;
	my $get_organismname_sql =
"SELECT description,dbxref_id from dbxref where accession=? and db_id =(SELECT db_id from db where name='NCBI_TAXONOMY');";
	my $get_organismname_prepare = $dbh->prepare($get_organismname_sql);
	my ( $abbreviation, $genus, $species, $organism_id, $id );
	my $check_organism_sql =
	  "SELECT organism_id from ORGANISM where genus=? AND species=?";
	my $check_organism_prepare = $dbh->prepare($check_organism_sql);
	my $insert_organism_sql =
	  "INSERT INTO ORGANISM (genus,species,abbreviation) VALUES (?,?,?)";
	my $insert_organism_prepare = $dbh->prepare($insert_organism_sql);
	my $check_organismdbxref_sql =
"SELECT organism_dbxref_id from organism_dbxref where organism_id=? AND dbxref_id=?";
	my $check_organismdbxref_prepare = $dbh->prepare($check_organismdbxref_sql);
	my $insert_organismdbxref_sql =
	  "INSERT INTO organism_dbxref (organism_id,dbxref_id) VALUES (?,?)";
	my $insert_organismdbxref_prepare =
	  $dbh->prepare($insert_organismdbxref_sql);
	$get_organismname_prepare->execute($ncbi_taxid);
	( $species_latin, $dbxref_id ) = $get_organismname_prepare->fetchrow_array;

	# if it doesn't exist, then we have to add it from the hash_lib
	if ( !$species_latin && !$dbxref_id ) {
		$species_latin = $org_name_from_external;
		#my $ref;
		if ( !$species_latin ) {
			#my ( $x, $y, $xml_file ) = &get_ncbi_taxid_lineage($ncbi_taxid);
			 my ( $taxid, $species, $ranks_ref) = &get_ncbi_taxid_lineage($ncbi_taxid);
			 $species_latin=$species;
			#( $ref, $species_latin ) = &parse_ncbi_tax_xml($xml_file);
		}
		#undef($ref);
		if ( !$species_latin ) { return; }
		&get_set_dbxref( $ncbi_taxid, 'NCBI_TAXONOMY', $species_latin );
	}
	# species latin will never be set if dbxref_id was not.
	 elsif ( $dbxref_id && !$species_latin ) {
		$species_latin = $org_name_from_external;
		$dbh->do(
"UPDATE dbxref SET description='$species_latin' where db_id=(SELECT db_id from db where name='NCBI_TAXONOMY') AND accession=$ncbi_taxid"
		);
	}
	$get_organismname_prepare->execute($ncbi_taxid);
	( $species_latin, $dbxref_id ) = $get_organismname_prepare->fetchrow_array;
	if ( $species_latin =~ /^([A-Z]{1}\w+)\s(.+)\s*$/ ) {
		$genus        = $1;
		$species      = $2;
		$abbreviation = substr( $genus, 0, 1 ) . '.' . $species;
	}
	if ( $abbreviation && $genus && $species ) {
		$check_organism_prepare->execute( $genus, $species );
		($organism_id) = $check_organism_prepare->fetchrow_array;
		unless ($organism_id) {
			$insert_organism_prepare->execute( $genus, $species,
											   $abbreviation );
			$check_organism_prepare->execute( $genus, $species );
			($organism_id) = $check_organism_prepare->fetchrow_array;
		}
	}
	if ( $organism_id && $dbxref_id ) {
		$check_organismdbxref_prepare->execute( $organism_id, $dbxref_id );
		($id) = $check_organismdbxref_prepare->fetchrow_array;
		unless ($id) {
			$insert_organismdbxref_prepare->execute( $organism_id, $dbxref_id );
			$check_organismdbxref_prepare->execute( $organism_id, $dbxref_id );
			($id) = $check_organismdbxref_prepare->fetchrow_array;
		}
	}
	unless ( $organism_id && $dbxref_id && $id ) {
		warn
"Something failed to insert (org_id=$organism_id dbxref_id=$dbxref_id organism_dbxref=$id)\n";
	}
	$check_organismdbxref_prepare->finish();
	$insert_organismdbxref_prepare->finish();
	$insert_organism_prepare->finish();
	$check_organism_prepare->finish();
	$get_organismname_prepare->finish();
	return ($organism_id);
}

sub connecttodb ($) {

# Reading the database connection settings from a seperate file.
# E.g. dbi:Pg:dbname=InsectaCentral;host=localhost;port=15432;user=postgres;password=123321
	my $dsn = shift || die("Need a DSN file or text\n");
	my $dbConnectionSettings;

	# make it one line
	if ( -s $dsn ) {
		open( IN, $dsn );
		$dsn = "";
		while ( my $line = <IN> ) {
			unless ( $line =~ /^\s*$/ ) { chomp($line); $dsn .= $line; }
		}
		close(IN);
		if ( $dsn =~ /^\s*(\S+)\s*$/ ) { $dbConnectionSettings = $1; }
	} else {
		$dbConnectionSettings = $dsn;
	}
	my $dbh = DBI->connect($dbConnectionSettings);
	die("Error: Unable to connect to database") if ( !$dbh );
	my $dbConnectionSettings_text = $dbConnectionSettings;
	$dbConnectionSettings_text =~ s/password.+/password=XXXX/;
	print "Connected to $dbConnectionSettings_text\n";
	return ($dbh);
}

sub disconnectdb ($) {
	my $dbh = shift;

	#Clean up
	$dbh->disconnect();
}

sub get_ncbi_dbest() {
	$counter = int(0);
	my $libs = scalar( keys %get_libs );
	$timer->attr( min => 0, max => $libs );
	$timer->restart;
	print "Fetching info for $libs libraries\n";
	foreach my $lib_dbest_id ( keys %get_libs ) {
		$counter++;
		if ( $counter =~ /0$/ ) {
			print $timer->report( "eta: %E min, %40b %p\r", $counter );
		}
		my $genbank = $get_libs{$lib_dbest_id};

		# if it has been skipped.
		if ( $skipped_libs{$lib_dbest_id} || !$genbank ) {
			next;
		}
		my $xmlfilename = 'ncbi_ests/' . $genbank . '.xml';
		my $link        = $xmlurl . $genbank;
		unless ( -s $xmlfilename ) {
			system("wget -q \'$link\' -O $xmlfilename");

			#allow for IO to catch up & avoiding hammering NCBI
			sleep(1);
		}
		if ( -s $xmlfilename ) {
			my $grep = `grep -m 1 "$lineage_user" $xmlfilename`;
			if ( !$grep ) {
				$skipped_libs{$lib_dbest_id} = 1;
				next;
			}
			my $xs = XML::Simple->new(
									   NormaliseSpace => 2,
									   KeyAttr        => [''],
									   ForceArray     => ['GBQualifier'],
			);
			my $ref = $xs->XMLin($xmlfilename);   # all data
			                                      # check if it's really an EST.
			my $division_check = $ref->{'GBSeq'}->{'GBSeq_division'};
			unless ( $division_check && $division_check eq "EST" ) {
				next;
			}
			my $lineage = $ref->{'GBSeq'}->{'GBSeq_taxonomy'};
			unless ( $lineage && $lineage =~ /$lineage_user/ ) {
				$skipped_libs{$lib_dbest_id} = 1;
				next;
			}
			my $organism = $ref->{'GBSeq'}->{'GBSeq_organism'};

			# allow for multiple identical names
			my @qualifiers =
			  @{ $ref->{'GBSeq'}->{'GBSeq_feature-table'}->{'GBFeature'}
				  ->{'GBFeature_quals'}->{'GBQualifier'} };
			my $ncbi_taxid;
			$hash_lib{$lib_dbest_id}{'dbEST_libid'} = $lib_dbest_id;
			$hash_lib{$lib_dbest_id}{'lineage'}     = $lineage;
			$hash_lib{$lib_dbest_id}{'organism'}    = $organism;
			foreach my $qualifier (@qualifiers) {
				my $name  = $qualifier->{'GBQualifier_name'};
				my $value = $qualifier->{'GBQualifier_value'};
				if ( $name eq 'db_xref' && $value =~ /taxon:(\d+)/ ) {
					$ncbi_taxid = $1;
					$hash_lib{$lib_dbest_id}{'ncbi_taxid'} = $ncbi_taxid;
				} elsif ( $name eq 'db_xref' ) {
					next;
				} else {
					$hash_lib{$lib_dbest_id}{$name} = $value;
				}
			}
			if ( !$hash_lib{$lib_dbest_id}{'ncbi_taxid'} ) {
				delete $hash_lib{$lib_dbest_id};
				next;
			}

# ok  all criteria have been passed, so copy over to $hash_est and empty temp hash
		}    # end if -s xml
	}
	print "\n";
}

sub _obsolete_get_ncbi_dbest() {

# we need to reduce memory here... process every est???
# we are using hash_est_temp. if it passes, we copy it to hash_est, otherwise we delete hash_est_temp
# actually it has only one key
	foreach my $dbest_id ( keys %hash_est_temp ) {
		my $lib_dbest_id = $hash_est_temp{$dbest_id}{'dbEST_libid'};

		# if already processed and lib passed our criteria.
		if ( $hash_lib{$lib_dbest_id} ) {
			$hash_est{$dbest_id} = $hash_est_temp{$dbest_id};
			next;
		}

		# if it has been skipped.
		elsif ( $skipped_libs{$lib_dbest_id} ) {
			delete( $hash_est_temp{$dbest_id} );
			next;
		}
		my $genbank     = $hash_est_temp{$dbest_id}{'genbank'};
		my $xmlfilename = 'ncbi_ests/' . $genbank . '.xml';
		my $link        = $xmlurl . $genbank;
		unless ( -s $xmlfilename ) {
			system("wget -q \'$link\' -O $xmlfilename");

			#allow for IO to catch up & avoiding hammering NCBI
			sleep(1);
		}
		if ( -s $xmlfilename ) {
			my $grep = `grep -m 1 "$lineage_user" $xmlfilename`;
			unless ($grep) {
				$skipped_libs{$lib_dbest_id} = 1;
				next;
			}
			my $xs = XML::Simple->new(
									   NormaliseSpace => 2,
									   KeyAttr        => [''],
									   ForceArray     => ['GBQualifier'],
			);
			my $ref = $xs->XMLin($xmlfilename);   # all data
			                                      # check if it's really an EST.
			my $division_check = $ref->{'GBSeq'}->{'GBSeq_division'};
			unless ( $division_check && $division_check eq "EST" ) {
				next;
			}
			my $lineage = $ref->{'GBSeq'}->{'GBSeq_taxonomy'};
			unless ( $lineage && $lineage =~ /$lineage_user/ ) {
				$skipped_libs{$lib_dbest_id} = 1;
				next;
			}

			# allow for multiple identical names
			my @qualifiers =
			  @{ $ref->{'GBSeq'}->{'GBSeq_feature-table'}->{'GBFeature'}
				  ->{'GBFeature_quals'}->{'GBQualifier'} };

# not done because for this we will need to download the xml for every est and there, really, is no benefit for this info
#my $date_create = $ref->{'GBSeq'}->{'GBSeq_create-date'};
#my @other_ids =  $ref->{'GBSeq'}->{'GBSeq_other-seqids'}->{'GBSeqid'};
			my $ncbi_taxid;
			$hash_lib{$lib_dbest_id} = {
										 'dbEST_libid' => $lib_dbest_id,
										 'lineage'     => $lineage
			};
			foreach my $qualifier (@qualifiers) {
				if (    $qualifier->{'GBQualifier_name'} eq 'db_xref'
					 && $qualifier->{'GBQualifier_value'} =~ /taxon:(\d+)/ )
				{
					$ncbi_taxid = $1;
					$hash_lib{$lib_dbest_id}{'ncbi_taxid'} = $ncbi_taxid;
				} elsif ( $qualifier->{'GBQualifier_name'} eq 'db_xref' ) {
					next;
				} else {
					$hash_lib{$lib_dbest_id}
					  { $qualifier->{'GBQualifier_name'} } =
					  $qualifier->{'GBQualifier_value'};
				}

#debug print "my name is ".$qualifier->{'GBQualifier_name'}." and my value is ".$qualifier->{'GBQualifier_value'}."\n";
			}
			unless ($ncbi_taxid) {
				warn(
"Failed to find an ncbi taxonomy id or lineage in XML $xmlfilename\n" );
				delete $hash_lib{$lib_dbest_id};
				next;
			}

# ok  all criteria have been passed, so copy over to $hash_est and empty temp hash
			$hash_est{$dbest_id} = $hash_est_temp{$dbest_id};
		}    # end if -s xml
		delete( $hash_est_temp{$dbest_id} );
	}    # end each est
}

sub get_ncbi_taxid_lineage ($) {
	my $species_latin = shift;
	my $ncbi_taxid=get_ncbi_taxid($species_latin);
	my $taxgrp = "Unknown";
	my %ranks;
	if ($ncbi_taxid) {
	   
	   my $taxon = $taxondb->get_taxon(-taxonid => $ncbi_taxid)  || die "Cannot find your query\n";;
	   # we want class, order, family, genus, species, common name
        my ($class,$order,$family,$common)=('unknown','unknown','unknown','unknown');
        my $ncbi_rank=$taxon->rank();
        $common=$taxon->common_names();
        my $species = $taxon->scientific_name();
        if (!$species){
            die "Cannot find a scientific name for your query\n";
            return;
        }
        $species =~s/^(\w+)\s//;
        my $sort_order=1;
        $ranks{$species}{'rank'}=$ncbi_rank;
        $ranks{$species}{'ncbi'}=$ncbi_taxid;
        $ranks{$species}{'sort_order'}=$sort_order;
        my $t=$taxon->ancestor();
        
        while ($t=$t->ancestor()){
        	$sort_order++;
            my $rank=$t->rank();
            my $name=$t->scientific_name();
            my $taxid=$t->ncbi_taxid();
            $ranks{$name}{'rank'}=$rank;
            $ranks{$name}{'ncbi'}=$taxid;
            $ranks{$name}{'sort_order'}=$sort_order;
        }
	   # Now we have to reverse the sort_order. sort desc
	   my %ranks_sorted;
	   my $new_sort=1;
	   foreach my $rank (sort {$ranks{$b}{'sort_order'}<=>$ranks{$a}{'sort_order'}} keys %ranks){
	   	 $ranks_sorted{$rank}=$ranks{$rank};
	   	 $ranks_sorted{$rank}{'sort_order'}=$new_sort;
	   	 $new_sort++;
	   }
	   undef(%ranks);
	   return ( $ncbi_taxid, $species, \%ranks_sorted);	
	}
}

sub clean_species_name($) {
	my $species_latin = shift;
	if ( $species_latin =~ /^\d+$/ ) { return ($species_latin); }
	else {

		#		$species_latin =~ s/[\/\\]//g;				# slashes
		#		$species_latin =~ s/[\'\"]//g;				# quotes
		$species_latin =~ s/\s+s?u?b?sp\.?.*$//g
		  ;    # subsp sp until end of line. dot is not required
		$species_latin =~ s/\s+str\.?.*$//g;      # strain str until end of line
		$species_latin =~ s/\s+serovar\.?.*$//g;  # serovar until end of line
		$species_latin =~
		  s/\s*\(.+\)\s*/ /g;    # remove brackets and their contents
		$species_latin =~ s/\s+[a-z]+\.\s+/ /g;    # remove abbreviations
		$species_latin =~ s/\W/ /g;                # non chars
		$species_latin =~ s/\s{2,}/ /g;            # extra whitespace
		$species_latin =~ s/^\s+|\s+$//;           # remove end whitespace
		return ($species_latin);
	}
}

sub create_taxonomy_cv () {

#TODO we now have to make sure that the taxonomy ranks have a parent-child relationship proper_part_of
	my $relationship_term = 'proper_part_of';
	my $relationship_cvterm =
	  &get_set_cv_check_dbxref( $relationship_term, 'relationship' );
	my @ranks = (
				  'no rank',       'superkingdom',
				  'kingdom',       'subkingdom',
				  'phylum',        'subphylum',
				  'superclass',    'class',
				  'subclass',      'infraclass',
				  'superorder',    'order',
				  'suborder',      'infraorder',
				  'parvorder',     'superfamily',
				  'family',        'subfamily',
				  'tribe',         'subtribe',
				  'genus',         'subgenus',
				  'species group', 'species subgroup',
				  'species',       'subspecies'
	);
	my $relationship_insert_sql =
"INSERT INTO cvterm_relationship (type_id,subject_id,object_id) VALUES ('$relationship_cvterm',?,?)";
	my $relationship_insert_prepare = $dbh->prepare($relationship_insert_sql);
	my $relationship_check_sql =
"SELECT cvterm_relationship_id from cvterm_relationship where type_id='$relationship_cvterm' AND subject_id=? AND object_id=? ";
	my $relationship_check_prepare = $dbh->prepare($relationship_check_sql);

	foreach my $rank (@ranks) {
		&get_set_cv_check_dbxref( $rank, 'NCBI_TAXONOMY_RANKS', \%rank_cache );
	}
	for ( my $i = 2 ; $i < scalar(@ranks) ; $i++ ) {
		my $subject_id = $rank_cache{ $ranks[$i] };
		my $object_id  = $rank_cache{ $ranks[ $i - 1 ] };
		my $existing =
		  $relationship_check_prepare->execute( $subject_id, $object_id );
		if ( !$existing || $existing == 0 ) {
			$relationship_insert_prepare->execute( $subject_id, $object_id );
		}
	}
	$relationship_insert_prepare->finish();
	$relationship_check_prepare->finish();
}

sub get_set_cv_check_dbxref($$) {
	my $cvtermname = shift || die;    #name in cvterm -> will mimic db.accession
	my $cvname     = shift || die;    # name in cv -> will mimic db.name
	                                  # HACK: stupid dbest.
	if ( $cvname eq "Isolate" ) { $cvname = 'isolate'; }
	my $cache_hash_ref = shift;
	my $termid_sql =
"select cvterm_id from cvterm where cv_id=(select cv_id from cv where name=?) and name=?";
	my $termid_prepare = $dbh->prepare($termid_sql);
	$termid_prepare->execute( $cvname, $cvtermname );
	my $insert_dbxrefid_sql =
	  "INSERT into dbxref (db_id,accession) VALUES (?,?)";
	my $insert_dbxrefid_prepare = $dbh->prepare($insert_dbxrefid_sql);
	my $insert_termid_sql =
	  "INSERT into cvterm (cv_id,dbxref_id,name) VALUES (?,?,?)";
	my $insert_termid_prepare = $dbh->prepare($insert_termid_sql);

	# to return:
	my ($termid) = $termid_prepare->fetchrow_array;
	if ( !$termid ) {
		my $check_cv_sql     = "select cv_id from cv where name=?";
		my $check_db_sql     = "select db_id from db where name=?";
		my $check_cv_prepare = $dbh->prepare($check_cv_sql);
		my $check_db_prepare = $dbh->prepare($check_db_sql);
		$check_cv_prepare->execute($cvname);
		my ($type_cvid) = $check_cv_prepare->fetchrow_array;
		$check_db_prepare->execute($cvname);
		my ($type_dbid) = $check_db_prepare->fetchrow_array;

		if ( !$type_cvid ) {
			$dbh->do("INSERT into cv (name) VALUES ('$cvname')");
			$check_cv_prepare->execute();
			($type_cvid) = $check_cv_prepare->fetchrow_array;
		}
		if ( !$type_dbid ) {
			$dbh->do("INSERT into db (name) VALUES ('$cvname')");
			$check_db_prepare->execute();
			($type_dbid) = $check_db_prepare->fetchrow_array;
		}
		if ( !$type_cvid ) { die "Cannot insert into CV table $cvname\n"; }
		if ( !$type_dbid ) { die "Cannot insert into DB table $cvname\n"; }
		my $type_dbxrefid_sql =
		  "SELECT dbxref_id from dbxref where db_id=? and accession=?";
		my $type_dbxrefid_prepare = $dbh->prepare($type_dbxrefid_sql);
		$type_dbxrefid_prepare->execute( $type_dbid, $cvtermname );
		my ($type_dbxrefid) = $type_dbxrefid_prepare->fetchrow_array;
		if ( !$type_dbxrefid ) {
			$insert_dbxrefid_prepare->execute( $type_dbid, $cvtermname );
			$type_dbxrefid_prepare->execute( $type_dbid, $cvtermname );
			($type_dbxrefid) = $type_dbxrefid_prepare->fetchrow_array;
		}
		if ( !$type_dbxrefid ) {
			die "Cannot insert $cvtermname in dbxref table";
		}
		$insert_termid_prepare->execute( $type_cvid, $type_dbxrefid,
										 $cvtermname );
		$termid_prepare->execute();
		($termid) = $termid_prepare->fetchrow_array;
		$check_cv_prepare->finish();
		$check_db_prepare->finish();
		$type_dbxrefid_prepare->finish();
		$termid_prepare->finish();
	}
	if ( !$termid ) { die "Cannot insert $cvtermname in cvterm table"; }
	else {
		if ($cache_hash_ref) { $cache_hash_ref->{$cvtermname} = $termid; }
		return ($termid);
	}
}

sub get_set_dbxref($$$) {
	my $dbaccession = shift || die;    # db.accession
	my $dbname      = shift || die;    # db.name
	my $dbdescription = shift;         # optional
	                                   # to return
	my $type_dbxrefid;
	my $check_db_sql     = "select db_id from db where name='$dbname'";
	my $check_db_prepare = $dbh->prepare($check_db_sql);
	my $insert_dbxref_sql1 =
	  "INSERT into dbxref (db_id,accession,description) VALUES (?,?,?)";
	my $insert_dbxref_sql2 =
	  "INSERT into dbxref (db_id,accession) VALUES (?,?)";
	my $insert_dbxref_prepare1 = $dbh->prepare($insert_dbxref_sql1);
	my $insert_dbxref_prepare2 = $dbh->prepare($insert_dbxref_sql2);
	$check_db_prepare->execute();
	my ($type_dbid) = $check_db_prepare->fetchrow_array;

	if ( !$type_dbid ) {
		$dbh->do("INSERT into db (name) VALUES ('$dbname')");
		$check_db_prepare->execute();
		($type_dbid) = $check_db_prepare->fetchrow_array;
	}
	if ( !$type_dbid ) { die "Cannot insert into DB table $dbname\n"; }
	my $type_dbxrefid_sql =
"SELECT dbxref_id from dbxref where db_id='$type_dbid' and accession='$dbaccession'";
	my $type_dbxrefid_prepare = $dbh->prepare($type_dbxrefid_sql);
	$type_dbxrefid_prepare->execute;
	($type_dbxrefid) = $type_dbxrefid_prepare->fetchrow_array;
	if ( !$type_dbxrefid ) {
		if ($dbdescription) {
			$insert_dbxref_prepare1->execute( $type_dbid, $dbaccession,
											  $dbdescription );
		} else {
			$insert_dbxref_prepare2->execute( $type_dbid, $dbaccession );
		}
		$type_dbxrefid_prepare->execute;
		($type_dbxrefid) = $type_dbxrefid_prepare->fetchrow_array;
	}
	$check_db_prepare->finish();
	$type_dbxrefid_prepare->finish();
	if ( !$type_dbxrefid ) { die "Cannot insert $dbaccession in dbxref table"; }
	else {
		return ($type_dbxrefid);
	}
}

sub get_set_pub($$$) {
	my $uniquename       = shift || die;
	my $title            = shift;
	my $term_id          = shift || die;
	my $pub_get          = "SELECT pub_id FROM pub where uniquename=?";
	my $prepared_pub_get = $dbh->prepare($pub_get);
	$prepared_pub_get->execute($uniquename);
	my ($pub_id) = $prepared_pub_get->fetchrow_array;
	if ( !$pub_id ) {
		my $code = &get_set_cv_check_dbxref( $term_id, 'publications' );
		$dbh->do(
"INSERT INTO pub (title,uniquename,type_id) VALUES ('$title','$uniquename','$code');"
		);
		$prepared_pub_get->execute($uniquename);
		($pub_id) = $prepared_pub_get->fetchrow_array;
		if ( !$pub_id ) {
			die("Still cannot get the publication ID for the IEA code.\n");
		}
	}
	return ($pub_id);
}

sub store_in_chado() {

# use globals of %hash_est && $hash_lib
# hash lib goes into library and library_dbxref
# hash_est goes into library_feature
# issues we cannot solve: the values of tissue, sex, stage are not CVs
# we are going to store them into a new CV so that we can later curate them (in due time).
# during curation, we will keep the original dbxref accession and only curate the cvterm
# we do not want to start using the publication record for this library because dbest is very unclean
# so we create a general pub record for all dbest libraries.
	my $pub_id_for_cvterm =
	  &get_set_pub( 'dbEST record', 'dbEST record', 'database record' );

#TODO the lineage (or phylum if XML was not used) values can be stored as phylonodes but wouldn't that be too complicated?
# get organism_id using organism_dbxref and ncbi_taxid
# library type will always be "cDNA library" 	#if it doesn't exist. add it
	my $get_organism_id_sql =
	    "SELECT organism_id from organism_dbxref where dbxref_id="
	  . "(SELECT dbxref_id from dbxref where accession=? and db_id="
	  . "(SELECT db_id from db where name='NCBI_TAXONOMY'))";
	my $get_organism_id_prepare = $dbh->prepare_cached($get_organism_id_sql);
	my ($lib_termid) =
	  &get_set_cv_check_dbxref( 'cDNA library', 'Library types' );

	# the uniquename will be the dbest dbxref as we are certain this is unique.
	my $insert_lib_sql =
"INSERT INTO library (type_id,organism_id,name,uniquename) VALUES ('$lib_termid',?,?,?)";
	my $insert_lib_prepare = $dbh->prepare($insert_lib_sql);
	my $insert_libdbxref_sql =
	  "INSERT INTO library_dbxref (library_id,dbxref_id) VALUES (?,?)";
	my $check_libdbxref_sql =
"SELECT library_dbxref_id from library_dbxref where library_id=? AND dbxref_id=?";
	my $insert_libdbxref_prepare = $dbh->prepare($insert_libdbxref_sql);
	my $check_libdbxref_prepare  = $dbh->prepare($check_libdbxref_sql);
	my $select_lib_sql = "SELECT library_id from library where uniquename=?";
	my $select_lib_prepare = $dbh->prepare($select_lib_sql);
	my $insert_libcvterm_sql =
"INSERT INTO library_cvterm (pub_id,library_id,cvterm_id) VALUES ('$pub_id_for_cvterm',?,?)";
	my $check_libcvterm_sql =
"select library_cvterm_id from library_cvterm where pub_id='$pub_id_for_cvterm' AND library_id=? AND cvterm_id=?";
	my $insert_libcvterm_prepare = $dbh->prepare($insert_libcvterm_sql);
	my $check_libcvterm_prepare  = $dbh->prepare($check_libcvterm_sql);
	my $insert_libprop_sql =
"INSERT INTO libraryprop (library_id,type_id,value,rank) VALUES (?,?,?,?)";
	my $check_libprop_sql =
"SELECT libraryprop_id from libraryprop where library_id=? AND type_id=? AND value=?";
	my $rank_libprop_sql =
	  "SELECT rank from libraryprop where library_id=? AND type_id=? ";
	my $insert_libprop_prepare = $dbh->prepare($insert_libprop_sql);
	my $check_libprop_prepare  = $dbh->prepare($check_libprop_sql);
	my $rank_libprop_prepare   = $dbh->prepare($rank_libprop_sql);
	my $select_est_featureid_sql =
	  "SELECT feature_id from feature where uniquename=?";
	my $select_est_featureid_prepare = $dbh->prepare($select_est_featureid_sql);
	my $insert_libfeature_sql =
	  "INSERT INTO library_feature (library_id,feature_id) VALUES (?,?)";
	my $insert_libfeature_prepare = $dbh->prepare($insert_libfeature_sql);
	my $check_libfeature_sql =
"SELECT library_feature_id from library_feature where library_id=? AND feature_id=?";
	my $check_libfeature_prepare = $dbh->prepare($check_libfeature_sql);
	my $check_feature_dbxref_sql =
"SELECT feature_dbxref_id from feature_dbxref where feature_id=? and dbxref_id=?";
	my $check_feature_dbxref_prepare = $dbh->prepare($check_feature_dbxref_sql);
	my $insert_feature_dbxref_sql =
	  "INSERT INTO feature_dbxref (feature_id,dbxref_id) VALUES (?,?)";
	my $insert_feature_dbxref_prepare =
	  $dbh->prepare($insert_feature_dbxref_sql);

	### LIBRARY
	if ( scalar( keys %hash_lib ) > 0 ) {
		foreach my $dbest_lib_id ( keys %hash_lib ) {
			my $ncbi_taxid = $hash_lib{$dbest_lib_id}{'ncbi_taxid'}
			  || die("Received no ncbi tax id for library $dbest_lib_id");
			my $name = $dbest_lib_id;
			$name = "dbEST id $dbest_lib_id" if ( $type eq 'dbEST' );
			$name = "SRA id $dbest_lib_id"   if ( $type eq 'SRA' );
			$name = $hash_lib{$dbest_lib_id}{'libname'}
			  if $hash_lib{$dbest_lib_id}{'libname'};
			$get_organism_id_prepare->execute($ncbi_taxid);
			my ($organism_id) = $get_organism_id_prepare->fetchrow_array;
			if ( !$organism_id ) {
				my $organism_name = $hash_lib{$dbest_lib_id}{'orgname'};
				if ($organism_name) {
					$organism_id =&insert_organism_chado( $ncbi_taxid, $organism_name );
				}
			}
			if ( !$organism_id ) {
				warn "Insert of organism id failed for $name and $ncbi_taxid TaxID.\n";
				die "Exiting: This error could point towards a serious error somewhere in your input/database.\n";
			}
			$select_lib_prepare->execute($dbest_lib_id);
			my ($lib_id) = $select_lib_prepare->fetchrow_array;
			if ( !$lib_id ) {
				$insert_lib_prepare->execute( $organism_id, $name,
											  $dbest_lib_id );
				$select_lib_prepare->execute($dbest_lib_id);
				($lib_id) = $select_lib_prepare->fetchrow_array;
			}
			if ( !$lib_id ) {
				die
"Can't find library and insertion failed: Library $name ($organism_id,$name,$dbest_lib_id).\n";
			}

			#dbxref link
			my $dbxref_id = &get_set_dbxref( $dbest_lib_id, 'dbEST' );
			if ( $lib_id && $dbxref_id ) {
				my $existing =
				  $check_libdbxref_prepare->execute( $lib_id, $dbxref_id );
				if ( !$existing ) {
					$insert_libdbxref_prepare->execute( $lib_id, $dbxref_id );
				}
			}
#cvterms
# if a sequencing technology has been applied to this library. for dbEST it's for sure at least Sanger
# notes go to libraryprop as Note type_id
			my $tech_cvterm =
			  &get_set_cv_check_dbxref( $seqtech, 'sequencing_technology' );
			$dbh->begin_work;
			foreach my $term ( keys %{ $hash_lib{$dbest_lib_id}{'tags'} } ) {
				my $tag = $hash_lib{$dbest_lib_id}{'tags'}{$term};
				if ( !$tag ) { next; }
				if ( $tag =~ /'/ ) { $tag = $dbh->quote($tag); }
				if ( $term eq 'note' ) {

					# is this cv ok?
					my $cvterm =
					  &get_set_cv_check_dbxref( 'Note', 'feature_property' );
					my $existing =
					  $check_libprop_prepare->execute( $lib_id, $cvterm, $tag );
					if ( !$existing || $existing == 0 ) {
						$rank_libprop_prepare->execute( $lib_id, $cvterm );
						my ($rank) = $rank_libprop_prepare->fetchrow_array;
						$rank++;
						$insert_libprop_prepare->execute( $lib_id, $cvterm,
														  $tag, $rank );
					}
					next;
				}
				my $cvterm = &get_set_cv_check_dbxref( $tag, $term );
				my $existing =
				  $check_libcvterm_prepare->execute( $lib_id, $cvterm );
				if ( !$existing || $existing == 0 ) {
					$insert_libcvterm_prepare->execute( $lib_id, $cvterm );
				}
			}
			if ($tech_cvterm) {
				my $existing =
				  $check_libcvterm_prepare->execute( $lib_id, $tech_cvterm );
				if ( !$existing || $existing == 0 ) {
					$insert_libcvterm_prepare->execute( $lib_id, $tech_cvterm );
				}
			}
			$dbh->commit;
		}
	}    # END LIBRARY
	$get_organism_id_prepare->finish();
	$insert_lib_prepare->finish();
	$check_libdbxref_prepare->finish();
	$insert_libdbxref_prepare->finish();
	$insert_libcvterm_prepare->finish();
	$check_libcvterm_prepare->finish();
	
	### ESTs
	if ( scalar( keys %hash_est ) > 0 ) {
		print "Checking/adding library information for "
		  . scalar( keys %hash_est )
		  . " features...\n";
		$counter = int(0);
		$timer->attr( min => 0, max => scalar( keys %hash_est ) );
		$timer->restart;
		$dbh->begin_work;
		foreach my $est_id ( keys %hash_est ) {
			if ( scalar( keys %hash_est ) > 0 ) {
				$counter++;
				if ( $counter =~ /0000$/ ) {
					print $timer->report( "eta: %E min, %40b %p\r", $counter );
				}
			}

			my $genbank = $hash_est{$est_id}{'genbank'};
			# if we have no genbank, the est_id is the official ID (e.g SRA)
			if (!$genbank){$genbank=$est_id;}
			$select_est_featureid_prepare->execute($genbank);
			my ($feature_id) = $select_est_featureid_prepare->fetchrow_array;
			if ($feature_id) {
				my $lib = $hash_est{$est_id}{'libid'} || die;
				$select_lib_prepare->execute($lib);
				my ($lib_id) = $select_lib_prepare->fetchrow_array;
				if ( !$lib_id ) {
					die "Can't find library of est $est_id.\n";
				}
				$check_libfeature_prepare->execute( $lib_id, $feature_id );
				my ($id) = $check_libfeature_prepare->fetchrow_array;
				unless ($id) {
					$insert_libfeature_prepare->execute( $lib_id, $feature_id );
				}
				undef($id);
				if ( $type eq "dbEST" ) {
					#TODO feature synonym for genbank record (if in dbEST)
					my $feature_dbxref_id = &get_set_dbxref( $est_id, 'dbEST' );
					$check_feature_dbxref_prepare->execute( $feature_id,
														   $feature_dbxref_id );
					($id) = $check_feature_dbxref_prepare->fetchrow_array;
					unless ($id) {
						$insert_feature_dbxref_prepare->execute( $feature_id,
														   $feature_dbxref_id );
					}
				}
			}
			
		}
		print "\nCommitting transaction... this may take a while\n";
		$dbh->commit;
	}    # END EST
	$check_libfeature_prepare->finish();
	$insert_libfeature_prepare->finish();
	$select_lib_prepare->finish();
	$select_est_featureid_prepare->finish();
	print "\n";
}

sub set_phylogeny() {
	print "Updating organism phylogeny table\n";

# use the phylonode?
# we are going to have one phylotree for the whole taxonomy...?!
# after a bit of thinking and playing around I decided that the
# phylomodule is not useful for a dynamic database content, it is of no use to have more than
# phylotree but all nodes must be known before a phylonodes can be built (due to the required indexing)
# in addition, the increased complexity will simply make it far more difficult (read expensive) to process
# therefore it is best to use organismprop where the key is a cvterm and value is a dbxref entry.
# i.e. explicitly set organismprop values where the type id is the rank cvterm
# and the values is the taxon ncbi id (which is indexed in my copy...).
# then we utilize the dbxref to link taxon ncbi ids to db.
# the organismprop.rank = 0 and allow only one entry as a species can belong to only one taxonomic entity at each level
# we will not use one cvterm and multiple ranks as that way we lose the ability to insert/remove (taxonomic) ranks at will
	my $get_organisms_id_sql =
	    "SELECT organism_id,accession from organism_dbxref natural join dbxref "
	  . "where dbxref.db_id=(SELECT db_id from db where name='NCBI_TAXONOMY');";
	my $get_organisms_id_prepare = $dbh->prepare_cached($get_organisms_id_sql);
	my $insert_phylogeny_sql =
"INSERT into organismprop (rank,organism_id,type_id,value) VALUES (?,?,?,?)";
	my $insert_phylogeny_prepare = $dbh->prepare($insert_phylogeny_sql);
	$get_organisms_id_prepare->execute();
	while ( my ( $organism_id, $ncbi_taxid ) =
			$get_organisms_id_prepare->fetchrow_array )
	{
		my $delete_sql ="DELETE from organismprop where organism_id=$organism_id and type_id IN (select cvterm_id from cvterm where cv_id = (SELECT cv_id from cv where name='NCBI_TAXONOMY_RANKS'))";
		my ( $taxid, $species, $ranks_ref) = &get_ncbi_taxid_lineage($ncbi_taxid);
		#my ( $x, $y, $xml_file ) = &get_ncbi_taxid_lineage($ncbi_taxid);
		#my ( $ref, $orgname ) = &parse_ncbi_tax_xml($xml_file);
		#undef($x);
		#undef($y);
		my %ranks = %$ranks_ref;
		$dbh->do($delete_sql);
		foreach my $name ( keys %ranks ) {
			my $taxon_rank     = $ranks{$name}{'rank'};
			my $rank_cvterm_id = $rank_cache{$taxon_rank};
			my $order          = $ranks{$name}{'sort_order'};
			my $ncbi_id        = $ranks{$name}{'ncbi'};
			my $taxon_dbxref_id =
			  &get_set_dbxref( $ncbi_id, 'NCBI_TAXONOMY', $name );
			$insert_phylogeny_prepare->execute( $order, $organism_id,
												$rank_cvterm_id, $name );
		  }
	   }

	$insert_phylogeny_prepare->finish();
	$get_organisms_id_prepare->finish();
}

sub get_ncbi_taxid ($){

    my $species_latin=shift;
    my $ncbi_taxid;
    if ( $species_latin =~ /^\d+$/ ) { 
        $ncbi_taxid = $species_latin; 
        return $ncbi_taxid;
    }
    
    #at least three letters
    elsif ( $species_latin =~ /^([A-Za-z]{2}[a-z]+)\s*/ ) {
        $ncbi_taxid = $taxondb->get_taxonid($species_latin);
        if ($ncbi_taxid){
            return $ncbi_taxid;
        }
        my $binomial =$1;
        if ($species_latin =~ /^[A-Za-z]{2}[a-z]+\s+(\w+)/ ) {
            $binomial.=' '.$1;
        }
        if ($binomial){
            $ncbi_taxid = $taxondb->get_taxonid($binomial);
            if (!$ncbi_taxid){
                warn "Could not find NCBI taxid for $species_latin.\n";
            }else{
                return $ncbi_taxid;
            }
        }
    }
} 
 
sub _obsolete_parse_ncbi_tax_xml($) {
	#no longer used
	my $xml_file = shift;
	unless ( -s $xml_file ) { return; }
	my $orgname;

	#parse xml as a tree to get taxa with ranks
	my $xs = XML::Simple->new(
							   NormaliseSpace => 2,
							   KeyAttr        => [],
							   ForceArray     => 0
	);
	my $ref         = $xs->XMLin($xml_file);          # all data
	my $order_str   = $ref->{'Taxon'}->{'Lineage'};
	my @order_array = split( ';', $order_str );
	foreach (@order_array) { $_ =~ s/^\s+|\s+$//; }
	$orgname = $ref->{'Taxon'}->{'ScientificName'};
	$ref     = $ref->{'Taxon'}->{'LineageEx'}->{'Taxon'};
	my %hash;

	for ( my $i = 0 ; $i < scalar(@$ref) ; $i++ ) {
		$hash{$i} = @$ref[$i];
	}
	my %ranks;    # to return
	              # rebuild xml to rank based.
	foreach my $key ( keys %hash ) {
		my $order;
		my $name        = $hash{$key}{'ScientificName'};
		my $taxon_taxid = $hash{$key}{'TaxId'};
		for ( my $i = 0 ; $i < scalar(@order_array) ; $i++ ) {
			if ( $name eq $order_array[$i] ) {
				$order = $i;
			}
		}
		if ( !$order ) { next; }
		my $taxon_rank = $hash{$key}{'Rank'};
		$ranks{$name}{'rank'}  = $taxon_rank;
		$ranks{$name}{'ncbi'}  = $taxon_taxid;
		$ranks{$name}{'sort_order'} = $order;
	}
	return ( \%ranks, $orgname );
}

sub store_in_xml ($) {
	my $type = shift;
	if ( $type eq "est" ) {
		my $results_est = new IO::File( ">$dbestfile" . "_parse.est.xml" );
		my %hash;
		$hash{'EST'} = \%hash_est;
		my $xml = XMLout(
						  \%hash,
						  AttrIndent => 1,
						  NoAttr     => 1,
						  KeyAttr    => ['dbEST_id'],
						  XMLDecl  => '<?xml version="1.0" encoding="utf-8"?>',
						  RootName => "dbEST_sequences",
						  OutputFile => $results_est
		);
		$results_est->close;
	}
	if ( $type eq "lib" ) {
		my $results_lib = new IO::File( ">$dbestfile" . "_parse.lib.xml" );

# how to parse this:
# my $xs = XML::Simple->new(NormaliseSpace => 2,ForceArray => 0,KeyAttr=>'id');$ref = $xs->XMLin($file);
		my %hash;
		$hash{'library'} = \%hash_lib;
		my $xml = XMLout(
						  \%hash,
						  AttrIndent => 1,
						  NoAttr     => 1,
						  KeyAttr    => ['dbEST_libid'],
						  XMLDecl  => '<?xml version="1.0" encoding="utf-8"?>',
						  RootName => "dbEST_libraries",
						  OutputFile => $results_lib
		);
		$results_lib->close;
	}
}
