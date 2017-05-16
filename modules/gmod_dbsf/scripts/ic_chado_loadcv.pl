#!/usr/bin/perl -w

#TODO physical properties
#TODO library info

#CHANGELOG
# 14Mar10: Added Iprscan subdatabases

=head1 NAME

ic_chado_loadcv.pl 

=head1 USAGE

Process annot8r GFF files to load CV terms into features and avoid using gmod_bulk_load

ic_chado_loadcv.pl [options] @[GFF or BLAST files]

	* 'dsn|database_connection:s'	=> of database to create terms from
	'cv_name:s'					=> EC KEGG GO, IPR or BLAST. 
	'db_name:s'			=> value in db table to use (defaults to cvterm above). 
	'chado|usechado' 		=> use chado instead of GFF files. Not implemented yet.
	'gff'          => File is in GFF format (defaults to yes for EC KEGG GO IPR, FALSE for BLAST)
	'csvfile:s'		=> CSV file with IprScan ids and names, in order to update the description in the dbxref table.
	'nohits'		=> For BLAST only: Don't process hits, just store that a query has been processed with a BLAST database.
    	'flat:ncbi_flatfiledir:s'   => Directory containing flatfiles of NCBI Taxonomy (def. /db/ncbi_taxonomy/) 

=head1 DESCRIPTION

 The -cv_name is optional as it can be part of the file name (if EC/KEGG/GO just before the .gff suffix).
 The NCBI Tax flatfile can be acquired from ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz


=head2 BLASTs

Databasing of BLAST files into chado does not follow the standard practice: we are not really interested in visualizing
 the alignments of hits and HSPs but rather utilize BLAST to offer IEA terms for searching. This is accomplished by 
 harvesting the description line of each hit. 
A special file needs to be given in order to speed up processing. This is very easy to produce using egrep:
  egrep -hA 1 "^Query=|^>|^Database: " blastfiles > special_file
 This special file can be from multiple files (therefore do use egrep -h).

Otherwise, you can also load a GFF file (use the -gff flagg)

=head1 AUTHORS

    Alexie Papanicolaou

    Centre for Ecology and Conservation, University of Exeter, UK
    * recipient of complaints:  alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

 This software is released under the GNU General Public License version 3 (GPLv3).
 It is provided "as is" without warranty of any kind.
 You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html. 
 Please note that incorporating the whole software or parts of its code in proprietary software 
 is prohibited under the current license.

=head1 BUGS & LIMITATIONS

 None known. Please report them.

 You will need a preloaded chado database in order to get get the correct CV terms for GO

An example of a DSN file:

dbi:Pg:dbname=chado;host=localhost;port=5433;user=chado_user;password=123

=cut

use strict;
use Text::CSV_XS;
use Bio::Taxon;
use Bio::DB::Taxonomy;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Bio::Tools::GFF;
use DBI;
use Data::Dumper;
use Time::Progress;
$| = 1;

my $flatfile_dir='/databases/ncbi_taxonomy/';

my ( $dsn, $cv_name, $usechado, $csvnames, $db_name, $noskip,$nohits,$gff_format,%toskip );
GetOptions(
	'dsn|database_connection:s' => \$dsn,     # of database to create terms from
	'cv_name:s'                 => \$cv_name, # Ec Kegg
	'db_name:s'      => \$db_name,  # value in db table to use (def. == cv_name)
	'chado|usechado' => \$usechado, # use chado instead of GFF files.
	'csvfile:s'      => \$csvnames,
	'nohits'		=> \$nohits,
	'flat|ncbi_flatfiledir:s' => \$flatfile_dir,
	'gff'          => \$gff_format,    #only useful for BLAST
	, # csv file with IprScan ids, names and other attributes, in order to update the description in the dbxref table.
	'noskipdb' => \$noskip
);
my @infiles = @ARGV;
unless ($dsn) { 
	warn "Providing a DSN is mandatory.\n";
	pod2usage; 
}
my $dbh = &connecttodb($dsn) || die("No database connection provided!\n");
unless ( @infiles || $usechado || $csvnames ) {
	&disconnectdb($dbh);
	pod2usage;
	die();
}
my $select_global_db_sql="SELECT db_id from db where name=? OR name='DB:'||? ";
my $insert_global_db_sql="INSERT INTO db (name) VALUES (?) ";
my $select_global_dbxref_sql = "SELECT dbxref_id from dbxref where accession=? AND db_id = ? ";
my $insert_global_dbxref_sql="INSERT INTO dbxref (db_id,accession) VALUES (?,?)";
my $insert_global_dbxref_prepared=$dbh->prepare($insert_global_dbxref_sql);
my $select_global_dbxref_prepared=$dbh->prepare_cached($select_global_dbxref_sql);
my $select_global_db_prepared=$dbh->prepare_cached($select_global_db_sql);
my $insert_global_db_prepared=$dbh->prepare($insert_global_db_sql);

my $taxondb = Bio::DB::Taxonomy->new(-source => 'entrez');
if ($flatfile_dir && -d $flatfile_dir && -s $flatfile_dir.'/nodes.dmp' && -s $flatfile_dir.'/names.dmp'){
  $taxondb = Bio::DB::Taxonomy->new(-source => 'flatfile',-nodesfile => $flatfile_dir.'/nodes.dmp',-namesfile => $flatfile_dir.'/names.dmp',-directory => $flatfile_dir);
}

# what we want to do is
# foreach file or for the chadodb
# identify which CV we are using (EC, GO, KEGG, IPR) and
# i) get the reference feature ID (feature's display_name) and any annotation IDs (Alias and/or dbxref)
# iia) get the cv_id and cvterm foreach (hash). store it in feature_cvterm
# iib) if using the dbxref, they should already exist in feature_dbxref. We need to ensure that a proper name exists in dbxref
#
# formally, we use the cvterm if there is an ontology (a multilevel hierachy) but in practice there are other issues to consider
# pros and cons for alias (cvterm table) and dbxref (db table)
# pros
# in cvterm the whole vocab is loaded, along with definitions
# in dbxref the links are already loaded by bulk_load as feature_dbxref -> worth replicating as cvterms?
# if cvterm table changes, then feature_cvterm must be reloaded. in dbxref, why don't have to do this.
# i guess because the cvterm contains the definitions and for GO, EC, KEGG is unlikely to change, then loaded them will be fine, especially since we can create a glossary
# for IPR there are no pre-loaded definitions or the vocabulary is huge and frequently changes, so keep them as dbxref.
# the cvterm_dbxref will replicate the dbxref that had the highest hit (UNIPROT)
# so we also have to update the dbxref table to add a description (descr field from UNIPROT?).
#
#
my %file_data;
my %failed;    # failed accession ids, so that we report failure once per accession.

print "Building HASHes\n";
if ($usechado) {
	die("Not implemented yet\n");
} else {
	my $count = int(0);
	foreach my $file (@infiles) {
		unless ( -s $file ) {
			warn("Cannot find $file! Skipping...\n");
			delete( $infiles[$count] );
			next;
		}
		my $filename = $file;
		$filename =~ s/\.gff$//;

		#validate
		unless ( -s $file ) {
			&disconnectdb($dbh);
			die("Couldn't find $file\n");
		}
		if ( !$cv_name ) {    # user didn't provide cvname
			if ( $filename =~ /EC$/ ) {
				$file_data{$file}{"cv_name"} = "EC";
				$file_data{$file}{"db_name"} = "EC" unless $db_name;
			} elsif ( $filename =~ /GO$/ ) {
				$file_data{$file}{"cv_name"} ="'biological_process','molecular_function','cellular_component'";
				$file_data{$file}{"db_name"} = "GO" unless $db_name;
			} elsif ( $filename =~ /KEGG$/ ) {
				$file_data{$file}{"cv_name"} = "KEGG_PATHWAY";
				$file_data{$file}{"db_name"} = "KEGG_PATHWAY" unless $db_name;
			} elsif (    $filename =~ /IPR$/i
					  || $filename =~ /interpro$/i
					  || $filename =~ /iprscan$/i )
			{
				$file_data{$file}{"cv_name"} = "InterPro";
				$file_data{$file}{"db_name"} = "InterPro" unless $db_name;
			} elsif ($filename =~ /blastall/i || $filename =~ /\.\w?blast\w?$/i
					  || ( $cv_name && $cv_name =~ /blast/i ) )
			{
				$file_data{$file}{"cv_name"} =
				  "BLAST";    # for symmetry, ignored really...
				$file_data{$file}{"db_name"} = "BLAST";
#TODO: physical properties databasing
#			} elsif ( $filename =~ /physprop\w$/i
#					  || ( $cv_name && $cv_name =~ /physprop/i ) )
#			{
#				$file_data{$file}{"cv_name"} ="physprop";
#				$file_data{$file}{"db_name"} ="physprop";
			} else {
				warn("Couldn't tell what type of CV you are trying to annotate with $file... Skipping.\n");
				next;
			}
		} else {

			#cv_name provided by user
			$cv_name = uc($cv_name);
			unless (
					 $cv_name    eq "EC"          # -> feature_cvterm
					 || $cv_name eq "GO"          # -> feature_cvterm
					 || $cv_name eq "KEGG"        # -> feature_cvterm
					 || $cv_name eq "InterPro"    # -> feature_dbxref
					 || $cv_name eq "BLAST"       # -> feature_dbxref
#					 || $cv_name eq "physprop"	  # -> featureprop
			  )
			{
				&disconnectdb($dbh);
				die("I don't know what kind of cvtype $cv_name is...\n");
			}
			if ( !$db_name ) {               # but no db_name
				$file_data{$file}{"db_name"} = $cv_name;
			}
		}
		$count++;
	}
	print "\nProcessing\n";
	foreach my $file (@infiles) {
		my %cvterm;      # global hash, contains cvterm id for an accession id
		my %features;    # for inserting into DB with one transaction
		
		my ( $db_name_local, $cv_name_local );
		$cv_name_local = $file_data{$file}{"cv_name"}
		  if $file_data{$file}{"cv_name"};
		$cv_name_local = $cv_name if $cv_name;
		if ( !$cv_name_local ) { warn "No cvname for $file. Skipping\n"; next; }
		$db_name_local = $file_data{$file}{"db_name"}
		  if $file_data{$file}{"db_name"};
		$db_name_local = $db_name if $db_name;
		if ( !$db_name_local ) { warn "No dbname for $file. Skipping"; next }

		#end validate
		if ( $db_name_local eq "BLAST" ) {
			&process_blast($file);
			next;
		}

		# Globals for performance reasons
		my $sql_terms;
		## EC KEGG
		if ( $db_name_local eq "EC" || $db_name_local eq "KEGG_PATHWAY" ) {
			$sql_terms =
"select cvterm_id from cvterm where cv_id IN (select cv_id from cv where name=? ) and dbxref_id IN (select dbxref_id from dbxref where accession=? and db_id =(select db_id from db where name=?))";
		}

#problem with GO and the fact that we lost the PFC term and only have the cvterm. As Chado has splitted the GO to PFC cvs we must give cv.name a LIKE search option. and as expected, it's circa 3 times slower
#       GO
		elsif ( $db_name_local eq "GO" ) {
			$sql_terms =
			    "select cvterm_id from cvterm where"
			  . " cv_id IN (select cv_id from cv where name IN ($cv_name_local) )"
			  . " and dbxref_id IN (select dbxref_id from dbxref "
			  . " where accession=? and db_id =(select db_id from db where name=?))";
		} elsif ( $db_name_local eq "InterPro" ) {
			$sql_terms = "SELECT dbxref_id from dbxref where "
			  . " accession=? AND db_id = (SELECT db_id from db where name =? )";
			  
		}
		if (!$sql_terms){
			die "I cannot understand which data you are loading...\n";
		}
		my $prepared_terms = $dbh->prepare_cached($sql_terms);
# Do we want to make use of BioPerl (robust) or shall we parse the text file manually (fast)? -> Let's go for BioPerl even if i grow to regret it
		print "Processing $file as $cv_name_local\n";
		my $gff_obj = Bio::Tools::GFF->new( -gff_version => 3, -file => $file )
		  || die();
		$gff_obj->ignore_sequence(1);

		# get Bio::SeqFeatureI
		while ( my $feature = $gff_obj->next_feature() ) {
			my $feat_name = $feature->seq_id();    # reference.
			if ( !$feat_name ) {
				&disconnectdb($dbh);
				die("Couldn't parse feature name\n");
			}
			unless ( $feat_name =~ /pep/ ) {
				warn
				  "$feat_name doesn't seem to be an est2assembly protein...\n";
			}
			#$features{$feat_name}{"dbxrefs"} = [];
			#$features{$feat_name}{"aliases"} = [];

			# get data
			my @aliases = $feature->get_tag_values('Alias')
			  if $feature->has_tag('Alias');
			my @dbxrefs = $feature->get_tag_values('Dbxref')
			  if $feature->has_tag('Dbxref');

			# alias: get cvterm id for this alias
			# only for GO EC
			if ( $db_name_local eq "GO" || $db_name_local eq "EC" ) {
				for ( my $i = 0 ; $i < @aliases ; $i++ ) {
					unless ( $cvterm{ $aliases[$i] } ) {
						$cvterm{ $aliases[$i] } =
						  &get_cvterm( $aliases[$i],   $cv_name_local,
									   $db_name_local, $prepared_terms );

						#only add those that have a cvterm
						if ( !$cvterm{ $aliases[$i] } ) {
							delete( $aliases[$i] );
						}
					}
				}
			}

			#only for KEGG and IPR
			if (    $db_name_local eq "KEGG_PATHWAY"
				 || $db_name_local eq "InterPro" ){
				for ( my $i = 0 ; $i < @dbxrefs ; $i++ ) {
					$dbxrefs[$i] =~ s/^\s*\://;
					$dbxrefs[$i] =~ /^(\w+)\:/;
					my $db = $1;
					if(!$db){next;}
                    elsif ($db eq 'IPR'){$db='IPR';}
                    elsif($db_name_local eq "KEGG_PATHWAY" && $db ne $db_name_local){
                    	next;
                    }
					unless ( $cvterm{ $dbxrefs[$i] } ) {
						$cvterm{ $dbxrefs[$i] } =
						  &get_cvterm( $dbxrefs[$i],   $cv_name_local,
									   $db, $prepared_terms );
					}if ( !$cvterm{ $dbxrefs[$i] } ) {
						#warn "No term for $dbxrefs[$i]...\n";
						 delete( $dbxrefs[$i] ); 
					}
				}
			}
            push( @{ $features{$feat_name}{"aliases"} }, @aliases );
            push( @{ $features{$feat_name}{"dbxrefs"} }, @dbxrefs );
		}
		$prepared_terms->finish;
        $insert_global_dbxref_prepared->finish;
        $select_global_db_prepared->finish;
        $insert_global_db_prepared->finish;
        $select_global_dbxref_prepared->finish; 
		# insert into feature_cvterm and/or feature_dbxref
		print "Inserting cvterms\n";
		sleep(1);
		&insert_cvterms( \%features, $db_name_local,\%cvterm );
	}    #end of file
}    #we have done the inserts
foreach my $fail ( keys %failed ) {
	warn $failed{$fail};
}
print "Done.\n";
if ($csvnames) {
	print "Inserting IPR names\n";
	&repair_dbxref('InterPro');
	if ( !@infiles && !$usechado ) { exit; }
}
print "\n";
&disconnectdb($dbh);
#########################################################################
sub build_blash_hashes_gff($) {
	my $gffblastfile = shift || die();
	print "Processing GFF file $gffblastfile\n";
	my %hash_id;   #   keeps counts: $hash_id_ref->{$id}{"hits"}{$hit}++ and $hash_id_ref->{$id}{"dbs"}{$db}++; 
	my %hash_hits; # builds: $hash_hits_ref->{$hit}{'descr'} and $hash_hits_ref->{$hit}{'db'}  
	my %skipped;   # skipped databases
	my $counter=0;
	my $gff_obj =
	  Bio::Tools::GFF->new( -gff_version => 3, -file => $gffblastfile )
	  || die();
	$gff_obj->ignore_sequence(1);

	# get Bio::SeqFeatureI
	while ( my $feature = $gff_obj->next_feature() ) {
		if ($feature->primary_tag eq 'match_part'){
			next;
		}
		my $id = $feature->seq_id();
		if ( !$id ) { next; }
		my $db= $feature->source_tag();
		$db=~s/^t?BLAST[xpn]_//i;
        $db='BLASTDB:'.uc($db);
        #skip some dbs as they are processed differently
        if ( !$noskip && ($db =~ /_EC|_KEGG|_GO/ || $gffblastfile=~ /_EC|_KEGG|_GO/)) {
            $skipped{$db} = 1 if !$skipped{$db};
            next; # Could be 'last' but maybe we've got a concatanated GFF with multiple databases.
        }
        my @notes = $feature->get_tag_values('Note')
          if $feature->has_tag('Note');
        my @hits = $feature->get_tag_values('Alias')
          if $feature->has_tag('Alias');
        my $hit=join(' ',@hits);
        my $descr=join(' ',@notes);
        &process_blast_hit( $id, $hit, $db, $descr, \%hash_id,\%hash_hits ) if $hit;
	}
	return ( \%hash_id, \%hash_hits, \%skipped );
	
}

sub build_blash_hashes_text($){
	my $blastfile=shift;
	print "Processing BLAST textfile $blastfile\n";
	my ( %hash_id, %hash_hits, %skipped,$counter );
	my $total_lines = `wc -l <$blastfile`;chomp($total_lines);
    print "Building HASH\n";
    open( BLAST, $blastfile ) || die;
    my $timer = new Time::Progress;
    $timer->attr( min => 0, max => $total_lines );
    $timer->restart;

    #used during parsing
    my ( $id,  $db );
    my ( $hit, $descr );
    while ( my $line = <BLAST> ) {
        $counter++;
        if ( $line =~ /^\s*$/ || $line =~ /^--/ ) { next; }
        if ( $line =~ /Length =|letters\)/ ) { next; }

        if ( $counter =~ /0000$/ ) {
            print $timer->report( "eta: %E min, %40b %p\r", $counter );
        }
        chomp($line);
        if ( $line =~ /^Query= (\S+)/ ) {

            #process last hit of previous query
            if ( $id ) {
                if (!$hit){$hit='None';}
                &process_blast_hit( $id, $hit, $db, $descr, \%hash_id,
                                    \%hash_hits );
            }
            $id = $1;
            undef($hit);
            undef($descr);
            undef($db);
            unless ($id=~/^[A-Z]{2}\d+[A-Z][a-z]\w+$/){
                # Not in ic_create_naming format
                undef($id);
                next;
            }
        } elsif ( $line =~ /^Database: (.+)$/ && $id) {
            $db = $1;
            chomp($db);
            $line=<BLAST>;
            if ($line=~/^\w+/){
                chomp($line);
                $db.=' '.$line;
            }
            $db =~s/\.fa?st?a.+//;
            $db=~s/\s+$//;
            $db=~s/^\s+//;
            $db =~s/^.+\///;
            $db=~s/^t?BLAST[xpn]_//i;
            $db='BLASTDB:'.uc($db);
            undef($hit);
            undef($descr);
        }
    unless ($nohits){
#get description. Subject's second line is always with a space but Query's second line does not.
        if ( $line =~ /^>|\s+/ && $id && $db) {

            #skip some dbs as they are processed differently
            if ( !$noskip && ($db =~ /_EC|_KEGG|_GO/ || $blastfile=~ /_EC|_KEGG|_GO/)) {
                $skipped{$db} = 1;
                next;
            }
            if ( $line =~ /^>(\S+)\s*(.*)\s*$/ ) {

                #new hit, process previous one
                if ($hit) {
                    &process_blast_hit( $id, $hit, $db, $descr, \%hash_id,
                                    \%hash_hits );
                }
                $hit   = $1;
                $descr = $2;
                $hit=~s/[\w\-_\.]+\|//g;
            }
        }
            #second descr line
        elsif ( $line =~ /^\s+(\S.+)/ && $hit && $db && $id) {
                $descr .= " " . $1;
        }
    }
    }
    close(BLAST);
	return ( \%hash_id, \%hash_hits, \%skipped );
}

sub process_blast_hit ($$$$$$) {
	my $id  = shift || next;
	my $hit = shift || next;
	my $db  = shift ||die;
	my $descr         = shift;
	my $hash_id_ref   = shift || die;
	my $hash_hits_ref = shift ||die;

	# if we have no hits, still store DB.

	unless ($hit eq "None"){
		
		
#The problem with Tax is that it is not delimited except by end of line or another tag
# we assume, thus, that Tags appear at the end of any other description.
# safest is to grab tax to end of line and then trim back.
	if ( $descr && $descr =~ /Tax=(.+)/ ) {
		my $tax = $1;
		$tax =~ s/Tax=([\w\s]+)/$1/;
		$tax =~ s/\w+=.+$//g;
		$tax =~ s/\s+\(.+$//g;
		$tax =~ s/[\'\"]//g;
		$tax =~ s/\s+s?u?b?sp\..*$//g;
		$tax =~ s/\s+str\..*$//g;
		$tax =~ s/\s+serovar\..*$//g;
		$tax =~ s/^\s+|\s+$//;
		$hash_id_ref->{$id}{'taxa'}{$tax}++ if ($tax);
		#$descr =~ s/Tax=.+//;
	} elsif ( $descr && $descr =~ /OS=(.+)/ ) {
		my $tax = $1;
		$tax =~ s/OS=([\w\s]+)/$1/g;
		$tax =~ s/\w+=.+$//g;
		$tax =~ s/\s+\(.+$//g;
		$tax =~ s/[\'\"]//g;
		$tax =~ s/\s+s?u?b?sp\..*$//g;
		$tax =~ s/\s+str\..*$//g;
		$tax =~ s/\s+serovar\..*$//g;
		$tax =~ s/^\s+|\s+$//;
		$hash_id_ref->{$id}{'taxa'}{$tax}++ if ($tax);
		#$descr =~ s/OS=.+//;
	}
	if ($descr) {
		$hash_hits_ref->{$hit}{'descr'} = $descr;
	}
	$hash_hits_ref->{$hit}{'db'} = $db;
	$hash_id_ref->{$id}{"hits"}{$hit}++;
	}
	$hash_id_ref->{$id}{"dbs"}{$db}++;
}
sub process_blast ($) {
#TODO: even if there is no hit, still add the blastdb in the dbxref so we know at least that we searched for it.
# We have a custom BLAST parser because we will process huge amounts of data and we don't want to fall asleep.
# the alternative would have been to process the .gff file we produce for Bio::SeqFeature::Store
	my $blastfile = shift || die();
	unless ( -s $blastfile ) { die("Cannot find $blastfile"); }

	# the following should go to a general section but it is fine for now as we expect only one input file.
	my $select_cvterm_sql ="select cvterm_id from cvterm where name=? and is_obsolete=0 limit 1";
	my $select_cvterm_prepare = $dbh->prepare($select_cvterm_sql);
	$select_cvterm_prepare ->execute('inferred from electronic annotation');
	my ($code_IEA) = $select_cvterm_prepare ->fetchrow_array;
	if ( !$code_IEA ) {
		die("Cannot get cvterm for inferred from electronic annotation part of the Evidence Codes CV. Have you loaded it?\n");
	}
	# we can add another CVTERM here for inserting the fact a database was used. a decent CV would be 'local' and thus db='null'
	my $blastdb_term='searched_against';
	my $blastdb_term_desc='Defines that feature has been processed with a similarity search using a specific database. Used in featureprop where the value defines which database was used.';
	my $blastdb_term_desc_escaped=$dbh->quote($blastdb_term_desc);
	# the following sqls to insert cvterms are better than other found in this script-> 
	#TODO update the others.

	my $select_test_dbxref_sql="SELECT dbxref_id from dbxref where db_id=(select db_id from db where name=?) and accession=? and description=?";
	my $select_test_cvterm_sql="SELECT cvterm_id from cvterm where cv_id=(select cv_id from cv where name=?) AND name=? AND definition=? AND dbxref_id=(select dbxref_id from dbxref where db_id=(select db_id from db where name='null') AND accession=?)";
	my $insert_dbxref_sql="INSERT into dbxref (db_id,accession,description) VALUES ((select db_id from db where name=?),?,?)";
	my $insert_cvterm_sql="INSERT into cvterm (cv_id,name,definition,dbxref_id) VALUES ((select cv_id from cv where name=?),?,?,(select dbxref_id from dbxref where db_id=(select db_id from db where name='null') AND accession=?))";

	my $select_test_dbxref_prepare=$dbh->prepare($select_test_dbxref_sql);
	my $select_test_cvterm_prepare=$dbh->prepare($select_test_cvterm_sql);
	my $insert_dbxref_prepare=$dbh->prepare($insert_dbxref_sql);
	my $insert_cvterm_prepare=$dbh->prepare($insert_cvterm_sql);


	$select_test_dbxref_prepare->execute('null',$blastdb_term,$blastdb_term_desc_escaped);
	$select_test_cvterm_prepare->execute('local',$blastdb_term,$blastdb_term_desc_escaped,$blastdb_term);
	my ($dbtest)= $select_test_dbxref_prepare ->fetchrow_array;
	my ($cvtest)=$select_test_cvterm_prepare->fetchrow_array;
	unless($dbtest){
		$insert_dbxref_prepare->execute('null',$blastdb_term,$blastdb_term_desc_escaped);
	}
	unless ($cvtest){
		$insert_cvterm_prepare->execute('local',$blastdb_term,$blastdb_term_desc_escaped,$blastdb_term);
	}

	$select_cvterm_prepare ->execute($blastdb_term);
	my ($blastdb_termid)= $select_cvterm_prepare ->fetchrow_array;
	if (!$blastdb_termid){die ("Failed to get the cvterm for the BLAST db term, $blastdb_term.\n");}

	# if a featureprop is already added by this script don't check for it.
	my %existence_hash;	#$feature_id,$blastdb_termid,$db
	
	
	print "Preprocessing...\n";
	if ($nohits){
		unless ($blastfile=~/_dbonly/||$gff_format){
			system ("egrep '^Query=|^Database: ' $blastfile >$blastfile.tmp");
			$blastfile.='.tmp';
		}
	}
    my ( $hash_id_ref, $hash_hits_ref, $skipped_ref );
    if (!$gff_format){
        ( $hash_id_ref, $hash_hits_ref, $skipped_ref )=&build_blash_hashes_text($blastfile) ;
    }else{
        ( $hash_id_ref, $hash_hits_ref, $skipped_ref )=&build_blash_hashes_gff($blastfile) ;
    }
    my %hash_id=%$hash_id_ref;
    my %hash_hits=%$hash_hits_ref;
    my %skipped=%$skipped_ref;

	my $skip;
	foreach my $db ( sort keys %skipped ) { $skip .= " " . $db; }
	if ($skip) { warn "These databases were skipped:$skip\n"; }
	$skip    = '';
	%skipped = ();
	print "Parsing HASH\n";
	my $get_feat_id_sql      = 'SELECT feature_id from feature where uniquename=?';
	my $select_featureprop_rank_sql='SELECT max(rank) from featureprop where feature_id=? AND type_id=?';
	my $select_featureprop_id_sql='SELECT featureprop_id from featureprop where feature_id=? AND type_id=? AND value=?';
	my $check_featureprop_sql='SELECT featureprop_id from featureprop where feature_id=? AND type_id=? AND rank=?';
	my $insert_featureprop_sql = 'INSERT into featureprop (type_id,feature_id,value,rank) VALUES (?,?,?,?)';
	my $insert_into_db_sql = 'INSERT into db (name) VALUES (?)';
	my $select_db_sql='SELECT db_id from db where name=?';
	my $select_dbxref_sql='SELECT dbxref_id from dbxref where db_id=? AND accession=?';
	my $insert_into_dbxref_sql='INSERT into dbxref (db_id,accession,description) VALUES (?,?,?)';
	my $select_feature_dbxref_sql='SELECT dbxref_id from feature_dbxref where feature_id=? AND dbxref_id=?';
	my $insert_into_feature_dbxref_sql='INSERT into feature_dbxref (feature_id,dbxref_id) VALUES (?,?)';

	my $prepared_get_feat_id = $dbh->prepare($get_feat_id_sql);
	my $prepared_select_featureprop_rank=$dbh->prepare($select_featureprop_rank_sql);
	my $prepared_select_featureprop_id=$dbh->prepare($select_featureprop_id_sql);
	my $prepared_check_featureprop=$dbh->prepare($check_featureprop_sql);
	my $prepared_insert_featureprop =$dbh-> prepare($insert_featureprop_sql);
	my $prepared_insert_into_db=$dbh-> prepare($insert_into_db_sql);
	my $prepared_select_db_check =$dbh-> prepare($select_db_sql);
	my $prepared_select_db_cached =$dbh-> prepare_cached($select_db_sql);
	my $prepared_select_dbxref_check =$dbh-> prepare ($select_dbxref_sql);
	my $prepared_select_dbxref_cached =$dbh-> prepare_cached ($select_dbxref_sql);
	my $prepared_insert_into_dbxref =$dbh-> prepare ($insert_into_dbxref_sql);
	my $prepared_select_feature_dbxref=$dbh-> prepare ($select_feature_dbxref_sql);
	my $prepared_insert_into_feature_dbxref=$dbh-> prepare ($insert_into_feature_dbxref_sql);
	

	# prepare db and dbxref
	my %dbs_hit;
	my %taxa_found;
	# build db array -> would be more efficient to have built it previously during file parsing
#	foreach my $hit (keys %hash_hits){
#		$dbs_hit{$hash_hits{$hit}{"db"}}{'exist'}=1;
#	}
	foreach my $id (keys %hash_id){
		foreach my $tax (keys %{$hash_id{$id}{"taxa"}}){
			$taxa_found{$tax}{'exist'}=1;
		}
		foreach my $db (keys %{$hash_id{$id}{"dbs"}}){
			$dbs_hit{$db}{'exist'}=1;
		}
	}

	#add dbs, store db_id
	$dbh->begin_work;

	foreach my $db (keys %dbs_hit) {
		$db=~s/^\s*|\s*$//g;
		$prepared_select_db_check->execute($db);
		my ($db_id) = $prepared_select_db_check->fetchrow_array;
		if (!$db_id){
			$prepared_insert_into_db->execute($db);
			$prepared_select_db_check->execute($db);
			($db_id) = $prepared_select_db_check->fetchrow_array;
			if (!$db_id){die ("The insert command failed for $db.\n");}
		}
		#cache it if it does not exist
		$dbs_hit{$db}{'db_id'}=$db_id;
	}
	$dbh->commit;
	# we shall add taxonomy as a dbxref of the NCBI taxonomy. Use a hash built previously.
	if (%taxa_found && keys %taxa_found>0){
	print "Preparing taxonomy reports. This may take a while...\n";
	$dbh->begin_work;
	foreach my $tax (sort keys %taxa_found){
		#add db
		my $db='NCBI_TAXONOMY';
		$prepared_select_db_check->execute($db);
		my ($db_id) = $prepared_select_db_check->fetchrow_array;
		if (!$db_id){
			$prepared_insert_into_db->execute($db);
			$prepared_select_db_check->execute($db);
			($db_id) = $prepared_select_db_check->fetchrow_array;
			if (!$db_id){die ("The insert command failed for $db\n");}
		}
		#cache it if it does not exist.
		$dbs_hit{$db}{'db_id'}=$db_id if !$dbs_hit{$db}{'db_id'};
		# add dbxref and cache it in hash.
		my $ncbi_taxid=get_ncbi_taxid($tax);
		if (!$ncbi_taxid){
			#warn "No NCBI TAXID for $tax. Skipping.\n";
			next;
		}
		$prepared_select_dbxref_check->execute($db_id,$ncbi_taxid);
		my ($tax_id) = $prepared_select_dbxref_check->fetchrow_array;
		if (!$tax_id){
			$prepared_insert_into_dbxref->execute($db_id,$ncbi_taxid,$tax);
			$prepared_select_dbxref_check->execute($db_id,$ncbi_taxid);
			($tax_id) = $prepared_select_dbxref_check->fetchrow_array;
			if (!$tax_id){die ("The insert command failed for taxon $tax ($ncbi_taxid) of db $db_id\n");}
		}
		#cache it if it does exist
		$taxa_found{$tax}{'dbxref_id'}=$tax_id;
	 }
	$dbh->commit;
	}
	print "Done. Parsing hits into dbxref...\n";
	$dbh->begin_work;
	#add HIT as a dbxref, use db_id from where it was found -> redundancy if two blastdbs are redundant
	foreach my $hit (keys %hash_hits){
		my $descr=$hash_hits{$hit}{"descr"};
		my $descr_escaped=$dbh->quote($descr);
		my $db=$hash_hits{$hit}{"db"};
		if (!$db){die "Cannot find db for $hit\n";}
		#$prepared_select_db_check->execute($db);
		$prepared_select_db_cached->execute($db);
		my ($db_id) = $dbs_hit{$db}{'db_id'};
		if (!$db_id){$db_id = $prepared_select_db_cached->fetchrow_array;}
		#my ($db_id) = $dbs_hit{$db}{'db_id'};
		if (!$db_id){die ("The insert command had failed for $db.\n");}
		# rely on unique index key? NO
		$prepared_select_dbxref_check->execute($db_id,$hit);
		my ($dbxref_id)=$prepared_select_dbxref_check->fetchrow_array;
		if (!$dbxref_id){$prepared_insert_into_dbxref->execute($db_id,$hit,$descr_escaped);}
		$prepared_select_dbxref_cached->execute($db_id,$hit);
		($dbxref_id)=$prepared_select_dbxref_cached->fetchrow_array;
		if (!$dbxref_id){die;}
	}
	$dbh->commit;
	my $id_number=scalar(keys %hash_id);
	print "Done. Storing hits from $id_number features into featureprop (and taxonomy hits into feature_dbxref)...\n";
	my $timer = new Time::Progress;
	$timer->attr( min => 0, max => scalar( keys %hash_id ) );
	$timer->restart;
	my $counter = int(0);
	foreach my $id ( keys %hash_id ) {
	#debugprint "processing $id\n";
		
		$counter++;
		if ( $counter =~ /000$/ ) {print $timer->report( "eta: %E min, %40b %p\r", $counter );}
		#get feature_id from db. if it doesn't exist skip and warn.
		$prepared_get_feat_id->execute($id);
		my ($feature_id) = $prepared_get_feat_id->fetchrow_array;
		if ( !$feature_id ) { $skipped{$id} = 1; next; }
		my @dbs_id  = sort ( keys %{ $hash_id{$id}{"dbs"} } );
		my @hits = sort ( keys %{ $hash_id{$id}{"hits"} } );
		my @taxa = sort ( keys %{ $hash_id{$id}{"taxa"} } );

		$dbh->begin_work;

		# add dbs into featureprop
		my $rank_dbs;
		foreach my $db (sort @dbs_id){
			#debug print $db;
			#check if term exists, if it does ->next
			my $exists;
			if (exists $existence_hash{$feature_id}{$blastdb_termid}{$db}){$exists=1;}
			else {
				$prepared_select_featureprop_id->execute($feature_id,$blastdb_termid,$db);
				($exists)=$prepared_select_featureprop_id->fetchrow_array;
			}
			if (!defined $exists){
#warn "Did not find featureprop entry for feature_id $feature_id term_id $blastdb_termid value $db \n";
				#get next available rank
				$prepared_select_featureprop_rank->execute($feature_id,$blastdb_termid);
				if (!defined $rank_dbs){
					($rank_dbs)=$prepared_select_featureprop_rank->fetchrow_array;
					if (!defined $rank_dbs){
						$rank_dbs=int(0);
					}
					else {$rank_dbs++;}
				}
#warn "Trying to insert into featureprop feature_id $feature_id term_id $blastdb_termid value $db rank $rank_dbs\n";
				my $success=$prepared_insert_featureprop->execute($blastdb_termid,$feature_id,$db,$rank_dbs);
				$rank_dbs++;
#warn "Result is $success and rank incremented to $rank_dbs\n";
			}
			$existence_hash{$feature_id}{$blastdb_termid}{$db}=1;
#die;
		}
		unless ($nohits){
		my $rank_term;
		
		foreach my $hit (@hits) {
			#get db_id & dbxref_id from cache; die otherwise
			my $db=$hash_hits{$hit}{"db"};
			my ($db_id) =$dbs_hit{$db}{'db_id'};
			$prepared_select_dbxref_cached->execute($db_id,$hit);
			my ($dbxref_id)=$prepared_select_dbxref_cached->fetchrow_array;
			unless ($db_id && $dbxref_id){die "Fatal error: this term ($hit) is not in any of the defined dbs ($db).\n";}
		
			#prepare descr
			my $descr=$hash_hits{$hit}{"descr"};

			if ($descr){
#debugprint "from $descr\n";
				#remove any tags \w+=\S+ but also 
				$descr =~ s/\w+=.+//g;

				#remove all non-words
				$descr=~s/[\W]+/ /g;
				
				#for description line we will remove common words |\b\b See short word removal above:
				# it's better to look for the few words that are short and should be preserved.
				#$descr =~ s/\bin\b|\band\b|\bthe\b|\ba\b|\bfor\b|\bat\b|\bof\b//g;
				#remove short words (0-3) written in lowercase capital could be keyword. Keep pol?
				$descr=~s/pol\s/Pol /g;
				$descr=~s/gag\s/Gag /g;
				$descr=~s/env\s/Env /g;
				$descr =~ s/\b[a-z]{0,3}\b//g;

				# make lower case for efficient searching
				$descr=lc($descr);
				# remove longer undesired words
				my @undesireds=qw(predicted similar protein putative uncharacterized fragment partial scaffold shotgun);
				foreach my $undesired (@undesireds){
					$descr=~s/\b$undesired\b//g;
				}
				#remove numbers
#				$descr=~s/\b\d+\b/ /g;
				#whitespace
				$descr =~ s/^\s+|\s+$//;
				$descr =~ s/\s\s+/ /g;
				my @descr_array=split(" ",$descr);
				foreach my $term (@descr_array){
					#check if term exists, if it does ->next
					my $exists;
					if (exists $existence_hash{$feature_id}{$code_IEA}{$term}){$exists=1;}
					else{
						$prepared_select_featureprop_id->execute($feature_id,$code_IEA,$term);
						($exists)=$prepared_select_featureprop_id->fetchrow_array;						
					}
					if (!defined $exists){
						#	rank_sql='SELECT rank from featureprop where feature_id=? AND type_id=?';
						#	id_sql='SELECT featureprop_id from featureprop where feature_id=? AND type_id=? AND value=?';
						#	$insert_featureprop_sql = 'INSERT into featureprop (type_id,feature_id,value,rank) VALUES (?,?,?,?)';
						
						#get next available rank
						$prepared_select_featureprop_rank->execute($feature_id,$code_IEA);
						if (!defined $rank_term){
							($rank_term)=$prepared_select_featureprop_rank->fetchrow_array;
							if (!defined $rank_term){
								$rank_term=int(0);
							}
							else {
								$rank_term++;
							}
						}
						#my $rc=$prepared_check_featureprop->execute($code_IEA,$feature_id,$rank_term);
						#if ($rc ne 'E0E' ){
                        #    $rank_term++;
                        #    $rc=$prepared_check_featureprop->execute($code_IEA,$feature_id,$rank_term);
                        #    if ($rc ne '0E0'){
                        #        warn "Bugtracking: Tried 10 times to increase the rank, but there is already a term for this term/feature/rank ($code_IEA,$feature_id,$rank_term). (for $term). Skipping\n";
                        #        next;
                        #    }
                        #}
    				    $prepared_insert_featureprop->execute($code_IEA,$feature_id,$term,$rank_term);
					    $rank_term++;
					}
					$existence_hash{$feature_id}{$code_IEA}{$term}=1;
				}

			}
		}


		foreach my $tax (@taxa){
			#debug print $tax;
			my $tax_dbxref=$taxa_found{$tax}{'dbxref_id'};
			if (!$tax_dbxref){next;}
			my $exists;
			if (exists $existence_hash{$feature_id}{'dbxref'}{$tax_dbxref}){$exists=1;}
			else{
				$prepared_select_feature_dbxref->execute($feature_id,$tax_dbxref);
				($exists)=$prepared_select_feature_dbxref->fetchrow_array;
			}
			unless (defined $exists){
#debugprint  "\twould insert dbxref $tax_dbxref ($tax) for feat $feature_id tax $tax\n";
				$prepared_insert_into_feature_dbxref->execute($feature_id,$tax_dbxref);
				$existence_hash{$feature_id}{'dbxref'}{$tax_dbxref}=1;
			}
		}
		
	}
	$dbh->commit;
	}
	$select_cvterm_prepare -> finish;
	$select_test_dbxref_prepare-> finish;
	$select_test_cvterm_prepare-> finish;
	$insert_dbxref_prepare-> finish;
	$insert_cvterm_prepare-> finish;
	$prepared_get_feat_id -> finish;
	$prepared_select_featureprop_rank-> finish;
	$prepared_select_featureprop_id-> finish;
	$prepared_insert_featureprop -> finish;
	$prepared_insert_into_db-> finish;
	$prepared_select_db_check -> finish;
	$prepared_select_db_cached -> finish;
	$prepared_select_dbxref_check -> finish;
	$prepared_select_dbxref_cached -> finish;
	$prepared_insert_into_dbxref -> finish;
	$prepared_select_feature_dbxref-> finish;
	$prepared_insert_into_feature_dbxref-> finish;
	
	foreach my $id ( sort keys %skipped ) { $skip .= " " . $id; }
	if ($skip) { warn "These Feature IDs were skipped:$skip\n"; }
	# delete tmp
	if ($nohits){unlink($blastfile);}
	else {"\nAll finished well.\n";}
}

sub insert_cvterms($$$) {

	#TODO fix KEGG ->get dbxref??
	my $features_ref       = shift ||die;
	my $db_name_local = shift || die;
	my $cvterm_ref = shift;
	my %features        = %$features_ref;
	my %cvterm=%$cvterm_ref;

	my ( $pub_id);

	#Make sure Pub id exists so we can cvterms
	if (    $db_name_local eq "GO"
		 || $db_name_local eq "EC"
		 || $db_name_local eq "KEGG_PATHWAY" ){
		my $pub_get          = "SELECT pub_id FROM pub where title=?";
		my $prepared_pub_get = $dbh->prepare($pub_get);
		$prepared_pub_get->execute('inferred from electronic annotation');
		($pub_id) = $prepared_pub_get->fetchrow_array;
		if ( !$pub_id ) {
			my $sql ="select cvterm_id from cvterm where name='inferred from electronic annotation';";
			my $code_IEA_prep = $dbh->prepare($sql);
			$code_IEA_prep->execute();
			my ($code_IEA) = $code_IEA_prep->fetchrow_array;
			if ( !$code_IEA ) {
				&disconnectdb($dbh);
				die("Cannot get cvterm for inferred from electronic annotation part of the evidence_codes CV. Have you loaded it?\n"				);
			}
			my $results = $dbh->do("INSERT INTO pub (title,uniquename,type_id) VALUES ('inferred from electronic annotation','IEA',$code_IEA);"			);
			if ( !$results ) {
				&disconnectdb($dbh);
				die("Could not insert into pub the IEA values\n");
			}
			$prepared_pub_get->execute('inferred from electronic annotation');
			($pub_id) = $prepared_pub_get->fetchrow_array;
			if ( !$pub_id ) {
				&disconnectdb($dbh);
				die("Still cannot get the publication ID for the IEA code.\n");
			}
		}
	}
	
	$dbh->begin_work;

	foreach my $feat_name ( keys %features) {
		if ($toskip{$feat_name}){next;}
		my $feature_id = get_feature_id($feat_name);
		if ( !$feature_id ) {
			warn "No feature ID for $feat_name, skipping...\n";
			$toskip{$feat_name}=1;
			sleep(5) if !%toskip; #if it is the first missing feature
			next;
		}

		#IPR goes to feature_dbxref, the others to feature_cvterm
		if ( $db_name_local eq "InterPro" ){
			my $insert ="INSERT INTO feature_dbxref (feature_id,dbxref_id) VALUES (?,?)";
			my $select ="SELECT feature_dbxref_id from feature_dbxref where feature_id=? and dbxref_id=?";
			my $prepared_insert = $dbh->prepare($insert);
			my $prepared_select = $dbh->prepare($select);
			my @dbxrefs = @{ $features{$feat_name}{"dbxrefs"} };
			foreach my $dbxref (@dbxrefs) {
				if ( !$dbxref ) { next; }
				#unless ( $dbxref =~ /^$db_name_local\:/ || $dbxref =~ /^IPR\:/ ){
				#	warn "Unrecognised dbxref $dbxref\n";
				#	next;
				#}
				if ( !$cvterm{$dbxref} ) {
#					warn "Cannot find cvterm for $dbxref... Skipping\n";
					next;
				}
				my $result = $prepared_select->execute( $feature_id, $cvterm{$dbxref} );
				if ( !$result || $result == 0 ) {
#					warn "Adding new term for $dbxref (".$cvterm{$dbxref}.")\n";
				$prepared_insert->execute( $feature_id, $cvterm{$dbxref});
#				}else{
#					warn "Term $dbxref (".$cvterm{$dbxref}.") already exists\n";
				}
			}
		}
		#KEGG goes to feature_cvterm
		if ( $db_name_local eq "KEGG_PATHWAY"){
			my $insert ="INSERT INTO feature_cvterm (feature_id,cvterm_id,pub_id) VALUES (?,?,$pub_id)";
			my $select ="SELECT feature_cvterm_id from feature_cvterm where feature_id=? and cvterm_id=?";
			my $prepared_insert = $dbh->prepare($insert);
			my $prepared_select = $dbh->prepare($select);
			my @dbxrefs = @{ $features{$feat_name}{"dbxrefs"} };
			
			foreach my $dbxref (@dbxrefs) {
				if ( !$dbxref ) { next; }
				unless ( $dbxref =~ /^$db_name_local\:/ ){next;}
				if ( !$cvterm{$dbxref} ) {
					warn "Cannot find cvterm for $dbxref... Skipping\n";
					sleep(1);
					next;
				}
				my $result = $prepared_select->execute( $feature_id, $cvterm{$dbxref} );
				unless ( $result && $result > 0 ) {
					$prepared_insert->execute( $feature_id, $cvterm{$dbxref} );
				}
			}
		}

		#GO EC goes to feature_cvterm
		if ( $db_name_local eq "GO" || $db_name_local eq "EC" ) {
			my $insert ="INSERT INTO feature_cvterm (feature_id,cvterm_id,pub_id) VALUES (?,?,$pub_id)";
			my $select ="SELECT feature_cvterm_id from feature_cvterm where feature_id=? and cvterm_id=?";
			my $prepared_insert = $dbh->prepare($insert);
			my $prepared_select = $dbh->prepare($select);
			my @aliases = @{ $features{$feat_name}{"aliases"} };
			foreach my $alias (@aliases) {
				if ( !$alias ) { next; }
				if ( !$cvterm{$alias} ) {
					warn "Cannot find cvterm for $alias... Skipping\n";
					next;
				}
				my $result = $prepared_select->execute( $feature_id, $cvterm{$alias} );
				unless ( $result && $result > 0 ) {
					$prepared_insert->execute( $feature_id, $cvterm{$alias} );
				}
			}
		}
	}
	$dbh->commit or die $dbh->errstr;
}

sub get_cvterm ($$$$$) {
	my $accession      = shift || die("No accession provided");
	my $cv_name_local  = shift || die("No cv name provided");
	my $db_name_local  = shift || die("No db name provided");
	my $prepared_terms = shift || die("No database connection passed");
	my $select_db_sql = "SELECT db_id,name from db where name=?";
	my $insert_sql_term = "INSERT INTO dbxref (accession,db_id) VALUES (?,?)";
    my $insert_term_prepared= $dbh->prepare($insert_sql_term);
    my $select_db_prepared= $dbh->prepare($select_db_sql); 

	if ( $accession eq "NA" ) { return; }
    
	#strip db name to get accession
	if ( $db_name_local eq "KEGG_PATHWAY" || $db_name_local eq "InterPro" ) {
		$accession =~ s/^\w+\://;
	}
	my $result;
	if ($db_name_local eq "EC" || $db_name_local eq "KEGG_PATHWAY" ) {
        $result =
          $prepared_terms->execute( $cv_name_local, $accession,
                                    $db_name_local );
    }else{
    	$accession=~s/^$db_name_local\://;
        $result = $prepared_terms->execute( $accession, $db_name_local );
    }
    my ($cvterm_id) = $prepared_terms->fetchrow_array;
	if ( !$failed{$accession} && !$cvterm_id ) {
		#$failed{$accession} ="Non-fatal Warning: Unable to query database for cv,acc,db: $cv_name_local $accession $db_name_local ($result)\n";
		$select_global_db_prepared->execute($db_name_local,$db_name_local);
		my ($db_id) = $select_global_db_prepared->fetchrow_array;
		if (!$db_id){
			$insert_global_db_prepared->execute($db_name_local);
            $select_global_db_prepared->execute($db_name_local,$db_name_local);
			($db_id) = $select_global_db_prepared->fetchrow_array;
			if (!$db_id){
				die "Failed to insert new database $db_name_local\n";
			}
		}
		my $res=$select_global_dbxref_prepared->execute($accession,$db_id);
		unless ($res || $res==0){
		  $insert_global_dbxref_prepared->execute($db_id,$accession);
		}
		if ($db_name_local eq "EC" || $db_name_local eq "KEGG_PATHWAY" ) {
            $result =$prepared_terms->execute( $cv_name_local, $accession,$db_name_local );
        }else{
        	#IPR & GO
            $result = $prepared_terms->execute( $accession, $db_name_local );
            ($cvterm_id) = $prepared_terms->fetchrow_array;
        }if (!$cvterm_id){
        	$result = $prepared_terms->execute( $accession, 'DB:'.$db_name_local );
            ($cvterm_id) = $prepared_terms->fetchrow_array;
            if (!$cvterm_id){
            	if ($insert_term_prepared){
            		$db_name_local=~s/^DB://;
            		my $r=$select_db_prepared->execute('DB:'.$db_name_local);
            		my ($db_id,$db_name_c)=$select_db_prepared->fetchrow_array;
            		if (!$db_id){
            			$r=$select_db_prepared->execute($db_name_local);
            			($db_id,$db_name_c)=$select_db_prepared->fetchrow_array;
            		}
            	    $insert_term_prepared->execute($accession,$db_id) if $db_id;
            	    $result = $prepared_terms->execute( $accession, $db_name_c);
                    ($cvterm_id) = $prepared_terms->fetchrow_array;
            	}
            	if (!$cvterm_id){
                    $failed{$accession} ="No id for $accession of $db_name_local ($cv_name_local)\n";
                    return;
            	}
            }
        }
	}
	if ( !$cvterm_id && !$failed{$accession} ) {
		$failed{$accession} ="No id for $accession of $db_name_local ($cv_name_local)\n";
	} else {
		return ($cvterm_id);
	}
}

sub get_feature_id($) {
	my $feature_name = shift;
	if ( !$feature_name ) { &disconnectdb($dbh); die("No id provided"); }
	my $sql          = "select feature_id from feature where uniquename=?";
	my $prepared_sql = $dbh->prepare($sql);
	$prepared_sql->execute($feature_name);
	my ($feature_id) = $prepared_sql->fetchrow_array;
	if   ( !$feature_id ) { return;}
	else                  { return $feature_id; }
}

sub repair_dbxref($) {
	my $db_name_local = shift;

	# using $csvnames;
	my (%update_dbxref);

# parse CSV or TABSV file e.g. with IPRID, ipr description (e.g. obtained from ebi's ipr's biomart)
# the domain and anything after the second field will be added as a property in dbxrefprop
# the cv of which will be hot-build
# e.g. biomart: InterPro Entry Accession,InterPro Entry Name,InterPro Entry Type
	if ( !$csvnames || !$db_name_local ) {
		&disconnectdb($dbh);
		die("No dbname or file given to update dbxref data\n");
	}

# produce a hash $update_dbxref{$id}{"descr"}=$descr; $update_dbxref{$id}{"domain"}=$domain; etc
# also a special hash for inserting into cvterms and then holding the cvterm_id  $domains{$domain}
	my $delimiter=",";
	my $tabs=`grep -c "\t"  $csvnames`;chomp($tabs);
	my $commas=`grep -c ","  $csvnames`;chomp($commas);
	if ($tabs && $tabs > $commas){
		$delimiter="\t";
	}
	open (CSV,$csvnames) or die ("Cannot open $csvnames: $!\n");
	while ( my $row = <CSV> ) {
		if ($row=~/^\s*#/){next;}
		chomp($row);
		my @data = split($delimiter,$row);
		my $id   = $data[0];
		my $descr = $data[1];
		my $domain = $data[2] if $data[2];
		if ($id) {
			if ( $id =~ /^#/ ) { next; }
			if ($descr) {
				$update_dbxref{$id}{"descr"} = $descr;
			}
			if ($domain) {
				$update_dbxref{$id}{"domain"} = $domain;
			}
		}
	}
	close CSV;

	# get db_id of db_name
	my $dbsql       = "select db_id from db where name=?";
	my $prepared_db = $dbh->prepare($dbsql);
	$prepared_db->execute($db_name_local);
	my ($db_id) = $prepared_db->fetchrow_array();
	$prepared_db->finish;
	unless ($db_id) { die("Failed to find db.db_id of $db_name_local\n"); }

# cvterm is cv.name='relationship' cvterm.name='is_a' reference???? -> lose structure if domain id is in value but keep logic :-(
	my $sql_cvterm_sel =
"select cvterm_id from cvterm where name='is_a' and cv_id IN (SELECT cv_id from cv where name='relationship')";
	my $sql_cvterm_sel_prep = $dbh->prepare($sql_cvterm_sel);
	$sql_cvterm_sel_prep->execute();
	my ($cvterm_is_a) = $sql_cvterm_sel_prep->fetchrow_array();

	# rest goes in single transaction
	$dbh->begin_work;
	my $sql_dbxref =
	  "update dbxref set description=? where accession=? and db_id=$db_id";
	my $sql_select =
	  "select dbxref_id from dbxref where accession=? and db_id=$db_id";
	my $sql_prop =
	  "INSERT into dbxrefprop (dbxref_id,type_id,value) VALUES (?,?,?)";
	my $sql_prop_select =
	  "SELECT value from dbxrefprop where dbxref_id=?";
	my $prepared_seq_select  = $dbh->prepare($sql_select);
	my $prepared_sql_dbxref  = $dbh->prepare($sql_dbxref);
	my $prepared_sql_prop    = $dbh->prepare($sql_prop);
	my $prepared_prop_select = $dbh->prepare($sql_prop_select);

# foreach IPR domain get domain type and descr, add it to dbxrefprop. Also update descr. at dbxref.
	foreach my $id ( keys %update_dbxref ) {
		my $descr = $update_dbxref{$id}{"descr"}
		  || warn("No description for $id\n");
		if ( !$descr ) { $descr = "No description"; }
		my $result = $prepared_seq_select->execute($id);
		if ( $result && $result > 0 ) {
			$prepared_sql_dbxref->execute( $descr, $id );
		} else {
			my $sql_dbxref_insert = "INSERT INTO dbxref (db_id,accession,description) VALUES ($db_id,?,?)";
			my $sql_dbxref_insert_prep = $dbh->prepare($sql_dbxref_insert);
			my $insert_result = $sql_dbxref_insert_prep->execute( $id, $descr );
			if ( !$insert_result ) {
				die "Insert to dbxref failed for $db_id $id $descr";
			}
			$prepared_seq_select->execute($id);
		}
		my ($dbxref_id) = $prepared_seq_select->fetchrow_array();
		my $dbprop_result = $prepared_prop_select->execute($dbxref_id);
		if ( $dbxref_id && ( !$dbprop_result || $dbprop_result == 0 ) ) {
			my $prop_result;
			if ( $update_dbxref{$id}{'domain'} ) {
				$prop_result = $prepared_sql_prop->execute( $dbxref_id, $cvterm_is_a, $update_dbxref{$id}{'domain'} );
			}
			unless ( $prop_result && $prop_result > 0 ) {
				warn "Failed to insert dbxrefprop ($dbxref_id,$cvterm_is_a,".$update_dbxref{$id}{'domain'}.")\n";
			}
		}
	}
	$dbh->commit;
	$prepared_seq_select->finish;
	$prepared_sql_dbxref->finish;
	$prepared_sql_prop->finish;
	print "$db_name_local updated!\n";
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
	my $dbConnectionSettings_text=$dbConnectionSettings;
	$dbConnectionSettings_text=~s/password.+/password=XXXX/;
	print "Connected to $dbConnectionSettings_text\n";
	return ($dbh);
}

sub disconnectdb ($) {
	my $dbh = shift;

	#Clean up
	$dbh->disconnect();
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
            	#my $taxon = $taxondb->get_taxon(-taxonid => $ncbi_taxid)  || die "Cannot find your query\n";;
            }
        }
	}
}
