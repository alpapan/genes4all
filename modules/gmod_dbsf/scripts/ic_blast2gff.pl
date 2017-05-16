#!/usr/bin/perl -w
#$Id$

#Considering: put hit descr into name so it is searchable? Probably not!
#Considering: pub dbxref in source? Trial.

=pod 

=head1  NAME 

 ic_blast2gff.pl - Program to convert a BLAST report into a GFF file. 

=head1 VERSION

 Version 1.0 Feb 2009

=head1 SYNOPSIS

 ic_blast2gff.pl [optional options e.g. -cs -ce -o -b -d -l etc] etfile(s)

	-cs|--cutoff_score     Specify bit-score cut-off for BLAST hits.
	-ce|--cutoff_eval      Specify e-value cut-off for BLAST hits.
	-o|organism        Organism value (abbrv or NCBI Taxonomy ID is fine). See perldoc
	-abbreviate     If -abbreviate is given then and the -organism is an integer then it is assumed it is an NCBI taxonomy ID and is translated to an abbreviation
	-check_org	Alternatively if your IDs have NCBI TaxIDs in the form ~/^[A-Z]{2}(\d+)/ then you can use this option to look it up.
	-biofeature        Prepare Bio::Seq::Feature compatible database (def. off -> chado)
	-l|--limit     Limit BLAST hits to the top -l only (sorted by score). Defaults to 10. Set 0 to disable
    
 See perldoc for all options. E.g. ic_blast2gff can be told to split the BLAST report, containing multiple queries, 
 into individual files. HTML and text based BLAST report can also be returned. In addition, if 
 a chado database is provided, ic_blast2gff will ensure that the GFF will be chado compatible.

Use perldoc ic_blast2gff.pl for an overview.

=head1 ORGANISMs

 Since v. 0.995 you can give a tab delimited file with multiple organisms. The first and last column are the only ones looked at.
 The first column is the NCBI Taxonomy ID and the last column is the abbreviation (e.g. H.sapiens). 
 The NCBI ID is matched with the contig ID and the abbreviation is inserted as a value for the organism tag

 It turns out, of course, that abbreviations are not unique; two species can share one. Therefore from est2assembly v.1.05 we
 forgoe abbreviation and allow the NCBI TaxId to be saved directly as the organism attribute. The switch -abbreviate forces the
 old behaviour. 
 
 The GMOD bulk loader, however, does not expect this, which is why you have to use the adaptor AP.pm during loading (--custom_ada AP). 
 Note that you will have to copy it where your GMOD Adaptors reside
  locate Wormbase.pm|grep GMOD 
 e.g. /usr/local/share/perl/5.10.0/Bio/GMOD/DB/Adapter/AP.pm
   
 This adaptor will use a chado VIEW called organism_ncbi_taxa to link the NCBI taxonomy ID with the organism ID. The view is created
 with 
 CREATE VIEW organism_ncbi_taxa AS SELECT organism_dbxref.organism_id, dbxref.description, dbxref.accession AS ncbi_taxid FROM organism_dbxref JOIN dbxref ON organism_dbxref.dbxref_id = dbxref.dbxref_id JOIN db ON db.db_id=dbxref.db_id WHERE db.name='NCBI_TAXONOMY' ORDER BY dbxref.description;
 
 NB Don't forget to grant ALL permissions to the view to the administrator/users

=head1 OPTIONS
 
	* ARGUMENT(S) 			SearchIO compatible file(s) to convert to GFF.

	  -cs|--cutoff_score 		Specify bit-score cut-off for BLAST hits.
	  -ce|--cutoff_eval 		Specify e-value cut-off for BLAST hits.
	  -l|--limit			Limit BLAST hits to the top -l only (sorted by score). Defaults to 10. Set to 0 to disable
	  -biofeature	Prepare Bio::Seq::Feature compatible database (def. off -> chado)
	  -overwrite	Produce a GFF file even if it already exists
	& -d|--dsn|--database_connection Path to dsn file. Only for Chado
				If a proper database connection is provided, the program can query
				the database to aim for unique identifiers. The dsn file should contain only a single line (+ any # comments) 
				E.g. dbi:Pg:dbname=chado;host=localhost;port=5432;user=postgres;password=123321
				Warning: this is using the IPv4 (not local) connection to your database.
	  -o|organism	Organism value (abbrv is fine).
        -abbreviate     If -abbreviate is given then and the -organism is an integer then it is assumed it is an NCBI taxonomy ID and is translated
        -check_org      Alternatively if your IDs have NCBI TaxIDs in the form ~/^[A-Z]{2}(\d+)/ then you can use this option to look it up.
	  -s|--split_blast        	Split input file into individual BLAST reports. If this option is specified, 
				individual BLAST reports will be extracted and written to file. 
				The reports will be written to a subdirectory called 'out'. It's location 
				is the same as the GFF report. Valid parameter are:
					blast      -  BLAST reports only
					html       -  HTML formatted BLAST reports
					blasthtml  -  normal and HTML formatted BLAST reports
	  -f|--format			SearchIO format of infile (defaults to BLAST).
	  -g|--gff_version        	GFF version output (default:3) Used to genetrate GFF reports in a given standard. 
				Currently, only GFF Version 2 and 3 are supported. 
	  -v|--verbose            	specify the verbosity (default:0)

=head1 DESCRIPTION

 Program to convert a BLAST report into a GFF file. The program can be
 told to split the BLAST report, containing multiple queries, into 
 individual files. HTML and text based BLAST report can be returned. If 
 a chado database is provided, ic_blast2gff will ensure that the GFF will be chado compatible.

=head1 TIPS

 a dsn file is simple text file which has one line like this:
 dbi:Pg:dbname=InsectaCentral;host=localhost;port=15432;user=postgres;password=123321

=head1 AUTHORS

	Remo Stierli 1 and Alexie Papanicolaou 2 3 *

	1 Dept. of Computer Science, University of Rhode Island, USA
	2 Dept. of Entomology, Max Planck Institute for Chemical Ecology, Germany
	3 Centre for Ecology and Conservation, University of Exeter, UK
	* recipient of complaints:	alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

 This software is released under the GNU General Public License version 3 (GPLv3).
 It is provided "as is" without warranty of any kind.
 You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html. 
 Please note that incorporating the whole software or parts of its code in proprietary software 
 is prohibited under the current license.

=head1 BUGS & LIMITATIONS

 BioPerl's tile_hsps is a very useful but buggy method. They can throw an exception of the type
  e.g Bio::Root::Exception MSG: Undefined sub-sequence (1,2). Valid range = 1 - 54

 We don't think that these are important to the GFF creating process.

 The timer function which is provided counts the number of contigs which have been processed from the total and estimates
 the remaining time. Since the more complete and complicated contigs are in the beginning, the estimated time is very likely to be off. 
 We get estimates of 200min but the program completes in ca 10-15min.

 No bugs are known. Please report them.

=cut

use strict;
use Getopt::Long;		#Extended processing of command line options
use Pod::Usage;			#Usage message from embedded pod documentation
use Bio::Tools::GFF;
use Bio::Taxon;
use Bio::DB::Taxonomy;
use Bio::SearchIO;
use Bio::SearchIO::Writer::HTMLResultWriter;
use Bio::SearchIO::Writer::TextResultWriter;
use DBI;
use File::Spec;
use Time::Progress;
$| = 1;

#Command line options
my $gffVersion = 3;
my ($help,$biofeature);
my $split="";
my $verbose = 0;
my $limit = 10;
my $blastformat="blast";
my $cutoff_score=1;
my $cutoff_eval=100;
my ($organism,$abbreviate,$overwrite,$dbh,$dsn,$outdir,$check_org,$count_seq);
my $flatfile_dir = '/db/ncbi_taxonomy/';
GetOptions(
	'count_seq:i'	=> \$count_seq,
	'biofeature'		=> \$biofeature,
	'dsn|database_connection:s'	=> \$dsn,
	'o|organism:s'			=> \$organism
	,'split_blast:s'			=> \$split
	,'gff_version:i'		=> \$gffVersion
	,'verbose:i'			=> \$verbose
	,'limit:i'			=> \$limit
	,'cs|cutoff_score:i'		=> \$cutoff_score
	,'ce|cutoff_eval:s'		=> \$cutoff_eval
	,'format:s'			=> \$blastformat
	,'overwrite'		=> \$overwrite,
	,'check_org'		=> \$check_org,
	'outdir:s'         =>\$outdir,
	'abbreviate'       => \$abbreviate,
        'f|taxonomy_flatfile_dir' => \$flatfile_dir,
);
my @infiles=@ARGV;
pod2usage if (!@infiles) ;
if ($limit && $limit == 0 ){undef($limit);}
my %organisms;
if (!$organism && !$check_org){
	warn "No -organism variable given. Are you sure you will not database your results?\n Provide one or give -check_org if your IDs contain NCBI tax IDs\n";
	sleep(3);
}elsif ($organism && -s $organism){	
	# then organism is retrieved from a tab delimited list. First and last element is taxid and abbreviation respectively
	open (ORGANISM,$organism);
	while (my $line=<ORGANISM>){
		if ($line=~/^#|^\s*$/){next;}
		chomp($line);
		my @data=split('\t',$line);
		if ($data[0] && $data[-1]){
			$organisms{$data[0]}=$data[-1];
		}
	}
}elsif ($organism && $organism=~/^\d+$/ && $abbreviate){
    $organism=&get_species_abbr($organism);
}elsif ($organism && $organism eq 'fromdata'){
	undef($organism);
	$check_org =1 ;
}
$verbose = 0 if $verbose < 0;
$verbose = 2 if $verbose > 2;

die "GFF version ".$gffVersion." is not supported\n"
    if ($gffVersion<2 or $gffVersion>3);
    
die "Bio::SearchIO method '".$split."' is not supported\n"
    if (lc($split) !~ /^|blast|html|blasthtml$/i);    

die "Bio::SearchIO method '".$blastformat."' is not supported\n"
    if (lc($blastformat) ne "blast" && lc($blastformat) ne "blastxml");    


foreach my $infile (@infiles){
	my $gff = `basename $infile`;chomp($gff);
	$gff.=".gff";
	if ($outdir){$gff=$outdir.'/'.$gff;}
	if ($biofeature){$gff.=".biofeature";}
	if ((-z $gff || -s $gff) && !$overwrite){warn "$gff already exists, skipping...\n";next;}
	print "Processing $infile as $gff\n";
	if ($dsn){
		$dbh=&connecttodb($dsn);
		&add_pg_functions($dbh);
	}
	&processBLASTFile($dbh,$infile, $gff, $split);
	if ($dsn){&disconnectdb($dbh);}
}

#################################################################
=head1 METHODS

=head3

sub processBLASTFile

 This function reads a .blast file and prepares the content 
 for GFF3.

=cut

sub processBLASTFile($$) {
	my $dbh = shift;
	my $blastFile = shift;
	my $gff = shift;
	my $split = shift;
	my $outdir = "out";

		
	#Preparing BLAST file 
	my $blastReader = Bio::SearchIO->new(
		-format 	=> $blastformat
		,-file 		=> $blastFile
		,-signif	=> $cutoff_eval
		,-score		=> $cutoff_score
	);
	
	#Preparing output files
	my($volume,$directory,$gffFile) = File::Spec->splitpath($gff);
	if($split =~ /blast|html/i){
		$outdir .= $directory."/" if($directory ne '');
		mkdir($outdir);
	}
	$outdir .= "/";
	my $writerhtml = Bio::SearchIO::Writer::HTMLResultWriter->new() if($split =~ /html/i);
	my $writertext = Bio::SearchIO::Writer::TextResultWriter->new() if($split =~ /blast/i);
	my $gffWriter = Bio::Tools::GFF->new(-file => ">$gff",-gff_version => $gffVersion);

	my $counter=int(0);
	my $total_queries=int(0);
	if ($count_seq){
		$total_queries=$count_seq;
		print "You asked for $total_queries queries to be processed, please wait...\n";
	}else{
		$total_queries=`grep -c "Query=" $blastFile` if $blastformat eq 'blast';
		$total_queries=`grep -c "<Iteration>" $blastFile` if $blastformat eq 'blastxml';
		chomp($total_queries);
		print "There are $total_queries queries to be processed, please wait...\n";
	}
	my $timer = new Time::Progress;
	$timer->attr( min => 0, max => $total_queries );
	$timer->restart;

	
	#Parsing BLAST file ========================================
	while(my $result = $blastReader->next_result){
		$counter++;
		if ($counter=~/000$/){print  $timer->report( "eta: %E min, %40b %p\r", $counter );}
		my $query = $result->query_name();
		#strip accession ids.
		if ($query=~/^gnl/){$query=~s/^gnl\|\w+\|//;}
		elsif ($query=~/^lc/){$query=~s/^lc\|//;}
		elsif ($query=~/^gi/){$query=~s/^gi\|\w+\|//;}
		
		#Reading BLAST result ----------------------------------
		
		# dbxref should nice and concise. Let's hope the blast report has it.
		# standardise to uc()
		my $dbname = uc($result->database_name());
		
		# sometimes there is a directory attached to the db. Get rid of it using basename
		# sometimes database name is a friendly name. In that case no need to call basename
		unless ($dbname=~/\s/){$dbname=`basename \"$dbname\"`;chomp($dbname);}
		
		#get rid of any suffix(es)
		$dbname =~s/\.\w+//g;
		$dbname =~s/\s+(\S)/_$1/g;
		$dbname =~s/^(\S+)\s*.*$/$1/;
		
		my $blastType = $result->algorithm();
	
		print $query." - ".$dbname."\n" if $verbose >= 1;;
		my $serial; $serial=int(0);
		while(my $hit = $result->next_hit()){
			if ($limit && $hit->rank()>$limit){last;}
			#Values
			my $seqid = $query;
			#expand here with dbname?
			my $source = $blastType."_".$dbname;
			
			my $type = "match";	#Using SOFA term 'match' as fail-safe
			if ($blastType=~/blastp/i || $blastType=~/tblastx/i || $blastType=~/blastx/i){
				$type = "protein_match";
			}elsif ($blastType=~/blastn/i || $blastType=~/tblastn/i){
				$type = "nucleotide_match";
			}

			my ($start,$end) = $hit->range('query');

			#skip too small hits/hsps.
			if (abs($start-$end)<3){next;}
		
			my ($targetStart,$targetEnd) = $hit->range('sbjct');
			my $score = $hit->significance();
			my $querystrand = $hit->strand('query');
			my $targetStrand = $hit->strand('sbjct'); 
			my $phase = $hit->frame();		# SearchIO frame is -2..+2 #GFF3 frame is 0,1,2
			$phase=~s/[-+]//;

			#Attributes	
			#Get unique id
			my $id;
			if ($dbh && $serial==0){
				$id = &calculateUniqueHitId($dbh,$seqid,$dbname);
				if($id =~ /Hit(\d+)$/){
					$serial = $1;
				}
			}else{
				$id = $seqid.":".$dbname.":Hit".$serial;	
			}
			
			my $parent = $id;
			my $accessionID = $hit->accession();
		
			#Reformating: BLAST -> GFF -----------------------------
			my $blast = Bio::SeqFeature::Generic->new(
				-seq_id		=> $seqid
				,-source_tag	=> $source
				,-primary	=> $type
				,-start		=> $start
				,-end		=> $end
				,-score		=> $score
				,-strand	=> $querystrand
				,-frame		=> "."			# frame doesn't work for HITs
				,-tag		=> {
					ID 	=> $id
					,Name	=> $id
					,Dbxref	=> $dbname.":".$accessionID
				});
		my $hit_length=$hit->length();

		#start < end for Target in chado as strand is no longer used
		if ($targetStart>$targetEnd){
			my $temp=$targetEnd;
			$targetEnd=$targetStart;
			$targetStart=$temp;
		}
		#Add 'Target' & co attributes to the feature. breaks in bio:seqfeature
		unless ($biofeature){
    		$blast->add_tag_value('Target', $accessionID);
    		$blast->add_tag_value('Target', $targetStart);
    		$blast->add_tag_value('Target', $targetEnd);
#    		$blast->add_tag_value('Target', $querystrand);
		}
		$blast->add_tag_value("Alias","$accessionID");
    		
    		$blast->add_tag_value('Note', $hit->description()) 		if($hit->description() && $hit->description() ne 'No definition line found');
    		$blast->add_tag_value('locus', $hit->locus()) 			if($hit->locus());
    		$blast->add_tag_value('length', $hit_length) 			if($hit_length);
    		$blast->add_tag_value('algorithm', $hit->algorithm()) 		if($hit->algorithm());
    		$blast->add_tag_value('total significance', $hit->significance()) 	if($hit->significance());
    		$blast->add_tag_value('top raw_score', $hit->raw_score()) 		if($hit->raw_score());
    		$blast->add_tag_value('top bit_score', $hit->bits()) 		if($hit->bits());
    		$blast->add_tag_value('logical_length', $hit->logical_length()) if($hit->logical_length());
    		$blast->add_tag_value('hit_rank', $hit->rank())			if($hit->rank());
	    	$blast->add_tag_value('fraction_of_identical_positions', sprintf("%.4f",($hit->matches("id")/$hit_length)))	if($hit->matches("id"));
	    	$blast->add_tag_value('fraction_of_conserved_positions', sprintf("%.4f",($hit->matches("cons")/$hit_length)))	if($hit->matches("cons"));
    		
			#Write GFF feature -------------------------------------
			print $gffWriter->_gff3_string($blast) if $verbose >= 2;
			$seqid=~/^[A-Z]{2}(\d+)/;
			my $taxid=$1 if $1;
			if ($organism && -f $organism && $taxid){
				$blast->add_tag_value('organism', $organisms{$taxid}) if $organisms{$taxid};
			}elsif ($check_org && $taxid){
				$organisms{$taxid}=&get_species_abbr($taxid) if !$organisms{$taxid};
				$blast->add_tag_value('organism', $organisms{$taxid}) if $organisms{$taxid};
			}elsif ($organism) {
		    		$blast->add_tag_value('organism', $organism);
			}
			$gffWriter->write_feature($blast);
			
			
			#Represent HSPs as 'match_part' ========================
			my $hspSerial; $hspSerial = 0;
			while(my $hsp = $hit->next_hsp()){
				
				#Skip HSP if alignment is too short
				if($hsp->length('sbjct') < 2){next;}
				
				#Values
				my ($hspStart,$hspEnd) = $hsp->range('query');
				my ($hspTargetStart,$hspTargetEnd) = $hsp->range('sbjct');
				
				my $hspScore = $hsp->bits();
				if (!$hspScore){$hspScore=$hsp->score();}				

				my $hspquerystrand = $hsp->strand('query');
				
				#todo : get $hsp frame into the hit?					
				my $hspPhase = $hsp->frame('query');	# SearchIO for HIT frame is -2..+2 SearchIO for HSP is 0,1,2! 
				#$hspPhase=~s/[-+]//;			#GFF3 frame is 0,1,2
				
				#Attributes				
				#Get unique id
				my $hspId;
				if ($dbh && $hspSerial == 0){
					$hspId = &calculateUniqueHspId($dbh,$id,$dbname);
					if($hspId =~ /Hsp([0-9]*)$/){
						$hspSerial = $1;
					}
				}else{
					$hspId = $id.":Hsp".$hspSerial;	
				}
				my $parent = $id;
				my $hspAlias = "Similar to $accessionID";
				
				#Reformating: BLAST-HSP -> GFF -----------------------------
				my $hspBlast = Bio::SeqFeature::Generic->new();
				$hspBlast->gff_format(Bio::Tools::GFF->new(-gff_version => 3));
				$hspBlast->set_attributes(
					-seq_id		=> $seqid
					,-source_tag	=> $source
					,-primary	=> "match_part"
					,-start		=> $hspStart
					,-end		=> $hspEnd
					,-score		=> $hspScore
					,-strand	=> $hspquerystrand
					,-frame		=> $hspPhase
					,-tag		=> {
						ID 		=> $hspId
						,Name	=> $hspId
#no searching for match_part
#						,Alias	=> $hspAlias
#						,Dbxref	=> $dbname.":".$accessionID
						,Parent	=> $id

					});
					
				#start must < the end for chado as strand is no longer used
				if ($hspTargetStart>$hspTargetEnd){
					my $temp=$hspTargetEnd;
					$hspTargetEnd=$hspTargetStart;
					$hspTargetStart=$temp;
				}
				
				#Add 'Target' & co attributes to the feature. breaks in biofeature
#not for match_part?
#				unless ($biofeature){
#	    		$hspBlast->add_tag_value('Target', $hspId);
#	    		$hspBlast->add_tag_value('Target', $hspTargetStart);
#	    		$hspBlast->add_tag_value('Target', $hspTargetEnd);
#	    		#$hspBlast->add_tag_value('Target', $hspquerystrand);
#				}
		        $hspBlast->add_tag_value('raw_score', $hsp->score())        	if($hsp->score());
				$hspBlast->add_tag_value('p_value', $hsp->pvalue())		if($hsp->pvalue());
	    		$hspBlast->add_tag_value('bit_score', $hsp->bits())		if($hsp->bits());
	    		$hspBlast->add_tag_value('significance', $hsp->evalue())	if($hsp->evalue());
	    		$hspBlast->add_tag_value('length', $hsp->length())		if($hsp->length());
	    		$hspBlast->add_tag_value('gaps', $hsp->gaps())			if($hsp->gaps());
	    		$hspBlast->add_tag_value('hsp_rank', $hsp->rank())			if($hsp->rank());
		
				#Write GFF feature -------------------------------------
	                        if ($organism && -f $organism && $taxid){
					$hspBlast->add_tag_value('organism', $organisms{$taxid}) if $organisms{$taxid};
				}elsif ($check_org && $taxid){
					$organisms{$taxid}=&get_species_abbr($taxid) if !$organisms{$taxid};
					$hspBlast->add_tag_value('organism', $organisms{$taxid}) if $organisms{$taxid};
				}
				elsif ($organism) {
					$hspBlast->add_tag_value('organism', $organism);
				}
				print $gffWriter->_gff3_string($hspBlast) if $verbose >= 2;
	    		$gffWriter->write_feature($hspBlast);
		    
				$hspSerial++;
			}
			
			$serial++;
		}
		# This has to go to the end because otherwise, it doesn't process the HITs.
		#Writing BLAST report ------------------------------
		if ($split){
			if($split =~ /html/i){
				my $outhtml = Bio::SearchIO->new(
					-writer => $writerhtml,
					-file => ">".$outdir.$query.".html");
				$outhtml->write_result($result);
				$outhtml->close();
			}
			if($split =~ /blast/i){
			        my $outtext = Bio::SearchIO->new(
			        	-writer => $writertext,
			        	-file => ">".$outdir.$query.".txt");
		  		$outtext->write_result($result);
		  		$outtext->close();
			}
		}
	}
	my $elapsed=$timer->report("%L");
	print "\nTotal time elapsed: $elapsed\n";
}

=head3

sub calculateUniqueSubfeatureId

 This function queries the database for the highest subfeature
 within the match group (e.g.match, protein_match, 
 nucleotide_match). The reslut is used to calculate an 
 incremented unique id.

=cut

sub connecttodb ($){
    my $dsn = shift;
    my $dbConnectionSettings;
    if (-s $dsn){
    # Reading the database connection settings from a seperate file.
    # E.g. dbi:Pg:dbname=InsectaCentral;host=localhost;port=15432;user=postgres;password=123321
    open(DSN, $dsn) || die("Error: Cannot open dsn file");
    my $dbConnectionSettings ;
    while (my $line=<DSN>){
        if ($line=~/^\s*#/){next;}
        elsif ($line=~/^\s*(\S+)\s*$/){$dbConnectionSettings = $1;}
    }
    close(DSN);
    }
    else {
        $dbConnectionSettings = $dsn;
    }
    my $dbh = DBI->connect($dbConnectionSettings);
    die("Error: Unable to connect to database") if(!$dbh);
    return ($dbh);
}

sub disconnectdb ($){
	my $dbh = shift;
	#Clean up
	$dbh->disconnect();
}
sub calculateUniqueSubfeatureId($$$$$) {
	my $dbh = shift;
	my $srcfeature = shift;
	my $matchgroup = shift;
	my $key = shift;
	my $dbname=shift;

	#get last subfeature index available and all the other ones continue from there (so only query db once)
	my $sql = "SELECT hit.uniquename AS id";
	$sql .= " FROM feature AS hit, feature AS src, featureloc AS loc";
	$sql .= " WHERE loc.feature_id = hit.feature_id";
	$sql .= " AND loc.srcfeature_id = src.feature_id";
	$sql .= " AND src.uniquename = '" . $srcfeature . "'";
	$sql .= " AND hit.type_id IN (".$matchgroup.")";
	$sql .= " ORDER BY ic_accession_serial(hit.uniquename) DESC";
	$sql .= " LIMIT 1";
	my $sth = $dbh->prepare($sql);
	
	$sth->execute or die("Error: Unable to query database");
	my $vals = $sth->fetchrow_hashref;
	#die ("val from $sql is $vals");

	my $newid = ":".$dbname.":".$key."0";
	if ($vals->{'id'}){
		if ($vals->{'id'} =~ /Hit([0-9*])$/){
			$newid = $1 + 1;
			$newid = ":".$key.$newid;
		}
	}
	$sth->finish;
	return $srcfeature.$newid;
}

=head3

sub calculateUniqueHitId

 This function queries the database for the highest hit id 
 associated with the source feature (subfeature thereof).
 The reslut is used to calculate an incremented unique id.

=cut

sub calculateUniqueHitId($$$){
	my $dbh = shift;
	my $srcfeature = shift;
	my $dbname=shift;
	#485 match
	#491 nucleotide_match
	#489 protein_match
	return &calculateUniqueSubfeatureId($dbh,$srcfeature,"491,489,485","Hit",$dbname);
}

=head3

sub calculateUniqueHspId

 This function queries the database for the highest hsp id 
 associated with the source feature (subfeature thereof).
 The reslut is used to calculate an incremented unique id.

=cut

sub calculateUniqueHspId($$$){
	my $dbh = shift;
	my $srcfeature = shift;
	my $dbname=shift;
	
	#182 match_part
	return &calculateUniqueSubfeatureId($dbh,$srcfeature,"182","Hsp",$dbname);
}


sub add_pg_functions($){
    my $dbh=shift;
    my @sqls;
    push (@sqls,
    'CREATE OR REPLACE FUNCTION ic_accession_serial(varchar) RETURNS int AS $$ SELECT substring($1 '
    ."from e'\\d+\$')::int;"
    .' $$ LANGUAGE SQL IMMUTABLE STRICT');
    push (@sqls,
    'CREATE OR REPLACE FUNCTION ic_accession_assembly(varchar) RETURNS char(2) AS $$ SELECT substring(substring($1, '
    ."e'^\\w{2}\\d+[A-Z][a-z]')"
    .' from length(substring($1,'
    ." e\'^\\w{2}\\d+[A-Z][a-z]\') )-1"
    .' for 2 ); $$ LANGUAGE SQL IMMUTABLE STRICT');
    foreach my $sql(@sqls){
       my $sth = $dbh ->prepare($sql);
       $sth->execute or die("Error: Unable to store function in database");
    }
}

sub get_species_abbr($) {
        my $ncbi_taxid = shift;
	my ($genus,$species,$abbr);
        my $taxondb = Bio::DB::Taxonomy->new( -source => 'entrez' );
        if (    $flatfile_dir
                 && -d $flatfile_dir
                 && -s $flatfile_dir . '/nodes.dmp'
                 && -s $flatfile_dir . '/names.dmp' )
        {
                $taxondb =
                  Bio::DB::Taxonomy->new(
                                                                  -source    => 'flatfile',
                                                                  -nodesfile => $flatfile_dir . '/nodes.dmp',
                                                                  -namesfile => $flatfile_dir . '/names.dmp',
                                                                  -directory => $flatfile_dir
                  );
        }
 my $taxon = $taxondb->get_taxon(-taxonid => $ncbi_taxid)  || die "Cannot find your species abbreviation\n";;
 my $species_rank=$taxon->rank();
 $species = $taxon->scientific_name() if ($species_rank eq 'species');
 if (!$species){
	my $t=$taxon;
	while ($t=$t->ancestor()){
		$species_rank=$t->rank();
		next if !$species_rank;
		$species = $taxon->scientific_name() if ($species_rank eq 'species'); 
	}
	if (!$species){
	 	die "Cannot find a scientific species name for your query at $species_rank\n";
	}
 }
 $species =~s/^(\w+)\s//;
 my $t=$taxon;
 while ($t=$t->ancestor()){
        my $rank=$t->rank();
        next if !$rank;
        my $name = $t->scientific_name();
        next if !$name;
        if ($rank eq 'genus'){
                $genus=$name;
                last;
        }
 }
 if (!$genus || !$species){die "Cannot find genus or species\n";}
 $genus=~/^(\w)/;
 $abbr=$1.'.'.$species;
 return $abbr;
}

