#!/usr/bin/perl -w
#$Id$

#CHANGELOG
#13JUL09
# added HMMPIR.
# Removed dependency on physprop if using only IPRSCAN
#04Aug09
# Dbxref IPR renamed to InterPro for Chado compat.


=head1 NAME 

 ic_annot8r2gff.pl

=head1 USAGE

		'physprop:s'			=> annot8r_physprop output. Necessary if loading GO EC KEGG
		'dsn|database_connection:s'	=> DSN connection file, needed to define the features unique IDs and to define the local feature types
		'nodb'				=> Disable the need (and benefits of the above). Feature IDs will start from ONE so make sure the db you will upload to will be empty for these datatype/refseq combination
		'id:s'				=> Force prefix for the unique id code; such as GO EC KEGG IPRSCAN. Defaults to autodetect
		'gff_version:i'			=> GFF version (2 or 3; defaults to 3)
		'source:s'			=> Value for Source field
		'cs|cutoff_score:i'		=> Cutoff bitscore for consideration
		'ce|cutoff_eval:s'		=> Cutoff e-value for consideration
        -o|organism        Organism value (abbrv or NCBI Taxonomy ID is fine; set 'fromdata' to get it from reference). See perldoc
        -abbreviate     If -abbreviate is given then and the -organism is an integer then it is assumed it is an NCBI taxonomy ID and is translated to an abbreviation
	-delimiter:s 	Delimiter of files (incl. physical properties file). Defaults to tab

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

=head1 AUTHORS

 Alexie Papanicolaou 1 2

	1 Max Planck Institute for Chemical Ecology, Germany
	2 Centre for Ecology and Conservation, University of Exeter, UK
	alexie@butterflybase.org

=head1 TIPS

 a dsn file is simple text file which has one line like this:
 dbi:Pg:dbname=InsectaCentral;host=localhost;port=15432;user=postgres;password=123321

=head1 DISCLAIMER & LICENSE

 This software is released under the GNU General Public License version 3 (GPLv3).
 It is provided "as is" without warranty of any kind. You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html. Please note that incorporating the whole software or parts of its code in proprietary software is prohibited under the current license.

=head1 BUGS & LIMITATIONS

 None known. Please report them.

 You will need a preloaded chado database in order to get get the correct CV terms for GO

 An example of a DSN file:

 dbi:Pg:dbname=chado;host=localhost;port=5433;user=chado_user;password=123

=cut

use strict;
use Getopt::Long;		
use Pod::Usage;
use Bio::Tools::GFF;
use Bio::Taxon;
use Bio::DB::Taxonomy;
use DBI;
use Time::Progress;
$| = 1;
use Data::Dumper;

#Command line options
my $gffVersion = 3;
my ($dsn,$id_prefix,$user_id_prefix);
my $blastformat="blast";
my $cutoff_score=1;
my $cutoff_eval=100;
my $delimiter="\t";
my $source="ANNOT8R";
my $dbname="local";
my $flatfile_dir = '/db/ncbi_taxonomy/';
my $physprop;
my (%links,%go_terms_chado,$nodb,$organism,$iprflag, $outdir,$abbreviate);
GetOptions(
	'nodb'					=> \$nodb,
	'o|organism:s'			=> \$organism,
	'dsn|database_connection:s'	=> \$dsn,
	'id:s'					=> \$user_id_prefix,		# for the unique id code; such as GO EC
	'physprop:s'			=> \$physprop,
	'gff_version:i'		=> \$gffVersion,
	'source:s'				=> \$source,
	'cs|cutoff_score:i'	=> \$cutoff_score,
	'ce|cutoff_eval:s'		=> \$cutoff_eval,
	'outdir:s'          => \$outdir,
	'abbreviate'       => \$abbreviate,
        'f|taxonomy_flatfile_dir' => \$flatfile_dir,
	'delimiter:s'	=> \$delimiter,

);
my @infiles=@ARGV;
if (!@infiles){pod2usage;die("No input"); };
unless ($nodb ||$dsn){die "You must specify a database connection or specifically give the -nodb option.\n"}
foreach my $file (@infiles){if ($file=~/iprscan/ || $file=~/\.raw$/){$iprflag=1;}}
unless ($iprflag){if (!$physprop || !-s $physprop ){die "Couldn't find physprop file\n";}}
if ($user_id_prefix){$id_prefix=$user_id_prefix;}
my %organisms;
if (!$organism){
    warn "No organism provided. Are you sure you will not database the results?\n";sleep(5);
}elsif (-s $organism ){
	# then organism is retrieved from a tab delimited list. First and last element is taxid and abbreviation respectively
	open( ORGANISM, $organism );
	while ( my $line = <ORGANISM> ) {
		if ( $line =~ /^#|^\s*$/ ) { next; }
		chomp($line);
		my @data = split( '\t', $line );
		$organisms{ $data[0] } = $data[-1];
	}
}elsif ($organism=~/^\d+$/ && $abbreviate){
    $organism=&get_species_abbr($organism);
}

#other globals
# We need a DSN to find the type which is needed to find the next unique feature. But no necessary if $nodb
my ($dbh,$prepared_id,$prepared_terms );
if ($dsn && -f $dsn){
	open (IN,$dsn);
	$dsn="";
	while (my $line=<IN>){unless ($line=~/^\s*$/){chomp($line);$dsn.=$line;}}
	close (IN);
}
unless ($nodb){
	$dbh=&connecttodb($dsn) || die();
	&add_pg_functions($dbh);
	my $sql_terms="select name,cvterm_id as id from cvterm where cv_id = (select cv_id from cv where name=? ) and dbxref_id = (select dbxref_id from dbxref where accession=? and db_id =(select db_id from db where name=?) order by version limit 1)";
# probably add below "and cv_id = (select cv_id from cv where name=? )"
	my $sql_id="select cvterm_id as id, name from cvterm where name=?";
	#my $prepared_terms = $dbh->prepare_cached($sql_terms);
	$prepared_terms = $dbh->prepare($sql_terms);
	$prepared_id = $dbh->prepare_cached($sql_id);
}

# build global %links
if ($physprop){
print "Processing physical properties\n";
    &process_physprop($physprop);
}else{
	warn "Warning physical properties file not processed. Annot8r processing may fail\n";
}


foreach my $infile (@infiles){
	unless (-s $infile ){die "Couldn't find $infile \n";}
	#get cvterm foreach unique go term reported. but only if we have a go file...
	my $gff =`basename $infile`;chomp($gff);
	$gff.=".gff";
	if ($outdir){$gff=$outdir.'/'.$gff;}
	if (-s $gff){
                warn "Outfile $gff already exists. Skipping\n";
                next;
        }

	print "\nProcessing $infile as $gff\n";
	if ($infile=~/go/i){
# removed because there is no need to get the cvterms while creating the GFF, we should get them prior loading to chado
#		unless ($nodb){
#			print "Finding names and internal IDs for GO terms\n";
#			%go_terms_chado=%{&get_goterms($infile)};
#		}
		print "Processing...\n";
		&processGO($infile, $gff);
	}
	elsif ($infile=~/ec/i){
		print "Processing...\n";
		&processEC($infile, $gff);
	}
	elsif ($infile=~/kegg/i){
		print "Processing...\n";
		&processKEGG($infile,$gff);
	}
	elsif ($infile=~/iprscan/i){
		print "Processing...\n";
		&process_iprscan($infile,$gff);
	}
}
unless ($nodb){
	$prepared_id->finish;
	$prepared_terms->finish;
	&disconnectdb($dbh);
}
print "\nCompleted.\n";

########################################################
sub get_type($$$){
	# we are using this to mainly get the next serial ID TODO: Check if necessary since we have the name
	# and also to get the cvterms of the GO terms. (necessary at this step? due to -> ic_chado_loadcv.pl)
	my $ontology_name=shift;
	my $ontology_abbrv=shift;
	my $ontology_accession=shift || die();
	my $results_hash;
	if ($ontology_accession=~/^\d+$/){
	#this is for GO only
		#print "$ontology_name,$ontology_accession,$ontology_abbrv\n";
		$prepared_terms->execute($ontology_name,$ontology_accession,$ontology_abbrv) || die("Error: Unable to query database for $ontology_name $ontology_abbrv,$ontology_accession\n");
		$results_hash = $prepared_terms->fetchrow_hashref;
	}
	else {
		$prepared_id->execute($ontology_accession);
		$results_hash = $prepared_id->fetchrow_hashref;
	}
	my $term_id=$results_hash->{'id'};
	my $term_name=$results_hash->{'name'};
	if (!$term_id || !$term_name){warn ("Failed to get ID and/or NAME for $ontology_name $ontology_abbrv $ontology_accession\nHave you looaded $ontology_name? Perhaps it has been made obsolete\n");}
	return ($term_id,$term_name);
}

sub process_iprscan ($$){
	my $infile = shift;
	my $gff = shift;
	
	my $serial; $serial=int(1);
	my $gffWriter = Bio::Tools::GFF->new(-file => ">$gff",-gff_version => $gffVersion);
	my $counter=int(0);
	my $total_queries=`wc -l $infile`;
	if($total_queries=~/(\d+)/){$total_queries=$1;}
	print "There are $total_queries queries to be processed, please wait...\n";
	my $timer = new Time::Progress;
	$timer->attr( min => 0, max => $total_queries );
	$timer->restart;
	open (ANNOT8RIPR,$infile);
	while (my $line=<ANNOT8RIPR>){
		$counter++;
		if ($counter=~/000$/){print $timer->report( "eta: %E min, %40b %p\r", $counter );}
		if($line!~/^#/ && $line!~/^\s*$/){
			chomp($line);
			my @data=split($delimiter,$line);
			#friendly names
			my $ref=$data[0];
			my $ref_length=$data[2];
			my $ipr_alg=$data[3];						#store : algorithm

			# store : dbxref http://www.ebi.ac.uk/interpro/ISignature?ac=
			my $ipr_domid=$data[4] if ($ipr_alg!~/coil|seg|SignalP|TMHMM/i);	
			my $ipr_domid_pure=$ipr_domid;
			if ($ipr_domid){
				#alg	exampleid	db id in chado	link (unless otherw stated, includes full id)
				#BlastProDom	PD021112	PRODOM	http://prodom.prabi.fr/prodom/current/cgi-bin/request.pl?question=SPTR&query=
				#superfamily	SSF55874	SUPERFAMILY	http://www.ebi.ac.uk/interpro/ISignature?ac= or supfam.org/SUPERFAMILY/cgi-bin/scop.cgi?sunid= (without SSF)
				#ProfileScan	PS50126	PROSITE (new)	http://www.expasy.ch/prosite/
				#PatternScan	PS00732	PROSITE	http://www.expasy.ch/prosite/
				#FPrintScan	PR00681	PRINTS	http://www.bioinf.manchester.ac.uk/cgi-bin/dbbrowser/sprint/searchprintss.cgi?display_opts=Prints&category=None&queryform=false&regexpr=off&prints_accn=
				#HMMSmart	SM00316	SMART	http://smart.embl-heidelberg.de/smart/do_annotation.pl?DOMAIN=
				#Gene3D	G3DSA:2.40.50.140	:CATH (NEW) http://www.cathdb.info/cathnode/2.40.50.140 (no g3dsa)
				#HMMTigr	TIGR00002	TIGR	http://cmr.jcvi.org/tigr-scripts/CMR/HmmReport.cgi?hmm_acc=
				#HMMPfam	PF00672	PFAM	www.sanger.ac.uk/cgi-bin/Pfam/getacc?
				#HMMPanther	PTHR23283	Panther (new)	http://www.pantherdb.org/panther/family.do?clsAccession=
				#HMMPIR	PIRSF002148	PIRSF (new)	http://pir.georgetown.edu/cgi-bin/ipcSF?id=
				if ($ipr_alg=~/BlastProDom/i ){$ipr_domid="PRODOM:".$ipr_domid;}
				if ($ipr_alg=~/superfamily/i ){$ipr_domid="SUPERFAMILY:".$ipr_domid;}
				if ($ipr_alg=~/ProfileScan/i ){$ipr_domid="PROSITE:".$ipr_domid;}
				if ($ipr_alg=~/PatternScan/i ){$ipr_domid="PROSITE:".$ipr_domid;}
				if ($ipr_alg=~/FPrintScan/i ){$ipr_domid="PRINTS:".$ipr_domid;}
				if ($ipr_alg=~/HMMSmart/i ){$ipr_domid="SMART:".$ipr_domid;}
				if ($ipr_alg=~/Gene3D/i ){
					$ipr_domid=~s/G3DSA://;
					$ipr_domid=":CATH:".$ipr_domid;
				}
				if ($ipr_alg=~/HMMTigr/i ){$ipr_domid="TIGR:".$ipr_domid;}
				if ($ipr_alg=~/HMMPfam/i ){$ipr_domid="PFAM:".$ipr_domid;}
				if ($ipr_alg=~/HMMPanther/i ){$ipr_domid="Panther:".$ipr_domid;}
				if ($ipr_alg=~/HMMPIR/i){$ipr_domid="PIRSF:".$ipr_domid;}
			}
			my ($ipr_best_hit,$ipr_best_hit_desc);
			# sometimes, there's a description.
			if ($data[5]!~/no description/i || $data[5]!~/transmembrane.+region/i || 
				$data[5] !~/coiled-coil/|| $data[5] ne "seg" || $data[5] !~/signal-peptide/i){
					$ipr_best_hit=$data[4];
					$ipr_best_hit_desc=$data[5];
					$ipr_best_hit_desc=~s/\;//;
					$ipr_best_hit_desc=ucfirst($ipr_best_hit);
				}

			my $ipr_qstart=$data[6];
			my $ipr_qend=$data[7];
			my $ipr_best_eval=$data[8] if ($data[8] ne "NA");
#			if (!$ipr_best_eval){$ipr_best_eval=int(0);}

			my $date=$data[10];
			my $ipr_id=$data[11] if ($data[11]&& $data[11] ne "NULL");		# where avail. db xref: http://www.ebi.ac.uk/interpro/ISignature?ac=
			my $ipr_desc=$data[12] if ($data[12] && $data[12] ne "NULL");	# quite long. has GO terms where available

			if ($ipr_best_eval && $ipr_best_eval > $cutoff_eval){next;}

			my $length=$ref_length;
			my $name;
			my $local_source="IPRSCAN";
			if ($ipr_alg=~/BlastProDom|HMMPfam|HMMPIR|PatternScan|ProfileScan|FPrintScan|superfamily|HMMTigr|HMMPanther|HMMSmart|Gene3D/i){
				$name="supported_by_sequence_similarity";
			}
			elsif ($ipr_alg=~/Coil|Seg|SignalPHMM|TMHMM/i){$name="predicted_by_ab_initio_computation";}
			else{$name="protein_match";}
			unless($user_id_prefix){$id_prefix="IPRSCAN";}
			if (!$nodb && $dbh && $serial==1){
				my ($type,$name)=&get_type("sequence","SO",$name);
				$serial = &calculateUniqueSubfeatureId($ref,$type);
			}
			my $gff_id=$data[0].":".$id_prefix.":".$serial;
					
			my $gff_ipr_obj = Bio::SeqFeature::Generic->new(
				-seq_id		=> $ref
				,-source_tag	=> $local_source
				,-primary	=> $name
				,-start		=> $ipr_qstart
				,-end		=> $ipr_qend
				,-strand	=> "."
				,-frame		=> "."
				,-tag		=> {
					ID 	=> $gff_id
					,Name	=> $gff_id
					}
				);
				$gff_ipr_obj->score($ipr_best_eval) if ($ipr_best_eval);
#no UNIPROT id... $gff_ipr_obj->add_tag_value('Dbxref','UNIPROT:'.$ipr_best_hit) if ($ipr_best_hit);
				$gff_ipr_obj->add_tag_value('Alias',$ipr_domid_pure) if ($ipr_domid_pure);
				$gff_ipr_obj->add_tag_value('Dbxref',"InterPro:".$ipr_id) if ($ipr_id);
				$gff_ipr_obj->add_tag_value('Dbxref',$ipr_domid) if ($ipr_domid); #if $ipr_alg!~/coil|seg|SignalP|TMHMM/i
				$gff_ipr_obj->add_tag_value('Note',$ipr_desc) if ($ipr_desc);
				$gff_ipr_obj->add_tag_value('algorithm',$ipr_alg) if ($ipr_alg);
				$gff_ipr_obj->add_tag_value('date',$date) if ($date);
			if (    $organism
				 && -f $organism
				 && $ref =~ /^[A-Z]{2}(\d+)/ )
			{
				my $taxid = $1;
				$gff_ipr_obj->add_tag_value( 'organism', $organisms{$taxid} )
				  if $organisms{$taxid};
			} elsif($organism eq 'fromdata'){
				$ref =~ /^[A-Z]{2}(\d+)/;
				$gff_ipr_obj->add_tag_value( 'organism',$1) if $1;
			} elsif ($organism) {
				$gff_ipr_obj->add_tag_value('organism', $organism );
			}
				$gffWriter->write_feature($gff_ipr_obj);
			$serial++;
		}
	}	
	close (ANNOT8RIPR);
}

sub get_goterms ($){
	# see ic_chado_loadcv.pl
	my $infile=shift;
	my %hash;
	my @go_results=`egrep -o \"\\bGO:[0-9]+\\s+[PFC]\\b\" $infile`;chomp(@go_results);
	foreach my $result (@go_results){
		my @data=split ($delimiter,$result);
		$data[0]=~s/GO://;
		$hash{$data[0]}{"occurance"}++;
		$hash{$data[0]}{"pfc"}=$data[1];
	}
	foreach my $go (keys %hash){
		if ($hash{$go}{"pfc"} eq "C"){
			($hash{$go}{"cvterm_id"},$hash{$go}{"cvterm_name"})=&get_type("cellular_component","GO",$go);
			if (!$hash{$go}{"cvterm_id"} || !$hash{$go}{"cvterm_name"}){($hash{$go}{"cvterm_id"},$hash{$go}{"cvterm_name"})=($go,"obsoleted GO term");}
		}
		elsif ($hash{$go}{"pfc"} eq "F"){
			($hash{$go}{"cvterm_id"},$hash{$go}{"cvterm_name"})=&get_type("molecular_function","GO",$go); 
			if (!$hash{$go}{"cvterm_id"} || !$hash{$go}{"cvterm_name"}){($hash{$go}{"cvterm_id"},$hash{$go}{"cvterm_name"})=($go,"obsoleted GO term");}
		}
		elsif ($hash{$go}{"pfc"} eq "P"){
			($hash{$go}{"cvterm_id"},$hash{$go}{"cvterm_name"})=&get_type("biological_process","GO",$go); 
			if (!$hash{$go}{"cvterm_id"} || !$hash{$go}{"cvterm_name"}){($hash{$go}{"cvterm_id"},$hash{$go}{"cvterm_name"})=($go,"obsoleted GO term");}
		}
		
	}
	return \%hash;
}
sub process_physprop ($){
	my $infile = shift;
	my $name="predicted_by_ab_initio_computation";
	
	my $local_source=$source."_physprop";	
	#id: SO:0000911
	#name: predicted_by_ab_initio_computation
	#def: "An attribute describing a feature that is predicted by a computer program that did not rely on sequence similarity." [SO:ke]

	my $gff=$infile.".gff";
	my $gffWriter = Bio::Tools::GFF->new(-file => ">$gff",-gff_version => $gffVersion);
	my $no_gff;
        if (-s $gff){warn "Outfile exists, skipping...\n";$no_gff=1;}
	open (CSV,$infile);
	while (my $line=<CSV>){
		if($line!~/^#/ && $line!~/^\s*$/){
			chomp($line);
			my @data=split($delimiter,$line);

			$links{$data[0]}={
				"length"		=>$data[1],
				"unknown_aa"	=>$data[2],
				"minwm"			=>$data[3],
				"maxmw"			=>$data[4],
				"charge_pos"	=>$data[5],
				"charge_neg"	=>$data[6],
				"pI"			=>$data[7],
				"charge_ph5"	=>$data[8],
				"charge_ph7"	=>$data[9],
				"charge_ph8"	=>$data[10]
			};
			next if $no_gff;
			my $gff_id=$data[0].":physprop";
			#build gff
			my $gff_phys_obj = Bio::SeqFeature::Generic->new(
				-seq_id			=> $data[0]
				,-source_tag	=> $local_source
				,-primary		=> $name
				,-start			=> "1"
				,-end			=> $data[1]		# because we don't really want to bother(index the blast report)
				,-strand		=> "."
				,-frame			=> "."
				,-tag			=> {
					ID 		=> $gff_id
					,Name	=> $gff_id
				});
			$gff_phys_obj->add_tag_value('molecular_weight',"$data[3],$data[4]");
			$gff_phys_obj->add_tag_value('charged',"+$data[5],-$data[6]");
			$gff_phys_obj->add_tag_value('pI',"$data[6]");
			$gff_phys_obj->add_tag_value('charge pH5',"$data[8]");
			$gff_phys_obj->add_tag_value('charge pH7',"$data[9]");
			$gff_phys_obj->add_tag_value('charge pH8',"$data[10]");
			if (    $organism
				 && -f $organism
				 && $data[0] =~ /^[A-Z]{2}(\d+)/ )
			{
				my $taxid = $1;
				$gff_phys_obj->add_tag_value( 'organism', $organisms{$taxid} )
				  if $organisms{$taxid};
			} elsif($organism eq 'fromdata'){
				$data[0] =~ /^[A-Z]{2}(\d+)/;
				$gff_phys_obj->add_tag_value( 'organism',$1) if $1;
			} elsif ($organism) {
				$gff_phys_obj->add_tag_value( 'organism', $organism );
			}
			$gffWriter->write_feature($gff_phys_obj);
		}
	}
	close (CSV);	
}

sub processGO ($$){
	my $infile = shift;
	my $gff = shift;
	my $local_source=$source."_GO";
	my $name="supported_by_sequence_similarity";
	unless($user_id_prefix){$id_prefix="GO";}

	my $serial; $serial=int(1);
	my $gffWriter = Bio::Tools::GFF->new(-file => ">$gff",-gff_version => $gffVersion);
	my $counter=int(0);
	my $total_queries=`wc -l $infile`;
	if($total_queries=~/(\d+)/){$total_queries=$1;}
	print "There are $total_queries queries to be processed, please wait...\n";
	my $timer = new Time::Progress;
	$timer->attr( min => 0, max => $total_queries );
	$timer->restart;
	open (ANNOT8RGO,$infile);
	while (my $line=<ANNOT8RGO>){
		$counter++;
		if ($counter=~/000$/){print $timer->report( "eta: %E min, %40b %p\r", $counter );}
		if($line!~/^#/ && $line!~/^\s*$/){
			chomp($line);
			my @data=split($delimiter,$line);
			#friendly names
			my $ref=$data[0];
			my $go_id=$data[1];
			$go_id=~s/GO://;
			my $go_pcf=$data[2];
			my $go_desc=$data[3];
			$go_desc=~s/\"|\.//g;
			my $go_slim=$data[4];
			$go_slim=~s/GO://;
			my $go_best_hit=$data[5];
			my $go_best_score=$data[6];
			my $go_best_eval=$data[7];
			my $go_no_hits=$data[8];
			my $go_no_ids=$data[9];
			my $go_stat=$data[10];
			if ($go_best_eval > $cutoff_eval || $go_best_score < $cutoff_score){next;}
			my $length=$links{$ref}{"length"} || die();


			if (!$nodb && $dbh && $serial==1){
				# we just want the type but nevermind
				my ($type,$name)=&get_type("sequence","SO",$name);
				$serial = &calculateUniqueSubfeatureId($ref,$type);
			}

			my $gff_id=$data[0].":".$id_prefix.":".$serial;
					
			my $gff_go_obj = Bio::SeqFeature::Generic->new(
				-seq_id		=> $ref
				,-source_tag	=> $local_source
				,-primary	=> $name
				,-start		=> "1"
				,-end		=> $length
				,-score		=> $go_best_eval
				,-strand	=> "."
				,-frame		=> "."
				,-tag		=> {
					ID 	=> $gff_id
					,Name	=> $gff_id
					}
				);
#				unless ($nodb){
#					my $cvterm=$go_terms_chado{$go_id}{"cvterm_id"} || warn ("Couldn't find cvterm for $go_id\n");
#					$gff_go_obj->add_tag_value('Dbxref',"local:".$cvterm) if $cvterm;
#				}
				$gff_go_obj->add_tag_value('Alias',$go_id) if $go_id;
				$gff_go_obj->add_tag_value('Alias',$go_slim) if ($go_slim && $go_slim ne $go_id);
				$gff_go_obj->add_tag_value('Note',$go_desc) if $go_desc;
			# the gene we took the annotation from
				$gff_go_obj->add_tag_value('Dbxref',"UNIPROT:".$go_best_hit) if $go_best_hit;
				$gff_go_obj->add_tag_value('Dbxref',"GO:".$go_id) if $go_id;
				$gff_go_obj->add_tag_value('Dbxref',"GO:".$go_slim) if $go_slim;
				if (    $organism
					 && -f $organism
					 && $ref =~ /^[A-Z]{2}(\d+)/ )
				{
					my $taxid = $1;
					$gff_go_obj->add_tag_value( 'organism', $organisms{$taxid} )
					  if $organisms{$taxid};
				} elsif($organism eq 'fromdata'){
					$ref =~ /^[A-Z]{2}(\d+)/;
					$gff_go_obj->add_tag_value( 'organism',$1) if $1;
				} elsif ($organism) {
					$gff_go_obj->add_tag_value('organism', $organism );
				}
				$gffWriter->write_feature($gff_go_obj);
			$serial++;
		}
	}	
	close ANNOT8RGO;
}

sub processKEGG ($$){
	# here we are not interested in using a preloaded CV as it will be autocreated for us.
	my $infile = shift;
	my $gff = shift;
	my $name="supported_by_sequence_similarity";
	my $local_source=$source."_KEGG";
	unless($user_id_prefix){$id_prefix="KEGG";}

	my $serial; $serial=int(1);
	my $gffWriter = Bio::Tools::GFF->new(-file => ">$gff",-gff_version => $gffVersion);
	my $counter=int(0);
	my $total_queries=`wc -l $infile`;
	if($total_queries=~/(\d+)/){$total_queries=$1;}
	print "There are $total_queries queries to be processed, please wait...\n";
	my $timer = new Time::Progress;
	$timer->attr( min => 0, max => $total_queries );
	$timer->restart;
	open (ANNOT8RKEGG,$infile);
	while (my $line=<ANNOT8RKEGG>){
		$counter++;
		if ($counter=~/000$/){print $timer->report( "eta: %E min, %40b %p\r", $counter );}
		if($line!~/^#/ && $line!~/^\s*$/){
			chomp($line);
			my @data=split($delimiter,$line);
			#friendly names
			my $ref=$data[0];
			my $kegg_id=$data[1];
			my $kegg_path=$data[2];
			my $kegg_desc=$data[3];
			$kegg_desc=~s/\"|\.//g;
			my $kegg_best_hit=$data[4];
			my $kegg_best_score=$data[5];
			my $kegg_best_eval=$data[6];
			my $kegg_no_hits=$data[7];
			my $kegg_no_ids=$data[8];
			my $kegg_stat=$data[9];
#debugif($data[10]){	#die "ref $ref kegg_id $kegg_id path $kegg_path desc: $kegg_desc hit: $kegg_best_hit eval: $kegg_best_eval score: $kegg_best_score hits: $kegg_no_hits ids: $kegg_no_ids stat: $kegg_stat\n";	next;}
			if ($kegg_best_eval > $cutoff_eval || $kegg_best_score < $cutoff_score){next;}
			my $length=$links{$ref}{"length"} || die();
			
			if (!$nodb && $dbh && $serial==1){
				my ($type,$name)=&get_type("sequence","SO",$name);
				$serial = &calculateUniqueSubfeatureId($ref,$type);
			}
			
			my $gff_id=$data[0].":".$id_prefix.":".$serial;
					
			my $gff_csv_obj = Bio::SeqFeature::Generic->new(
				-seq_id		=> $ref
				,-source_tag	=> $local_source
				,-primary	=> $name
				,-start		=> "1"
				,-end		=> $length
				,-score		=> $kegg_best_eval
				,-strand	=> "."
				,-frame		=> "."
				,-tag		=> {
					ID 	=> $gff_id
					,Name	=> $gff_id
					}
				);
				$gff_csv_obj->add_tag_value('Dbxref',"KEGG_ORTHOLOGY:".$kegg_id) if $kegg_id;
				$gff_csv_obj->add_tag_value('Dbxref',"KEGG_PATHWAY:".$kegg_path) if $kegg_path;
			# the gene we took the annotation from
				$gff_csv_obj->add_tag_value('Dbxref',"UNIPROT:".$kegg_best_hit) if $kegg_best_hit;
				$gff_csv_obj->add_tag_value('Note',$kegg_desc) if $kegg_desc;
				$gff_csv_obj->add_tag_value('Alias',$kegg_id) if $kegg_id;
				if (    $organism
					 && -f $organism
					 && $ref =~ /^[A-Z]{2}(\d+)/ )
				{
					my $taxid = $1;
					$gff_csv_obj->add_tag_value( 'organism', $organisms{$taxid} )
					  if $organisms{$taxid};
				} elsif($organism eq 'fromdata'){
					$ref =~ /^[A-Z]{2}(\d+)/;
					$gff_csv_obj->add_tag_value( 'organism',$1) if $1;
				} elsif ($organism) {
					$gff_csv_obj->add_tag_value('organism', $organism );
				}
				$gffWriter->write_feature($gff_csv_obj);
			$serial++;
		}
	}
	close ANNOT8RKEGG;
}

sub processEC ($$){
	# here we are not interested in using a preloaded CV as it will be autocreated for us.

	my $infile = shift;
	my $gff = shift;
	my $name="supported_by_sequence_similarity";

	my $local_source=$source."_EC";
	unless($user_id_prefix){$id_prefix="EC";}

	my $serial; $serial=int(1);
	my $gffWriter = Bio::Tools::GFF->new(-file => ">$gff",-gff_version => $gffVersion);
	my $counter=int(0);
	my $total_queries=`wc -l $infile`;
	if($total_queries=~/(\d+)/){$total_queries=$1;}
	print "There are $total_queries queries to be processed, please wait...\n";
	my $timer = new Time::Progress;
	$timer->attr( min => 0, max => $total_queries );
	$timer->restart;
	open (ANNOT8REC,$infile);
	while (my $line=<ANNOT8REC>){
		$counter++;
		if ($counter=~/000$/){print $timer->report( "eta: %E min, %40b %p\r", $counter );}
		if($line!~/^#/ && $line!~/^\s*$/){
			chomp($line);
			my @data=split($delimiter,$line);
			#friendly names
			my $ref=$data[0];
			my $ec_id=$data[1];
			my $ec_best_desc=$data[2];
			$ec_best_desc=~s/\"|\.//g;
			my $ec_best_hit=$data[3];
			my $ec_best_score=$data[4];
			my $ec_best_eval=$data[5];
			my $ec_no_hits=$data[6];
			my $ec_no_ids=$data[7];
			my $ec_stat=$data[8];
		
			if ($ec_best_eval > $cutoff_eval || $ec_best_score < $cutoff_score){next;}
			my $length=$links{$ref}{"length"} || die("Cannot find length of reference $ref\n");

			if (!$nodb && $dbh && $serial==1){
					my ($type,$name)=&get_type("sequence","SO",$name);
					$serial = &calculateUniqueSubfeatureId($ref,$type);
			}
			my $gff_id=$data[0].":".$id_prefix.":".$serial;
					
			my $gff_csv_obj = Bio::SeqFeature::Generic->new(
				-seq_id		=> $ref
				,-source_tag	=> $local_source
				,-primary	=> $name
				,-start		=> "1"
				,-end		=> $length
				,-score		=> $ec_best_eval
				,-strand	=> "."
				,-frame		=> "."
				,-tag		=> {
					ID 	=> $gff_id
					,Name	=> $gff_id
					}
				);
				$gff_csv_obj->add_tag_value('Note',$ec_best_desc) if $ec_best_desc;
			# the gene we took the annotation from
				$gff_csv_obj->add_tag_value('Dbxref',"UNIPROT:".$ec_best_hit) if $ec_best_hit;
				$gff_csv_obj->add_tag_value('Dbxref',"EC:".$ec_id) if $ec_id;
				$gff_csv_obj->add_tag_value('Alias',$ec_id) if $ec_id;
				if ($organism
					 && -f $organism
					 && $ref =~ /^[A-Z]{2}(\d+)/ ){
					my $taxid = $1;
					$gff_csv_obj->add_tag_value( 'organism', $organisms{$taxid} )
					if $organisms{$taxid};
				} elsif($organism eq 'fromdata'){
					$ref =~ /^[A-Z]{2}(\d+)/;
					$gff_csv_obj->add_tag_value( 'organism',$1) if $1;
				} elsif ($organism) {
					$gff_csv_obj->add_tag_value('organism',$organism);
				}
				$gffWriter->write_feature($gff_csv_obj);
			$serial++;
		}
	}
	close (ANNOT8REC);
}

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

sub calculateUniqueSubfeatureId($$) {
	my $srcfeature = shift;
	my $matchgroup = shift;

	#get last subfeature index available and all the other ones continue from there (so only query db once)
	my $sql = "SELECT hit.uniquename AS id";
	$sql .= " FROM feature AS hit, feature AS src, featureloc AS loc";
	$sql .= " WHERE loc.feature_id = hit.feature_id";
	$sql .= " AND loc.srcfeature_id = src.feature_id";
	$sql .= " AND src.uniquename = '".$srcfeature."'";
	$sql .= " AND hit.type_id = '".$matchgroup."'";
	$sql .= " ORDER BY ic_accession_serial(hit.uniquename) DESC";
	$sql .= " LIMIT 1";
	my $sth = $dbh->prepare($sql);
	
	$sth->execute or die("Error: Unable to query database");
	my $vals = $sth->fetchrow_hashref;
	#die ("val from $sql is $vals");
	my $next_id=int(1);
	if ($vals->{'id'}){
		if ($vals->{'id'} =~ /([0-9+])$/){
			$next_id=$1+1;
		}
	}
	
	$sth->finish;
	return $next_id ;
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

