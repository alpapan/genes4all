#!/usr/bin/perl -w
#$Id$

=pod

=head1 NAME

 extract_experiment.pl
 
=head1 USAGE

 -dsn => DSN connection string. e.g. dbi:Pg:dbname=drupal;host=localhost;port=5432;user=www-db;password=123
 -passkey => Extract passkey as well

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


=cut

use strict;
use DBI;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
my $dbprefix='';
my ($unfinished,$get_passkey,$help);
my $dsn =  'dbi:Pg:dbname=drupal;host=localhost;port=5432;user=www-db;password=123';
GetOptions(
    'dsn:s' => \$dsn,
    'passkey'=>\$get_passkey,
    'help'  => \$help,
    'prefix'=> \$dbprefix,
    'unifinished'=> \$unfinished,
);
if (!$dsn || $help){pod2usage}
my $dbh = &connecttodb($dsn);
## prepare SQL commands
#features. removed tax_class
my $get_feature_sql =
'SELECT feature.uniquename,feature_id,feature.dbxref_id,dbxref.accession,db.name as db_name,residues,feature.organism_id,tax_order,tax_family,'
  . 'genus,species,ncbi_taxid,cvterm.name as type_name from '.$dbprefix.'gmod_dbsf_feature as feature JOIN '
  . ' '.$dbprefix.'gmod_dbsf_dbxref as dbxref ON dbxref.dbxref_id=feature.dbxref_id JOIN '.$dbprefix.'gmod_dbsf_db as db ON db.db_id=dbxref.db_id '
  . ' JOIN '.$dbprefix.'gmod_dbsf_organism as organism ON organism.organism_id=feature.organism_id '
  . ' JOIN '.$dbprefix.'gmod_dbsf_cvterm as cvterm ON cvterm.cvterm_id=feature.type_id where feature_id=?';
my $get_featurepro_sql =
'SELECT cvterm.name as cvterm_name,value from '.$dbprefix.'gmod_dbsf_featureprop as featureprop '
  . ' JOIN '.$dbprefix.'gmod_dbsf_cvterm as cvterm ON featureprop.type_id=cvterm.cvterm_id'
  . ' WHERE featureprop.feature_id=?';
my $get_feature_cvterms_sql =
'SELECT cvterm.name as cvterm_name,cv.name as cv_name from '.$dbprefix.'gmod_dbsf_feature_cvterm as feature_cvterm '
  . ' JOIN '.$dbprefix.'gmod_dbsf_cvterm as cvterm ON feature_cvterm.cvterm_id=cvterm.cvterm_id'
  . ' JOIN '.$dbprefix.'gmod_dbsf_cv as cv ON cv.cv_id=cvterm.cv_id '
  . ' WHERE feature_cvterm.feature_id=?';
my $get_feature_prepare         = $dbh->prepare($get_feature_sql);
my $get_feature_cvterms_prepare = $dbh->prepare($get_feature_cvterms_sql);
my $get_featurepro_prepare      = $dbh->prepare($get_featurepro_sql);


#resources, removed tax_class
my $get_resource_sql =
'SELECT resource.uniquename,resource_id,resource.dbxref_id,dbxref.accession,db.name as db_name,resource.organism_id,tax_order,tax_family,'
  . 'genus,species,ncbi_taxid,cvterm.name as type_name from '.$dbprefix.'gmod_dbsf_resource as resource JOIN '
  . ' '.$dbprefix.'gmod_dbsf_dbxref as dbxref ON dbxref.dbxref_id=resource.dbxref_id JOIN '.$dbprefix.'gmod_dbsf_db as db ON db.db_id=dbxref.db_id '
  . ' JOIN '.$dbprefix.'gmod_dbsf_organism as organism ON organism.organism_id=resource.organism_id '
  . ' JOIN '.$dbprefix.'gmod_dbsf_cvterm as cvterm ON cvterm.cvterm_id=resource.type_id'
  . ' where resource_id=?';
my $get_resource_cvterms_sql =
'SELECT cvterm.name as cvterm_name,cv.name as cv_name from '.$dbprefix.'gmod_dbsf_resource_cvterm as resource_cvterm '
  . ' JOIN '.$dbprefix.'gmod_dbsf_cvterm as cvterm ON resource_cvterm.cvterm_id=cvterm.cvterm_id'
  . ' JOIN '.$dbprefix.'gmod_dbsf_cv as cv ON cv.cv_id=cvterm.cv_id '
  . ' WHERE resource_cvterm.resource_id=?';
my $get_resource_prop_sql =
'SELECT cvterm.name as cvterm_name,value from '.$dbprefix.'gmod_dbsf_resourceprop as resourceprop '
  . ' JOIN '.$dbprefix.'gmod_dbsf_cvterm as cvterm ON resourceprop.type_id=cvterm.cvterm_id'
  . ' WHERE resourceprop.resource_id=?';
  
my $get_resource_prepare=$dbh->prepare($get_resource_sql);
my $get_resource_cvterms_prepare=$dbh->prepare($get_resource_cvterms_sql);
my $get_resourcepro_prepare =$dbh->prepare($get_resource_prop_sql);

#pub
my $get_pub_sql= 'SELECT uniquename,dbxref.accession as dbxref,db.name as dbname from '.$dbprefix.'gmod_dbsf_pub as pub '
 . ' join '.$dbprefix.'gmod_dbsf_dbxref as dbxref ON dbxref.dbxref_id=pub.dbxref_id JOIN '.$dbprefix.'gmod_dbsf_db as db ON db.db_id=dbxref.db_id'
 . " where pub_id=?";
 my $author_sql    = 'SELECT first_names,last_names,email,rank from '.$dbprefix.'gmod_dbsf_author author JOIN '
 . ' '.$dbprefix.'gmod_dbsf_pub_author pa ON pa.author_id=author.author_id'. " where pub_id=?";
 my $author_prepare=$dbh->prepare($author_sql);
my $get_pub_prepare=$dbh->prepare($get_pub_sql);

# study
my $get_study_sql='SELECT passkey,uniquename,description,cvterm.name as type_name,passkey,submitter_email from '.$dbprefix.'gmod_dbsf_study as study '
.'JOIN '.$dbprefix.'gmod_dbsf_cvterm as cvterm on cvterm.cvterm_id=study.type_id where study_id=?';
my $get_all_studies_sql=  'SELECT study_id as id from '.$dbprefix.'gmod_dbsf_study as study ';
$get_all_studies_sql.= " WHERE type_id!=(select cvterm_id from ".$dbprefix."gmod_dbsf_cvterm where name='unfinished')" if !$unfinished;
my $pub_select = 'SELECT pub_id as id from '.$dbprefix.'gmod_dbsf_study as study where study_id=?';
my $resource_select='select sr.resource_id as id from  '.$dbprefix.'gmod_dbsf_study_resource sr join '.$dbprefix.'gmod_dbsf_resource r ON '
    .' sr.resource_id=r.resource_id join '.$dbprefix.'gmod_dbsf_cvterm cvterm ON type_id=cvterm_id  where cvterm.name=? AND  study_id=?';
my $feature_select='select sr.feature_id as id from  '.$dbprefix.'gmod_dbsf_study_feature sr join '.$dbprefix.'gmod_dbsf_feature r ON sr.feature_id=r.feature_id'
    .' join '.$dbprefix.'gmod_dbsf_cvterm cvterm ON type_id=cvterm_id where cvterm.name=? AND study_id=?';

my $get_study_prepare=$dbh->prepare($get_study_sql);
my $get_all_studies_prepare=$dbh->prepare($get_all_studies_sql);
my $pub_select_prepare=$dbh->prepare($pub_select);
my $resource_select_prepare=$dbh->prepare($resource_select);
my $feature_select_prepare=$dbh->prepare($feature_select);


### end prepare SQL
## CSV
open (CSV, ">experiment_results.csv");
my @headers;
push (@headers,'Passkey') if $get_passkey;
push (@headers,qw/Study_ID Study_Name Study_Submitter_email Study_Type Study_Description/);
push (@headers,qw/Publication_Name Publication_DB_name Publication_DB_accession Publication_Author_Last Publication_Author_First Publication_Author_email/);
push (@headers,qw/Target_Name Target_sequence Target_Tax_order Target_Tax_family Target_Tax_genus Target_Tax_species Target_Tax_NCBI_ID Target_alignment Target_Repository_Name Target_Repository_Accession/);
push (@headers,qw/Construct_Name Construct_Sequence Construct_type Construct_in-vitro Construct_purification_method Construct_annealing_method Construct_protocol Construct_Repository_Name Construct_Repository_ID/);
push (@headers,qw/Animal_Name Animal_Tax_Order Animal_Tax_family Animal_Tax_genus Animal_Tax_species Animal_Tax_NCBI_id Animal_Stock_Center Animal_Stock_ID Animal_Dev_stage Colony_Infection Indiv_Infection Colony_Origin/);
push (@headers,qw/Delivery_Name Delivery_Construct_Concentration Delivery_method Delivery_buffer Delivery_adjuvant_name Delivery_adjuvant_amnt Delivery_Control_method Delivery_replicates Delivery_organ Control_method_detail/);
push (@headers,qw/Assay_Name Assay_dev_stage Assay_detection_method Assay_tissue Assay_processing_detection Assay_detection_time Assay_average_silencing Assay_accurate_silencing/);
#my $csv = Text::CSV_XS->new ({eol => "\n" });
#$csv->print (\*CSV, \@headers);
print CSV join("\t",@headers)."\n";

$get_all_studies_prepare->execute();
my $count=int(0);
my %hash;
while (my $studies_res=$get_all_studies_prepare->fetchrow_hashref){
	$count++;
	my $experiment_id=$studies_res->{'id'};
	my @print_data=($experiment_id);

	#Study
	#Study_Name Study_Type Study_Description
	$get_study_prepare->execute($experiment_id);
	my $study_res=$get_study_prepare->fetchrow_hashref();
	push(@print_data,$study_res->{'passkey'}) if $get_passkey;
	push(@print_data,($study_res->{'uniquename'},$study_res->{'submitter_email'},$study_res->{'type_name'},$study_res->{'description'}));
    $hash{$study_res->{'uniquename'}}=$experiment_id;	
	#Pub
	#Publication_Name Publication_DB_name Publication_DB_accession Publication_Author_Last Publication_Author_First 
	#Publication_Author_email/);
	$pub_select_prepare->execute($experiment_id);
    my $res=$pub_select_prepare->fetchrow_hashref();
    my $pub_id=$res->{'id'};
    my $pub_data_ref=&get_pub($pub_id);
    push (@print_data,($pub_data_ref->{'uname'},$pub_data_ref->{'dbname'},$pub_data_ref->{'accession'},
    $pub_data_ref->{'author_last'},$pub_data_ref->{'author_first'},$pub_data_ref->{'author_email'}));
    
    
    #target
    #Target_Name Target_sequence Target_Tax_order Target_Tax_family 
    #Target_Tax_genus Target_Tax_species Target_Tax_NCBI_ID Target_alignment 
    #Target_Repository_Name Target_Repository_Accession
    $feature_select_prepare->execute('gene_target',$experiment_id);
    $res=$feature_select_prepare->fetchrow_hashref();
   my $target_id=$res->{'id'};
   my $target_data_ref=&get_feature($target_id);
   push(@print_data,($target_data_ref->{'general'}->{'uname'},$target_data_ref->{'general'}->{'sequence'},$target_data_ref->{'general'}->{'tax_order'},$target_data_ref->{'general'}->{'tax_family'},
   $target_data_ref->{'general'}->{'genus'},$target_data_ref->{'general'}->{'species'},$target_data_ref->{'general'}->{'ncbi_taxid'},
   $target_data_ref->{'cvterm'}->{'aln_region'},$target_data_ref->{'general'}->{'db_name'},$target_data_ref->{'general'}->{'accession'}));
   
   
   
   #construct
   #Construct_Name Construct_Sequence Construct_type Construct_in-vitro 
   #Construct_purification_method Construct_annealing_method Construct_protocol
   # Construct_Repository_Name Construct_Repository_ID
    $feature_select_prepare->execute('rnai_construct',$experiment_id);
    $res=$feature_select_prepare->fetchrow_hashref();
    my $construct_id=$res->{'id'};
    my $construct_data_ref=&get_feature($construct_id);
    push(@print_data,($construct_data_ref->{'general'}->{'uname'},$construct_data_ref->{'general'}->{'sequence'},$construct_data_ref->{'cvterm'}->{'rna_probe_type'},
    $construct_data_ref->{'prop'}->{'construct_protocol_kit'},$construct_data_ref->{'cvterm'}->{'purification'},$construct_data_ref->{'cvterm'}->{'annealing_method'},
    $construct_data_ref->{'prop'}->{'construct_protocol_details'},$construct_data_ref->{'general'}->{'db_name'},$construct_data_ref->{'general'}->{'accession'}));
    
   #animals
   #Animal_Name Animal_Tax_Order Animal_Tax_family Animal_Tax_genus 
   #Animal_Tax_species Animal_Tax_NCBI_id Animal_Stock_Center Animal_Stock_ID
   #Animal_Dev_stage Colony_Infection Indiv_infection Colony_Origin 
    $resource_select_prepare->execute('experimental animals',$experiment_id);
    $res=$resource_select_prepare->fetchrow_hashref();
   my $animals_id=$res->{'id'};
   my $animal_data_ref=&get_resource($animals_id);
    push(@print_data,($animal_data_ref->{'general'}->{'uname'},$animal_data_ref->{'general'}->{'tax_order'},$animal_data_ref->{'general'}->{'tax_family'},
    $animal_data_ref->{'general'}->{'genus'},$animal_data_ref->{'general'}->{'species'},$animal_data_ref->{'general'}->{'ncbi_taxid'},
    $animal_data_ref->{'general'}->{'db_name'},$animal_data_ref->{'general'}->{'accession'},$animal_data_ref->{'cvterm'}->{'dev_stage'},
    $animal_data_ref->{'cvterm'}->{'colony_infection'},$animal_data_ref->{'cvterm'}->{'indiv_infection'},
    $animal_data_ref->{'cvterm'}->{'origin'}));

   #delivery
   #Delivery_Name Delivery_Construct_Concentration Delivery_method Delivery_buffer
   # Delivery_adjuvant_name Delivery_adjuvant_amnt Delivery_Control_method Delivery_replicates
    $resource_select_prepare->execute('delivery protocol',$experiment_id);
    $res=$resource_select_prepare->fetchrow_hashref();
   my $delivery_id=$res->{'id'};
   my $delivery_data_ref=&get_resource($delivery_id);
   push(@print_data,($delivery_data_ref->{'general'}->{'uname'},$delivery_data_ref->{'prop'}->{'construct_concentration'},
   $delivery_data_ref->{'cvterm'}->{'delivery_method'},$delivery_data_ref->{'prop'}->{'delivery_buffer'},
   $delivery_data_ref->{'prop'}->{'adjuvant_name'},$delivery_data_ref->{'prop'}->{'adjuvant_amount'},
   $delivery_data_ref->{'cvterm'}->{'control_method'},$delivery_data_ref->{'prop'}->{'replicates'},
   $delivery_data_ref->{'prop'}->{'delivery_organ'},$delivery_data_ref->{'prop'}->{'control_method_detail'}));
   
   
   #assay
   #Assay_Name Assay_dev_stage Assay_detection_method Assay_tissue Assay_processing_detection
   # Assay_detection_time Assay_average_silencing Assay_accurate_silencing
    $resource_select_prepare->execute('assay protocol',$experiment_id);
    $res=$resource_select_prepare->fetchrow_hashref();
    my $assay_id=$res->{'id'};
    my $assay_data_ref=&get_resource($assay_id);
    my $general=$assay_data_ref->{'general'};
    my $cvterm=$assay_data_ref->{'cvterm'};
    my $prop=$assay_data_ref->{'prop'};
    push(@print_data,($general->{'uname'},
    $cvterm->{'dev_stage'},
    $cvterm->{'detection_method'},$prop->{'assay_organ'},
    $cvterm->{'detect_dsrna_processing'},$prop->{'detect_time'},
    $cvterm->{'silencing_level'},$prop->{'silencing_accurate_level'},
    ));
    
    for (my $i=0;$i<scalar(@print_data);$i++){
    	if (ref($print_data[$i]) eq 'ARRAY'){
    		$print_data[$i]=join (' ',@{$print_data[$i]});
    	}
    	$print_data[$i]=~s/\r|\n|\t/ /g if $print_data[$i];
    	if (!$print_data[$i]){$print_data[$i]=' ';}
    	#else{$print_data[$i]="'".$print_data[$i]."'";}
    }
    print CSV join("\t",@print_data)."\n";
}
close CSV;
print "Done. Found $count studies\n";
#print Dumper %hash;

$resource_select_prepare->finish();
$get_resource_prepare->finish();
$get_resource_cvterms_prepare->finish();
$get_resourcepro_prepare ->finish();
$author_prepare->finish();
$get_pub_prepare->finish();
$get_study_prepare->finish();
$get_all_studies_prepare->finish();
$pub_select_prepare->finish();
$feature_select_prepare->finish();
$get_feature_cvterms_prepare->finish();
$get_featurepro_prepare->finish();
$get_feature_prepare->finish();
&disconnectdb($dbh);
######################################################
sub get_feature($) {
	my %data;
	my $feature_id = shift;
	$get_feature_prepare->execute($feature_id);
	
	my $results = $get_feature_prepare->fetchrow_hashref;
	$data{'general'} = {
		'uname'       => $results->{'uniquename'},
		'feature_id'  => $results->{'feature_id'},
		'dbxref_id'   => $results->{'dbxref_id'},
		'accession'   => $results->{'accession'},
		'db_name'     => $results->{'db_name'},
		'sequence'    => $results->{'residues'},
		'organism_id' => $results->{'organism_id'},

		#'tax_class' => $results->{'tax_class'},
		'tax_order'  => $results->{'tax_order'},
		'tax_family' => $results->{'tax_family'},
		'genus'      => $results->{'genus'},
		'species'    => $results->{'species'},
		'ncbi_taxid' => $results->{'ncbi_taxid'},
		'type_name'  => $results->{'type_name'},
	};
	$get_feature_cvterms_prepare->execute($feature_id);
	while ( my $res = $get_feature_cvterms_prepare->fetchrow_hashref ) {
		push( @{ $data{'cvterm'}{ $res->{'cv_name'} } },
			  $res->{'cvterm_name'} );
	}
	$get_featurepro_prepare->execute($feature_id);
	while ( my $res = $get_featurepro_prepare->fetchrow_hashref ) {
		push( @{ $data{'prop'}{ $res->{'cvterm_name'} } }, $res->{'value'} );
	}
	return \%data;
}

sub get_resource($) {
	my %data;
	my $resource_id = shift;
	   $get_resource_prepare->execute($resource_id);
    my $results = $get_resource_prepare->fetchrow_hashref;
    $data{'general'} = {
        'uname'       => $results->{'uniquename'},
        'resource_id'  => $results->{'resource_id'},
        'dbxref_id'   => $results->{'dbxref_id'},
        'accession'   => $results->{'accession'},
        'db_name'     => $results->{'db_name'},
        'sequence'    => $results->{'residues'},
        'organism_id' => $results->{'organism_id'},

        #'tax_class' => $results->{'tax_class'},
        'tax_order'  => $results->{'tax_order'},
        'tax_family' => $results->{'tax_family'},
        'genus'      => $results->{'genus'},
        'species'    => $results->{'species'},
        'ncbi_taxid' => $results->{'ncbi_taxid'},
        'type_name'  => $results->{'type_name'},
    };
    $get_resource_cvterms_prepare->execute($resource_id);
    while ( my $res = $get_resource_cvterms_prepare->fetchrow_hashref ) {
        push( @{ $data{'cvterm'}{ $res->{'cv_name'} } },
              $res->{'cvterm_name'} );
    }
    $get_resourcepro_prepare->execute($resource_id);
    while ( my $res = $get_resourcepro_prepare->fetchrow_hashref ) {
        push( @{ $data{'prop'}{ $res->{'cvterm_name'} } }, $res->{'value'} );
    }
	
	
	
	return \%data;
}


sub get_pub($){
    my $pub_id=shift;
    my %data;
    $get_pub_prepare->execute($pub_id);
    my $pub_res=$get_pub_prepare->fetchrow_hashref;
    $author_prepare->execute($pub_id);
    my $author_res=$author_prepare->fetchrow_hashref;
    %data=(
       'uname'=>$pub_res->{'uniquename'},
       'accession'=>$pub_res->{'dbxref'},
       'dbname'=>$pub_res->{'dbname'},
       'author_last'=>$author_res->{'last_names'},
       'author_first'=>$author_res->{'first_names'},
       'author_email'=>$author_res->{'email'},
#       'author_address'=>$author_res->{'address'},
    );
    return \%data;
}

sub connecttodb() {
	if ( $dsn && -f $dsn ) {
		open( IN, $dsn );
		$dsn = "";
		while ( my $line = <IN> ) {
			unless ( $line =~ /^\s*$/ ) { chomp($line); $dsn .= $line; }
		}
		close(IN);
	}
	if ( !$dsn ) {
		die "No connection to database provided.\n";
	}
	my $dbh = DBI->connect($dsn);
	die("Error: Unable to connect to database") if ( !$dbh );
	return ($dbh);
}

sub disconnectdb ($) { 
	my $dbh = shift;
	$dbh->disconnect();
}
