#!/usr/bin/perl -w
#Changelog
#December 2011: first attempt; works with Rabbit calicivirus
#April 2012: added support for non-genbank features

=pod

=head1 TODO

 # TODO AP low until we get some relevant data Process metadata add_geolocation_metadata() add_activity_metadata();
 # TODO AP high make web interface
 
=head1 NAME

 ic_upload_geo_experiment.pl
 
=head1 USAGE

 ic_upload_geo_activity.pl [.csv]
 Options:
            
            'chado:s'                          => Chado connection string (e.g. dbi:Pg:dbname=...)
            'profile:s'                        => GMOD DB profile
            'quote:s'                      => \$quote_char,
            'delim:s'                          => Delimiter for text files (def tab \t)
            'detic:s'                          => GeoDETIC code (def WGS84)
            'flat|ncbi_flatfiledir:s'          => Directory with NCBI taxonomy flatfiles
            'mainfiles:s{1,}'     => Collection site 
            'geolocation_metadata_files:s{1,}' => Metadata for collection sites (e.g. location metadata) NOT USED
            'activity_metadata_files:s{1,}' => Metadata for activities (e.g. event metadata such as year)
            'genotype_data_file:s{1,}'     => Genotype data
            'phenotype_data_file:s{1,}'         => \@phenotype_data_files,
            'sequence_metadata_file:s{1,}'     => Sequence metadata
            'fasta_file:s'                    => FASTA with sequence data if sequences are not from GenBank
            'gff_file:s'                      => \$gff_file,
            'organism:s'                       => Organism name/taxID
            'phylogeny_files:s{,}'             => \@phylogeny_files,
            'collection_project_name:s'    => \$collection_project,
            'project_name:s'    => \$main_project,

            
=head1 AUTHORS

    Alexie Papanicolaou * 
    and authors of bulk_load_gff (Scott Cain et al?)
    
    CSIRO Ecosystem Sciences, GPO 1700, Canberra, Australia
    * recipient of complaints:  alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

 This software is released under the GNU General Public License version 3 (GPLv3).
 It is provided "as is" without warranty of any kind.
 You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html. 
 Please note that incorporating the whole software or parts of its code in proprietary software 
 is prohibited under the current license.

=head1 Controlled Vocabularies

CVs are derived from two sources: obofoundry (http://www.obofoundry.org/) for genetics and Darwin Core for ecological data (http://rs.tdwg.org/dwc/). 
Some custom CVs have been made for nd_experiment.type and stock.type  

  Geolocation Metadata CVs:
* Locality: http://purl.org/dc/terms/Location http://code.google.com/p/darwincore/wiki/Location
* GeologicalContext
  
  activity Metadata CVs:
* Event
    
  Stock Metadata CVs:
* Occurrence
* Identification 

  Organism Metadata CVs
* Taxon 
    
  ?? cvterm Metadata CVs:
* MeasurementOrFact

See GMOD's description of the Chado database or docs/darwin_core.cvs.txt docs/nd_types.cvs.txt for more details

=head1 Implementation of ND

=head2 Introduction

 (writing in progress)

The superior entity is I<project>. A project brings together activity I<experiments> so for easy querying. Such activitys can be collections, genotyping, phenotyping
or any kind of activity a researcher is undertaking. 

A tier 1 project houses all experiments of a particular activity type. 

A tier 2 project houses all experiments of a particular research topic (e.g. paper) regardless of experiment type.

Meta-projects can be created to link experiments of particular interest.

The project.name (unique) is our central entity when preparing and loading CSV files. 

=head3 Collection Experiments/activities

Collection experiments involve the collection of a specimen (called I<stock>) from an I<environment> such as a habitat, a laboratory setting or a field trial.

We define environment as the collection of biotic and abiotic factors that define a geolocation at a particular point in time.
As environments are not static, our implementation created a link between environment and I<geolocation> and demands a I<date>. 
Further, one geolocation can have many environments. For example, different seasons.

Thus the basic line of thought when adding data to a database is (and hence how this script works):
  PROJECT -> ENVIRONMENT -> GEOLOCATION + time -> collection EXPERIMENT -> genotyping EXPERIMENT -> STOCK

 in other words: For a project, at a particular location and time, a particular environment is investigated collections of specimens are taken and genotyped
  
The linkage of tables in ND are however:
  PROJECT -> EXPERIMENT -> GEOLOCATION -> ENVIRONMENT 
                        -> STOCK
                        
 in other words, experiment is the central entity.

=head3 Genotyping and phenotyping experiments

A collection results in a specimen (stock). That specimen (or its descendants) can get genotyped by a genotyping experiment. A feature is I<always> created.
It can also be phenotyped by a phenotyping experiment.

phenotype:
 A phenotype describes a trait, something that can be measured and quantified or classified. It has one or more genotypes, linked via
 phenstatement (also requires an environment which is poorly defined and used by flybase)
  
genotype:
 A genotype is the variant of a genomic region that influences a phenotype. It has one or more features (usually one) 
 Each such feature is linked via feature_genotype
 The genotype /available/ for a given loci is stored in genotype
 The genotype observed for an individual is stored BOTH in stock_genotype (flybase) and in nd_experiment_stock (ND)
 
 "The sample stock is linked to the observed genotype through the nd_experiment table since a genotypic assay was performed in order to determine what the observed genotype for a given sample is.
 Thus, as with phenotypic data, a single experiment in the nd_experiment table identifies the genotype of a single genetic locus in a single sample."

 Per locus: the 'feature_genotype' individual holds an entry for each allelic sequence feature (i.e. AGTCcGG and AGTCgGG for diploids) and uses the same genotype_id

=head2 Notes

On stock-to-stock transversal:
 Our implementation is KISS: it does not aspire to be a LIMS and therefore we do not track stock relationships: A collected specimen stock
is the same stock that is used in a genotyping experiment even after a series of extractions, amplification and library preparations.

Next-Gen Sequencing:
 We do not worry about database inflation when storing v. large numbers of sequence features. Chado is our store, not a warehouse that drives a software
(we will let mongodb deal with that).

=cut
use strict;
use Data::Dumper;
use DBI;
use Getopt::Long;
use Pod::Usage;
use Bio::Taxon;
use Bio::DB::Taxonomy;
use Bio::GMOD::Config;
use Bio::GMOD::DB::Config;
use Bio::Chado::Schema;
use Bio::GMOD::DB::Adapter;

#loading sequence data to chado from genbank2gff
#TODO: reduce and remove
use Bio::SeqFeature::Tools::Unflattener;
use Bio::SeqFeature::Tools::TypeMapper;
use Bio::SeqFeature::Tools::IDHandler;
use Bio::DB::GenBank;
use Bio::FeatureIO;
use URI::Escape;
use Module::Load;
my $flatfile_dir = '/databases/ncbi_taxonomy/';
my $delimiter    = "\t";
my $geodetic     = 'WGS84';
my $dbprefix     = '';
my $quote_char   = '"';

#my $dsn_drupal = 'dbi:Pg:dbname=drupal-test;host=localhost;port=5432;user=www-db;password=123';
my $dsn_chado =
'dbi:Pg:dbname=biosecurity_chado;host=localhost;port=5432;user=chado_dbadmin;password=123';

#read only my $dsn_chado = 'dbi:Pg:dbname=chado;host=localhost;port=5432;user=chado_dbuser;password=123';
my @collection_sites_files
  ; # a list of collection sites; site.uniquename,long,lat,altitude,description; &extract_store_geolocation()
my @geolocation_metadata_files
  ; # data describing collection sites as location; site.uniquename,cv.name,cvterm.name,value

#my @activity_data_files
; # description of collection activity; activity.uniquename,activity_type_name,site.uniquename &extract_store_activity()
my @activity_metadata_files
  ; # data derived from a collection activity such as dates;  activity.uniquename,cv.name,cvterm.name,value
my ( $ncbi_taxid, @genotype_data_files );
my $chado_profile = 'default';
my (
     %genbank_script_method, $genbank_script_idh,
     %genbank_script_seen,   %genbank_script_proteinfa,
     @genbank_script_return, $sequence_metadata_file,
     $chado_lock,            $chado_obj,
     @phylogeny_files,       $sequence_file,
     @phenotype_data_files,  $collection_project,
     $main_project,          $gff_file
);
GetOptions(
  'quote:s'                 => \$quote_char,
  'chado:s'                 => \$dsn_chado,
  'prefix:s'                => \$dbprefix,
  'delim:s'                 => \$delimiter,
  'detic:s'                 => \$geodetic,
  'flat|ncbi_flatfiledir:s' => \$flatfile_dir,
  'mainfiles:s{1,}'         => \@collection_sites_files,
  'activity_metadata_files:s{1,}' =>
    \@activity_metadata_files,         # metadata about the collection activity
  'geolocation_metadata_files:s{1,}' => \@geolocation_metadata_files,
  'genotype_data_file:s{1,}'         => \@genotype_data_files,
  'phenotype_data_file:s{1,}'        => \@phenotype_data_files,
  'sequence_metadata_files:s'        => \$sequence_metadata_file,
  'fasta_file:s'                     => \$sequence_file,
  'gff_file:s'                       => \$gff_file,
  'organism:s'                       => \$ncbi_taxid,
  'phylogeny_files:s{,}'             => \@phylogeny_files,
  'profile:s'                        => \$chado_profile,
  'collection_project_name:s'        => \$collection_project,
  'project_name:s'                   => \$main_project,
);

# checks
pod2usage "Need main project name\n" unless $main_project;
pod2usage "You need to provide collection activities as CSV file(s)\n"
  unless (@collection_sites_files);
pod2usage "In order to database genotypes, you need sequence CSV file(s)\n"
  unless (@genotype_data_files);
$ncbi_taxid = &get_ncbi_taxid($ncbi_taxid) if $ncbi_taxid;
my $db_conf = &prep_chado($chado_profile);
my $schema  = Bio::Chado::Schema->connect( $db_conf->dsn, $db_conf->user,
                                          $db_conf->password || "",
                                          { AutoCommit => 1 } );

# Taxonomy database
my $taxondb = &get_taxondb();
my $dbh     = &connecttodb($dsn_chado);
my ( $organism_id, $species_name ) = &add_organism();

# this will be populated later on
my ($main_project_id);
my (
     $check_project_prepared,
     $get_project_prepared,
     $insert_project_prepared,
     $get_projectrel_prepared,
     $insert_projectrel_prepared,
     $select_experiment_project_prepared,
     $insert_experiment_project_prepared,
     $select_stock_genotype_prepared,
     $insert_stock_genotype_prepared,
     $get_feat_name_prepared,
     $select_stock_prepared,
     $next_stock_prepared,
     $insert_stock_prepared,
     $select_featureloc_prepared,
     $insert_featureloc_prepared,
     $check_genotype_prepared,
     $insert_genotype_prepared,
     $check_phenotype_prepared,
     $insert_phenotype_prepared,
     $select_environment_prepared,
     $insert_environment_prepared,
     $update_environment_descr_prepared,
     $select_geolocation_id_prepared,
     $insert_geolocation_prepared,
     $update_geolocation_description_prepared,
     $select_experiment_id_prepared,
     $insert_experiment_prepared,
     $select_experiment_genotype_prepared,
     $insert_experiment_genotype_prepared,
     $select_experiment_phenotype_prepared,
     $insert_experiment_phenotype_prepared,
     $select_feature_genotype_prepared,
     $insert_feature_genotype_prepared,
     $select_cv_prepared,
     $insert_cv_prepared,
     $select_cvterm_via_dbxref_prepared,
     $select_test_cvterm_prepared,
     $select_cvterm_prepared,
     $insert_cvterm_prepared,
     $select_featureprop_rank_prepared,
     $select_featureprop_id_prepared,
     $check_featureprop_prepared,
     $insert_featureprop_prepared,
     $insert_into_db_prepared,
     $select_dbxref_check_prepared,
     $insert_into_dbxref_prepared,
     $select_feature_dbxref_prepared,
     $insert_into_feature_dbxref_prepared,
     $select_global_dbxref_prepared,
     $select_global_db_prepared,
     $cvterms_list_prepared,
     $dbxref_list_prepared,
     $get_feat_id_prepared,
     $select_dbxref_prepared_cached,
     $select_experimentprop_rank_prepared,
     $select_experimentprop_id_prepared,
     $check_experimentprop_prepared,
     $insert_experimentprop_prepared,
     $select_experiments_of_environment_prepared,
     $select_experiment_stock_prepared,
     $insert_experiment_stock_prepared
) = &prepare_sqls($dbh);

# PRELUDE
# prepare database
#load CVs
my (
     $geolocation_metadata, $activity_metadata, $stock_metadata,
     $organism_metadata,    $cvterm_metadata
) = &create_darwin_core_cvs();
my ( $experiment_cv, $stock_cv, $genotype_cv, $relationship_cv ) =
  &other_nd_cvs();
sleep(0.1);
# Ends preamble
## Go alpha leader
my ($experiments_hash_ref);    # keys are the geo-experiment u.name
## The following is not making use of the feature group of tables. it is simpler
if (@collection_sites_files) {
## The following is for data that contain sequence data, such as genotypes etc
#geo experiment data; $experiments_hash_ref will have environement_id, geolocation id and type_id as key-value pairs
  foreach my $collection_site (@collection_sites_files) {
    unless ($collection_project) {
      open( IN, $collection_site ) || die();
      my $header = <IN>;
      close IN;
      $header =~ /^#([^\s\|]+)\|/;
      $collection_project = $1;
    }
    print "Processing experiment $collection_project\n";
    ( $main_project_id, my $collection_project_id, my $projectrel_id ) =
      &add_projects( $main_project, $collection_project );
    $experiments_hash_ref =
      &extract_store_geolocation( $collection_site, $collection_project,
                                  $collection_project_id );
  }

#create/link experiment with geolocation; adds $experiments_hash_ref value 'experiment_id'
  &add_geo_experiment();
  sleep(0.1);
  &link_projects();
  sleep(0.1);

  # add stocks using $experiments_hash_ref
  &add_stocks();
  sleep(0.1);
  &add_features();
  sleep(0.1);

#feature data! let's make some genotypes... ($experiments_hash_ref key 'genotypes')
  &process_genotypes($collection_project) if @genotype_data_files;
  sleep(0.1);
  &process_phenotypes($collection_project) if @phenotype_data_files;
  sleep(0.1);
}

#for debug  warn Dumper $experiments_hash_ref;

# add metadata
&add_activity_metadata() if scalar(@activity_metadata_files) > 0;
sleep(0.1);
&add_geolocation_metadata() if scalar(@geolocation_metadata_files) > 0;
sleep(0.1);
&disconnectdb($dbh);
########################################################
sub add_organism() {
  my $taxon = $taxondb->get_taxon( -taxonid => $ncbi_taxid )
    || die "Cannot find your query\n";
  my $species_name = $taxon->scientific_name();
  $species_name =~ /^(\S+)\s+(\S+)/;
  my $genus   = $1;
  my $species = $2;
  my $comment = "";
  my $common_name =
    $taxon->common_names() ? $taxon->common_names() : $genus . ' ' . $species
    if $genus ne 'Unknown';
  my $abbr = substr( $genus, 0, 1 ) . '.' . $species;
  warn "Adding $genus $species ($abbr)\n";
  my $result = $schema->resultset("Organism::Organism")->find_or_new(
                                                         {
                                                           genus   => $genus,
                                                           species => $species,
                                                           abbreviation => $abbr
                                                         }
  );
  $result->insert;
  my $organism_id = $result->organism_id();
  return ( $organism_id, $abbr );
}

sub add_features() {
  if ($gff_file) {
    system(
"gmod_bulk_load_gff3.pl -org $species_name -dbprofile $chado_profile -gff $gff_file"
    );
  }
  if ($sequence_file) {
    system(
"gmod_bulk_load_gff3.pl -org $species_name -dbprofile $chado_profile --fastafile $sequence_file"
    );
  }
}

sub add_projects() {

  # add/get main project
  # add/get collection project
  my $main_project       = shift;
  my $collection_project = shift;
  my ($main_project_id);
  my $partof_cvterm = &get_cvterm( 'relationship', 'part_of' );
  if ( $main_project =~ /^\d+$/ ) {
    $main_project_id = $main_project;
    $check_project_prepared->execute($main_project_id);
    ($main_project_id) = $check_project_prepared->fetchrow_array();
    if ( !$main_project_id ) {
      die "ERROR: Cannot find project in database.\n";
    }
  } else {
    $get_project_prepared->execute($main_project);
    ($main_project_id) = $get_project_prepared->fetchrow_array();
    if ( !$main_project_id ) {
      $insert_project_prepared->execute($main_project);
      $get_project_prepared->execute($main_project);
      ($main_project_id) = $get_project_prepared->fetchrow_array();
      if ( !$main_project_id ) {
        die "ERROR: Cannot add project to database. Check for permissions\n";
      }
    }
  }
  $get_project_prepared->execute($collection_project);
  my ($collection_project_id) = $get_project_prepared->fetchrow_array();
  if ( !$collection_project_id ) {
    $insert_project_prepared->execute($collection_project);
    $get_project_prepared->execute($collection_project);
    ($collection_project_id) = $get_project_prepared->fetchrow_array();
    if ( !$collection_project_id ) {
      die "ERROR: Cannot add project to database. Check for permissions\n";
    }
  }

  #link the two
  $get_projectrel_prepared->execute( $main_project_id, $collection_project_id,
                                     $partof_cvterm );
  my ($projectrel_id) = $get_projectrel_prepared->fetchrow_array();
  if ( !$projectrel_id ) {
    $insert_projectrel_prepared->execute( $main_project_id,
                                       $collection_project_id, $partof_cvterm );
    $get_projectrel_prepared->execute( $main_project_id, $collection_project_id,
                                       $partof_cvterm );
    ($projectrel_id) = $get_projectrel_prepared->fetchrow_array();
    if ( !$projectrel_id ) {
      die "ERROR: Cannot link projects in database. Check for permissions\n";
    }
  }
  return ( $main_project_id, $collection_project_id, $projectrel_id );
}

sub link_projects() {

  # with nd_experiment_project,
  # $collection_project_id
  #$select_experiment_project_prepared ,$insert_experiment_project_prepared
  foreach my $experiment ( keys %{$experiments_hash_ref} ) {
    my $experiment_id = $experiments_hash_ref->{$experiment}->{'experiment_id'};
    my $collection_project_id =
      $experiments_hash_ref->{$experiment}->{'collection_project_id'};
    $select_experiment_project_prepared->execute( $collection_project_id,
                                                  $experiment_id );
    my ($id) = $select_experiment_project_prepared->fetchrow_array();
    if ( !$id ) {
      $insert_experiment_project_prepared->execute( $collection_project_id,
                                                    $experiment_id );
      $select_experiment_project_prepared->execute( $collection_project_id,
                                                    $experiment_id );
      ($id) = $select_experiment_project_prepared->fetchrow_array();
      if ( !$id ) {
        die
"ERROR: Failed to link project ($collection_project_id) with experiment ($experiment_id). Check for permissions\n";
      }
    }
  }
}

sub prep_chado() {
  $chado_profile = shift;
  my $gmod_conf =
    ( $ENV{'GMOD_ROOT'} || "/usr/local/share/gmod" )
    ? Bio::GMOD::Config->new( $ENV{'GMOD_ROOT'} || "/usr/local/share/gmod" )
    : Bio::GMOD::Config->new();
  system(
"gmod_bulk_load_gff3.pl --recreate_cache -org fromdata -dbprofile $chado_profile </dev/null"
  );
  return Bio::GMOD::DB::Config->new( $gmod_conf, $chado_profile );
}

sub prepare_sqls() {
  my $dbh = shift;

  # GLOBAL SQLs
  # DBXREF
  my $select_global_db_sql = "SELECT db_id,name from db where name=?";
  my $select_global_dbxref_sql =
"SELECT dbxref_id from dbxref where accession=? AND db_id = (SELECT db_id from db where name =? ) ";
  my $dbxref_list_sql =
"SELECT dbxref_id from dbxref where accession=? AND db_id = (SELECT db_id from db where name =? )";
  my $insert_into_db_sql = 'INSERT into db (name) VALUES (?)';
  my $insert_into_dbxref_sql =
    'INSERT into dbxref (accession,db_id,description) VALUES (?,?,?)';

  #CVTERMs
  my $select_cvterm_sql =
"SELECT cvterm_id from cvterm where name=? AND cv_id=(select cv_id from cv where name=?) AND is_obsolete=0 limit 1";
  my $select_test_cvterm_sql =
"SELECT cvterm_id from cvterm where cv_id=(select cv_id from cv where name=?) AND name=? AND definition=? AND dbxref_id=(select dbxref_id from dbxref where db_id=(select db_id from db where name='null') AND accession=?)";
  my $select_cvterm_via_dbxref_sql =
    "SELECT cvterm_id from cvterm where dbxref_id=?";
  my $insert_cvterm_sql =
"INSERT into cvterm (cv_id,name,definition,dbxref_id) VALUES ((select cv_id from cv where name=?),?,?,?)";
  my $cvterms_list_sql =
"SELECT cvterm_id from cvterm where cv_id IN (select cv_id from cv where name = ? ) and dbxref_id IN (select dbxref_id from dbxref where accession=? and db_id =(select db_id from db where name=?))";
  my $select_cv_sql = "SELECT cv_id from cv where name = ?";
  my $insert_cv_sql = "INSERT INTO cv (name) VALUES (?)";

  #feature
  my $get_feat_id_sql   = 'SELECT feature_id from feature where uniquename=?';
  my $get_feat_name_sql = 'SELECT uniquename from feature where feature_id=?';
  my $select_featureprop_rank_sql =
    'SELECT max(rank) from featureprop where feature_id=? AND type_id=?';
  my $select_featureprop_id_sql =
'SELECT featureprop_id from featureprop where feature_id=? AND type_id=? AND value=?';
  my $check_featureprop_sql =
'SELECT featureprop_id from featureprop where feature_id=? AND type_id=? AND rank=?';
  my $insert_featureprop_sql =
    'INSERT into featureprop (type_id,feature_id,value,rank) VALUES (?,?,?,?)';
  my $select_dbxref_sql =
    'SELECT dbxref_id from dbxref where db_id=? AND accession=?';
  my $select_feature_dbxref_sql =
    'SELECT dbxref_id from feature_dbxref where feature_id=? AND dbxref_id=?';
  my $insert_into_feature_dbxref_sql =
    'INSERT into feature_dbxref (feature_id,dbxref_id) VALUES (?,?)';

  #featureloc
  my $select_featureloc_sql =
"SELECT featureloc_id from featureloc where srcfeature_id = (SELECT feature_id from feature where uniquename=?) and feature_id= (SELECT feature_id from feature where uniquename=?)";
  my $insert_featureloc_sql =
"INSERT INTO featureloc (srcfeature_id,feature_id,fmin,fmax,strand,rank) VALUES ((SELECT feature_id from feature where uniquename=?,(SELECT feature_id from feature where uniquename=?),"
    . ",?,?,?,(SELECT max(rank)+1 from featureloc where srcfeature_id=? AND feature_id=? AND locgroup=0) )";

  # genotype
  my $check_genotype_sql =
    "SELECT genotype_id from genotype where uniquename=?";
  my $insert_genotype_sql =
    "INSERT INTO genotype (uniquename,type_id) VALUES (?,?)";

  # phenotype
  my $check_phenotype_sql =
    "SELECT phenotype_id from phenotype where uniquename=?";
  my $insert_phenotype_sql =
    "INSERT INTO phenotype (uniquename,observable_id,value) VALUES (?,?,?)";

  #environment
  my $select_environment_sql =
    "SELECT environment_id from environment WHERE uniquename=?";
  my $insert_environment_sql =
    "INSERT INTO environment (uniquename) VALUES (?)";
  my $update_environment_descr_sql =
    "UPDATE environment SET description=? where uniquename=?";

  #geolocation
  my $select_geolocation_id_sql =
"SELECT max(nd_geolocation_id) from nd_geolocation where latitude=? AND longitude=? AND geodetic_datum=?";
  my $insert_geolocation_sql =
"INSERT INTO nd_geolocation (latitude,longitude,geodetic_datum,altitude) VALUES (?,?,?,?)";
  my $update_geolocation_description_sql =
"UPDATE nd_geolocation set description=? where nd_geolocation_id=? AND description IS NOT NULL";

  #experiment
  my $select_experiment_id_sql =
"SELECT nd_experiment_id from nd_experiment where nd_geolocation_id=? and type_id=?";
  my $insert_experiment_sql =
    "INSERT INTO nd_experiment (nd_geolocation_id,type_id) VALUES (?,?)";

  #experiment_genotype
  my $select_experiment_genotype_sql =
"SELECT nd_experiment_genotype_id from nd_experiment_genotype where nd_experiment_id=? and genotype_id=?";
  my $insert_experiment_genotype_sql =
"INSERT INTO nd_experiment_genotype (nd_experiment_id,genotype_id) VALUES (?,?)";

  #experiment_phenotype
  my $select_experiment_phenotype_sql =
"SELECT nd_experiment_phenotype_id from nd_experiment_phenotype where nd_experiment_id=? and phenotype_id=?";
  my $insert_experiment_phenotype_sql =
"INSERT INTO nd_experiment_phenotype (nd_experiment_id,phenotype_id) VALUES (?,?)";

  #experimentprop
  my $select_experimentprop_rank_sql =
'SELECT max(rank) from nd_experimentprop where nd_experiment_id=? AND type_id=?';
  my $select_experimentprop_id_sql =
'SELECT nd_experimentprop_id from nd_experimentprop where nd_experiment_id=? AND type_id=? AND value=?';
  my $check_experimentprop_sql =
'SELECT nd_experimentprop_id from nd_experimentprop where nd_experiment_id=? AND type_id=? AND rank=?';
  my $insert_experimentprop_sql =
'INSERT into nd_experimentprop (type_id,nd_experiment_id,value,rank) VALUES (?,?,?,?)';

  #experiments of environment
  my $select_experiments_of_environment_sql =
      "SELECT nd_experiment_id as id from nd_experiment WHERE "
    . "nd_geolocation_id IN (SELECT nd_geolocation_id FROM nd_geolocation_environment WHERE "
    . "environment_id=(SELECT environment_id from environment WHERE uniquename=? ))";
  my $select_feature_genotype_sql =
"SELECT feature_genotype_id from feature_genotype where feature_id=? AND genotype_id=? AND rank=0 "
    . "AND cgroup=0 AND cvterm_id=(SELECT cvterm_id from cvterm where name='null' and cv_id =(SELECT cv_id from cv where name='null'))";
  my $insert_feature_genotype_sql =
"INSERT INTO feature_genotype (feature_id,genotype_id,rank,cgroup,cvterm_id) VALUES (?,?,0,0,"
    . "(SELECT cvterm_id from cvterm where name='null' and cv_id =(SELECT cv_id from cv where name='null')))";
### Stocks; organism is derived from the feature; relationship cvterm used is 'produces'
  my $select_stock_sql =
"SELECT stock_id from stock where organism_id=$organism_id AND type_id=? AND uniquename=? ";
  my $insert_stock_sql =
"INSERT INTO stock (organism_id,type_id,uniquename,name) VALUES ($organism_id,?,?,?)";
  my $next_stock_sql = "SELECT nextval('stock_stock_id_seq')";
### stock_genotpye
  my $select_stock_genotype_sql =
"SELECT stock_genotype_id from stock_genotype where stock_id=? AND genotype_id=?";
  my $insert_stock_genotype_sql =
    "INSERT INTO stock_genotype (stock_id,genotype_id) VALUES (?,?)";
  my $select_experiment_stock_sql =
"SELECT nd_experiment_stock_id from nd_experiment_stock where nd_experiment_id=? AND stock_id=? AND type_id="
    . "(SELECT cvterm_id from cvterm where name='produces' and cv_id=(SELECT cv_id from cv where name='relationship') )";
  my $insert_experiment_stock_sql =
"INSERT INTO nd_experiment_stock (nd_experiment_id,stock_id,type_id) VALUES (?,?,"
    . "(SELECT cvterm_id from cvterm where name='produces' and cv_id=(SELECT cv_id from cv where name='relationship') ) )";
  ### project
  my $check_project_sql = "SELECT project_id from project WHERE project_id=?";
  my $get_project_sql   = "SELECT project_id from project WHERE name=?";
  my $insert_project_sql =
    "INSERT INTO project (name,description) VALUES (?,'None')";
  my $get_projectrel_sql =
"SELECT project_relationship_id from project_relationship WHERE subject_project_id=? AND object_project_id=? AND type_id=?";
  my $insert_projectrel_sql =
"INSERT INTO project_relationship (subject_project_id,object_project_id,type_id) VALUES (?,?,?)";
  ##experiment_project
  my $select_experiment_project_sql =
"SELECT nd_experiment_project_id as id from nd_experiment_project where project_id=? AND nd_experiment_id=?";
  my $insert_experiment_project_sql =
"INSERT INTO nd_experiment_project (project_id,nd_experiment_id) VALUES (?,?)";
#####
##    prepare statements (some selects are cached)
#####
  my $check_project_prepared     = $dbh->prepare($check_project_sql);
  my $get_project_prepared       = $dbh->prepare($get_project_sql);
  my $insert_project_prepared    = $dbh->prepare($insert_project_sql);
  my $get_projectrel_prepared    = $dbh->prepare($get_projectrel_sql);
  my $insert_projectrel_prepared = $dbh->prepare($insert_projectrel_sql);
  my $select_experiment_project_prepared =
    $dbh->prepare($select_experiment_project_sql);
  my $insert_experiment_project_prepared =
    $dbh->prepare($insert_experiment_project_sql);
  my $select_stock_genotype_prepared =
    $dbh->prepare($select_stock_genotype_sql);
  my $insert_stock_genotype_prepared =
    $dbh->prepare($insert_stock_genotype_sql);
  my $get_feat_name_prepared      = $dbh->prepare($get_feat_name_sql);
  my $select_stock_prepared       = $dbh->prepare($select_stock_sql);
  my $next_stock_prepared         = $dbh->prepare($next_stock_sql);
  my $insert_stock_prepared       = $dbh->prepare($insert_stock_sql);
  my $select_featureloc_prepared  = $dbh->prepare($select_featureloc_sql);
  my $insert_featureloc_prepared  = $dbh->prepare($insert_featureloc_sql);
  my $check_genotype_prepared     = $dbh->prepare($check_genotype_sql);
  my $insert_genotype_prepared    = $dbh->prepare($insert_genotype_sql);
  my $check_phenotype_prepared    = $dbh->prepare($check_phenotype_sql);
  my $insert_phenotype_prepared   = $dbh->prepare($insert_phenotype_sql);
  my $select_environment_prepared = $dbh->prepare($select_environment_sql);
  my $insert_environment_prepared = $dbh->prepare($insert_environment_sql);
  my $update_environment_descr_prepared =
    $dbh->prepare($update_environment_descr_sql);
  my $select_geolocation_id_prepared =
    $dbh->prepare($select_geolocation_id_sql);
  my $insert_geolocation_prepared = $dbh->prepare($insert_geolocation_sql);
  my $update_geolocation_description_prepared =
    $dbh->prepare($update_geolocation_description_sql);
  my $select_experiment_id_prepared =
    $dbh->prepare($select_experiment_id_sql);    # no need for cache
  my $insert_experiment_prepared = $dbh->prepare($insert_experiment_sql);
  my $select_experiment_genotype_prepared =
    $dbh->prepare($select_experiment_genotype_sql);
  my $insert_experiment_genotype_prepared =
    $dbh->prepare($insert_experiment_genotype_sql);
  my $select_experiment_phenotype_prepared =
    $dbh->prepare($select_experiment_phenotype_sql);
  my $insert_experiment_phenotype_prepared =
    $dbh->prepare($insert_experiment_phenotype_sql);
  my $select_feature_genotype_prepared =
    $dbh->prepare($select_feature_genotype_sql);
  my $insert_feature_genotype_prepared =
    $dbh->prepare($insert_feature_genotype_sql);
  my $select_cv_prepared = $dbh->prepare_cached($select_cv_sql);
  my $insert_cv_prepared = $dbh->prepare_cached($insert_cv_sql);
  my $select_cvterm_via_dbxref_prepared =
    $dbh->prepare_cached($select_cvterm_via_dbxref_sql);
  my $select_test_cvterm_prepared = $dbh->prepare($select_test_cvterm_sql);
  my $select_cvterm_prepared      = $dbh->prepare($select_cvterm_sql);
  my $insert_cvterm_prepared      = $dbh->prepare($insert_cvterm_sql);
  my $select_featureprop_rank_prepared =
    $dbh->prepare($select_featureprop_rank_sql);
  my $select_featureprop_id_prepared =
    $dbh->prepare($select_featureprop_id_sql);
  my $check_featureprop_prepared   = $dbh->prepare($check_featureprop_sql);
  my $insert_featureprop_prepared  = $dbh->prepare($insert_featureprop_sql);
  my $insert_into_db_prepared      = $dbh->prepare($insert_into_db_sql);
  my $select_dbxref_check_prepared = $dbh->prepare($select_dbxref_sql);
  my $insert_into_dbxref_prepared  = $dbh->prepare($insert_into_dbxref_sql);
  my $select_feature_dbxref_prepared =
    $dbh->prepare($select_feature_dbxref_sql);
  my $insert_into_feature_dbxref_prepared =
    $dbh->prepare($insert_into_feature_dbxref_sql);
  my $select_global_dbxref_prepared =
    $dbh->prepare_cached($select_global_dbxref_sql);
  my $select_global_db_prepared = $dbh->prepare_cached($select_global_db_sql);
  my $cvterms_list_prepared     = $dbh->prepare_cached($cvterms_list_sql);
  my $dbxref_list_prepared      = $dbh->prepare($dbxref_list_sql);
  my $get_feat_id_prepared      = $dbh->prepare_cached($get_feat_id_sql);
  my $select_dbxref_prepared_cached = $dbh->prepare_cached($select_dbxref_sql);
  my $select_experimentprop_rank_prepared =
    $dbh->prepare($select_experimentprop_rank_sql);
  my $select_experimentprop_id_prepared =
    $dbh->prepare($select_experimentprop_id_sql);
  my $check_experimentprop_prepared = $dbh->prepare($check_experimentprop_sql);
  my $insert_experimentprop_prepared =
    $dbh->prepare($insert_experimentprop_sql);
  my $select_experiments_of_environment_prepared =
    $dbh->prepare_cached($select_experiments_of_environment_sql);
  my $select_experiment_stock_prepared =
    $dbh->prepare($select_experiment_stock_sql);
  my $insert_experiment_stock_prepared =
    $dbh->prepare($insert_experiment_stock_sql);
  return (
           $check_project_prepared,
           $get_project_prepared,
           $insert_project_prepared,
           $get_projectrel_prepared,
           $insert_projectrel_prepared,
           $select_experiment_project_prepared,
           $insert_experiment_project_prepared,
           $select_stock_genotype_prepared,
           $insert_stock_genotype_prepared,
           $get_feat_name_prepared,
           $select_stock_prepared,
           $next_stock_prepared,
           $insert_stock_prepared,
           $select_featureloc_prepared,
           $insert_featureloc_prepared,
           $check_genotype_prepared,
           $insert_genotype_prepared,
           $check_phenotype_prepared,
           $insert_phenotype_prepared,
           $select_environment_prepared,
           $insert_environment_prepared,
           $update_environment_descr_prepared,
           $select_geolocation_id_prepared,
           $insert_geolocation_prepared,
           $update_geolocation_description_prepared,
           $select_experiment_id_prepared,
           $insert_experiment_prepared,
           $select_experiment_genotype_prepared,
           $insert_experiment_genotype_prepared,
           $select_experiment_phenotype_prepared,
           $insert_experiment_phenotype_prepared,
           $select_feature_genotype_prepared,
           $insert_feature_genotype_prepared,
           $select_cv_prepared,
           $insert_cv_prepared,
           $select_cvterm_via_dbxref_prepared,
           $select_test_cvterm_prepared,
           $select_cvterm_prepared,
           $insert_cvterm_prepared,
           $select_featureprop_rank_prepared,
           $select_featureprop_id_prepared,
           $check_featureprop_prepared,
           $insert_featureprop_prepared,
           $insert_into_db_prepared,
           $select_dbxref_check_prepared,
           $insert_into_dbxref_prepared,
           $select_feature_dbxref_prepared,
           $insert_into_feature_dbxref_prepared,
           $select_global_dbxref_prepared,
           $select_global_db_prepared,
           $cvterms_list_prepared,
           $dbxref_list_prepared,
           $get_feat_id_prepared,
           $select_dbxref_prepared_cached,
           $select_experimentprop_rank_prepared,
           $select_experimentprop_id_prepared,
           $check_experimentprop_prepared,
           $insert_experimentprop_prepared,
           $select_experiments_of_environment_prepared,
           $select_experiment_stock_prepared,
           $insert_experiment_stock_prepared
  );
}

sub get_taxondb() {
  my $taxondb;
  Bio::DB::Taxonomy->new( -source => 'entrez' );
  if (    $flatfile_dir
       && -d $flatfile_dir
       && -s $flatfile_dir . '/nodes.dmp'
       && -s $flatfile_dir . '/names.dmp' )
  {
    $taxondb = Bio::DB::Taxonomy->new(
                                     -source    => 'flatfile',
                                     -nodesfile => $flatfile_dir . '/nodes.dmp',
                                     -namesfile => $flatfile_dir . '/names.dmp',
                                     -directory => $flatfile_dir
    );
  }
  return $taxondb;
}

sub lat2dec() {
  my $n = shift;
  die "Latitude column does not exist\n" unless $n;
  return $n if $n =~ /^[+\-]?\d+\.?\d*$/;
  my $dec;
  my @points = split( '/', $n );
  die "Latitude column does not seem to be in degrees\n" unless $points[2];
  $dec = sprintf( "%.4f",
                  ( $points[0] + ( $points[1] / 60 ) + ( $points[2] / 3600 ) )
  );

  #TODO australia
  $dec *= -1;
  return $dec;
}

sub long2dec() {
  my $n = shift;
  die "Longitude column does not exist\n" unless $n;
  return $n if $n =~ /^[+\-]?\d+\.?\d*$/;
  my $dec;
  my @points = split( '/', $n );
  die "Longitude column does not seem to be in degrees\n" unless $points[2];
  $dec = sprintf( "%.4f",
                  ( $points[0] + ( $points[1] / 60 ) + ( $points[2] / 3600 ) )
  );

  #TODO australia
  return $dec;
}

=pod

=head1 Term definitions within (our) ND

TODO cleanup this


environment:
 Luckily we can define the environment very well here, it is taken from the geolocation but there is no usage/link within the ND.
So we can make use of it within our implementation of the ND by using it as - what I call - project specific geo-experiment uniquename
This works rather well because environment only requires a name and we had nowhere else to store it.

I've proposed to the group that we create a environment_geolocation with (environment_id FK,geolocation_id FK, date) + PK to host the snapshot of the geolocation's
environment (incl. biotic and abiotic factors). We could store the layers here as properties (we will need to create an environmentprop)

phendesc and phenstatement:
We can use this for population genetic/phylogenetic statements later on.  

Conclusion:
 For viruses, we will create a pathovar 'phenotype' which will be equal to one 'genotype' (a haplotype). This genotype has one feature. 
For diploids, at 1 particular locus which controls 1 trait, we have 1 phenotype for 1 genotype for 2 alleles (features e.g. 'sequence_variant'; etc for n-ploids)
  
also @see http://gmod.oicr.on.ca/wiki/Talk:Chado_Natural_Diversity_Module_Working_Group#Genotype

=cut 

=pod

=head1 ROUTINES

=head2 add_stocks

 TODO write details here.

=cut

sub add_stocks() {
  foreach my $experiment ( keys %{$experiments_hash_ref} ) {
    my $experiment_id = $experiments_hash_ref->{$experiment}->{'experiment_id'};
    my $stock_type_id = $experiments_hash_ref->{$experiment}->{'stock_type_id'};
    my $stock_type_name =
      $experiments_hash_ref->{$experiment}->{'stock_type_name'};
    die "Cannot find stock type/name\n"
      unless $stock_type_id && $stock_type_name;
    $next_stock_prepared->execute();
    my ($next_stock) = $next_stock_prepared->fetchrow_array();
    my $stock_uname =
      'urn:lsid:biosecurity.cbf.csiro.au:geospatial_genetics_chado_stock:'
      . $next_stock;
    my $stock_name = $experiment;
    $stock_name =~ s/^\S+:://;
    $select_stock_prepared->execute( $stock_type_id, $stock_uname );
    my ($stock_id) = $select_stock_prepared->fetchrow_array();

    if ( !$stock_id ) {
      $insert_stock_prepared->execute( $stock_type_id, $stock_uname,
                                       $stock_name );
      $select_stock_prepared->execute( $stock_type_id, $stock_uname );
      ($stock_id) = $select_stock_prepared->fetchrow_array();
    }
    if ( !$stock_id ) {
      die
"ERROR: Cannot add to database: stock with uniquename=$stock_uname and type=$stock_type_name for experiment $experiment. Are you sure you have write permissions?\n";
    }
    $experiments_hash_ref->{$experiment}->{'stock_id'} = $stock_id;

    # link with nd_experiment_stock
    $select_experiment_stock_prepared->execute( $experiment_id, $stock_id );
    my ($check) = $select_experiment_stock_prepared->fetchrow_array();
    if ( !$check ) {
      $insert_experiment_stock_prepared->execute( $experiment_id, $stock_id );
      $select_experiment_stock_prepared->execute( $experiment_id, $stock_id );
      ($check) = $select_experiment_stock_prepared->fetchrow_array();
    }
    if ( !$check ) {
      die
"ERROR: Cannot add to database: link experiment with stock with stock_id=$stock_id for experiment $experiment. Are you sure you have write permissions?\n";
    }
  }
}

sub link_stock_genotypes() {
  foreach my $experiment ( keys %{$experiments_hash_ref} ) {
    my $stock_id = $experiments_hash_ref->{$experiment}->{'stock_id'};
    foreach my $feature (
                  keys %{ $experiments_hash_ref->{$experiment}->{'genotypes'} } )
    {
      my $genotypes_ref_array =
        $experiments_hash_ref->{$experiment}->{'genotypes'}->{$feature}
        ->{'genotype_ids'}
        if $experiments_hash_ref->{$experiment}->{'genotypes'}->{$feature}
          ->{'genotype_ids'};
      foreach my $genotype_id (@$genotypes_ref_array) {

        #link with stock_genotype
        $select_stock_genotype_prepared->execute( $stock_id, $genotype_id );
        my ($check) = $select_stock_genotype_prepared->fetchrow_array();
        if ( !$check ) {
          $insert_stock_genotype_prepared->execute( $stock_id, $genotype_id );
          $select_stock_genotype_prepared->execute( $stock_id, $genotype_id );
          ($check) = $select_stock_genotype_prepared->fetchrow_array();
        }
        if ( !$check ) {
          die "ERROR: Cannot add to database: link genotype with stock with stock_id=$stock_id and genotype_id=$genotype_id for experiment $experiment. Are you sure you have write permissions?\n";
        }
      }
    }
  }
}

=pod

=head2 METADATA

 TODO write details here.

=cut

sub add_geolocation_metadata() {

  # these are for environment_cvterm
  foreach my $file (@activity_metadata_files) {
    open( IN, $file );
    while ( my $ln = <IN> ) {
      my @data = split( "\t", $ln );
      next unless $data[3];
      $select_experiments_of_environment_prepared->execute( $data[0] );
      while ( my $experiment_id =
              $select_experiments_of_environment_prepared->fetchrow_array() ) 
      {
      }
    }
  }
}

sub add_activity_metadata() {

  # these are experimentprop
}

=pod 

=head2 geo_experiments

 TODO WRITE details here
 
=cut

sub add_geo_experiment() {
  foreach my $experiment ( keys %{$experiments_hash_ref} ) {
    $select_experiment_id_prepared->execute(
                       $experiments_hash_ref->{$experiment}->{'geolocation_id'},
                       $experiments_hash_ref->{$experiment}->{'type_id'} );
    my $experiment_id = $select_experiment_id_prepared->fetchrow_array();
    if ( !$experiment_id ) {
      $insert_experiment_prepared->execute(
                       $experiments_hash_ref->{$experiment}->{'geolocation_id'},
                       $experiments_hash_ref->{$experiment}->{'type_id'} );
      $select_experiment_id_prepared->execute(
                       $experiments_hash_ref->{$experiment}->{'geolocation_id'},
                       $experiments_hash_ref->{$experiment}->{'type_id'} );
      $experiment_id = $select_experiment_id_prepared->fetchrow_array();
    }
    if ( !$experiment_id ) {
      die "ERROR: Cannot add new experiment with nd_geolocation_id="
        . $experiments_hash_ref->{$experiment}->{'geolocation_id'}
        . " and type_id="
        . $experiments_hash_ref->{$experiment}->{'type_id'}
        . " to database. Are you sure you have write permissions?\n";
    }
    $experiments_hash_ref->{$experiment}->{'experiment_id'} = $experiment_id;
  }
}

=pod

=pod 

=head2 feature_genotype

links to feature table. this is quite complicated 
 the complication is mostly unnecessary and even flybase has not been using it properly (http://gmod.oicr.on.ca/wiki/Chado_Phenotype_Module_at_FlyBase#feature_genotype) 

so: 
 let's start with 1 genotype per feature (aka haplotype)
 we assign cvterm as null.null since it is undocumented (and unnecessary)
 cgroup is set to 0 since we are not using anchors to any groups of features (featureloc does this well enough; it is unnecessary)
 that way, we have ONE feature = ONE genotype which works for genotype.type haplotype
 
if we have >1 haplotype (e.g. SNP which doesn't have it's own feature) then we have to set rank to something or simply implement a feature
for each marker genotype (see previous rant). We can't really not implement a feature per haplotype because otherwise we will not know what the sequence
of that genotype.

The conclusion is always set one feature (e.g. of type sequence_variant if we have >1) per haplotype and rank=0.
We need at least one genotype entry
For diploids and recessive markers, we have 2 haplotypes, 2 features, 1 genotype, 2 feature_genotype
For the user: they can represent this by having multiple lines in the feature file(s) but then we have to refactor experiment.uname (for it will no longer be unique)

=cut

=head2 process_genotypes

Expect a -delimiter separated file with mandatory four columns: experiment.uname, stock.type (the type of DNA used in the genotyping experiment) ,the GenBank accession (or other global accession) and the type of genotype.
the other columns can be use for featureloc (in this order): Reference SeqMin  SeqMax

if SeqMax < SeqMin then it is the antisense strand, otherwise it is the sense strand

=cut

sub process_phenotypes() {
  my $collection_project = shift;
  my %data;
  foreach my $file (@phenotype_data_files) {

    #testingAnders|stock.uname      phenotype.uname allele
    open( IN, $file );
    while ( my $ln = <IN> ) {
      next if $ln =~ /^#/ || $ln =~ /^\s*$/;
      chomp($ln);
      my @ln_data = split( $delimiter, $ln );
      my $experiment_name =
        $collection_project . '::' . $ln_data[0];    # also the stock name
      my $cv                             = $ln_data[1];
      my $marker                         = $ln_data[2];
      my $value                          = $ln_data[3];
      my $phenotype_observable_cvterm_id = get_set_cvterm( $cv, $marker );
      die "No cvterm found for CV $cv term $marker\n"
        unless $phenotype_observable_cvterm_id;

      # add to phenotype, experiment_phenotype
      $data{'experiment'}{$experiment_name}{$marker}{'observable_id'} =
        $phenotype_observable_cvterm_id;
      push(
            @{ $data{'experiment'}{$experiment_name}{$marker}{'alleles'} },
            $ln_data[1] . '::' . $ln_data[2] . '::' . $ln_data[3]
      );
      push(
            @{ $data{'experiment'}{$experiment_name}{$marker}{'values'} },
            $ln_data[3]
      );
    }
    close IN;
  }

  # add data to database
  foreach my $experiment ( keys %{ $data{'experiment'} } ) {
    my $experiment_id = $experiments_hash_ref->{$experiment}->{'experiment_id'};
    if ( !$experiment_id ) {
      warn Dumper $experiments_hash_ref->{$experiment};
      warn
"Experiment ID for $experiment was not found. You are providing a new experiment that we have no geolocation for! Skipping....\n";
      next;
    }
    foreach my $marker ( keys %{ $data{'experiment'}{$experiment} } ) {
      my @alleles = @{ $data{'experiment'}{$experiment}{$marker}{'alleles'} };
      my @values  = @{ $data{'experiment'}{$experiment}{$marker}{'values'} };
      for ( my $i = 0 ; $i < scalar(@alleles) ; $i++ ) {
        my $allele = $alleles[$i];
        my $value  = $values[$i];
        $check_phenotype_prepared->execute($allele);
        my $phenotype_id = $check_phenotype_prepared->fetchrow_array();
        if ( !$phenotype_id ) {

       # uniquename of phenotype is now equal to feature::allele name...
       #"INSERT INTO phenotype (uniquename,observable_id,value) VALUES (?,?,?)";
          $insert_phenotype_prepared->execute( $allele,
                     $data{'experiment'}{$experiment}{$marker}{'observable_id'},
                     $value );
          $check_phenotype_prepared->execute($allele);
          $phenotype_id = $check_phenotype_prepared->fetchrow_array();
        }
        if ( !$phenotype_id ) {
          die
"ERROR: Cannot add to database: phenotype with uniquename=$allele and observable_id="
            . $data{'experiment'}{$experiment}{$marker}{'observable_id'}
            . " for experiment $experiment. Are you sure you have write permissions?\n";
        }

#link to experiment via experiment_phenotype with experiment.uniquename %{$data{'experiment'}}
        $select_experiment_phenotype_prepared->execute( $experiment_id,
                                                        $phenotype_id );
        my ($link_id) = $select_experiment_phenotype_prepared->fetchrow_array();
        if ( !$link_id ) {
          $insert_experiment_phenotype_prepared->execute( $experiment_id,
                                                          $phenotype_id );
          $select_experiment_phenotype_prepared->execute( $experiment_id,
                                                         $phenotype_id );
          ($link_id) = $select_experiment_phenotype_prepared->fetchrow_array();
          if ( !$link_id ) {
            die "ERROR: Cannot add to database: experiment_phenotype with phenotype_id= $phenotype_id and experiment= $experiment_id "
              . "for experiment $experiment. Are you sure you have write permissions?\n";
          }
        }
                # here we are populating the array with the genotype info
        push(
              @{
                $experiments_hash_ref->{$experiment}->{'phenotypes'}->{$marker}->{'phenotype_ids'}
                },
              $phenotype_id
        );
      }
    }
  }
}

sub process_genotypes() {
  my $collection_project = shift;
  my %data;
  my @genbank_to_add;
  foreach my $file (@genotype_data_files) {
    open( IN, $file );
    while ( my $ln = <IN> ) {
      next if $ln =~ /^#/ || $ln =~ /^\s*$/;
      chomp($ln);

      #stock.uname      marker.type     marker  allele
      my @ln_data          = split( $delimiter, $ln );
      my $experiment_name  = $collection_project . '::' . $ln_data[0];
      my $genotype_type_id = $genotype_cv->{'genotype_type'}->{ $ln_data[1] };
      die "ERROR: Genotype type '" . $ln_data[1] . "' is not supported\n"
        if !$genotype_type_id;
      my $accession = $ln_data[2];
      $data{'feature'}{$accession}{'exists'} = 1;

      # to be used later for genotypes
      $data{'experiment'}{$experiment_name}{$accession}{'type_id'} =
        $genotype_type_id;
      push(
            @{ $data{'experiment'}{$experiment_name}{$accession}{'alleles'} },
            $ln_data[2] . '::' . $ln_data[3]
      );

# addition to the feature table
# the following ensure it only happens if we have featureloc data for this accession
      if ( $ln_data[4] && $ln_data[6] ) {
        my $ref = $ln_data[4];
        $data{'featureloc'}{$ref}{$accession}{'exists'} = 1;
        my ( $start, $end, $strand );
        if ( $ln_data[5] <= $ln_data[6] ) {
          my $start  = $ln_data[5];
          my $end    = $ln_data[6];
          my $strand = 1;
        } else {
          my $start  = $ln_data[6];
          my $end    = $ln_data[5];
          my $strand = -1;
        }
        $data{'featureloc'}{$ref}{$accession}{'start'}  = $start;
        $data{'featureloc'}{$ref}{$accession}{'end'}    = $end;
        $data{'featureloc'}{$ref}{$accession}{'strand'} = $strand;
      }
    }
    close IN;

    # add references first
    foreach my $id ( keys %{ $data{'featureloc'} } ) {
      push( @genbank_to_add, $id );
    }
    foreach my $id ( keys %{ $data{'feature'} } ) {
      push( @genbank_to_add, $id ) unless $data{'featureloc'}{$id};
    }
  }

  #add all sequences
  &load_sequences( \@genbank_to_add );

  # update feature with feature_ids
  foreach my $feature ( keys %{ $data{'feature'} } ) {
    $get_feat_id_prepared->execute($feature);
    my ($feature_id) = $get_feat_id_prepared->fetchrow_array();
    $data{'feature'}{$feature}{'feature_id'} = $feature_id;
  }
  ## add featureloc
  foreach my $ref ( keys %{ $data{'referenced'} } ) {
    foreach my $location ( keys %{ $data{'referenced'}{$ref} } ) {
      my $start  = $data{'featureloc'}{$ref}{$location}{'start'};
      my $end    = $data{'featureloc'}{$ref}{$location}{'start'};
      my $strand = $data{'featureloc'}{$ref}{$location}{'start'};
      $select_featureloc_prepared->execute( $ref, $location );
      my $loc_id = $select_featureloc_prepared->fetchrow_array();
      if ( !$loc_id ) {
        $insert_featureloc_prepared->execute( $ref, $location, $start, $end,
                                              $strand );
        $select_featureloc_prepared->execute( $ref, $location );
      }
      if ( !$loc_id ) {
        $select_featureloc_prepared->finish();
        &disconnectdb($dbh);
        die(
"ERROR: Cannot find nor add featureloc for srcfeature $ref and feature $location... Do you have write permissions?\n"
        );
      }
    }
  }

# create a genotype from
# @$data{'experiment'}{$experiment_name}{$accession}{'alleles'} = $data[3];
# each feature + genotype -> feature_genotype
# -> experiment_genotype
# genotype.type_id comes from $data{'experiment'}{ $experiment_name }{ $global_accession }{'type_id'}
# feature_id comes from $data{'feature'}{$feat}
  foreach my $experiment ( keys %{ $data{'experiment'} } ) {
    my $experiment_id = $experiments_hash_ref->{$experiment}->{'experiment_id'};
    if ( !$experiment_id ) {
      warn Dumper $experiments_hash_ref->{$experiment};
      warn
"Experiment ID for $experiment was not found. You are providing a new experiment that we have no geolocation for! Skipping....\n";
      next;
    }
    foreach my $marker ( keys %{ $data{'experiment'}{$experiment} } ) {
      my $feature_id = $data{'feature'}{$marker}{'feature_id'};
      my @alleles = @{ $data{'experiment'}{$experiment}{$marker}{'alleles'} };
      foreach my $allele (@alleles) {
        $check_genotype_prepared->execute($allele);
        my $genotype_id = $check_genotype_prepared->fetchrow_array();
        if ( !$genotype_id ) {

          # uniquename of genotype is now equal to feature::allele name...
          $insert_genotype_prepared->execute( $allele,
                         $data{'experiment'}{$experiment}{$marker}{'type_id'} );
          $check_genotype_prepared->execute($allele);
          $genotype_id = $check_genotype_prepared->fetchrow_array();
        }
        if ( !$genotype_id ) {
          die
"ERROR: Cannot add to database: genotype with uniquename=$allele and type="
            . $data{'experiment'}{$experiment}{$marker}{'type_id'}
            . " for experiment $experiment. Are you sure you have write permissions?\n";
        }
        $select_feature_genotype_prepared->execute( $feature_id, $genotype_id );
        my $res = $select_feature_genotype_prepared->fetchrow_array();
        if ( !$res ) {
          $insert_feature_genotype_prepared->execute( $feature_id,
                                                      $genotype_id );
          $select_feature_genotype_prepared->execute( $feature_id,
                                                      $genotype_id );
          $res = $select_feature_genotype_prepared->fetchrow_array();
        }
        if ( !$res ) {
          die
"ERROR: Cannot link feature and genotype in database: genotype with genotype_id=$genotype_id and feature_id=$feature_id for experiment $experiment. Are you sure you have write permissions?\n";
        }

#link to experiment via experiment_genotype with experiment.uniquename %{$data{'experiment'}}
        $select_experiment_genotype_prepared->execute( $experiment_id,
                                                       $genotype_id );
        my ($link_id) = $select_experiment_genotype_prepared->fetchrow_array();
        if ( !$link_id ) {
          $insert_experiment_genotype_prepared->execute( $experiment_id,
                                                         $genotype_id );
          $select_experiment_genotype_prepared->execute( $experiment_id,
                                                         $genotype_id );
          ($link_id) = $select_experiment_genotype_prepared->fetchrow_array();
          if ( !$link_id ) {
            die "ERROR: Cannot add to database: experiment_genotype with genotype_id= $genotype_id and experiment= $experiment_id "
              . "for experiment $experiment. Are you sure you have write permissions?\n";
          }
        }
        
        # here we are populating the array with the genotype info
        push(
              @{
                $experiments_hash_ref->{$experiment}->{'genotypes'}->{$marker}
                  ->{'genotype_ids'}
                },
              $genotype_id
        );
        push(
              @{
                $experiments_hash_ref->{$experiment}->{'genotypes'}->{$marker}
                  ->{'feature_ids'}
                },
              $feature_id
        );
      }
    }
  }

 # foreach experiment u.name update $experiments_hash_ref with 'genotypes' array. this is redundant as nd_experiment_genotype is used as well
  &link_stock_genotypes();
}

sub create_darwin_core_cvs() {

# this still needs a bit of curation. when it doubt, always choose a chado solution over a DC. The plan would be to create a connector to allow for DC queries.
# http://darwincore.googlecode.com/svn/trunk/terms/index.htm
  my $geolocation_metadata = {
    'Location' => {
      'continent' =>
'The name of the continent in which the Location occurs. Recommended best practice is to use a controlled vocabulary such as the Getty Thesaurus of Geographic Names or the ISO 3166 Continent code.',
      'county' =>
'The full, unabbreviated name of the next smaller administrative region than stateProvince (county, shire, department, etc.) in which the Location occurs.',
      'waterBody' =>
'The name of the water body in which the Location occurs. Recommended best practice is to use a controlled vocabulary such as the Getty Thesaurus of Geographic Names.',
      'islandGroup' =>
'The name of the island group in which the Location occurs. Recommended best practice is to use a controlled vocabulary such as the Getty Thesaurus of Geographic Names.',
      'island' =>
'The name of the island on or near which the Location occurs. Recommended best practice is to use a controlled vocabulary such as the Getty Thesaurus of Geographic Names.',
      'stateProvince' =>
'The name of the next smaller administrative region than country (state, province, canton, department, region, etc.) in which the Location occurs.',
      'country' =>
'The name of the country or major administrative unit in which the Location occurs. Recommended best practice is to use a controlled vocabulary such as the Getty Thesaurus of Geographic Names.',
      'municipality' =>
'The full, unabbreviated name of the next smaller administrative region than county (city, municipality, etc.) in which the Location occurs. Do not use this term for a nearby named place that does not contain the actual location.',
      'locality' =>
'The specific description of the place. Less specific geographic information can be provided in other geographic terms (higherGeography, continent, country, stateProvince, county, municipality, waterBody, island, islandGroup). This term may contain information modified from the original to correct perceived errors or standardize the description.',
      'verbatimLocality' =>
'The original textual specific description of the place. Less specific geographic information can be provided in other geographic terms (higherGeography, continent, country, stateProvince, county, municipality, waterBody, island, islandGroup).'
    },
    'GeologicalContext' => {
      'earliestEonOrLowestEonothem' =>
'The full name of the earliest possible geochronologic eon or lowest chrono-stratigraphic eonothem or the informal name ("Precambrian") attributable to the stratigraphic horizon from which the cataloged item was collected.',
      'latestEonOrHighestEonothem' =>
'The full name of the latest possible geochronologic eon or highest chrono-stratigraphic eonothem or the informal name ("Precambrian") attributable to the stratigraphic horizon from which the cataloged item was collected.',
      'earliestEraOrLowestErathem' =>
'The full name of the earliest possible geochronologic era or lowest chronostratigraphic erathem attributable to the stratigraphic horizon from which the cataloged item was collected.',
      'latestEraOrHighestErathem' =>
'The full name of the latest possible geochronologic era or highest chronostratigraphic erathem attributable to the stratigraphic horizon from which the cataloged item was collected.',
      'earliestPeriodOrLowestSystem' =>
'The full name of the earliest possible geochronologic period or lowest chronostratigraphic system attributable to the stratigraphic horizon from which the cataloged item was collected.',
      'latestPeriodOrHighestSystem' =>
'The full name of the latest possible geochronologic period or highest chronostratigraphic system attributable to the stratigraphic horizon from which the cataloged item was collected.',
      'earliestEpochOrLowestSeries' =>
'The full name of the earliest possible geochronologic epoch or lowest chronostratigraphic series attributable to the stratigraphic horizon from which the cataloged item was collected.',
      'latestEpochOrHighestSeries' =>
'The full name of the latest possible geochronologic epoch or highest chronostratigraphic series attributable to the stratigraphic horizon from which the cataloged item was collected.',
      'earliestAgeOrLowestStage' =>
'The full name of the earliest possible geochronologic age or lowest chronostratigraphic stage attributable to the stratigraphic horizon from which the cataloged item was collected.',
      'latestAgeOrHighestStage' =>
'The full name of the latest possible geochronologic age or highest chronostratigraphic stage attributable to the stratigraphic horizon from which the cataloged item was collected.',
      'lowestBiostratigraphicZone' =>
'The full name of the lowest possible geological biostratigraphic zone of the stratigraphic horizon from which the cataloged item was collected.',
      'highestBiostratigraphicZone' =>
'The full name of the highest possible geological biostratigraphic zone of the stratigraphic horizon from which the cataloged item was collected.',
      'lithostratigraphicTerms' =>
'The combination of all litho-stratigraphic names for the rock from which the cataloged item was collected.',
      'group' =>
'The full name of the lithostratigraphic group from which the cataloged item was collected.',
      'formation' =>
'The full name of the lithostratigraphic formation from which the cataloged item was collected.',
      'member' =>
'The full name of the lithostratigraphic member from which the cataloged item was collected.',
      'bed' =>
'The full name of the lithostratigraphic bed from which the cataloged item was collected.'
    }
  };
  my $activity_metadata = {
    'Event' => {
      'samplingProtocol' =>
'The name of, reference to, or description of the method or protocol used during an Event.',
      'samplingEffort' => 'The amount of effort expended during an Event.',
      'eventDate' =>
'The date-time or interval during which an Event occurred. For occurrences, this is the date-time when the event was recorded. Not suitable for a time in a geological context. Recommended best practice is to use an encoding scheme, such as ISO 8601:2004(E).',
      'eventTime' =>
'The time or interval during which an Event occurred. Recommended best practice is to use an encoding scheme, such as ISO 8601:2004(E).',
      'startDayOfYear' =>
'The earliest ordinal day of the year on which the Event occurred (1 for January 1, 365 for December 31, except in a leap year, in which case it is 366).',
      'endDayOfYear' =>
'The latest ordinal day of the year on which the Event occurred (1 for January 1, 365 for December 31, except in a leap year, in which case it is 366).',
      'year' =>
'The four-digit year in which the Event occurred, according to the Common Era Calendar.',
      'month' => 'The ordinal month in which the Event occurred.',
      'day'   => 'The integer day of the month on which the Event occurred.',
      'verbatimEventDate' =>
'The verbatim original representation of the date and time information for an Event.',
      'habitat' =>
        'A category or description of the habitat in which the Event occurred.',
      'fieldNumber' =>
'An identifier given to the event in the field. Often serves as a link between field notes and the Event.',
      'fieldNotes' =>
'One of a) an indicator of the existence of, b) a reference to (publication, URI), or c) the text of notes taken in the field about the Event.',
      'eventRemarks' => 'Comments or notes about the Event.',
    }
  };
  my $stock_metadata = {
    'Occurrence' => {
      'occurrenceID' =>
'An identifier for the Occurrence (as opposed to a particular digital record of the occurrence). In the absence of a persistent global unique identifier, construct one from a combination of identifiers in the record that will most closely make the occurrenceID globally unique.',
      'occurrenceDetails' =>
'A reference (publication, URI) to the most detailed information available about the Occurrence.',
      'occurrenceRemarks' => 'Comments or notes about the Occurrence.',
      'recordNumber' =>
"An identifier given to the Occurrence at the time it was recorded. Often serves as a link between field notes and an Occurrence record, such as a specimen collector's number.",
      'recordedBy' =>
'A list (concatenated and separated) of names of people, groups, or organizations responsible for recording the original Occurrence. The primary collector or observer, especially one who applies a personal identifier (recordNumber), should be listed first.',
      'individualID' =>
'An identifier for an individual or named group of individual organisms represented in the Occurrence. Meant to accommodate resampling of the same individual or group for monitoring purposes. May be a global unique identifier or an identifier specific to a data set.',
      'individualCount' =>
'The number of individuals represented present at the time of the Occurrence.',
      'sex' =>
'The sex of the biological individual(s) represented in the Occurrence. Recommended best practice is to use a controlled vocabulary.',
      'lifeStage' =>
'The age class or life stage of the biological individual(s) at the time the Occurrence was recorded. Recommended best practice is to use a controlled vocabulary.',
      'reproductiveCondition' =>
'The reproductive condition of the biological individual(s) represented in the Occurrence. Recommended best practice is to use a controlled vocabulary.',
      'behavior' =>
'A description of the behavior shown by the subject at the time the Occurrence was recorded. Recommended best practice is to use a controlled vocabulary.',
      'establishmentMeans' =>
'The process by which the biological individual(s) represented in the Occurrence became established at the location. Recommended best practice is to use a controlled vocabulary.',
      'occurrenceStatus' =>
'A statement about the presence or absence of a Taxon at a Location. Recommended best practice is to use a controlled vocabulary.',
      'preparations' =>
'A list (concatenated and separated) of preparations and preservation methods for a specimen.',
      'disposition' =>
'The current state of a specimen with respect to the collection identified in collectionCode or collectionID. Recommended best practice is to use a controlled vocabulary.',
      'otherCatalogNumbers' =>
'A list (concatenated and separated) of previous or alternate fully qualified catalog numbers or other human-used identifiers for the same Occurrence, whether in the current or any other data set or collection.',
      'previousIdentifications' =>
'A list (concatenated and separated) of previous assignments of names to the Occurrence.',
      'associatedMedia' =>
'A list (concatenated and separated) of identifiers (publication, global unique identifier, URI) of media associated with the Occurrence.',
      'associatedReferences' =>
'A list (concatenated and separated) of identifiers (publication, bibliographic reference, global unique identifier, URI) of literature associated with the Occurrence.',
      'associatedOccurrences' =>
'A list (concatenated and separated) of identifiers of other Occurrence records and their associations to this Occurrence.',
      'associatedSequences' =>
'A list (concatenated and separated) of identifiers (publication, global unique identifier, URI) of genetic sequence information associated with the Occurrence.',
      'associatedTaxa' =>
'A list (concatenated and separated) of identifiers or names of taxa and their associations with the Occurrence.'
    },
    'Identification' => {
      'identifiedBy' =>
'A list (concatenated and separated) of names of people, groups, or organizations who assigned the Taxon to the subject.',
      'dateIdentified' =>
'The date on which the subject was identified as representing the Taxon. Recommended best practice is to use an encoding scheme, such as ISO 8601:2004(E).',
      'identificationReferences' =>
'A list (concatenated and separated) of references (publication, global unique identifier, URI) used in the Identification.',
      'identificationRemarks' => 'Comments or notes about the Identification.',
      'identificationQualifier' =>
"A brief phrase or a standard term (\"cf.\", \"aff.\") to express the determiner's doubts about the Identification.",
      'typeStatus' =>
'A list (concatenated and separated) of nomenclatural types (type status, typified scientific name, publication) applied to the subject.'
    }
  };
  my $organism_metadata = {
    'Taxon' => {
      'scientificNameID' =>
'An identifier for the nomenclatural (not taxonomic) details of a scientific name.',
      'acceptedNameUsageID' =>
'An identifier for the name usage (documented meaning of the name according to a source) of the currently valid (zoological) or accepted (botanical) taxon.',
      'parentNameUsageID' =>
'An identifier for the name usage (documented meaning of the name according to a source) of the direct, most proximate higher-rank parent taxon (in a classification) of the most specific element of the scientificName.',
      'originalNameUsageID' =>
'An identifier for the name usage (documented meaning of the name according to a source) in which the terminal element of the scientificName was originally established under the rules of the associated nomenclaturalCode.',
      'nameAccordingToID' =>
'An identifier for the source in which the specific taxon concept circumscription is defined or implied. See nameAccordingTo.',
      'namePublishedInID' =>
'An identifier for the publication in which the scientificName was originally established under the rules of the associated nomenclaturalCode.',
      'taxonConceptID' =>
'An identifier for the taxonomic concept to which the record refers - not for the nomenclatural details of a taxon.',
      'scientificName' =>
'The full scientific name, with authorship and date information if known. When forming part of an Identification, this should be the name in lowest level taxonomic rank that can be determined. This term should not contain identification qualifications, which should instead be supplied in the IdentificationQualifier term.',
      'acceptedNameUsage' =>
'The full name, with authorship and date information if known, of the currently valid (zoological) or accepted (botanical) taxon.',
      'parentNameUsage' =>
'The full name, with authorship and date information if known, of the direct, most proximate higher-rank parent taxon (in a classification) of the most specific element of the scientificName.',
      'originalNameUsage' =>
'The taxon name, with authorship and date information if known, as it originally appeared when first established under the rules of the associated nomenclaturalCode. The basionym (botany) or basonym (bacteriology) of the scientificName or the senior/earlier homonym for replaced names.',
      'nameAccordingTo' =>
'The reference to the source in which the specific taxon concept circumscription is defined or implied - traditionally signified by the Latin "sensu" or "sec." (from secundum, meaning "according to"). For taxa that result from identifications, a reference to the keys, monographs, experts and other sources should be given.',
      'namePublishedIn' =>
'A reference for the publication in which the scientificName was originally established under the rules of the associated nomenclaturalCode.',
      'higherClassification' =>
'A list (concatenated and separated) of taxa names terminating at the rank immediately superior to the taxon referenced in the taxon record. Recommended best practice is to order the list starting with the highest rank and separating the names for each rank with a semi-colon (";").',
      'kingdom' =>
'The full scientific name of the kingdom in which the taxon is classified.',
      'phylum' =>
'The full scientific name of the phylum or division in which the taxon is classified.',
      'class' =>
'The full scientific name of the class in which the taxon is classified.',
      'order' =>
'The full scientific name of the order in which the taxon is classified.',
      'family' =>
'The full scientific name of the family in which the taxon is classified.',
      'genus' =>
'The full scientific name of the genus in which the taxon is classified.',
      'subgenus' =>
'The full scientific name of the subgenus in which the taxon is classified. Values should include the genus to avoid homonym confusion.',
      'specificEpithet' =>
        'The name of the first or species epithet of the scientificName.',
      'infraspecificEpithet' =>
'The name of the lowest or terminal infraspecific epithet of the scientificName, excluding any rank designation.',
      'taxonRank' =>
'The taxonomic rank of the most specific name in the scientificName. Recommended best practice is to use a controlled vocabulary.',
      'verbatimTaxonRank' =>
'The taxonomic rank of the most specific name in the scientificName as it appears in the original record.',
      'scientificNameAuthorship' =>
'The authorship information for the scientificName formatted according to the conventions of the applicable nomenclaturalCode.',
      'vernacularName' => 'A common or vernacular name.',
      'nomenclaturalCode' =>
'The nomenclatural code (or codes in the case of an ambiregnal name) under which the scientificName is constructed. Recommended best practice is to use a controlled vocabulary.',
      'taxonomicStatus' =>
'The status of the use of the scientificName as a label for a taxon. Requires taxonomic opinion to define the scope of a taxon. Rules of priority then are used to define the taxonomic status of the nomenclature contained in that scope, combined with the experts opinion. It must be linked to a specific taxonomic reference that defines the concept. Recommended best practice is to use a controlled vocabulary.',
      'nomenclaturalStatus' =>
'The status related to the original publication of the name and its conformance to the relevant rules of nomenclature. It is based essentially on an algorithm according to the business rules of the code. It requires no taxonomic opinion.',
      'taxonRemarks' => 'Comments or notes about the taxon or name.'
    },
  };
  my $cvterm_metadata = {
    'MeasurementOrFact' => {
      'measurementType' =>
'The nature of the measurement, fact, characteristic, or assertion. Recommended best practice is to use a controlled vocabulary.',
      'measurementValue' =>
        'The value of the measurement, fact, characteristic, or assertion.',
      'measurementAccuracy' =>
'The description of the potential error associated with the measurementValue.',
      'measurementUnit' =>
'The units associated with the measurementValue. Recommended best practice is to use the International System of Units (SI).',
      'measurementDeterminedDate' =>
'The date on which the MeasurementOrFact was made. Recommended best practice is to use an encoding scheme, such as ISO 8601:2004(E).',
      'measurementDeterminedBy' =>
'A list (concatenated and separated) of names of people, groups, or organizations who determined the value of the MeasurementOrFact.',
      'measurementMethod' =>
'A description of or reference to (publication, URI) the method or protocol used to determine the measurement, fact, characteristic, or assertion.',
      'measurementRemarks' =>
        'Comments or notes accompanying the MeasurementOrFact.'
    }
  };

#die Dumper  (           $geolocation_metadata, $activity_metadata,           $stock_metadata,       $organism_metadata,           $cvterm_metadata  );
# foreach CV add it to chado and rewrite description with cvterm_id
  foreach my $cv ( keys %{$geolocation_metadata} ) {
    foreach my $cvterm ( keys %{ $geolocation_metadata->{$cv} } ) {
      $geolocation_metadata->{$cv}->{$cvterm} =
        &get_set_cvterm( $cv, $cvterm, 'dcterms#' . $cv,
                         $cvterm, $geolocation_metadata->{$cv}->{$cvterm} );
    }
  }
  foreach my $cv ( keys %{$activity_metadata} ) {
    foreach my $cvterm ( keys %{ $activity_metadata->{$cv} } ) {
      $activity_metadata->{$cv}->{$cvterm} =
        &get_set_cvterm( $cv, $cvterm, 'dcterms#' . $cv,
                         $cvterm, $activity_metadata->{$cv}->{$cvterm} );
    }
  }
  foreach my $cv ( keys %{$stock_metadata} ) {
    foreach my $cvterm ( keys %{ $stock_metadata->{$cv} } ) {
      $stock_metadata->{$cv}->{$cvterm} =
        &get_set_cvterm( $cv, $cvterm, 'dcterms#' . $cv,
                         $cvterm, $stock_metadata->{$cv}->{$cvterm} );
    }
  }
  foreach my $cv ( keys %{$organism_metadata} ) {
    foreach my $cvterm ( keys %{ $organism_metadata->{$cv} } ) {
      $organism_metadata->{$cv}->{$cvterm} =
        &get_set_cvterm( $cv, $cvterm, 'dcterms#' . $cv,
                         $cvterm, $organism_metadata->{$cv}->{$cvterm} );
    }
  }
  foreach my $cv ( keys %{$cvterm_metadata} ) {
    foreach my $cvterm ( keys %{ $cvterm_metadata->{$cv} } ) {
      $cvterm_metadata->{$cv}->{$cvterm} =
        &get_set_cvterm( $cv, $cvterm, 'dcterms#' . $cv,
                         $cvterm, $cvterm_metadata->{$cv}->{$cvterm} );
    }
  }
  return (
           $geolocation_metadata, $activity_metadata, $stock_metadata,
           $organism_metadata,    $cvterm_metadata
  );
}

sub other_nd_cvs() {
  my $experiment_cv = {
    'experiment_type' => {
      'natural_collection' =>
'A collection from the natural environment of the organism, i.e. where this species is expected to occur.',
      'lab_collection' =>
'A collection from a laboratory stock of the organism, i.e. not in this species\' natural environment.',
      'genotyping'  => 'Characterization of a haplotype.',
      'phenotyping' => 'Characterization of a phenotype.',
      'generic' => 'Placeholder for experiments whose type must be considered.',
    }
  };
  my $stock_cv = {
    'stock_type' => {

      #flybase:
      'living stock' => 'A living, reproducing culture of organisms.',
      'molecular extract' =>
        'Molecular material extracted from a living stock.',
      'genomic DNA' => 'Total DNA extracted from a living stock.',
      'preserved specimen' =>
        'A dead organism or tissue treated to prevent decay.',
      'desiccated specimen' =>
        'An organism or part of an organism preserved by desiccation.',
      'ethanol-preserved specimen' =>
        'An organism or part of an organism preserved in ethanol.',
      'frozen specimen' =>
        'An organism or part of an organism preserved at low temperature.',

      #own:
      'lab_stock_reared' =>
'Single specimen is derived from a standard laboratory stock that is available from a stock center.',
      'lab_custom_reared' =>
'Single specimen is derived from a laboratory stock not available from a stock center.',
      'natural_history_museum' =>
'Single specimen is derived from a museum collection and is therefore considered unique (i.e. is a voucher).',
      'field_collected_once' =>
'Single specimen is derived from the field but was not reared/cultured in the lab. It is therefore considered unique but it may or may not exist as a voucher (it may have been destroyed).',
      'field_collected_reared' =>
'Single specimen is derived from the field and was reared/cultured in the lab. The specimen is unique but it may or may not exist as a voucher (it may have been destroyed). Its progeny are available (and will be type lab_custom_reared).',
      'host_mass_collected_once' =>
'A number of specimens were mass derived from a host (e.g. pathogens). Specimens were processed without separation. No cultivation occured and specimen may or may not exist as a voucher (it may have been destroyed).',
      'host_collected_reared' =>
'A number of specimens were mass derived from a host (e.g. pathogens). Specimens were processed without separation. Specimens may or may not exist as a voucher (it may have been destroyed) but its a laboratory colony of the progeny exists.',
      'host_mass_collected_separated_once' =>
'A number of specimens were mass derived from a host (e.g. pathogens). Specimens were separated to - or select as - haplotypes before prior to phenotype/genotype experiment. No cultivation occured and specimen may or may not exist as a voucher (it may have been destroyed).',
      'host_collected_separated_reared' =>
'A number of specimens were mass derived from a host (e.g. pathogens). Specimens were separated to - or select as - haplotypes before prior to phenotype/genotype experiment. Specimens may or may not exist as a voucher (it may have been destroyed) but its a laboratory colony of the progeny exists.',
      'generic' => 'Placeholder for stocks whose type must be considered.',
    }
  };
  my $relationship_cv = {
    'relationship' => {
      'produces' =>
'Subject creates the object. For example a collection experiment produces a collected specimen.',
      'produced_by' =>
'Subject is created by the object. For example a collection specimen was produced during a collection experiment.',
    }
  };
  my $genotype_cv = {
    ## select the most informative of these. you may add more as genotypeprop e.g. with type_id='homozygote',value=1
    'genotype_type' => {
                         'homozygote'   => '',
                         'heterozygote' => '',
                         'resistant'    => '',
                         'sensitive'    => '',
                         'dominant'     => '',
                         'recessive'    => '',
                         'co_dominant'  => '',
                         'haplotype'    => '',
    }
  };

  #die Dumper ($experiment_cv,$stock_cv);
  foreach my $cv ( keys %{$experiment_cv} ) {
    foreach my $cvterm ( keys %{ $experiment_cv->{$cv} } ) {
      $experiment_cv->{$cv}->{$cvterm} =
        &get_set_cvterm( $cv, $cvterm, $cv, $cvterm,
                         $experiment_cv->{$cv}->{$cvterm} );
    }
  }
  foreach my $cv ( keys %{$stock_cv} ) {
    foreach my $cvterm ( keys %{ $stock_cv->{$cv} } ) {
      $stock_cv->{$cv}->{$cvterm} =
        &get_set_cvterm( $cv, $cvterm, $cv, $cvterm,
                         $stock_cv->{$cv}->{$cvterm} );
    }
  }
  foreach my $cv ( keys %{$genotype_cv} ) {
    foreach my $cvterm ( keys %{ $genotype_cv->{$cv} } ) {
      $genotype_cv->{$cv}->{$cvterm} =
        &get_set_cvterm( $cv, $cvterm, $cv, $cvterm,
                         $genotype_cv->{$cv}->{$cvterm} );
    }
  }
  foreach my $cv ( keys %{$relationship_cv} ) {
    foreach my $cvterm ( keys %{ $relationship_cv->{$cv} } ) {
      $relationship_cv->{$cv}->{$cvterm} =
        &get_set_cvterm( $cv, $cvterm, $cv, $cvterm,
                         $relationship_cv->{$cv}->{$cvterm} );
    }
  }

  # create new tables: environmentprop and nd_geolocation_environment
  my $check_table_environmentprop =
    $dbh->table_info( "", 'public', 'environmentprop', "TABLE" );
  if ( $check_table_environmentprop->fetch ) { }
  else {
    $dbh->do(
"SELECT nd_geolocation_environment_id from nd_geolocation_environment limit 1"
    );
    $dbh->do(
"CREATE TABLE environmentprop (environmentprop_id serial primary key,environment_id integer NOT NULL, type_id integer NOT NULL,  value text,  rank integer DEFAULT 0 NOT NULL )"
    );
    $dbh->do(
"ALTER TABLE ONLY environmentprop ADD CONSTRAINT environmentprop_c1 UNIQUE (environment_id, type_id, rank)"
    );
    $dbh->do(
"CREATE INDEX environmentprop_idx1 ON environmentprop USING btree (environment_id)"
    );
    $dbh->do(
"CREATE INDEX environmentprop_idx2 ON environmentprop USING btree (type_id)" );
    $dbh->do(
"ALTER TABLE ONLY environmentprop ADD CONSTRAINT environmentprop_environment_id_fkey FOREIGN KEY (environment_id) REFERENCES environment(environment_id) ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED"
    );
    $dbh->do(
"ALTER TABLE ONLY environmentprop ADD CONSTRAINT environmentprop_type_id_fkey FOREIGN KEY (type_id) REFERENCES cvterm(cvterm_id) ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED"
    );
  }
  my $create_table_nd_geolocation_environment =
    $dbh->table_info( "", 'public', 'nd_geolocation_environment', "TABLE" );
  if ( $create_table_nd_geolocation_environment->fetch ) { }
  else {
    $dbh->do(
"CREATE TABLE nd_geolocation_environment (nd_geolocation_environment_id serial primary key,nd_geolocation_id integer not NULL,environment_id integer not NULL,datesampled date)"
    );
    $dbh->do(
"ALTER TABLE ONLY nd_geolocation_environment ADD CONSTRAINT nd_geolocation_environment_nd_geolocation_id_fkey FOREIGN KEY (nd_geolocation_id) REFERENCES nd_geolocation(nd_geolocation_id) ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED"
    );
    $dbh->do(
"CREATE UNIQUE INDEX nd_geolocation_environment_c1 ON nd_geolocation_environment USING btree (nd_geolocation_id,environment_id)"
    );
    $dbh->do(
"ALTER TABLE ONLY nd_geolocation_environment ADD CONSTRAINT nd_geolocation_environment_environment_id_fkey FOREIGN KEY (environment_id) REFERENCES environment(environment_id) ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED"
    );
  }
  return ( $experiment_cv, $stock_cv, $genotype_cv, $relationship_cv );
}
####

=pod

=head2 extract_geo

 IN: We expect a CSV with these columns a unique ID, experiment_type, Long and Lat.
     If the unique id has a space, then it is split to name + description. uname and description (if it exists) are stored in the environment table
     If there is a 5th column, then it is an altitude
     If there is a 5th column, then it is a description
     If there is no description, the unique ID is used as 'Location of X'
     Unique ID is enforced to force people to provide location identifiers
     lines starting with # are ignored (e.g headers)
 Process: Store in nd_geolocation
 
=cut
sub extract_store_geolocation() {
  my $infile                = shift;
  my $collection_project    = shift;
  my $collection_project_id = shift;
  my %geolocations
    ; # key is uname of experiment, key-value pair will be nd_geolocation_id and type_id
  die "ERROR: Cannot find/open file $infile\n" unless $infile && -s $infile;
  open( IN, $infile );

  # check if there is
  while ( my $ln = <IN> ) {
    next if ( $ln =~ /^#/ || $ln =~ /^\s*$/ );
    chomp($ln);
    my @data = split( $delimiter, $ln );
    next if !$data[4];
    my $experiment_name = $collection_project . '::' . $data[0];
    if ( $geolocations{$experiment_name} ) {
      warn "Warning: ID "
        . $experiment_name
        . " already exists but it is supposed to be a unique ID. This line will be ignored:\n$ln\n";
      sleep(0.001);
    }
    my $experimental_stock_type_name = $data[2];
    my $experimental_stock_type_id =
      $stock_cv->{'stock_type'}->{$experimental_stock_type_name};
    die "Stock type '" . $experimental_stock_type_name . "' is not supported\n"
      unless $experimental_stock_type_id;
    $geolocations{$experiment_name}{'stock_type_id'} =
      $experimental_stock_type_id;
    $geolocations{$experiment_name}{'stock_type_name'} =
      $experimental_stock_type_name;
    $geolocations{$experiment_name}{'collection_project'} = $collection_project;
    $geolocations{$experiment_name}{'collection_project_id'} =
      $collection_project_id;
    my $description = $1 if $1;
    my $long        = $data[3];
    my $lat         = $data[4];

    if ( $long =~ /[^\d\.\-\+WE]/ || $lat =~ /[^\d\.\-\+NS]/ ) {
      die
"ERROR: Longitude and latitude columns should only contain numbers and + - . NSWE characters\n"
        . $long . "\n"
        . $lat . "\n";
    }
    my $geoexperiment_type_id = $experiment_cv->{'experiment_type'}{ $data[1] };
    if ( !$geoexperiment_type_id ) {
      warn( "Error: Experiment type " . $data[1] . " is not supported.\n" );
      die Dumper $experiment_cv->{'experiment_type'};
    }
    $select_environment_prepared->execute($experiment_name);
    my $enviro_id = $select_environment_prepared->fetchrow_array();
    if ( !$enviro_id ) {
      $insert_environment_prepared->execute($experiment_name);
      $select_environment_prepared->execute($experiment_name);
      $enviro_id = $select_environment_prepared->fetchrow_array();
    }
    if ( !$enviro_id ) {
      die
"ERROR: Cannot add to db the new environment with uniquename=$experiment_name. Do you have write permissions?\n";
    }
    $geolocations{$experiment_name}{'environment_id'} = $enviro_id;
    my $alt  = $data[4] ? $data[5] : undef;
    my $desc = $data[5] ? $data[6] : 'Location of ' . $experiment_name;
    $geolocations{$experiment_name}{'geolocation_id'} =
      &_add_geolocation( $experiment_name, $long, $lat, $alt, $desc );
    $geolocations{$experiment_name}{'type_id'} = $geoexperiment_type_id;
  }
  close(IN);
  return \%geolocations;
}

sub _add_geolocation() {
  my $id   = shift;
  my $long = shift;
  my $lat  = shift;
  my $alt  = shift;
  my $desc = shift;
  my $geolocation_id;    #returned;
  $select_geolocation_id_prepared->execute( $lat, $long, $geodetic );
  $geolocation_id = $select_geolocation_id_prepared->fetchrow_array();

  if ($geolocation_id) {

#warn "Warn: Geolocation with lat $lat long $long geodetic $geodetic already exists in database\n";
  } else {
    $insert_geolocation_prepared->execute( $lat, $long, $geodetic, $alt );
    $select_geolocation_id_prepared->execute( $lat, $long, $geodetic );
    $geolocation_id = $select_geolocation_id_prepared->fetchrow_array();
  }
  if ( !$geolocation_id ) {
    &disconnectdb($dbh);
    die(
"ERROR: Cannot find nor add GEOLOCATION with lat= $lat long= $long geodetic= $geodetic to database... Do you have write permissions?"
    );
  }
  if ($desc) {
    $update_geolocation_description_prepared->execute( $desc, $geolocation_id );
  }
  return $geolocation_id;
}
###
sub connecttodb() {
  my $dsn = shift;
  if ( $dsn && -f $dsn ) {
    open( IN, $dsn );
    $dsn = "";
    while ( my $line = <IN> ) {
      unless ( $line =~ /^\s*$/ ) { chomp($line); $dsn .= $line; }
    }
    close(IN);
  }
  if ( !$dsn ) {
    die "ERROR: No connection to database provided.\n";
  }
  my $dbh = DBI->connect($dsn);
  die("ERROR: Unable to connect to database") if ( !$dbh );
  return ($dbh);
}
###
sub disconnectdb ($) {
  my $dbh = shift;
  $check_project_prepared->finish();
  $get_project_prepared->finish();
  $insert_project_prepared->finish();
  $get_projectrel_prepared->finish();
  $insert_projectrel_prepared->finish();
  $select_experiment_project_prepared->finish();
  $insert_experiment_project_prepared->finish();
  $select_stock_genotype_prepared->finish();
  $get_feat_name_prepared->finish();
  $select_experiment_stock_prepared->finish();
  $select_featureloc_prepared->finish();
  $check_genotype_prepared->finish();
  $check_phenotype_prepared->finish();
  $select_environment_prepared->finish();
  $select_geolocation_id_prepared->finish();
  $select_experiment_id_prepared->finish();
  $select_experiment_genotype_prepared->finish();
  $select_experiment_phenotype_prepared->finish();
  $select_feature_genotype_prepared->finish();
  $select_cvterm_via_dbxref_prepared->finish();
  $select_stock_prepared->finish();
  $next_stock_prepared->finish();
  $select_cv_prepared->finish();
  $select_geolocation_id_prepared->finish();
  $select_test_cvterm_prepared->finish();
  $select_cvterm_prepared->finish();
  $get_feat_id_prepared->finish();
  $select_featureprop_rank_prepared->finish();
  $select_featureprop_id_prepared->finish();
  $check_featureprop_prepared->finish();
  $select_dbxref_check_prepared->finish();
  $select_dbxref_prepared_cached->finish();
  $select_feature_dbxref_prepared->finish();
  $cvterms_list_prepared->finish();
  $dbxref_list_prepared->finish();
  $select_global_db_prepared->finish();
  $select_global_dbxref_prepared->finish();
  $dbh->disconnect();
}
###
sub get_ncbi_taxid ($) {
  my $species_latin = shift;
  my $ncbi_taxid;
  if ( $species_latin =~ /^\d+$/ ) {
    $ncbi_taxid = $species_latin;
    return $ncbi_taxid;
  }

  #at least three letters
  elsif ( $species_latin =~ /^([A-Za-z]{2}[a-z]+)\s*/ ) {
    $ncbi_taxid = $taxondb->get_taxonid($species_latin);
    if ($ncbi_taxid) {
      return $ncbi_taxid;
    }
    my $binomial = $1;
    if ( $species_latin =~ /^[A-Za-z]{2}[a-z]+\s+(\w+)/ ) {
      $binomial .= ' ' . $1;
    }
    if ($binomial) {
      $ncbi_taxid = $taxondb->get_taxonid($binomial);
      if ( !$ncbi_taxid ) {
        die "ERROR: Could not find NCBI taxid for $species_latin.\n";
      } else {
        return $ncbi_taxid;

#my $taxon = $taxondb->get_taxon(-taxonid => $ncbi_taxid)  || die "Cannot find your query\n";;
      }
    }
  }
}
###
sub get_set_dbxref() {
  my $dbxref_accession = shift || die("No accession provided");
  my $db_name          = shift || die("No db name provided");
  my ($dbxref_id);
  $db_name =~ s/^DB://;

  # get db_id
  my $r = $select_global_db_prepared->execute($db_name);
  my ( $db_id, $db_name_c ) = $select_global_db_prepared->fetchrow_array();
  undef($r);
  if ( !$db_id ) {

    # insert
    $insert_into_db_prepared->execute($db_name);
    $r = $select_global_db_prepared->execute($db_name);
    ( $db_id, $db_name_c ) = $select_global_db_prepared->fetchrow_array();
  }
  if ( !$db_id ) {
    &disconnectdb($dbh);
    die(
"ERROR: Cannot find nor add DB with db.name= $db_name to database... Do you have write permissions?"
    );
  }

  # Get dbxref_id
  undef($r);
  $select_global_dbxref_prepared->execute( $dbxref_accession, $db_name_c );
  $r = $select_global_dbxref_prepared->execute( $dbxref_accession, $db_name_c );
  ($dbxref_id) = $select_global_dbxref_prepared->fetchrow_array();
  if ( !$dbxref_id ) {

    # insert
    $insert_into_dbxref_prepared->execute( $dbxref_accession, $db_id, "NULL" );
    $select_global_dbxref_prepared->execute( $dbxref_accession, $db_name_c );
    $r =
      $select_global_dbxref_prepared->execute( $dbxref_accession, $db_name_c );
    ($dbxref_id) = $select_global_dbxref_prepared->fetchrow_array();
  }
  if ( !$dbxref_id ) {
    &disconnectdb($dbh);
    die(
"ERROR: Cannot find nor add DBXREF with dbxref_accession=$dbxref_accession db= $db_name to database... Do you have write permissions?"
    );
  }
  return $dbxref_id;
}

sub get_cvterm() {
  my $cv_name     = shift || die("No cv name provided");
  my $cvterm_name = shift || die("No cvterm name provided");
  $select_cvterm_prepared->execute( $cvterm_name, $cv_name );
  my $cvterm_id = $select_cvterm_prepared->fetchrow_array();
  return $cvterm_id if $cvterm_id;
}

sub get_set_cvterm () {
  my $cv_name     = shift || die("No cv name provided");
  my $cvterm_name = shift || die("No accession provided");
  my $db_name     = shift;
  my $dbxref_accession = shift;
  my $definition       = shift;
  $definition       = 'NULL'       if !$definition;
  $dbxref_accession = $cvterm_name if !$dbxref_accession;
  $db_name          = $cv_name     if !$db_name;
  my ($cvterm_id);

  # does it exist in the cvterm table?
  $select_cvterm_prepared->execute( $cvterm_name, $cv_name );
  $cvterm_id = $select_cvterm_prepared->fetchrow_array();
  return $cvterm_id if $cvterm_id;

  # Do we have a dbxref? if not, add it
  my $dbxref_id = &get_set_dbxref( $dbxref_accession, $db_name );
  die
"ERROR: cannot find dbxref_id for accession=$dbxref_accession of db.name=$db_name\n"
    unless $dbxref_id;

  #Do we have a CV? If not add it
  my $r     = $select_cv_prepared->execute($cv_name);
  my $cv_id = $select_cv_prepared->fetchrow_array();
  if ( !$cv_id ) {
    $insert_cv_prepared->execute($cv_name);
    $select_cv_prepared->execute($cv_name);
    $cv_id = $select_cv_prepared->fetchrow_array();
  }
  if ( !$cv_id ) {
    &disconnectdb($dbh);
    die(
"ERROR: Cannot find nor add CV with name= $cv_name database... Do you have write permissions?"
    );
  }

  # Do we have a cvterm with that dbxref?
  $r         = $select_cvterm_via_dbxref_prepared->execute($dbxref_id);
  $cvterm_id = $select_cvterm_via_dbxref_prepared->fetchrow_array();
  if ( !$cvterm_id ) {

    #insert it
    $insert_cvterm_prepared->execute( $cv_name,    $cvterm_name,
                                      $definition, $dbxref_id );
    $r         = $select_cvterm_via_dbxref_prepared->execute($dbxref_id);
    $cvterm_id = $select_cvterm_via_dbxref_prepared->fetchrow_array();
  }
  if ( !$cvterm_id ) {
    &disconnectdb($dbh);
    die(
"ERROR: Cannot find nor add CVTERM with dbxref_accession= $dbxref_accession cv= $cv_name db= $db_name to database... Do you have write permissions?"
    );
  }
  return $cvterm_id;
}

sub prepare_genbank() {
  my $ref = shift;
  return if !$ref;
  my $outfile    = "features_to_load.gff";
  my @accessions = @$ref;
  my $FORMAT     = "GenBank";
  my $SOURCEID   = "genbank";
  my %TAG_MAP = (
                  db_xref => 'Dbxref',
                  name    => 'Name',
                  note    => 'Note',
                  synonym => 'Alias',
                  symbol  => 'Alias',
  );
  $FORMAT = "swiss" if $FORMAT =~ /UniProt|trembl/;
  my $unflattener = Bio::SeqFeature::Tools::Unflattener->new;
  my $tm          = Bio::SeqFeature::Tools::TypeMapper->new;
  $genbank_script_idh = Bio::SeqFeature::Tools::IDHandler->new;
  my $source_type ||= "region";
  my ( %terms, %syn );
  my $get_so_terms_sql =
"SELECT name,cvterm_id from cvterm where cv_id=(SELECT cv_id from cv where name='sequence')";
  my $get_so_synonyms_sql =
    "SELECT synonym from cvtermsynonym where cvterm_id=?";
  my $get_so_terms_prepared    = $dbh->prepare($get_so_terms_sql);
  my $get_so_synonyms_prepared = $dbh->prepare($get_so_synonyms_sql);
  $get_so_terms_prepared->execute();

  while ( my ( $cvterm, $cvterm_id ) = $get_so_terms_prepared->fetchrow_array )
  {
    $terms{$cvterm} = $cvterm;
    $get_so_synonyms_prepared->execute($cvterm_id);
    while ( my ($synonym) = $get_so_synonyms_prepared->fetchrow_array ) {
      $syn{$synonym} = $cvterm;
    }
  }
  my $FTSOmap      = \%terms;
  my $FTSOsynonyms = \%syn;
  my %hardTerms    = %{ $tm->FT_SO_map() };
  map { $FTSOmap->{$_} ||= $hardTerms{$_} } keys %hardTerms;
  my $TYPE_MAP = $FTSOmap;
  my $SYN_MAP  = $FTSOsynonyms;
  open my $out, ">$outfile";
  print $out "##gff-version 3\n";
  open my $lumpfa_fh, ">$outfile.fsa";
  my $gb = Bio::DB::GenBank->new();

  foreach my $accession (@accessions) {
    $get_feat_id_prepared->execute($accession);
    my ($feature_id) = $get_feat_id_prepared->fetchrow_array();
    if ($feature_id) {
      warn "Feature $accession already exists in the database\n";
      next;
    } else {
      print
"Processing $accession after 3 seconds of sleep (being nice to NCBI webservice)...\n";
      sleep(3);    # being nice to NCBI
    }
    my $seq      = $gb->get_Seq_by_id($accession);
    my $gffio    = Bio::Tools::GFF->new( -noparse => 1, -gff_version => 3 );
    my $seq_name = $seq->accession_number;
    my $end      = $seq->length;
    my ( @to_print, @GFF_LINE_FEAT );
    my $source_feat = undef;
    my $organism    = undef;
    my @source      = _filter($seq);
    $source_feat = $source[0];
    ( $source_type, $source_feat ) =
      _getSourceInfo( $seq, $source_type, $source_feat );
    warn "$seq_name has no features, skipping\n" and next
      if !$seq->all_SeqFeatures;
    $FTSOmap->{'source'} = $source_type;
    _unflatten_seq( $seq, $unflattener, $FTSOmap, $FTSOsynonyms, $tm );

    for my $feature ( _get_all_SeqFeatures($seq) ) {
      my $method = $feature->primary_tag;
      next unless $method eq 'region';
      $feature->seq_id( $seq->id ) unless ( $feature->seq_id );
      $feature->source_tag($SOURCEID);
      _maptags2gff( $feature, \%TAG_MAP );
      my ( $_gene_name, $gene_id );
      $_gene_name = _gene_name($feature);
      _convert_to_name($feature);
      ($organism) = $feature->get_tag_values('organism')
        if $feature->has_tag('organism');
      _add_generic_id( $feature, $_gene_name, "" );
      my $gff = $gffio->gff_string($feature);
      push @GFF_LINE_FEAT, $feature;
    }
    @to_print = _print_held( $out, $gffio, \@to_print, @GFF_LINE_FEAT );
    _gff_validate(@GFF_LINE_FEAT);
    for my $feature (@GFF_LINE_FEAT) {
      my $gff = $gffio->gff_string($feature);
      print $out "$gff\n" if $gff;
    }
    my $dna = $seq->seq;
    if ( $dna || %genbank_script_proteinfa ) {
      $genbank_script_method{'RESIDUES'} += length($dna);
      $dna =~ s/(\S{60})/$1\n/g;
      $dna .= "\n";
      print $lumpfa_fh ">$seq_name\n$dna" if $dna;
      foreach my $aid ( sort keys %genbank_script_proteinfa ) {
        my $aa = delete $genbank_script_proteinfa{$aid};
        $genbank_script_method{'RESIDUES(tr)'} += length($aa);
        $aa =~ s/(\S{60})/$1\n/g;
        print $lumpfa_fh ">$aid\n$aa\n";
      }
      %genbank_script_proteinfa = ();
    }
  }
  close($lumpfa_fh) if $lumpfa_fh;
  print $out "##FASTA\n";
  open( IN, "$outfile.fsa" );
  while ( my $ln = <IN> ) {
    print $out $ln;
  }
  unlink("$outfile.fsa");
  close IN;
  close $out;
  return $outfile;
}

=pod 

=head2 load_sequences

We have two options: get data from genbank and create a FASTA file or have an existing FASTA file.

=cut

sub load_sequences() {
  my $ref = shift;
  warn "ERROR: No data provided for sequence loading...\n" if !$ref;
  return if !$ref;
  my $outfile;
  if ( $sequence_file && -s $sequence_file ) {
    $outfile = &_fasta2gff( $sequence_file, $sequence_metadata_file );
  } else {
    $outfile = &prepare_genbank($ref);
  }
  die "There is no sequence file to load!\n" if !$outfile || !-s $outfile;
  ######################################
  ###### CHADO load from bulk_load #####
  ######################################
  $SIG{__DIE__} = $SIG{INT} = 'cleanup_handler';
  my $adapter  = 'Bio::GMOD::DB::Adapter';
  my $ANALYSIS = 0;
  my $DBNAME   = $db_conf->name();
  my $DBUSER   = $db_conf->user();
  my $DBPASS   = $db_conf->password();
  my $DBHOST   = $db_conf->host();
  my $DBPORT   = $db_conf->port();
  my $VALIDATE    ||= 0;
  my $NOTRANSACT  ||= 0;
  my $INSERTS     ||= 0;
  my $SCORE_COL   ||= 'significance';
  my $FP_CV       ||= 'feature_property';
  my $NO_ADDFP_CV ||= 0;
  my $SKIP_VACUUM           = 1;
  my $RECREATE_CACHE        = 0;
  my $REMOVE_LOCK           = 1;
  my $ORGANISM              = $species_name;
  my $ORGANISM_FROM_CMDLINE = $ORGANISM;
  my ( $DELETE_CONFIRM, $NOEXON );
  my %argv;
  $argv{gfffile}               = $outfile;
  $argv{is_analysis}           = $ANALYSIS;
  $argv{dbprofile}             = $chado_profile;
  $argv{dbname}                = $DBNAME;
  $argv{dbuser}                = $DBUSER;
  $argv{dbpass}                = $DBPASS;
  $argv{dbhost}                = $DBHOST;
  $argv{dbport}                = $DBPORT;
  $argv{validate}              = $VALIDATE;
  $argv{notransact}            = $NOTRANSACT;
  $argv{score_col}             = $SCORE_COL;
  $argv{recreate_cache}        = $RECREATE_CACHE;
  $argv{skip_vacuum}           = $SKIP_VACUUM;
  $argv{inserts}               = $INSERTS;
  $argv{recreate_cache}        = $RECREATE_CACHE;
  $argv{fp_cv}                 = $FP_CV;
  $argv{addpropertycv}         = !$NO_ADDFP_CV;     #dgg
  $argv{allow_external_parent} = 1;
  my %locgroup  = ();
  my $chado_obj = $adapter->new(%argv);
  $chado_obj->remove_lock( force => 1 ) if $REMOVE_LOCK;
  $chado_obj->place_lock();
  my $chado_lock = 1;
  $dbh->do("TRUNCATE tmp_gff_load_cache ");
  my $select_uniquename =
    "SELECT feature_id from tmp_gff_load_cache where uniquename = ?";
  my $select_uniquename_prepared = $dbh->prepare($select_uniquename);

  #if we need custom ontology mapping, cache them here
  $chado_obj->file_handles();
  my $gffio = Bio::FeatureIO->new(
                                   -file              => $argv{gfffile},
                                   -format            => 'gff',
                                   -ignore_seq_region => 1,
                                   -validate_terms    => $VALIDATE
  );
  warn "Preparing data for inserting into the $DBNAME database\n";
  warn "(This may take a while ...)\n";
  my $seen_cds = my $seen_exon = my $seen_bad_cds = 0;
  my $ORGANISM_FROMDATA = ( $ORGANISM =~ /fromdata/ );
  if ($ORGANISM_FROMDATA) {
    $ORGANISM = "null";
  } elsif ($ORGANISM_FROM_CMDLINE) {
    $ORGANISM = $ORGANISM_FROM_CMDLINE;
    $chado_obj->organism($ORGANISM);
    $chado_obj->organism_id($ORGANISM)
      or die "ERROR: $ORGANISM organism not found in the database";
  } elsif ( defined $gffio && $gffio->organism ) {
    $ORGANISM = $gffio->organism;
    $chado_obj->organism($ORGANISM);
    $chado_obj->organism_id($ORGANISM)
      or die "ERROR: $ORGANISM organism not found in the database";
  } else {
    $chado_obj->organism($ORGANISM);
    $chado_obj->organism_id($ORGANISM)
      or die "ERROR: $ORGANISM organism not found in the database";
  }
  my ( $analysis_err_str, $cds_err_str );
  my $processing_CDS = 0;
  my $feature_iterator;
  my $itern = 0;
  my %seen_organism;
FEATURE:
  while ( my $feature =
             ( defined $feature_iterator && $feature_iterator->next_feature )
          || ( defined $gffio && $gffio->next_feature ) )
  {
    $chado_obj->primary_dbxref('');
    my $featuretype  = $feature->type->name;
    my $gff_organism = ( $feature->get_tag_values('organism') )[0];
    $gff_organism =
      defined( $feature->annotation->get_Annotations('organism') )
      ? ( $feature->annotation->get_Annotations('organism') )[0]
      : ''
      if !$gff_organism;
    if ( !$chado_obj->organism()
         || ( $gff_organism && $gff_organism ne $chado_obj->organism() ) )
    {
      warn "Organism $gff_organism from data\n"
        unless ( $seen_organism{$ORGANISM}++ );
      $chado_obj->organism($gff_organism);
      $chado_obj->organism_id($gff_organism)
        or die "ERROR: $ORGANISM organism not found in the database";
    }
    $itern++;
    $seen_exon = 1 if $featuretype =~ /exon$/ and !$processing_CDS;
    if ( $featuretype =~ /(CDS|UTR)/ ) {
      my $continue_on = $chado_obj->handle_CDS($feature);
      $seen_cds = 1 if !$seen_cds && $featuretype =~ /CDS/;
      next FEATURE unless ( $continue_on == 0 );
      if ( !$seen_bad_cds ) {
        warn <<END;

This GFF file has CDS and/or UTR features that do not belong to a 
'central dogma' gene (ie, gene/transcript/CDS).  The features of 
this type are being stored in the database as is.

END
        $seen_bad_cds = 1;
      }
    }
    if (   !$cds_err_str
         && $seen_cds
         && $seen_exon
         && !$NOEXON
         && !$processing_CDS )
    {
      $cds_err_str =
          "\n\nThere are both CDS and exon features in this file, but\n"
        . "you did not set the --noexon option, which you probably want.\n"
        . "Please see `perldoc gmod_bulk_load_gff3.pl` for more information.\n\n";
      warn $cds_err_str;
    }
    my $nextfeature    = $chado_obj->nextfeature();
    my $nextfeatureloc = $chado_obj->nextfeatureloc();
    my $type           = $chado_obj->get_type($featuretype);
    my ( $src, $seqlen ) = $chado_obj->get_src_seqlen($feature);
    if ( !$src ) {
      $src = $chado_obj->src_second_chance($feature);
    }
    die "ERROR: no feature for " . $feature->seq_id unless $src;
    if ( $feature->annotation->get_Annotations('Parent') ) {
      $chado_obj->handle_parent($feature);
    }
    if ( $feature->annotation->get_Annotations('Derives_from') ) {
      $chado_obj->handle_derives_from($feature);
    }
    my $source =
      defined( $feature->source )
      ? $feature->source->value
      : '.';
    my ($uniquename) =
      defined( ( $feature->annotation->get_Annotations('ID') )[0] )
      ? ( $feature->annotation->get_Annotations('ID') )[0]->value
      : "auto" . $nextfeature;
    $uniquename = $uniquename->value if ref($uniquename);
    $select_uniquename_prepared->execute($uniquename);
    my $exists = $select_uniquename_prepared->fetchrow_array();
    if ($exists) {
      warn "Sequence $uniquename already exists in database. Skipping...\n";
      next;
    }
    if ( defined( ( $feature->annotation->get_Annotations('ID') )[0] )
         && ( $feature->annotation->get_Annotations('ID') )[0]->value ne
         $uniquename )
    {

      #need to keep a temporary map of modified uniquenames
      $chado_obj->modified_uniquename(
           orig_id => ( $feature->annotation->get_Annotations('ID') )[0]->value,
           modified_id => $uniquename,
           organism_id => $chado_obj->organism_id
      );
    }
    my ($name) =
      defined( ( $feature->annotation->get_Annotations('Name') )[0] )
      ? ( $feature->annotation->get_Annotations('Name') )[0]->value
      : defined( ( $feature->annotation->get_Annotations('ID') )[0] )
      ? ( $feature->annotation->get_Annotations('ID') )[0]->value
      : "$featuretype-$uniquename";
    $name = $name->value if ref($name);
    if (
         $uniquename eq $feature->seq_id
         or (
              defined(
                       $chado_obj->modified_uniquename(
                                          modified_id => $uniquename,
                                          organism_id => $chado_obj->organism_id
                       )
              )
              and $chado_obj->modified_uniquename(
                                          modified_id => $uniquename,
                                          organism_id => $chado_obj->organism_id
              ) eq $feature->seq_id
         )
      )
    {
      $chado_obj->reftype_property( $featuretype, $type )
        ;    #  a reference sequence?? yes
    }
    my $current_feature_id = 0;
    if ( $chado_obj->cache( 'feature', $uniquename ) ) {

      #seen this feature before
      $locgroup{$uniquename}++;
    } else {
      $chado_obj->cache( 'feature', $uniquename, $nextfeature );
      $locgroup{$uniquename} = 0;
      $current_feature_id = $nextfeature;
    }

    #if there are Targets, match types or scores and this is not ANALYSIS,
    #there is a problem.
    #
    if (
           !$analysis_err_str
        && !$ANALYSIS
        && (
          (
            ( scalar( $feature->annotation->get_Annotations('Target') ) > 0 )
            and (
              (
                ( $feature->annotation->get_Annotations('Target') )[0]
                ->can('value')
                && ( $feature->annotation->get_Annotations('Target') )[0]->value
              )
              or ( ( $feature->annotation->get_Annotations('Target') )[0]
                   ->can('display_text')
                   && ( $feature->annotation->get_Annotations('Target') )[0]
                   ->display_text )
            )
          )
          or ( defined( $feature->score ) and $feature->score ne '.' )
          or $featuretype =~ /match/
        )
      )
    {
      my @errs;
      push @errs, '* Has Target attributes'
        if ( scalar( $feature->annotation->get_Annotations('Target') ) > 0
          and ( $feature->annotation->get_Annotations('Target') )[0]->as_text );
      push @errs, '* Has scores'
        if ( defined( $feature->score ) and $feature->score ne '.' );
      push @errs, '* Has a match feature type'
        if ( $featuretype =~ /match/ );
      $analysis_err_str = join( "\n", @errs );
      warn
"\nThis file was not declared as analysis results (with the --analysis flag,\nbut this file contains attributes that imply that it is:\n$analysis_err_str\nYou can safely ignore this message if you don't need to access scores\nassociated with these features.\n\n";
    }
    if (    $ANALYSIS
         && $featuretype =~ /match/
         && !defined( $feature->annotation->get_Annotations('Target') ) )
    {
      if ( ( $feature->annotation->get_Annotations('ID') )[0]->can('value') ) {
        $chado_obj->cache( 'feature',
                      ( $feature->annotation->get_Annotations('ID') )[0]->value,
                      $nextfeature );
      } elsif (
         ( $feature->annotation->get_Annotations('ID') )[0]->can('display_text')
        )
      {
        $chado_obj->cache( 'feature',
               ( $feature->annotation->get_Annotations('ID') )[0]->display_text,
               $nextfeature );
      }
    }

    #don't write a featureloc entry for srcfeatures
    unless ( $src eq '\N' or $src == $nextfeature ) {

      #need to convert from base to interbase coords
      my $start = $feature->start eq '.' ? '\N' : ( $feature->start - 1 );
      my $end =
          $feature->end eq '.' ? '\N'
        : defined( $feature->end ) ? $feature->end
        :                            '\N';
      my $phase =
        ( $feature->phase eq '.' or $feature->phase eq '' )
        ? '\N'
        : $feature->phase;
      $chado_obj->print_floc(
                              $nextfeatureloc,
                              $chado_obj->cache( 'feature', $uniquename ),
                              $src,
                              $start,
                              $end,
                              $feature->strand,
                              $phase,
                              '0',
                              $locgroup{$uniquename}
      );
    }
    if ( $feature->annotation->get_Annotations('Gap') ) {
      $chado_obj->handle_gap( $feature, $uniquename );
    }
    if ( $feature->annotation->get_Annotations('Note') ) {
      $chado_obj->handle_note( $feature, $uniquename );
    }

    #try to put unreserved tags in featureprop
    #this requires that the terms exist in cvterm (and therefore that they
    #have a dbxref)
    my @unreserved_tags =
      grep { /^[a-z]/ } $feature->annotation->get_all_annotation_keys();
    if ( @unreserved_tags > 0 ) {
      $chado_obj->handle_unreserved_tags( $feature, $uniquename,
                                          @unreserved_tags );
    }
    if ( $chado_obj->{const}{source_success} && $source && $source ne '.' ) {
      $chado_obj->handle_source( $feature, $uniquename, $source );
    }
    if ( $feature->annotation->get_Annotations('Ontology_term') ) {
      $chado_obj->handle_ontology_term( $feature, $uniquename );
    }
    if (    $feature->annotation->get_Annotations('Dbxref')
         or $feature->annotation->get_Annotations('dbxref') )
    {
      $chado_obj->handle_dbxref( $feature, $uniquename );
    }
    my @aliases;
    if ( $feature->annotation->get_Annotations('Alias') ) {
      @aliases =
        map { $_->value } $feature->annotation->get_Annotations('Alias');
    }

    #if the uniquename was modified, put the orig ID in the alias list
    push @aliases,
      $chado_obj->modified_uniquename( modified_id => $uniquename,
                                       organism_id => $chado_obj->organism_id )
      if $chado_obj->modified_uniquename( modified_id => $uniquename,
                                          organism_id => $chado_obj->organism_id
      );
    my %count;
    my @ualiases = grep { ++$count{$_} < 2 } @aliases;
    foreach my $alias (@ualiases) {
      $chado_obj->synonyms( $alias,
                            $chado_obj->cache( 'feature', $uniquename ) );
    }
    if (    $current_feature_id
         or $chado_obj->cache( 'srcfeature', $nextfeature ) )
    {
      $chado_obj->print_f( $nextfeature, $chado_obj->organism_id, $name,
                           $uniquename, $type, $seqlen,
                           $chado_obj->primary_dbxref );
    }
    if ( $ANALYSIS
         && !defined( ( $feature->annotation->get_Annotations('Target') )[0] ) )
    {
      $chado_obj->handle_nontarget_analysis( $feature, $uniquename );
    }
    $chado_obj->nextfeatureloc('++');

    #now deal with creating another feature for targets
    if ( defined( ( $feature->annotation->get_Annotations('Target') )[0] ) ) {
      $chado_obj->handle_target( $feature, $uniquename, $name, $featuretype,
                                 $type );
    }
    $chado_obj->nextfeature('++');
  }
  if ( $feature_iterator = $chado_obj->process_CDS() ) {
    $processing_CDS = 1;
    goto FEATURE;    #?
  }
  $chado_obj->end_files();
  my $fh = $chado_obj->file_handles('sequence');
  while ( my $seq = $gffio->next_seq ) {
    print $fh "UPDATE feature set residues='"
      . $seq->seq()
      . "', seqlen="
      . length( $seq->seq() )
      . " WHERE feature_id=(SELECT feature_id from feature where uniquename = '"
      . $seq->display_id() . "');\n";
  }
  $chado_obj->flush_caches();
  $chado_obj->load_data();
  $chado_obj->remove_lock();
}

sub cleanup_handler {
  warn "@_\nAbnormal termination, trying to clean up...\n\n"
    if @_;    #gets the message that the die signal sent if there is one
  if ( $chado_obj && $chado_obj->dbh->ping ) {
    $chado_obj->cleanup_tmp_table;
    if ($chado_lock) {
      warn
"Trying to remove the run lock (so that --remove_lock won't be needed)...\n";
      $chado_obj->remove_lock;    #remove the lock only if we've set it
    }
    print STDERR "Exiting...\n";
  }
  exit(1);
}
##########################################################
##### SUBROUTINES FOR prepare_genbank FOLLOW #############
##########################################################
sub _typeorder {
  return 1 if ( $_[0] =~ /gene/ );
  return 2 if ( $_[0] =~ /RNA|transcript/ );
  return 3 if ( $_[0] =~ /protein|peptide/ );
  return 4 if ( $_[0] =~ /exon|CDS/ );
  return 3;
}

sub _sort_by_feattype {
  my ( $at, $bt ) = ( $a->primary_tag, $b->primary_tag );
  return ( _typeorder($at) <=> _typeorder($bt) )
    or ( $at cmp $bt );
}

sub _print_held {
  my ( $out, $gffio, $to_print, @GFF_LINE_FEAT ) = @_;
  return unless (@$to_print);
  @$to_print = sort _sort_by_feattype @$to_print;
  while ( my $feature = shift @$to_print ) {
    my $gff = $gffio->gff_string($feature);
    push @GFF_LINE_FEAT, $feature;
  }
  return ();
}

sub _maptags2gff {
  my $f       = shift;
  my $TAG_MAP = shift;
  foreach my $tag ( keys %{$TAG_MAP} ) {
    if ( $f->has_tag($tag) ) {
      my $newtag = $TAG_MAP->{$tag};
      my @v      = $f->get_tag_values($tag);
      $f->remove_tag($tag);
      $f->add_tag_value( $newtag, @v );
      if ( $tag eq 'note' ) {
        map { s/\[goid (\d+)/\[goid GO:$1/g; } @v;
        my @go = map { m/(GO:\d+)/g } @v;
        $f->add_tag_value( 'Ontology_term', @go ) if (@go);
      }
    }
  }
}

sub _getSourceInfo {
  my ( $seq, $source_type, $sf ) = @_;
  my $PROTEIN_TYPE = 'polypeptide';
  my $SOURCEID     = "GenBank";
  my $is_swiss     = ( $SOURCEID =~ /UniProt|swiss|trembl/i );
  my $is_gene      = ( $SOURCEID =~ /entrezgene/i );
  my $is_rich      = ( ref($seq) =~ /RichSeq/ );
  my $seq_name     = $seq->accession_number();
  unless ($sf) {
    $source_type =
        $is_swiss
      ? $PROTEIN_TYPE
      : $is_gene ? "eneg"    # "gene"  # "region" #
      :            $is_rich ? $seq->molecule : $source_type;
    $sf = Bio::SeqFeature::Generic->direct_new();
    my $len = $seq->length();
    $len = 1 if ( $len < 1 );
    my $start = 1;           ##$start= $len if ($len<1);
    my $loc =
        $seq->can('location')
      ? $seq->location()
      : new Bio::Location::Simple( -start => $start, -end => $len );
    $sf->location($loc);
    $sf->primary_tag($source_type);
    $sf->source_tag($SOURCEID);
    $sf->seq_id($seq_name);
    $sf->add_tag_value( Alias => $seq->id() );    # unless id == accession
    $seq->add_SeqFeature($sf);
  }
  if ( $sf->has_tag("chromosome") ) {
    $source_type = "chromosome";
    my ($chrname) = $sf->get_tag_values("chromosome");
    $sf->add_tag_value( Alias => $chrname );
  }
  my %AnnotTagMap = (
                      '_gene_name'    => 'Alias',
                      'ALIAS_SYMBOL'  => 'Alias',                   # Entrezgene
                      'LOCUS_SYNONYM' => 'Alias',                   #?
                      'symbol'        => 'Alias',
                      'synonym'       => 'Alias',
                      'dblink'        => 'Dbxref',
                      'product'       => 'product',
                      'Reference'     => 'reference',
                      'OntologyTerm'  => 'Ontology_term',
                      'comment'       => 'Note',
                      'comment1'      => 'Note',
  );
  my ($desc) = $seq->annotation->get_Annotations("desc")
    || ( $seq->desc() );
  my ($date) =
       $seq->annotation->get_Annotations("dates")
    || $seq->annotation->get_Annotations("update-date")
    || $is_rich ? $seq->get_dates() : ();
  my ($comment) = $seq->annotation->get_Annotations("comment");
  my ($species) = $seq->annotation->get_Annotations("species");
  if (   !$species
       && $seq->can('species')
       && defined $seq->species()
       && $seq->species()->can('binomial') )
  {
    $species = $seq->species()->binomial();
  }
  $sf->add_tag_value( ID => $seq_name ) unless $sf->has_tag('ID');
  $sf->add_tag_value( Note => $desc ) if ( $desc && !$sf->has_tag('Note') );
  $sf->add_tag_value( organism => $species )
    if ( $species && !$sf->has_tag('organism') );
  $sf->add_tag_value( comment1 => $comment )
    if ( !$is_swiss && $comment && !$sf->has_tag('comment1') );
  $sf->add_tag_value( date => $date ) if ( $date && !$sf->has_tag('date') );
  $sf->add_tag_value( Dbxref => $SOURCEID . ':' . $seq_name );
  foreach my $atag ( sort keys %AnnotTagMap ) {
    my $gtag = $AnnotTagMap{$atag};
    next unless ($gtag);
    my @anno = map {
      if ( ref $_ && $_->can('get_all_values') )
      {
        split( /[,;] */, join ";", $_->get_all_values );
      } elsif ( ref $_ && $_->can('display_text') ) {
        split( /[,;] */, $_->display_text );
      } elsif ( ref $_ && $_->can('value') ) {
        split( /[,;] */, $_->value );
      } else {
        ();
      }
    } $seq->annotation->get_Annotations($atag);
    foreach (@anno) { $sf->add_tag_value( $gtag => $_ ); }
  }
  return ( $source_type, $sf ) if ($sf);
  return $source_type;
}

sub _add_generic_id {
  my ( $f, $ft_name, $flags ) = @_;
  my $method = $f->primary_tag;
  $genbank_script_method{$method}++
    unless ( $flags =~ /nocount/ );
  if ( $f->has_tag('ID') ) {
  } elsif ( $f->has_tag($method) ) {
    my ($name) = $f->get_tag_values($method);
    $f->add_tag_value( ID => "$method:$name" );
  } elsif ($ft_name) {
    $f->add_tag_value( ID => $ft_name );
  } else {
    $genbank_script_idh->generate_unique_persistent_id($f);
  }
  _move_translation_fasta( $f, ( $f->get_tag_values("ID") )[0] )
    if ( $method =~ /CDS/ );
}

sub _move_translation_fasta {
  my ( $f, $ft_id ) = @_;
  if ( $ft_id && $f->has_tag('translation') ) {
    my ($aa) = $f->get_tag_values("translation");
    if ( $aa && $aa !~ /^length/ ) {
      $genbank_script_proteinfa{$ft_id} = $aa;
      $f->remove_tag("translation");
      $f->add_tag_value( "translation", "length." . length($aa) );
    }
  }
}

sub _unflatten_seq {
  my $seq          = shift;
  my $unflattener  = shift;
  my $FTSOmap      = shift;
  my $FTSOsynonyms = shift;
  my $tm           = shift;
  my $uh_oh =
      "Possible gene unflattening error with"
    . $seq->accession_number
    . ": consult STDERR\n";
  eval {
    $unflattener->_unflatten_seq(
                                  -seq       => $seq,
                                  -noinfer   => 0,
                                  -use_magic => 1
    );
  };

  if ($@) {
    warn $seq->accession_number . " Unflattening error:\n";
    warn "Details: $@\n";
    print "# " . $uh_oh;
  }
  return 0 if !$seq || !$seq->all_SeqFeatures;
  _map_types(
              $tm,
              -seq       => $seq,
              -type_map  => $FTSOmap,
              -syn_map   => $FTSOsynonyms,
              -undefined => "region"
  );
}

sub _filter {
  my $seq = shift;
  my @feats;
  my @sources;
  for my $f ( $seq->get_SeqFeatures ) {
    my $m = $f->primary_tag;
    push @sources, $f if ( $m eq 'source' );
  }
  return @sources;
}

sub _get_all_SeqFeatures {
  my $seq = shift;
  my @flatarr;
  foreach my $feat ( $seq->get_SeqFeatures ) {
    push( @flatarr, $feat );
    _add_flattened_SeqFeatures( \@flatarr, $feat );
  }
  return @flatarr;
}

sub _gene_name {
  my $g       = shift;
  my $gene_id = '';
  if ( $g->has_tag('locus_tag') ) {
    ($gene_id) = $g->get_tag_values('locus_tag');
  }
  if ( !$gene_id && $g->has_tag('gene') ) {
    ($gene_id) = $g->get_tag_values('gene');
  }
  if ( !$gene_id && $g->has_tag('ID') ) {    # for non-Genbank > Entrezgene
    ($gene_id) = $g->get_tag_values('ID');
  }
  if ( !$gene_id && $g->has_tag('product') ) {
    my ($name) = $g->get_tag_values('product');
    ($gene_id) = $name unless ( $name =~ / / );    # a description not name
  }
  if ( !$gene_id && $g->has_tag('protein_id') ) {
    my ($name) = $g->get_tag_values('protein_id');
    ($gene_id) = $name;
  }
  if ( !$gene_id && $g->has_tag('transposon') ) {
    my ($name) = $g->get_tag_values('transposon');
    ($gene_id) = $name unless ( $name =~ / / );    # a description not name
  }
  return $gene_id;
}

sub _convert_to_name {
  my $g       = shift;
  my $gene_id = '';                                # zero it;
  if ( $g->has_tag('gene') ) {
    ($gene_id) = $g->get_tag_values('gene');
    $g->remove_tag('gene');
    $g->add_tag_value( 'Name', $gene_id ) if $gene_id;
  }
  if ( !$gene_id && $g->has_tag('protein_id') ) {
    ($gene_id) = $g->get_tag_values('protein_id');
    $g->add_tag_value( 'Name', $gene_id ) if $gene_id;
  }
  if ( !$gene_id && $g->has_tag('locus_tag') ) {
    ($gene_id) = $g->get_tag_values('locus_tag');
    $g->remove_tag('locus_tag');
    $g->add_tag_value( 'Name', $gene_id ) if $gene_id;
  }
  if ( !$gene_id && $g->has_tag('product') ) {
    my ($name) = $g->get_tag_values('product');
    ($gene_id) = $name unless ( $name =~ / / );    # a description not name
    $g->add_tag_value( 'Name', $gene_id ) if $gene_id;
  }
  if ( !$gene_id && $g->has_tag('transposon') ) {
    my ($name) = $g->get_tag_values('transposon');
    ($gene_id) = $name unless ( $name =~ / / );    # a description not name
    $g->add_tag_value( 'Name', $gene_id ) if $gene_id;
  }
  if ( !$gene_id && $g->has_tag('ID') ) {
    my ($name) = $g->get_tag_values('ID');
    $g->add_tag_value( 'Name', $name );
  }
  return $gene_id;
}

sub _add_flattened_SeqFeatures {
  my ( $arrayref, $feat ) = @_;
  my @subs = ();
  if ( $feat->isa("Bio::FeatureHolderI") ) {
    @subs = $feat->get_SeqFeatures;
  } elsif ( $feat->isa("Bio::SeqFeatureI") ) {
    @subs = $feat->sub_SeqFeature;
  } else {
    warn ref($feat)
      . " is neither a FeatureHolderI nor a SeqFeatureI. "
      . "Don't know how to flatten.";
  }
  for my $sub (@subs) {
    push( @$arrayref, $sub );
    _add_flattened_SeqFeatures( $arrayref, $sub );
  }
}

sub _map_types {
  my ( $self, @args ) = @_;
  my ( $sf, $seq, $type_map, $syn_map, $undefmap ) = $self->_rearrange(
    [
      qw(FEATURE
        SEQ
        TYPE_MAP
        SYN_MAP
        UNDEFINED
        )
    ],
    @args
  );
  if ( !$sf && !$seq ) {
    $self->throw("you need to pass in either -feature or -seq");
  }
  my @sfs = ($sf);
  if ($seq) {
    $seq->isa("Bio::SeqI") || $self->throw("$seq NOT A SeqI");
    @sfs = $seq->_get_all_SeqFeatures;
  }
  $type_map = $type_map || $self->typemap;    # dgg: was type_map;
  foreach my $feat (@sfs) {
    $feat->isa("Bio::SeqFeatureI")
      || $self->throw("$feat NOT A SeqFeatureI");
    $feat->isa("Bio::FeatureHolderI")
      || $self->throw("$feat NOT A FeatureHolderI");
    my $primary_tag = $feat->primary_tag;
    my $mtype       = $type_map->{$primary_tag};
    if ($mtype) {
      if ( ref($mtype) ) {
        if ( ref($mtype) eq 'ARRAY' ) {
          my $soID;
          ( $mtype, $soID ) = @$mtype;
          if ( $undefmap && $mtype eq 'undefined' ) {    # dgg
            $mtype = $undefmap;
          }
          $type_map->{$primary_tag} = $mtype if $mtype;
        } elsif ( ref($mtype) eq 'CODE' ) {
          $mtype = $mtype->($feat);
        } else {
          $self->throw('must be scalar or CODE ref');
        }
      } elsif ( $undefmap && $mtype eq 'undefined' ) {    # dgg
        $mtype = $undefmap;
      }
      $feat->primary_tag($mtype);
    }
    if ( !$mtype && $syn_map ) {
      if ( $feat->has_tag('note') ) {
        my @all_matches;
        my @note = $feat->each_tag_value('note');
        for my $k ( keys %$syn_map ) {
          if ( $k =~ /"(.+)"/ ) {
            my $syn = $1;
            for my $note (@note) {
              if ( $syn eq $note ) {
                my @map        = @{ $syn_map->{$k} };
                my $best_guess = $map[0];
                unshift @{ $all_matches[-1] }, [$best_guess];
                $mtype = $best_guess;
                print '#' x 78
                  . "\nGuessing the proper SO term for GenBank"
                  . " entry:\n\n"
                  . _GenBank_entry($feat)
                  . "\nis:\t$mtype\n"
                  . '#' x 78 . "\n\n";
              } else {
                _SO_fuzzy_match( $k, $primary_tag, $note, $syn, \@all_matches );
              }
            }
          }
        }

        #unless ($mtype) {
        for my $note (@note) {
          for my $name ( values %$type_map ) {
            _SO_fuzzy_match( $name, $primary_tag, $note, $name, \@all_matches );
          }
        }

        #}
        if ( scalar(@all_matches) && !$mtype ) {
          my $top_matches = first { defined $_ } @{ $all_matches[-1] };
          my $best_guess = $top_matches->[0];
          if ( $best_guess =~ /"(.+)"/ ) {
            $best_guess = $syn_map->{$best_guess}->[0];
          }
          @genbank_script_return = @all_matches;
          $mtype                 = $best_guess;
          print '#' x 78
            . "\nGuessing the proper SO term for GenBank"
            . " entry:\n\n"
            . _GenBank_entry($feat)
            . "\nis:\t$mtype\n"
            . '#' x 78 . "\n\n";
        }
      }
      $mtype ||= $undefmap;
      $feat->primary_tag($mtype);
    }
  }
}

sub _SO_fuzzy_match {
  my $candidate        = shift;
  my $primary_tag      = shift;
  my $note             = shift;
  my $SO_terms         = shift;
  my $best_matches_ref = shift;
  my $modifier         = shift;
  $modifier ||= '';
  my @feat_terms;

  for ( split( " |_", $primary_tag ) ) {
    my @camelCase = /(?:[A-Z]|[a-z])(?:[A-Z]+|[a-z]*)(?=$|[A-Z]|[;:.,])/g;
    push @feat_terms, @camelCase;
  }
  for ( split( " |_", $note ) ) {
    ( my $word = $_ ) =~ s/[;:.,]//g;
    push @feat_terms, $word;
  }
  my @SO_terms = split( " |_", $SO_terms );
  my ( $plus, $minus ) = ( 0, 0 );
  my %feat_terms;
  my %SO_terms;
  map { $feat_terms{$_} = 1 } @feat_terms;
  map { $SO_terms{$_}   = 1 } @SO_terms;
  for my $st ( keys %SO_terms ) {

    for my $ft ( keys %feat_terms ) {
      ( $st =~ m/$modifier\Q$ft\E/ ) ? $plus++ : $minus++;
    }
  }
  push @{ $$best_matches_ref[$plus][$minus] }, $candidate if $plus;
}

sub _GenBank_entry {
  my ( $f, $delimiter, $num ) = @_;
  $delimiter ||= "\n";
  my $entry =
      ( $num ? ' [1] ' : ' ' x 5 ) 
    . $f->primary_tag
    . (
        $num
        ? ' ' x ( 12 - length $f->primary_tag ) . ' [2] '
        : ' ' x ( 15 - length $f->primary_tag )
    )
    . $f->start . '..'
    . $f->end
    . "$delimiter";
  if ($num) {
    _words_tag( $f, \$entry );
  } else {
    for my $tag ( $f->all_tags ) {
      for my $val ( $f->each_tag_value($tag) ) {
        $entry .= ' ' x 20;
        $entry .=
          $val eq '_no_value'
          ? "/$tag$delimiter"
          : "/$tag=\"$val\"$delimiter";
      }
    }
  }
  return $entry;
}

sub _gff_validate {
  my @feat = @_;
  my ( %parent2child, %all_ids, %descendants, %reserved );
  for my $f (@feat) {
    for my $aTags ( [ 'Parent', \%parent2child ], [ 'ID', \%all_ids ] ) {
      map { push @{ $$aTags[1]->{$_} }, $f } $f->get_tag_values( $$aTags[0] )
        if $f->has_tag( $$aTags[0] );
    }

    #unless ($f->has_tag('organism')){
    #  warn Dumper $f;
    #}
  }
  _id_validate( \%all_ids, \%reserved );
}

sub _id_validate {
  my ( $hAll, $hReserved ) = @_;
  for my $id ( keys %$hAll ) {
    shift @{ $hAll->{$id} } unless $hReserved->{$id};
    for my $feat ( @{ $hAll->{$id} } ) {
      _id_uniquify( 0, $id, $feat, $hAll );
    }
  }
  for my $parentID ( keys %$hReserved ) {
    my @keys = keys %{ $hReserved->{$parentID} };
    shift @keys;
    for my $k (@keys) {
      my $parent    = $hReserved->{$parentID}{$k}{'parent'};
      my $aChildren = $hReserved->{$parentID}{$k}{'children'};
      my $value     = _id_uniquify( 0, $parentID, $parent, $hAll );
      for my $child (@$aChildren) {
        my %parents;
        map { $parents{$_}++ } $child->get_tag_values('Parent');
        $child->remove_tag('Parent');
        delete $parents{$parentID};
        $parents{$value}++;
        my @parents = keys %parents;
        $child->add_tag_value( 'Parent', @parents );
      }
    }
  }
}

sub _id_uniquify {
  my ( $count, $value, $feat, $hAll ) = @_;
  if ( !$count ) {
    $feat->add_tag_value( Alias => $value );
    $value .= ( '.' . $feat->primary_tag );
  } elsif ( $count == 1 ) {
    $value .= ".$count";
  } else {
    chop $value;
    $value .= $count;
  }
  $count++;
  if ( $hAll->{$value} ) {
    _id_uniquify( $count, $value, $feat, $hAll );
  } else {
    $feat->remove_tag('ID');
    $feat->add_tag_value( 'ID', $value );
    push @{ $hAll->{$value} }, $value;
  }
  $value;
}

sub _words_tag {
  my ( $feat, $entry ) = @_;
  use constant ALPHABET_DIVISOR => 26;
  use constant ALPHABET =>
    [qw(a b c d e f g h i j k l m n o p q r s t u v w x y z)];
  my @tags;
  @tags[ 1, 2 ] = (
                { 'all' => [ 'primary_tag', $feat->primary_tag ] },
                { 'all' => [ 'location',    $feat->start . '..' . $feat->end ] }
  );
  my $i = 3;
  foreach my $tag ( $feat->all_tags ) {
    foreach my $value ( $feat->each_tag_value($tag) ) {
      my ( $string, $tagged_string );
      my @words = split( /(?=\w+?)/, $value );
      my $pos = 0;
      foreach my $word (@words) {
        ( my $sanitized_word = $word ) =~ s/\W+?//g;
        $string .= $word;
        my $lead = int( $pos / ALPHABET_DIVISOR );
        my $lag  = $pos % ALPHABET_DIVISOR;
        my $a    = $lead ? ${ (ALPHABET) }[ $lead - 1 ] : '';
        $a .= $lag ? ${ (ALPHABET) }[$lag] : 'a';
        $tagged_string .= " ($a) $word";
        $tags[$i]{$a} = [ $tag, $sanitized_word ];
        $pos++;
      }
      $value = $tagged_string if scalar(@words) > 1;
      $$entry .= "[$i] /$tag=\"$value\"\r";
      $tags[$i]{'all'} = [ $tag, $string ];
    }
    $i++;
  }
  return @tags;
}

sub _fasta2gff() {
  my $fasta_file    = shift;
  my $metadata_file = shift;
  my %metadata;
  if ( -s $metadata_file ) {
    open( IN, $metadata_file );
    while ( my $ln = <IN> ) {
      chomp($ln);
      my @data = split( "\t", $ln );
      next unless $data[3];
      push( @{ $metadata{ $data[0] }{ $data[2] } }, $data[3] );
    }
    close IN;
    foreach my $id ( keys %metadata ) {
      my @array;
      foreach my $key ( keys %{ $metadata{$id} } ) {
        push( @array, $key . '=' . join( ',', @{ $metadata{$id}{$key} } ) );
      }
      $metadata{$id} = join( ';', @array );
    }
  }
  open( IN, $fasta_file );
  open my $out, ">$fasta_file.gff";
  print $out "##gff-version 3\n";
  my ( $id, $seq );
  while ( my $ln = <IN> ) {
    chomp($ln);
    if ( $ln =~ /^>(\S+)/ ) {
      my $newid = $1;
      $newid =~ s/\|.+//;
      if ( $id && $seq ) {
        print $out $id
          . "\tuser\tregion\t1\t"
          . length($seq)
          . "\t.\t+\t.\tID=$id";
        print $out ';' . $metadata{$id} if $metadata{$id};
        print $out "\n";
      }
      $seq = '';
      $id  = $newid;
    } else {
      $seq .= $ln;
    }
  }
  if ( $id && $seq ) {
    print $out $id . "\tuser\tregion\t1\t" . length($seq) . "\t.\t+\t.\tID=$id";
    print $out ';' . $metadata{$id} if $metadata{$id};
    print $out "\n";
  }
  close IN;
  print $out "##FASTA\n";
  open( IN, $fasta_file );
  while ( my $ln = <IN> ) {
    print $out $ln;
  }
  close IN;
  close $out;
  return "$fasta_file.gff";
}
