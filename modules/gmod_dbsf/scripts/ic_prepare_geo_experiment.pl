#!/usr/bin/perl -w

=pod

=head1 USAGE

  -data_files :s{1,}' => \@experiment_data_files,
  -longitude_column :i' =>    \$experiment_lon_column,           # column number with longitude data
  -latitude_column :i' =>    \$experiment_lat_column,           # column number with latitude data
  -habitat_columns :i{,}'  => \@experiment_habitat_columns,
  -genotype_columns :i{,}' => \@experiment_genotype_columns,
  -phenotype_columns :i{,}' => \@experiment_phenotype_columns,
  -uname_column :i'        => \$experiment_uname_column,
  -type :s'               => \$experiment_type,
  -project :s' => \$experiment_project,
  
  
=cut
use strict;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Text::CSV_XS;
my $delimiter  = "\t";
my $quote_char = '"';
my (
     @experiment_data_files,   $experiment_lon_column,
     $experiment_lat_column,   @experiment_habitat_columns,
     $experiment_uname_column, @experiment_genotype_columns,
     $experiment_project,      @experiment_phenotype_columns,@experiment_metadata_columns
);

#data needed:
my ( $ncbi_taxid, $marker_type, $experiment_type, $genotype_type, $stock_type );
GetOptions(
  'quote:s'            => \$quote_char,
  'delim:s'            => \$delimiter,
  'i|data_files:s{1,}' => \@experiment_data_files,
  'longitude_column:i' =>
    \$experiment_lon_column,    # column number with longitude data
  'latitude_column:i' =>
    \$experiment_lat_column,    # column number with latitude data
  'habitat_columns:i{,}'   => \@experiment_habitat_columns,
  'metadata_columns:i{,}'   => \@experiment_metadata_columns,
  'genotype_columns:i{,}'  => \@experiment_genotype_columns,
  'phenotype_columns:i{,}' => \@experiment_phenotype_columns,
  'uname_column:i'         => \$experiment_uname_column,
  'experiment_type:s'      => \$experiment_type,
  'project:s'              => \$experiment_project,
  'stock_type:s'           => \$stock_type,
  'marker_type:s'          => \$marker_type,
  'genotype_type:s'        => \$genotype_type,
  'organism|ncbi_taxid:i'  => \$ncbi_taxid,
);
pod2usage "No files provided\n" unless (@experiment_data_files);
pod2usage "NCBI taxid is required\n" unless $ncbi_taxid;
pod2usage "Marker type is required (e.g. microsatellite)\n" unless $marker_type;
pod2usage "Genotype is required (e.g. haplotype)\n" unless $genotype_type;
pod2usage "Valid experiment type is required (e.g. natural_collection)\n"
  unless $experiment_type;
pod2usage "Valid stock type is required (e.g. host_mass_collected_once)\n"
  unless $stock_type;
pod2usage "You must provide the name of your project\n"
  unless $experiment_project;
&preprocess_single_experiment();
#######################################################
sub preprocess_single_experiment() {
  my (%hash_names);
  my %habitat_cols   = columns2hash( \@experiment_habitat_columns );
  my %experiment_cols   = columns2hash( \@experiment_metadata_columns );
  my %genotype_cols  = columns2hash( \@experiment_genotype_columns );
  my %phenotype_cols = columns2hash( \@experiment_phenotype_columns );
  my (
       $experiment_cv,     $stock_cv,             $relationship_cv,
       $genotype_cv,       $event_cv,             $occurrence_cv,
       $identification_cv, $MeasurementOrFact_cv, $location_cv,
       $GeologicalContext_cv
  ) = &see_core_cvs();
  foreach my $file (@experiment_data_files) {
    next unless $file && -s $file;
    my $basename = basename($file);
    open my ($fh), $file;
    my $header_str = <$fh>;
    if ( $header_str =~ /\015$/ ) {
      close $fh;
      $file .= '.unix';

      #rewrite file
      $header_str =~ tr/\015/\012/;
      open( OUT, ">" . $file );
      print OUT $header_str;
      close OUT;
      open $fh, $file;
      $header_str = <$fh>;
    }
    $header_str =~ s/[\r\n]+$//g;
    my @headers = split( $delimiter, $header_str );

    # check which habitat headers, if any, are actually in the Darwin core
    # those that are we can use as is. the others have to be prefixed
    # with the project name
    my %unverified_headers;
    for ( my $i = 0 ; $i < scalar(@headers) ; $i++ ) {
      next unless $habitat_cols{$i} || $phenotype_cols{$i} || $experiment_cols{$i};
      my $head = lc( $headers[$i] );
      if ( $experiment_cv->{$head} ) {
        $head = 'experiment_type' . "\t" . $head;
      } elsif ( $stock_cv->{$head} ) {
        $head = 'stock_type' . "\t" . $head;
      } elsif ( $relationship_cv->{$head} ) {
        $head = 'relationship' . "\t" . $head;
      } elsif ( $genotype_cv->{$head} ) {
        $head = 'genotype_type' . "\t" . $head;
      } elsif ( $event_cv->{$head} ) {
        $head = 'Event' . "\t" . $head;
      } elsif ( $occurrence_cv->{$head} ) {
        $head = 'Occurrence' . "\t" . $head;
      } elsif ( $identification_cv->{$head} ) {
        $head = 'Identification' . "\t" . $head;
      } elsif ( $MeasurementOrFact_cv->{$head} ) {
        $head = 'MeasurementOrFact' . "\t" . $head;
      } elsif ( $location_cv->{$head} ) {
        $head = 'Location' . "\t" . $head;
      } elsif ( $GeologicalContext_cv->{$head} ) {
        $head = 'GeologicalContext' . "\t" . $head;
      } else {
        $unverified_headers{$head} = 1;
        $head = $experiment_project . '::habitat' . "\t" . $head;
      }
      $headers[$i] = $head;
    }
    warn "Custom CVs will be created for these headers:\n"
      . join( ", ", keys %unverified_headers ) . "\n";
    my $eol = ( $header_str =~ /\r\n$/ ) ? "\015\012" : "\012";
    my $csv = Text::CSV_XS->new(
                                 {
                                   binary           => 1,
                                   eol              => $eol,
                                   sep_char         => $delimiter,
                                   allow_whitespace => 0,
                                   empty_is_undef   => 1,
                                   quote_char       => $quote_char
                                 }
    ) or die "Cannot use CSV: " . Text::CSV_XS->error_diag();
    my ($row);
    open( COLLECT,   ">$file.collect" );
    open( GENOTYPE,  ">$file.genotypes" );
    open( PHENOTYPE, ">$file.phenotypes" );
    open( HABITAT,  ">$file.geo.metadata" );
    open( EXPERIMENT,  ">$file.exp.metadata" );
    print COLLECT "#$experiment_project|stock.uname\texperiment.type\tstock.type\tLong\tLat\taltitude\tdescription\n";
    print HABITAT "#$experiment_project|stock.uname\tCV\tcvterm\tvalue\n";
    print EXPERIMENT "#$experiment_project|stock.uname\tCV\tcvterm\tvalue\n";
    print GENOTYPE
      "#$experiment_project|stock.uname\tmarker.type\tmarker\tallele\n";
    print PHENOTYPE "#$experiment_project|stock.uname\tphenotype.uname\tallele\n";

    # create marker data
    #TODO somehow link this to a sequence file....
    # store all possible alleles by parsing the files my %alleles;
    # each allele is a genotype
    open( MARKER, ">$file.marker" );
    open( MARKERGFF, ">$file.marker.gff" );
    print MARKERGFF "##gff-version 3\n";
    print MARKER "#ID\tType\tNCBI TAXID\n";
    foreach my $col ( keys %genotype_cols ) {
      my $marker = $headers[$col] || next;
      print MARKER "$marker\t$marker_type\t$ncbi_taxid\n";
      print MARKERGFF "$marker\tmarker\t$marker_type\t1\t2\t.\t.\t.\tID=$marker\n";
    }
    print MARKERGFF "##FASTA\n";
    close MARKER;
    close MARKERGFF;
    while ( my $ln_obj = $csv->getline($fh) ) {
      $row++;
      my $column_number = scalar(@$ln_obj);
      next if ( !$column_number || $column_number < 1 );
      my $name =
          $experiment_uname_column
        ? $ln_obj->[ $experiment_uname_column - 1 ]
        : $row;
      $hash_names{$name}++;
      die "Name $name appears more than once!" if $hash_names{$name} > 1;
      my $lat  = &lat2dec( $ln_obj->[ $experiment_lat_column - 1 ] );
      my $long = &lat2dec( $ln_obj->[ $experiment_lon_column - 1 ] );
      die "Latitude is not a number:" . $lat
        unless $lat && $lat =~ /^[\d\.\-\+]+$/;
      die "Longitude is not a number:" . $long
        unless $long && $long =~ /^[\d\.\-\+]+$/;

      #print out long/lat
      print COLLECT "$name\t$experiment_type\t$stock_type\t$long\t$lat\n";

      #metadata
      foreach my $hcol ( keys %habitat_cols ) {
        my $cv = $headers[$hcol] || next;
        my $value = $ln_obj->[$hcol];
        print HABITAT "$name\t$cv\t$value\n";
      }

      #genotypes
      foreach my $col ( keys %genotype_cols ) {
        my $marker = $headers[$col] || next;
        my $allele = $ln_obj->[$col];
        print GENOTYPE "$name\t" . $genotype_type . "\t$marker\t$allele\n";
      }

      #phenotypes
      foreach my $col ( keys %phenotype_cols ) {
        my $marker = $headers[$col];
        my $allele = $ln_obj->[$col];
        print PHENOTYPE "$name\t$marker\t$allele\n";
      }

    }
    close COLLECT;
    close GENOTYPE;
    close PHENOTYPE;
    close HABITAT;
    close EXPERIMENT;
  }
}

sub columns2hash() {

  # also convert column numbering from 0
  my $ref = shift;
  return if !$ref;
  my %hash;
  foreach my $r ( sort { $a <=> $b } @$ref ) {
    $r--;
    $hash{$r} = 1;
  }
  return %hash;
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

sub see_core_cvs() {
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

#darwin core:
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
  my $event_cv             = $activity_metadata->{'Event'};
  my $occurrence_cv        = $stock_metadata->{'Occurrence'};
  my $identification_cv    = $stock_metadata->{'Identification'};
  my $MeasurementOrFact_cv = $cvterm_metadata->{'MeasurementOrFact'};
  my $location_cv          = $geolocation_metadata->{'Location'};
  my $GeologicalContext_cv = $activity_metadata->{'GeologicalContext'};
  return (
           $experiment_cv,     $stock_cv,             $relationship_cv,
           $genotype_cv,       $event_cv,             $occurrence_cv,
           $identification_cv, $MeasurementOrFact_cv, $location_cv,
           $GeologicalContext_cv
  );
}
