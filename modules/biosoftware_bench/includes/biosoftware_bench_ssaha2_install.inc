<?php
// $Id$

/**
 * @file
 * Installs SSAHA2
 */
function biosoftware_bench_add_software_ssaha2() {
  // get cvterm for software option
  $cvterm_id = gmod_dbsf_get_add_cv_withnames('resource_type', 'software option');
  $link_sql = 'INSERT INTO {gmod_dbsf_software_resource} (software_id,resource_id) VALUES (' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ",(SELECT resource_id from {gmod_dbsf_resource} where uniquename='%s')" . ')';


  /*PLUGIN_EDIT_HERE
   * Add CVterms which might be needed
   * $dbxref_id -> who made this option/can provide help? For SSAHA2 it's the Sanger institute
   * $options_array contains option names and descriptions. Prefix name with $software.
   */

  $check_sql  = "SELECT software_id as id from {gmod_dbsf_software} where uniquename='%s'";
  $insert_sql = "INSERT INTO {gmod_dbsf_software} (uniquename,description) VALUES ('%s','%s')";
  $software   = 'ssaha2';
  gmod_dbsf_add_cv('ssaha_hash_type', array('sanger', 'solexa', '454'));
  // Name and description
  $check = db_fetch_array(db_query($check_sql, 'ssaha2Build'));
  if (empty($check['id'])) {
    db_query($insert_sql, 'ssaha2Build', 'Format reference sequences as a SSAHA2 database. Needed for creating your own SSAHA2 databases.');
  }
  $dbxref_id = gmod_dbsf_get_add_dbxref_withnames('Sanger', 'SSAHA2');
  // Name and descriptions appearing on user

  $options_array = array(
    //  'ssaha_seeds' => 'Number of kmer matches required to flag a hit',
    'ssaha2_query_type' => 'Optimize for query size/sequencing technology',
    'ssaha2_best' => 'Only report best alignment per query',
    'ssaha2_score' => 'Minimum score for match to be reported',
    'ssaha2_identity' => 'Minimum identity for match to be reported',
    'ssaha2_cmatch' => 'Minimum match length for cross_match matching',
    // lower than 13
    'ssaha2_ckmer' => 'Word size for cross_match matching',
    // max 100 billion
    'ssaha2_cut' => 'Number of repeats allowed before this kmer is ignored',
    'ssaha2_depth' => 'Number of hits to consider for alignment',
  );


  foreach ($options_array as $option => $description) {
    biosoftware_bench_add_software_option($option, $dbxref_id, $cvterm_id, $description);

    switch ($option) {
      /* PLUGIN_EDIT_HERE
       * define values for each option. first appears as the default
       */

      case 'ssaha2_score':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            '30', '5', '12', '15', '18', '20', '25', '30', '50', '100', '200', '500', '1000', '1500', '2000',
          )
        );
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('select'));
        break;

      case 'ssaha2_identity':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            '50', '20', '30', '40', '45', '55', '60', '65', '70', '75', '80', '85', '90', '95', '98', '99', '100',
          )
        );
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('select'));
        break;

      case 'ssaha2_cmatch':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            '14', '5', '9', '16', '18', '20', '22', '24', '28', '32', '40', '50', '100',
          )
        );
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('select'));
        break;

      case 'ssaha2_ckmer':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            '10', '4', '5', '6', '8', '11', '12',
          )
        );
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('select'));
        break;

      case 'ssaha2_cut':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            '10000', '3000', '5000', '8000', '15000', '20000', '40000', '80000', '100000', '500000', '1000000',
          )
        );
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('select'));
        break;

      case 'ssaha2_depth':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            '50',
          )
        );
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('textfield'));
        break;

      case 'ssaha2_best':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            '0',
          )
        );
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('checkbox'));
        break;

      case 'ssaha2_query_type':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            'Sanger ABI', '454', 'Illumina-Solexa',
          )
        );
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('radios'));
        break;
    }
    db_query($link_sql, $software, $option);
  }
  /*PLUGIN_EDIT_HERE
   * Store allowed dataset suffixes. Rank 0 is the switch specifying if ANY or ALL of the suffixes must be present for a dataset to be allowed
   * If you don't want to use suffixes but still want to link datasets, then give simply array('ANY')
   * Some software e.g. see iprscan, do not need any datasets
   * remember to include the "." if there is one before the suffix
   */


  //index is for the indexed file using index_fasta.pl (in script directory, a bioperl script)
  // empty string '' is for requiring the base (FASTA) file. index_fasta.pl uses it to get the index sequence
  $suffix_cvterm_id = gmod_dbsf_get_add_cv_withnames('software_setting', 'dataset_suffix');
  $dataset_suffixes = gmod_dbsf_populate_prop('software', $software, $suffix_cvterm_id,
    array('ALL', '', '.index', '_454.head.gz', '_454.body', '_454.base',
      '_454.name', '_454.size', '_sanger.head.gz', '_sanger.body', '_sanger.base',
      '_sanger.name', '_sanger.size', '_solexa.head.gz', '_solexa.body', '_solexa.base',
      '_solexa.name', '_solexa.size',
    )
  );

  return TRUE;
}

