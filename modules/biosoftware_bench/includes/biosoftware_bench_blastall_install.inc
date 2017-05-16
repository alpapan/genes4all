<?php
// $Id$

/**
 * @file
 * Installs BLASTALL
 */

/**
 * Helper function to add BLAST software
 *
 * Use it as a template for adding plugins
 *
 */
function biosoftware_bench_add_software_blastall() {
  $software   = 'blastall';
  $dbxref_id  = gmod_dbsf_get_add_dbxref_withnames('NCBI', 'BLASTALL');
  $cvterm_id  = gmod_dbsf_get_add_cv_withnames('resource_type', 'software option');
  $check_sql  = "SELECT software_id as id from {gmod_dbsf_software} where uniquename='%s'";
  $insert_sql = "INSERT INTO {gmod_dbsf_software} (uniquename,description) VALUES ('%s','%s')";
  $link_sql   = 'INSERT INTO {gmod_dbsf_software_resource} (software_id,resource_id) VALUES (' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='blastall')" . ",(SELECT resource_id from {gmod_dbsf_resource} where uniquename='%s')" . ')';
  if (empty($dbxref_id) || empty($cvterm_id)) {
    drupal_set_message(t('There has been a problem getting the dbxref_id and/or cvterm_id for BLASTALL / "software option" from the database (!dbxref_id,!cvterm_id).', array('!dbxref_id' => $dbxref_id, '!cvterm_id' => $cvterm_id)));
    return FALSE;
  }
  $check = db_fetch_array(db_query($check_sql, 'formatdb'));
  if (empty($check['id'])) {
    db_query($insert_sql, 'formatdb', 'Format reference sequence as a BLAST database - NCBI C toolkit. Required for blastall to work properly.');
  }
  $check = db_fetch_array(db_query($check_sql, 'fastacmd'));
  if (empty($check['id'])) {
    db_query($insert_sql, 'fastacmd', 'Retrieve data from formatted BLAST databases - NCBI C toolkit. Required for blastall to work properly.');
  }

  $options_array = array(
    'Scoring matrix' => '',
    'Filtering' => '',
    'E-value cutoff' => '',
    'Advanced BLAST options' => '',
    'query genetic code' => '',
    'DB genetic code' => '',
    'Frame shift penalty' => '',
    'Number of results' => '',
  );
  foreach ($options_array as $option => $description) {
    biosoftware_bench_add_software_option($option, $dbxref_id, $cvterm_id, $description);
    switch ($option) {
      case 'Scoring matrix':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            'BLOSUM62', 'PAM30', 'PAM70', 'BLOSUM80', 'BLOSUM45',
          ));
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('select'));
        break;

      case 'Filtering':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            'Low complexity regions', 'Lower-case letters',
          ));
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('checkboxes'));
        break;

      case 'E-value cutoff':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            '1e-3', '1e-30', '1e-10', '1e-5', '0.1', '0', '10',
          ));
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('select'));
        break;

      case 'Advanced BLAST options':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            'ungapped', 'megablast',
          ));
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('checkboxes'));
        break;

      case 'query genetic code':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            'Standard', 'Vertebrate Mitochondrial', 'Yeast Mitochondrial', 'Mold, Protozoan, and Coelocoel Mitochondrial'
            , 'Invertebrate Mitochondrial', 'Ciliate Nuclear', 'Echinoderm Mitochondrial', 'Euplotid Nuclear'
            , 'Bacterial', 'Alternative Yeast Nuclear', 'Ascidian Mitochondrial', 'Flatworm Mitochondrial'
            , 'Blepharisma Macronuclear',
          ));
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('select'));
        break;

      case 'DB genetic code':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            'Standard', 'Vertebrate Mitochondrial', 'Yeast Mitochondrial', 'Mold, Protozoan, and Coelocoel Mitochondrial'
            , 'Invertebrate Mitochondrial', 'Ciliate Nuclear', 'Echinoderm Mitochondrial', 'Euplotid Nuclear'
            , 'Bacterial', 'Alternative Yeast Nuclear', 'Ascidian Mitochondrial', 'Flatworm Mitochondrial'
            , 'Blepharisma Macronuclear',
          ));
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('select'));
        break;

      case 'Frame shift penalty':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            '0',
          ));
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('textfield'));
        break;

      case 'Number of results':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            '10', '100', '1', '50',
          )
        );
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('select'));
        break;
    }
    //Now link them to the software
    db_query($link_sql, $option);
  }
  //Store allowed dataset suffixes. Rank 0 is the switch specifying if ANY or ALL of
  //the suffixes must be present for a dataset to be allowed
  //remember to include the "." if there is one before the suffix
  $suffix_cvterm_id = gmod_dbsf_get_add_cv_withnames('software_setting', 'dataset_suffix');
  $dataset_suffixes = gmod_dbsf_populate_prop('software', 'blastall', $suffix_cvterm_id, array('ANY', '.pin', '.nin', '.pal', '.nal'));


  return TRUE;
}

