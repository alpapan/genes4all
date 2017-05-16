<?php
// $Id$

/**
 * @file
 * Installs annot8r
 */
function biosoftware_bench_add_software_annot8r() {
  gmod_dbsf_add_cv('software_setting', array('executable', 'executable physprop', 'executable annot8r', 'executable annot8r2gff', 'credential port', 'credential host', 'credential user', 'credential pass', 'GOdb', 'ECdb', 'KEGGdb'));
  $dbxref_id = gmod_dbsf_get_add_dbxref_withnames('Edinburgh', 'annot8r');
  $cvterm_id = gmod_dbsf_get_add_cv_withnames('resource_type', 'software option');
  //$check_sql= "SELECT software_id as id from {gmod_dbsf_software} where uniquename='%s'";
  //$insert_sql="INSERT INTO {gmod_dbsf_software} (uniquename,description) VALUES ('%s','%s')";
  $link_sql = 'INSERT INTO {gmod_dbsf_software_resource} (software_id,resource_id) VALUES (' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ",(SELECT resource_id from {gmod_dbsf_resource} where uniquename='%s')" . ')';
  $options_array = array(
    'annot8r_applications' => 'Applications to run',
  );
  $software = 'annot8r';
  foreach ($options_array as $option => $description) {
    biosoftware_bench_add_software_option($option, $dbxref_id, $cvterm_id, $description);

    switch ($option) {
      case 'annot8r_applications':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            'GO', 'KEGG', 'EC',
          )
        );
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('checkboxes'));
        break;
    }
    db_query($link_sql, 'annot8r', $option);
  }
  $suffix_cvterm_id = gmod_dbsf_get_add_cv_withnames('software_setting', 'dataset_suffix');
  $dataset_suffixes = gmod_dbsf_populate_prop('software', $software, $suffix_cvterm_id, array('ANY'));
}

