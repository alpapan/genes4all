<?php
// $Id$

/**
 * @file
 * Installs rnaseq
 */



/**
 * Helper function to add rnaseq software
 *
 * Use it as a template for adding plugins
 *
 */
function biosoftware_bench_add_software_rnaseq() {
  $software   = 'rnaseq';
  $dbxref_id  = gmod_dbsf_get_add_dbxref_withnames('alexie@butterflybase.org', $software);
  $check_sql  = "SELECT software_id as id from {gmod_dbsf_software} where uniquename='%s'";
  $cvterm_id  = gmod_dbsf_get_add_cv_withnames('resource_type', 'software option');
  $insert_sql = "INSERT INTO {gmod_dbsf_software} (uniquename,description) VALUES ('%s','%s')";
  $link_sql   = 'INSERT INTO {gmod_dbsf_software_resource} (software_id,resource_id) VALUES (' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='$software')" . ",(SELECT resource_id from {gmod_dbsf_resource} where uniquename='%s')" . ')';
  if (empty($dbxref_id) || empty($cvterm_id)) {
    drupal_set_message(t('There has been a problem getting the dbxref_id and/or cvterm_id for %software / "software option" from the database (!dbxref_id,!cvterm_id).', array('!dbxref_id' => $dbxref_id, '!cvterm_id' => $cvterm_id,'%software'=> $software)));
    return FALSE;
  }
  $suffix_cvterm_id = gmod_dbsf_get_add_cv_withnames('software_setting', 'dataset_suffix');
  $dataset_suffixes = gmod_dbsf_populate_prop('software', $software, $suffix_cvterm_id, array('ALL', '.bam'));
  return TRUE;
}

