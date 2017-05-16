<?php
// $Id$

/**
 * @file
 * Installs InterProScan
 */
function biosoftware_bench_add_software_iprscan() {
  $dbxref_id  = gmod_dbsf_get_add_dbxref_withnames('Sanger', 'iprscan');
  $cvterm_id  = gmod_dbsf_get_add_cv_withnames('resource_type', 'software option');
  $software   = 'iprscan';
  $link_sql   = 'INSERT INTO {gmod_dbsf_software_resource} (software_id,resource_id) VALUES (' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ",(SELECT resource_id from {gmod_dbsf_resource} where uniquename='%s')" . ')';
  $check_sql  = "SELECT software_id as id from {gmod_dbsf_software} where uniquename='%s'";
  $insert_sql = "INSERT INTO {gmod_dbsf_software} (uniquename,description) VALUES ('%s','%s')";
  $check      = db_fetch_array(db_query($check_sql, 'ipr_convert'));
  if (empty($check['id'])) {
    db_query($insert_sql, 'ipr_convert', 'Formats InterProScan output files.');
  }
  $options_array = array(
    'iprscan_applications' => 'Applications to run (TMHMM and signalp may not be available due to licensing)',
    'iprscan_crc' => 'Use the CRC algorithm to determine whether results already exist in DB (useful for existing InterPro proteins)',
    'iprscan_iprlookup' => 'Iprlookup: give IPR domain names',
    'iprscan_goterms' => 'Search for GO terms if iprlookup is active',
  );
  foreach ($options_array as $option => $description) {
    biosoftware_bench_add_software_option($option, $dbxref_id, $cvterm_id, $description);

    switch ($option) {
      case 'iprscan_applications':
        // in the form, anything before a - will be title and after will become description
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            'blastprodom - ',
            'fprintscan - ',
            'hmmpfam - ',
            'hmmpir - ',
            'hmmpanther - ',
            'hmmtigr - ',
            'hmmsmart - ',
            'superfamily - ',
            'gene3d - ',
            'patternscan - ',
            'profilescan - ',
            'seg - ',
            'coils - ',
            //these need a license
            //'tmhmm - ',
            //'signalp - ',
          )
        );
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('checkboxes'));
        break;

      case 'iprscan_crc':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            '0',
          )
        );
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('checkbox'));
        break;

      case 'iprscan_iprlookup':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            '0',
          )
        );
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('checkbox'));
        break;

      case 'iprscan_goterms':
        gmod_dbsf_populate_prop('resource', $option, $cvterm_id, array(
            '0',
          )
        );
        gmod_dbsf_populate_cvterms('resource', $option, 'input type', array('checkbox'));
        break;
    }
    db_query($link_sql, 'iprscan', $option);
  }
  $suffix_cvterm_id = gmod_dbsf_get_add_cv_withnames('software_setting', 'dataset_suffix');
  $dataset_suffixes = gmod_dbsf_populate_prop('software', $software, $suffix_cvterm_id, array('ANY'));
}

