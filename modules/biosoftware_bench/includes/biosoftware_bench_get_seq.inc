<?php
// $Id$

/**
 * @file
 * Fetch sequences from various programs using a indexed database
 */

/**
 * Retrieve the sequence from a formatdb-ed dataset
 *
 * TODO: support other programs? Or use another method for indexing
 * and storing data?
 * Used when retrieving hits of BLAST reports
 *
 * @param $uid
 *   the UID of the BLAST report, used to grab the database
 * @param $database
 *   Explicitly set the database to get the sequence from.
 *   It can be a string or a list array.
 *
 * @return
 *   HTML text of sequences in FASTA format plus a link to download a text file
 */
function biosoftware_bench_get_sequence($uid = NULL, $database = NULL) {
  if (empty($uid) && empty($database)) {
    return 'No method to acquire the database was provided';
  }
  if ($database === $uid) {
    $database = '';
  }
  if ($database == 'iprscan') {
    return 'InterProScan is currently not supported, please click on the graph to see the HTML report.';
  }
  $data = $_POST;

  // to avoid duplicates
  $seqs_array = array();
  foreach ($data as $key => $value) {
    if (is_numeric($key) && !is_array($value)) {
      $seqs_array[$value] = escapeshellarg(trim($value));
    }
  }
  if (empty($seqs_array)) {
    return 'No sequences were requested for retrieval';
  }
  $seqs_to_get = gmod_dbsf_array2csv($seqs_array);
  // currently only FASTA
  //$operation = isset($_POST['op']) ? $_POST['op'] : 'FASTA';
  //Alter below when we fetch the file.
  $software = drupal_get_path('module', 'biosoftware_bench') .'/scripts/index_fasta.pl';
  $database_array = array();
  //first get database using $uid.
  if (!empty($uid)) {
    // get databases with which to search for ID.
    $dirpath = file_create_path() .'/bench/';
    $f = file_scan_directory($dirpath, "^$uid.+done$", array('.', '..', 'CVS'), NULL, FALSE, 'filename');
    foreach ($f as $infile => $obj) {
      //PLUGIN_EDIT_HERE -> if want another program not fastacmd to act for fetching.
      if (strpos($infile, 'blastall') !== FALSE) {
        $software = 'fastacmd';
      }
      if ($inhandle = fopen($infile, 'rb')) {
        while (!feof($inhandle)) {
          $line = trim(fgets($inhandle));
          if (!empty($line) && preg_match('/^database:/', $line)) {
            $line_data = explode(':', $line);
            //set as key so that we avoid duplicates
            $database_array[trim(escapeshellarg($line_data[1]))] = 1;
          }
        }
      }
      else {
        drupal_set_message('There is a problem reading a file. Perhaps a permissions error. Please '. l('contact', 'contact') . t(' the administrator citing %t', array('%t' => $uid .' '. $database)), 'error');
        return FALSE;
      }
    }
  }
  // We can have both a UID and a database search at the same time.
  if (!empty($database)) {
    if (is_array($database)) {
      foreach ($database as $db) {
        $database_array[$db] = 1;
      }
    }
    else {
      $database_array[$database] = 1;
    }
  }
  if (empty($database_array)) {
    return 'No database was found to search with.';
  }

  $software_exec = '';
  if (file_exists($software)) {
    $software_exec = $software;
  }
  else {
    $select_software_setting = 'SELECT value FROM {gmod_dbsf_softwareprop} as softwareprop where software_id='
     . "(SELECT software_id from {gmod_dbsf_software} as software where uniquename='%s')"
     . ' AND rank=0 AND type_id='
     . "(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id "
     . "where cv.name='software_setting' AND cvterm.name='%s')";
    $res = db_fetch_array(db_query($select_software_setting, $software, 'executable'));
    $software_exec = $res['value'];
    if (empty($software_exec)) {
      drupal_set_message(t('Failed to find an executable for %software', array('%software' => $software)), 'error');
      return ' ';
    }
  }

  // run $software
  $exec_output = array();
  $html_output = '';
  $fasta_data  = '';
  foreach ($database_array as $database_str => $dummy) {
    $database_str=trim($database_str);
    if (empty($database_str)) {
      continue;
    }
    $args = " -d $database_str -s $seqs_to_get";
    $result = exec($software_exec . $args, $exec_output);

//dpm($software_exec . $args);
  if (empty($exec_output)) {
    drupal_set_message(t('Sorry, cannot find your sequence entries (%seqs) in any database. If you want, you can go back and click on the graph to see the HTML report or ',array('%seqs' => $seqs_to_get)). l('contact', 'contact') .t(' the administrator.'),'error');
  }
  foreach ($exec_output as $line) {
    if (strpos($line, '>') !== FALSE) {
      $line = str_replace(' No definition line found', '', $line);
      $fasta_data .= $line ."\n";
      $html_output .= '</p><p>'. $line .'<br>';
    }
    else {
      $fasta_data .= $line ."\n";
      $html_output .= $line .'<br>';
    }
  }
  }
  $return_text = '';
  $outfile = $infile . time() .'.get';
  // OS compatibility
  ini_set('auto_detect_line_endings', TRUE);
  if (!$outhandle = fopen($outfile, 'wt')) {
    drupal_set_message(t("Could not create output FASTA file %outfile.", array('%outfile' => $outfile)), 'error');
  }
  else {
    fwrite($outhandle, $fasta_data);
    fclose($outhandle);
    $return_text .= '<p>You can download a '. l('FASTA file', base_path() . $outfile,
      array('external' => TRUE, 'attributes' => array('target' => '_blank'))
    ) .' or if small enough, it is available for copy-pasting below.</p>';
  }
  if (strlen($html_output) < 100000) {
    $return_text .= '<pre>'. $html_output .'</pre>';
  }
  return $return_text;
}

