<?php
// $Id$

/**
 * @file
 * Universal Search Server Module plugin
 * A Drupal module to implement iprscan (InterProScan) usage via the Universal Search Server
 * Developed by Alexie Papanicolaou. Licensed under GPLv3
 * University of Exeter. Licensed under GPLv3
 * @see http://insectacentral.org/
 *
 * This module was built using iprscan V 4.5
 *
 * by
 * Zdobnov E.M. and Apweiler R.
 *    "InterProScan - an integration platform for the signature-recognition methods in InterPro."
 *     Bioinformatics, 2001, 17(9): 847-8.
 *
 *  [1] The InterPro Consortium (*R.Apweiler, T.K.Attwood, A.Bairoch, A.Bateman,
 *  E.Birney, M.Biswas, P.Bucher, L.Cerutti, F.Corpet, M.D.R.Croning, R.Durbin,
 *  L.Falquet, W.Fleischmann, J.Gouzy, H.Hermjakob, N.Hulo, I.Jonassen, D.Kahn,
 *  A.Kanapin, Y.Karavidopoulou, R.Lopez, B.Marx, N.J.Mulder, T.M.Oinn, M.Pagni,
 *  F.Servant, C.J.A.Sigrist, E.M.Zdobnov).
 *  "The InterPro database, an integrated documentation resource for protein
 *  families, domains and functional sites."
 *  Nucleic Acids Research, 2001, 29(1): 37-40.
 *
 */

/*
 * We assume the iprscan has been already setup and is working well.
 *
 * For server usage we will implement
 *  iprscan (wrapper for all software)
 *
 * The following functions are needed in this include file for iprscan to operate
 *  biosoftware_bench_admin_iprscan_page
 *  biosoftware_bench_iprscan_page
 *  biosoftware_bench_add_software_iprscan
 *
 * The following can be done conditionally via the main module
 *  administration-> where the executable is
 *
 * The following must be done conditionally via the main .install file
 *  add iprscan as a software in ontology
 *  add options as CV
 *
 **/
function biosoftware_bench_iprscan_menu() {
  $items = array();
  $items['admin/bench/iprscan'] = array(
    'file' => 'includes/biosoftware_bench_iprscan.inc',
    'title' => 'iprscan Server administration',
    'page callback' => 'biosoftware_bench_admin_iprscan_page',
    'access arguments' => array('Administrate iprscan Server'),
    'description' => 'Configure iprscan specific settings',
    'type' => MENU_CALLBACK,
  );
  $items['bench/iprscan'] = array(
    'file' => 'includes/biosoftware_bench_iprscan.inc',
    'title' => 'biosoftware_bench iprscan Server',
    'page callback' => 'biosoftware_bench_iprscan_page',
    'access arguments' => array('Access iprscan Server'),
    'type' => MENU_CALLBACK,
  );
  $items['bench/iprscan_result'] = array(
    'file' => 'includes/biosoftware_bench_iprscan.inc',
    'title' => 'biosoftware_bench iprscan Results',
    'page callback' => 'biosoftware_bench_iprscan_result_page',
    'page arguments' => array(2, 3),
    'access arguments' => array('Access iprscan Server'),
    'type' => MENU_CALLBACK,
  );
  return $items;
}

function biosoftware_bench_software_iprscan_core_settings_form($form_state) {
  $software = 'iprscan';
  $form     = array();
  $res      = db_fetch_array(db_query("SELECT active from {gmod_dbsf_software} where uniquename='%s'", $software));
  if ($res['active'] == 'f') {
    return $form;
  }
  $select_software_setting = 'SELECT value FROM {gmod_dbsf_softwareprop} where software_id=' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ' AND rank=0 AND type_id=' . "(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')";
  $core_exec = db_fetch_array(db_query($select_software_setting, $software, 'executable'));
  $format_exec = db_fetch_array(db_query($select_software_setting, 'ipr_convert', 'executable'));
    $external_server = gmod_dbsf_get_add_var('biosoftware_bench_external_server');

  $form['settings'] = array(
    '#type' => 'fieldset',
    '#title' => $software . t(' core settings'),
    '#description' => t('Please set some important settings for %s', array('%s' => $sofware)),
    'submit' => array(
      '#type' => 'submit',
      '#value' => "Save core $software settings",
      '#weight' => 5,
    ),
    'core_executable' => array(
      '#type' => 'textfield',
      '#title' => t('Core executable path'),
      '#description' => t('Please provide the <strong>full path</strong> to the <strong>%software</strong> executable on the server, e.g. /usr/bin/%software .', array('%software' => $software)),
      '#required' => TRUE,
      '#default_value' => $core_exec['value'],
    ), 'format_executable' => array(
      '#type' => 'textfield',
      '#title' => t('InterProScan format converter path'),
      '#description' => t('Please provide the <strong>full path</strong> to the <strong>converter.pl</strong> executable on the server. It is part of the iprscan package in the bin directory, e.g. /software/iprscan/bin/converter.pl .'),
      '#required' => TRUE,
      '#default_value' => $format_exec['value'],
    ), 'info' => array(
      '#type' => 'fieldset',
      '#title' => t('Information'),
      '#description' => t("If you are new to iprscan and you're having problems to make it work (e.g. with newer processors), it helps if you install EMBOSS locally and specify the correct " . "environmental variables EMBOSS_ROOT, EMBOSS_ACDROOT and EMBOSS_DATA (e.g. /usr/share/EMBOSS/, /usr/share/EMBOSS/acd, /usr/share/EMBOSS/data/ respectively) " . "in the files iprscan/conf/seqret.sh and sixpack.sh. See iprscan/conf/iprscan.conf for more info."
      ),
    ),
    'external_server'=> array(
      '#type'=>'checkbox',
      '#title'=> t('I\' using an external server'),
      '#description'=>t('Check this if you do not want Drupal to check if the file paths actually exist on the web-server (because you are running it on an external machine'),
      '#weight'=> 40,
      '#default_value'=>$external_server,
    ),
  );
  return $form;
}

function biosoftware_bench_software_iprscan_core_settings_form_validate($form, &$form_state) {
  $software    = 'iprscan';
  $core_exec   = escapeshellcmd(trim($form_state['values']['core_executable']));
  $format_exec = escapeshellcmd(trim($form_state['values']['format_executable']));
  if (!empty($form_state['values']['external_server'])){
    gmod_dbsf_get_add_var('biosoftware_bench_external_server','TRUE');
  }elseif (gmod_dbsf_get_add_var('biosoftware_bench_external_server')){
    gmod_dbsf_get_add_var('biosoftware_bench_external_server','DELETE');
  }
  if (!empty($core_exec)) {
    if (!file_exists($core_exec) && empty($form_state['values']['external_server'])) {
      form_set_error('core_executable', t('Cannot find/access executable %exec on the server.', array('%exec' => $core_exec)));
    }
    elseif (!is_executable($core_exec) && empty($form_state['values']['external_server'])) {
      form_set_error('core_executable', t('Cannot execute the %f program.', array('%f' => $core_exec)));
      return FALSE;
    }
    $form_state['values']['core_executable'] = $core_exec;
  }
  if (!empty($format_exec)) {
    if (!file_exists($format_exec) && empty($form_state['values']['external_server'])) {
      form_set_error('format_executable', t('Cannot find/access executable %exec on the server.', array('%exec' => $format_exec)));
    }
    elseif (!is_executable($format_exec) && empty($form_state['values']['external_server'])) {
      form_set_error('format_executable', t('Cannot execute the %f program.', array('%f' => $format_exec)));
      return FALSE;
    }
    $form_state['values']['format_executable'] = $format_exec;
  }
}

function biosoftware_bench_software_iprscan_core_settings_form_submit($form, &$form_state) {
  $software = 'iprscan';
  $store_software_setting = 'INSERT INTO {gmod_dbsf_softwareprop} (software_id,type_id,rank,value) VALUES (' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ",(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')" . ",0,'%s')";
  $delete_software_setting = 'DELETE FROM {gmod_dbsf_softwareprop} where software_id=' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ' AND rank=0 AND type_id=' . "(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')";
  $select_software_setting = 'SELECT value FROM {gmod_dbsf_softwareprop} where software_id=' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ' AND rank=0 AND type_id=' . "(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')";

  $core_exec        = db_fetch_array(db_query($select_software_setting, $software, 'executable'));
  $format_exec      = db_fetch_array(db_query($select_software_setting, 'ipr_convert', 'executable'));
  $user_core_exec   = $form_state['values']['core_executable'];
  $user_format_exec = $form_state['values']['format_executable'];
  if (!empty($user_core_exec) && $user_core_exec != $core_exec['value']) {
    db_query($delete_software_setting, $software, 'executable');
    db_query($store_software_setting, $software, 'executable', $user_core_exec);
    drupal_set_message(t('%exec is now the %software executable.', array('%exec' => $user_core_exec, '%software' => $software)), 'warning');
  }
  if (!empty($user_format_exec) && $user_format_exec != $format_exec['value']) {
    db_query($delete_software_setting, 'ipr_convert', 'executable');
    db_query($store_software_setting, 'ipr_convert', 'executable', $user_format_exec);
    drupal_set_message(t('%exec is now the InterProScan converter.pl executable.', array('%exec' => $user_format_exec)), 'warning');
  }
}

function biosoftware_bench_admin_iprscan_page($software = 'iprscan') {
  // Include biosoftware_bench admin function file
  require_once(drupal_get_path('module', 'biosoftware_bench') .'/includes/biosoftware_bench_admin.inc');
  $check_active = biosoftware_bench_check_software_active($software);
  if (empty($check_active)) {
    drupal_set_message(t("I'm sorry, but %software does not seem to have been activated yet.<br>", array('%software' => $software)) .
      l(t('See the software settings page'), 'admin/bench/software'), 'error'
    );
    return FALSE;
  }
  $path_form = drupal_get_form('biosoftware_bench_software_iprscan_core_settings_form');

  $page_tabs = array(
    'settings' => array(
      '#type' => 'tabset',
      'core' => array(
        '#type' => 'tabpage',
        '#title' => 'Core settings',
        '#content' => $path_form,
        '#weight' => 0,
      ),
    ),
  );
  $return_text = tabs_render($page_tabs);

  return $return_text;
}

function biosoftware_bench_check_iprscan($software = 'iprscan') {
  $check_active = biosoftware_bench_check_software_active($software);
  if (empty($check_active)) {
    return 'active';
  }
  $check = 0;
  $software_prop = gmod_dbsf_get_softwareprop(NULL, TRUE);
  if (!empty($software_prop)) {
    foreach ($software_prop as $name => $data) {
      if (in_array($name, array($software))) {
        foreach ($data as $term_name => $v) {
          if ($term_name == 'executable') {
            $check++;
          }
        }
      }
    }
  }
  if ($check < 1) {
    return 'variables';
  }

  return TRUE;
}

function biosoftware_bench_iprscan_page() {
  $setup = biosoftware_bench_check_iprscan('iprscan');
  if ($setup === TRUE) {
    return drupal_get_form('biosoftware_bench_iprscan_form');
  }
  else {
    $return_text = "<p>I'm sorry but your iprscan software has not been setup properly.</p><ul>";
    switch ($setup) {
      case 'active':
        $return_text .= '<li>iprscan has not been activated. Please activate it at the '. l('software settings page', 'admin/bench/software');
        break;

      case 'variables':
        $return_text .= '<li>iprscan variable paths have not been defined. Please define them at at the '. l('iprscan settings page', 'admin/bench/software', array('fragment' => 'iprscan'));
        break;
    }
    drupal_set_message(t('!t', array('!t' => $return_text .'</ul>')), 'error');
    return ' ';
  }
}

function biosoftware_bench_iprscan_form($form_state) {
  $form             = array();
  $software         = 'iprscan';
  $software_options = biosoftware_bench_generate_options($software);
  $option_form      = array();
  if (!empty($software_options)) {

    foreach ($software_options as $name => $data) {
      $description = !empty($data['description']) ? t('@d', array('@d' => $data['description'])) : '';
      $title       = str_replace($software, '', $name);
      $title       = ucfirst(trim(str_replace('_', ' ', $title)));
      if ($data['type'] == 'checkbox') {
        unset($data['options']);
      }
      $option_form[$name] = array(
        '#type' => $data['type'],
        '#title' => $title,
        '#description' => $description,
      );
      if ($data['type'] !== 'textfield' && $data['type'] !== 'checkbox' && !empty($data['options'])) {
        $option_form[$name]['#options'] = $data['options'];
        if ($data['type'] !== 'checkboxes') {
          $option_form[$name]['#default_value'] = current($data['options']);
        }
        //interproscan specific to set all checkboxes as checked.
        else {
          $option_form[$name]['#default_value'] = $data['options'];
        }
      }
    }
  }
  $option_form['iprlookup']['#default_value'] = 1;
  $option_form['goterms']['#default_value'] = 1;
  $software_description = '';
  $description_file = drupal_get_path('module', 'biosoftware_bench') .'/includes/biosoftware_bench_'. $software .'_description.html';
  if (file_exists($description_file) && $inhandle = fopen($description_file, 'rb')) {
    while (!feof($inhandle)) {
      $software_description .= (fgets($inhandle));
    }
  }
  if (!empty($software_description)) {
    $form['description'] = array(
      '#type' => 'fieldset',
      '#title' => 'Description',
      '#collapsible' => TRUE,
      '#collapsed' => TRUE,
      '#description' => $software_description,
    );
  }
  $form['program'] = array(
    '#type' => 'fieldset',
    '#title' => strtoupper($software),
    'query_sequence' => array(
      '#type' => 'textarea',
      '#title' => t('Enter query sequence in simple text, FASTA format'),
    ),
    'query_file_'. $software => array(
      '#type' => 'file',
      '#title' => t('or upload query sequence in FASTA format'),
      '#description' => t('Please give a text file, not a MS-Word or other document, you can upload up to %m Mb.', array('%m' => gmod_dbsf_get_add_var('biosoftware_bench_upload_size'))),
      // Needed because of drupal bug
      '#tree' => FALSE,
    ),
  );
  if (!empty($option_form)) {
    $form['program']['advanced parameters'] = array(
      '#type' => 'fieldset',
      '#title' => 'Advance Search Parameters',
      '#collapsible' => TRUE,
      '#collapsed' => FALSE,
      $option_form,
    );
  }
  $form['buttons'] = array(
    '#weight' => 10,
    'db_submit' => array(
      '#type' => 'submit',
      '#value' => t('Run '. $software),
    ),
    'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'submit',
      '#value' => 'Reset data',
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 2,
    ),
  );
  $form['program']['seqtype'] = array(
    '#type' => 'value',
    '#value' => 'p',
  );
  $form['#attributes']['enctype'] = "multipart/form-data";
  return $form;
}

function biosoftware_bench_iprscan_form_validate($form, &$form_state) {
  $software = 'iprscan';
  if (empty($form_state['values']['query_sequence']) && empty($_FILES['files']['tmp_name']['query_file_'. $software])) {
    form_set_error("query_sequence", t('It seems you neither gave a query sequence nor uploaded one.'));
    return FALSE;
  }
  if (!empty($form_state['values']['query_sequence'])) {
    if (!empty($_FILES['files']['tmp_name']['query_file_'. $software])) {
      form_set_error("query_sequence", t('It seems you both typed a query sequence and uploaded one.'));
      return FALSE;
    }
    // First Query (i.e. check for query first; query must exist if subject exists)
    $result = gmod_dbsf_validate_seq_protein($form_state['values']['query_sequence']);
    if ($result === FALSE || is_numeric($result)) {
      $result = gmod_dbsf_validate_seq_dna($form_state['values']['query_sequence']);
      if ($result === FALSE || is_numeric($result)) {
        form_set_error("query_sequence", t('Sorry your sequence does not seem to be a valid DNA or protein sequence (%result errors).', array('%result' => $result)));
        return FALSE;
      }
      else {
        $form_state['values']['seqtype'] = 'n';
        $form_state['values']['query_sequence'] = $result;
      }
    }
    else {
      $form_state['values']['query_sequence'] = $result;
    }
  }
  if (empty($form_state['values']['iprscan_applications'])) {
    form_set_error('ipr_applications', t('You must select at least one application.'));
  }
  if (!empty($form_state['values']['iprscan_goterms']) && empty($form_state['values']['iprscan_iprlookup'])) {
    $form_state['values']['iprscan_iprlookup'] = 1;
  }
}

function biosoftware_bench_iprscan_form_submit($form, &$form_state) {
  // -cli -i seq -o seq.ipr_out -altjobs -nocrc -seqtype p -goterms -iprlookup -format raw -appl as list (no commas, quotes etc)
  $seqtype   = check_plain($form_state['values']['seqtype']);
  $crc       = $form_state['values']['iprscan_crc'];
  $iprlookup = $form_state['values']['iprscan_iprlookup'];
  $goterms   = $form_state['values']['iprscan_goterms'];
  $par       = " -cli -format raw -seqtype $seqtype ";
  foreach ($form_state['values']['iprscan_applications'] as $value) {
    if (!empty($value)) {
      $par .= ' -appl '. check_plain($value);
    }
  }
  if (empty($par)) {
    form_set_error('ipr_applications', 'No programs selected');
  }
  if (empty($form_state['values']['crc'])) {
    $par .= ' -nocrc';
  }
  if (!empty($form_state['values']['iprlookup'])) {
    $par .= ' -iprlookup';
  }
  if (!empty($form_state['values']['goterms'])) {
    $par .= ' -goterms';
  }

  $software = 'iprscan';
  $file_size_allowed = gmod_dbsf_get_add_var('biosoftware_bench_upload_size');

  $validators_file = array('file_validate_size' => array($file_size_allowed * 1000));
  $tmpdir = file_directory_temp();
  // all files, same timestamp
  $timestamp = time();
  $sessionid = session_id();
  $dirpath   = file_create_path() .'/bench';
  $data      = $form_state['values'];
  // BATCH API
  $operations = array();
  // for future use
  $save_array = array();
  // the form_uid will allow users to have multiple windows open.
  $form_uid = gmod_dbsf_create_uid($sessionid, $timestamp, $software);
  $batch_file_data = array();

  $file_type = 'protein';
  $verify    = array();
  $uid       = gmod_dbsf_create_uid($sessionid, $timestamp, $software);
  if (isset($_FILES['files']) && !empty($_FILES['files']['tmp_name']['query_file_'. $software])) {
    $file = file_save_upload('query_file_'. $software, $validators_file, $tmpdir, FILE_EXISTS_RENAME);
    if (empty($file)) {
      form_set_error('query_file_'. $software, t('Sorry your file for <em>iprscan</em> was not saved. Maybe it is too large (>%file_size_allowed Mb)? Otherwise, '. l('contact', 'contact') .' the administrator  (quote %uid).', array('%uid' => $uid, '%file_size_allowed' => $file_size_allowed)));
    }
    else {
      $verify['runq']++;
      file_set_status($file, FILE_STATUS_TEMPORARY);
      $batch_file_data['infile'][$uid] = $file->filepath;
      $batch_file_data['outfile'][$uid] = $dirpath .'/'. $uid .'.query';
      $batch_file_data['filetype'][$uid] = $file_type;
      $batch_file_data['format'][$uid] = FALSE;
    }
  }
  elseif (!empty($data['query_sequence'])) {
    $verify['runq']++;
    // textfield
    $batch_file_data['infile'][$uid] = $data['query_sequence'];
    $batch_file_data['outfile'][$uid] = $dirpath .'/'. $uid .'.query';
    $batch_file_data['filetype'][$uid] = $file_type;
    $batch_file_data['format'][$uid] = FALSE;
  }
  else {
    form_set_error('', t('No data given!'));
    return FALSE;
  }
  //if it is going to run (both query and subject are ok)
  if (!empty($verify['runq']) && $verify['runq'] == 1) {
    $save_array[$form_uid][$uid]['par'] = $par;
    $save_array[$form_uid][$uid]['algorithm'] = $software;
  }

  // ended each algorithm
  $operations[] = array('gmod_dbsf_batch_upload_fasta', array($batch_file_data));
  $operations[] = array('gmod_dbsf_batch_save_data', array($save_array));
  $batch        = array(
    'title' => t('Preparing data needed for iprscan...'),
    'operations' => $operations,
    'init_message' => t('Starting iprscan submission...'),
    'progress_message' => t('@remaining operations remaining...'),
    'error_message' => t('Your iprscan submission encountered an error.'),
    'finished' => 'biosoftware_bench_iprscan_batch_finished',
    'file' => drupal_get_path('module', 'biosoftware_bench') .'/includes/biosoftware_bench_iprscan.inc',
  );
  batch_set($batch);
  $form_state['redirect'] = array('bench/iprscan_result', "submission_uid=$form_uid");
  return "bench/iprscan_result?submission_uid=$form_uid";
}



/* Implements hook_batch_finish
 *
 * Handles the final operations of the BATCH API
 * @see http://api.drupal.org/api/group/batch
 *
 * @param $success
 *   See Drupal's BATCH API
 * @param $results
 *   See Drupal's BATCH API
 * @param $operations
 *   See Drupal's BATCH API
 *
 * @return
 *   FALSE on failure
 */
function biosoftware_bench_iprscan_batch_finished($success, $results, $operations) {
  $software = 'iprscan';
  // store uids to report in case of error.
  $uid_array      = array();
  $submission_uid = '';
  $data           = array();
  $to_store       = "software:biosoftware_bench\n";
  $dirpath        = file_create_path() .'/bench';
  $message        = '';
  //one form uid
  foreach ($results['data'] as $form_uid => $batch_data) {
    //for use later
    $submission_uid = $form_uid;
    $to_store .= "submission_uid:$submission_uid\n";
    $data = $batch_data;
    foreach ($batch_data as $uid => $values) {
      // this is not a tag value but tag;tag;tag... = value;value;value...
      $to_store .= 'algorithm;id='. $values['algorithm'] .";$uid\n";
      $uid_array[] = $uid;
    }
  }

  if (!empty($success)) {
    // store what the users is performing in this submission, so that the results
    // can be accurately retrieved.
    $outfile = $dirpath .'/'. $submission_uid .'.submission';
    if (!$outhandle = fopen($outfile, 'wb')) {
      drupal_set_message(t('Could not create %outfile.', array('%outfile' => $outfile)), 'error');
      return FALSE;
    }
    fwrite($outhandle, $to_store);
    fclose($outhandle);
    foreach ($data as $uid => $values) {
      gmod_dbsf_parameter_daemon($uid, $dirpath, $values['par'],
        $software, array(), array()
      );
    }
    $message = t('Your request has been submitted to the queue for processing...');
  }
  else {
    $error_operation = reset($operations);
    $uids            = implode(' ', $uid_array);
    $message         = t('An error occurred while processing your request. Please '. l('contact', 'contact') ." the administrator reporting: $uids.");
    //dpm($error_operation);
    drupal_set_message($message, 'warning');
    return FALSE;
  }
  drupal_set_message($message, 'warning');
}

/**
 * Present software results by providing the UID
 *
 * This function must be customized
 */
function biosoftware_bench_iprscan_result_page($submission_uid = NULL, $level = 1) {
  if (empty($submission_uid)) {
    $submission_uid = $_GET['submission_uid'];
    if (empty($submission_uid)) {
      return 'No submission data provided';
    }
  }
  if (empty($level)) {
    $level = 1;
  }
  $basepath = file_create_path();
  $dirpath  = $basepath .'/bench/';
  $data     = array();
  // one algorithm type for each submission.
  // we get this from a file.
  $infile = $dirpath . $submission_uid .'.submission';
  if (file_exists($infile)) {
    if ($inhandle = fopen($infile, 'rb')) {
      while (!feof($inhandle)) {
        $line = trim(fgets($inhandle));
        if (!empty($line) && preg_match('/^algorithm;id=/', $line)) {
          $line_data          = explode('=', $line);
          $alg_data           = explode(';', $line_data[1]);
          $data[$alg_data[0]] = $alg_data[1];
        }
      }
    }
  }
  else {
    //dpm("Did not find $infile");
    return 'Submission ID not found. Maybe it has expired or is invalid?';
  }
  if (empty($data)) {
    //dpm("$infile has no data");
    return 'Not a valid submission ID or ID has expired.';
  }
  $content = '';
  foreach ($data as $algorithm => $uid) {
    $outfile   = $dirpath . $uid .".$algorithm.output";
    $graph     = $outfile .'.graph.png';
    $errorfile = $dirpath . $uid .".$algorithm.error";
    /*PLUGIN_EDIT_HERE
     * link name
     */


    $content = t('Your report is not ready yet. Please wait or ') . l('try again', 'bench/iprscan_result/'. $submission_uid
    ) . t(' in a few moments or bookmark this page.');
    if (file_exists($graph)) {
      $results = drupal_get_form('biosoftware_bench_iprscan_result_form', $uid, $outfile .'.gff', $level, $submission_uid);
      if (strlen($results) < 600) {
        $results = "No (more) hits found for any of your queries using $algorithm. Perhaps you can try a different algorithm or database?";
      }
      $content = '<p>'. l(
        theme_image($graph, 'IPRSCAN results of first ten queries', 'IPRSCAN results graphical overview', array(), FALSE)
        , base_path() . $outfile .'.html', array('html' => TRUE, 'external' => TRUE, 'attributes' => array('target' => '_blank'))
      ) .'</p>' . '<p>All results as :';
      if (file_exists($outfile .".html")) {
        $content .= l(
          theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/html.png', 'Report as HTML', 'Report as HTML'), base_path() . $outfile .".html",
          array(
            'html' => TRUE,
            'external' => TRUE,
            'attributes' => array('target' => '_blank'),
          )
        );
      }
      if (file_exists($outfile .".xml")) {
        $content .= l(
          theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/xml.png', 'Report as InterProScan XML', 'Report as InterProScan XML (recommended for parsing with other software)'), base_path() . $outfile .".xml",
          array(
            'html' => TRUE,
            'external' => TRUE,
            'attributes' => array('target' => '_blank'),
          )
        );
      }
      if (file_exists($outfile)) {
        $content .= l(theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/iprscan.png', 'Report as InterProScan', 'Report as InterProScan'), base_path() . $outfile,
          array(
            'html' => TRUE,
            'external' => TRUE,
            'attributes' => array('target' => '_blank'),
          )
        );
      }
      if (file_exists($outfile .".gff")) {
        $content .= l(theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/gff.png', 'Report as GFF3', 'Report as GFF3'), base_path() . $outfile .".gff",
          array(
            'html' => TRUE,
            'external' => TRUE,
            'attributes' => array('target' => '_blank'),
          )
        );
      }
      if (file_exists($outfile .".txt")) {
        $content .= l(
          theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/txt.png', 'Report as text', 'Report as text'), base_path() . $outfile .".txt",
          array(
            'html' => TRUE,
            'external' => TRUE,
            'attributes' => array('target' => '_blank'),
          )
        );
      }

      $content .= '</p>'. $results;
    }
    elseif (file_exists($errorfile) && filesize($errorfile) > 0) {
      // There was an error
      $errorhandle = fopen($errorfile, 'r');
      if (isset($errorhandle)) {
        $content = '<p>An error has been encountered with your submission. Please notify an administrator:</p></pre>';
        $content .= fread($errorhandle, 100) .'...';
        fclose($errorhandle);
      }
    }
    else {
      global $base_url;
      drupal_set_html_head('<meta http-equiv="refresh" content="15;url='. $base_url .'/bench/iprscan_result/'. $submission_uid .'" />');
    }
  }
  return $content;
}

/**
 * Implements hook_form()
 *
 * Used to produce a table of results with checkboxes which we can then
 * select for processing (e.g. downloading). This form is themed.
 * @see biosoftware_bench_theme_blast_result_form()
 * @see biosoftware_bench_blast_result_page()
 * @ingroup forms
 *
 * @param $form_state
 *   Form API variable
 * @param $uid
 *   Unique Identifier of algorithm's BLAST search
 * @param $blastfile
 *   Full path to BLAST results as an XML file
 */
function biosoftware_bench_iprscan_result_form($form_state, $uid, $resultfile, $level = 1, $suid = NULL) {
  $form = array();
  if (empty($resultfile)) {
    return $form;
  }
  if (empty($level)) {
    $level = 1;
  }
  $outfile           = str_replace('.gff', '', $resultfile);
  $array_from_result = array();
  $array_from_result = gmod_dbsf_parse_gfffile($resultfile);
  // all data
  $table_row = array();
  // for checkboxes to get hit
  $row_check_data = array();

  if (empty($array_from_result) || is_array($array_from_result) === FALSE) {
    return $form;
  }
  $weight_series = 0;
  $query_rank    = 0;
  $total_queries = count($array_from_result);
  foreach ($array_from_result as $query_id => $data_array) {
    $query_rank++;
    if ($query_rank != $level) {
      continue;
    }
    foreach ($data_array as $key => $gff_data) {
      if ($gff_data['type'] !== 'match_set') {
        continue;
      }
      $hit_id = $gff_data['attributes']['Name'][0];
      $other_info = implode(', ', (array)$gff_data['attributes']['Note']);
      $other_info .= '. '. implode(', ', (array)$gff_data['attributes']['Ontology_term']);
      $length_of_align = abs($gff_data['end'] - $gff_data['start']);
      $match_data = $data_array[$key + 1];
      if ($match_data['type'] == 'match_part') {
        $score = $match_data['score'];
      }
      $links = '';
      if (file_exists($outfile .".Query$query_rank.html")) {
        $links .= l(
          theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/html.png', 'Report as HTML', 'Report as HTML'), base_path() . $outfile .".Query$query_rank.html",
          array(
            'html' => TRUE,
            'external' => TRUE,
            'attributes' => array('target' => '_blank'),
          )
        );
      }
      if (file_exists($outfile .".Query$query_rank.raw")) {
        $links .= l(
          theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/xml.png', 'Report as InterProScan XML', 'Report as InterProScan XML (recommended for parsing with other software)'), base_path() . $outfile .".Query$query_rank.xml",
          array(
            'html' => TRUE,
            'external' => TRUE,
            'attributes' => array('target' => '_blank'),
          )
        );
      }
      if (file_exists($outfile .".Query$query_rank.raw")) {
        $links .= l(theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/iprscan.png', 'Report as InterProScan', 'Report as InterProScan'), base_path() . $outfile .".Query$query_rank.raw",
          array(
            'html' => TRUE,
            'external' => TRUE,
            'attributes' => array('target' => '_blank'),
          )
        );
      }
      if (file_exists($outfile .".Query$query_rank.gff")) {
        $links .= l(theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/gff.png', 'Report as GFF3', 'Report as GFF3'), base_path() . $outfile .".Query$query_rank.gff",
          array(
            'html' => TRUE,
            'external' => TRUE,
            'attributes' => array('target' => '_blank'),
          )
        );
      }
      if (file_exists($outfile .".Query$query_rank.txt")) {
        $links .= l(
          theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/txt.png', 'Report as text', 'Report as text'), base_path() . $outfile .".Query$query_rank.txt",
          array(
            'html' => TRUE,
            'external' => TRUE,
            'attributes' => array('target' => '_blank'),
          )
        );
      }

      $table_row[] = array($query_id, $hit_id, $score, $length_of_align, $other_info, $links);
      $form['features'][$weight_series] = array(
        '#type' => 'checkbox',
        '#title' => $query_id,
        '#return_value' => $hit_id,
        '#weight' => $weight_series,
      );
      $weight_series++;
    }
  }
  $pager = array();
  global $base_url;
  $dirname           = drupal_get_path('module', 'gmod_dbsf');
  $base              = $base_url .'/bench/iprscan_result/'. $uid;
  $pager['bookmark'] = l(theme_image($dirname .'/images/bookmark_page.png', 'Bookmark this page', 'Bookmark this page'), $base .'/'. $level,
    array(
      'html' => TRUE,
      // needed
      'external' => TRUE,
    )
  );;
  if ($level > 1) {
    $pager['previous'] = l(theme_image($dirname .'/images/previous-button.png', 'Previous Query', 'Previous Query'), $base .'/'. ($level - 1),
      array(
        'html' => TRUE,
        // needed
        'external' => TRUE,
        'attributes' => array('title' => 'Get next Query'),
      )
    );
  }
  if ($total_queries > $level) {
    $pager['next'] = l(theme_image($dirname .'/images/next-button.png', 'Next Query', 'Next Query'), $base .'/'. ($level + 1),
      array(
        'html' => TRUE,
        // needed
        'external' => TRUE,
        'attributes' => array('title' => 'Get next Query'),
      )
    );
  }
  $table_header     = array('Query ID', 'Hit ID', 'Score', 'Aln length', 'Notes', 'Links');
  $table_attributes = array('class' => 'padded-table sortable');
  $table_caption    = 'Overview of significant results';
  $table_array      = array(
    '#links' => $pager,
    'header' => $table_header,
    'caption' => $table_caption,
    'attributes' => $table_attributes,
    'data' => $table_row,
  );

  $form['resultfile'] = array('#type' => 'value', '#value' => $outfile);
  $form['data']       = array('#type' => 'value', '#value' => $table_array);
  $form['buttons']    = array(
    '#type' => 'fieldset',
    '#title' => t('Download hits of selected results as'),
    '#description' => t('If the reference database has been indexed then you can download the hits.'),
    '#collapsible' => FALSE,
    '#weight' => 20 + $weight_series,
    'FASTA' => array(
      '#type' => 'submit',
      '#value' => t('FASTA'),
      '#weight' => 1,
      /*'GFF' => array(   '#type' => 'submit',   '#value' => t('GFF'),   '#weight' => 16   )*/
    ),
  );
  $form['#action'] = url("bench/get_sequences/$uid/iprscan");
  return $form;
}

