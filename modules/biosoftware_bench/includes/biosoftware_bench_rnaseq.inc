<?php
// $Id$

/**
 * @file
 * rnaseq include file
 */
function biosoftware_bench_rnaseq_menu() {
  $items = array();
  $items['bench/rnaseq_result'] = array(
    'file' => 'includes/biosoftware_bench_rnaseq.inc',
    'title' => 'biosoftware_bench rnaseq Results',
    'page callback' => 'biosoftware_bench_rnaseq_result_page',
    'page arguments' => array(2, 3),
    'access arguments' => array('Access rnaseq Server'),
    'type' => MENU_CALLBACK,
  );
  $items['bench/rnaseq'] = array(
    'file' => 'includes/biosoftware_bench_rnaseq.inc',
    'title' => 'biosoftware_bench rnaseq Server',
    'page callback' => 'biosoftware_bench_rnaseq_page',
    'page arguments' => array(2),
    'access arguments' => array('Access rnaseq Server'),
    'type' => MENU_CALLBACK,
  );
  $items['admin/bench/rnaseq'] = array(
    'file' => 'includes/biosoftware_bench_rnaseq.inc',
    'title' => 'rnaseq Server administration',
    'page callback' => 'biosoftware_bench_admin_rnaseq_page',
    'access arguments' => array('Administrate rnaseq Server'),
    'description' => 'Configure rnaseq specific settings',
    'type' => MENU_CALLBACK,
  );
  return $items;
}


/**
 * Create HTML for the administration rnaseq page
 *
 * Menu as admin/bench/rnaseq
 */
function biosoftware_bench_admin_rnaseq_page() {
  $software = 'rnaseq';
  $return_text='';
  require_once(drupal_get_path('module', 'biosoftware_bench') .'/includes/biosoftware_bench_admin.inc');
  $check_active = biosoftware_bench_check_software_active($software);
  if (empty($check_active)) {
    drupal_set_message(t("I'm sorry, but %software does not seem to have been activated yet.<br>", array('%software' => $software)) .
      l('See the software settings page', 'admin/bench/software'), 'error'
    );
    return FALSE;
  }
  $database_link_form = drupal_get_form('biosoftware_bench_admin_link_dataset_software_form', $software);
  if (strlen($database_link_form) < 600) {
    $database_link_form = "You have no datasets defined to link them to $software. Please see ". l('the dataset administration page', 'admin/bench/dataset');
  }
  $table = biosoftware_bench_admin_database_table($software);

  $page_tabs = array(
    'settings' => array(
      '#type' => 'tabset',
      'variables' => array(
        '#type' => 'tabpage',
        '#title' => 'Link datasets to '. strtoupper($software),
        '#content' => $database_link_form . $table,
        '#weight' => 0,
      ),
    ),
  );
  $setup = biosoftware_bench_check_rnaseq('rnaseq');
  if ($setup===TRUE){
      $return_text.='<p>rnaseq is setup correctly.</p>';
  }else{
    $return_text.='<p>rnaseq is <b>not setup</b> correctly.'
    ." For rnaseq to be activated, at least 1 protein and 1 nucleotide dataset must be available.</p>";
  }
  $return_text .= "<p>New datasets can be registered ".l('here','admin/bench/dataset').'.</p>'.tabs_render($page_tabs);
  return $return_text;
}


/**
 * Implements hook_form()
 *
 * @see biosoftware_bench_software_rnaseq_core_settings_form_validate()
 * @see biosoftware_bench_software_rnaseq_core_settings_form_submit()
 * @ingroup forms
 *
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_software_rnaseq_core_settings_form($form_state) {
  $form     = array();
  $software = 'rnaseq';
  $active   = biosoftware_bench_check_software_active($software);
  if ($active === FALSE) {
    return $form;
  }
  $select_software_setting = 'SELECT value FROM {gmod_dbsf_softwareprop} where software_id='
  . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')"
  . ' AND rank=0 AND type_id='
  . "(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')";
  $external_server = gmod_dbsf_get_add_var('biosoftware_bench_external_server');

  $bwa_exec        = db_fetch_array(db_query($select_software_setting, 'bwa', 'executable'));
  $samtools_exec   = db_fetch_array(db_query($select_software_setting, 'samtools', 'executable'));

  $form['settings'] = array(
    '#type' => 'fieldset',
    '#title' => t('rnaseq core settings'),
    '#description' => t('Please set some important settings for the rnaseq program.'),
    'submit' => array(
      '#type' => 'submit',
      '#value' => 'Save core rnaseq settings',
      '#weight' => 5,
    ),
    'bwa_executable' => array(
      '#type' => 'textfield',
      '#title' => t('BWA executable path'),
      '#description' => t('Please provide the <strong>full path</strong> to the <strong>BWA</strong> executable on the server, e.g. /usr/bin/bwa .'),
      '#required' => TRUE,
      '#default_value' => $bwa_exec['value'],
    ),
    'samtools_executable' => array(
      '#type' => 'textfield',
      '#title' => t('Samtools executable path'),
      '#description' => t('Please provide the <strong>full path</strong> to the <strong>samtools</strong> executable on the server, e.g. /usr/bin/samtools .'),
      '#required' => TRUE,
      '#default_value' => $samtools_exec['value'],
    ),
    'external_server'=> array(
      '#type'=>'checkbox',
      '#default_value'=>$external_server,
      '#title'=> t('I\'m using an external server'),
      '#description'=>t('Check this if you do not want Drupal to check if the file paths actually exist on the web-server (because you are running it on an external machine'),
      '#weight'=> 40,
    ),
  );
  return $form;
}

/**
 * Implements hook_form_validate()
 *
 * @see biosoftware_bench_software_rnaseq_core_settings_form()
 * @see biosoftware_bench_software_rnaseq_core_settings_form_submit()
 *
 * @param $form
 *   hook_form data
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_software_rnaseq_core_settings_form_validate($form, &$form_state) {
  $software = 'rnaseq';

  $bwa_exec     = escapeshellcmd(trim($form_state['values']['bwa_executable']));
  $samtools_exec   = escapeshellcmd(trim($form_state['values']['samtools_executable']));
  if (!empty($form_state['values']['external_server'])){
    gmod_dbsf_get_add_var('biosoftware_bench_external_server','TRUE');
  }elseif (gmod_dbsf_get_add_var('biosoftware_bench_external_server')){
    gmod_dbsf_get_add_var('biosoftware_bench_external_server','DELETE');
  }
  if (!empty($data) && substr($data, -1) !== '/') {
    $data .= '/';
  }
  if (!empty($bwa_exec)) {
    if (!file_exists($bwa_exec) && empty($form_state['values']['external_server'])) {
      form_set_error('bwa_executable', t('Cannot find executable %bwa_exec on the server.', array('bwa_exec' => $bwa_exec)));
    }
    elseif (!is_executable($bwa_exec) && empty($form_state['values']['external_server'])) {
      form_set_error('rnaseq_executable', t('Cannot execute the %f program.', array('%f' => $bwa_exec)));
      return FALSE;
    }
    $form_state['values']['rnaseq_executable'] = $bwa_exec;
  }
}

/**
 * Implements hook_form_submit()
 *
 * @see biosoftware_bench_software_rnaseq_core_settings_form_validate()
 * @see biosoftware_bench_software_rnaseq_core_settings_form()
 *
 * @param $form
 *   hook_form data
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_software_rnaseq_core_settings_form_submit($form, &$form_state) {
  $software = 'rnaseq';

  $store_software_setting = 'INSERT INTO {gmod_dbsf_softwareprop} (software_id,type_id,rank,value) VALUES (' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ",(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')" . ",0,'%s')";
  $delete_software_setting = 'DELETE FROM {gmod_dbsf_softwareprop} where software_id=' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ' AND rank=0 AND type_id=' . "(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')";
  $select_software_setting = 'SELECT value FROM {gmod_dbsf_softwareprop} where software_id=' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ' AND rank=0 AND type_id=' . "(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')";

  $bwa_exec     = db_fetch_array(db_query($select_software_setting, 'bwa', 'executable'));
  $samtools_exec   = db_fetch_array(db_query($select_software_setting, 'samtools', 'executable'));

  $user_bexec         = $form_state['values']['bwa_executable'];
  $user_fexec         = $form_state['values']['samtools_executable'];

  if (!empty($user_bexec) && $user_bexec != $bwa_exec['value']) {
    db_query($delete_software_setting, $software, 'executable');
    db_query($store_software_setting, $software, 'executable', $user_bexec);
    drupal_set_message(t('%user_bexec is now the BWA executable.', array('%user_bexec' => $user_bexec)), 'warning');
  }
  if (!empty($user_fexec) && $user_fexec != $samtools_exec['value']) {
    db_query($delete_software_setting, 'formatdb', 'executable');
    db_query($store_software_setting, 'formatdb', 'executable', $user_fexec);
    drupal_set_message(t('%user_fexec is now the samtools executable.', array('%user_fexec' => $user_fexec)), 'warning');
  }
}

/**
 * Checks if rnaseq is setup correctly
 *
 * Use this function as a template for plugins.
 *
 * @param $software
 *   The name of software to check
 *
 * @return
 *   TRUE if software is correctly setup. Otherwise it return a string
 *   depending what is wrong:
 *   - variables: One of the variables such as execuutable path
 *   - datasets: Not enough datasets have been linked to the software
 *   - active: It has not been activated in the main settings page
 */
function biosoftware_bench_check_rnaseq($software = 'rnaseq') {
  $check_active = biosoftware_bench_check_software_active($software);
  if (empty($check_active)) {
    return 'active';
  }
  $check = 0;
  $software_prop = gmod_dbsf_get_softwareprop(NULL, TRUE);
  if (!empty($software_prop)) {
    foreach ($software_prop as $name => $data) {
      if (!empty($data) && in_array($name, array($software, 'formatdb', 'fastacmd'))) {
        foreach ($data as $term_name => $v) {
          if ($term_name == 'executable' || $term_name == 'data') {
            $check++;
          }
        }
      }
    }
  }
  if ($check < 2) {
    return 'variables';
  }
  $sql = 'SELECT count(distinct rcv.cvterm_id) from {gmod_dbsf_resource_cvterm} rcv ' . ' JOIN {gmod_dbsf_software_resource} sr ON sr.resource_id=rcv.resource_id ' . " WHERE sr.software_id=(SELECT software_id from {gmod_dbsf_software} where uniquename='rnaseq') " . " AND rcv.cvterm_id IN (SELECT cvterm_id FROM {gmod_dbsf_cvterm} " . " WHERE cv_id=(SELECT cv_id from {gmod_dbsf_cv} where name='dataset_type'))";
  $count = db_fetch_array(db_query($sql));
  if ($count['count'] < 2) {
    return 'datasets';
  }

  return TRUE;
}


/**
 * Creates HTML for the rnaseq page
 *
 * Menu as bench/rnaseq
 *
 * @param $sample_data
 *   If true then the form will be prefilled with sample data
 *
 * @return
 *   HTML of rnaseq page
 */
function biosoftware_bench_rnaseq_page() {
  $setup = biosoftware_bench_check_rnaseq('rnaseq');
  if ($setup === TRUE) {
    $return = '';
    $link = base_path() . drupal_get_path('module', 'biosoftware_bench') .'/examples/';
    $return .= drupal_get_form('biosoftware_bench_rnaseq_form',$switch);
    return $return;
  }
  else {
    $return_text = "<p>I'm sorry but your rnaseq software has not been setup properly.</p><ul>";
    switch ($setup) {
      case 'datasets':
        $return_text .= '<li>Not enough datasets have been defined. Please define at least one protein and one nucleotide dataset at the '. l('dataset registration page', 'admin/bench/dataset', array('fragment' => 'available-datasets')) .' and then register them with rnaseq at the '. l('rnaseq database settings page', 'admin/bench/rnaseq');
        break;

      case 'active':
        $return_text .= '<li>rnaseq has not been activated. Please activate it at the '. l('software settings page', 'admin/bench/software');
        break;

      case 'variables':
        $return_text .= '<li>rnaseq variable paths have not been defined. Please define them at at the '. l('rnaseq settings page', 'admin/bench/software', array('fragment' => 'rnaseq'));
        break;
    }
    drupal_set_message(t('!t', array('!t' => $return_text .'</ul>')), 'error');
    return ' ';
  }
}

/**
 * Implements hook_form()
 *
 * @see biosoftware_bench_rnaseq_form_validate()
 * @see biosoftware_bench_rnaseq_form_submit()
 * @ingroup forms
 *
 * @param $form_state
 *   hook_form data
 * @param $sample_data
 *   Boolean, if TRUE then fill form with sample data
 */
function biosoftware_bench_rnaseq_form($form_state) {
  $software_options = biosoftware_bench_generate_options('rnaseq');
  $form             = array();
  $dbs              = biosoftware_bench_get_user_datasets('rnaseq');
  //returns in the form of ['dbtype']['group_name']['id']]=['name'];
  $option_form = array();
  if (empty($dbs)) {
    global $user;
    $roles_str = implode(', ', $user->roles);
    drupal_set_message(t("I'm sorry, but the system administrator has not authorized your role (%roles_str) to have access to at least 1 dataset. Please let them know.", array('%roles_str' => $roles_str)), 'error');
    return FALSE;
  }

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

    if ($data['type'] !== 'textfield' && !empty($data['options'])) {
      $option_form[$name]['#options'] = $data['options'];
      if ($data['type'] !== 'checkboxes') {
        $option_form[$name]['#default_value'] = current($data['options']);
      }
    }
  }
  $file_upload_desc = t('Please give a text file, not a MS-Word or other document, you can upload up to %m Mb.', array('%m' => gmod_dbsf_get_add_var('biosoftware_bench_upload_size')));
  $form['program'] = array(
    '#type' => 'fieldset',
    '#title' => 'RNA-seq expression query',
      'query_sequence' => array(
        '#type' => 'textarea',
        '#title' => t('Enter query sequence in simple text or FASTA format'),
        '#default_value' => '',
      ),
      'query_file' => array(
        '#type' => 'file',
        '#title' => t('or upload query sequence in FASTA format'),
        '#description' => $file_upload_desc,
        '#tree' => FALSE,
      ),
      'subject_db' => array(
        '#type' => 'select',
        '#title' => t('Dataset(s)'),
        '#description' => t('multiples allowed (use option/control keys; likewise to deselect).'),
        '#size' => 10,
        '#multiple' => '1',
        '#options' => $dbs,
      ),
      'advanced parameters' => array(
        '#type' => 'fieldset',
        '#title' => 'Advance Search Parameters',
        '#collapsible' => TRUE,
        '#collapsed' => TRUE,
        $option_form,
      ),
  );
  $form['buttons'] = array(
    '#weight' => 10,
    'db_submit' => array(
      '#type' => 'submit',
      '#value' => t('Run search'),
    ),
    'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'button',
      '#value' => 'Reset data',
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 2,
    ),
  );
  $form['#attributes']['enctype'] = "multipart/form-data";
  return $form;
}

/**
 * Implements hook_form_validate()
 *
 * @see biosoftware_bench_rnaseq_form()
 * @see biosoftware_bench_rnaseq_form_submit()
 *
 * @param $form
 *   hook_form data
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_rnaseq_form_validate($form, &$form_state) {
  foreach ($form_state['values'] as $key => $value) {
    if (is_array($value)) {
      continue;
    }
    elseif (empty($value)) {
      unset($form_state['values'][$key]);
    }
    elseif (preg_match('/^(\w+)-dataset-(\d+)/', $value, $matches)) {
      //dpm($matches);
      $form_state['values']['program'][$matches[1]]['subject_db'][] = $matches[2];
      unset($form_state['values'][$key]);
    }
  }

  $data = $form_state['values']['program'];
    if (!empty($data['query_sequence'])) {
      if (!empty($_FILES['files']['tmp_name']['query_file'])) {
        form_set_error("program][query_sequence", t('It seems you both typed a query sequence and uploaded one.'));
        return FALSE;
      }
      if (empty($data['subject_db'])) {
        form_set_error("program][query_sequence", t('It seems you provided a query sequence but no subject database.'));
        return FALSE;
      }
        $result = gmod_dbsf_validate_seq_dna($data['query_sequence']);
        if ($result === FALSE || is_numeric($result)) {
          form_set_error("program][query_sequence", t('Sorry your sequence does not seem to be a valid DNA sequence (%result errors).', array('%result' => $result)));
          return FALSE;
        }
        else {
          $form_state['values']['program']['query_sequence'] = $result;
        }

    }elseif (!empty($_FILES['files']['tmp_name']['query_file'])) {
      //report an error if no database has been selected. and query is uploaded
      if (empty($data['subject_db'])) {
        form_set_error("program][query_sequence", t('It seems you provided a query but no database database (or file).'));
        return FALSE;
      }
    }
    // now check for subject and make sure that if a query is selected, then subject is too.
    if (!empty($data['subject_db'])) {
      if (empty($data['query_sequence']) && empty($_FILES['files']['tmp_name']['query_file'])) {
        form_set_error("program][subject_db", t('It seems you provided a subject database but no query sequence (or file).'));
        return FALSE;
      }
    }
}

/**
 * Implements hook_form_submit()
 *
 * @see biosoftware_bench_rnaseq_form()
 * @see biosoftware_bench_rnaseq_form_validate()
 *
 * @param $form
 *   hook_form data
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_rnaseq_form_submit($form, &$form_state) {
  // Maximum Size in Mb of uploaded (subject or query) files.
  $file_size_allowed = gmod_dbsf_get_add_var('biosoftware_bench_upload_size');
  $validators_file   = array('file_validate_size' => array($file_size_allowed * 1048576));
  $tmpdir            = file_directory_temp();
  // all files, same timestamp
  $timestamp      = time();
  $sessionid      = session_id();
  $dirpath        = file_create_path() .'/bench';
  $data = $form_state['values']['program'];
  // BATCH API
  $operations = array();
  // for future use
  $save_array = array();
  // the form_uid will allow users to have multiple windows open.
  $form_uid = gmod_dbsf_create_uid($sessionid, $timestamp, 'rnaseq');
  $batch_file_data = array();

    $subject_dbs = array();
    $file_type = 'nucleotide';
    $verify = array();
    $uid = gmod_dbsf_create_uid($sessionid, $timestamp, 'rnaseq');
    // process user data
    // first, the query. check if file or text
    if (isset($_FILES['files']) && !empty($_FILES['files']['tmp_name']['query_file'])) {
      $file = file_save_upload('query_file', $validators_file, $tmpdir, FILE_EXISTS_RENAME);
      if (empty($file)) {
        form_set_error('query_file', t('Sorry your file was not saved. Maybe it is too large (>%file_size_allowed Mb)? Otherwise, '. l('contact', 'contact') .' the administrator  (quote %uid).', array('%uid' => $uid, '%file_size_allowed' => $file_size_allowed)));
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
    // now subject which will build the -d variable

    $file_type = 'nucleotide';
    if (!empty($data['subject_db'])) {
      $result = array();
      if (is_array($data['subject_db'])) {
        foreach ($data['subject_db'] as $db) {
          $verify['rundb']++;
          $result[$db] = biosoftware_bench_get_dataset_path($db);
        }
      }
      else {
        $verify['rundb']++;
        $result[$db] = biosoftware_bench_get_dataset_path($db);
      }
      foreach ($result as $db_data) {
        $subject_dbs[] = $db_data;
      }
    }
    //if it is going to run (both query and subject are ok)
    if (!empty($verify['runq']) && $verify['runq'] == 1 && !empty($verify['rundb'])) {
      $par .= biosoftware_bench_rnaseq_software_pars($data['advanced parameters'][0]);
      $save_array[$form_uid][$uid]['par'] = $par;
      $save_array[$form_uid][$uid]['sub'] = $subject_dbs;
    }

  // ended each algorithm
  $operations[] = array('gmod_dbsf_batch_upload_fasta', array($batch_file_data));
  $operations[] = array('gmod_dbsf_batch_save_data', array($save_array));
  $batch        = array(
    'title' => t('Preparing data needed for rnaseq jobs...'),
    'operations' => $operations,
    'init_message' => t('Starting rnaseq submission...'),
    'progress_message' => t('@remaining operations remaining...'),
    'error_message' => t('Your rnaseq submission encountered an error.'),
    'finished' => 'biosoftware_bench_rnaseq_batch_finished',
    'file' => drupal_get_path('module', 'biosoftware_bench') .'/includes/biosoftware_bench_rnaseq.inc',
    'progressive'=> TRUE,
  );
  $redirect_url=url('bench/rnaseq_result',
  array(
  'query'=>array('submission_uid'=>$form_uid),
  'alias'=>TRUE,'external'=>TRUE,
  )
  );
  $_REQUEST['destination']=$redirect_url;
  batch_set($batch);
  // We MUST do a batch_process, otherwise the form api attempts to load the entire
  // $form into $batch (as a variable and as a db) resulting in a 450 Mb variable
  batch_process();
}


/**
 * Prepares the rnaseq software cmdline parameters
 *
 * They will then be return to a file for parsing by the daemon
 * @see gmod_dbsf_parameter_daemon
 * Use this as a template for plugins
 *
 * @param $data
 *   Associative array with rnaseq option key/values for parsing into a
 *   parameter string
 * @param $algorithm
 *   The rnaseq algorithm used
 *   TODO: Integrate it into $data
 *
 * @return
 *   A string with all the parameters, suitable to  be used with cmdline
 *   interface. FALSE on failure.
 */
function biosoftware_bench_rnaseq_software_pars($data) {
  if (empty($data)) {
    return FALSE;
  }
  $par = '';
  foreach ($data as $name => $option) {
    if (empty($option)) {
      continue;
    }
  return $par;
}

/**
 * Implements hook_batch_finish
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
function biosoftware_bench_rnaseq_batch_finished($success, $results, $operations) {
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
      $to_store .= "algorithm;id=rnaseq;$uid\n";
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
      $par = $values['par'];
      $more_data = array();
      if (!empty($values['sub']) || !empty($results[$uid]['subject'])) {
        $subject_dbs = array_merge((array)$values['sub'], (array)$results[$uid]['subject']);
        $subject_dbs_str = trim(implode(' ', $subject_dbs));
        $par .= " -d '". $subject_dbs_str ."'";
        $more_data = array('database' => $subject_dbs);
      }
      gmod_dbsf_parameter_daemon($uid, $dirpath, $par, 'rnaseq', $more_data);
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
 * Parses a rnaseq XML file.
 *
 * @param $resultfile
 *   Full path to rnaseq XML file
 *
 * @return
 *   An associative array of the parsed XML, suitable for downstream processing
 */
function biosoftware_bench_rnaseq_result_xml($outfile = NULL) {
  $result_array = array();
  if (!file_exists($outfile)) {
    drupal_set_message(t('Failed to find %resultfile.', array('%resultfile' => $outfile)), 'error');
    return FALSE;
  }
  $xml_in    = simplexml_load_file($outfile, 'SimpleXMLElement');
  $algorithm = (string)$xml_in->rnaseqOutput_program;
  $xml       = $xml_in->rnaseqOutput_iterations;
  //is either an object (one query) or an array. so convert to array.
  $array_iterations = get_object_vars($xml);
  $array_iterations = $array_iterations['Iteration'];
  if (is_object($array_iterations)) {
    $array_iterations = array(0 => $array_iterations);
  }
  $query_rank;
  foreach ($array_iterations as $iteration_obj) {
    $query_rank++;
    $query_id           = '';
    $name_of_hit        = '';
    $accession_of_hit   = '';
    $average_bit_score  = '';
    $top_bit_score      = '';
    $average_evalue     = '';
    $lowest_evalue      = '';
    $average_identity   = '';
    $average_similarity = '';
    // we need array-ize to get query_id
    $iteration_data_array = get_object_vars($iteration_obj);

    $query_id = $iteration_data_array['Iteration_query-ID'] .' '. $iteration_data_array['Iteration_query-def'];
    $query_id = preg_replace('/\s*lcl\|\d+_\d+\s*/', '', $query_id);
    $query_id = str_replace('No definition line found', '', $query_id);
    if (empty($query_id)) {
      $query_id = "User query $query_rank";
    }
    if (!empty($iteration_data_array['Iteration_message'])) {
      drupal_set_message(t('%a : %b for %query_id', array('%query_id' => $query_id, '%a' => $algorithm, '%b' => $iteration_data_array['Iteration_message'])), 'warning');
    }
    $hit_iteration_hits = $iteration_data_array['Iteration_hits'];
    if (empty($hit_iteration_hits)) {
      continue;
    }

    $hits_array = get_object_vars($hit_iteration_hits);
    //dvm($hits_array);
    $hits_array = $hits_array['Hit'];
    if (is_object($hits_array)) {
      $hits_array = array(0 => $hits_array);
    }
    //dvm($hits_array);

    foreach ($hits_array as $rank0 => $hit_object) {

      // cast objects as strings to avoid extra lines of code
      $name_of_hit        = (string)$hit_object->Hit_id;
      $accession_of_hit   = (string)$hit_object->Hit_accession;
      $description_of_hit = (string)$hit_object->Hit_def;
      if ($description_of_hit == 'No definition line found') {
        $description_of_hit = '';
      }
      $hsps_array = get_object_vars($hit_object->Hit_hsps);
      // cannot get Hit_hsps->Hsp directly, because if HSP is an array, then get_object_vars only returns first element of array
      // so get Hit_hsps only. Then retrieve second element of resulting array (first is the XML key)
      $hsps_array = $hsps_array['Hsp'];
      // result is sometimes an array (multiple HSPs) and sometimes an object (one HSP)
      // So check
      if (is_object($hsps_array)) {
        // then one HSP and we got first HSP object directly. Convert to an array with one key.
        $hsps_array = array(0 => $hsps_array);
      }
      //$hsps_array is now an array of objects, each object is an hsp, key is hsp_rank.
      $total_bit_score  = 0;
      $total_evalue     = 0;
      $hsp_number       = 0;
      $total_identity   = 0;
      $total_similarity = 0;
      foreach ($hsps_array as $hsp_rank => $hsp_obj) {
        $hsp_data = get_object_vars($hsp_obj);
        //dpm($hsp_data);
        $bit_score  = $hsp_data['Hsp_bit-score'];
        $aln_length = $hsp_data['Hsp_align-len'];
        $evalue     = $hsp_data['Hsp_evalue'];
        if (!empty($aln_length)) {
          $identity = $hsp_data['Hsp_identity'] / $aln_length;
          $similarity = $hsp_data['Hsp_positive'] / $aln_length;
        }
        $total_bit_score += $bit_score;
        $total_evalue += $evalue;
        $total_identity += $identity;
        $total_similarity += $similarity;
        $hsp_number++;
        if (empty($lowest_evalue) || $evalue < $lowest_evalue) {
          $lowest_evalue = $evalue;
        }
        if (empty($top_bit_score) || $bit_score > $top_bit_score) {
          $top_bit_score = $bit_score;
        }
      }
      $lowest_evalue = sprintf('%.2e', $lowest_evalue);
      $average_bit_score = $hsp_number ? sprintf('%.2f', $total_bit_score / $hsp_number) : 'NA';
      $average_evalue = $hsp_number ? sprintf('%.2e', $total_evalue / $hsp_number) : 'NA';
      $top_bit_score = sprintf('%.2f', $top_bit_score);
      $average_identity = $hsp_number ? sprintf('%.2f %%', $total_identity * 100 / $hsp_number) : 'NA';
      $average_similarity = $hsp_number ? sprintf('%.2f %%', $total_similarity * 100 / $hsp_number) : 'NA';
      $result_array[$query_id][$rank0 + 1] = array(
        'name_of_hit' => $name_of_hit,
        'accession_of_hit' => $accession_of_hit,
        'average_bit_score' => $average_bit_score,
        'top_bit_score' => $top_bit_score,
        'average_evalue' => $average_evalue,
        'lowest_evalue' => $lowest_evalue,
        'average_identity' => $average_identity,
        'average_similarity' => $average_similarity,
      );
    }
  }
  return $result_array;
}

/**
 * Produces HTML results for a rnaseq
 *
 * @see biosoftware_bench_rnaseq_result_form()
 *
 * @param $submission_uid
 *   The UID relating to the submission. This is not the same as the individual
 *   rnaseq algorithm UID as users can perform multiple algorithm searches from
 *   one form
 *
 * @return
 *   HTML of results, include a table and graph of hits.
 */
function biosoftware_bench_rnaseq_result_page($submission_uid = NULL, $level = 1) {
  $software = 'rnaseq';
  if (empty($submission_uid)) {
    $submission_uid = $_GET['submission_uid'];
    if (empty($submission_uid)) {
      return 'No submission data provided';
    }
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

  $page_tabs = array(
    'rnaseq_results' => array(
      '#type' => 'tabset',
    ),
  );
  $return_html = '';
  foreach ($data as $algorithm => $uid) {
    // add base_path to allow downloading but must not be used for anything else...
    $resultfile = $dirpath . $uid .'.rnaseq.output';
    $graph      = $resultfile .'.graph.png';
    $errorfile  = $dirpath . $uid .'.rnaseq.error';

    $content = t('Your report is not ready yet. Please wait or ') . l('try again', 'bench/rnaseq_result/'. $submission_uid
    ) . t(' in a few moments or bookmark this page.');
    if (file_exists($graph)) {
      $results = drupal_get_form('biosoftware_bench_rnaseq_result_form', $uid, $resultfile, $level, $submission_uid);
      if (strlen($results) < 600) {
        $results = "No (more) hits found for any of your queries using $algorithm. Perhaps you can try a different algorithm or database?";
      }
      $content = '<p>'. l(
        theme_image($graph, 'rnaseq results of first ten queries', 'rnaseq results graphical overview', array(), FALSE)
        , base_path() . $resultfile .'.html', array('html' => TRUE, 'external' => TRUE, 'attributes' => array('target' => '_blank'))
      ) .'</p>';

      $content .= $results;
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
      drupal_set_html_head('<meta http-equiv="refresh" content="15;url='. $base_url .'?q=/bench/rnaseq_result/'. $submission_uid .'" />');
    }
    $page_tabs['rnaseq_results'][$algorithm] = array(
      '#type' => 'tabpage',
      '#title' => 'Results for '. $algorithm,
      '#content' => $content,
    );
  }
  $return_text = '<p>Do a '. l('new search', "bench/$software") ." or use your browser's back button</p>";
  $return_text .= '<br><p>If you want, you can bookmark this page: your results will be stored on the server at least until '. date("F j, Y, g:i a", time() + 604800) .'. Let us '. l('know', 'contact') .' if you would like to keep them for longer.</p>';
  $return_text .= tabs_render($page_tabs);
  return $return_text;
}

/**
 * Implements hook_form()
 *
 * Used to produce a table of results with checkboxes which we can then
 * select for processing (e.g. downloading). This form is themed.
 * @see biosoftware_bench_software_result_form()
 * @see biosoftware_bench_rnaseq_result_page()
 * @ingroup forms
 *
 * @param $form_state
 *   Form API variable
 * @param $uid
 *   Unique Identifier of algorithm's rnaseq search
 * @param $resultfile
 *   Full path to rnaseq results as an XML file
 */
function biosoftware_bench_rnaseq_result_form($form_state, $uid, $resultfile, $level = 1, $suid = NULL) {
  $form = array();
  if (empty($resultfile)) {
    return $form;
  }
  if (empty($level)) {
    $level = 1;
  }
  $array_from_xml = array();
  $array_from_xml = biosoftware_bench_rnaseq_result_xml($resultfile);
  //$xml=biosoftware_bench_rnaseq_result_xml($resultfile);dpm($xml);
  // all data
  $table_row = array();
  // for checkboxes to get hit
  $row_check_data = array();
  if (empty($array_from_xml) || is_array($array_from_xml) === FALSE) {
    return $form;
  }
  $weight_series = 0;
  $query_rank    = 0;
  $total_queries = count($array_from_xml);
  $table_caption = '';
  foreach ($array_from_xml as $query_id => $hit_array) {
    $query_rank++;
    if ($query_rank != $level) {
      continue;
    }
    $table_caption = "Overview of significant results for $query_id (query $level)";
    foreach ($hit_array as $hit_rank => $data) {
      $hit_id    = $data['accession_of_hit'] ? $data['accession_of_hit'] : $data['name_of_hit'];
      $bit_score = $data['average_bit_score'] .' ('. $data['top_bit_score'] .')';
      $evalue    = $data['average_evalue'] .' ('. $data['lowest_evalue'] .')';
      $identity  = $data['average_identity'] .' ('. $data['average_similarity'] .')';
      $links     = l(
        theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/xml.png', 'Report as XML', 'Report as XML (recommended for parsing with other software)'), base_path() . $resultfile,
        array(
          'html' => TRUE,
          'external' => TRUE,
          'attributes' => array('target' => '_blank'),
        )
      ) . l(
        theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/html.png', 'Report as HTML', 'Report as HTML'), base_path() ."$resultfile.Query$query_rank.html",
        array(
          'html' => TRUE,
          'external' => TRUE,
          'attributes' => array('target' => '_blank'),
        )
      ) . l(
        theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/txt.png', 'Report as text', 'Report as text'), base_path() ."$resultfile.Query$query_rank.txt",
        array(
          'html' => TRUE,
          'external' => TRUE,
          'attributes' => array('target' => '_blank'),
        )
      );
      $table_row[] = array($query_id, $hit_id, $bit_score, $evalue, $identity, $links);
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
  $base              = $base_url .'/bench/rnaseq_result/'. $suid;
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


  $table_header = array('Query ID', 'Hit ID', 'Average bit score (top)', 'Average e-value (lowest)', 'Average identity (av. similarity)', 'Links');
  $table_attributes = array('class' => 'padded-table sortable');

  $table_array = array(
    '#links' => $pager,
    'header' => $table_header,
    'caption' => $table_caption,
    'attributes' => $table_attributes,
    'data' => $table_row,
  );
  $form['resultfile'] = array('#type' => 'value', '#value' => $resultfile);
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
  $form['#action'] = url("bench/get_sequences/$uid");
  return $form;
}

