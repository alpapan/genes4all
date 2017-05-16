<?php
// $Id$

/**
 * @file
 *
 * Please use this file as a template to make new modules
 * You will need to search for lines which have
 * PLUGIN_EDIT_HERE and follow the instructions
 * you will also need to rename the functions. I suggest you
 * pick a name for your software (e.g. iprscan) and
 * search and replace ssaha2 with e.g. iprscan. then continue editing.
 */

/**
 * @file
 * Universal Search Server Module plugin
 * A Drupal module to implement SSAHA2 usage via the Universal Search Server
 *
 * A SSAHA dataset must have the following to be allowed
 * 1) The actual raw fasta file (eg dmel_genome_5.16)
 * 2) formatted with index_fasta.pl (in scripts; a BioPerl script)
 *    e.g. index_fasta.pl -d dmel_genome_5.16
 * 3) hash table build with the ssaha2Build -rtype abi as -save dmel_genome_5.16_sanger
 * 4) hash table build with the ssaha2Build -rtype 454 as -save dmel_genome_5.16_454
 * 5) hash table build with the ssaha2Build -rtype solexa as -save dmel_genome_5.16_solexa
 * 6) ALL the *.head files produced by ssaha2Build MUST be compressed; use gzip -1 as it is fast.
 * so that for a dataset 'dmel' we have the files
 * - dmel, dmel.index
 * - dmel_sanger.base,dmel_sanger.body,dmel_sanger.head.gz,dmel_sanger.name,dmel_sanger.size
 * - dmel_454.base,dmel_454.body,dmel_454.head.gz,dmel_454.name,dmel_454.size,
 * - dmel_solexa.basedmel_solexa.body,dmel_solexa.head.gz,dmel_solexa.name,dmel_solexa.size
 * Symbolic links are also allowed but must be set in your Apache configuration file as allowed ("follow symlinks" etc)
 * NB The run-user of biosoftware_bench_daemon.pl must have permission to write in the directory where the .head.gz file resides
 *
 * Drupal module developed by Alexie Papanicolaou.
 * University of Exeter
 * @see http://insectacentral.org/
 *
 * This module was built using SSAHA2 v2.4.1
 *
 * SSAHA2 is a package combining SSAHA with cross_match developed
 *   by Phil Green at the University of Washington.
 *   Reference: Ning Z, Cox AJ, Mullikin JC.
 *          SSAHA: a fast search method for large DNAdatabases.
 *          Genome Res. 2001 Oct;11(10):1725-9.
 */
function biosoftware_bench_ssaha2_menu() {
  $plugin_path = drupal_get_path('module', 'biosoftware_bench') .'/includes/';
  $items = array();

  /*PLUGIN_EDIT_HERE
   * We are making a menu e.g. admin/bench/ssaha2 -> http://localhost/admin/bench/ssaha2
   * need to change title, description
   * access arguments and page_callback
   * The access arguments defines user permissions to access it this menu
   * the page_callback is the name of the function called when this menu is visited
   *
   */


  $items['admin/bench/ssaha2'] = array(
    'file' => 'includes/biosoftware_bench_ssaha2.inc',
    'title' => 'SSAHA2 Server administration',
    'page callback' => 'biosoftware_bench_admin_ssaha2_page',
    'access arguments' => array('Administrate SSAHA2 Server'),
    'description' => 'Configure SSAHA2 specific settings',
    'type' => MENU_CALLBACK,
  );
  $items['bench/ssaha2'] = array(
    'file' => 'includes/biosoftware_bench_ssaha2.inc',
    'title' => 'biosoftware_bench SSAHA2 Server',
    'page callback' => 'biosoftware_bench_ssaha2_page',
    'access arguments' => array('Access SSAHA2 Server'),
    'type' => MENU_CALLBACK,
  );
  $items['bench/ssaha2_result'] = array(
    'file' => 'includes/biosoftware_bench_ssaha2.inc',
    'title' => 'biosoftware_bench SSAHA2 Results',
    'page callback' => 'biosoftware_bench_ssaha2_result_page',
    'page arguments' => array(2, 3),
    'access arguments' => array('Access SSAHA2 Server'),
    'type' => MENU_CALLBACK,
  );
  return $items;
}

function biosoftware_bench_software_ssaha2_core_settings_form($form_state) {
  /*PLUGIN_EDIT_HERE
   * $software is the name of your software
   */
  $software = 'ssaha2';
  $form = array();
  $res = db_fetch_array(db_query("SELECT active from {gmod_dbsf_software} where uniquename='%s'", $software));
  $select_software_setting = 'SELECT value FROM {gmod_dbsf_softwareprop} where software_id=' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ' AND rank=0 AND type_id=' . "(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')";
  if ($res['active'] == 'f') {
    return $form;
  }

  /*PLUGIN_EDIT_HERE
   * SSAHA has two executables, one for formatting databases, and one for executing
   * In reality we are not really using the formating executable but we might in the future
   * your software can have multiple executables (e.g. see annot8r)
   */


  $core_exec = db_fetch_array(db_query($select_software_setting, $software, 'executable'));
  $format_exec = db_fetch_array(db_query($select_software_setting, 'ssaha2Build', 'executable'));
  $external_server = gmod_dbsf_get_add_var('biosoftware_bench_external_server');

  $form['settings'] = array(
    '#type' => 'fieldset',
    '#title' => $software . t(' core settings'),
    '#description' => t('Please set some important settings for %soft', array('%soft' => $software)),
    'submit' => array(
      '#type' => 'submit',
      '#value' => "Save core $software settings",
      '#weight' => 5,
    ),
    'core_executable' => array(
      '#type' => 'textfield',
      '#title' => t('Core executable path'),
      '#description' => t('Please provide the <strong>full path</strong> to the <strong>%soft</strong> executable on the server, e.g. /usr/bin/%soft .', array('%soft' => $software)),
      '#required' => TRUE,
      '#default_value' => $core_exec['value'],
    ),
    /*PLUGIN_EDIT_HERE
   * For each executable we make an addition to the form array
   * the core executable is above and we don't need to change it.
   * You might want to edit/delete the format_executable array below
   */
    'format_executable' => array(
      '#type' => 'textfield',
      '#title' => t('Executable path for database formating'),
      '#description' => t('Please provide the <strong>full path</strong> to the <strong>ssaha2Build</strong> executable on the server, e.g. /usr/bin/ssaha2Build .'),
      '#required' => TRUE,
      '#default_value' => $format_exec['value'],
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

function biosoftware_bench_software_ssaha2_core_settings_form_validate($form, &$form_state) {
  $core_exec = escapeshellcmd(trim($form_state['values']['core_executable']));
  $format_exec = escapeshellcmd(trim($form_state['values']['format_executable']));
  if (!empty($form_state['values']['external_server'])){
    gmod_dbsf_get_add_var('biosoftware_bench_external_server','TRUE');
  }elseif (gmod_dbsf_get_add_var('biosoftware_bench_external_server')){
    gmod_dbsf_get_add_var('biosoftware_bench_external_server','DELETE');
  }
  if (!empty($core_exec)) {
    if (!file_exists($core_exec) && empty($form_state['values']['external_server'])) {
      form_set_error('core_executable', t('Cannot find/access executable %core_exec on the server.', array('%core_exec' => $core_exec)));
    }
    elseif (!is_executable($core_exec) && empty($form_state['values']['external_server'])) {
      form_set_error('core_executable', t('Cannot execute the %f program.', array('%f' => $core_exec)));
    }
    $form_state['values']['core_executable'] = $core_exec;
  }

  if (!empty($format_exec)) {
    if (!file_exists($format_exec) && empty($form_state['values']['external_server'])) {
      form_set_error('format_executable', t('Cannot find/access executable %format_exec on the server.', array('%format_exec' => $format_exec)));
    }
    elseif (!is_executable($format_exec) && empty($form_state['values']['external_server'])) {
      form_set_error('format_executable', t('Cannot execute the %f program.', array('%f' => $format_exec)));
    }
    $form_state['values']['format_executable'] = $format_exec;
  }
}

function biosoftware_bench_software_ssaha2_core_settings_form_submit($form, &$form_state) {
  /*PLUGIN_EDIT_HERE
   *edit software name
   */


  $software = 'ssaha2';

  $store_software_setting = 'INSERT INTO {gmod_dbsf_softwareprop} (software_id,type_id,rank,value) VALUES (' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ",(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')" . ",0,'%s')";
  $delete_software_setting = 'DELETE FROM {gmod_dbsf_softwareprop} where software_id=' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ' AND rank=0 AND type_id=' . "(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')";
  $select_software_setting = 'SELECT value FROM {gmod_dbsf_softwareprop} where software_id=' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ' AND rank=0 AND type_id=' . "(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')";

  $core_exec = db_fetch_array(db_query($select_software_setting, $software, 'executable'));
  $format_exec = db_fetch_array(db_query($select_software_setting, 'ssaha2Build', 'executable'));

  $user_core_exec = $form_state['values']['core_executable'];
  $user_format_exec = $form_state['values']['format_executable'];

  if (!empty($user_core_exec) && $user_core_exec != $core_exec['value']) {
    db_query($delete_software_setting, $software, 'executable');
    db_query($store_software_setting, $software, 'executable', $user_core_exec);
    drupal_set_message(t('%user_core_exec is now the %software executable.', array('%user_core_exec' => $user_core_exec, '%software' => $software)), 'warning');
  }
  if (!empty($user_format_exec) && $user_format_exec != $format_exec['value']) {
    db_query($delete_software_setting, 'ssaha2Build', 'executable');
    db_query($store_software_setting, 'ssaha2Build', 'executable', $user_format_exec);
    drupal_set_message(t('%user_format_exec is now the ssaha2Build executable.', array('%user_format_exec' => $user_format_exec)), 'warning');
  }
}

function biosoftware_bench_admin_ssaha2_page($software = 'ssaha2') {
  /*PLUGIN_EDIT_HERE
   * edit software name above in function
   */



  // Include biosoftware_bench admin function file
  require_once(drupal_get_path('module', 'biosoftware_bench') .'/includes/biosoftware_bench_admin.inc');
  $check_active = biosoftware_bench_check_software_active($software);
  if (empty($check_active)) {
    drupal_set_message(t("I'm sorry, but %software does not seem to have been activated yet.<br>", array('%software' => $software)) .
      l(t('See the software settings page'), 'admin/bench/software'), 'error'
    );
    return FALSE;
  }
  $path_form = drupal_get_form('biosoftware_bench_software_ssaha2_core_settings_form');
  $database_link_form = drupal_get_form('biosoftware_bench_admin_link_dataset_software_form', $software);
  if (strlen($database_link_form) < 600) {
    $database_link_form = "You have no datasets defined to link them to $software. Please see ". l('the dataset administration page', 'admin/bench/software', array('fragment' => 'available-datasets'));
  }
  $table = biosoftware_bench_admin_database_table($software);

  $page_tabs = array(
    'settings' => array(
      '#type' => 'tabset',
      'core' => array(
        '#type' => 'tabpage',
        '#title' => 'Core settings',
        '#content' => $path_form,
        '#weight' => 0,
      ),
      'variables' => array(
        '#type' => 'tabpage',
        '#title' => 'Link datasets to '. strtoupper($software),
        '#content' => $database_link_form . $table,
        '#weight' => 1,
      ),
    ),
  );
  $return_text = tabs_render($page_tabs);
  return $return_text;
}

/**
 * Check if the ssaha2 software is setup correctly
 *
 * @param $software
 * software name
 *
 * @return boolean
 * TRUE if everything is ok, a string signifying what is wrong otherwise
 */
function biosoftware_bench_check_ssaha2($software) {
  $sql = 'SELECT count(distinct rcv.cvterm_id) from {gmod_dbsf_resource_cvterm} as rcv ' . ' JOIN {gmod_dbsf_software_resource} sr ON sr.resource_id=rcv.resource_id ' . " WHERE sr.software_id=(SELECT software_id from {gmod_dbsf_software} where uniquename='%s') " . " AND rcv.cvterm_id IN (SELECT cvterm_id FROM {gmod_dbsf_cvterm} where cv_id=(SELECT cv_id from {gmod_dbsf_cv} where name='dataset_type'))";

  /*PLUGIN_EDIT_HERE
   * This function is generalized but can be customized
   * you probably want to change (2) and leave (1) and (3) as they are
   * currently it checks if:
   */



  // 1. is the software activated?
  $check_active = biosoftware_bench_check_software_active($software);
  if (empty($check_active)) {
    return 'active';
  }

  //2. are the two executable variables set?
  $check = 0;
  $software_prop = gmod_dbsf_get_softwareprop(NULL, TRUE);
  if (!empty($software_prop)) {
    foreach ($software_prop as $name => $data) {
      if (in_array($name, array('ssaha2', 'ssaha2Build'))) {
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

  //3. is there are at least one dataset defined?
  $count = db_fetch_array(db_query($sql, $software));
  if ($count['count'] < 1) {
    return 'datasets';
  }
  return TRUE;
}

function biosoftware_bench_ssaha2_page() {
  /*PLUGIN_EDIT_HERE
   * software name and rename functions below
   */


  $software = 'ssaha2';
  $setup = biosoftware_bench_check_ssaha2($software);
  if ($setup === TRUE) {
    return drupal_get_form('biosoftware_bench_ssaha2_form');
  }

  // Check has failed. Report to the user/admin what the error is
  else {
    $return_text = t("<p>I'm sorry but your %software software has not been setup properly.</p><ul>", array('%software' => $software));
    switch ($setup) {
      case 'datasets':
        /*PLUGIN_EDIT_HERE
         * edit menu link below
         */

        $return_text .= '<li>Not enough datasets have been defined. Please define at least one dataset at the '. l('dataset registration page', 'admin/bench/dataset', array('fragment' => 'available-datasets')) . t(' and then register them with %software at the ', array('%software' => $software)) . l("$software settings page", 'admin/bench/ssaha2');
        break;

      case 'active':
        $return_text .= "<li>$software has not been activated. Please activate it at the ". l('software settings page', 'admin/bench/software');
        break;

      case 'variables':
        /*PLUGIN_EDIT_HERE
         * edit fragment link below
         */


        $return_text .= "<li>$software variable paths have not been defined. Please define them at at the ". l("$software settings page", 'admin/bench/software', array('fragment' => 'ssaha2'));
        break;
    }
    drupal_set_message(t('!t', array('!t' => $return_text .'</ul>')), 'error');
    return ' ';
  }
}

function biosoftware_bench_ssaha2_form($form_state) {
  $form = array();

  /*PLUGIN_EDIT_HERE
   * set software name
   */


  $software         = 'ssaha2';
  $software_options = biosoftware_bench_generate_options($software);
  $dbs              = biosoftware_bench_get_user_datasets($software);

  /*PLUGIN_EDIT_HERE
   * here are getting the datasets. For ssaha they are all nucleotide
   */


  //returns in the form of ['dbtype']['group_name']['id']]=['name'];
  $nucleotide_dbs = (array)$dbs['nucleotide'];
  $option_form = array();
  if (empty($nucleotide_dbs)) {
    global $user;
    $roles_str = implode(', ', $user->roles);
    drupal_set_message(t("I'm sorry, but the system administrator has not authorized your role (%roles_str) to have access to at least 1 dataset with %software. Please let them know.", array('%roles_str' => $roles_str, '%software' => $software)), 'error');
    return FALSE;
  }
  if (!empty($software_options)) {
    foreach ($software_options as $name => $data) {
      $description = !empty($data['description']) ? t('@d', array('@d' => $data['description'])) : '';
      $title = str_replace($software, '', $name);
      //if name has a " - " then put anything after - into description
      if (strpos($title, ' - ')) {
        $array = explode(' - ', $title);
        $title = $array[0];
        $description .= $title .': '. $array[1];
      }
      $title = ucfirst(trim(str_replace('_', ' ', $title)));
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
    //PLUGIN_EDIT_HERE
    // some special cases go here
    $option_form['ssaha2_query_type']['#weight'] = -1;
    $option_form['ssaha2_query_type']['#default_value'] = 'Sanger ABI';
    $option_form['ssaha2_score']['#options'] = array('0' => 'Default') + $option_form['ssaha2_score']['#options'];
    $option_form['ssaha2_score']['#default_value'] = array('0' => 'Default');
    $option_form['ssaha2_cmatch']['#options'] = array('0' => 'Default') + $option_form['ssaha2_cmatch']['#options'];
    $option_form['ssaha2_cmatch']['#default_value'] = array('0' => 'Default');
    $option_form['ssaha2_ckmer']['#options'] = array('0' => 'Default') + $option_form['ssaha2_ckmer']['#options'];
    $option_form['ssaha2_ckmer']['#default_value'] = array('0' => 'Default');
  }
  //PLUGIN_EDIT_HERE
  // You may prepare a description HTML file. It will be parsed here
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


  /*PLUGIN_EDIT_HERE
   * This is the form. Customize title/description as you wish but probably leave
   * the array keys as they are
   */


  $form['program'] = array(
    '#type' => 'fieldset',
    '#title' => strtoupper($software),
    'query_sequence' => array(
      '#type' => 'textarea',
      '#title' => t('Enter query nucleotide sequence in simple text, FASTA or FASTQ format'),
    ),
    'query_file_'. $software => array(
      '#type' => 'file',
      '#title' => t('or upload query nucleotide sequence in FASTA or FASTQ format'),
      '#description' => t('Please give a text file, not a MS-Word or other document, you can upload up to %m Mb.', array('%m' => gmod_dbsf_get_add_var('biosoftware_bench_upload_size'))),
      // Needed because of drupal bug
      '#tree' => FALSE,
    ),
    'subject_db' => array(
      '#type' => 'select',
      '#title' => t('Dataset(s)'),
      '#description' => t('Please select one dataset.'),
      //PLUGIN_EDIT_HERE
      // SSAHA allows only one database. you can change it here
      '#size' => 1,
      '#multiple' => FALSE,
      '#options' => $nucleotide_dbs + array(0 => 'Use Upload'),
    ),
    'subject_file_'. $software => array(
      '#type' => 'file',
      '#title' => t('or upload your own database as a FASTA'),
      '#description' => t('Please give a FASTA text file, not a MS-Word or other document'),
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
      //'#tree'=>TRUE,
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
  $form['#attributes']['enctype'] = "multipart/form-data";
  return $form;
}

function biosoftware_bench_ssaha2_form_validate($form, &$form_state) {
  /*PLUGIN_EDIT_HERE
   * give software name
   */


  $software = 'ssaha2';
  if (empty($form_state['values']['query_sequence']) && empty($_FILES['files']['tmp_name']['query_file_'. $software])) {
    form_set_error("query_sequence", t('It seems you neither gave a query sequence nor uploaded one.'));
    return FALSE;
  }
  if (!empty($form_state['values']['query_sequence'])) {
    if (!empty($_FILES['files']['tmp_name']['query_file_'. $software])) {
      form_set_error("query_sequence", t('It seems you both typed a query sequence and uploaded one.'));
      return FALSE;
    }
    if (empty($form_state['values']['subject_db']) && empty($_FILES['files']['tmp_name']['subject_file_'. $software])) {
      form_set_error("query_sequence", t('It seems you provided a query sequence but no subject database (or file).'));
      return FALSE;
    }
  }
  elseif (!empty($_FILES['files']['tmp_name']['query_file_'. $software])) {
    //report an error if no database has been selected. and query is uploaded
    if (empty($form_state['values']['subject_db']) && empty($_FILES['files']['tmp_name']['subject_file_'. $software])) {
      form_set_error("query_sequence", t('It seems you provided a query but no database database (or file).'));
      return FALSE;
    }
  }
  // now check for subject and make sure that if a query is selected, then subject is too.
  if (!empty($form_state['values']['subject_db'])) {
    if (empty($form_state['values']['query_sequence']) && empty($_FILES['files']['tmp_name']['query_file_'. $software])) {
      form_set_error("subject_db", t('It seems you provided a subject database but no query sequence (or file).'));
      return FALSE;
    }
    /*PLUGIN_EDIT_HERE
     *You can add additional validation steps
     *e.g. here we preventing both upload and selection of subject databases.
     */


    if (!empty($_FILES['files']['tmp_name']['subject_file_'. $software])) {
      form_set_error("subject_db", t('It seems you are uploading and selecting a subject database. The %software software can only search one database at a time.', array('%software' => $software)));
    }
    //  No need to check if it is nucleootide or protein as it is a pre-formated db
  }
  elseif (!empty($_FILES['files']['tmp_name']['subject_file_'. $software])) {
    //report an error if no query has been selected.
    if (empty($form_state['values']['query_sequence']) && empty($_FILES['files']['tmp_name']['query_file_'. $software])) {
      form_set_error("subject_db", t('It seems you provided a subject database but no query sequence (or file).'));
      return FALSE;
    }
  }
}

/**
 * Submit the search form
 *
 * This form does all the work, including the batch processing
 */
function biosoftware_bench_ssaha2_form_submit($form, &$form_state) {
  /*PLUGIN_EDIT_HERE
   * give software name. you may want to edit file_validate_size below,
   * if you expect uploaded files to have a different size
   */


  $software = 'ssaha2';
  // Maximum Size in Mb of uploaded (subject or query) files.
  $file_size_allowed = gmod_dbsf_get_add_var('biosoftware_bench_upload_size');;

  $validators_file = array('file_validate_size' => array($file_size_allowed * 1000));
  $tmpdir = file_directory_temp();
  // all files, same timestamp
  $timestamp = time();
  $sessionid = session_id();
  $dirpath   = file_create_path() .'/bench';
  $data      = $form_state['values'];
  // BATCH API
  $operations = array();
  $save_array = array();
  // the form_uid will allow users to have multiple windows open.
  $form_uid        = gmod_dbsf_create_uid($sessionid, $timestamp, $software);
  $batch_file_data = array();
  $file_type       = 'nucleotide';
  $verify          = array();
  // in blast we have separate uids for each algorithm. not needed here
  $uid = $form_uid;

  // if you are making a loop (like in BLAST) then this needs to go inside
  $subject_dbs = array();
  if (isset($_FILES['files']) && !empty($_FILES['files']['tmp_name']['query_file_'. $software])) {
    $file = file_save_upload('query_file_'. $software, $validators_file, $tmpdir, FILE_EXISTS_RENAME);
    if (empty($file)) {
      form_set_error('query_file_'. $software, t('Sorry your file for <em>%software</em> was not saved. Maybe it is too large (>%file_size_allowed Mb)? Otherwise, '. l('contact', 'contact') .' the administrator  (quote %uid).', array('%uid' => $uid, '%file_size_allowed' => $file_size_allowed, '%software' => $software)));
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
  // now subject which will decide the database
  $file_type = 'nucleotide';
  // if subject has been provided by user, either format it and add it to -d or reject it.
  if (isset($_FILES['files']) && !empty($_FILES['files']['tmp_name']['subject_file_'. $software])) {
    $file = file_save_upload('subject_file_'. $software, $validators_file, $tmpdir, FILE_EXISTS_RENAME);
    if (empty($file)) {
      form_set_error('subject_file_'. $software, t('Sorry your file for <em>ssaha2</em> was not saved. Maybe it is too large (>%file_size_allowed Mb)? Otherwise, '. l('contact', 'contact') .' the administrator  (quote %uid).', array('%uid' => $uid, '%file_size_allowed' => $file_size_allowed)));
    }
    else {
      $verify['rundb']++;
      file_set_status($file, FILE_STATUS_TEMPORARY);
      $batch_file_data['infile'][$uid] = $file->filepath;
      $batch_file_data['outfile'][$uid] = $dirpath .'/'. $uid .'.subject';
      $batch_file_data['filetype'][$uid] = $file_type;
      $batch_file_data['format'][$uid] = TRUE;
    }
  }
  if (!empty($data['subject_db'])) {
    if (is_array($data['subject_db'])) {
      foreach ($data['subject_db'] as $db) {
        $verify['rundb']++;
        $subject_dbs[] = biosoftware_bench_get_dataset_path($db);
      }
    }
    else {
      $verify['rundb']++;
      $subject_dbs[] = biosoftware_bench_get_dataset_path($data['subject_db']);
    }
  }
  //check if it is going to run (i.e. both query and subject are ok)
  if (!empty($verify['runq']) && $verify['runq'] == 1 && !empty($verify['rundb'])) {
    /*PLUGIN_EDIT_HERE
     * a small function to convert option in the form to cmdline arguments
     */


    $par = gmod_dbsf_ssaha_parameter($data);
    $save_array[$form_uid][$uid]['par'] = $par;
    $save_array[$form_uid][$uid]['rtype'] = $data['ssaha2_query_type'];
    $save_array[$form_uid][$uid]['sub'] = $subject_dbs;
    $save_array[$form_uid][$uid]['algorithm'] = $software;
  }

  // ended each algorithm
  $operations[] = array('gmod_dbsf_batch_upload_fasta', array($batch_file_data));
  $operations[] = array('gmod_dbsf_batch_save_data', array($save_array));
  $batch        = array(
    'title' => t('Preparing data needed for %software...', array('%software' => strtoupper($software))),
    'operations' => $operations,
    'init_message' => t('Starting %software submission...', array('%software' => strtoupper($software))),
    'progress_message' => t('@remaining operations remaining...'),
    'error_message' => t('Your %software submission encountered an error.', array('%software' => strtoupper($software))),
    //PLUGIN_EDIT_HERE
    // Function name and file below,
    'finished' => 'biosoftware_bench_ssaha2_batch_finished',
    'file' => drupal_get_path('module', 'biosoftware_bench') .'/includes/biosoftware_bench_ssaha2.inc',
  );
  batch_set($batch);
  /*PLUGIN_EDIT_HERE
   * menu link below
   */


  $form_state['redirect'] = array('bench/ssaha2_result', "submission_uid=$form_uid");
  return "bench/ssaha2_result?submission_uid=$form_uid";
}

function gmod_dbsf_ssaha_parameter($data) {
  if (empty($data)) {
    return;
  }
  //PLUGIN_EDIT_HERE
  //edit this function to convert form options to cmd line arguments
  // to return
  $cmds = ' ';
  switch ($data['ssaha2_query_type']) {
    case 'Sanger ABI':
      $cmds .= ' -rtype abi';
      break;

    case '454';
    //454 -score 30 -sense 1 -cmatch 10 -ckmer 6
    $cmds .= ' -rtype 454';
    break;

  case 'Illumina-Solexa';
  //solexa -kmer 13 -skip 2 -seeds 2 -score 12 -cmatch 9 -ckmer 6
  $cmds .= ' -rtype solexa';
  break;
}

if (!empty($data['ssaha2_identity'])) {
$cmds .= ' -identity '. $data['ssaha2_identity'];
}
if (!empty($data['ssaha2_score'])) {
$cmds .= ' -score '. $data['ssaha2_score'];
}
if (!empty($data['ssaha2_cmatch'])) {
$cmds .= ' -cmatch '. $data['ssaha2_cmatch'];
}
if (!empty($data['ssaha2_ckmer'])) {
$cmds .= ' -ckmer '. $data['ssaha2_ckmer'];
}
if (!empty($data['ssaha2_cut'])) {
$cmds .= ' -cut '. $data['ssaha2_cut'];
}
if (!empty($data['ssaha2_depth'])) {
$cmds .= ' -depth '. $data['ssaha2_depth'];
}
if (!empty($data['ssaha2_best'])) {
$cmds .= ' -best '. $data['ssaha2_best'];
}

return $cmds;
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
function biosoftware_bench_ssaha2_batch_finished($success, $results, $operations) {
/*PLUGIN_EDIT_HERE
   * give software name
   */


$software = 'ssaha2';
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
  $par = $values['par'];
  $more_data = array();
  if (!empty($values['sub']) || !empty($results[$uid]['subject'])) {
    /*PLUGIN_EDIT_HERE
         * we now have to put the database in the cmd line. the commented out code
         * below is for BLAST. For ssaha we have only one database, which if it is
         * uploaded, it has to be hash_built (formatted) on the fly
         //$subject_dbs     = array_merge((array)$values['sub'], (array)$results[$uid]['subject']);
         //$subject_dbs_str = trim(implode(' ', $subject_dbs));
         //$par         .= " -d '". $subject_dbs_str ."'";
         *
         */


    $local_db = trim(implode(' ', (array)$values['sub']));
    $user_db = trim(implode(' ', (array)$results[$uid]['subject']));

    if (!empty($local_db) && !empty($user_db)) {
      //not supposed to happen!
      drupal_set_message(t('This is not supposed to happen: more than one subject database was requested.'), 'warning');
      return FALSE;
    }
    $db_used = !empty($local_db) ? $local_db : $user_db;
    //if user provided, then give as is
    //if hash table exists in server, then build it depending on sequencing technology
    $subject_dbs = $user_db;
    if (empty($subject_dbs)) {
      switch ($values['rtype']) {
        case 'Sanger ABI':
          $subject_dbs = ' -save '. $local_db .'_sanger';
          break;

        case 'Illumina-Solexa':
          $subject_dbs = ' -save '. $local_db .'_solexa';
          break;

        case '454':
          $subject_dbs = ' -save '. $local_db .'_454';
          break;
      }
    }
    //dpm(array($local_db,$user_db,$subject_dbs));
    //$more_data = array('database' => $subject_dbs);
  }

  gmod_dbsf_parameter_daemon($uid, $dirpath, $par, $software, array('database' => $db_used), array('in' => " $subject_dbs ", 'out' => '>'));
}
$message = t('Your request has been submitted to the queue for processing...');
}
else {
$error_operation = reset($operations);
$uids            = implode(' ', $uid_array);
$message         = t('An error occurred while processing your request. Please '. l('contact', 'contact') ." the administrator reporting: $uids."
);
//dpm($error_operation);
drupal_set_message($message, 'warning');
return FALSE;
}
drupal_set_message($message, 'warning');
}

/**
 * Parse a ssaha result of a software run
 *
 * @param $file
 * full path to file
 *
 * @return
 * array with data
 */
function biosoftware_bench_ssaha2_parse_result($outfile) {
/*PLUGIN_EDIT_HERE
   * You will need to edit this entire function
   * to work with the output of your software.
   * This function works with GFF files
   * See the biosoftware_bench.module for an XML solution used for BLAST
   */



/*
   * In BLAST the array looks like this:
   *       $result_array[$query_id][$hit_rank] = array(
   'name_of_hit' => $name_of_hit,
   'accession_of_hit' => $accession_of_hit,
   'average_bit_score' => $average_bit_score,
   'top_bit_score' => $top_bit_score,
   'average_evalue' => $average_evalue,
   'lowest_evalue' => $lowest_evalue,
   'average_identity' => $average_identity,
   'average_similarity' => $average_similarity,
   );

   So we suggest you keep it similar
   */



//SSAHA uses the format name in the beginning of every output line
//SSAHA uses GFF version 2 which is tab delimited
// 0.query name 1. SSAHA2 2. similarity 3. query start 4. query end
// 5. score Smith-Waterman score obtained from the cross_match alignment.
// 6. strand of query 7. no phase
// 8. Subject "subject name" start end
// 8 Align <local_query_start> <local_subject_start> [length];
/* SSAHA FORMAT
  $VAR0  = 'ALIGNMENT::50';    tag+score cutoff
  $VAR1  = 131         score
  $VAR2  = 'query';   name of query
  $VAR3  = '2L';      name of hit
  $VAR4  = '1';       start of query
  $VAR5  = '1279';    end of query
  $VAR6  = '64311';   start of hit
  $VAR7  = '65589';   end of hit
  $VAR8  = 'F';       strand
  $VAR9  = '1279';    bases matching
  $VAR10 = '100.00';  % identity
  $VAR11 = '1279';   length of query
  */


if (!file_exists($outfile)) {
drupal_set_message(t('Failed to find %outfile.', array('%outfile' => $outfile)), 'error');
return FALSE;
}
// return an associative array
$result_array = array();
// how to know the result line
$format_grep = 'ALIGNMENT::';


// OS compatibility
ini_set('auto_detect_line_endings', TRUE);
if ($inhandle = fopen($outfile, 'rb')) {
$rank = 0;
$previous_query = '';
while (!feof($inhandle)) {
  $line = trim(fgets($inhandle));
  if (!empty($line) && strpos($line, $format_grep) !== FALSE) {
    $line = preg_replace('/\s+/', ' ', $line);
    $data = explode(" ", $line);
    if (!empty($data[11])) {
      $result_array[$data[2]][$rank] = array(
        'accession_of_hit' => $data[3],
        'length_of_alignment' => $data[11],
        'top_bit_score' => $data[1],
        'orientation' => $data[8],
        'average_identity' => $data[10],
      );
      if ($data[2] == $previous_query) {
        $rank++;
      }
      else {
        $rank = 0;
      }
      $previous_query = $data[2];
    }
  }
}
fclose($inhandle);
}
if (empty($result_array)) {
// No data or error
return FALSE;
}
else {
return $result_array;
}
}

function biosoftware_bench_ssaha2_result_page($submission_uid = NULL, $level = 1) {
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
dpm("Did not find $infile");
return 'Submission ID not found. Maybe it has expired or is invalid?';
}
if (empty($data)) {
dpm("$infile has no data");
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


$content = t('Your report is not ready yet. Please wait or ') . l('try again', 'bench/ssaha2_result/'. $submission_uid
) . t(' in a few moments or bookmark this page.');
if (file_exists($graph)) {
  $results = drupal_get_form('biosoftware_bench_ssaha2_result_form', $uid, $outfile, $level, $submission_uid);
  if (strlen($results) < 600) {
    $results = "No (more) hits found for any of your queries using $algorithm. Perhaps you can try a different algorithm or database?";
  }
  $content = '<p>'. l(
    /*PLUGIN_EDIT_HERE
       * link filename is $outfile.txt here but for blast it is $outfile.html, set it to what your program does.
       */
    theme_image($graph, 'SSAHA2 results of first ten queries', 'SSAHA2 results graphical overview', array(), FALSE)
    , base_path() . $outfile .'.txt', array('html' => TRUE, 'external' => TRUE, 'attributes' => array('target' => '_blank'))
  ) .'</p>'. $results;
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
  drupal_set_html_head('<meta http-equiv="refresh" content="15;url='. $base_url .'/bench/ssaha2_result/'. $submission_uid .'" />');
}
}
return $content;
}

/**
 * Implements hook_form()
 *
 * Used to produce a table of results with checkboxes which we can then
 * select for processing (e.g. downloading). This form is themed.
 * @see biosoftware_bench_software_result_form()
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
function biosoftware_bench_ssaha2_result_form($form_state, $uid, $outfile, $level = 1, $suid = NULL) {
$form = array();
if (empty($outfile)) {
return $form;
}
if (empty($level)) {
$level = 1;
}
$array_from_result = array();
$array_from_result = biosoftware_bench_ssaha2_parse_result($outfile);
// all data
$table_row = array();
// for checkboxes to get hit
$row_check_data = array();

if (empty($array_from_result) || is_array($array_from_result) === FALSE) {
return $form;
}
/* PLUGIN_EDIT_HERE
   * we are going to customize the result form
   * from the function biosoftware_bench_ssaha2_parse_result we get
   *        'accession_of_hit'
   'length_of_alignment'
   'top_bit_score'
   'orientation'
   'average_identity'
   */


$weight_series = 0;
$query_rank    = 0;
$total_queries = count($array_from_result);
foreach ($array_from_result as $query_id => $hit_array) {
$query_rank++;
if ($query_rank != $level) {
  continue;
}
foreach ($hit_array as $hit_rank => $data) {
  $hit_id          = $data['accession_of_hit'];
  $length_of_align = $data['length_of_alignment'];
  $bit_score       = $data['top_bit_score'];
  $orientation     = $data['orientation'];
  $identity        = $data['average_identity'];
  $links           = ''. l(
    theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/ssaha2.png', 'Report as SSAHA2', 'Report as SSAHA2 (recommended for parsing with other software)'), base_path() . $outfile,
    array(
      'html' => TRUE,
      'external' => TRUE,
      'attributes' => array('target' => '_blank'),
    )
  ) . l(
    theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/txt.png', 'Report as text', 'Report as text'), base_path() ."$outfile.Query$query_rank.txt",
    array(
      'html' => TRUE,
      'external' => TRUE,
      'attributes' => array('target' => '_blank'),
    )
  );
  $table_row[] = array($query_id, $hit_id, $bit_score, $identity, "Strand ". $orientation ."; Aln length ". $length_of_align ." bp", $links);
  // checkboxes will retrieve hit_id but will be anchored on query_id
  // key to be passed for checkboxes and value to show to user.
  #   $row_check_data[$hit_id] = $query_id;
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
$base              = $base_url .'/bench/ssaha2_result/'. $uid;
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
$table_header     = array('Query ID', 'Hit ID', 'SW Score', '% identity ', 'Notes', 'Links');
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
$form['#action'] = url("bench/get_sequences/$uid");
return $form;
}

