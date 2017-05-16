<?php
// $Id$

/**
 * @file
 * Universal Search Server Module plugin
 * A Drupal module to implement annot8r  usage via the Universal Search Server
 * Developed by Alexie Papanicolaou. Licensed under GPLv3
 * University of Exeter. Licensed under GPLv3
 * @see http://insectacentral.org/
 *
 * This module requires annot8r to be properly installed
 * annot8r: GO, EC and KEGG annotation of EST datasets
 * by Ralf Schmid and Mark L Blaxter
 * BMC Bioinformatics 2008, 9:180doi:10.1186/1471-2105-9-180
 * http://www.biomedcentral.com/1471-2105/9/180
 *
 *
 *
 */

/*
 * We assume the annot8r has been already setup and is working well.
 *
 * For server usage we will implement
 *  annot8r (wrapper for all software)
 *
 * The following functions are needed in this include file for annot8r to operate
 *  biosoftware_bench_admin_annot8r_page
 *  biosoftware_bench_annot8r_page
 *  biosoftware_bench_add_software_annot8r
 *
 * The following can be done conditionally via the main module
 *  administration-> where the executable is
 *
 * The following must be done conditionally via the main .install file
 *  add annot8r as a software in ontology
 *  add options as CV
 *
 **/
function biosoftware_bench_annot8r_menu() {
  $items = array();
  $items['admin/bench/annot8r'] = array(
    'file' => 'includes/biosoftware_bench_annot8r.inc',
    'title' => 'annot8r Server administration',
    'page callback' => 'biosoftware_bench_admin_annot8r_page',
    'access arguments' => array('Administrate annot8r Server'),
    'description' => 'Configure annot8r specific settings',
    'type' => MENU_CALLBACK,
  );
  $items['bench/annot8r'] = array(
    'file' => 'includes/biosoftware_bench_annot8r.inc',
    'title' => 'biosoftware_bench annot8r Server',
    'page callback' => 'biosoftware_bench_annot8r_page',
    'access arguments' => array('Access annot8r Server'),
    'type' => MENU_CALLBACK,
  );
  $items['bench/annot8r_result'] = array(
    'file' => 'includes/biosoftware_bench_annot8r.inc',
    'title' => 'biosoftware_bench annot8r Results',
    'page callback' => 'biosoftware_bench_annot8r_result_page',
    'page arguments' => array(2, 3),
    'access arguments' => array('Access annot8r Server'),
    'type' => MENU_CALLBACK,
  );
  return $items;
}

function biosoftware_bench_software_annot8r_core_settings_form($form_state) {
  $software     = 'annot8r';
  $form         = array();
  $check_active = biosoftware_bench_check_software_active('annot8r');
  $external_server = gmod_dbsf_get_add_var('biosoftware_bench_external_server');
  if (empty($check_active)) {
    return $form;
  }
  $plugin_path = drupal_get_path('module', 'biosoftware_bench') .'/includes';
  if (file_exists($plugin_path .'/biosoftware_bench_blastall.inc')) {
    require_once($plugin_path .'/biosoftware_bench_blastall.inc');
  }
  if (biosoftware_bench_check_blast() !== TRUE) {
    drupal_set_message(t("Will not show settings for %software: the BLASTALL dependency does not seem to be fully setup yet.<br>", array('%software' => $software)) .
      l('See the core software settings page', 'admin/bench/software') .' and '. l('BLASTALL settings page', 'admin/bench/blastall'), 'error'
    );
    return $form;
  }

  $select_software_setting = 'SELECT value FROM {gmod_dbsf_softwareprop} where software_id=' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ' AND rank=0 AND type_id=' . "(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')";
  $core_exec = db_fetch_array(db_query($select_software_setting, $software, 'executable'));
  $physprop_exec = db_fetch_array(db_query($select_software_setting, $software, 'executable physprop'));
  $annot8r_executable = db_fetch_array(db_query($select_software_setting, $software, 'executable annot8r'));
  $annot8r2gff_exec = db_fetch_array(db_query($select_software_setting, $software, 'executable annot8r2gff'));
  $uniprot_GO = db_fetch_array(db_query($select_software_setting, $software, 'GOdb'));
  $uniprot_EC = db_fetch_array(db_query($select_software_setting, $software, 'ECdb'));
  $uniprot_KEGG = db_fetch_array(db_query($select_software_setting, $software, 'KEGGdb'));
  /*
 * $annot8r_exec= db_fetch_array(db_query($select_software_setting, $software, 'executable annot8r'));
  $physprop_exec= db_fetch_array(db_query($select_software_setting, $software, 'executable physprop'));
  $annot8r2gff_exec= db_fetch_array(db_query($select_software_setting, $software, 'executable annot8r2gff'));
 */


  $dbs = array(0 => 'None');
  //$dataset_array = gmod_dbsf_get_add_resource('dataset', NULL, FALSE, 'BOTH');
  $dataset_array = biosoftware_bench_get_user_datasets();
  if (!empty($dataset_array['protein'])) {
    $dbs = $dbs + $dataset_array['protein'];
  }
  $dbcredentials    = array();
  $res              = db_fetch_array(db_query($select_software_setting, $software, 'credential user'));
  $dbcredentials[1] = $res['value'];
  $res              = db_fetch_array(db_query($select_software_setting, $software, 'credential host'));
  $dbcredentials[2] = $res['value'];
  $res              = db_fetch_array(db_query($select_software_setting, $software, 'credential port'));
  $dbcredentials[3] = $res['value'];
  $res              = db_fetch_array(db_query($select_software_setting, $software, 'credential pass'));
  $dbcredentials[4] = $res['value'];
  $cpus             = gmod_dbsf_get_add_var('biosoftware_bench_cpus');
  if (empty($dbcredentials[1])) {
    global $db_url;
    $db = !empty($db_url['default']) ? $db_url['default'] : $db_url;
    preg_match('/:\/\/(\S+):\S+@(\S+):(\d+)\//', $db, $dbcredentials);
  }
  $form['settings'] = array(
    '#type' => 'fieldset',
    '#title' => $software . t(' core settings'),
    '#description' => t('Please set some important settings for %software', array('%software' => $software)),
    'submit' => array(
      '#type' => 'submit',
      '#value' => "Save core $software settings",
      '#weight' => 5,
    ), 'annot8r_wrapper_executable' => array(
      '#type' => 'textfield',
      '#title' => t('Core executable path'),
      '#description' => t('Please provide the <strong>full path</strong> to the annot8r_wrapper.pl executable on the server, e.g. /usr/bin/annot8r_wrapper.pl .'),
      '#required' => TRUE,
      '#default_value' => $core_exec['value'],
    ), 'annot8r_executable' => array(
      '#type' => 'textfield',
      '#title' => t('Annot8r (est2assembly) path'),
      '#description' => t('Please provide the <strong>full path</strong> pointing to the annot8r script from est2assembly, e.g. /usr/bin/annot8rAP-gff.pl .'),
      '#required' => TRUE,
      '#default_value' => $annot8r_executable['value'],
    ), 'physprop_executable' => array(
      '#type' => 'textfield',
      '#title' => t('Annot8r physical properties executable path'),
      '#description' => t('Please give the full path to the annot8r_physprop.pl (or similar) script.'),
      '#required' => TRUE,
      '#default_value' => $physprop_exec['value'],
    ), 'annot8r2gff_executable' => array(
      '#type' => 'textfield',
      '#title' => t('Annot8r to GFF conversion script path'),
      '#description' => t('Please provide the full path to the ic_annot8r2gff.pl script.'),
      '#required' => TRUE,
      '#default_value' => $annot8r2gff_exec['value'],
    ), 'cpus' => array(
      '#type' => 'textfield',
      '#title' => t('Number of CPUs'),
      '#description' => t('Please provide how many CPUs you would like to use for BLAST searches.'),
      '#required' => TRUE,
      '#default_value' => !empty($cpus) ? $cpus : '2',
    ),'external_server'=> array(
      '#type'=>'checkbox',
      '#title'=> t('I\' using an external server'),
      '#description'=>t('Check this if you do not want Drupal to check if the file paths actually exist on the web-server (because you are running it on an external machine'),
      '#default_value'=>$external_server,
    ),
    'blastdbs' => array(
      '#type' => 'fieldset',
      '#title' => 'BLAST databases',
      '#description' => t('These are set throught the main biosoftware_bench '. l('dataset settings', 'admin/bench/dataset') .' page.'),
      'a8r_gobase' => array(
        '#type' => 'select',
        '#title' => t('Gene Ontology Uniprot BLAST Database'),
        '#description' => t('Please select the Uniprot derived BLAST database with Gene Ontology (GO) terms, as formatted with annot8r or distributed with est2assembly.'),
        '#required' => TRUE,
        '#options' => $dbs,
        '#default_value' => $uniprot_GO['value'],
      ), 'a8r_ecbase' => array(
        '#type' => 'select',
        '#title' => t('Enzyme Classification Uniprot BLAST Database'),
        '#description' => t('Please select the Uniprot derived BLAST database with Enzyme Classification (EC) terms, as formatted with annot8r or distributed est2assembly.'),
        '#required' => TRUE,
        '#options' => $dbs,
        '#default_value' => $uniprot_EC['value'],
      ), 'a8r_keggbase' => array(
        '#type' => 'select',
        '#title' => t('Kyoto Encyclopedia of Genes and Genomes Uniprot BLAST Database'),
        '#description' => t('Please select the Uniprot derived BLAST database with KEGG terms, as formatted with annot8r or distributed with est2assembly.'),
        '#required' => TRUE,
        '#options' => $dbs,
        '#default_value' => $uniprot_KEGG['value'],
      ),
    ), 'pgdata' => array(
      '#weight' => 2,
      '#type' => 'fieldset',
      '#tree' => TRUE,
      '#title' => 'Database server settings',
      'pgport' => array(
        '#weight' => 0,
        '#type' => 'textfield',
        '#title' => t('Postgres database port'),
        '#description' => t("Please provide the port number of the database server which contains the postgres databases created by annot8r (a8r_gobase, a8r_ecbase, a8r_keggbase)."),
        '#required' => TRUE,
        '#default_value' => $dbcredentials[3],
        '#size' => 8,
      ), 'pghost' => array(
        '#weight' => 1,
        '#type' => 'textfield',
        '#title' => t('Postgres database host'),
        '#description' => t('Please provide the hostname of the above database server.'),
        '#required' => TRUE,
        '#default_value' => $dbcredentials[2],
        '#size' => 20,
      ), 'pguser' => array(
        '#weight' => 2,
        '#type' => 'textfield',
        '#title' => t('Postgres database username'),
        '#description' => t('Please provide the username to conenct to the above database server. Note that read ("select") permissions must have been set for this user.'),
        '#required' => TRUE,
        '#default_value' => $dbcredentials[1],
        '#size' => 20,
      ), 'pgpass' => array(
        '#weight' => 3,
        '#type' => 'textfield',
        '#title' => t('Postgres username password'),
        '#description' => t('Please provide the database server password of the above user. Note that some database server configurations may not require the "correct" password (if password method has been set to TRUST in the pg_hba file).'),
        '#size' => 20,
        // never taken from db_url
        '#default_value' => $dbcredentials[4],
      ),
    ),
  );
  return $form;
}

function biosoftware_bench_software_annot8r_core_settings_form_validate($form, &$form_state) {
  $software         = 'annot8r';
  $core_exec        = escapeshellcmd(trim($form_state['values']['annot8r_wrapper_executable']));
  $physprop_exec    = escapeshellcmd(trim($form_state['values']['physprop_executable']));
  $annot8r2gff_exec = escapeshellcmd(trim($form_state['values']['annot8r2gff_executable']));
  $annot8r_exec     = escapeshellcmd(trim($form_state['values']['annot8r_executable']));
  $pgport           = check_plain(trim($form_state['values']['pgdata']['pgport']));
  $pghost           = check_plain(trim($form_state['values']['pgdata']['pghost']));
  $pguser           = check_plain(trim($form_state['values']['pgdata']['pguser']));
  $pgpass           = $form_state['values']['pgdata']['pgpass'];
  $dbnames          = array('a8r_ecbase', 'a8r_gobase', 'a8r_keggbase');
  $cpus             = check_plain(trim($form_state['values']['cpus']));
  if (!empty($form_state['values']['external_server'])){
    gmod_dbsf_get_add_var('biosoftware_bench_external_server','TRUE');
  }elseif (gmod_dbsf_get_add_var('biosoftware_bench_external_server')){
    gmod_dbsf_get_add_var('biosoftware_bench_external_server','DELETE');
  }
  if (!empty($cpus)) {
    if (!is_numeric($cpus)) {
      form_set_error('cpus', t('Number of CPUs must be a number.'));
    }
  }
  if (!empty($core_exec)) {
    if (!file_exists($core_exec) && !empty($form_state['values']['external_server'])) {
      form_set_error('annot8r_wrapper_executable', t('Cannot find/access executable %core_exec on the server. Perhaps a permissions error or a symbolic link was used?', array('%core_exec' => $core_exec)));
    }
    elseif (!is_executable($core_exec) && !empty($form_state['values']['external_server'])) {
      form_set_error('annot8r_wrapper_executable', t('Cannot execute the %f program.', array('%f' => $core_exec)));
    }
    $form_state['values']['annot8r_wrapper_executable'] = $core_exec;
  }
  if (!empty($physprop_exec)) {
    if (!file_exists($physprop_exec) && !empty($form_state['values']['external_server'])) {
      form_set_error('physprop_executable', t('Cannot find/access executable %physprop_exec on the server. Perhaps a permissions error or a symbolic link was used?', array('%physprop_exec' => $physprop_exec)));
    }
    elseif (!is_executable($physprop_exec) && !empty($form_state['values']['external_server'])) {
      form_set_error('physprop_executable', t('Cannot execute the %f program.', array('%f' => $physprop_exec)));
    }
  }
  if (!empty($annot8r_exec)) {
    if (!file_exists($annot8r_exec) && !empty($form_state['values']['external_server'])) {
      form_set_error('annot8r_executable', t('Cannot find/access executable %annot8r_exec on the server. Perhaps a permissions error or a symbolic link was used?', array('%annot8r_exec' => $annot8r_exec)));
    }
    elseif (!is_executable($annot8r_exec) && !empty($form_state['values']['external_server'])) {
      form_set_error('annot8r_executable', t('Cannot execute the %f program.', array('%f' => $annot8r_exec)));
    }
  }
  if (!empty($annot8r2gff_exec)) {
    if (!file_exists($annot8r2gff_exec) && !empty($form_state['values']['external_server'])) {
      form_set_error('annot8r2gff_executable', t('Cannot find/access executable %annot8r2gff_exec on the server. Perhaps a permissions error or a symbolic link was used?', array('%annot8r2gff_exec' => $annot8r2gff_exec)));
    }
    elseif (!is_executable($annot8r2gff_exec) && !empty($form_state['values']['external_server'])) {
      form_set_error('annot8r2gff_executable', t('Cannot execute the %f program.', array('%f' => $annot8r2gff_exec)));
    }
  }
  if (!empty($pgport) && !empty($pghost) && !empty($pguser)) {
    if (!is_numeric($pgport)) {
      form_set_error('pgport', t('Port must be a number!'));
      return;
    }
    foreach ($dbnames as $dbname) {
      if (empty($form_state['values'][$dbname])) {
        form_set_error("$dbname", t('You need to speficy all of the UniProt databases.'));
        return;
      }
      $test_connection = pg_connect("dbname=$dbname host=$pghost port=$pgport user=$pguser password='$pgpass'");
      if ($test_connection === FALSE) {
        form_set_error("pgdata", t('Database %dbname either does not exists or cannot connect to it with the credentials you provided.', array('%dbname' => $dbname)));
        return;
      }
      else {
        pg_close($test_connection);
      }
    }
  }
  $form_state['values']['blastall_executable'] = $core_exec;
  $form_state['values']['physprop_executable'] = $physprop_exec;
  $form_state['values']['annot8r2gff_executable'] = $annot8r2gff_exec;
  $form_state['values']['annot8r_executable'] = $annot8r_exec;
  $form_state['values']['pgdata']['pgport'] = $pgport;
  $form_state['values']['pgdata']['pghost'] = $pghost;
  $form_state['values']['pgdata']['pguser'] = $pguser;
  $form_state['values']['cpus'] = $cpus;
}

function biosoftware_bench_software_annot8r_core_settings_form_submit($form, &$form_state) {
  $software = 'annot8r';
  $store_software_setting = 'INSERT INTO {gmod_dbsf_softwareprop} (software_id,type_id,rank,value) VALUES (' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ",(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')" . ",0,'%s')";
  $delete_software_setting = 'DELETE FROM {gmod_dbsf_softwareprop} where software_id=' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ' AND rank=0 AND type_id=' . "(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')";
  $select_software_setting = 'SELECT value FROM {gmod_dbsf_softwareprop} where software_id=' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ' AND rank=0 AND type_id=' . "(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id where cv.name='software_setting' AND cvterm.name='%s')";

  $core_exec = db_fetch_array(db_query($select_software_setting, $software, 'executable'));
  $user_core_exec = $form_state['values']['annot8r_wrapper_executable'];
  if (!empty($user_core_exec) && $user_core_exec != $core_exec['value']) {
    db_query($delete_software_setting, $software, 'executable');
    db_query($store_software_setting, $software, 'executable', $user_core_exec);
    drupal_set_message(t('%user_core_exec is now the %software executable.', array('%user_core_exec' => $user_core_exec, '%software' => $software)), 'warning');
  }
  gmod_dbsf_get_add_var('biosoftware_bench_cpus', $form_state['values']['cpus']);

  $user_annot8r_exec = $form_state['values']['annot8r_executable'];
  $user_physprop_exec = $form_state['values']['physprop_executable'];
  $user_annot8r2gff_exec = $form_state['values']['annot8r2gff_executable'];
  $user_pgport = $form_state['values']['pgdata']['pgport'];
  $user_pghost = $form_state['values']['pgdata']['pghost'];
  $user_pguser = $form_state['values']['pgdata']['pguser'];
  $user_pgpass = $form_state['values']['pgdata']['pgpass'];
  $user_uniprot_GO = $form_state['values']['a8r_gobase'];
  $user_uniprot_EC = $form_state['values']['a8r_ecbase'];
  $user_uniprot_KEGG = $form_state['values']['a8r_keggbase'];


  $annot8r_exec     = db_fetch_array(db_query($select_software_setting, $software, 'executable annot8r'));
  $physprop_exec    = db_fetch_array(db_query($select_software_setting, $software, 'executable physprop'));
  $annot8r2gff_exec = db_fetch_array(db_query($select_software_setting, $software, 'executable annot8r2gff'));
  $pgport           = db_fetch_array(db_query($select_software_setting, $software, 'credential port'));
  $pghost           = db_fetch_array(db_query($select_software_setting, $software, 'credential host'));
  $pguser           = db_fetch_array(db_query($select_software_setting, $software, 'credential user'));
  $pgpass           = db_fetch_array(db_query($select_software_setting, $software, 'credential pass'));
  $uniprot_GO       = db_fetch_array(db_query($select_software_setting, $software, 'GOdb'));
  $uniprot_EC       = db_fetch_array(db_query($select_software_setting, $software, 'ECdb'));
  $uniprot_KEGG     = db_fetch_array(db_query($select_software_setting, $software, 'KEGGdb'));

  if (!empty($user_physprop_exec) && $user_physprop_exec != $physprop_exec['value']) {
    db_query($delete_software_setting, $software, 'executable physprop');
    db_query($store_software_setting, $software, 'executable physprop', $user_physprop_exec);
    drupal_set_message(t('%user_physprop_exec is now the physprop executable.', array('%user_physprop_exec' => $user_physprop_exec)), 'warning');
  }
  if (!empty($user_annot8r_exec) && $user_annot8r_exec != $annot8r_exec['value']) {
    db_query($delete_software_setting, $software, 'executable annot8r');
    db_query($store_software_setting, $software, 'executable annot8r', $user_annot8r_exec);
    drupal_set_message(t('%user_annot8r_exec is now the annot8r executable.', array('%user_annot8r_exec' => $user_annot8r_exec)), 'warning');
  }
  if (!empty($user_annot8r2gff_exec) && $user_annot8r2gff_exec != $annot8r2gff_exec['value']) {
    db_query($delete_software_setting, $software, 'executable annot8r2gff');
    db_query($store_software_setting, $software, 'executable annot8r2gff', $user_annot8r2gff_exec);
    drupal_set_message(t('%user_annot8r2gff_exec is now the annot8r2gff executable.', array('%user_annot8r2gff_exec' => $user_annot8r2gff_exec)), 'warning');
  }
  if (!empty($user_pgport) && $user_pgport !== $pgport['value']) {
    db_query($delete_software_setting, $software, 'credential port');
    db_query($store_software_setting, $software, 'credential port', $user_pgport);
    drupal_set_message(t('%user_pgport is now the PostGres server port.', array('%user_pgport' => $user_pgport)), 'warning');
  }
  if (!empty($user_pghost) && $user_pghost !== $pghost['value']) {
    db_query($delete_software_setting, $software, 'credential host');
    db_query($store_software_setting, $software, 'credential host', $user_pghost);
    drupal_set_message(t('%user_pghost is now the PostGres hostname.', array('%user_pghost' => $user_pghost)), 'warning');
  }
  if (!empty($user_pguser) && $user_pguser !== $pguser['value']) {
    db_query($delete_software_setting, $software, 'credential user');
    db_query($store_software_setting, $software, 'credential user', $user_pguser);
    drupal_set_message(t('%user_pguser is now the PostGres username.', array('%user_pguser' => $user_pguser)), 'warning');
  }
  if (!empty($user_pgpass) && $user_pgpass !== $pgpass['value']) {
    db_query($delete_software_setting, $software, 'credential pass');
    db_query($store_software_setting, $software, 'credential pass', $user_pgpass);
    drupal_set_message(t('%user_pgpass is now the PostGres username password.', array('%user_pgpass' => $user_pgpass)), 'warning');
  }

  if (!empty($user_uniprot_GO) && $user_uniprot_GO != $uniprot_GO['value']) {
    db_query($delete_software_setting, $software, 'GOdb');
    db_query($store_software_setting, $software, 'GOdb', $user_uniprot_GO);
    drupal_set_message(t('Set the UniProt GO dataset.'), 'warning');
  }
  if (!empty($user_uniprot_EC) && $user_uniprot_EC != $uniprot_EC['value']) {
    db_query($delete_software_setting, $software, 'ECdb');
    db_query($store_software_setting, $software, 'ECdb', $user_uniprot_EC);
    drupal_set_message(t('Set the UniProt EC dataset.'), 'warning');
  }
  if (!empty($user_uniprot_KEGG) && $user_uniprot_KEGG != $uniprot_KEGG['value']) {
    db_query($delete_software_setting, $software, 'KEGGdb');
    db_query($store_software_setting, $software, 'KEGGdb', $user_uniprot_KEGG);
    drupal_set_message(t('Set the UniProt KEGG dataset.'), 'warning');
  }
}

function biosoftware_bench_admin_annot8r_page($software = 'annot8r') {
  // Include biosoftware_bench admin function file
  require_once(drupal_get_path('module', 'biosoftware_bench') .'/includes/biosoftware_bench_admin.inc');
  $check_active = biosoftware_bench_check_software_active($software);
  if (empty($check_active)) {
    drupal_set_message(t("I'm sorry, but %software does not seem to have been activated yet.<br>", array('%software' => $software)) .
      l(t('See the software settings page'), 'admin/bench/software'), 'error'
    );
    return FALSE;
  }
  $path_form = drupal_get_form('biosoftware_bench_software_annot8r_core_settings_form');

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

function biosoftware_bench_check_annot8r($software = 'annot8r') {
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
          if (strpos($term_name, 'executable') !== FALSE || strpos($term_name, 'credential') !== FALSE || strpos($term_name, 'db') !== FALSE) {
            $check++;
          }
        }
      }
    }
  }
  if ($check < 10) {
    return 'variables';
  }

  return TRUE;
}

function biosoftware_bench_annot8r_page() {
  $setup = biosoftware_bench_check_annot8r('annot8r');
  if ($setup === TRUE) {
    return drupal_get_form('biosoftware_bench_annot8r_form');
  }
  else {
    $return_text = "<p>I'm sorry but your annot8r software has not been setup properly.</p><ul>";
    switch ($setup) {
      case 'active':
        $return_text .= '<li>annot8r has not been activated. Please activate it at the '. l('software settings page', 'admin/bench/software');
        break;

      case 'variables':
        $return_text .= '<li>annot8r variable paths have not been defined. Please define them at at the '. l('annot8r settings page', 'admin/bench/software', array('fragment' => 'annot8r'));
        break;
    }
    drupal_set_message(t('!return_text', array('!return_text' => $return_text .'</ul>')), 'error');
    return ' ';
  }
}

function biosoftware_bench_annot8r_form($form_state) {
  $form             = array();
  $software         = 'annot8r';
  $software_options = biosoftware_bench_generate_options($software);

  $option_form = array();
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
  }

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
      '#title' => t('Enter query <strong>protein</strong> sequence in simple text, FASTA format'),
    ),
    'query_file_'. $software => array(
      '#type' => 'file',
      '#title' => t('or upload <strong>protein</strong> query sequences in FASTA format'),
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
      '#tree' => TRUE,
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

function biosoftware_bench_annot8r_form_validate($form, &$form_state) {
  $data = $form_state['values'];

  if (empty($data['query_sequence']) && empty($_FILES['files']['tmp_name']['query_file_' . 'annot8r'])) {
    form_set_error("query_sequence", t('It seems you neither gave a query sequence nor uploaded one.'));
    return FALSE;
  }
  if (!empty($data['query_sequence'])) {
    if (!empty($_FILES['files']['tmp_name']['query_file_' . 'annot8r'])) {
      form_set_error("query_sequence", t('It seems you both typed a query sequence and uploaded one.'));
      return FALSE;
    }
    // First Query (i.e. check for query first; query must exist if subject exists)
    $result = gmod_dbsf_validate_seq_protein($data['query_sequence']);
    if ($result === FALSE || is_numeric($result)) {
      form_set_error("query_sequence", t('Sorry your sequence does not seem to be a valid protein sequence (%result errors).', array('%result' => $result)));
      return FALSE;
    }
    else {
      $form_state['query_sequence'] = $result;
    }
  }
}

function biosoftware_bench_annot8r_form_submit($form, &$form_state) {
  $software          = 'annot8r';
  $file_size_allowed = gmod_dbsf_get_add_var('biosoftware_bench_upload_size');
  $cpus              = gmod_dbsf_get_add_var('biosoftware_bench_cpus');
  $validators_file   = array('file_validate_size' => array($file_size_allowed * 1000));
  $tmpdir            = file_directory_temp();
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
  $form_uid        = gmod_dbsf_create_uid($sessionid, $timestamp, $software);
  $uid             = gmod_dbsf_create_uid($sessionid, $timestamp, 'annot8r');
  $batch_file_data = array();

  $file_type = 'protein';
  $verify    = array();
  $par       = ' -cpu '. $cpus;
  $par .= biosoftware_bench_annot8r_software_pars($data['applications']['GO'], $data['applications']['EC'], $data['applications']['KEGG']);

  if (isset($_FILES['files']) && !empty($_FILES['files']['tmp_name']['query_file_' . 'annot8r'])) {
    $file = file_save_upload('query_file_' . 'annot8r', $validators_file, $tmpdir, FILE_EXISTS_RENAME);
    if (empty($file)) {
      form_set_error('query_file_' . 'annot8r', t('Sorry your file for <em>annot8r</em> was not saved. Maybe it is too large (>%file_size_allowed Mb)? Otherwise, '. l('contact', 'contact') .' the administrator  (quote %uid).', array('%uid' => $uid, '%file_size_allowed' => $file_size_allowed)
        ));
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
    $save_array[$form_uid][$uid]['algorithm'] = 'annot8r';
  }

  // ended each algorithm
  $operations[] = array('gmod_dbsf_batch_upload_fasta', array($batch_file_data));
  $operations[] = array('gmod_dbsf_batch_save_data', array($save_array));
  $batch        = array(
    'title' => t('Preparing data needed for annot8r...'),
    'operations' => $operations,
    'init_message' => t('Starting annot8r submission...'),
    'progress_message' => t('@remaining operations remaining...'),
    'error_message' => t('Your annot8r submission encountered an error.'),
    'finished' => 'biosoftware_bench_annot8r_batch_finished',
    'file' => drupal_get_path('module', 'biosoftware_bench') .'/includes/biosoftware_bench_annot8r.inc',
  );
  batch_set($batch);
  $form_state['redirect'] = array('bench/annot8r_result', "submission_uid=$uid");
  return "bench/annot8r_result?submission_uid=$uid";
}

function biosoftware_bench_annot8r_software_pars($go = FALSE, $ec = FALSE, $kegg = FALSE) {
  $software = 'annot8r';
  $select_software_setting = 'SELECT value FROM {gmod_dbsf_softwareprop} where software_id=' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='%s')" . ' AND rank=0 AND type_id=' . "(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv on cv.cv_id=cvterm.cv_id " . " where cv.name='software_setting' AND cvterm.name='%s')";

  $blast_exec       = db_fetch_array(db_query($select_software_setting, 'blastall', 'executable'));
  $blast_exec       = $blast_exec['value'];
  $blast_data       = db_fetch_array(db_query($select_software_setting, 'blastall', 'data'));
  $blast_data       = $blast_data['value'];
  $annot8r_exec     = db_fetch_array(db_query($select_software_setting, $software, 'executable annot8r'));
  $annot8r_exec     = $annot8r_exec['value'];
  $physprop_exec    = db_fetch_array(db_query($select_software_setting, $software, 'executable physprop'));
  $physprop_exec    = $physprop_exec['value'];
  $annot8r2gff_exec = db_fetch_array(db_query($select_software_setting, $software, 'executable annot8r2gff'));
  $annot8r2gff_exec = $annot8r2gff_exec['value'];
  $uniprot_GO       = db_fetch_array(db_query($select_software_setting, $software, 'GOdb'));
  $uniprot_GO       = biosoftware_bench_get_dataset_path($uniprot_GO['value']);
  $uniprot_EC       = db_fetch_array(db_query($select_software_setting, $software, 'ECdb'));
  $uniprot_EC       = biosoftware_bench_get_dataset_path($uniprot_EC['value']);
  $uniprot_KEGG     = db_fetch_array(db_query($select_software_setting, $software, 'KEGGdb'));
  $uniprot_KEGG     = biosoftware_bench_get_dataset_path($uniprot_KEGG['value']);
  $dbcredentials    = array();
  $pguser           = db_fetch_array(db_query($select_software_setting, $software, 'credential user'));
  $pguser           = $pguser['value'];
  $pghost           = db_fetch_array(db_query($select_software_setting, $software, 'credential host'));
  $pghost           = $pghost['value'];
  $pgport           = db_fetch_array(db_query($select_software_setting, $software, 'credential port'));
  $pgport           = $pgport['value'];
  $pgpass           = db_fetch_array(db_query($select_software_setting, $software, 'credential pass'));
  $pgpass           = $pgpass['value'];

  $par = " -physprop $physprop_exec -annot $annot8r_exec -annot2gff $annot8r2gff_exec -blast $blast_exec -matrix $blast_data ". " -dbgo $uniprot_GO -dbec $uniprot_EC -dbkegg $uniprot_KEGG -host $pghost -user $pguser -port $pgport ";
  if (!empty($pgpass)) {
    $par .= " -pass $pgpass ";
  }
  if (!empty($go)) {
    $par .= ' -go';
  }
  if (!empty($ec)) {
    $par .= ' -ec';
  }
  if (!empty($kegg)) {
    $par .= ' -kegg';
  }
  // Now we just need to add the infile (-i) and perhaps organism (-o)
  return $par;
}

function biosoftware_bench_annot8r_batch_finished($success, $results, $operations) {
  $software = 'annot8r';
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
function biosoftware_bench_annot8r_result_page($submission_uid = NULL, $level = 1) {
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
    $content   = t('Your report is not ready yet. Please wait or ') . l('tbench/tware_bench/annot8r_result/'. $submission_uid
    ) . t(' in a few moments or bookmark this page.');
    if (file_exists($graph)) {
      $results = drupal_get_form('biosoftware_bench_annot8r_result_form', $uid, $outfile, $level, $submission_uid);
      if (strlen($results) < 600) {
        $results = "No (more) hits found for any of your queries using $algorithm. Perhaps you can try a different algorithm or database?";
      }
      $content = '<p>'. l(
        theme_image($graph, 'Annot8r results of first ten queries', 'annot8r results graphical overview', array(), FALSE)
        , base_path() . $outfile, array('html' => TRUE, 'external' => TRUE, 'attributes' => array('target' => '_blank'))
      ) .'</p>' . '<p>All results as :';
      if (file_exists($outfile)) {
        $content .= l(
          theme_image(drupal_get_path('module', 'biosoftware_bench') .'/images/gff.png', 'Report as GFF3', 'Report as GFF3'), base_path() . $outfile,
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
      drupal_set_html_head('<meta http-equiv="refresh" content="15;url='. $base_url .'/bench/annot8r_result/'. $submission_uid .'" />');
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
function biosoftware_bench_annot8r_result_form($form_state, $uid, $resultfile, $level = 1, $suid = NULL) {
  $form = array();
  if (empty($resultfile)) {
    return $form;
  }
  if (empty($level)) {
    $level = 1;
  }
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
      $hit_id     = '';
      $other_info = '';
      $score      = '';
      if ($gff_data['source'] == 'ANNOT8R_physprop') {
        $hit_id = 'Physical properties';
        $other_info = 'Estimations: ';
        $other_info .= 'charge at pH5: '. implode(', ', (array)$gff_data['attributes']['charge pH5']);
        $other_info .= '; at pH7: '. implode(', ', (array)$gff_data['attributes']['charge pH7']);
        $other_info .= '; at pH8: '. implode(', ', (array)$gff_data['attributes']['charge pH8']);
        $other_info .= '; charged residues: '. implode(', ', (array)$gff_data['attributes']['charged']);
        $other_info .= '; molecular weight: '. implode(', ', (array)$gff_data['attributes']['molecular_weight']);
        $other_info .= '; pI: '. implode(', ', (array)$gff_data['attributes']['pI']);
        $other_info .= implode(', ', (array)$gff_data['attributes']['Note']);
      }
      else {
        foreach ((array)$gff_data['attributes']['Note'] as $key => $note) {
          $gff_data['attributes']['Note'][$key] = l($gff_data['attributes']['Note'][$key],
            'http://amigo.geneontology.org/cgi-bin/amigo/search.cgi',
            array(
              'query' => array(
                'query' => $note,
                'term' => 'term',
              ),
              'absolute' => TRUE, 'external' => TRUE,
              'attributes' => array('target' => '_blank'),
            )
          );
        }
        $hit_id     = implode(', ', $gff_data['attributes']['Dbxref']);
        $other_info = implode(', ', (array)$gff_data['attributes']['Note']);
        $score      = $gff_data['score'];
      }
      //$table_row[] = array($query_id, $hit_id, $score, $other_info);
      $table_row[] = array($query_id, $score, $other_info);
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
  //$table_header     = array('Query ID', 'Hit ID', 'Score', 'Notes');
  $table_header     = array('Query ID', 'Score', 'Notes');
  $table_attributes = array('class' => 'padded-table sortable');
  $table_caption    = 'Overview of significant results';
  $table_array      = array(
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
    ),
  );
  $form['#action'] = url("bench/get_sequences/$uid/annot8r");
  return $form;
}

