<?php
// $Id$

/**
 * @file
 * include file for curating libraries
 */

/**
 * Curate Library
 *
 * @param $lib_id
 *
 * @return
 * HTML
 */
function genes4all_curate_page_library_curate($lib_id = NULL) {
  $return = "This page allows you to curate a library";
  if (empty($lib_id)) {
    $lib_id = $_GET['library_id'];
    if (empty($lib_id)){
      drupal_set_message(t('You did not specify an identifier.'), 'error');
      return drupal_get_form('genes4all_curate_library_view_form');
    }
  }
  $lib_id = check_plain(trim($lib_id));
  $access = gmod_dbsf_chadoapi_check_restriction('select', 'library', $lib_id);
  if ($access != TRUE) {
    return "You are not allowed to access this resource.";
  }
  $lib_table = cache_get($lib_id.'_info','cache_genes4all_explore_library');
  if (empty($lib_table)){
    $lib_table = gmod_dbsf_db_execute('chado', 'genes4all_get_lib_info', $lib_id);
    cache_set($lib_id.'_info', $lib_table, 'cache_genes4all_explore_library',time()+7*3600*24);
  }else{
    $lib_table =$lib_table->data;
  }
  $return .= '<p>'.call_user_func_array('theme_table',$lib_table) . '</p><p>';
  $return .= drupal_get_form('genes4all_curate_library_curate_form', $lib_id);
  $return .= drupal_get_form('genes4all_curate_library_link_features_form', $lib_id);
  $return .= '</p>';
  return $return;
}

function genes4all_curate_library_link_features_form($form_state, $lib_id) {
  $form = array();
  if (empty($lib_id)) {
    return $form;
  }
  $form['data'] = array(
    '#type' => 'fieldset',
    '#title' => 'Link features to this library',
    '#description'=>t('Linking features to a library allows you to search for features using the library properties, such as developmental stage or tissue.'),
    'lib_id' => array(
      '#type' => 'value',
      '#value' => $lib_id,
  ),
    'feature_file' => array(
      '#type' => 'file',
      '#title' => 'EST file',
      '#description' => t('A file containing a list of IDs (one per line). A ">" character is allowed in the beginning of the line so that you could "grep" the contents of a FASTA file. In this case only the first word after the ">" symbol is captured.'),
      '#tree' => FALSE,
  ),
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Upload data'),
      '#weight' => 10,
  ),
  );
  $form['#attributes']['enctype'] = "multipart/form-data";
  $form['#redirect'] = 'genes4all/library/curate/'. $lib_id;
  return $form;
}

function genes4all_curate_library_link_features_form_submit($form, &$form_state) {
  $lib_id = check_plain($form_state['values']['lib_id']);
  if (empty($_FILES['files']['tmp_name']['feature_file'])) {
    form_set_error('feature_file', t('Upload failed'));
    return;
  }
  // 100 Mb
  $validators_file = array('file_validate_size' => array(100000));
  $tmpdir          = file_directory_temp();
  $file            = file_save_upload('feature_file', $validators_file, $tmpdir, FILE_EXISTS_RENAME);
  if (empty($file)) {
    form_set_error('feature_file'. $algorithm, t('Sorry your file was not saved. Maybe it is too large (>100 Mb)? Otherwise, '. l('contact', 'contact') .' the administrator.'));
    return;
  }
  file_set_status($file, FILE_STATUS_TEMPORARY);
  $batch_file_data['infile'] = $file->filepath;
  $operations[] = array('gmod_dbsf_batch_save_data', array(array('lib_id' => $lib_id)));
  $operations[] = array('genes4all_curate_batch_upload_features', array($batch_file_data));
  $operations[] = array('genes4all_curate_batch_link_features', array());
  $batch = array(
    'title' => t('Preparing data for submission...'),
    'operations' => $operations,
    'init_message' => t('Starting upload...'),
    'progress_message' => t('@remaining operations remaining...'),
    'error_message' => t('Your submission encountered an error.'),
    'finished' => 'genes4all_curate_batch_upload_features_finished',
    'file' => drupal_get_path('module', 'blingblob') .'/includes/genes4all_curate_library.inc',
  );
  batch_set($batch);
  $form_state['redirect'] = array('genes4all/library/curate/'. $lib_id);
  return 'genes4all/library/curate/'. $lib_id;
}

function genes4all_curate_batch_upload_features($batch_file_data, &$context) {
  $infile = $batch_file_data['infile'];
  if (empty($infile)) {
    return;
  }
  ini_set('auto_detect_line_endings', TRUE);
  $features = array();
  if (file_exists($infile)) {
    if ($inhandle = fopen($infile, 'rb')) {
      while (!feof($inhandle)) {
        $line = trim(fgets($inhandle));
        if (!empty($line)) {
          $defline = strpos($line, '>');
          $matches = preg_match('/^>?\s*(\S+)/', $line, $match);
          if (!empty($matches)) {
            if (!empty($match[1])) {
              $feat = check_plain($match[1]);
              $features[$feat] = 1;
            }
          }
        }
      }
      fclose($inhandle);
    }
    $context['results']['features'] = $features;
  }
  else {
    drupal_set_message(t('Cannot find uploaded file. Perhaps a permissions error?'), 'error');
  }
}

function genes4all_curate_batch_link_features(&$context) {
  if (empty($context['results']['features'])) {
    drupal_set_message(t('No features found in your uploaded file.'), 'error');
    return FALSE;
  }
  $lib_id = $context['results']['data']['lib_id'];
  $context['results']['linked'] = (int)0;
  // The following are Chado SQLs so no {} fields
  $sql_link_library_feature = "INSERT INTO library_feature (library_id,feature_id) VALUES ($lib_id,%d)";
  $sql_check_feature = "SELECT feature_id as id from feature where uniquename='%s'";
  $sql_check_library_feature = "SELECT library_feature_id as id from library_feature where feature_id=%d and library_id=$lib_id";
  gmod_dbsf_db_execute('chado_edit');
  foreach ($context['results']['features'] as $feature => $dum) {
    $feat_data = db_fetch_array(db_query($sql_check_feature, $feature));
    $feature_id = $feat_data['id'];
    if (!empty($feature_id)) {
      $check = db_fetch_array(db_query($sql_check_library_feature, $feature_id));
      if (empty($check)) {
        db_query($sql_link_library_feature, $feature_id);
        $check = db_fetch_array(db_query($sql_check_library_feature, $feature_id));
      }
      if (!empty($check)) {
        $context['results']['linked']++;
      }
    }
  }
  gmod_dbsf_db_execute();
}

function genes4all_curate_batch_upload_features_finished($success, $results, $operations) {
  $message = '';
  if (!empty($success)) {
    $features = count($results['features']);
    $linked = !empty($results['linked']) ? $results['linked'] : 0;
    if (!empty($results['features'])) {
      $message .= "$features were uploaded.";
    }
    $message .= " $linked features have been linked.";
  }
  else {
    $error_operation = reset($operations);
    $message = t('An error occurred while processing your request. Please '. l('contact', 'contact') .' the administrator.');
    dpm($error_operation);
  }
  drupal_set_message($message, 'warning');
}

function genes4all_curate_library_curate_form($form_state, $lib_id = NULL) {
  $form = array();
  if (empty($lib_id)) {
    return $form;
  }
  $cv_options = array(0 => 'None') + gmod_dbsf_db_execute('chado',genes4all_curate_cv_list,'library');
  $property_types = array(0 => 'None') + gmod_dbsf_db_execute('chado',genes4all_curate_cv_list,'library_property');

  $form['data'] = array(
    '#type' => 'fieldset',
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Submit data'),
      '#weight' => 10,
  ), 'lib_id' => array(
      '#type' => 'value',
      '#value' => $lib_id,
  ),
    'ontologies' => array(
      '#type' => 'fieldset',
      '#title' => 'Ontologies',
      '#weight' => 1,
      '#tree' => TRUE,
      '#description'=>t('These selected ontologies allow you to set descriptive tags from a formal vocabulary.'),
      'cv' => array(
        '#type' => 'select',
        '#title' => 'Controlled Vocabulary (CV)',
        '#options' => $cv_options,
        '#weight' => 1,
        '#description' => t('These CVs are derived from dbEST and aim to describe a cDNA library. If you have suggestions for other CVs, please '. l('let us know', 'contact') .'.'),
  // TODO list: add an interface so that you can explore the AMIGO terms interactively)
  //      '#description'=>t('Select one CV, if you do not know what they are, please check the relevant publications/websites (e.g. '.l('AMIGO','http://amigo.geneontology.org/',array('attributes'=>array('target'=>'_blank'))))
  //.' for GO instead of potentially adding the wrong term.',
  ),
      'cvterm' => array(
        '#type' => 'textfield',
        '#title' => 'Value for this Controlled Vocabulary',
        '#description' => t('Enter here the information you want to store the particular CV. Please store only the information for that particular CV: e.g. enter "female" for Sex and not "female adult" since "adult" is the stage. Existing entries (which are derived from the original authors and thus may be erroneous) are shown during autocomplete to help you but you are not obliged to pick one of them.'),
        '#weight' => 2,
        '#size' => 30,
        '#autocomplete_path' => 'genes4all/curate/library_cvterm/autocomplete',
  ),
  ), 'dbxref' => array(
      '#type' => 'fieldset',
      '#tree' => TRUE,
      '#title' => 'Database cross-references',
      '#weight' => 8,
      '#description' => t('You can store a reference to another database (e.g. FlyBase) by providing its formal name (FlyBase) and the unique identifier from that database which you wish to store (e.g Fbgn001291). Please make sure you provide the correct spelling for the database name (the autocomplete may help you). Please note that this cross-reference can be from your own private database as well. In that case it would be beneficial if your database had an online interface or existed in other published forms (e.g. a book) so that users can explore what information this cross-reference provides.'),
      'dbname' => array(
        '#type' => 'textfield',
        '#title' => 'Database name',
        '#description' => t('The formal name of your database.'),
        '#size' => 15,
        '#autocomplete_path' => 'genes4all/curate/library_dbname/autocomplete',
  ),
      'accession' => array(
        '#type' => 'textfield',
        '#title' => 'Database ID/accession',
        '#description' => t('The unique identifier from this database.'),
        '#size' => 15,
  ),
  ),
    'properties' => array(
      '#type' => 'fieldset',
      '#title' => 'Extended properties',
      '#tree' => TRUE,
      '#weight' => 5,
      '#description' => t("It is recommended that you store extended information about the library you wish to curate. Please provide as much information as you'd like to appear in a publication."),
      'property_type' => array(
        '#title' => 'Property type',
        '#options' => $property_types,
        '#description' => t('Choose one of these terms. A protocol is the complete protocol to generate the library; a description give a short outline regarding the nature of the library and the study it may relate to. A note is a statement, usually a single sentence; you may add multiple notes.'),
        '#type' => 'select',
  ), 'property' => array(
        '#title' => 'Information to store',
        '#description' => t('The information you wish to link to the above term.'),
        '#type' => 'textarea',
  ),
  ),
  );

  $form['data']['ontologies']['evidence'] = genes4all_curation_cv_evidence('inferred', TRUE, TRUE);

  $form['#redirect'] = 'genes4all/library/curate/'. $lib_id;
  return $form;
}

function genes4all_curate_library_curate_form_validate($form, &$form_state) {
  $form_state['values']['ontologies']['cv'] = check_plain(trim($form_state['values']['ontologies']['cv']));
  $form_state['values']['ontologies']['cvterm'] = check_plain(trim($form_state['values']['ontologies']['cvterm']));
  $form_state['values']['ontologies']['evidence'] = check_plain(trim($form_state['values']['ontologies']['evidence']));

  $form_state['values']['dbxref']['dbname'] = check_plain(trim($form_state['values']['dbxref']['dbname']));
  $form_state['values']['dbxref']['accession'] = check_plain(trim($form_state['values']['dbxref']['accession']));
  $form_state['values']['lib_id'] = check_plain(trim($form_state['values']['lib_id']));
  $done=0;
  if (!empty($form_state['values']['ontologies']['cv']) || !empty($form_state['values']['ontologies']['cvterm'])) {
    if (empty($form_state['values']['ontologies']['cv']) || empty($form_state['values']['ontologies']['cvterm']) || empty($form_state['values']['ontologies']['evidence'])) {
      form_set_error('ontologies', t('You must provide all the highlighted options.'));
    }
    $done++;
  }
  if (!empty($form_state['values']['dbxref']['dbname']) || !empty($form_state['values']['dbxref']['accession'])) {
    if (empty($form_state['values']['dbxref']['dbname']) || empty($form_state['values']['dbxref']['accession'])) {
      form_set_error('dbxref', t('You must provide all the highlighted options.'));
    }
    $done++;
  }
  if (!empty($form_state['values']['properties']['property_type']) || !empty($form_state['values']['properties']['property'])) {
    if (empty($form_state['values']['properties']['property_type']) || empty($form_state['values']['properties']['property'])) {
      form_set_error('properties', t('You must provide all the highlighted options.'));
    }
    $done++;
  }

  if (empty($form_state['values']['dbxref']['dbname']) && empty($form_state['values']['dbxref']['accession'])
  && empty($form_state['values']['ontologies']['cv']) && empty($form_state['values']['ontologies']['cvterm']) && empty($form_state['values']['ontologies']['evidence'])
  && empty($form_state['values']['properties']['property_type']) && empty($form_state['values']['properties']['property'])
  ) {
    form_set_error('', t("You didn't submit anything!"));
  }
  if (empty($done)){
    form_set_error('', t("You didn't submit anything!"));
  }
}

function genes4all_curate_library_curate_form_submit($form, &$form_state) {
  $lib_id        = $form_state['values']['lib_id'];
  $cv            = $form_state['values']['ontologies']['cv'];
  $cvterm        = $form_state['values']['ontologies']['cvterm'];
  $dbname        = $form_state['values']['dbxref']['dbname'];
  $accession     = $form_state['values']['dbxref']['accession'];
  $evidence      = $form_state['values']['ontologies']['evidence'];
  $cvterm_id     = '';
  $dbxref_id     = '';
  $property_type = check_plain($form_state['values']['properties']['property_type']);
  $property      = check_plain($form_state['values']['properties']['property']);
  $done          = 0;
  gmod_dbsf_db_execute('chado_edit');
  if (!empty($cv) && !empty($cvterm)) {
    $cvterm_id = gmod_dbsf_get_add_cv_withnames($cv, $cvterm);
    $pub_id    = gmod_dbsf_get_add_pub('annotation', $evidence, NULL, 'ADD', TRUE);
    $pub_id    = $pub_id[$evidence];
  }
  if (!empty($dbname) && !empty($accession)) {
    $dbxref_id = gmod_dbsf_get_add_dbxref_withnames($dbname, $accession);
  }
  if (!empty($property_type) && !empty($property)) {
    $success = gmod_dbsf_populate_prop('library', $lib_id, $property_type, array($property));
    if (!empty($success)) {
      $done++;
    }
  }
  if (!empty($cvterm_id) && !empty($pub_id)) {
    $success = gmod_dbsf_chado_add_library_cvterm($lib_id, $cvterm_id, $pub_id);
    if (!empty($success)) {
      $done++;
    }
  }
  if (!empty($dbxref_id)) {
    $success = gmod_dbsf_chado_add_library_dbxref($lib_id, $dbxref_id);
    if (!empty($success)) {
      $done++;
    }
  }
  gmod_dbsf_db_execute();
  if (empty($done)) {
    form_set_error('', t('Could not complete request. Perhaps a permissions error?'));
  }
  else {
    //empty library cache
    cache_clear_all($lib_id.'_info','cache_genes4all_explore_library');
    drupal_set_message(t('@done request(s) stored.', array('@done' => $done)), 'warning');
  }
  $form_state['redirect'] = array('genes4all/library/curate/'. $lib_id);
  return 'genes4all/library/curate/'. $lib_id;
}

function genes4all_curate_page_library_add() {
  $return = "This page will allow you to create a new library. You then be able to curate it via ". l('the curation link', 'genes4all/library/curate');

  $return .= drupal_get_form(genes4all_curate_library_add_form);

  return $return;
}

function genes4all_curate_library_add_form($form_state) {
  // name, db and db.accession,
  $type_options    = gmod_dbsf_db_execute('chado', 'gmod_dbsf_get_add_cv_withnames', 'Library types');
  $org_array       = gmod_dbsf_db_execute('chado', 'gmod_dbsf_getspecies', NULL, TRUE);
  $species_options = $org_array['phylogeny_array'];
  if (empty($species_options)){
    $species_options =$org_array['species_array'];
    if (empty($species_options)){
      drupal_set_message(t('Sorry, this function is not available as the organism table has not been properly populated.'),'error');
      return '';
    }
  }
  unset($org_array);
  $db_options   = array(0 => 'None', 'dbEST' => 'dbEST', 'SRA' => 'Short Read Archive');
  $form         = array();
  $form['data'] = array(
    '#title' => 'New library data',
    '#type' => 'fieldset',
    '#weight' => 0,
    'name' => array(
      '#title' => 'Library name',
      '#type' => 'textfield',
      '#size' => 25,
      '#required' => TRUE,
      '#description' => t('Please provide a name which will be used for your library. It is preferred that it is both succinct and informative: e.g. the name "fdg4" is not informative and "cDNA from male and female Fall Armyworm 3rd instar larvae wing discs" is not succinct, "S. frugiperda wing disc cDNA" is probably better.'),
      '#weight' => 0,
  ), 'db_field' => array(
      '#type' => 'fieldset',
      '#title' => 'Database/repository cross-reference',
      '#description' => t('You can store a reference to another database (e.g. FlyBase) by providing its formal name (FlyBase) and the unique identifier from that database which you wish to store (e.g Fbgn001291). Please make sure you provide the correct spelling for the database name (the autocomplete may help you). Please note that this cross-reference can be from your own private database as well. In that case it would be beneficial if your database had an online interface or existed in other published forms (e.g. a book) so that users can explore what information this cross-reference provides.'),
      '#weight' => 10,
      'db_name' => array(
        '#title' => 'Database name',
        '#type' => 'select',
        '#size' => 1,
        '#description' => t('The formal name of the cross-referenced database.'),
        '#weight' => 1,
        '#options' => $db_options,
  ), 'accession' => array(
        '#title' => 'Accession',
        '#type' => 'textfield',
        '#size' => 12,
        '#description' => t('The unique accession from this database'),
        '#weight' => 2,
  ),
  ),
    'type' => array(
      '#title' => 'Library type',
      '#type' => 'select',
      '#size' => 1,
      '#weight' => 3,
      '#required' => TRUE,
      '#description' => t("What kind of library you're generating"),
      '#options' => $type_options,
  ), 'species' => array(
      '#title' => 'Species',
      '#type' => 'select',
      '#size' => 1,
      '#weight' => 1,
      '#required' => TRUE,
      '#options' => $species_options,
      '#description' => t('To which species this library belongs to. To select faster, the select box will respond to key strokes: try typing the name. If your species does not appear on the list, then you need to '. l('add', 'genes4all/organism/add') .' it first.'),
  ),
  );
  $form['button'] = array(
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Submit data'),
      '#weight' => 10,
  ),
  );
  return $form;
}

function genes4all_curate_library_add_form_validate($form, &$form_state) {
  $needed_array = array('species', 'type', 'name');
  foreach ($form_state['values'] as $key => $value) {
    if (is_array($value)) {
      continue;
    }
    $form_state['values'][$key] = check_plain(trim($value));
  }
  foreach ($needed_array as $needed) {
    if (empty($form_state['values'][$needed])) {
      form_set_error($needed, t('Some necessary information is missing...'));
    }
  }
  $db_name   = check_plain(trim($form_state['values']['db_name']));
  $accession = check_plain(trim($form_state['values']['accession']));
  $name      = check_plain(trim($form_state['values']['name']));

  $uniquename = '';
  if (!empty($db_name) && !empty($accession)) {
    $uniquename = $db_name .':'. $accession;
  }
  else {
    $db_name    = '';
    $accession  = '';
    $uniquename = $name;
  }
  $serial         = 1;
  $original_uname = $uniquename;
  $check          = gmod_dbsf_db_execute('chado', 'gmod_dbsf_decypher_library_id', $uniquename, 'uniquename', TRUE);
  while (!empty($check['library_id'])) {
    $serial++;
    $uniquename = $original_uname .' V'. $serial;
    $check = gmod_dbsf_db_execute('chado', 'gmod_dbsf_decypher_library_id', $uniquename, 'uniquename', TRUE);
  }
  if (!empty($uniquename)) {
    $form_state['values']['uniquename'] = $uniquename;
  }
  else {
    form_set_error('name', t('Cannot generate a uniquename for you submission. Please contact an administrator.'));
  }
}

function genes4all_curate_library_add_form_submit($form, &$form_state) {
  gmod_dbsf_db_execute('chado_edit');

  $lib_id = gmod_dbsf_get_add_library(array(
      'name' => $form_state['values']['name'],
      'type_id' => $form_state['values']['type'],
      'organism_id' => $form_state['values']['species'],
      'uniquename' => $form_state['values']['uniquename'],
      'db_name' => $form_state['values']['db_name'],
      'accession' => $form_state['values']['accession'],
  ));
  gmod_dbsf_db_execute();

  $form_state['redirect'] = array('genes4all/library/curate/'. $lib_id);
  return "genes4all/library/curate/$lib_id";
}

function genes4all_curate_view_library_page() {
  $ret = genes4all_curate_library_view_form();
  return drupal_get_form($ret);
}

function genes4all_curate_library_view_form($form_state) {
  $lib_data   = gmod_dbsf_db_execute('chado', 'gmod_dbsf_get_add_library');
  if (empty($lib_data)){
    drupal_set_message('There are currently no libraries stored in the database.','error');
    return FALSE;
  }
  $type_array = array();
  $name_array = array();
  $none       = array(NULL => 'None');
  foreach ($lib_data as $lib_id => $data) {
    $type_array[$data['type_id']] = $data['cvterm.name'];
    $name_array[$lib_id] = $data['library.name'];
  }
  $form = array();
  $form['data'] = array(
    '#tree' => TRUE,
    '#type' => 'fieldset',
    '#title' => 'Find a library',
    '#description' => t('Here, you can search for a library. The search will include all the limits you set below and the results will be shown on the next page.'),
    'type' => array(
      '#type' => 'select',
      '#title' => 'Type',
      '#weight' => 0,
      '#options' => $none + $type_array,
      '#description'=> t('These are the types of libraries, e.g. cDNA or genomic, that already exist in the database.'),
  ), 'name' => array(
      '#type' => 'select',
      '#title' => 'Name',
      '#weight' => 1,
      '#options' => $none + $name_array,
      '#description'=> t('Available library names.'),
  ),
  );
  $org_array = gmod_dbsf_db_execute('chado', 'gmod_dbsf_getspecies');
  $species_options = $org_array['phylogeny_array'];
  if (!empty($species_options)){
    $form['data']['species'] = gmod_dbsf_phylogeny_tree($species_options);
  }
  $form['button'] = array(
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Submit data'),
      '#weight' => 10,
  ),
  );
  return $form;
}

function genes4all_curate_library_view_form_validate($form, &$form_state) {
  // new: check that submitted species is indeed allowed.
  $org_array = gmod_dbsf_db_execute('chado', 'gmod_dbsf_getspecies', NULL, 'all');

  foreach ($form_state['values'] as $key => $value) {
    //HACK HACK HACK
    if (is_array($value)) {
      continue;
    }
    if (empty($value)) {
      unset($form_state['values'][$key]);
      continue;
    }
    if (preg_match('/^species-(\d+)/', $value, $matches)) {
      if (empty($org_array['species_array'][$matches[1]])) {
        // attempt to DOM hack: stop dead in your tracks
        form_set_error('organism', t('You requested a species which you are not allowed access to.'));
        return FALSE;
      }

      $form_state['values']['species_select'][$matches[1]] = $matches[1];
      unset($form_state['values'][$key]);
    }
  }
  if (empty($form_state['values']['data']['type']) && empty($form_state['values']['data']['name']) && empty($form_state['values']['species_select'])) {
    form_set_error('data', t('You need to select at least one limit.'));
  }
}

function genes4all_curate_library_view_form_submit($form, &$form_state) {
  $types   = $form_state['values']['data']['type'];
  $species = $form_state['values']['species_select'];
  $names   = $form_state['values']['data']['name'];
  // The following are Chado SQLs so no {} fields
  $sql = "SELECT library_id from library where library is not null ";
  if (!empty($types)) {
    $type_str = gmod_dbsf_array2csv($types);
    $sql .= " AND type_id IN ($type_str)";
  }
  if (!empty($species)) {
    $species_str = gmod_dbsf_array2csv($species);
    $sql .= " AND organism_id IN ($species_str)";
  }
  if (!empty($names)) {
    $names_str = gmod_dbsf_array2csv($names);
    $sql .= " AND names IN ($names_str)";
  }
  $sql = check_plain($sql);
  $results_array = array();
  gmod_dbsf_db_execute('chado');
  $dbdata = db_query($sql);
  while ($row = db_fetch_array($dbdata)) {
    $results_array[] = $row['library_id'];
  }
  gmod_dbsf_db_execute();
  if (empty($results_array)) {
    form_set_error('data', t('Your search returned no results.'));
    return;
  }
  $enc = implode('-', $results_array);
  $form_state['redirect'] = array('genes4all/library/result/'. $enc);
  return "genes4all/library/result/$enc";
}

function genes4all_curate_page_library_find_results($encoded_data) {
  // generate table with links which redirect you to genes4all/info/library?library_id=$uniquename
  // Name, type, species.
  if (empty($encoded_data)) {
    $encoded_data = $_POST['query'];
    drupal_set_message(t('Your search returned no results.'), 'error');
    return '';
  }
  $data             = explode('-', $encoded_data);
  $table_rows       = array();
  $table_header     = array('Task', 'Name', 'Type', 'Species');
  $table_attributes = array('class' => 'genes4all_explore-library-table sortable', 'width' => '200');
  $table_caption    = 'Overview of libraries found';
  gmod_dbsf_db_execute('chado');
  //library_id,type_id,cvterm.name as type_name,library.name as lib_name,library.organism_id
  foreach ($data as $lib_id) {
    $lib_data = gmod_dbsf_get_add_library(array('library_id' => $lib_id));
    if (empty($lib_data)) {
      continue;
    }
    $org_data = gmod_dbsf_get_add_organism(array('organism_id' => $lib_data['organism_id'], 'pgtype' => 'chado'));
    $species_name = $org_data['genus'] .' '. $org_data['species'];
    if (!empty($org_data['common_name'])) {
      $species_name .= ' ('. $org_data['common_name'] .')';
    }
    if (strlen($lib_data['lib_name']) > 55) {
      $lib_data['lib_name'] = substr($lib_data['lib_name'], 0, 50) .'...';
    }
    $table_rows[] = array(
    l('Edit', 'genes4all/library/curate', array('alias' => TRUE, 'query' => array('library_id' => $lib_id))) .'<br>'.
    l('View', 'genes4all/info/library', array('alias' => TRUE, 'query' => array('library_id' => $lib_id))),
    $lib_data['lib_name'], $lib_data['type_name'], $species_name,
    );
  }
  gmod_dbsf_db_execute();
  $table = theme_table($table_header, $table_rows, $table_attributes, $table_caption);
  return $table;
}

