<?php
// $Id$

/**
 * @file
 * include file for curating organisms
 */

/**
 * Curate an organism
 *
 * @param $organism_id
 *
 * @return
 * HTML
 */
function genes4all_curate_page_organism_curate($id = NULL) {
  if (empty($id)) {
    $id = $_POST['id'];
    if (empty($id)) {
      drupal_set_message(t('You did not specify an identifier'), 'error');
      return "This page will allow you to add or edit an organism";
    }
  }
  // genus species organism_id
  $org_ids = gmod_dbsf_db_execute('chado', 'gmod_dbsf_decypher_organism_id', $id);
  if (empty($org_ids)) {
    return "Organism $id not found.";
  }
  $access = gmod_dbsf_chadoapi_check_restriction('select', 'organism', $org_ids['organism_id']);
  if ($access != TRUE) {
    return "You are not allowed to access this resource.";
  }
  //class, order, family, genus, species, organism_id, common_name, ncbi_taxid, description (i.e. binomial)
  $org_data = gmod_dbsf_db_execute('chado', 'gmod_dbsf_get_add_organism', array('organism_id' => $org_ids['organism_id'], 'pgtype' => 'chado'));
  return drupal_get_form('genes4all_curate_organism_curate_form', $org_data);
}

function genes4all_curate_organism_curate_form($form_state, $org) {
  $form = array();
  $form['organism_id'] = array(
    '#type' => 'value',
    '#value' => $org['organism_id'],
  );
  $form['data'] = array(
    '#type' => 'fieldset',
    '#title' => 'Phylogeny',
    '#description' => t('The following data are derived from NCBI and stored for your organism. You may change only the common name.'),
    '#tree' => TRUE,
    'class' => array(
      '#type' => 'textfield',
      '#title' => 'Class',
      '#weight' => 0,
      '#default_value' => $org['class'],
      '#disabled' => TRUE,
    ), 'order' => array(
      '#type' => 'textfield',
      '#title' => 'Order',
      '#weight' => 1,
      '#default_value' => $org['order'],
      '#disabled' => TRUE,
    ), 'family' => array(
      '#type' => 'textfield',
      '#title' => 'family',
      '#weight' => 1,
      '#default_value' => $org['family'],
      '#disabled' => TRUE,
    ), 'genus' => array(
      '#type' => 'textfield',
      '#title' => 'Genus',
      '#weight' => 1,
      '#default_value' => $org['genus'],
      '#disabled' => TRUE,
    ), 'species' => array(
      '#type' => 'textfield',
      '#title' => 'Species',
      '#weight' => 1,
      '#default_value' => $org['species'],
      '#disabled' => TRUE,
    ), 'common_name' => array(
      '#type' => 'textfield',
      '#title' => 'Common name',
      '#weight' => 1,
      '#default_value' => $org['common_name'],
    ), 'description' => array(
      '#type' => 'textfield',
      '#title' => 'Latin binomial name',
      '#weight' => 1,
      '#default_value' => $org['description'],
      '#disabled' => TRUE,
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

function genes4all_curate_organism_curate_form_validate($form, &$form_state) {
  foreach ($form_state['values']['data'] as $key => $value) {
    $form_state['values']['data'][$key] = check_plain(trim($value));
  }
}

function genes4all_curate_organism_curate_form_submit($form, &$form_state) {
  $org_id = $form_state['values']['organism_id'];
  $org_data = gmod_dbsf_db_execute('chado', 'gmod_dbsf_get_add_organism', array('organism_id' => $org_id, 'pgtype' => 'chado'));
  foreach ($form_state['values']['data'] as $key => $value) {
    if ($value !== $org_data[$key] && $key == 'common_name') {
      gmod_dbsf_db_execute('chado_edit', 'gmod_dbsf_get_add_organism',
        array('organism_id' => $org_id, 'common_name' => $value, 'pgtype' => 'chado', 'update' => TRUE)
      );
    }
  }
}

function genes4all_curate_page_organism_add() {
  $exec = drupal_get_path('module', 'gmod_dbsf') .'/scripts/get_species_ncbi.pl';
  if (!file_exists($exec)) {
    drupal_set_message(t('Cannot find the get_species_ncbi.pl program.'), 'error');
    return FALSE;
  }
  elseif (!is_executable($exec)) {
    drupal_set_message(t('Failed to execute %exec', array('%exec' => $exec)), 'error');
    return FALSE;
  }
  $return_text = "This page allows you to create a new organism.";
  $add_form = drupal_get_form('genes4all_curate_add_organism_form');
  if (strlen($add_form) < 400) {
    $add_form = '';
  }
  $return_text .= $add_form;
  return $return_text;
}

function genes4all_curate_view_organism_page() {
  return drupal_get_form('genes4all_curate_view_organism_form');
}

function genes4all_curate_view_organism_form($form_state) {
  $form = array();

  $form['data'] = array(
    '#type' => 'fieldset',
    '#title' => 'View an organism',
    'genus' => array(
      '#type' => 'textfield',
      '#title' => 'Genus name',
      '#weight' => 0,
      '#required' => TRUE,
      '#autocomplete_path' => 'genes4all/curate/genus/autocomplete',
    ), 'species' => array(
      '#type' => 'textfield',
      '#title' => 'Species name',
      '#weight' => 1,
      '#required' => TRUE,
      '#autocomplete_path' => 'genes4all/curate/species/autocomplete',
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

function genes4all_curate_view_organism_form_validate($form, &$form_state) {
  $genus   = trim($form_state['values']['genus']);
  $species = trim($form_state['values']['species']);
  $org_ids = gmod_dbsf_db_execute('chado', 'gmod_dbsf_decypher_organism_id', $genus .'_'. $species);
  if (empty($org_ids)) {
    form_set_error('data', t('Organism %genus %species does not exists, has no data or you have no access to.', array('%genus' => $genus, '%species' => $species)));
    return;
  }
  $form_state['values']['org_id'] = $org_ids['organism_id'];
  $form_state['values']['genus'] = $genus;
  $form_state['values']['species'] = $species;
}

function genes4all_curate_view_organism_form_submit($form, &$form_state) {
  $org_id = $form_state['values']['org_id'];
  $form_state['redirect'] = array('genes4all/organism/curate/'. $form_state['values']['genus'] .'_'. $form_state['values']['species']);
  return "genes4all/organism/curate/". $form_state['values']['genus'] .'_'. $form_state['values']['species'];
}

function genes4all_curate_add_organism_form($form_state) {
  $form = array();

  $form['data'] = array(
    '#type' => 'fieldset',
    '#title' => 'Add new organism',
    'ncbi_taxid' => array(
      '#type' => 'textfield',
      '#title' => 'NCBI taxonomy ID',
      '#description' => t('Please give the '. l('NCBI taxonomy', 'http://www.ncbi.nlm.nih.gov/Taxonomy/', array('attributes' => array('target' => '_blank'))) .' ID. Information about the names of the species, family, order and class will then be populated automatically.' . ' If your species is not in NCBI yet, then you must contact NCBI to acquire an ID.'
      ),
      '#size' => 6,
      '#weight' => 0,
      '#required' => TRUE,
    ),
  );

  $form['button'] = array(
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Submit data'),
      '#weight' => 10,
    ), 'clear' => array(
      '#type' => 'button',
      '#value' => t('Clear form'),
      '#weight' => 11,
      '#validate' => array('gmod_dbsf_form_clear'),
    ),
  );
  return $form;
}

function genes4all_curate_add_organism_form_validate($form, &$form_state) {
  $ncbi_taxid = check_plain(trim($form_state['values']['ncbi_taxid']));
  if (empty($ncbi_taxid) || !is_numeric($ncbi_taxid)) {
    form_set_error('ncbi_taxid', t('You must provide a number for the NCBI Taxonomy ID.'));
  }

  $form_state['values']['ncbi_taxid'] = $ncbi_taxid;
}

function genes4all_curate_add_organism_form_submit($form, &$form_state) {
  $ncbi_taxid = $form_state['values']['ncbi_taxid'];
  //to add to admin menu.
  $operations[] = array('genes4all_curate_add_organism_batch', array($ncbi_taxid));

  $batch = array(
    'title' => t('Fetching data for your submission...'),
    'operations' => $operations,
    'init_message' => t('Fetching data for your submission...'),
    'progress_message' => t('@remaining operations remaining...'),
    'error_message' => t('Your submission encountered an error.'),
    'finished' => 'genes4all_curate_add_organism_batch_finished',
    'file' => drupal_get_path('module', 'genes4all_curate') .'/includes/genes4all_curate_organism.inc',
  );
  batch_set($batch);
}

function genes4all_curate_add_organism_batch($ncbi_taxid, &$context) {
  if (empty($ncbi_taxid)) {
    return;
  }
  $tax_db = gmod_dbsf_get_add_var('gmod_dbsf_taxonomy_directory');
  $exec   = drupal_get_path('module', 'gmod_dbsf') .'/scripts/get_species_ncbi.pl';
  $args   = " -n $ncbi_taxid";
  if (!empty($tax_db) && file_exists($tax_db .'/names.dmp')) {
    $args .= " -flat $tax_db";
  }
  $output   = array();
  $taxonomy = array();
  $result   = exec($exec . $args .' 2>&1', $output, $return_var);
  if (empty($result)) {
    //dpm(array($exec . $args .' 2>&1', $output));
    drupal_set_message(t('Failed to find any entry with ID %ncbi_taxid.', array('%ncbi_taxid' => $ncbi_taxid)), 'error');
    return FALSE;
  }
  if (preg_match('/class:/', $result)) {
    $result = trim($result);
    // TODO i'm sure there is a php function for this.
    $ln_data = explode(';', $result);
    foreach ($ln_data as $data) {
      $values = explode(':', $data);
      $taxonomy[$values[0]] = $values[1];
    }
  }
  else {
    //dpm(array($exec . $args .' 2>&1', $output));
    drupal_set_message(t('Failed to find any entry with ID %ncbi_taxid.', array('%ncbi_taxid' => $ncbi_taxid)), 'error');
    return FALSE;
  }
  if (empty($taxonomy['class']) || empty($taxonomy['order']) || empty($taxonomy['family'])
    || empty($taxonomy['genus']) || empty($taxonomy['species'])
  ) {
    drupal_set_message(t('Failed to find a complete taxonomy ID %ncbi_taxid.', array('%ncbi_taxid' => $ncbi_taxid)));
    //dpm($taxonomy);
    return FALSE;
  }
  $org_id = gmod_dbsf_db_execute('chado_edit', 'gmod_dbsf_get_add_organism', array(
      'pgtype' => 'chado',
      'common_name' => $taxonomy['common'],
      'species' => $taxonomy['species'],
      'genus' => $taxonomy['genus'],
      'class' => $taxonomy['class'],
      'order' => $taxonomy['order'],
      'family' => $taxonomy['family'],
      'ncbi_taxid' => $taxonomy['ncbi'],
    )
  );
  $context['results']['org_id'] = $org_id;
  $context['results']['taxonomy'] = $taxonomy;
}

function genes4all_curate_add_organism_batch_finished($success, $results, $operations) {
  if (!empty($success)) {
    $org_id   = $results['org_id'];
    $taxonomy = $results['taxonomy'];
    $tax      = $taxonomy['genus'] .' '. $taxonomy['species'] .' (class: '. $taxonomy['class'] .', order: '. $taxonomy['order'] .', family: '. $taxonomy['family'] .')';
    if (!empty($org_id)) {
      $curate_link = 'gmod_dbsf/organism/curate/NCBI'. $taxonomy['ncbi'];
      drupal_set_message(t('Organism %tax is in the Chado database. You can ', array('%tax' => $tax)) . l(t('curate it here'), $curate_link), 'warning');
    }
    else {
      // has the NCBI_TAXONOMY been populated properly?
      drupal_set_message(t('Could not add organism %gen %spe to the database.', array('%gen' => $taxonomy['genus'], '%spe' => $taxonomy['species'])), 'error');
  //    dpm($taxonomy);
    }
  }
  else {
    $error_operation = reset($operations);
    $message = t('An error occurred while processing your request. Please '. l('contact', 'contact') .' the administrator.');
    //dpm($error_operation);
    drupal_set_message($message, 'warning');
    return FALSE;
  }
}

