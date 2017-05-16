<?php
// $Id$

/**
 * @file
 * Administration pages for gmod_dbsf
 */

/**
 * Prepare administration page
 */
function gmod_dbsf_admin_page() {
  return drupal_get_form('gmod_dbsf_admin_settings_form');
}

/**
 * Implement hook_form
 */
function gmod_dbsf_admin_settings_form($form_state) {
  $form             = array();
  $tax_db           = gmod_dbsf_get_add_var('gmod_dbsf_taxonomy_directory');
  $form['settings'] = array(
    '#type' => 'fieldset',
    '#title' => t('Basic settings'),
    '#tree' => FALSE,
    '#weight' => 1,
    'ncbi_taxdb' => array(
      '#type' => 'textfield',
      '#default_value' => $tax_db,
      '#title' => t('NCBI Taxonomy flatfiles'),
      '#description' => t('It is not necessary but in order to fetch NCBI Taxonomies faster, you can opt to install the database locally. Please give the directory where the NCBI taxonomy flatfiles (nodes.dmp, names.dmp etc) reside. They can be retrieved ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz .'),
    ),
  );
  $form['buttons'] = array(
    '#type' => 'submit',
    '#value' => t('Update'),
    '#weight' => 3,
  );
  return $form;
}

/**
 * Implement hook_form_validate
 */
function gmod_dbsf_admin_settings_form_validate($form, &$form_state) {
  if (!empty($form_state['values']['ncbi_taxdb'])) {
    $db = check_plain(trim($form_state['values']['ncbi_taxdb']));
    if (!file_exists($db)) {
      form_set_error('ncbi_taxdb', t('Directory %d does not exist!', array('%d' => $db)));
    }
    elseif (!is_writable($db)) {
      form_set_error('ncbi_taxdb', t('Cannot write in directory %d', array('%d' => $db)));
    }
    else {
      $form_state['values']['ncbi_taxdb'] = $db;
    }
  }
}

/**
 * Implement hook_form_submit
 */
function gmod_dbsf_admin_settings_form_submit($form, &$form_state) {
  if (!empty($form_state['values']['ncbi_taxdb'])) {
    $tax_db = gmod_dbsf_get_add_var('gmod_dbsf_taxonomy_directory', $form_state['values']['ncbi_taxdb']);
  }
}

