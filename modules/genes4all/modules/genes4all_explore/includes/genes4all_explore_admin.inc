<?php
// $Id$

/**
 * @file
 * include file for administating the exploration module
 */
function genes4all_explore_admin_page() {
  $rendered_form = drupal_get_form('genes4all_explore_admin_page_form');
  return $rendered_form;
}

function genes4all_explore_admin_page_form($form_state) {
  $form = array();
  $external_link_name = gmod_dbsf_get_add_var('genes4all_explore_external_search_name');
  $external_link      = gmod_dbsf_get_add_var('genes4all_explore_external_search_link');
  $gbrowse_url        = gmod_dbsf_get_add_var('genes4all_explore_gbrowse_link');
  $form['data']       = array(
    '#type' => 'fieldset',
    '#title' => t('Links'),
    '#tree' => FALSE,
    '#description' => t('Please specify links used by the module.'),
    'existing_external' => array(
      '#type' => 'value',
      '#value' => $external_link,
    ), 'existing_gbrowse' => array(
      '#type' => 'value',
      '#value' => $gbrowse_url,
    ),
    'external_search_link' => array(
      '#type' => 'textfield',
      '#title' => t('External Search Page'),
      '#description' => t('Optionally, you may specify an external link that is offered in the ') . l('search page', 'genes4all/explore'),
      '#weight' => 0,
      '#default_value' => $external_link,
    ), 'external_search_name' => array(
      '#type' => 'textfield',
      '#title' => t('Name for your external search page'),
      '#description' => t('Provide a name for your external link. It will be shown to users in the ') . l('search page', 'genes4all/explore'),
      '#weight' => 1,
      '#default_value' => $external_link_name,
    ), 'gbrowse_link' => array(
      '#type' => 'textfield',
      '#title' => t('URL for your GBrowse script'),
      '#description' => t('If you want to use GBrowse, please give the complete URL to your GBrowse CGI script, e.g. http://insectacentral.org/cgi-bin/gbrowse/ . Otherwise, leave this field empty.'),
      '#weight' => 1,
      '#default_value' => $gbrowse_url,
    ),
  );
  $form['buttons'] = array(
    '#type' => 'submit',
    '#value' => t('Submit'),
  );
  return $form;
}

function genes4all_explore_admin_page_form_validate($form, &$form_state) {
  $external_name = check_plain(trim($form_state['values']['external_search_name']));
  $external_link = check_url(trim($form_state['values']['external_search_link']));
  $gbrowse_url = check_url(trim($form_state['values']['gbrowse_link']));
  $form_state['values']['existing_external'] = check_url(trim($form_state['values']['existing_external']));
  $form_state['values']['existing_gbrowse'] = check_url(trim($form_state['values']['existing_gbrowse']));
  if (empty($gbrowse_url) && !module_exists('jbrowse')) {
    drupal_set_message('error', t('JBrowse not detected in your Drupal installation. In that case, you should provide a GBrowse link.'));
  }
  elseif (!empty($gbrowse_url) && substr($gbrowse_url, -1, 1) !== '/') {
    $gbrowse_url .= '/';
  }
  if (!empty($external_link) && empty($external_name)) {
    form_set_error('external_search_name', t('You must provide a name for your link.'));
    return FALSE;
  }
  elseif (empty($external_link) && !empty($external_name)) {
    form_set_error('external_search_link', t('You provided a link name but no link!'));
    return FALSE;
  }
  if (!empty($gbrowse_url) && strpos($gbrowse_url, 'http://') !== 0) {
    form_set_error('gbrowse_link', t('Links must be full, i.e. contain http:// in the beginning.'));
    return FALSE;
  }
  if (!empty($external_link) && strpos($external_link, 'http://') !== 0) {
    form_set_error('external_search_link', t('Links must be full, i.e. contain http:// in the beginning.'));
    return FALSE;
  }

  if (!empty($gbrowse_url)) {
    if ($gbrowse_url === $form_state['values']['existing_gbrowse']) {
      unset($gbrowse_url);
    }
    elseif (get_headers($gbrowse_url) === FALSE) {
      form_set_error('gbrowse_link', t('%l does not seem to be a valid link.', array('%l' => $gbrowse_url)));
      return FALSE;
    }
  }
  if (!empty($external_link)) {
    if ($external_link === $form_state['values']['existing_external']) {
      unset($external_link);
    }
    elseif (get_headers($external_link) === FALSE) {
      form_set_error('external_search_link', t('%l does not seem to be a valid link.', array('%l' => $external_url)));
      return FALSE;
    }
  }

  $form_state['values']['external_search_link'] = $external_link;
  $form_state['values']['external_search_name'] = $external_name;
  $form_state['values']['gbrowse_link'] = $gbrowse_url;
}

function genes4all_explore_admin_page_form_submit($form, &$form_state) {
  $external_name = $form_state['values']['external_search_name'];
  $external_link = $form_state['values']['external_search_link'];
  $gbrowse_url   = $form_state['values']['gbrowse_link'];
  if (!empty($external_name) && !empty($external_link)) {
    gmod_dbsf_get_add_var('genes4all_explore_external_search_name', $external_name);
    $success = gmod_dbsf_get_add_var('genes4all_explore_external_search_link', $external_link);
    if (!empty($success)) {
      drupal_set_message(t('Added an external link to %s.', array('%s' => $external_link)), 'warning');
    }
  }
  elseif (!empty($form_state['values']['existing_external'])) {
    gmod_dbsf_get_add_var('genes4all_explore_external_search_name', 'DELETE');
    $success = gmod_dbsf_get_add_var('genes4all_explore_external_search_link', 'DELETE');
    if (!empty($success)) {
      drupal_set_message('Deleted the external link.', 'warning');
    }
  }
  if (!empty($gbrowse_url) && $gbrowse_url !== $form_state['values']['existing_gbrowse']) {
    $success = gmod_dbsf_get_add_var('genes4all_explore_gbrowse_link', $gbrowse_url);
    if (!empty($success)) {
      drupal_set_message(t('Updated the GBrowse link to %s.', array('%s' => $gbrowse_url)), 'warning');
    }
  }elseif (empty($gbrowse_url) && !empty($form_state['values']['existing_gbrowse'])){
    $success = gmod_dbsf_get_add_var('genes4all_explore_gbrowse_link', 'DELETE');
    if (!empty($success)) {
     drupal_set_message(t('Deleted existing GBrowse link'), 'warning');
    }
  }
}

