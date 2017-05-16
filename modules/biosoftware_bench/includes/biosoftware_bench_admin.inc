<?php
// $Id$

/**
 * @file
 * Administration interface for biosoftware_bench
 */

/**
 * Implements of hook_form()
 *
 * Used when selecting which software to make available
 * @see biosoftware_bench_software_form_submit()
 * @ingroup forms
 *
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_software_form($form_state) {
  $software_array  = array();
  $software_select = 'SELECT uniquename,description,active from {gmod_dbsf_software} as software';
  $res             = db_query($software_select);
  while ($row = db_fetch_array($res)) {
    $software_array[$row['uniquename']]['description'] = $row['description'];
    $software_array[$row['uniquename']]['active'] = $row['active'];
  }

  $form = array();
  $form['software'] = array(
    '#type' => 'fieldset',
    '#title' => t('Available software'),
    '#description' => t('Please specify which software will be active.'),
    '#tree' => TRUE,
    '#weight' => 1,
  );
  foreach ($software_array as $software => $value) {
    if ($value['active'] == 'f') {
      $value['active'] = FALSE;
    }
    else {
      $value['active'] = TRUE;
    }
    $description = t('@d', array('@d' => $value['description']));
    $form['software'][$software] = array(
      '#type' => 'checkbox',
      '#title' => $software,
      '#description' => $description,
      '#default_value' => $value['active'],
    );
  }
  $form['other_settings'] = array(
    '#type' => 'fieldset',
    '#title' => t('Other options'),
    '#tree' => FALSE,
    '#weight' => 2,
    'cron' => array(
      '#type' => 'textfield',
      '#size' => 3,
      '#default_value' => 7,
      '#required' => TRUE,
      '#title' => 'Cron interval (days)',
      '#description' => t('Every how many days should the module purge old files?'),
    ), 'upload_size' => array(
      '#type' => 'textfield',
      '#size' => 5,
      '#default_value' => 10,
      '#required' => TRUE,
      '#title' => 'Upload limit (Mb)',
      '#description' => t('The maximum file size (megabytes per file) that users can upload for similarity searches?'),
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
 * Implements hook_validate
 */
function biosoftware_bench_software_form_validate($form, &$form_state) {
  $cron = check_plain(trim($form_state['values']['cron']));
  if (!empty($cron) && !is_numeric($cron)) {
    form_set_error('cron', t('Cron days must be numeric.'));
  }
  $form_state['values']['cron'] = $cron;
  $upload_size = check_plain(trim($form_state['values']['upload_size']));
  if (!empty($upload_size) && !is_numeric($upload_size)) {
    form_set_error('upload_size', t('Upload limit must be numeric.'));
  }
  $form_state['values']['upload_size'] = $upload_size;
}

/**
 * Implements hook_form_submit()
 *
 * @see biosoftware_bench_software_form()
 *
 * @param $form
 *   hook_form data
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_software_form_submit($form, &$form_state) {
  //alter 1: for each option, activate into the _software table.
  //alter 2: If a software is present in _software table but not present in database, inactivate it
  $options        = $form_state['values']['software'];
  $cron           = $form_state['values']['cron'];
  $upload_size    = $form_state['values']['upload_size'];
  $sql_activate   = "UPDATE {gmod_dbsf_software} set active=TRUE where uniquename='%s'";
  $sql_deactivate = "UPDATE {gmod_dbsf_software} set active=FALSE where uniquename='%s'";
  foreach ($options as $name => $active) {
    if ($active === 1) {
      db_query($sql_activate, $name);
      drupal_set_message(t('Software %name is active.', array('%name' => $name)), 'warning');
    }
    else {
      db_query($sql_deactivate, $name);
      drupal_set_message(t('Software %name is not active.', array('%name' => $name)), 'warning');
    }
  }
  if (!empty($cron)) {
    gmod_dbsf_get_add_var('biosoftware_bench_cron_cleanup', $cron);
    drupal_set_message(t('Old files will be purged every %cron days.', array('%cron' => $cron)), 'warning');
  }
  if (!empty($upload_size)) {
    gmod_dbsf_get_add_var('biosoftware_bench_upload_size', $upload_size);
    drupal_set_message(t('Users will be able to upload files up to %mb megabytes.', array('%mb' => $upload_size)), 'warning');
  }
}



/**
 * Builds a table with all files which could be added as datasets
 *
 * @return
 *   HTML of a table
 */
function biosoftware_bench_admin_database_files_table() {
  $software_avail_array = biosoftware_bench_get_software(array('active' => TRUE), TRUE);
  $directory_array = gmod_dbsf_get_add_resource('directory');
  if (empty($directory_array)) {
    drupal_set_message(t('It seems that no directory paths have been registered yet. Will show no databases.'), 'error');
    return FALSE;
  }
  //array with [software]['dataset_suffix'][rank]=value where rank==0 is switch ANY or ALL
  $dataset_suffixes = gmod_dbsf_get_softwareprop(array('cvterm.name' => 'dataset_suffix', 'active' => TRUE), TRUE);
  if (empty($dataset_suffixes)) {
      return FALSE;
  }

  $table_caption    = 'Files in your registered directories which could act as databases';
  $table_header     = array('Directory', 'Filename', 'Modified date', 'Formatted for', 'Suffixes found');
  $table_rows       = array();
  foreach ($directory_array as $id => $dirname) {
    $dirname = rtrim($dirname, '/');
    foreach ($dataset_suffixes as $software_name => $data) {
      $suffixes_needed = $data['dataset_suffix'][0] === 'ALL' ? count($data['dataset_suffix']) - 1 : 1;
      $suffixes_found  = 0;
      $files_array     = array();
      foreach ($data['dataset_suffix'] as $key => $suffix) {
        if ($key === 0) {
          continue;
          //rank 0 is switch
        }
        $f = file_scan_directory($dirname, '.+'. $suffix .'$', array('.', '..', 'CVS'), NULL, FALSE, 'basename');
        foreach ($f as $fname => $obj) {
          $fname = preg_replace("/$suffix$/", '', $fname, 1);
          $files_array[$fname][] = $suffix;
        }
      }
      foreach ($files_array as $basename => $f_data) {
        $suffixes_found = count($f_data);
        if ($suffixes_found == $suffixes_needed) {
          $suffix_str   = implode(', ', $f_data);
          $mtime        = filemtime($dirname .'/'. $basename . $f_data[0]);
          $mtime        = date("F j, Y, g:i a", $mtime);
          $table_rows[] = array($dirname, $basename, $mtime, $software_name, $suffix_str);
        }
      }
    }
    //end software
  }
  //end dir
  usort($table_rows, create_function('$a,$b', 'return strcmp($a[1], $b[1]);'));
  $table = theme_table($table_header, $table_rows, array('width' => '750px', 'class' => 'vertical'), $table_caption);
  return $table;
}

/**
 * Builds a table containing all datasets that have been uploaded
 *
 * @param $software
 *  Optionally, software to create the table for. Defaults to NULL i.e. all
 *
 * @return
 *   HTML of a table
 */
function biosoftware_bench_admin_database_table($software = NULL) {
  $software_avail_array = biosoftware_bench_get_software(array('active' => TRUE), TRUE);
  if (!empty($software) && empty($software_avail_array[$software])) {
    drupal_set_message(t('It seems that %software is not active yet. Will not show registered databases for it.', array('%software' => $software)), 'error');
    return FALSE;
  }
  elseif (empty($software_avail_array)) {
    drupal_set_message(t('It seems that no software is not active yet. Will not show any registered databases.'), 'error', FALSE);
    return FALSE;
  }
  $table_caption = 'Databases currently registered';
  $table_header  = array('Linked with', 'Role permitted', 'File', 'Name shown to users', 'Grouping', 'Database type', 'NCBI TaxID');
  $table_rows    = array();
  // get all datasets
  if (empty($software)) {
    $dataset_software = gmod_dbsf_get_software_resources(array('rcvterm.name' => 'dataset'), 'RESOURCE');
    $dataset_array = gmod_dbsf_get_add_resource('dataset', FALSE, FALSE, 'ALL');
    if (!empty($dataset_array)) {
      foreach ($dataset_array as $dataset_id => $dataset_data) {
        $dataset_name = $dataset_data['uniquename'];
        $ncbi_tax     = $dataset_data['organism_id'] > 1 ? gmod_dbsf_getspecies($dataset_data['organism_id']) : 0;
        $ncbi_tax     = !empty($ncbi_tax['species_array']) ? implode('', array_values($ncbi_tax['species_array'])) : 'None';
        $linked_with  = '';
        if (!empty($dataset_software) && !empty($dataset_software[$dataset_name])) {
          foreach ($dataset_software[$dataset_name] as $software_name => $d) {
            $linked_with .= $software_name .' ';
          }
        }
        $linked_with = trim($linked_with);
        if (empty($linked_with)) {
          $linked_with = 'none';
        }
        $to_add    = array();
        $role_str  = 'none';
        $role_data = biosoftware_bench_get_resource_roles($dataset_id);
        if (!empty($role_data)) {
          $role_str = implode(', ', $role_data);
        }
        $dataset_extra   = gmod_dbsf_get_resource_cvterm_and_prop($dataset_id);
        $dbtype          = implode(', ', $dataset_extra['cvterm']['dataset_type']);
        $linked_res_data = gmod_dbsf_get_linked_resources($dataset_id, NULL, TRUE);
        foreach ($linked_res_data as $object_name => $object_data) {
          $object_type = $object_data['type_name'];
          $to_add[$object_type] = $object_name;
        }

        $table_rows[] = array($linked_with, $role_str, $to_add['directory'] . $to_add['filename'], $dataset_name, $to_add['dataset_group'], $dbtype, $ncbi_tax);
      }
    }
  }
  else {
    $table_caption .= " to $software";
    $table_header = array('Role permitted', 'File', 'Name shown to users', 'Grouping', 'Database type');
    foreach ($software_avail_array as $software_name => $software_id) {
      if (!empty($software) && $software_name !== $software) {
        continue;
      }
      $software_data = gmod_dbsf_get_software_resources(
        array('sr.software_id' => $software_id, 'rcvterm.name' => 'dataset'),
        FALSE
      );
      if (empty($software_data[$software_id])) {
        continue;
      }
      foreach ($software_data[$software_id] as $dataset_id => $dataset_data) {
        $to_add       = array();
        $dataset_name = $dataset_data['uniquename'];
        $role_data    = biosoftware_bench_get_resource_roles($dataset_id);
        if (!empty($role_data)) {
          $role_str = implode(', ', $role_data);
        }
        $dataset_extra   = gmod_dbsf_get_resource_cvterm_and_prop($dataset_id);
        $dbtype          = implode(', ', $dataset_extra['cvterm']['dataset_type']);
        $linked_res_data = gmod_dbsf_get_linked_resources($dataset_id, NULL, TRUE);
        foreach ($linked_res_data as $object_name => $object_data) {
          $object_type = $object_data['type_name'];
          $to_add[$object_type] = $object_name;
        }
        $table_rows[] = array($role_str, $to_add['directory'] . $to_add['filename'], $dataset_name, $to_add['dataset_group'], $dbtype);
      }
    }
  }
  //$table_rows[] = array($linked_with, $role_str, $to_add['directory'] . $to_add['filename'], $dataset_name, $to_add['dataset_group'], $dbtype);
  usort($table_rows, create_function('$a,$b', 'return strcmp($a[4], $b[4]);'));

  $table = theme_table($table_header, $table_rows, array('width' => '750px', 'class' => 'vertical'), $table_caption);
  return $table;
}

/**
 * Produce HTML for the 'administrate dataset page'.
 *
 * Menu admin/bench/dataset
 */
function biosoftware_bench_admin_dataset_page() {
  $add_dataset_variables = drupal_get_form('biosoftware_bench_admin_dataset_variables_form');
  $add_dataset_data = drupal_get_form('biosoftware_bench_admin_newdataset_form');
  $dataset_perm = drupal_get_form('biosoftware_bench_admin_dataset_perms_form');
  $link_dataset = drupal_get_form('biosoftware_bench_admin_link_dataset_software_form');
  $edit_dataset = drupal_get_form('biosoftware_bench_admin_edit_dataset_form');
  $db_table = biosoftware_bench_admin_database_table();
  $file_table = biosoftware_bench_admin_database_files_table();


  if (strlen($link_dataset) < 600) {
    $link_dataset = '';
  }
  if (strlen($edit_dataset) < 600) {
    $edit_dataset = '';
  }
  if (strlen($add_dataset_data) < 600) {
    $add_dataset_data = '';
    drupal_set_message(t('You will need to define at least one group and one directory in <strong>Dataset variables</strong> before registering new datasets.'), 'error');
  }
  if (strlen($dataset_perm) < 600) {
    $dataset_perm = '';
  }
  if (!empty($dataset_perm)) {
    $dataset_perm .= $db_table;
  }
  if (!empty($link_dataset)) {
    $link_dataset .= $db_table;
  }
  if (!empty($add_dataset_data)) {
    $add_dataset_data .= $file_table;
  }
  if (!empty($edit_dataset)) {
    $edit_dataset .= $db_table;
  }
  $page_tabs = array(
    'dataset' => array(
      '#type' => 'tabset',
      'variables' => array(
        '#type' => 'tabpage',
        '#title' => 'Dataset variables',
        '#content' => $add_dataset_variables,
        '#weight' => 4,
      ),
      'add' => array(
        '#type' => 'tabpage',
        '#title' => 'Available datasets',
        '#content' => $add_dataset_data,
        '#weight' => 1,
      ),
      'edit' => array(
        '#type' => 'tabpage',
        '#title' => 'Edit dataset',
        '#content' => $edit_dataset,
        '#weight' => 2,
      ),
      'perm' => array(
        '#type' => 'tabpage',
        '#title' => 'User permissions',
        '#content' => $dataset_perm,
        '#weight' => 3,
      ),
      'link' => array(
        '#type' => 'tabpage',
        '#title' => 'Link datasets',
        '#content' => $link_dataset,
        '#weight' => 5,
      ),
    ),
  );
  $return_text = tabs_render($page_tabs);
  return $return_text;
}

/**
 * Implements of hook_form()
 *
 * Used to edit dataset properties
 * @see biosoftware_bench_admin_edit_dataset_form_validate()
 * @see biosoftware_bench_admin_edit_dataset_form_submit()
 * @ingroup forms
 *
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_admin_edit_dataset_form($form_state) {
  $form = array();
  // select of all dataset friendly names.
  $dataset_array = gmod_dbsf_get_add_resource('dataset', NULL, FALSE, 'BOTH');
  //allow to change:
  // select of new group
  $group_array = array(0 => 'No Change') + (array)gmod_dbsf_get_add_resource('dataset_group');
  if (empty($dataset_array) || empty($group_array)) {
    return $form;
  }
  $form['data'] = array(
    '#type' => 'fieldset',
    '#title' => 'Edit dataset',
    '#description' => t('Here you can change the group membership of a  dataset or delete a database.'),
    'current_name' => array(
      '#type' => 'select',
      '#size' => 10,
      '#title' => 'Dataset to change',
      '#options' => $dataset_array,
      '#weight' => 0,
      '#required' => TRUE,
      '#multiple' => TRUE,
    ),
    'friendly_name' => array(
      '#type' => 'textfield',
      '#size' => 50,
      '#title' => 'New friendly name',
      '#description' => t('This is shown to users.'),
      '#weight' => 1,
    ),
    'group_id' => array(
      '#weight' => 3,
      '#type' => 'select',
      '#description' => t('Select the new group.'),
      '#size' => 1,
      '#title' => 'Group',
      '#options' => $group_array,
    ), 'ncbi_id' => array(
      '#weight' => 4,
      '#type' => 'textfield',
      '#description' => t('Optionally, assign an NCBI Tax ID.'),
      '#size' => 8,
      '#title' => 'NCBI Tax ID',
    ),
  );

  $form['db_change'] = array(
    '#type' => 'submit',
    '#value' => t('Alter DB'),
    '#weight' => 1,
  );
  $form['db_remove'] = array(
    '#type' => 'submit',
    '#value' => 'Remove DB',
    '#weight' => 2,
  );
  return $form;
}

/**
 * Implements of hook_form()
 *
 * Used to edit dataset properties
 * @see biosoftware_bench_admin_edit_dataset_form()
 * @see biosoftware_bench_admin_edit_dataset_form_submit()
 * @ingroup forms
 *
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_admin_edit_dataset_form_validate($form, &$form_state) {
  $fullname = check_plain(trim($form_state['values']['friendly_name']));
  $form_state['values']['friendly_name'] = $fullname;
  $group_id = check_plain($form_state['values']['group_id']);
  $form_state['values']['group_id'] = check_plain($group_id);
  $current_name_array = $form_state['values']['current_name'];
  $current_name_array2 = array();
  $current_name = '';
  foreach ($current_name_array as $data) {
    $current_name_array2[] = check_plain($data);
  }
  $current_name_array = $current_name_array2;
  unset($current_name_array2);
  $form_state['values']['current_name'] = $current_name_array;
  if (count($current_name_array) == 1) {
    $current_name = implode('', $current_name_array);
  }

  $ncbi_id = check_plain(trim($form_state['values']['ncbi_id']));
  switch ($form_state['clicked_button']['#value']) {
    case 'Alter DB':
      if (!empty($fullname) && empty($current_name)) {
        form_set_error('current_name', t('You can only rename one database at a time!.'));
      }
      elseif (!empty($fullname)) {
        $res = gmod_dbsf_get_add_resource('dataset', array('uniquename' => $fullname));
        if (!empty($res)) {
          form_set_error('friendly_name', t('Name %fullname has been registered already.', array('%fullname' => $fullname)));
        }
      }
      if (!empty($ncbi_id) && !is_numeric($ncbi_id)) {
        form_set_error('ncbi_id', t('NCBI ID needs to be a numeric identifier from NCBI.'));
      }
      elseif (!empty($ncbi_id)) {
        //get taxonomy for ncbi_id
        $org_data = gmod_dbsf_get_taxonomy_from_ncbi($ncbi_id);
        $org_id = $org_data['organism_id'];
        if (empty($org_id)) {
          form_set_error('ncbi_id', t('There has been a problem with getting the organism data. Please report the error.'));
        }
        $form_state['values']['org_id'] = $org_id;
      }
      break;

    case 'Remove DB':
      foreach ($current_name_array as $current_name) {
        $res = gmod_dbsf_get_add_resource('dataset', array('uniquename' => $current_name));
        if (empty($res)) {
          form_set_error('name', t('Database %current_name does not exist.', array('%current_name' => $current_name)));
        }
      }
      break;
  }
}

/**
 * Implements of hook_form()
 *
 * Used to edit dataset properties
 * @see biosoftware_bench_admin_edit_dataset_form()
 * @see biosoftware_bench_admin_edit_dataset_form_validate()
 * @ingroup forms
 *
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_admin_edit_dataset_form_submit($form, &$form_state) {
  $fullname           = check_plain(trim($form_state['values']['friendly_name']));
  $group_id           = check_plain($form_state['values']['group_id']);
  $org_id             = $form_state['values']['org_id'];
  $current_group_id   = '';
  $filename           = '';
  $current_name_array = $form_state['values']['current_name'];
  foreach ($current_name_array as $current_name) {

    $linked_res_data = gmod_dbsf_get_linked_resources($current_name);
    foreach ($linked_res_data as $object_id => $object_data) {
      if ($object_data['type_name'] == 'dataset_group') {
        $current_group_id = $object_id;
      }
      elseif ($object_data['type_name'] == 'filename') {
        $filename = $object_data['uniquename'];
      }
    }
    $change_org_sql = "UPDATE {gmod_dbsf_resource} as resource set organism_id=%d WHERE uniquename='%s'";
    switch ($form_state['clicked_button']['#value']) {
      case 'Alter DB':
        if (!empty($group_id) && $group_id != $current_group_id) {
          $success = gmod_dbsf_link_unlink_resources($current_name, 'uses', $group_id, 'REPLACE');
          if (!empty($success)) {
            drupal_set_message(t('Dataset %current_name group membership changed.', array('%current_name' => $current_name)), 'warning');
          }
          else {
            drupal_set_message(t('For some reason, group membership of dataset %current_name could not be changed (%group_id).', array('%current_name' => $current_name, '%group_id' => $group_id)), 'error');
          }
        }
        if (!empty($fullname)) {
          $success = gmod_dbsf_rename_resource($current_name, $fullname);
          if (!empty($success)) {
            drupal_set_message(t('Dataset %current_name renamed to %fullname.', array('%current_name' => $current_name, '%fullname' => $fullname)), 'warning');
          }
          else {
            drupal_set_message(t('For some reason, name of dataset %current_name could not be changed.', array('%current_name' => $current_name)), 'error');
          }
        }
        if (!empty($org_id)) {
          db_query($change_org_sql, $org_id, $current_name);
          db_query($change_org_sql, $org_id, $filename);
          drupal_set_message(t("Dataset $current_name assigned to a new organism."), 'warning');
        }
        if (empty($org_id) && empty($fullname) && (empty($group_id) || $group_id == $current_group_id)) {
          drupal_set_message(t('No change was performed on dataset %current_name.', array('%current_name' => $current_name)), 'warning');
        }
        break;

      case 'Remove DB':
        $dataset_data = gmod_dbsf_get_add_resource('dataset', array('uniquename' => $current_name), 'DELETE');
        $filename_data = gmod_dbsf_get_add_resource('filename', array('uniquename' => $filename), 'DELETE');
        if ($dataset_data === TRUE) {
          drupal_set_message(t('Dataset %current_name removed.', array('%current_name' => $current_name)), 'warning');
        }
        break;
    }
  }
}

/**
 * Produce HTML for the 'administrate software page'.
 *
 * Menu admin/bench/software
 * Grabs available plugins as well.
 * @see biosoftware_bench_software_form()
 * @see biosoftware_bench_software_blastall_core_settings_form()
 */
function biosoftware_bench_admin_software_page() {
  $software_data = drupal_get_form('biosoftware_bench_software_form');
  $page_tabs = array(
    'software' => array(
      '#type' => 'tabset',
      'available' => array(
        '#type' => 'tabpage',
        '#title' => 'Available software',
        '#content' => $software_data,
        '#weight' => 0,
      ),
    ),
  );

  // plugins
  $plugin_path = drupal_get_path('module', 'biosoftware_bench') .'/includes/';
  $plugins = array();
  //PLUGIN_EDIT_HERE
  if (file_exists($plugin_path .'biosoftware_bench_blastall.inc')) {
    $plugins['blastall']['file'] = $plugin_path .'biosoftware_bench_blastall.inc';
  }
  if (file_exists($plugin_path .'biosoftware_bench_ssaha2.inc')) {
    $plugins['ssaha2']['file'] = $plugin_path .'biosoftware_bench_ssaha2.inc';
  }
  if (file_exists($plugin_path .'biosoftware_bench_iprscan.inc')) {
    $plugins['iprscan']['file'] = $plugin_path .'biosoftware_bench_iprscan.inc';
  }
  if (file_exists($plugin_path .'biosoftware_bench_annot8r.inc')) {
    $plugins['annot8r']['file'] = $plugin_path .'biosoftware_bench_annot8r.inc';
  }
  foreach ($plugins as $plugin_name => $plugin_data) {
    require_once($plugin_data['file']);
    $data = '';
    if (function_exists('biosoftware_bench_software_'. $plugin_name .'_core_settings_form')) {
      $data = drupal_get_form('biosoftware_bench_software_'. $plugin_name .'_core_settings_form');
      if (strlen($data) < 300) {
        $data = '';
      }
      if (!empty($data)) {
        $page_tabs['software'][$plugin_name] = array(
          '#type' => 'tabpage',
          '#title' => strtoupper($plugin_name),
          '#content' => $data,
          '#weight' => 5,
        );
      }
    }
  }
  $return_text = tabs_render($page_tabs);
  return $return_text;
}

/**
 * Implements hook_form()
 *
 * @see biosoftware_bench_admin_dataset_perms_form_validate()
 * @see biosoftware_bench_admin_dataset_perms_form_submit()
 * @ingroup forms
 *
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_admin_dataset_perms_form($form_state) {
  $form = array();
  //keys are resource_ids and values are the user-friendly names
  $dataset_array = gmod_dbsf_get_add_resource('dataset');
  $role_array = biosoftware_bench_get_role();
  if (empty($dataset_array)) {
    return $form;
  }
  if (empty($role_array)) {
    drupal_set_message(t('This is very weird, it seems I cannot get the roles from your Drupal installation. Please file a bug report.'));
    return $form;
  }
  $form['access'] = array(
    '#type' => 'fieldset',
    '#title' => "Data access settings",
    '#description' => t('If you wish, you can protect datasets (e.g. databases for BLAST or other software) by asking that they are accessible only to users belonging to specific user roles.' . ' By default, all new datasets are accessible only to logged in users (authenticated). If you wish for  all users to be able to see the database, add permissions for the anonymous' . ' group too; alternatively, you may remove permissions from the authenticated group and add them to a new private one you created via Drupal.'),
    'role' => array(
      '#type' => 'select',
      '#size' => 1,
      '#title' => 'Role',
      '#description' => t('These roles are defined by '. l("Drupal's role management", 'admin/user/roles') .'.'),
      '#options' => $role_array,
    ),
    'dataset' => array(
      '#type' => 'select',
      '#size' => 10,
      '#multiple' => TRUE,
      '#title' => 'Dataset',
      '#description' => t('Which dataset(s) to protect.'),
      '#options' => $dataset_array,
    ),
    'grant_submit' => array(
      '#type' => 'submit',
      '#weight' => 10,
      '#value' => t('Grant permission'),
    ), 'revoke_submit' => array(
      '#type' => 'submit',
      '#weight' => 11,
      '#value' => t('Revoke permission'),
    ),
  );
  return $form;
}

/**
 * Implements hook_form_submit()
 *
 * @see biosoftware_bench_admin_dataset_perms_form()
 * @see biosoftware_bench_admin_dataset_perms_form_validate()
 *
 * @param $form
 *   hook_form data
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_admin_dataset_perms_form_submit($form, &$form_state) {
  $role_id = $form_state['values']['role'];
  $dataset_id_array = $form_state['values']['dataset'];
  if (empty($dataset_id_array) || empty($role_id)) {
    return FALSE;
  }
  switch ($form_state['values']['op']) {
    case 'Grant permission':
      foreach ($dataset_id_array as $dataset_id) {
        $result = gmod_dbsf_link_unlink_resource_role($dataset_id, $role_id, 'GRANT');
        switch ($result) {
          case FALSE:
            drupal_set_message(t('Failed to change dataset permissions'), 'error');
            break;

          case 'GRANTED':
            drupal_set_message(t('Permission granted.'), 'warning');
            break;
        }
      }
      break;

    case 'Revoke permission':
      foreach ($dataset_id_array as $dataset_id) {
        $result = gmod_dbsf_link_unlink_resource_role($dataset_id, $role_id, 'REVOKE');
        switch ($result) {
          case FALSE:
            drupal_set_message(t('Failed to change dataset permissions'), 'error');
            break;

          case 'REVOKED':
            drupal_set_message(t('Permission revoked.'), 'warning');
            break;
        }
      }
      break;
  }
}

/**
 * Implements hook_form()
 *
 * @see biosoftware_bench_admin_dataset_variables_form_validate()
 * @see biosoftware_bench_admin_dataset_variables_form_submit()
 *
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_admin_dataset_variables_form($form_state) {
  $form = array();
  $access = user_access('Add new dataset paths');
  if (!empty($access)) {
    $form['directories'] = array(
      '#type' => 'fieldset',
      '#title' => "Dataset/DB paths",
      'directory_name' => array(
        '#type' => 'textfield',
        '#description' => t('Here you can register a directory on the server that contains one or more databases for BLAST or other software.'),
        '#size' => 50,
        '#title' => "Full path to directory containing your datasets",
      ),
      'directory_submit' => array(
        '#type' => 'submit',
        '#value' => t('Register directory'),
      ),
    );
  }
  else {
    $form['directories'] = array(
      '#type' => 'fieldset',
      '#title' => "Dataset/DB paths",
      '#description' => t('For security reasons, your role does not permit adding new dataset paths. You will have to use the existing ones or ') . l('contact an administrator', 'contact') .'.',
    );
  }
  $form['groups'] = array(
    '#type' => 'fieldset',
    '#title' => "File Groups",
    'group_name' => array(
      '#type' => 'textfield',
      '#description' => t('In order to help users when viewing a large number of datasets, you can create groups for them.'),
      '#size' => 20,
      '#title' => "Group name",
    ),
    'group_submit' => array(
      '#type' => 'submit',
      '#value' => t('Create group'),
    ),
  );
  return $form;
}

/**
 * Implements hook_form_validate()
 *
 * @see biosoftware_bench_admin_dataset_variables_form()
 * @see biosoftware_bench_admin_dataset_variables_form_submit()
 *
 * @param $form
 *   hook_form data
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_admin_dataset_variables_form_validate($form, &$form_state) {
  switch ($form_state['clicked_button']['#value']) {
    case 'Register directory':
      $directory_name = check_plain(trim($form_state['values']['directory_name']));
      if (empty($directory_name)) {
        form_set_error('directory_name', t('If you want to register a new directory, you must provide a name!'));
      }
      if (!file_exists($directory_name)) {
        form_set_error('directory_name', t('Cannot find directory %directory_name.', array('%directory_name' => $directory_name)));
      }
      elseif (substr($directory_name, -1) !== '/') {
        $directory_name .= '/';
        $form_state['values']['directory_name'] = $directory_name;
      }
      break;

    case 'Create group':
      $group_name = check_plain(trim($form_state['values']['group_name']));
      if (empty($group_name)) {
        form_set_error('group_name', t('If you want to create a new group, you must provide a name!'));
      }
      else {
        $form_state['values']['group_name'] = $group_name;
      }
      break;
  }
}

/**
 * Implements hook_form_submit()
 *
 * @see biosoftware_bench_admin_dataset_variables_form_validate()
 * @see biosoftware_bench_admin_dataset_variables_form()
 *
 * @param $form
 *   hook_form data
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_admin_dataset_variables_form_submit($form, &$form_state) {
  switch ($form_state['clicked_button']['#value']) {
    case 'Register directory':
      $access = user_access('Add new dataset paths');
      if (empty($access)) {
        return FALSE;
      }
      $directory_name = $form_state['values']['directory_name'];
      if (empty($directory_name)) {
        return FALSE;
      }
      $res = gmod_dbsf_get_add_resource('directory', array('uniquename' => $directory_name, 'dbxref_name' => 'local resource'), 'ADD');
      if (!empty($res)) {
        drupal_set_message(t('Directory %directory_name saved.', array('%directory_name' => $directory_name)), 'warning');
      }
      break;

    case 'Create group':
      $group_name = $form_state['values']['group_name'];
      if (empty($group_name)) {
        return FALSE;
      }
      gmod_dbsf_get_add_resource('dataset_group', array('uniquename' => $group_name, 'dbxref_name' => 'local resource'), 'ADD');
      if (!empty($res)) {
        drupal_set_message(t('Group %group_name saved.', array('%group_name' => $group_name)), 'warning');
      }
      break;
  }
}

/**
 * Implements hook_form()
 *
 * @see biosoftware_bench_admin_newdataset_form_validate()
 * @see biosoftware_bench_admin_newdataset_form_submit()
 * @ingroup forms
 *
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_admin_newdataset_form($form_state) {
  $form            = array();
  $directory_array = gmod_dbsf_get_add_resource('directory');
  $group_array     = gmod_dbsf_get_add_resource('dataset_group');
  $dbtypes_array   = gmod_dbsf_get_add_cv_withnames('dataset_type', NULL, 'BOTH');
  if (empty($directory_array)) {
    return $form;
  }
  if (empty($group_array)) {
    return $form;
  }

  $form['data'] = array(
    '#type' => 'fieldset',
    '#title' => 'Database settings',
    '#description' => t('You can register new datasets to the system by providing the filename and the directory on the server. In order to remove datasets, just provide the <strong>friendly name</strong> you gave before and click <strong>remove DB</strong>.'),
    'friendly_name' => array(
      '#type' => 'textfield',
      '#size' => 50,
      '#title' => 'Friendly name',
      '#description' => t('This will be shown to users instead of the filename'),
      '#weight' => 2,
      '#required' => TRUE,
    ),
    'filename' => array(
      '#weight' => 1,
      '#type' => 'textfield',
      '#size' => 100,
      '#required' => TRUE,
      '#title' => t('Filename'),
      '#description' => t('Please provide the basename (no directory) of the database (e.g. as you would pass it to the -d option of BLASTALL).' . ' <strong>Avoid any suffix created by a particular software</strong> (e.g. for BLAST, .pal/.nal can be created).'
      ),
    ),
    'db_type' => array(
      '#weight' => 3,
      '#type' => 'checkboxes',
      '#title' => 'Molecule type',
      '#description' => t('Select what molecule types are supported by this basename. You can select <i>both</i> if, for example, a nucleotide and protein BLAST database share the same basename.'),
      '#options' => $dbtypes_array,
      '#required' => TRUE,
    ),
    'directory_id' => array(
      '#weight' => -1,
      '#type' => 'select',
      '#size' => 1,
      '#title' => 'Directory',
      '#options' => $directory_array,
    ),
    'group_id' => array(
      '#weight' => 4,
      '#type' => 'select',
      '#size' => 1,
      '#title' => 'Group',
      '#options' => $group_array,
      '#required' => TRUE,
    ),
    'ncbi_id' => array(
      '#description' => t('Optionally, you may link an NCBI description to your dataset. Please provide the '. l('NCBI Taxonomy', 'http://www.ncbi.nlm.nih.gov/Taxonomy/') .' ID.'),
      '#weight' => 10,
      '#type' => 'textfield',
      '#size' => 8,
      '#title' => 'NCBI ID',
    ),
  );

  $form['db_add'] = array(
    '#type' => 'submit',
    '#value' => t('Add DB'),
    '#weight' => 1,
  );

  return $form;
}

/**
 * Implements hook_form_validate()
 *
 * @see biosoftware_bench_admin_newdataset_form()
 * @see biosoftware_bench_admin_newdataset_form_submit()
 *
 * @param $form
 *   hook_form data
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_admin_newdataset_form_validate($form, &$form_state) {
  // security checks
  $db_types = $form_state['values']['db_type'];
  foreach ($db_types as $key => $type) {
    if (empty($type)) {
      unset($db_types[$key]);
    }
  }
  $db_type = implode(',', $db_types);
  $filename = check_plain(trim($form_state['values']['filename']));
  $form_state['values']['filename'] = $filename;
  $fullname = check_plain(trim($form_state['values']['friendly_name']));
  $form_state['values']['friendly_name'] = $fullname;
  $directory_id = check_plain($form_state['values']['directory_id']);
  $form_state['values']['directory_id'] = $directory_id;
  $group_id = check_plain($form_state['values']['group_id']);
  $form_state['values']['group_id'] = $group_id;
  $blast_id = check_plain($form_state['values']['dataset']);
  $form_state['values']['dataset'] = $blast_id;
  $ncbi_id = check_plain(trim($form_state['values']['ncbi_id']));

  //one cmd to rule them all
  $sql_select = 'SELECT uniquename from {gmod_dbsf_resource} as resource WHERE resource_id=%d';
  //$dbtypes_array=gmod_dbsf_get_add_cv_withnames('dataset_type');
  if (!empty($filename)) {
    if (empty($fullname)) {
      form_set_error('friendly_name', t('Please give a friendly name for your database.'));
    }
    $res = gmod_dbsf_get_add_resource('dataset', array('uniquename' => $fullname));
    if (!empty($res)) {
      form_set_error('friendly_name', t('Name %fullname has been registered already.', array('%fullname' => $fullname)));
    }
    if (empty($group_id)) {
      form_set_error('group_id', t('You must specify a group.'));
    }
    if (empty($directory_id)) {
      form_set_error('directory_id', t('You must specify a directory.'));
    }
    if (empty($db_types{'nucleotide'}) && empty($db_types{'protein'})) {
      form_set_error('db_type', t('You must specify a type of molecule.'));
      return;
    }
    // check for file existence
    $dbdata = db_fetch_array(db_query($sql_select, $directory_id));
    $filename = $dbdata['uniquename'] . $filename;
    //array with [software]['dataset_suffix'][rank]=value where rank==0 is switch ANY or ALL
    $dataset_suffixes = gmod_dbsf_get_softwareprop(array('cvterm.name' => 'dataset_suffix', 'active' => TRUE), TRUE);
    //TODO consider removing the following and going through $dataset_suffixes
    if ((file_exists($filename .'.nin') || file_exists($filename .'.nal'))
      && (file_exists($filename .'.pin') || file_exists($filename .'.pal'))
    ) {
      if (empty($db_types{'nucleotide'}) || empty($db_types{'protein'})) {
        form_set_error('filename', t('You asked for %db_type as molecule type but I found it formatted as both nucleotide and protein.', array('%db_type' => $db_type)));
        return FALSE;
      }
    }
    elseif (file_exists($filename .'.nin') || file_exists($filename .'.nal')) {
      if (empty($db_types{'nucleotide'})) {
        form_set_error('filename', t('You asked for %db_type as molecule type but I found it formatted as nucleotide.', array('%db_type' => $db_type)));
        return FALSE;
      }
    }
    elseif (file_exists($filename .'.pin') || file_exists($filename .'.pal')) {
      if (empty($db_types{'protein'})) {
        form_set_error('filename', t('You asked for %db_type as molecule type but I found it formatted as protein.', array('%db_type' => $db_type)));
        return FALSE;
      }
    }
    elseif (file_exists($filename) && (strpos($filename, '.pin') !== FALSE || strpos($filename, '.pal') !== FALSE || strpos($filename, '.nin') !== FALSE || strpos($filename, '.nal') !== FALSE)) {
      form_set_error('filename', t('I found %filename but you gave the formatdb suffix as well. Please give the base name without nin/pin/nal/pal suffices.', array('%filename' => $filename)));
      return FALSE;
    }
    /*PLUGIN_EDIT_HERE
     * Optionally you may include checks here to detect filename suffixes
     */
    else {
      $found = '';
      foreach ($dataset_suffixes as $software_name => $suffix_data) {
        $suffixes_needed = $suffix_data['dataset_suffix'][0] === 'ALL' ? count($suffix_data['dataset_suffix']) - 1 : 1;
        $suffixes_found = 0;
        foreach ($suffix_data['dataset_suffix'] as $rank => $suffix) {
          if ($rank == 0) {
            continue;
          }
          if (file_exists($filename . $suffix)) {
            $suffixes_found++;
            if ($suffixes_found == $suffixes_needed) {
              $found = $software_name;
              break2;
            }
          }
        }
      }
      if (!empty($found)) {
        drupal_set_message(t('I found %filename formatted for %found software.', array('%filename' => $filename, '%found' => $found)), 'warning');
      }
      elseif (file_exists($filename)) {
        drupal_set_message(t('I found %filename but I do not know if it is formatted for any particular software (if it is found in the table below, then it must be).', array('%filename' => $filename)), 'warning');
      }
      else {
        form_set_error('filename', t('Sorry, cannot find file %filename in your system.', array('%filename' => $filename)));
        return FALSE;
      }
    }
  }
  else {
    form_set_error('filename', t("Filename cannot be empty when you're trying to add a database."));
    return FALSE;
  }
  if (!empty($ncbi_id) && !is_numeric($ncbi_id)) {
    form_set_error('ncbi_id', t('NCBI ID needs to be a numeric identifier from NCBI.'));
  }
  elseif (!empty($ncbi_id)) {
    $org_data = gmod_dbsf_get_taxonomy_from_ncbi($ncbi_id);
    $org_id = $org_data['organism_id'];
    if (empty($org_id)) {
      form_set_error('ncbi_id', t('Did not find a species using %id as ID. Currently only species-level phylogenies are supported.', array('%id' => $ncbi_id)));
      return FALSE;
    }
    else {
      $form_state['values']['org_id'] = $org_id;
    }
  }
}

/**
 * Implements hook_form_submit()
 *
 * @see biosoftware_bench_admin_newdataset_form()
 * @see biosoftware_bench_admin_newdataset_form_validate()
 *
 * @param $form
 *   hook_form data
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_admin_newdataset_form_submit($form, &$form_state) {
  // security checks already performed during validate
  $filename     = $form_state['values']['filename'];
  $fullname     = $form_state['values']['friendly_name'];
  $directory_id = $form_state['values']['directory_id'];
  $group_id     = $form_state['values']['group_id'];
  $db_types     = $form_state['values']['db_type'];
  $org_id       = $form_state['values']['org_id'];
  //$dbtypes_array=gmod_dbsf_get_add_cv_withnames('dataset_type');
  if (empty($filename)) {
    return;
  }
  // add database name as the core resource
  $dataset_data = array();
  if (!empty($org_id)) {
    $dataset_data = gmod_dbsf_get_add_resource('dataset', array('organism_id' => $org_id, 'uniquename' => $fullname, 'dbxref_name' => 'similarity search database'), 'REPLACE', TRUE);
    $filename_data = gmod_dbsf_get_add_resource('filename', array('organism_id' => $org_id, 'uniquename' => $filename, 'dbxref_name' => 'similarity search database'), 'REPLACE', TRUE);
  }
  else {
    $dataset_data = gmod_dbsf_get_add_resource('dataset', array('uniquename' => $fullname, 'dbxref_name' => 'similarity search database'), 'REPLACE', TRUE);
    //add database filename, link it (use resource_ids as a security check
    $filename_data = gmod_dbsf_get_add_resource('filename', array('uniquename' => $filename, 'dbxref_name' => 'similarity search database'), 'REPLACE', TRUE);
  }
  gmod_dbsf_link_unlink_resources($dataset_data[$fullname], 'uses', $filename_data[$filename]);
  // directory and group already exist as resources, link them.
  gmod_dbsf_link_unlink_resources($dataset_data[$fullname], 'uses', $group_id);
  gmod_dbsf_link_unlink_resources($dataset_data[$fullname], 'uses', $directory_id);
  // set dbtype into resource_cvterm
  if (!empty($db_types{'protein'})) {
    gmod_dbsf_populate_cvterms('resource', $fullname, 'dataset_type', array('protein'));
  }
  if (!empty($db_types{'nucleotide'})) {
    gmod_dbsf_populate_cvterms('resource', $fullname, 'dataset_type', array('nucleotide'));
  }
  // link role
  $role_link = gmod_dbsf_link_unlink_resource_role($fullname, 'authenticated user');
  drupal_set_message(t('Dataset %filename saved.', array('%filename' => $filename)), 'warning');
}

/**
 * Implements hook_form()
 *
 * @see biosoftware_bench_admin_link_dataset_software_form_validate()
 * @see biosoftware_bench_admin_link_dataset_software_form_submit()
 * @ingroup forms
 *
 * @param $form_state
 *   hook_form data
 * @param $software
 *   Which software to produce form for
 */
function biosoftware_bench_admin_link_dataset_software_form($form_state, $software = NULL) {
  $form = array();
  $software_array = array();
  if (empty($software)) {
    $software_array = biosoftware_bench_get_software(array('active' => TRUE), 'BOTH');
    if (empty($software_array)) {
      drupal_set_message(t('Sorry, but you have not activated any software.'), 'error', FALSE);
      return $form;
    }
  }
  $dataset_array = gmod_dbsf_get_add_resource('dataset', NULL, FALSE, 'BOTH');
  if (empty($dataset_array)) {
    return $form;
  }
  $form['dataset_link'] = array(
    '#type' => 'fieldset',
    '#title' => "Link datasets with software",
    '#description' => t('Here you can link datasets you have submitted earlier with activated software'),
    'dataset' => array(
      '#options' => $dataset_array,
      '#size' => 5,
      '#title' => 'Dataset',
      '#description' => t('Which dataset/database to link.'),
      '#required' => TRUE,
    ),
    'submit_link' => array(
      '#type' => 'submit',
      '#value' => t('Link'),
      '#weight' => 10,
    ), 'submit_unlink' => array(
      '#type' => 'submit',
      '#value' => t('Unlink'),
      '#weight' => 11,
    ),
  );

  if (empty($software)) {
    $form['dataset_link']['dataset']['#type'] = 'select';
    $form['dataset_link']['dataset']['#multiple'] = TRUE;
    $form['dataset_link']['software'] = array(
      '#type' => 'select',
      '#size' => 1,
      '#title' => 'Software',
      '#description' => t('These are the activated software. Note that not all software require linking with datasets.'),
      '#options' => $software_array,
    );
  }
  else {
    $form['dataset_link']['dataset']['#type'] = 'checkboxes';
    $form['dataset_link']['software'] = array(
      '#type' => 'value',
      '#value' => $software,
    );
  }

  return $form;
}

/**
 * Implements hook_form_validate()
 *
 * @see biosoftware_bench_admin_link_dataset_software_form()
 * @see biosoftware_bench_admin_link_dataset_software_form_submit()
 *
 * @param $form
 *   hook_form data
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_admin_link_dataset_software_form_validate($form, &$form_state) {
  $dataset_array = $form_state['values']['dataset'];
  $software = check_plain($form_state['values']['software']);
  if (empty($dataset_array) || empty($software)) {
    form_set_error('', t('Something went wrong.'));
  }
  //array with [software]['dataset_suffix'][rank]=value where rank==0 is switch ANY or ALL
  $dataset_suffixes = gmod_dbsf_get_softwareprop(array('cvterm.name' => 'dataset_suffix', 'active' => TRUE), TRUE);
  if (empty($dataset_suffixes[$software])) {
    form_set_error('software', t('Software %software does not expect/allow dataset linking.', array('%software' => $software)));
    return FALSE;
  }
  foreach ($dataset_array as $key => $dataset) {
    $dataset = check_plain($dataset);
    $dataset_array[$key] = $dataset;
    if (empty($dataset)) {
      continue;
    }
    $format_check = 0;
    $fullpath     = biosoftware_bench_get_dataset_path($dataset);
    $basename     = basename($fullpath);
    $dirname      = dirname($fullpath);
    $f            = file_scan_directory($dirname, $basename .'.*$', array('.', '..', 'CVS'), 0, FALSE, 'basename');
    foreach ($f as $fname => $obj) {
      $suffix = str_replace($basename, '', $fname);
      if (in_array($suffix, $dataset_suffixes[$software]['dataset_suffix'])) {
        $format_check++;
      }
    }
    switch ($dataset_suffixes[$software]['dataset_suffix'][0]) {
      case 'ANY':
        if (empty($format_check)) {
          form_set_error('dataset', t("It does not seem that %dataset is properly formatted for %software.", array('%dataset' => $dataset, '%software' => $software)));
        }
        break;

      case 'ALL':
        if ($format_check < count($dataset_suffixes[$software]['dataset_suffix']) - 1) {
          form_set_error('dataset', t("It doesn't seem that %dataset is properly formatted for %software.", array('%dataset' => $dataset, '%software' => $software)));
        }
        break;
    }
  }
  $form_state['values']['software'] = $software;
  $form_state['values']['dataset'] = $dataset_array;
}

/**
 * Implements hook_form_submit()
 *
 * @see biosoftware_bench_admin_link_dataset_software_form()
 * @see biosoftware_bench_admin_link_dataset_software_form_validate()
 *
 * @param $form
 *   hook_form data
 * @param $form_state
 *   hook_form data
 */
function biosoftware_bench_admin_link_dataset_software_form_submit($form, &$form_state) {
  $dataset_array = $form_state['values']['dataset'];
  $software = $form_state['values']['software'];
  if (empty($dataset_array) || empty($software)) {
    return FALSE;
  }
  switch ($form_state['values']['op']) {
    case 'Link':
      foreach ($dataset_array as $dataset) {
        if (empty($dataset)) {
          continue;
        }
        $result = gmod_dbsf_link_unlink_software_resource($software, $dataset, 'GRANT');
        switch ($result) {
          case FALSE:
            drupal_set_message(t('Failed to link of %dataset with %software', array('%dataset' => $dataset, '%software' => $software)), 'error');
            break;

          case 'GRANTED':
            drupal_set_message(t('Dataset %dataset linked with %software.', array('%dataset' => $dataset, '%software' => $software)), 'warning');
            break;
        }
      }
      break;

    case 'Unlink':
      foreach ($dataset_array as $dataset) {
        if (empty($dataset)) {
          continue;
        }
        $result = gmod_dbsf_link_unlink_software_resource($software, $dataset, 'REVOKE');
        switch ($result) {
          case FALSE:
            drupal_set_message(t('Failed to unlink of %dataset with %software', array('%dataset' => $dataset, '%software' => $software)), 'error');
            break;

          case 'REVOKED':
            drupal_set_message(t('Dataset %dataset unlinked from %software.', array('%dataset' => $dataset, '%software' => $software)), 'warning');
            break;
        }
      }
      break;
  }
}

/**
 * Produce HTML for the 'administrate software page'.
 *
 * Menu admin/settings/bench
 */
function biosoftware_bench_admin_page() {
  $return_text = '<dl><dt>'. l('Setup the software available', 'admin/bench/software', array('attributes' => array('target' => '_blank'))) .'</dt><dd>Specify which software to activate for your search server</dd>';
  $return_text .= '<dt>'. l('Setup available reference datasets/databases', 'admin/bench/dataset', array('attributes' => array('target' => '_blank'))) .'</dt><dd>Load databases and datasets, specify which software and which users will have access to them.';
  $return_text .= '<dt>'. l('Setup the BLAST server', 'admin/bench/blastall', array('attributes' => array('target' => '_blank'))) ."</dt><dd>Specify the necessary variables and databases to use with NCBI's BLAST</dd>";
  //PLUGIN_EDIT_HERE
  if (file_exists(drupal_get_path('module', 'biosoftware_bench') .'/includes/biosoftware_bench_ssaha2.inc')) {
    $return_text .= '<dt>'. l('Setup the SSAHA2 server', 'admin/bench/ssaha2', array('attributes' => array('target' => '_blank'))) ."</dt><dd>Specify the necessary variables and databases to use with Sangers's SSAHA2</dd>";
  }
  if (file_exists(drupal_get_path('module', 'biosoftware_bench') .'/includes/biosoftware_bench_iprscan.inc')) {
    $return_text .= '<dt>'. l('Setup the InterProScan server', 'admin/bench/iprscan', array('attributes' => array('target' => '_blank'))) ."</dt><dd>Specify the necessary variables to use EBI's InterProScan</dd>";
  }
  if (file_exists(drupal_get_path('module', 'biosoftware_bench') .'/includes/biosoftware_bench_annot8r.inc')) {
    $return_text .= '<dt>'. l('Setup the annot8r server', 'admin/bench/annot8r', array('attributes' => array('target' => '_blank'))) ."</dt><dd>Specify the necessary variables and databases to use with Ralf Schmid's annot8r</dd>";
  }
  $return_text .= '</dl>';
  $return_text .= '<p>Once setup, your software will be available at from bench/ + "software short name"'
  .'. See this '.l('page','bench/select').' or '.' '.l('BLAST','bench/blastall')
  .', '.l('annot8r','bench/annot8r').', '.l('SSAHA2','bench/ssaha2').' and '.l('InterProScan','bench/iprscan')
  .'. We recommend you use the '.l('administration menu','admin/build/menu').' to setup up links to your software servers.</p>';
  return $return_text;
}

