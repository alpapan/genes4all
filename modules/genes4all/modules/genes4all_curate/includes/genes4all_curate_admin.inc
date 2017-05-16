<?php
// $Id$

/**
 * @file
 * include file for admin of curation module
 */
function genes4all_curate_admin_page() {
  $return_text = 'The genes4all curation centre allows users with curation permissions to add and/or edit new features in your Chado database. ' . 'The submissions are first stored in your local Drupal database until a senior curator approves them. Further, it' . " utilizes the Drupal Bioinformatic Server Framework and the module biosoftware_bench to handle automatic annotation with InterProScan and NCBI's BLASTALL. " . 'Please use the '. l('BioSoftware Bench administration page', 'admin/settings/bench') .' to activate BLAST datasets and choose here which ones will be used for automatic annotation.';

  $db_data = drupal_get_form('genes4all_curate_admin_page_form');
  if (strlen($db_data) < 600) {
    $db_data = '';
  }

  $page_tabs = array(
    'datasets' => array(
      '#type' => 'tabset',
      'available' => array(
        '#type' => 'tabpage',
        '#title' => 'Basic settings',
        '#content' => $db_data,
        '#weight' => 0,
      ),
    ),
  );
  $return_text .= tabs_render($page_tabs);
  $already_linked = gmod_dbsf_get_linked_resources('Dataset for auto-annotation');
  if (!empty($already_linked)) {
    $links = array();
    foreach ($already_linked as $link_data) {
      $links[] = '<strong>'. $link_data['uniquename'] .'</strong>';
    }
    $return_text .= '<h3>Linked datasets</h3><p>These datasets have been already assigned for auto-annotation: '. implode(', ', $links) .'</p>';
  }
  return $return_text;
}

function genes4all_curate_admin_page_form($form_state) {
  $form = array();
  $database_array     =gmod_dbsf_get_add_var('genes4all_my_central_federation');
  $dbcode = array_pop(array_keys($database_array));
  $dbname = array_pop(array_values($database_array));
  $form['federation']=array(
    '#type' => 'fieldset',
    '#title' => t('Central Federation'),
    '#description' => t('Please register your Central database.'),
    '#weight' => 0,
    'dbcode' => array(
    '#type'=>'textfield',
    '#default_value'=>$dbcode,
    '#size'=> 5,
    '#maxlength'=>2,
    '#title'=>t('Code'),
    '#description'=>t('A 2 letter code to define your database according to the est2assembly format (e.g. IC74343Econ1 where IC is the code.'),
    '#weight'=>2,
    ),
    'dbname'=>array(
    '#type'=>'textfield',
    '#default_value'=>$dbname,
    '#title'=>t('Database name'),
    '#description'=>t('Please name your central database, eg. InsectaCentral.'),
    '#weight'=>0,
  ),
  'federate' => array(
    '#type' => 'submit',
    '#value' => t('Register Central'),
    '#weight' => 5,
  ),

    );
  if (empty($dbs)) {
    return $form;
  }
  $dbs = biosoftware_bench_get_user_datasets('blastall');
  $size = count($dbs) < 10 ? count($dbs) : 10;
    $form['data'] = array(
    '#type' => 'fieldset',
    '#title' => t('Availabe datasets'),
    '#tree' => TRUE,
    '#description' => t('Please specify which datasets will be used for curated feature annotation. You can use your Ctrl/option key to select multiple datasets or deselect them.'),
    '#weight' => 0,
    'proteins' => array(
      '#type' => 'select',
      '#multiple' => TRUE,
      '#size' => $size,
      '#title' => 'Protein datasets',
      '#options' => $dbs['protein'],
    ), 'nucleotide' => array(
      '#type' => 'select',
      '#multiple' => TRUE,
      '#size' => $size,
      '#title' => 'Nucleotide datasets',
      '#options' => $dbs['nucleotide'],
    ),
  'db_change' => array(
    '#type' => 'submit',
    '#value' => t('Add dataset'),
    '#weight' => 1,
  ),
  'db_remove' => array(
    '#type' => 'submit',
    '#value' => 'Remove dataset',
    '#weight' => 2,
  ));
  return $form;
}

function genes4all_curate_admin_page_form_validate($form, &$form_state) {
  switch($form_state['values']['op']){
    case 'Register Central':
      if (empty($form_state['values']['dbname'])){
        form_set_error('name',t('You must provide a name for your Central'));
      }if (empty($form_state['values']['dbcode'])){
        form_set_error('name',t('You must provide a code for your Central'));
      }
      break;
    case 'Remove dataset':
    $linked = gmod_dbsf_get_linked_resources('Dataset for auto-annotation');
    if (empty($linked)) {
      form_set_error('data', t('You cannot unlink datasets when you have not linked any yet!'));
      return;
    }
    foreach ($form_state['values']['data'] as $data_type => $value_array) {
      if (empty($value_array)) {
        continue;
      }
      foreach ($value_array as $key => $value) {
        if (empty($linked[$key])) {
          form_set_error($data_type, t('You asked to unlink a dataset which is not linked.'));
        }
      }
    }
    break;
    case 'Add dataset':
    $linked = gmod_dbsf_get_linked_resources('Dataset for auto-annotation');
    foreach ($form_state['values']['data'] as $data_type => $value_array) {
      if (empty($value_array)) {
        continue 1;
      }
      foreach ($value_array as $key => $value) {
        if (!empty($linked[$key])) {
          form_set_error($data_type, t('You asked to link a dataset which is already linked.'));
        }
      }
    }
    break;
  }
}

function genes4all_curate_admin_page_form_submit($form, &$form_state) {
  $res_id   = gmod_dbsf_get_add_resource('Dataset for auto-annotation', FALSE, FALSE, TRUE);
  $res_name = implode(',', array_keys($res_id));
  $res_id   = implode(',', $res_id);
  switch($form_state['values']['op']){
      case 'Register Central':
        $dbname = check_plain($form_state['values']['dbname']);
        $dbcode = check_plain($form_state['values']['dbcode']);
        gmod_dbsf_get_add_var('genes4all_my_central_federation',array($dbcode=>$dbname));
        drupal_set_message(t('Code set to %c and name set to %n',array('%c'=>$dbcode,'%n'=>$dbname)),'warning');
        break;
  case 'Remove dataset':
    foreach ($form_state['values']['data'] as $data_type => $value_array) {
      if (empty($value_array)) {
        continue;
      }
      foreach ($value_array as $id => $value) {
        $res = gmod_dbsf_link_unlink_resources($res_id, 'linked_with', $id, 'DELETE');
        if ($res === TRUE) {
          drupal_set_message(t('Dataset has been removed from %res_name', array('%res_name' => $res_name)), 'warning');
        }
      }
    }
    break;
  case 'Add dataset':
    foreach ($form_state['values']['data'] as $data_type => $value_array) {
      if (empty($value_array)) {
        continue;
      }
      foreach ($value_array as $id => $value) {
        $res = gmod_dbsf_link_unlink_resources($res_id, 'linked_with', $id);
        if ($res === TRUE) {
          drupal_set_message(t('Dataset has been set as %res_name', array('%res_name' => $res_name)), 'warning');
        }
      }
    }
  break;
  }
}

function genes4all_curate_access_restrictions_page($show_all = NULL) {

  $return = "<p>The genes4all software allows you to restrict certain datasets so that only certain users and groups have access to them." . " These restrictions can be anchored to a library, an organism or just a particular feature. If a dataset has no restrictions, then it will be viewable to all." . " These restrictions are cumulative: access to a feature is restricted if the entire organism is restricted or if the library it belongs to. If you restrict an organism" . ", therefore, no unpriviliged users will be allowed to see any of the data for that organism.</p>";

  if (empty($show_all)) {
    $db_data = drupal_get_form('genes4all_curate_access_restrictions_form');
    if (strlen($db_data) > 600) {
      $return .= $db_data;
    }
    $return .= genes4all_curate_access_restrictions_table();
  }
  else {
    // just table, no form.
    $return .= genes4all_curate_access_restrictions_table('all');
  }
  return $return;
}

function genes4all_curate_access_restrictions_form($form_state) {
  static $ic_phylogeny_loaded = FALSE;

  $form = array();
  $restr_chosen = $form_state['storage']['restriction_level'];

  $form['data'] = array(
    '#type' => 'fieldset',
    '#title' => 'Access restrictions',
  );
  if (empty($restr_chosen)) {

    $restr_level = array('Organism' => 'Organism', 'Library' => 'Library', 'Feature' => 'Feature');
    $form['data']['restriction_level'] = array(
      '#weight' => -1,
      '#title' => 'Restriction level',
      '#description' => t('Please select at which level the restrictions will apply.'),
      '#type' => 'select',
      '#size' => 1,
      '#options' => $restr_level,
    );
    $form['data']['choose'] =array(
        '#type' => 'submit',
        '#value' => t('Choose'),
        '#weight' => 1,
    );
  }
  else {
    $restr_types = gmod_dbsf_db_execute('chado', 'gmod_dbsf_get_add_cv_withnames', 'ic_data_restrictions');
    switch ($restr_chosen) {
      case 'Organism':
        $org_array = gmod_dbsf_db_execute('chado', 'gmod_dbsf_getspecies');
        $form['data']['species'] = gmod_dbsf_db_execute('chado', 'gmod_dbsf_phylogeny_tree', $org_array['phylogeny_array']);
        break;

      case 'Library':
        $lib_data = gmod_dbsf_db_execute('chado', 'gmod_dbsf_get_add_library');
        $libs = array();
        foreach ($lib_data as $lib_id => $data) {
          $libs[$data['org_abbr']][$lib_id] = $data['library.name'];
        }
        unset($lib_data);
        $form['data']['libraries'] = array(
          '#type' => 'select',
          '#multiple' => TRUE,
          '#size' => 10,
          '#title' => 'Libraries',
          '#options' => $libs,
        );
        break;

      case 'Feature':
        $form['data']['features'] = array(
          '#type' => 'textfield',
          '#title' => 'Features',
          '#autocomplete_path' => 'genes4all/feature/autocomplete',
        );
        break;
    }
    $user_roles = gmod_dbsf_get_users_roles('BOTH');
    $form['data']['user'] = array(
      '#weight' => 5,
      '#title' => 'Which users to restrict',
      '#description' => t('Please select one or more user to allow access. There is no need to specifically select a user if you allow access to one group they belongs to.'),
      '#type' => 'select',
      '#size' => 10,
      '#multiple' => TRUE,
      '#options' => $user_roles['users'],
    );
    $form['data']['role'] = array(
      '#weight' => 4,
      '#title' => 'Which roles to restrict',
      '#description' => t('Please select one or more roles to allow access. You can '. l('add more roles', 'admin/user/roles') .' to your Drupal installation.'),
      '#type' => 'select',
      '#size' => 7,
      '#multiple' => TRUE,
      '#options' => $user_roles['roles'],
    );

    $form['data']['restriction_type'] = array(
      '#weight' => 6,
      '#title' => 'Restriction type',
      '#description' => t('Please select a type of access restriction. The user/role will be allowed to perform these operations and everyone else will be prevented.'),
      '#type' => 'select',
      '#options' => $restr_types,
    );
    $form['data']['restriction_level'] = array(
      '#disabled' => TRUE,
      '#weight' => -1,
      '#title' => 'Restriction level',
      '#description' => t('Please select at which level the restrictions will apply.'),
      '#type' => 'select',
      '#options' => array($restr_chosen => $restr_chosen),
    );

    $form[] = array(
      '#weight' => 10,
      '#type' => 'fieldset',
      '#title' => 'Submit',
      '#description' => t('You can grant or revoke permissions for the settings you have selected above. It is okey to attempt to revoke permissions from a dataset which had none, or to grant to a user which already has them.'),
      'grant' => array(
        '#type' => 'submit',
        '#value' => t('Grant'),
        '#weight' => 1,
      ),
      'revoke' => array(
        '#type' => 'submit',
        '#value' => 'Revoke',
        '#weight' => 2,
      ),
    );
  }

  return $form;
}

function genes4all_curate_access_restrictions_form_validate($form, &$form_state) {

  if (empty($form_state['storage'])) {
    $form_state['rebuild'] = TRUE;
    $form_state['storage']['restriction_level'] = $form_state['values']['restriction_level'];
  }
  else {
    $form_state['values']['restriction_level'] = $form_state['storage']['restriction_level'];
    if ($form_state['storage']['restriction_level'] == 'Organism') {
      $org_array = gmod_dbsf_db_execute('chado', 'gmod_dbsf_getspecies');
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
    }
    $form_state['rebuild'] = FALSE;
    unset($form_state['storage']);
    if (empty($form_state['values']['species_select']) && empty($form_state['values']['libraries']) && empty($form_state['values']['features'])) {
      form_set_error('data', t('You must give some input!'));
    }
  }
}

function genes4all_curate_access_restrictions_form_submit($form, &$form_state) {
  $data_array  = array_merge($form_state['values']['user'], $form_state['values']['role']);
  $restr_type  = $form_state['values']['restriction_type'];
  $restr_level = $form_state['values']['restriction_level'];
  $is_chado    = gmod_dbsf_is_chado();
  switch ($form_state['values']['op']) {
    case 'Grant':
      $data_str = gmod_dbsf_array2csv($data_array, TRUE);
      $str = '';
      switch ($form_state['values']['restriction_level']) {
        case 'Organism':
          //;
          foreach ($form_state['values']['species_select'] as $id) {
              gmod_dbsf_db_execute('chado_edit', 'gmod_dbsf_populate_prop', 'organism', $id, $restr_type, $data_array);
            }
            break;

          case 'Library':
            foreach ($form_state['values']['libraries'] as $id) {
                gmod_dbsf_db_execute('chado_edit', 'gmod_dbsf_populate_prop', 'library', $id, $restr_type, $data_array);
              }
              break;

            case 'Feature':
              gmod_dbsf_db_execute('chado_edit', 'gmod_dbsf_populate_prop', 'feature', $form_state['values']['features'], $restr_type, $data_array);
              break;
          }
          drupal_set_message(t('Granted %restr_level access for %data_str', array('%restr_level' => $restr_level, '%data_str' => $data_str)), 'warning');
          break;

        case 'Revoke':
          $data_str   = gmod_dbsf_array2csv($data_array, TRUE);
          $str        = '';
          //$sql_revoke = "DELETE from %s WHERE type_id=$restr_type AND value IN ($data_str) AND ";
          gmod_dbsf_db_execute('chado_edit');
          switch ($form_state['values']['restriction_level']) {
            case 'Organism':
              $str = gmod_dbsf_array2csv($form_state['values']['species_select']);
              db_query($sql_revoke ." organism_id IN ($str)", 'organismprop');
              break;

            case 'Library':
              $str = gmod_dbsf_array2csv($form_state['values']['libraries']);
              db_query($sql_revoke ." library_id IN ($str)", 'libraryprop');
              break;

            case 'Feature':
              $str = gmod_dbsf_array2csv($form_state['values']['features']);
              db_query($sql_revoke ." feature_id IN ($str)", 'featureprop');
              break;
          }
          gmod_dbsf_db_execute();
          drupal_set_message(t('Revoked %restr_level access for %data_str with %str', array('%restr_level' => $restr_level, '%data_str' => $data_str, '%str' => $str)), 'warning');
          break;
      }
    }

function genes4all_curate_access_restrictions_table($all = NULL) {
      // The following are CHADO SQLs so no {} fields
      $select_org_restriction_sql = "select o.genus,o.species,o.abbreviation as uniquename,o.common_name as name,op.value as has_access,cvterm.name as access_type " . " from organismprop op JOIN organism o ON o.organism_id=op.organism_id JOIN cvterm ON cvterm.cvterm_id=op.type_id " . " where cvterm.cv_id=(select cv_id from cv where name='ic_data_restrictions') ORDER by uniquename,access_type,has_access";
      $select_lib_restriction_sql = "select l.uniquename,l.name,lp.value as has_access,cvterm.name as access_type " . " from libraryprop lp JOIN library l ON l.library_id=lp.library_id JOIN cvterm ON cvterm.cvterm_id=lp.type_id " . " where cvterm.cv_id=(select cv_id from cv where name='ic_data_restrictions') ORDER by uniquename,access_type,has_access";
      $select_feature_restriction_sql = "select f.name,f.uniquename, fp.value as has_access,cvterm.name as access_type" . " from featureprop fp JOIN feature f ON f.feature_id=fp.feature_id JOIN cvterm ON cvterm.cvterm_id=fp.type_id " . " where cvterm.cv_id=(select cv_id from cv where name='ic_data_restrictions') ORDER by uniquename,access_type,has_access";
      $table_rows = array();
      gmod_dbsf_db_execute('chado');
      $res = db_query($select_org_restriction_sql);
      while ($row = db_fetch_array($res)) {
        $uname = $row['uniquename'];
        if (empty($uname)) {
          $uname = $row['genus'] .' '. $row['species'];
        }
        if (!empty($row['name'])) {
          $uname .= ' ('. $row['name'] .')';
        }
        $table_rows[] = array('Organism', $uname, $row['access_type'], $row['has_access']);
      }
      $res = db_query($select_lib_restriction_sql);
      while ($row = db_fetch_array($res)) {
        $uname = $row['uniquename'];
        if (!empty($row['name'])) {
          $uname .= ' ('. $row['name'] .')';
        }
        $table_rows[] = array('Library', $uname, $row['access_type'], $row['has_access']);
      }
      if (!empty($all)) {
        $res = db_query($select_feature_restriction_sql);
        while ($row = db_fetch_array($res)) {
          $uname = $row['uniquename'];
          if (!empty($row['name'])) {
            $uname .= ' ('. $row['name'] .')';
          }
          $table_rows[] = array('Feature', $uname, $row['access_type'], $row['has_access']);
        }
      }
      unset($res);
      gmod_dbsf_db_execute();

      $table_headers = array('Table', 'Object name (alias)', 'Restriction type', 'User/Role',
      );
      $table = theme_table($table_headers, $table_rows, array('class' => 'genes4all_explore-library-table'), 'Existing restrictions <font size=-5>(except feature due to speed; see '. l('here', 'admin/genes4all_curate/access/all') .')</font>');
      return $table;
    }

