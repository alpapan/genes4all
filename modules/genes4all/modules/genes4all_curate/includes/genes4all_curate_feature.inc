<?php
// $Id$

/**
 * @file
 * include file for curating features
 */


///////////////////////////////////////////////////////////////////////
// Feature
///////////////////////////////////////////////////////////////////////
function genes4all_curate_approve_features_form($form_state, $data) {
  $form = array();
  if (empty($data)) {
    return $form;
  }
  // $data contains all elements, we only want to show the wants the approval field is null.
  $sql_trim     = 'select uniquename,annotated from {genes4all_approval_queue} where approved is null';
  $not_approved = array();
  $res          = db_query($sql_trim);
  while ($row = db_fetch_array($res)) {
    $not_approved[$row['uniquename']] = $row['annotated'];
  }
  $form             = array('#type' => 'fieldset', '#title' => 'Features', '#tree' => TRUE);
  $img_status_done  = file_create_path() .'/genes4all_curate/status_done.gif';
  $img_status_wait  = file_create_path() .'/genes4all_curate/status_wait.gif';
  $img_status_error = file_create_path() .'/genes4all_curate/status_error.gif';
  $img_status_done  = theme_image($img_status_done, 'Annotation done', 'Annotation done', NULL, TRUE);
  $img_status_wait  = theme_image($img_status_wait, 'Waiting for annotation', 'Waiting for annotation', NULL, TRUE);
  $img_status_error = theme_image($img_status_error, 'Error with auto-annotation', 'Error with auto-annotation', NULL, TRUE);
  foreach ($data as $type => $features) {
    $tables = array();
    $table_helper = array();
    if (empty($features)) {
      continue;
    }
    foreach ($features as $feature_id => $feature_name) {
      $status = $not_approved[$feature_name];
      if (isset($status)) {
        $img = $img_status_error;
        switch ($status) {
          case '1':
            $img = $img_status_done;
            break;

          case '0':
            $img = $img_status_wait;
            break;
        }
        $table_helper['table_row'][$feature_id] = array($feature_name, l('Feature page', 'genes4all/curate/feature/view/'. $feature_name, array('attributes' => array('target' => '_blank'))), $img);
        //  key to be passed for checkboxes and value to show to user.
        $table_helper['row_check_data'][$feature_id] = $feature_name;
      }
    }
    if (!empty($table_helper)) {
      $table_header     = array('Name', 'Links', 'Status');
      $table_attributes = array('class' => 'padded-table sortable');
      $table_caption    = 'Overview of features in approval queue';
      $table_array      = array(
      //  '#links' => array('previous'=>$previous_link,'next'=>$next_link,'bookmark'=>$bookmark),
        'header' => $table_header,
        'caption' => $table_caption,
        'attributes' => $table_attributes,
        'data' => $table_helper['table_row'],
      );
      $form['data'][$type] = array('#value' => $table_array);
      $form['features'][$type] = array(
        '#type' => 'checkboxes',
        '#options' => $table_helper['row_check_data'],
      );
    }
  }
  $form['buttons'] = array(
    'accept' => array(
      '#type' => 'submit',
      '#value' => t('Approve'),
      '#weight' => 15,
  ), 'reject' => array(
      '#type' => 'submit',
      '#value' => t('Reject'),
      '#weight' => 16,
  ),
  );
  return $form;
}

function genes4all_curate_approve_features_form_validate($form, &$form_state) {
  // Approve or Reject
  $op = check_plain(trim($form_state['clicked_button']['#value']));
  global $user;
  $curator_name = check_plain($user->name);
  if (empty($curator_name)) {
    form_set_error('form', t('You must be logged in so we can track your curator actions.'));
    return;
  }
  $form_state['values']['curator_name'] = $curator_name;
  //make sure user has selected boxes with annotations done
  //this form will only check that the user has selected the entire annotation group.
  // if they have not, then give an error.
  //The following is a Drupal SQL
  $sql_check = "SELECT annotated from {genes4all_approval_queue} where uniquename=(SELECT uniquename from {gmod_dbsf_feature} where feature_id=%d)";
  $sql_parent_select_loc = "SELECT f.uniquename as parent_name,f.feature_id as parent_id from {gmod_dbsf_featureloc} fl JOIN " . " {gmod_dbsf_feature} f ON f.feature_id=fl.srcfeature_id where fl.feature_id=%d";
  $sql_children_select_loc = "SELECT f.uniquename as child_name,f.feature_id as child_id from {gmod_dbsf_featureloc} fl JOIN " . " {gmod_dbsf_feature} f ON f.feature_id=fl.feature_id where fl.srcfeature_id=%d";

  foreach ($form_state['values']['features'] as $features) {
    foreach ($features as $fid) {
      if (!empty($fid)) {
        $feature_list[] = $fid;
      }
    }
  }
  //we can get the srcfeature for each feature and check if it exists in Chado's feature or if it exists in the
  // feature list submitted.
  foreach ($feature_list as $fid) {
    $res = db_fetch_array(db_query($sql_check, $fid));
    if ($res['annotated'] < 1) {
      form_set_error('features', t('You are trying to operate on a feature whose annotations have not finished yet!'));
      return;
    }
    if ($res['annotated'] == 2 && $op == 'Approve') {
      form_set_error('features', t('You are trying to approve a feature whose annotations failed!'));
      return;
    }
    $parent_data_res = db_query($sql_parent_select_loc, $fid);
    while ($parent_data = db_fetch_array($parent_data_res)) {
      if (!empty($parent_data)) {
        //check feature_list
        if (!in_array($parent_data['parent_id'], $feature_list)) {
          //check chado
          $check = gmod_dbsf_db_execute('chado', 'gmod_dbsf_get_feature_name_id', $parent_data['parent_name'], TRUE);
          if (empty($check[$parent_data['parent_name']])) {
            //drupal_set_message("I also added the parents of a chosen feature.",'warning');
            $feature_list[] = $parent_data['parent_id'];
          }
        }
      }
    }
    $children_data_res = db_query($sql_children_select_loc, $fid);
    while ($children_data = db_fetch_array($children_data_res)) {
      if (!empty($children_data)) {
        //check feature_list
        if (!in_array($children_data['child_id'], $feature_list)) {
          //check chado
          $check = gmod_dbsf_db_execute('chado', 'gmod_dbsf_get_feature_name_id', $children_data['child_name'], TRUE);
          if (empty($check[$children_data['child_name']])) {
            $feature_list[] = $children_data['child_id'];
          }
        }
      }
    }
  }
}

function genes4all_curate_approve_features_form_submit($form, &$form_state) {
  // Approve or Reject
  $op = check_plain(trim($form_state['clicked_button']['#value']));

  $features         = $form_state['values']['features'];
  $dirpath          = file_create_path() .'/genes4all_curate/';
  $fullpath         = base_path() . $dirpath;
  $sql_select       = "SELECT annotation_id,approved,annotated from {genes4all_approval_queue} where uniquename='%s'";
  $gff_files        = array();
  $message          = '';
  $features_ordered = array('CDS' => (array)$features['CDS'], 'mRNA' => (array)$features['mRNA'], 'ORF' => (array)$features['ORF'],
    'polypeptide' => (array)$features['polypeptide'],
  );
  $features = $features_ordered;

  foreach ($features as $type => $features_array) {
    if (empty($features_array)) {
      continue;
    }
    foreach ($features_array as $fid) {
      $fid = check_plain(trim($fid));
      if (empty($fid)) {
        continue;
      }
      $sql_queue_update = '';
      $feature_data     = gmod_dbsf_get_feature_data($fid);
      $org_data         = gmod_dbsf_get_add_organism(array('organism_id' => $feature_data['organism_id']));
      $ncbi_taxid       = $org_data['ncbi_taxid'];
      if (empty($feature_data)) {
        drupal_set_message(t('Cannot find the feature data for this feature. Unexpected error, please contact the administrator (%fid).', array('%fid' => $fid)), 'error');
        continue;
      }
      unset($feature_data['feature_id']);
      unset($feature_data['type_id']);
      $feature_data['pgdb'] = 'chado_edit';
      // org_id relates to drupal.
      $feature_data['organism_db'] = 'default';
      $feature_name = $feature_data['uniquename'];
      $res = db_fetch_array(db_query($sql_select, $feature_name));
      $uid = $res['annotation_id'];

      //update timestamp on drupal db
      db_query('UPDATE {gmod_dbsf_feature} set timelastmodified=now() where feature_id=%d', $fid);
      switch ($op) {
        case 'Approve':
          $message .= t('Feature %feature_name is approved.<br>', array('%feature_name' => $feature_name));
          //add to chado. delete from drupal? not yet...
          $sql_queue_update = "UPDATE {genes4all_approval_queue} set approved = TRUE where uniquename='%s'";
          // store feature and featureloc
          $feature_data['add_history'] = array(
            'actor' => $form_state['values']['curator_name'],
            'action' => 'approve',
          );
          $res = gmod_dbsf_store_feature($feature_data);
          if ($res === FALSE) {
            return;
          }

          $files = file_scan_directory($fullpath, $uid .'.+\.output\.gff', array('.', '..', 'CVS'), 0, FALSE, 'filename');
          if (!empty($files)) {
            foreach ($files as $fname => $obj) {
              $gff_files[$ncbi_taxid][$feature_name][] = $fname;
            }
          }
          //gmod_dbsf_get_change_feature($feature_data['type_id'], $feature_data, 'DELETE');
          break;

        case 'Reject':
          $message .= t('Feature %feature_name is rejected.<br>', array('%feature_name' => $feature_name));
          $sql_queue_update = "UPDATE {genes4all_approval_queue} set approved = FALSE where uniquename='%s'";
          $update_cvterm_id = gmod_dbsf_get_add_cv_withnames('curator_actions', 'rejected');
          gmod_dbsf_populate_prop('feature', $fid, $update_cvterm_id, $form_state['values']['curator_name']);
          // keep it in history #gmod_dbsf_get_change_feature($feature_data['type_id'], $feature_data, 'DELETE');
          break;
      }
      db_query($sql_queue_update, $feature_name);
      db_query("UPDATE {genes4all_approval_queue} set date_curated=now() where uniquename='%s'", $feature_name);
    }
  }
  if ($op == 'Approve' && !empty($gff_files)) {
    genes4all_curate_load_approved($gff_files);
  }

  drupal_set_message($message, 'warning');
}

/**
 * Get all features that yet to be curated
 *
 * It queries the default database (Drupal) for any ORF, polypeptide or other features
 *
 * @return unknown_type
 */
function genes4all_curate_approve_features_page($features = NULL) {
  $return_text = '';
  genes4all_curate_check4annotated_features();
  $data['ORF']         = gmod_dbsf_get_change_feature('ORF');
  $data['mRNA']        = gmod_dbsf_get_change_feature('mRNA');
  $data['polypeptide'] = gmod_dbsf_get_change_feature('polypeptide');
  $data['CDS']         = gmod_dbsf_get_change_feature('CDS');
  $return_text .= drupal_get_form('genes4all_curate_approve_features_form', $data);
  return $return_text;
}

function genes4all_curate_page_feature_annotations_form($form,$table_array){
  if (empty($table_array['data'])) {
    drupal_set_message('This feature has no annotations stored currently.','warning');
    return ;
  }
  $row_check_data = array();
  $form = array();
  foreach ($table_array['data'] as $row => $value_array) {
    $cvterm_name = $value_array[0]['data'];
    $cvterm_id = $value_array[0]['id'];
    $row_check_data[$cvterm_id] = $cvterm_name;
  }
  $form             = array();
  $form['data']     = array(
  '#value' => $table_array
  );
  $form['checkboxes'] = array(
    '#type' => 'checkboxes',
    '#options' => $row_check_data,
  );
  $form['feature_id'] = array('#type' => 'value','#value' => $table_array['feature_id']);
  $form['buttons'] = array(
      'disagree' => array(
        '#type' => 'submit',
        '#value' => t('Disagree'),
        '#weight' => 15,
  ),'agree' => array(
        '#type' => 'submit',
        '#value' => t('Agree'),
        '#weight' => 10,
  ),
  );
  return $form;
}

function genes4all_curate_page_feature_annotations_form_submit($form,&$form_state){
  $cvterm_id_array = array();
  foreach ($form_state['values']['checkboxes'] as $cvterm_id){
    if (empty($cvterm_id)){continue;}
    $cvterm_id_array[]=check_plain($cvterm_id);
  }
  if (empty($cvterm_id_array)){
    drupal_set_message(t('No action was taken as no valid terms were found.'),'error');
    return;
  }
  switch($form_state['values']['op']){
    case 'Agree':
      genes4all_curation_store_user_curation_action('feature',$form_state['values']['feature_id'],$cvterm_id_array,FALSE);
      break;
    case 'Disagree':
      genes4all_curation_store_user_curation_action('feature',$form_state['values']['feature_id'],$cvterm_id_array,TRUE);
      break;
  }
  drupal_set_message(t('Request stored.'),'warning');
}

function genes4all_curate_page_feature_annotations_table($feature_name){
  if (empty($feature_name)) {return;}
  $feature_id = gmod_dbsf_chado_name2id($feature_name, 'feature');
  $local_annotations =genes4all_curation_get_all_curation_actions('feature',$feature_id);
  $table_rows       = array();
  $table_header     = array('Ontology', 'Term','1st submitter (date)','Agreements','Disagreements');
  if (!empty($local_annotations['cvterm'])) {
    foreach ($local_annotations['cvterm'] as $type => $data_array) {
      foreach ($data_array as $id=>$term) {
        $text =$type;
        $term_text = $term;
        if ($type ==='free_curation_terms'){$text = 'Free';}
        if (!empty($local_annotations['cvterm_def'][$type][$id][1])){
          $term_text = '<span>'.$term.'</span><span class="hovertip">'
          .$local_annotations['cvterm_def'][$type][$id][1].'</span>';
        }if (!empty($local_annotations['cvterm_def'][$type][$id][0])){
          $text= '<span>'.$type.'</span><span class="hovertip">'
          .$local_annotations['cvterm_def'][$type][$id][0].'</span>';
        }
        // get curator data
        $curator_data = genes4all_curation_get_user_curation_action('feature',$feature_id,$type,$term);
        $agreements = !empty($curator_data['agreements']) ? $curator_data['agreements'] :'0';
        $disagreements = !empty($curator_data['disagreements']) ? $curator_data['disagreements'] :'0';
        $first_submitter = !empty($curator_data['submitter_uid'][0])
        ? l($curator_data['submitter'][0],'user/'.$curator_data['submitter_uid'][0],
        array('attributes' => array('target' => '_blank')))
        : $curator_data['submitter'][0];
        // add data
        $table_rows[] = array(
        array('data' => $text,'id'=>$id),
        array('data' => $term_text),
        array('data' => $first_submitter.' ('.$curator_data['date'][0].')'),
        array('data' => $agreements),
        array('data' => $disagreements)
        );
      }
    }
  }
  if (!empty($local_annotations['dbxref'])) {
    foreach ($local_annotations['dbxref'] as $type => $data_array) {
      foreach ($data_array as $id=>$term) {
        $text =$type;
        $term_text = $term;
        if (!empty($local_annotations['dbxref_def'][$type][$id][1])){
          $term_text = '<span>'.$term.'</span><span class="hovertip">'
          .$local_annotations['dbxref_def'][$type][$id][1].'</span>';
        }if (!empty($local_annotations['dbxref_def'][$type][$id][0])){
          $text= '<span>'.$type.'</span><span class="hovertip">'
          .$local_annotations['dbxref_def'][$type][$id][0].'</span>';
        }
        // get curator data
        $curator_data = genes4all_curation_get_user_curation_action('feature',$feature_id,$type,$term);
        $agreements = !empty($curator_data['agreements']) ? $curator_data['agreements'] :'0';
        $disagreements = !empty($curator_data['disagreements']) ? $curator_data['disagreements'] :'0';
        $first_submitter = !empty($curator_data['submitter_uid'][0])
        ? l($curator_data['submitter'][0],'user/'.$curator_data['submitter_uid'][0],
        array('attributes' => array('target' => '_blank')))
        : $curator_data['submitter'][0];
        // add data
        $table_rows[] = array(
        array('data' => $text,'id'=>$id),
        array('data' => $term_text),
        array('data' => $first_submitter.' ('.$curator_data['date'][0].')'),
        array('data' => $agreements),
        array('data' => $disagreements)
        );
      }
    }
  }
  return array(
  'feature_id'=>$feature_id,
  'data'=>$table_rows,
  'attributes'=>array('class' => 'sortable', 'width' => '750px')
  ,'caption'=>t('Overview of existing annotations')
  ,'header'=>$table_header);
}

/**
 * Curate feature
 *
 * @param $feature_id
 *
 * @return
 * HTML
 */
function genes4all_curate_page_feature_curate($feature_name = NULL) {
  $return = "This page allows you to assign annotations to ";
  if (empty($feature_name)) {
    drupal_set_message(t('You did not specify an identifier'), 'error');
    return $return.' a feature';
  }
  $access = gmod_dbsf_chadoapi_check_restriction('select', 'feature', $feature_name);
  if ($access !== TRUE) {
    return $return. t('You are not allowed to access this resource.');
  }
  $return .=' feature '.l($feature_name,'genes4all/info/feature',array('query'=>array('feature_id'=>$feature_name)));
  $existing_annotations = drupal_get_form('genes4all_curate_page_feature_annotations_form', genes4all_curate_page_feature_annotations_table($feature_name));
  $basic_onto_form = drupal_get_form('genes4all_curate_curate_basic_ontologies_feature_form',$feature_name);
  $free_onto_form = drupal_get_form('genes4all_curate_curate_free_ontologies_feature_form',$feature_name);
  $user_onto_form = drupal_get_form('genes4all_curate_curate_user_ontologies_feature_form',$feature_name);
  $user_dbxref_form = drupal_get_form('genes4all_curate_curate_user_dbxref_feature_form',$feature_name);
  $control_form  = drupal_get_form('genes4all_curate_curate_control_feature_form', $feature_name);
  $page_tabs = array(
    'curation' => array(
      '#type' => 'tabset',
      'existing' => array(
        '#type' => 'tabpage',
        '#title' => 'Existing',
        '#content' => $existing_annotations,
        '#weight' => -1,
  ),'basic_ontologies' => array(
        '#type' => 'tabpage',
        '#title' => 'Basic',
        '#content' => $basic_onto_form,
        '#weight' => 1,
  ),'free_tags' => array(
        '#type' => 'tabpage',
        '#title' => 'Free tags',
        '#content' => $free_onto_form,
        '#weight' => 3,
  ),'control' => array(
        '#type' => 'tabpage',
        '#title' => 'Control',
        '#content' => $control_form,
        '#weight' => 5,
  ),
  ),
  );
  if (strlen($user_onto_form) > 700) {
    $page_tabs['curation']['user_ontologies'] = array(
        '#type' => 'tabpage',
        '#title' => 'Users\' CV' ,
        '#content' => $user_onto_form,
        '#weight' => 2,
    );
  }
  if (strlen($user_dbxref_form) > 700) {
    $page_tabs['curation']['user_dbxref'] = array(
        '#type' => 'tabpage',
        '#title' => 'Cross-refs',
        '#content' => $user_dbxref_form,
        '#weight' => 4,
    );
  }
  $return .= tabs_render($page_tabs);
  return $return;
}

/**
 * Get (drupal) feature details for approval
 *
 * Only senior curators can see this page.
 *
 * @param $feature_name
 *
 * @return unknown_type
 */
function genes4all_curate_view_feature_page($feature_name) {
  if (empty($feature_name)) {
    $feature_name = $_GET['feature'];
    if (empty($feature_name)) {
      drupal_set_message(t('You did not specify an identifier'), 'error');
      return '';
    }
  }
  $access = gmod_dbsf_chadoapi_check_restriction('select', 'feature', $feature_name);
  if ($access != TRUE) {
    return "You are not allowed to access this resource.";
  }
  genes4all_curate_check4annotated_features();
  $dirpath = file_create_path() .'/genes4all_curate/';

  $return_text = "You are looking at feature $feature_name.";

  $sql_select = "SELECT annotation_id,approved,annotated from {genes4all_approval_queue} where uniquename='%s'";
  $res        = db_fetch_array(db_query($sql_select, $feature_name));
  $uid        = $res['annotation_id'];
  $appr       = $res['approved'];
  if ($appr === FALSE) {
    $return_text .= ' It has been rejected.';
  }
  elseif ($appr === TRUE) {
    $return_text .= ' It has been approved.';
  }
  else {
    $return_text .= ' It has not been considered for approval.';
  }
  $anno = $res['annotated'];
  /*
   * annotated=0: not finished, not checked
   * annotated=1 finished
   * annotated=2 failed/non-existing uid.
   */


  if ($anno == 1) {
    $return_text .= ' Annotation runs for this feature have been completed.';
  }
  elseif ($anno == 0) {
    $return_text .= ' Annotation runs for this feature have yet to finish.';
  }
  elseif ($anno == 2) {
    $return_text .= ' Annotation runs for this feature encountered an error.';
  }

  $table_rows = array();
  $gff_files = file_scan_directory(base_path() .'/'. $dirpath, $uid .'.+\.gff', array('.', '..', 'CVS'), 0, FALSE, 'filename');
  if (!empty($gff_files)) {
    foreach ($gff_files as $fname => $obj) {
      $base = $obj->basename;
      $data = gmod_dbsf_parse_gfffile($fname);
      if (empty($data)) {
        $table_rows[] = array($base, 'No results');
        continue;
      }
      foreach ($data[$feature_name] as $software => $type_array) {
        foreach ($type_array as $type => $count) {
          if (strpos($type, 'part') === FALSE) {
            $table_rows[] = array($base, $software, $type, $count);
          }
        }
      }
    }
  }
  $table_header  = array('File', 'Software', 'Type', 'Count');
  $table_caption = 'Overview of automatic annotations';
  $existing_annotations = theme_table($table_header, $table_rows, array('class' => 'genes4all_explore-library-table'), $table_caption);
  return $return_text . $existing_annotations;
}

/**
 * Add new feature
 *
 * @param $feature_id
 *
 * @return
 * HTML
 */
function genes4all_curate_page_feature_add() {
  $return = "This page allows you to create a new sequence feature (e.g. a gene) and submit it to the automatic annotation queue. You will then be able to curate it via ". l('the curation link', 'genes4all/feature/curate');

  $add_form = drupal_get_form('genes4all_curate_add_feature_form');
  if (strlen($add_form) < 400) {
    $add_form = '';
  }
  $page_tabs = array(
    'curation' => array(
      '#type' => 'tabset',
      'ontologies' => array(
        '#type' => 'tabpage',
        '#title' => 'Basic properties',
        '#content' => $add_form,
        '#weight' => 0,
  ),
  ),
  );
  $return .= tabs_render($page_tabs);
  return $return;
}

/**
 * Add a new feature. Implements hook_form()
 *
 * @param $form_state
 *
 * @return
 * HTML
 */
function genes4all_curate_add_feature_form($form_state) {
  $form = array();
  $database_array     = gmod_dbsf_get_add_var('genes4all_my_central_federation');
  if (empty($database_array)){
    drupal_set_message(t('The module genes4all_curate has not been setup: please ask the system administrator to setup the federation using the module\'s ').l('administrative setting page','admin/settings/genes4all_curates').' .','error');
    return;
  }
  $org_array          = gmod_dbsf_db_execute('chado', 'gmod_dbsf_getspecies', NULL, TRUE);
  $type_array         = array('mRNA' => 'mRNA', 'ORF' => 'ORF', 'CDS' => 'partial CDS', 'polypeptide' => 'protein');
  $genetic_code_array = array('Standard', 'Vertebrate Mitochondrial', 'Yeast Mitochondrial', 'Mold, Protozoan, and Coelocoel Mitochondrial'
  , 'Invertebrate Mitochondrial', 'Ciliate Nuclear', 'Echinoderm Mitochondrial', 'Euplotid Nuclear'
  , 'Bacterial', 'Alternative Yeast Nuclear', 'Ascidian Mitochondrial', 'Flatworm Mitochondrial'
  , 'Blepharisma Macronuclear',
  );
  /*
   $type_array=gmod_dbsf_db_execute('chado','gmod_dbsf_get_add_cv_withnames','sequence');
   foreach ($type_array as $id=>$name){
   $name=str_replace('_',' ',$name);
   $type_array[$id]=$name;
   if (strlen($name)>85){
   $type_array[$id]=substr($name,0,80).' ...';
   }
   }
   */
  $species_options = $org_array['phylogeny_array'];
  if (empty($species_options)){
    $species_options =$org_array['species_array'];
    if (empty($species_options)){
      drupal_set_message(t('Sorry, this function is not available as the organism table has not been properly populated.'),'error');
      return '';
    }
  }

  $form['data'] = array(
    '#title' => 'Basic properties for new sequence-based feature',
    '#description' => t('Chado defines a feature to be a region of a biological polymer (typically a DNA, RNA, or a polypeptide molecule) or an aggregate of regions on this polymer. ' . 'Use this form to ask that your sequence is included in the database and a "feature page" (i.e. gene page) is constructed for it.'
    ),
    '#type' => 'fieldset',
    'database' => array(
      '#type' => 'select',
      '#size' => 1,
      '#title' => 'Database',
      '#description' => t('Which database references your new sequence-based feature. This is used in federated systems.'),
      '#options' => $database_array,
      '#weight' => -1,
      '#required' => TRUE,
    ), 'organism' => array(
      '#type' => 'select',
      '#size' => 1,
      '#title' => 'Species',
      '#description' => t('To select faster, the select box will respond to key strokes: try typing the name. If your species does not appear on the list, then you need to '. l('add', 'genes4all/organism/add') .' it first.'),
      '#options' => $species_array,
      '#weight' => 0,
      '#required' => TRUE,
    ), 'name' => array(
      '#type' => 'textfield',
      '#title' => 'Friendly name/alias',
      '#description' => t('You can optionally include a friendly name. This will be shown on the feature page (a.k.a. "gene page").'),
      '#maxlength' => 200,
      '#weight' => 2,
    ), 'residues' => array(
      '#type' => 'textarea',
      '#title' => 'Sequence',
      '#description' => t('Please provide the sequence you wish to store for this feature in plain text format. ' . 'It is case insensitive but please provide it in DNA or protein alphabet. ' . 'Degeneracy is allowed but take care not to include any non-IUPAC characters.' . ' If you are storing a <i>mRNA</i> or an <i>ORF</i> then include both the start and stop codons ' . '(if they are unavailable, then provide a <i>CDS</i>). If submitting a CDS, ensure you are submitting in the ' . 'correct frame (first base will be the 1st position of a codon, last base will be the 3rd position of a codon).' . ' Further, partial proteins are allowed but providing a partial CDS is prefered.'
      ),
      '#weight' => 4,
      '#required' => TRUE,
      ), 'type_name' => array(
      '#type' => 'select',
      '#size' => 1,
      '#title' => 'Sequence type',
      '#options' => $type_array,
      '#weight' => 3,
      '#required' => TRUE,
      '#description' => t('Please select whether your submission is a mRNA (transcribed sequence, can include UTR but no introns), a full Open Reading Frame (ORF), a partial CoDing Sequence (CDS) or a protein. It is prefered to provide the ORF rather than a protein as the latter can be acquired via a translation. We are currently not utilizing the UTR but may do so in the future.'),
      ), 'genetic_code' => array(
      '#type' => 'select',
      '#size' => 1,
      '#title' => 'Genetic code',
      '#description' => t('What genetic code is being used for this coding sequence?'),
      '#options' => $genetic_code_array,
      '#required' => TRUE,
      '#weight' => 5,
      ), 'gene_start' => array(
      '#type' => 'textfield',
      '#size' => 5,
      '#title' => 'First base of start codon',
      '#description' => t('If you selected <strong>mRNA</strong> as your sequence type, provide the position where translation starts (numbering starts from 1). With <i>ORF</i> it is assumed that the first three letters are the start codon.'),
      '#weight' => 6,
      ), 'gene_stop' => array(
      '#type' => 'textfield',
      '#size' => 5,
      '#title' => 'Last base of stop codon',
      '#description' => t('If you selected <strong>mRNA</strong> as your sequence type, provide the position where translation stops (numbering starts from 1 but you can also give a negative number to indicate how many residues from the end: e.g. -3 is three residues from the end (thus the first position of the stop codon must be =< -3 or >=1). With <i>ORF</i> it is assumed that the last three letters are the stop codon.'),
      '#weight' => 7,
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

/**
 * Implements hook_form_validate
 *
 * @param unknown_type $form
 * @param unknown_type $form_state
 *
 * @return unknown_type
 */
function genes4all_curate_add_feature_form_validate($form, &$form_state) {
  $organism_id = $form_state['values']['organism'];
  // new: check that submitted species is indeed allowed.
  $org_data = gmod_dbsf_db_execute('chado', 'gmod_dbsf_getspecies', NULL, 'all');
  if (empty($org_data['species_array'][$organism_id])) {
    // attempt to DOM hack: stop dead in your tracks
    form_set_error('organism', t('You requested a species which you are not allowed access to.'));
    return FALSE;
  }
  unset($org_data);
  // Synchronise organisms between drupal and Chado
  $org_data   = gmod_dbsf_db_execute('chado', 'gmod_dbsf_get_add_organism', array('organism_id' => $organism_id, 'pgtype' => 'chado'));
  $drupal_org = gmod_dbsf_db_execute('default', 'gmod_dbsf_get_add_organism', $org_data + array('pgtype' => 'drupal'));
  $drupal_org = $drupal_org['organism_id'];
  //TODO add to DB
  $available_start_codons = array('ATG');
  $available_stop_codons = array('TAG', 'TGA', 'TAA');

  // set genetic code to start from 0
  $form_state['values']['genetic_code']++;
  $form_state['values']['residues'] = strtoupper($form_state['values']['residues']);
  $check_seq = '';
  $type = 'DNA';
  $start_codon = '';
  $stop_codon = '';
  $prefix = '';

  if ($form_state['values']['type_name'] == 'polypeptide') {
    $check_seq = gmod_dbsf_validate_seq_protein($form_state['values']['residues'], TRUE);
    $type = 'protein';
  }
  else {
    $check_seq = gmod_dbsf_validate_seq_dna($form_state['values']['residues'], TRUE);
  }
  if ($check_seq == FALSE || is_numeric($check_seq)) {
    form_set_error('residues', t('Your sequence does not seem to be a valid %type sequence (%check_seq errors).', array('%type' => $type, '%check_seq' => $check_seq)));
  }
  else {
    if (substr_count($check_seq, '*') > 1) {
      form_set_error('residues', t('Your sequence contains more than one stop codon! You may only submit coding sequence, perhaps you are not using the correct genetic code?'));
    }
    $form_state['values']['residues'] = $check_seq;
  }


  $form_state['values']['alphabet'] = $type;
  $dna2orf = '';
  switch ($form_state['values']['type_name']) {
    case 'mRNA':
      if (empty($form_state['values']['gene_start'])) {
        form_set_error('gene_start', t('A start codon position is needed for mRNA sequences.'));
      }
      if (empty($form_state['values']['gene_stop'])) {
        form_set_error('gene_stop', t('A stop codon position is needed for mRNA sequences.'));
        return;
      }
      $gene_length = strlen($form_state['values']['residues']);
      $orf_length = $gene_length - ($form_state['values']['gene_start'] - 1);

      $gene_stop = $form_state['values']['gene_stop'];
      // user numbering starts from 1 or is neg.
      if ($gene_stop > 0) {
        // user gives first position of codon
        $orf_length -= ($gene_length - $gene_stop) + 2;
        // for $stop_codon
        $gene_stop--;
      }
      else {
        $orf_length += ($gene_stop + 1) + 2;
      }
      $dna2orf = substr($form_state['values']['residues'], $form_state['values']['gene_start'] - 1, $orf_length);
      $form_state['values']['dna2orf'] = $dna2orf;
      if ($orf_length % 3 !== 0) {
        form_set_error('residues', t('You gave a coding sequence which has partial codons!'));
      }
      $start_codon = substr($form_state['values']['residues'], $form_state['values']['gene_start'] - 1, 3);
      $stop_codon = substr($form_state['values']['residues'], $gene_stop, 3);

      if (empty($start_codon) || ($start_codon && !in_array($start_codon, $available_start_codons))) {
        form_set_error('residues', t('There does not seem to be a start codon where you specified in your sequence (%start_codon).', array('%start_codon' => $start_codon)));
      }
      if (empty($stop_codon) || ($stop_codon && !in_array($stop_codon, $available_stop_codons))) {
        form_set_error('residues', t('There does not seem to be a stop codon where you specified in your sequence (%stop_codon).', array('%stop_codon' => $stop_codon)));
      }
      $prefix = 'Cgen';
      break;

    case 'ORF':
      if (strlen($form_state['values']['residues']) % 3 !== 0) {
        form_set_error('residues', t('You gave a coding sequence which has partial codons!'));
      }
      $start_codon = substr($form_state['values']['residues'], 0, 3);
      $stop_codon = substr($form_state['values']['residues'], -3, 3);
      if (empty($start_codon) || ($start_codon && !in_array($start_codon, $available_start_codons))) {
        form_set_error('residues', t('There does not seem to be a start codon at the beginning of your sequence (%start_codon).', array('%start_codon' => $start_codon)));
      }
      if (empty($stop_codon) || ($stop_codon && !in_array($stop_codon, $available_stop_codons))) {
        form_set_error('residues', t('There does not seem to be a stop codon at the end of your sequence ($stop_codon).', array('%stop_codon' => $stop_codon)));
      }
      $dna2orf = $form_state['values']['residues'];
      $form_state['values']['dna2orf'] = $dna2orf;
      $prefix = 'Corf';
      break;

    case 'CDS':
      if (strlen($form_state['values']['residues']) % 3 !== 0) {
        form_set_error('residues', t('You gave a coding sequence which has partial codons!'));
      }
      $dna2orf = $form_state['values']['residues'];
      $form_state['values']['dna2orf'] = $dna2orf;

      $prefix = 'Ccds';
      break;

    case 'polypeptide':
      $prefix = 'Cpep';
      break;
  }

  // translate if needed
  $protein_seq = '';
  if (!empty($dna2orf)) {
    $protein_seq = gmod_dbsf_translate_DNA_to_protein($dna2orf, $form_state['values']['genetic_code']);
    if (substr_count($protein_seq, '*') > 1) {
      form_set_error('residues', t('Your sequence contains more than one stop codon! You may only submit coding sequence, perhaps you are not using the correct genetic code?'));
    }
    $form_state['values']['protein']['seq'] = $protein_seq;
    // decide which is latest: drupal or chado
    $array        = array();
    $array_drupal = gmod_dbsf_db_execute('default', 'gmod_dbsf_chado_latest_assembly_serial', array('drupal' => TRUE, 'org' => $drupal_org, 'prefix' => 'Cpep', 'db' => $form_state['values']['database']));
    $array        = gmod_dbsf_db_execute('chado', 'gmod_dbsf_chado_latest_assembly_serial', array('org' => $organism_id, 'prefix' => 'Cpep', 'db' => $form_state['values']['database']));
    if (empty($array['latest_assembly']) && !empty($array_drupal['latest_assembly'])) {
      $array = $array_drupal;
    }
    elseif (!empty($array_drupal['latest_assembly']) && !empty($array['latest_assembly'])) {
      $array['latest_assembly'] = (strcmp($array_drupal['latest_assembly'], $array['latest_assembly']) < 0) ? $array['latest_assembly'] : $array_drupal['latest_assembly'];
      $array['next_serial'] = (strcmp($array_drupal['next_serial'], $array['next_serial']) < 0) ? $array['next_serial'] : $array_drupal['next_serial'];
    }
    if (empty($array['latest_assembly'])) {
      $form_state['values']['protein']['latest_assembly'] = 'Aa';
      $form_state['values']['protein']['next_serial'] = 1;
    }
    else {
      $form_state['values']['protein']['latest_assembly'] = $array['latest_assembly'];
      $form_state['values']['protein']['next_serial'] = $array['next_serial'];
    }
    // if mRNA then extra ORF
    if ($form_state['values']['type_name'] == 'mRNA') {
      $array        = array();
      $array_drupal = gmod_dbsf_db_execute('default', 'gmod_dbsf_chado_latest_assembly_serial', array('drupal' => TRUE, 'org' => $drupal_org, 'prefix' => 'Corf', 'db' => $form_state['values']['database']));
      $array        = gmod_dbsf_db_execute('chado', 'gmod_dbsf_chado_latest_assembly_serial', array('org' => $organism_id, 'prefix' => 'Corf', 'db' => $form_state['values']['database']));
      if (empty($array['latest_assembly']) && !empty($array_drupal['latest_assembly'])) {
        $array = $array_drupal;
      }
      elseif (!empty($array_drupal['latest_assembly']) && !empty($array['latest_assembly'])) {
        $array['latest_assembly'] = (strcmp($array_drupal['latest_assembly'], $array['latest_assembly']) < 0) ? $array['latest_assembly'] : $array_drupal['latest_assembly'];
        $array['next_serial'] = (strcmp($array_drupal['next_serial'], $array['next_serial']) < 0) ? $array['next_serial'] : $array_drupal['next_serial'];
      }

      if (empty($array['latest_assembly'])) {
        $form_state['values']['extra_orf']['latest_assembly'] = 'Aa';
        $form_state['values']['extra_orf']['next_serial'] = 1;
      }
      else {
        $form_state['values']['extra_orf']['latest_assembly'] = $array['latest_assembly'];
        $form_state['values']['extra_orf']['next_serial'] = $array['next_serial'];
      }
    }
  }

  $form_state['values']['prefix'] = $prefix;
  // get org ncbi tax id
  $form_state['values']['ncbi_taxid'] = gmod_dbsf_db_execute('chado', 'gmod_dbsf_chado_get_ncbi_taxid', $organism_id);
  if (empty($form_state['values']['ncbi_taxid'])) {
    form_set_error('organism', t('Cannot find the NCBI taxonomy ID in your Chado database. This is a fault the administrator has to address.'));
    return;
  }
  //get latest assembly version for this organism
  $array        = array();
  $array_drupal = gmod_dbsf_db_execute('default', 'gmod_dbsf_chado_latest_assembly_serial', array('drupal' => TRUE, 'org' => $drupal_org, 'prefix' => $prefix, 'db' => $form_state['values']['database']));
  $array        = gmod_dbsf_db_execute('chado', 'gmod_dbsf_chado_latest_assembly_serial', array('org' => $organism_id, 'prefix' => $prefix, 'db' => $form_state['values']['database']));
  if (empty($array['latest_assembly']) && !empty($array_drupal['latest_assembly'])) {
    $array = $array_drupal;
  }
  elseif (!empty($array_drupal['latest_assembly']) && !empty($array['latest_assembly'])) {
    $array['latest_assembly'] = (strcmp($array_drupal['latest_assembly'], $array['latest_assembly']) < 0) ? $array['latest_assembly'] : $array_drupal['latest_assembly'];
    $array['next_serial'] = (strcmp($array_drupal['next_serial'], $array['next_serial']) < 0) ? $array['next_serial'] : $array_drupal['next_serial'];
  }

  if (empty($array['latest_assembly'])) {
    drupal_set_message(t('Found no stored sequences for this organism, setting latest assembly and serial to 1.'), 'warning');
    $form_state['values']['latest_assembly'] = 'Aa';
    $form_state['values']['next_serial'] = 1;
  }
  else {
    $form_state['values']['latest_assembly'] = $array['latest_assembly'];
    $form_state['values']['next_serial'] = $array['next_serial'];
  }
}

/**
 * Implements hook_form_submit()
 *
 * @param $form
 * @param $form_state
 *
 * @return unknown_type
 */
function genes4all_curate_add_feature_form_submit($form, &$form_state) {
  $db_id              = $form_state['values']['database'];
  $protein_data_array = $form_state['values']['protein'];
  $orf_data_array     = $form_state['values']['extra_orf'];
  $alphabet           = $form_state['values']['alphabet'];
  $name               = $form_state['values']['name'];
  $organism_id        = $form_state['values']['organism'];
  $genetic_code       = $form_state['values']['genetic_code'];
  $sequence           = $form_state['values']['residues'];
  $orf_seq            = $form_state['values']['dna2orf'];
  $orf_length         = strlen($orf_seq);
  $gene_start         = $form_state['values']['gene_start'];
  $gene_stop          = $form_state['values']['gene_stop'];
  if (!empty($gene_stop) && $gene_stop < 0) {
    $gene_stop = strlen($form_state['values']['residues']) + $gene_stop;
  }
  $data_type = $form_state['values']['type_name'];
  $ncbi_taxid = $form_state['values']['ncbi_taxid'];
  // e.g. Cgen Corf Cpep etc
  $prefix          = $form_state['values']['prefix'];;
  $latest_assembly = $form_state['values']['latest_assembly'];
  $next_serial     = $form_state['values']['next_serial'];
  // whole serial, including org, assembly, prefix and suffix
  $uniquename         = '';
  $orf_uniquename     = '';
  $pep_uniquename     = '';
  $submitted_features = array();
  if (empty($db_id) || empty($ncbi_taxid) || empty($latest_assembly) || empty($prefix) || empty($next_serial)) {
    form_set_error('all', t('Sorry but I do not have enough data to create new identifier. Please let the administrator know what you tried to do.'));
    return FALSE;
  }

  $dirpath = file_create_path() .'/genes4all_curate';
  $submission_uid = gmod_dbsf_create_uid(NULL, NULL, 'genes4all_curate');
  //-seqtype $seqtype -trtable $codon_table

  // prepare operations
  $batch_file_data = array();
  $operations      = array();
  $save_array      = array();

  //now we want to add it to the drupal database so that the admin can approve it.
  //Once approved, it will be created in the chado database.
  // the feature table will be used, remember to ensure that the time field is transferred.
  $uniquename = $db_id . $ncbi_taxid . $latest_assembly . $prefix . $next_serial;
  $submitted_features[] = $uniquename;
  //check if feature exists in drupal (but not chado). if it does, set it to last available one.

  if (empty($name)) {
    $name = $uniquename;
  }
  $uid = gmod_dbsf_create_uid(NULL, NULL, $uniquename);
  $feature_id = gmod_dbsf_db_execute('default', 'gmod_dbsf_store_feature', array(
      'pgdb' => 'default',
      'residues' => $sequence,
      'uniquename' => $uniquename,
      'name' => $name,
      'organism_id' => $organism_id,
      'organism_db' => 'chado',
      'type' => $data_type,
  ));
  if (empty($feature_id)) {
    drupal_set_message(t('There has been an error synchronizing the feature table. Unexpected: please report the bug.'), 'error');
    return;
  }
  $success = genes4all_curate_add_approval_queue($uniquename, $uid, $submission_uid);
  if (empty($success)) {
    drupal_set_message(t('Submission halted: unknown user.'), 'error');
    return;
  }

  // nucleotide sequence
  $batch_file_data['name'][$uid] = $uniquename;
  $batch_file_data['infile'][$uid] = $sequence;
  $batch_file_data['outfile'][$uid] = $dirpath .'/'. $uid .'.query';
  $batch_file_data['filetype'][$uid] = 'nucleotide';
  $batch_file_data['format'][$uid] = FALSE;
  $save_array[$submission_uid][$uid] = array('seqtype' => 'n', 'gencode' => $genetic_code);
  // if mRNA then extra ORF
  if ($prefix == 'Cgen') {
    $orf_uniquename = $db_id . $ncbi_taxid . $orf_data_array['latest_assembly'] .'Corf'. $orf_data_array['next_serial'];
    $submitted_features[] = $orf_uniquename;
    $orf_uid = gmod_dbsf_create_uid(NULL, NULL, $orf_uniquename);
    $orf_feature_id = gmod_dbsf_db_execute('default', 'gmod_dbsf_store_feature', array(
        'pgdb' => 'default',
        'residues' => $orf_seq,
        'uniquename' => $orf_uniquename,
        'name' => 'ORF in '. $name,
        'organism_id' => $organism_id,
        'organism_db' => 'chado',
        'type' => 'ORF',
        'parent' => array('parent_id' => $feature_id, 'start' => $gene_start, 'end' => $gene_stop, 'phase' => 0, 'strand' => 1),
    )
    );
    $orf_placeholder_feature_id = gmod_dbsf_db_execute('default', 'gmod_dbsf_store_feature', array(
        'pgdb' => 'default',
        'residues' => $orf_seq,
        'uniquename' => $orf_uniquename .'_0',
        'name' => 'Placeholder for ORF in '. $name,
        'organism_id' => $organism_id,
        'organism_db' => 'chado',
        'type' => 'ORF',
        'parent' => array('parent_id' => $feature_id, 'start' => $gene_start, 'end' => $gene_stop, 'phase' => 0, 'strand' => 1),
    )
    );
    if (empty($orf_feature_id)) {
      drupal_set_message(t('There has been an error synchronizing the feature table. Unexpected: please report the bug.'), 'error');
      return;
    }
    $success = genes4all_curate_add_approval_queue($orf_uniquename, $orf_uid, $submission_uid);
    if (empty($success)) {
      drupal_set_message(t('Submission halted: unknown user.'), 'error');
      return;
    }
    $batch_file_data['name'][$orf_uid] = $orf_uniquename;
    $batch_file_data['infile'][$orf_uid] = $orf_seq;
    $batch_file_data['outfile'][$orf_uid] = $dirpath .'/'. $orf_uid .'.query';
    $batch_file_data['filetype'][$orf_uid] = 'nucleotide';
    $batch_file_data['format'][$orf_uid] = FALSE;
    $save_array[$submission_uid][$orf_uid] = array('seqtype' => 'n', 'gencode' => $genetic_code);
  }
  else {
    //ensure polypeptide is linked directly to ORF and not mRNA.
    $orf_length = strlen($sequence);
    $orf_feature_id = $feature_id;
  }
  // repeat for translated protein
  $pep_feature_id = 0;
  if (!empty($protein_data_array)) {
    $pep_uniquename = $db_id . $ncbi_taxid . $protein_data_array['latest_assembly'] .'Cpep'. $protein_data_array['next_serial'];
    $submitted_features[] = $pep_uniquename;
    $pep_uid = gmod_dbsf_create_uid(NULL, NULL, $pep_uniquename);

    $pep_feature_id = gmod_dbsf_db_execute('default', 'gmod_dbsf_store_feature', array(
        'pgdb' => 'default',
        'residues' => $form_state['values']['protein']['seq'],
        'uniquename' => $pep_uniquename,
        'name' => 'Translation of '. $name,
        'organism_id' => $organism_id,
        'organism_db' => 'chado',
        'type' => 'polypeptide',
        'parent' => array('parent_id' => $orf_feature_id, 'start' => 1, 'end' => $orf_length, 'phase' => 0, 'strand' => 1),
    )
    );
    $pep_placeholder_feature_id = gmod_dbsf_db_execute('default', 'gmod_dbsf_store_feature', array(
        'pgdb' => 'default',
        'residues' => $form_state['values']['protein']['seq'],
        'uniquename' => $pep_uniquename .'_0',
        'name' => 'Placeholder for translation of '. $name,
        'organism_id' => $organism_id,
        'organism_db' => 'chado',
        'type' => 'polypeptide',
        'parent' => array('parent_id' => $orf_feature_id, 'start' => 1, 'end' => $orf_length, 'phase' => 0, 'strand' => 1),
    )
    );
    if (empty($pep_feature_id)) {
      drupal_set_message(t('There has been an error synchronizing the feature table. Unexpected: please report the bug.'), 'error');
      return;
    }
    $success = genes4all_curate_add_approval_queue($pep_uniquename, $pep_uid, $submission_uid);
    if (empty($success)) {
      drupal_set_message(t('Submission halted: unknown user.'), 'error');
      return;
    }
    $batch_file_data['name'][$pep_uid] = $pep_uniquename;
    $batch_file_data['infile'][$pep_uid] = $protein_data_array['seq'];
    $batch_file_data['outfile'][$pep_uid] = $dirpath .'/'. $pep_uid .'.query';
    $batch_file_data['filetype'][$pep_uid] = 'protein';
    $batch_file_data['format'][$pep_uid] = FALSE;
    $save_array[$submission_uid][$pep_uid] = array('seqtype' => 'p');
  }

  $operations[] = array('gmod_dbsf_batch_upload_fasta', array($batch_file_data));
  $operations[] = array('gmod_dbsf_batch_save_data', array($save_array));

  $batch = array(
    'title' => t('Preparing data needed your submission...'),
    'operations' => $operations,
    'init_message' => t('Processing your submission...'),
    'progress_message' => t('@remaining operations remaining...'),
    'error_message' => t('Your submission encountered an error.'),
    'finished' => 'genes4all_curate_batch_finished',
    'file' => drupal_get_path('module', 'genes4all_curate') .'/includes/genes4all_curate_feature.inc',
  );
  batch_set($batch);
  $links = '';
  foreach ($submitted_features as $feature) {
    $links .= ' '. l($feature, 'genes4all/curate/feature/view/', array('query' => array('feature' => $feature)));
  }
  drupal_set_message(t('If successful, your submissions will be assigned these IDs (one for each sequence type): %links.'), array('%links' => $links));
  $form_state['redirect'] = array('genes4all/feature/added', "submission_uid=$submission_uid");
  return "genes4all/feature/added?submission_uid=$submission_uid";
}

/**
 * Process the annoation batch for submission to condor/local processing
 *
 * Adapted from dbsf_blast_batch_finished(), an extra depth has been added
 * to the $results array to allow for multiple software per sequence.
 * @see dbsf_blast_batch_finished()
 * @see http://api.drupal.org/api/group/batch
 *
 * @param unknown_type $success
 * @param unknown_type $results
 * @param unknown_type $operations
 *
 * @return unknown_type
 */
function genes4all_curate_batch_finished($success, $results, $operations) {
  // store uids to report in case of error.
  $uid_array = array();
  $dirpath   = file_create_path() .'/genes4all_curate';
  $message   = '';


  if (!empty($success)) {
    $run_ipr = dbsf_check_software_active('iprscan');
    //TODO : parse interpro xml file
    $ipr_par = ' -cli -nocrc -goterms -iprlookup -format raw -seqtype p';
    $blast_par = ' -m 7 ';
    // Get the blast databases the administrator has decided to use for
    // curated features
    $blastdbs = gmod_dbsf_get_linked_resources('Dataset for auto-annotation', FALSE);
    $subject_dbs = array();
    if (empty($blastdbs)) {
      drupal_set_message(t('No BLAST databases have been defined for annotation. No annotations will run'), 'error');
      return;
    }
    foreach ($blastdbs as $dataset_id => $dataset_values) {
      $dataset_name    = $dataset_values['uniquename'];
      $dataset_extra   = gmod_dbsf_get_resource_cvterm_and_prop($dataset_id);
      $dbtype          = implode(', ', $dataset_extra['cvterm']['dataset_type']);
      $linked_res_data = gmod_dbsf_get_linked_resources($dataset_id, NULL, TRUE);
      foreach ($linked_res_data as $object_name => $object_data) {
        $object_type = $object_data['type_name'];
        $to_add[$object_type] = $object_name;
      }
      $subject_dbs[$dbtype][$dataset_name] = $to_add['directory'] . $to_add['filename'];
    }

    $submission_uid      = '';
    $submission_to_store = "software:genes4all_curate\n";
    $jobs_to_do          = array();

    foreach ($results['data'] as $submission_uid => $batch_data) {
      $submission_to_store .= "submission_uid:$submission_uid\n";
      foreach ($batch_data as $uid => $seq_data) {
        $jobs_to_do[$uid] = $seq_data;
      }
    }
    if (!empty($submission_uid)) {
      foreach ($jobs_to_do as $uid => $seq_data) {
        $seq_type  = $seq_data['seqtype'];
        $gen_code  = $seq_data['gencode'];
        $par_array = array();
        switch ($seq_type) {
          case 'p':
            $more_data = array();
            //blast
            foreach ($subject_dbs as $dbtype => $db_data) {
              if (empty($db_data)) {
                continue;
              }
              $algorithm = '';
              switch ($dbtype) {
                case 'protein':
                  $algorithm = 'blastp';
                  break;

                case 'nucleotide':
                  $algorithm = 'tblastn';
                  break;
              }
              foreach ($db_data as $db_name => $db_file) {
                $submission_to_store .= 'algorithm;id='. $algorithm .";$uid\n";
                $more_data['database'][] = $db_name;
                $par_str = $blast_par ." -p $algorithm -d '$db_file'";
                if (!empty($gen_code)) {
                  $par_str .= " -Q $gen_code ";
                }
                $par_array[] = $par_str;
              }
            }
            gmod_dbsf_parameter_daemon($uid, $dirpath, $par_array, 'blastall', $more_data);
            if (!empty($run_ipr)) {
              $submission_to_store .= 'algorithm;id=' . 'iprscan'. ";$uid\n";
              gmod_dbsf_parameter_daemon($uid, $dirpath, $ipr_par, 'iprscan');
            }
            break;

          case 'n':
            $more_data = array();

            //blast
            foreach ($subject_dbs as $dbtype => $db_data) {
              if (empty($db_data)) {
                continue;
              }
              $algorithm = '';
              switch ($dbtype) {
                case 'protein':
                  $algorithm = 'blastx';
                  break;

                case 'nucleotide':
                  $algorithm = 'tblastx';
                  break;
              }
              foreach ($db_data as $db_name => $db_file) {
                $submission_to_store .= 'algorithm;id='. $algorithm .";$uid\n";
                $more_data['database'][] = $db_name;
                $par_str = $blast_par ." -p $algorithm -d '$db_file'";
                if (!empty($gen_code)) {
                  $par_str .= " -Q $gen_code ";
                }
                $par_array[] = $par_str;
              }
            }
            gmod_dbsf_parameter_daemon($uid, $dirpath, $par_array, 'blastall', $more_data);
            break;
        }
      }

      // finally store all data in a submission file for later linked retrieval.

      $outfile = $dirpath .'/'. $submission_uid .'.submission';
      if (!$outhandle = fopen($outfile, 'wt')) {
        drupal_set_message(t('Could not create %outfile.', array('%outfile' => $outfile)), 'error');
        return FALSE;
      }
      fwrite($outhandle, $submission_to_store);
      fclose($outhandle);
      $message = t('Your submission has been included in the annotation queue for processing...');
    }
    else {
      $message = t('There has been a problem submitting your annotation. Please '. l('contact', 'contact') .' the administrator.');
    }
  }
  else {
    $error_operation = reset($operations);
    $uids            = implode(' ', $uid_array);
    $message         = t('An error occurred while processing your request. Please '. l('contact', 'contact') .' the administrator reporting: %uids.', array('%uids' => $uids));
    drupal_set_message($message, 'warning');
    return FALSE;
  }
  drupal_set_message($message, 'warning');
}

function genes4all_curate_page_feature_added($submission_uid = NULL) {
  if (empty($submission_uid)) {
    $submission_uid = $_GET['submission_uid'];
    if (empty($submission_uid)) {
      return 'No submission data provided';
    }
  }
  // there are two IDs: the one relating to the feature (annotation_id)
  // and submission IDs associated with the user's submission. The user has access only to the latter.

  $sql_select = "SELECT approved,annotated,annotation_id,uniquename from {genes4all_approval_queue} where submission_id='%s'";
  $res = db_query($sql_select, $submission_uid);

  $table_rows = array();
  while ($db_res = db_fetch_array($res)) {
    $annotation_id = $db_res['annotation_id'];
    $approved      = $db_res['approved'];
    $annotated     = $db_res['annotated'];
    $uname         = $db_res['uniquename'];
    if ($approved === TRUE) {
      $approved = 'Approved';
    }
    elseif ($approved === FALSE) {
      $approved = 'Rejected';
    }
    else {
      $approved = 'Not considered';
    }
    if (empty($annotated)) {
      $annotated = 'In queue';
    }
    elseif ($annotated == 1) {
      $annotated = 'Annotated';
    }
    else {
      $annotated = 'Encountered an error';
    }
    $table_rows[] = array($uname, $annotated, $approved);
  }

  if (empty($table_rows)) {
    return 'Submission ID not found. Either it is wrong, it has been manually deleted or it never made it to the database. If you are ' . 'certain this is a correct SUID, then '. l('contact', 'contact') .' the administrator.';
  }
  $table = theme_table(array('Feature name', 'Automatic annotation', 'Approval'), $table_rows, '', 'Features submitted for inclusion to database.');
  return "<p>Your submission has been submitted to the annotation queue." . " A senior curator will check them and approve or reject them via the ". l('annotation queue', 'genes4all/curate/features') .'</p>'. $table;
}

/**
 * Form to control (e.g. mark for deletion) an existing feature
 *
 * @param $form_state
 * @param $feature_id
 *
 * @return unknown_type
 */
function genes4all_curate_curate_control_feature_form($form_state, $feature_name = NULL) {
  $form = array();
  $feature_id = gmod_dbsf_chado_name2id($feature_name, 'feature');
  if (empty($feature_id)) {
    drupal_set_message(t('Resource %feature_name was not found', array('%feature_name' => $feature_name)), 'error');
    return $form;
  }
  $curator_suggestions = gmod_dbsf_db_execute('chado','gmod_dbsf_get_add_cv_withnames','curator_suggestions');

  $form = array(
    '#type' => 'fieldset',
    '#title' => "Decide upon a feature's fate",
    '#description' => t('You can request that certain actions are undertaken for a feature. None of these are immediate and must be approved by a main curator.'),
    'suggestions' => array(
      '#type' => 'checkboxes',
      '#options'=> $curator_suggestions
  ),
  );
  $form['button'] = array(
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Submit request'),
      '#weight' => 10,
  ), 'clear' => array(
      '#type' => 'button',
      '#value' => t('Clear form'),
      '#weight' => 11,
      '#validate' => array('gmod_dbsf_form_clear'),
  ),
  );
  $form['data']['#title'] =l($feature_name,'genes4all/info/feature',array('query'=>array('feature_id'=>$feature_name)));
  $form['data']['feature_id'] = array('#type' => 'value','#value' => $feature_id);
  return $form;
}

function genes4all_curate_curate_control_feature_form_submit($form, &$form_state) {
  $fid = check_plain($form_state['values']['feature_id']);
  $cvterm_array=array();
  foreach ($form_state['values']['suggestions'] as $cvterm_id){
    if (empty($cvterm_id)){continue;}
    $cvterm_array[]= check_plain($cvterm_id);
  }
  if (!empty($cvterm_array)){
    genes4all_curation_store_user_curation_action('feature',$fid,$cvterm_array);
  }else{
    drupal_set_message('No actions were requested!','error');
  }
}

/**
 * Form to add ontologies to an existing feature
 *
 * Adds entries to the feature_cvterm table of the chado db
 *
 * @param $form_state
 * @param $feature_id
 *
 * @return unknown_type
 */
function genes4all_curate_curate_basic_ontologies_feature_form($form_state,$feature_name=NULL) {
  $form = array();
  $feature_id = gmod_dbsf_chado_name2id($feature_name, 'feature');
  if (empty($feature_id)){return;}
  /*
   * We will use the following CVs
   * sequence,biological_process,molecular_function,cellular_component,cell,EC,KEGG_PATHWAY
   *
   * Also need the CV Evidence Codes
   */
  $basic_cv = array('sequence','biological_process','molecular_function','cellular_component','EC','KEGG_PATHWAY');
  $evidence_data = genes4all_curation_cv_evidence(NULL, TRUE, TRUE);
  $cv_data = cache_get('BASE_ONTOLOGIES','cache_genes4all_explore');
  $cv_descr =array();
  $desc_sql = "SELECT definition as descr from {cv} where name='%s'";
  if (empty($cv_data)){
    $cv_data = array();
    gmod_dbsf_db_execute('chado');
    foreach ($basic_cv as $cv_name){
      $cv_data[$cv_name] = array(0 => 'None') + gmod_dbsf_get_add_cv_withnames($cv_name,NULL,'FRIENDLY');
      if (!empty($cv_data[$cv_name])){
        $d = db_fetch_array(db_query($desc_sql,$cv_name));
        $cv_descr[$cv_name]=$d['descr'];
      }
    }
    gmod_dbsf_db_execute();
    if (!empty($cv_data)){
      foreach ($cv_data as $key => $cv) {
        foreach ($cv as $id => $name) {
          $name = str_replace('_', ' ', $name);
          $cv_data[$key][$id] = $name;
          if (strlen($name) > 85) {
            $cv_data[$key][$id] = substr($name, 0, 80) .' ...';
          }
        }
      }
    }
    cache_set('BASE_ONTOLOGIES', $cv_data , 'cache_genes4all_explore');
  } else{
    $cv_data =$cv_data->data;
  }

  $form[] = array(
  	'#title' =>l($feature_name,'genes4all/info/feature',array('query'=>array('feature_id'=>$feature_name))),
    '#type' => 'fieldset',
  	'feature_id' => array('#type' => 'value','#value' => $feature_id),
    'ontologies' => array(
      '#title' => 'Standard ontologies / Controlled Vocabularies (CVs)',
      '#type' => 'fieldset',
      '#tree' => TRUE,
      '#weight' => 2,
      'cv_sequence' => array(
        '#title' => 'Sequence Ontology',
        '#description' => $cv_descr['sequence'].t(' This CV can further define the type of sequence data.'),
        '#type' => 'select',
        '#size' => 5,
        '#options' => $cv_data['sequence'],
  ),
      'cv_go_bio' => array(
        '#title' => 'Gene Ontology: Biological Process',
        '#description' => $cv_descr['biological_process'].t(' This CV allows you to specify which process the sequence is part of/involved in.'),
        '#type' => 'select',
        '#size' => 5,
        '#options' => $cv_data['biological_process'],
  ),
      'cv_go_mol' => array(
        '#title' => 'Gene Ontology: Molecular Function',
        '#description' => $cv_descr['molecular_function'].t(' This CV allows you to specify what function this sequence performs'),
        '#type' => 'select',
        '#size' => 5,
        '#options' => $cv_data['molecular_function'],
  ),
      'cv_go_comp' => array(
        '#title' => 'Gene Ontology: Cellular Compartment',
        '#description' => $cv_descr['cellular_component'].t(' This CV can specify location(s) that the sequence feature is associated with.'),
        '#type' => 'select',
        '#size' => 5,
        '#options' => $cv_data['cellular_component'],
  ),
      'cv_ec' => array(
        '#title' => 'Enzyme Classification',
        '#description' => $cv_descr['EC'].t(' If you know what enzymatic reaction (if any) the feature catalyses, select an EC term'),
        '#type' => 'select',
        '#size' => 5,
        '#options' => $cv_data['EC'],
  ),
      'cv_kegg' => array(
        '#title' => 'KEGG Pathway term',
        '#description' => $cv_descr['KEGG_PATHWAY'].t(' The Kyoto Encyclopedia of Genes and Genomes provides more structured pathway maps the GO Biological process. If you know that this genes is involved in a specific pathway, then select the relevant term.'),
        '#type' => 'select',
        '#size' => 5,
        '#options' => $cv_data['KEGG_PATHWAY'],
  ),
  ),'evidence' => $evidence_data,
  'button' => array(
    '#weight'=>20,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Submit annotation'),
      '#weight' => 10,
  ), 'clear' => array(
      '#type' => 'button',
      '#value' => t('Clear form'),
      '#weight' => 11,
      '#validate' => array('gmod_dbsf_form_clear'),
  ),
  )
  );
  return $form;
}

function genes4all_curate_curate_basic_ontologies_feature_form_validate($form, &$form_state) {
  // text
  $evidence = check_plain($form_state['values']['evidence']);
  $pub_id = gmod_dbsf_db_execute('chado_edit', 'gmod_dbsf_get_add_pub', 'annotation', $evidence, NULL, 'ADD', TRUE);
  if (empty($pub_id)) {
    form_set_error('evidence', t('Cannot annotate without giving an evidence id.'));
  }
  $form_state['values']['pub_id'] = $pub_id[$evidence];
}

function genes4all_curate_curate_basic_ontologies_feature_form_submit($form, &$form_state) {
  $fid = check_plain($form_state['values']['feature_id']);
  $evidence = check_plain($form_state['values']['evidence']);
  $pub_id = check_plain($form_state['values']['pub_id']);

  //only add desc if it does not exist. i.e. set it as an update.
  $cvterm_array = array();
  foreach ( $form_state['values']['ontologies'] as $cvterm_id) {
    if (empty($cvterm_id)) {
      continue;
    }
    $cvterm_array[] = check_plain($cvterm_id);
  }
  if (!empty($cvterm_array)){
    gmod_dbsf_db_execute('chado_edit', 'gmod_dbsf_populate_cvterms', 'feature', $fid, $cvterm_array, NULL, $pub_id);
    genes4all_curation_store_user_curation_action('feature',$fid,$cvterm_array);
  }else{
    drupal_set_message('No actions were requested!','error');
  }
}

/**
 * Form to add ontologies to an existing feature
 *
 * Adds entries to the feature_cvterm table of the chado db
 *
 * @param $form_state
 * @param $feature_id
 *
 * @return unknown_type
 */
function genes4all_curate_curate_user_ontologies_feature_form($form_state,$feature_name=NULL) {
  $form = array();
  $feature_id = gmod_dbsf_chado_name2id($feature_name, 'feature');
  if (empty($feature_id)){return $form;}
  $evidence_data       = genes4all_curation_cv_evidence(NULL, TRUE, TRUE);
  $current_custom_cv = gmod_dbsf_get_add_var('genes4all_custom_cv_curate');
  $current_custom_cv = array_diff($current_custom_cv, array('EC','KEGG_PATHWAY'));
  if (empty($current_custom_cv)){return $form;}
  $cv_data = cache_get('CUSTOM_ONTOLOGIES','cache_genes4all_explore');
  $cv_descr =array();
  $desc_sql = "SELECT definition as descr from {cv} where name='%s'";

  if (empty($cv_data)){
    $cv_data = array();
    gmod_dbsf_db_execute('chado');
    foreach ($current_custom_cv as $cv_name){
      $cv_data[$cv_name] = array(0 => 'None') + gmod_dbsf_get_add_cv_withnames($cv_name);
      if (!empty($cv_data[$cv_name])){
        $d = db_fetch_array(db_query($desc_sql,$cv_name));
        $cv_descr[$cv_name]=$d['descr'];
      }
    }
    gmod_dbsf_db_execute();
  }
  if (!empty($cv_data)){
    foreach ($cv_data as $key => $cv) {
      foreach ($cv as $id => $name) {
        $name = str_replace('_', ' ', $name);
        $cv_data[$key][$id] = $name;
        if (strlen($name) > 85) {
          $cv_data[$key][$id] = substr($name, 0, 80) .' ...';
        }
      }
    }
    if (empty($cv_data)){return $form;}
    $day1 = time()+86400;
  }else{
    $cv_data =$cv_data->data;
    if (empty($cv_data)){return $form;}
  }
  cache_set('CUSTOM_ONTOLOGIES', $cv_data , 'cache_genes4all_explore',$day1);
  $form[] = array(
    '#type' => 'fieldset',
    'ontologies' => array(
      '#title' => 'User-provided ontologies / Controlled Vocabularies (CVs)',
      '#type' => 'fieldset',
      '#tree' => TRUE,
      '#weight' => 2,
  ),
  );
  foreach ($cv_data as $cv_name=>$data){
    $form['ontologies'][$cv_name] = array(
        '#type' => 'select',
        '#title' => $cv_name,
        '#description'=> $cv_descr[$cv_name],
        '#size' => 5,
        '#options' => $cv_data[$cv_name],
    );
  }


  $form['evidence'] = $evidence_data;
  $form['button'] = array(
  	'#weight' => 20,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Submit annotation'),
      '#weight' => 10,
  ), 'clear' => array(
      '#type' => 'button',
      '#value' => t('Clear form'),
      '#weight' => 11,
      '#validate' => array('gmod_dbsf_form_clear'),
  ),
  );
  $form['#title'] =l($feature_name,'genes4all/info/feature',array('query'=>array('feature_id'=>$feature_name)));
  $form['feature_id'] = array('#type' => 'value','#value' => $feature_id);
  return $form;
}

function genes4all_curate_curate_user_ontologies_feature_form_validate($form, &$form_state) {
  // text
  $evidence = check_plain($form_state['values']['evidence']);
  $pub_id = gmod_dbsf_db_execute('chado_edit', 'gmod_dbsf_get_add_pub', 'annotation', $evidence, NULL, 'ADD', TRUE);
  if (empty($pub_id)) {
    form_set_error('evidence', t('Cannot annotate without giving an evidence id.'));
  }
  $form_state['values']['pub_id'] = $pub_id[$evidence];
}

function genes4all_curate_curate_user_ontologies_feature_form_submit($form, &$form_state) {
  $fid = check_plain($form_state['values']['feature_id']);
  $evidence = check_plain($form_state['values']['evidence']);
  $pub_id = check_plain($form_state['values']['pub_id']);

  //only add desc if it does not exist. i.e. set it as an update.
  $cvterm_array = array();
  foreach ( $form_state['values']['ontologies'] as $cvterm_id) {
    if (empty($cvterm_id)) {
      continue;
    }
    $cvterm_array[] = check_plain($cvterm_id);
  }
  if (!empty($cvterm_array)){
    gmod_dbsf_db_execute('chado_edit', 'gmod_dbsf_populate_cvterms', 'feature', $fid, $cvterm_array, NULL, $pub_id);
    genes4all_curation_store_user_curation_action('feature',$fid,$cvterm_array);
  }else{
    drupal_set_message('No actions were requested!','error');
  }
}

/**
 * Form to add ontologies to an existing feature
 *
 * Adds entries to the feature_cvterm table of the chado db
 *
 * @param $form_state
 * @param $feature_id
 *
 * @return unknown_type
 */
function genes4all_curate_curate_user_dbxref_feature_form($form_state,$feature_name=NULL) {
  $form = array();
  $feature_id = gmod_dbsf_chado_name2id($feature_name, 'feature');
  if (empty($feature_id)){return $form;}
  $current_custom_db = gmod_dbsf_get_add_var('genes4all_custom_db_curate');
  if (empty($current_custom_db)){return $form;}
  $db_data = cache_get('CUSTOM_DBXREF','cache_genes4all_explore');

  if (empty($db_data)){
    $db_data = array();
    $desc_sql = "SELECT description as descr from {db} where name='%s'";
    gmod_dbsf_db_execute('chado');
    foreach ($current_custom_db as $db_name){
      $db_data[$db_name]['data'] = array(0 => 'None') + gmod_dbsf_get_add_dbxref_withnames($db_name,NULL,NULL,'FRIENDLY');
      if (!empty($db_data[$db_name])){
        $d = db_fetch_array(db_query($desc_sql,$db_name));
        $db_data[$db_name]['descr']=$d['descr'];
      }
    }
    gmod_dbsf_db_execute();
    if (!empty($db_data)){
      foreach ($db_data as $key => $db) {
        foreach ($db as $id => $name) {
          $name = str_replace('_', ' ', $name);
          $db_data[$key][$id] = $name;
          if (strlen($name) > 85) {
            $db_data[$key][$id] = substr($name, 0, 80) .' ...';
          }
        }
      }
    }
    $day1 = time()+86400;
    if (empty($db_data)){return $form;}
    cache_set('CUSTOM_DBXREF', $db_data , 'cache_genes4all_explore',$day1);
  } else{
    $db_data =$db_data->data;
    if (empty($db_data)){return $form;}
  }

  $form['data'] = array(
    '#type' => 'fieldset',
    'dbxref' => array(
      '#title' => t('User-provided DBXREFs'),
      '#description'=>t('external database cross-references'),
      '#type' => 'fieldset',
      '#tree' => TRUE,
      '#weight' => 2,
  ),
  );
  foreach ($db_data as $db_name=>$data){
    $form['data']['dbxref'][$db_name] = array(
        '#type' => 'select',
        '#title' => $db_name,
        '#description'=> $db_data[$db_name]['descr'],
        '#size' => 25,
        '#options' => $db_data[$db_name]['data'],
    );
  }


  $form['data']['button'] = array(
  '#weight' => 20,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Submit annotation'),
      '#weight' => 10,
  ), 'clear' => array(
      '#type' => 'button',
      '#value' => t('Clear form'),
      '#weight' => 11,
      '#validate' => array('gmod_dbsf_form_clear'),
  ),
  );
  $form['#title'] =l($feature_name,'genes4all/info/feature',array('query'=>array('feature_id'=>$feature_name)));
  $form['feature_id'] = array('#type' => 'value','#value' => $feature_id);
  return $form;
}

function genes4all_curate_curate_user_dbxref_feature_form_submit($form, &$form_state) {
  $fid = check_plain($form_state['values']['feature_id']);
  //only add desc if it does not exist. i.e. set it as an update.
  $dbxref_array = array();
  foreach ( $form_state['values']['dbxref'] as $dbxref_id) {
    if (empty($dbxref_id)) {
      continue;
    }
    $dbxref_array[] = check_plain($dbxref_id);
  }
  if (!empty($dbxref_array)){
    gmod_dbsf_db_execute('chado_edit', 'gmod_dbsf_populate_dbxrefs', 'feature', $fid, $dbxref_array);
    genes4all_curation_store_user_curation_action('feature',$fid,$dbxref_array,FALSE,TRUE);
  }else{
    drupal_set_message('No actions were requested!','error');
  }
}


/**
 * Form to add ontologies to an existing feature
 *
 * Adds entries to the feature_cvterm table of the chado db
 *
 * @param $form_state
 * @param $feature_id
 *
 * @return unknown_type
 */
function genes4all_curate_curate_free_ontologies_feature_form($form_state,$feature_name=NULL) {
  $form = array();
  $feature_id = gmod_dbsf_chado_name2id($feature_name, 'feature');
  if (empty($feature_id)){return;}
  /*
   * We will use the following CVs
   * sequence,biological_process,molecular_function,cellular_component,cell,EC,KEGG_PATHWAY
   *
   * Also need the CV Evidence Codes
   */
  $evidence_data       = genes4all_curation_cv_evidence(NULL, TRUE, TRUE);
  $form[]= array(
    '#title' =>l($feature_name,'genes4all/info/feature',array('query'=>array('feature_id'=>$feature_name))),
    '#type' => 'fieldset',
	'feature_id' => array('#type' => 'value','#value' => $feature_id),
  	'free_terms' => array(
      '#title' => 'Free terms',
      '#type' => 'fieldset',
      '#description' => t('Free terms allow you to give your own tags in order to drive future searches. If the tag is novel, then it will not appear in the autocomplete below; in that case we highly recommend you provide a description if the term is not self-explanatory.'),
      '#tree' => TRUE,
      '#weight' => 5,
      'free_term' => array(
        '#title' => 'Free term',
        '#description' => t('You can give any term which you think this feature can be associated with. Terms that have been already defined in the database (not necessarily for this feature) will appear in the autocomplete.'),
        '#type' => 'textfield',
        '#autocomplete_path' => 'genes4all/curate/free_term/autocomplete',
        '#weight' => 2,
  ),
      'free_term_def' => array(
        '#title' => 'Term definition',
        '#description' => t('Optionally, you may add a definition for your term (recommended if name is not self-explanatory).'),
        '#type' => 'textarea',
        '#weight' => 4,
  ), 'free_term_def_ow' => array(
        '#title' => 'Overwrite definition if it already exists?',
        '#description' => t('Please check this only if you understand that it will delete any previous definition someone else might have given for that particular term (i.e. is applied across all features that may already contain this term). If there is no definition provided for this term, and you do not check this box, then your definition will still be stored. If a term has not appeared in the autocomplete, then it does not exists and this checkbox is irrelevant.'),
        '#type' => 'checkbox',
        '#weight' => 5,
  ),
  ),
  'evidence' => $evidence_data,
  'button' => array(
    '#weight' => 20,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Submit annotation'),
      '#weight' => 10,
  ), 'clear' => array(
      '#type' => 'button',
      '#value' => t('Clear form'),
      '#weight' => 11,
      '#validate' => array('gmod_dbsf_form_clear'),
  ),
  )
  );
  return $form;
}

function genes4all_curate_curate_free_ontologies_feature_form_validate($form, &$form_state) {
  // text
  $evidence = check_plain($form_state['values']['evidence']);
  $pub_id = gmod_dbsf_db_execute('chado_edit', 'gmod_dbsf_get_add_pub', 'annotation', $evidence, NULL, 'ADD', TRUE);
  if (empty($pub_id)) {
    form_set_error('evidence', t('Cannot annotate without giving an evidence id.'));
  }
  $form_state['values']['pub_id'] = $pub_id[$evidence];
  $form_state['values']['free_terms']['free_term_def'] = check_plain(trim($form_state['values']['free_terms']['free_term_def']));
  $form_state['values']['free_terms']['free_term'] = check_plain(trim($form_state['values']['free_terms']['free_term']));
}

function genes4all_curate_curate_free_ontologies_feature_form_submit($form, &$form_state) {
  $fid = check_plain($form_state['values']['feature_id']);
  // text
  $evidence = check_plain($form_state['values']['evidence']);
  // return array of cvterm_id as values (key is unimportant)
  $free_cvterm      = $form_state['values']['free_terms']['free_term'];
  $free_cvterm_desc = $form_state['values']['free_terms']['free_term_def'];
  // 1 if overwrite
  $free_desc_overwrite = check_plain($form_state['values']['free_terms']['free_term_def_ow']);
  $pub_id = check_plain($form_state['values']['pub_id']);

  //only add desc if it does not exist. i.e. set it as an update.
  $cvterm_array = array();
  // add free cv
  if (!empty($free_cvterm)) {
    $free_cvterm_id = 0;
    if (!empty($free_cvterm_desc)) {
      if (!empty($free_desc_overwrite)) {
        $free_cvterm_id = gmod_dbsf_db_execute('chado_edit', 'gmod_dbsf_get_add_cv_withnames', 'free_curation_terms', $free_cvterm, FALSE, $free_cvterm_desc, TRUE);
      }
      else {
        $free_cvterm_id = gmod_dbsf_db_execute('chado_edit', 'gmod_dbsf_get_add_cv_withnames', 'free_curation_terms', $free_cvterm, FALSE, $free_cvterm_desc);
      }
    }
    else {
      $free_cvterm_id = gmod_dbsf_db_execute('chado_edit', 'gmod_dbsf_get_add_cv_withnames', 'free_curation_terms', $free_cvterm);
    }
    if (!empty($free_cvterm_id)) {
      $cvterm_array[] = check_plain($free_cvterm_id);
    }
  }
  gmod_dbsf_db_execute('chado_edit', 'gmod_dbsf_populate_cvterms', 'feature', $fid, $cvterm_array, NULL, $pub_id);
  genes4all_curation_store_user_curation_action('feature',$fid,$cvterm_array);
}


function genes4all_curate_load_approved($gff_files = NULL) {
  //TODO:: Add to administratable variables.
  $dbprefix = 'ic';
  global $db_url;
  preg_match('/^pgsql:\/\/(.+):(.+)\@(.+):(\d+)\/(.+)/', $db_url['chado_edit'], $match);
  $dsn = "'" . 'dbi:Pg:dbname='. $match[5] .';host='. $match[3] .';port='. $match[4] .';user='. $match[1] .';password='. $match[2] ."'";
  $db_array = array(
    'dbnames' => array(),
    'host' => $match[3],
    'port' => $match[4],
    'user' => $match[1],
    'password' => $match[2],
  );
  if (empty($dsn)) {
    drupal_set_message(t('I cannot get the Chado database connection settings in order to construct a DSN entry for the GMOD scripts.'), 'error');
    return;
  }
  //Drupal specific SQLs
  $sql_parent_select_loc = "SELECT f.uniquename as uname,f.feature_id as parent_id from {gmod_dbsf_featureloc} fl JOIN {gmod_dbsf_feature} f " . " ON f.feature_id=fl.srcfeature_id where fl.feature_id=(SELECT feature_id from {gmod_dbsf_feature} where uniquename='%s')";
  $sql_children_select_loc = "SELECT f.uniquename as uname,f.feature_id as child_id from {gmod_dbsf_featureloc} fl JOIN {gmod_dbsf_feature} f " . " ON f.feature_id=fl.feature_id where fl.srcfeature_id=(SELECT feature_id from {gmod_dbsf_feature} where uniquename='%s')";

  if (!empty($gff_files)) {
    $dirpath = base_path() . file_create_path() .'/genes4all_curate/';
    $tuid = gmod_dbsf_create_uid(NULL, NULL, 'gff_load');
    foreach ($gff_files as $ncbi_taxid => $features) {
      $db_array['dbnames'][] = $dbprefix .'_'. $ncbi_taxid;
      $outfile[$ncbi_taxid] = $dirpath .'/'. $tuid .".$ncbi_taxid.temp";
      $towrite_gff[$ncbi_taxid] = "##gff-version 3\n";
      $towrite_fsa[$ncbi_taxid] = "\n##FASTA\n";
    }
    foreach ($gff_files as $ncbi_taxid => $features) {
      foreach ($features as $feature_name => $infile_array) {
        $parent_data_res = db_query($sql_parent_select_loc, $feature_name);
        $children_data_res = db_query($sql_children_select_loc, $feature_name);
        //add feature parents
        while ($parent_data = db_fetch_array($parent_data_res)) {
          $check = gmod_dbsf_db_execute('chado', 'gmod_dbsf_get_feature_name_id', $parent_data['uname']);
          if (empty($check)) {
            $towrite_array = gmod_dbsf_db_execute('chado', 'gmod_dbsf_chado_feature2gff', $parent_data['uname']);
            $towrite_gff[$ncbi_taxid] .= $towrite_array['gff'];
            $towrite_fsa[$ncbi_taxid] .= $towrite_array['fsa'];
          }
        }
        // add features if not existing (shouldn't exist!)
        $check = gmod_dbsf_db_execute('chado', 'gmod_dbsf_get_feature_name_id', $feature_name);
        if (empty($check)) {
          $towrite_array = gmod_dbsf_db_execute('chado', 'gmod_dbsf_chado_feature2gff', $feature_name);
          $towrite_gff[$ncbi_taxid] .= $towrite_array['gff'];
          $towrite_fsa[$ncbi_taxid] .= $towrite_array['fsa'];
        }
        // add their children if not existing
        while ($children_data = db_fetch_array($children_data_res)) {
          $check = gmod_dbsf_db_execute('chado', 'gmod_dbsf_get_feature_name_id', $children_data['uname']);
          if (empty($check)) {
            $towrite_array = gmod_dbsf_db_execute('chado', 'gmod_dbsf_chado_feature2gff', $children_data['uname']);
            $towrite_gff[$ncbi_taxid] .= $towrite_array['gff'];
            $towrite_fsa[$ncbi_taxid] .= $towrite_array['fsa'];
          }
        }
      }
      if (!$outhandle[$ncbi_taxid] = fopen($outfile[$ncbi_taxid], 'ab')) {
        drupal_set_message(t('Could not create %outfile.', array('%outfile' => $outfile)), 'error');
      }
      fwrite($outhandle[$ncbi_taxid], $towrite_gff[$ncbi_taxid]);
      $towrite_gff[$ncbi_taxid] = '';
    }
    $operations = array();
    foreach ($gff_files as $ncbi_taxid => $features) {
      foreach ($features as $feature_name => $infile_array) {
        foreach ($infile_array as $infile) {
          if (file_exists($infile)) {
            if ($inhandle = fopen($infile, 'rb')) {
              while (!feof($inhandle)) {
                $line = fgets($inhandle);
                if (strpos($line, '#') !== 0) {
                  fwrite($outhandle[$ncbi_taxid], $line);
                }
              }
            }
            fclose($inhandle);
          }
        }
      }
      fwrite($outhandle[$ncbi_taxid], $towrite_fsa[$ncbi_taxid]);
      $towrite_fsa[$ncbi_taxid] = '';
      fclose($outhandle[$ncbi_taxid]);
    }
    $operations[] = array('gmod_dbsf_chado_gff2seqfeature', array($db_array, $outfile));
    $operations[] = array('gmod_dbsf_chado_gff2chado', array($dsn, $gff_files));
    $batch        = array(
      'title' => t('Checking for finished annotation...'),
      'operations' => $operations,
      'init_message' => t('Uploading annotation to Chado...'),
      'progress_message' => t('@remaining operations remaining...'),
      'error_message' => t('Batch system encountered an error.'),
      'finished' => 'genes4all_curate_batch_annot8_finished',
      'progressive' => FALSE,
      'file' => drupal_get_path('module', 'genes4all_curate') .'/includes/genes4all_curate_feature.inc',
    );
    batch_set($batch);
    batch_process();
  }
}

function genes4all_curate_batch_annot8_finished($success, $results, $operations) {
  // store uids to report in case of error.
  $uid_array = array();
  $dirpath   = file_create_path() .'/genes4all_curate';
  $message   = '';
  $status    = 'warning';
  if (!empty($success)) {
    $message = t('Annotations uploaded.');
  }
  else {
    $message = t('Something went wrong, please '. l('contact', 'contact') .' an administrator.');
  }
  drupal_set_message($message, $status);
}

function genes4all_curate_check4annotated_features() {
  /*
   * annotated=0: not finished, not checked
   * annotated=1 finished
   * annotated=2 failed/non-existing uid.
   *
   */


  global $db_url;
  preg_match('/^pgsql:\/\/(.+):(.+)\@(.+):(\d+)\/(.+)/', $db_url['chado_edit'], $match);
  $dsn = "'" . 'dbi:Pg:dbname='. $match[5] .';host='. $match[3] .';port='. $match[4] .';user='. $match[1] .';password='. $match[2] ."'";

  $dirpath  = file_create_path() .'/genes4all_curate/';
  $fullpath = base_path() . $dirpath;
  $execdir  = './'. drupal_get_path('module', 'genes4all') .'/scripts/';
  // this is a parent page, therefore we check if this exists, so that an error is reported. we do not use it in this page though.
  // must specify organism explicitly.
  //$gmod_bulk_load = $execdir.'/gmod_bulk_load_gff3-AP.pl';
  $ic_annot8r2gff    = $execdir .'ic_annot8r2gff.pl';
  $ic_blast2gff      = $execdir .'ic_blast2gff.pl';
  $seqfeat_load_exec = $execdir .'bp_seqfeature_load.pl';
  $chado_load_exec   = $execdir .'ic_chado_loadcv.pl';
  if (!file_exists($seqfeat_load_exec)) {
    drupal_set_message(t('Failed to find executable %exec', array('%exec' => $seqfeat_load_exec)), 'error');
    return FALSE;
  }
  elseif (!is_executable($seqfeat_load_exec)) {
    drupal_set_message(t('Failed to execute %exec', array('%exec' => $seqfeat_load_exec)), 'error');
    return FALSE;
  }
  if (!file_exists($chado_load_exec)) {
    drupal_set_message(t('Failed to find executable %exec', array('%exec' => $chado_load_exec)), 'error');
    //return FALSE;
  }
  elseif (!is_executable($chado_load_exec)) {
    drupal_set_message(t('Failed to execute %exec', array('%exec' => $chado_load_exec)), 'error');
    //return FALSE;
  }
  if (!file_exists($ic_annot8r2gff)) {
    drupal_set_message(t('Cannot find program ic_annot8r2gff in module directory. Will not parse annotation results.'), 'error');
    return FALSE;
  }
  elseif (!is_executable($ic_annot8r2gff)) {
    drupal_set_message(t('Failed to execute %exec', array('%exec' => $ic_annot8r2gff)), 'error');
    return FALSE;
  }
  if (!file_exists($ic_blast2gff)) {
    drupal_set_message(t('Cannot find program ic_blast2gff in module directory. Will not parse annotation results.'), 'error');
    return FALSE;
  }
  elseif (!is_executable($ic_blast2gff)) {
    drupal_set_message(t('Failed to execute %exec', array('%exec' => $ic_blast2gff)), 'error');
    return FALSE;
  }
  /*if (!file_exists($gmod_bulk_load)){
   drupal_set_message('Cannot find program gmod_bulk_load_gff3-AP.pl in module directory. Will not parse annotation results.','error');
   return;
   }*/


  $sql_org_select = "SELECT abbreviation from {gmod_dbsf_organism} where organism_id=(SELECT organism_id from {gmod_dbsf_feature} where uniquename='%s')";
  $sql_select     = "SELECT * from {genes4all_approval_queue} where approved is NULL and annotated=0";
  $sql_update     = "UPDATE {genes4all_approval_queue} set annotated=%d where approval_queue_id=%d";
  $db_res         = db_query($sql_select);
  // All the new GFF files
  $gff_files = array();
  while ($row = db_fetch_array($db_res)) {
    $feature_name = $row['uniquename'];
    $uid          = $row['annotation_id'];
    $id           = $row['approval_queue_id'];
    //cehck for non-existing files, mark them so don't check again
    $infile = $dirpath . $uid .'.query';
    if (!file_exists($infile)) {
      drupal_set_message(t('Found no input files for searches set to run for %feature_name', array('%feature_name' => $feature_name)), 'error');
      db_query($sql_update, 2, $id);
      continue;
    }
    //check for errors
    $f = file_scan_directory($dirpath, $uid .'.+\.error$', array('.', '..', 'CVS'), 0, FALSE, 'filename');
    foreach ($f as $infile => $obj) {
      if (filesize($infile) > 1) {
        db_query($sql_update, 2, $id);
        continue 2;
      }
    }
    $outfile_number = 0;
    //done files
    $f = file_scan_directory($dirpath, $uid .'.+\.done$', array('.', '..', 'CVS'), 0, FALSE, 'filename');
    foreach ($f as $infile => $obj) {
      $inhandle = fopen($infile, 'rb');
      while (!feof($inhandle)) {
        $line = trim(fgets($inhandle));
        if (!empty($line) && preg_match('/^par/', $line)) {
          $outfile_number++;
        }
      }
    }
    //wait files
    $f = file_scan_directory($dirpath, $uid .'.+\.wait$', array('.', '..', 'CVS'), 0, FALSE, 'filename');
    foreach ($f as $infile => $obj) {
      $inhandle = fopen($infile, 'rb');
      while (!feof($inhandle)) {
        $line = trim(fgets($inhandle));
        if (!empty($line) && preg_match('/^par/', $line)) {
          $outfile_number++;
        }
      }
    }
    unset($f);
    if (empty($outfile_number)) {
      drupal_set_message(t('Found no parameter file for searches set to run for %feature_name.', array('%feature_name' => $feature_name)), 'error');
      db_query($sql_update, 2, $id);
      continue;
    }
    $finished_outs = 0;
    //check for (xml & raw) outfile(s)
    $outfiles = file_scan_directory($fullpath, $uid .'.+\.output$', array('.', '..', 'CVS'), 0, FALSE, 'filename');

    foreach ($outfiles as $fname => $obj) {
      if (filesize($fname) < 1 && strpos($fname, 'iprscan') === FALSE) {
        continue;
      }
      $finished_outs++;
    }
    if ($finished_outs === $outfile_number) {
      // then feature is ready.
      // first find organism abbreviation (Genus species=>G.species) from Drupal db
      $org_abbr = db_fetch_array(db_query($sql_org_select, $feature_name));
      $org_abbr = $org_abbr['abbreviation'];
      if (empty($org_abbr)) {
        drupal_set_message(t('Cannot find species abbreviation for organism of %feature_name.', array('%feature_name' => $feature_name)), 'error');
        db_query($sql_update, 2, $id);
        continue;
      }
      // first convert to GFF.
      $error = array();
      foreach ($outfiles as $fname => $obj) {
        if (strpos($fname, 'iprscan') !== FALSE) {
          $exec = $ic_annot8r2gff ." -dsn $dsn -org $org_abbr -outdir $fullpath $fname ";
          exec($exec ." >/dev/null 2>&1", $error);
        }
        elseif (strpos($fname, 'blastall') !== FALSE) {
          $exec = $ic_blast2gff ." -dsn $dsn -org $org_abbr -cs 60 -ce 1e-4 -l 10 -format blastxml -outdir $fullpath  -b $fname ";
          exec($exec ." >/dev/null 2>&1", $error);
        }
      }
      // if error set to black
      if (!empty($error)) {
        drupal_set_message(t('An error occuring while attempting to load %feature_name into Chado.', array('%feature_name' => $feature_name)), 'error');
        db_query($sql_update, 2, $id);
        continue;
      }
      //set to green
      db_query($sql_update, 1, $id);
    }
  }
}

function genes4all_curate_add_approval_queue($uniquename, $annot_id, $subid) {
  global $user;
  $subname = check_plain($user->name);
  if (empty($subname)) {
    drupal_set_message(t('Warning, guest users are not be allowed to submit new features: we will not know who did it! Please login.'), 'error');
    return FALSE;
  }
  $subid      = check_plain($subid);
  $annot_id   = check_plain($annot_id);
  $uniquename = check_plain($uniquename);
  if (empty($uniquename) || empty($subname) || empty($subid) || empty($annot_id)) {
    return;
  }
  $sql_insert = "INSERT INTO {genes4all_approval_queue} (uniquename,submitter_name,annotation_id,submission_id) VALUES ('%s','%s','%s','%s')";
  db_query($sql_insert, $uniquename, $subname, $annot_id, $subid);
  return TRUE;
}

