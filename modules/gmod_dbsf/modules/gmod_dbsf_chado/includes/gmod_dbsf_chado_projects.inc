<?php

/**
 * @file
 * Administration of gmod_dbsf_chado
 * @author Alexie Papanicolaou <Alexie.Papanicolaou@csiro.au>
 * @see http://insectacentral.org
 */

function gmod_dbsf_chado_page_project_admin() {
  /*
   * we want to be able to do the following:
   * 1 add a project
   * 2 add a contact
   * 3 edit a project
   * 4 edit a contact
   * 5 link a contact with a project
   */
  drupal_add_js(drupal_get_path('module', 'gmod_dbsf') .'/includes/gmod_dbsf_sorttable.js');
  $html = 'Administer geolocation-enabled genetics projects';
  $add_project_form = drupal_get_form('gmod_dbsf_chado_project_new_form');
  $edit_project_form = drupal_get_form('gmod_dbsf_chado_project_edit_form');
  $add_contact_form = drupal_get_form('gmod_dbsf_chado_contact_new_form');
  $edit_contact_form = drupal_get_form('gmod_dbsf_chado_contact_edit_form');
  $link_project_contact_form = drupal_get_form('gmod_dbsf_chado_link_project_contact_form');
  $add_project_data_form = drupal_get_form('gmod_dbsf_chado_add_project_data_form');

  $edit_project_form .= gmod_dbsf_chado_project_table();
  $edit_contact_form .= gmod_dbsf_chado_contact_table();
  $link_project_contact_form .= gmod_dbsf_chado_project_contact_table();

  $page_tabs = array(
    'administration' => array(
      '#type' => 'tabset',
   'new_project' => array(
        '#type' => 'tabpage',
        '#title' => 'New project',
        '#content' => $add_project_form,
        '#weight' => 1,
  ),
  'edit_project'  => array(
        '#type' => 'tabpage',
        '#title' => 'Edit project',
        '#content' => $edit_project_form,
        '#weight' => 2,
  ),
  'add_contact'  => array(
        '#type' => 'tabpage',
        '#title' => 'New contact',
        '#content' => $add_contact_form,
        '#weight' => -2,
  ),
  'edit_contact'  => array(
        '#type' => 'tabpage',
        '#title' => 'Edit contact',
        '#content' => $edit_contact_form,
        '#weight' => -1,
  ),
  'contact_project'  => array(
        '#type' => 'tabpage',
        '#title' => 'Link contacts with projects',
        '#content' => $link_project_contact_form,
        '#weight' => 3,
  ),
  'project_data'  => array(
        '#type' => 'tabpage',
        '#title' => 'Add data to project',
        '#content' => $add_project_data_form,
        '#weight' => 3,
  ),
  ),
  );
  $html .= tabs_render($page_tabs);
  return $html;
}

function gmod_dbsf_chado_project_new_form($form_state){
  $form = array();

  $form['new_project']=array(
    '#type' => 'fieldset',
    '#tree'=>TRUE,
    '#title' => t('New Chado Project'),
    '#description' => t('Register a new project in Chado, e.g. one involving collecting, genotyping and/or phenotyping experiments.'),
      'name' => array(
      '#type'=>'textfield',
      '#maxlength'=>255,
      '#title'=>t('Name'),
      '#required'=>TRUE,
      '#description'=>t('A short name to specify the overall activity of the project. It will be used across the community so try to be informative yet succinct.'),
      '#weight'=>0,
  ),
    'description'=>array(
    	'#type'=>'textfield',
    	'#title'=>t('Description of project'),
    	'#maxlength'=>255,
    	'#required'=>TRUE,
    	'#description'=>t('Describe your project.'),
    	'#weight'=>1,
  ),
  );
  //buttons
  $form['buttons']=array(
  	'#weight' => 15,
  	'submit' => array(
       '#type' => 'submit',
       '#value' => t('Submit data'),
       '#weight' => 1,
  ),
  	'clear' => array(
        '#attributes' => array('title' => t('Clear the form')),
        '#type' => 'submit',
        '#value' => 'Reset this form',
        '#submit' => array('gmod_dbsf_form_clear'),
        '#weight' => 2,
  )
  );

  return $form;
}

function gmod_dbsf_chado_project_new_form_validate($form,&$form_state){
  //data is stored at $form_state['values']
  // property tree in the fieldset above allows us to store all data in a structure tree format under $form_state['values']['new_project'] etc. otherwise it would be flat (under 'values').
  // for every user input, sanitize (hence passing form_state by reference above)
  $form_state['values']['new_project']['name'] = trim(check_plain($form_state['values']['new_project']['name']));
  $form_state['values']['new_project']['description'] = trim(check_plain($form_state['values']['new_project']['description']));

}

function gmod_dbsf_chado_project_new_form_submit($form,&$form_state){
  gmod_dbsf_db_execute('chado_edit','gmod_dbsf_chado_get_add_project',$form_state['values']['new_project']['name'],$form_state['values']['new_project']['description']);
}

function gmod_dbsf_chado_project_edit_form($form_state){
  $form =array();
  //get all available projects
  $projects = gmod_dbsf_db_execute('chado','gmod_dbsf_chado_get_projects');
  if (empty($projects)){
    drupal_set_message('No projects have been defined!','warning',FALSE);
    return $form;
  }
  // format arrays into options
  $projects_options = array();
  foreach ($projects as $id=>$data){
    $projects_options[$id] = $data['name'];
  }
  $projects = $projects_options; unset($projects_options);
  //buttons
  $form['delete'] = array(
      '#type'=>'select',
      '#title'=>t('Available projects'),
      '#options' => $projects,
      '#size' => 1,
      '#weight'=>0,
  );
  $form['buttons']=array(
  	'#weight' => 15,
  	'submit' => array(
       '#type' => 'submit',
       '#value' => t('Delete'),
       '#weight' => 1,
  ));

  return $form;
}

function gmod_dbsf_chado_project_edit_form_submit($form,&$form_state){
  $delete = check_plain($form_state['values']['delete']);
  if (empty($delete)){return;}
  gmod_dbsf_db_execute('chado_edit','gmod_dbsf_chado_get_add_project',$delete,NULL,TRUE);
}

function gmod_dbsf_chado_contact_new_form($form_state){
  $form =array();

  //@TODO AP: move to .install
  $contact_type_array = array('Principal Investigator','Communicating Author','Primary Author','Curator','Data Entryist');
  $contact_type_array = gmod_dbsf_db_execute('chado_edit','gmod_dbsf_get_add_cv_withnames','project_contact_types',$contact_type_array);

  //end move

  $form['contact_details']=array(
    '#type' => 'fieldset',
    '#tree'=>TRUE,
    '#title' => t('New Contact'),
  	  '#description'=>t('Details for a person associated with a project.'),
        'contact_type'=>array(
  	        '#type'=>'select',
  	        '#title'=>t('Type'),
  	        '#required'=>TRUE,
            '#options'=>$contact_type_array,
  	        '#weight'=>0,
  ),
  	    'prefix'=>array(
  	        '#type'=>'textfield',
  	        '#title'=>t('Prefix'),
  			'#description'=>t('Give any name prefix.'),
  	        '#weight'=>1,
  ),
  		'first_name'=>array(
  	        '#type'=>'textfield',
  	        '#title'=>t('First name'),
  			'#description'=>t('Write in full your proper first name. For example if your name is %a %b %c but you go by the name %b, then give %a.',array('%a'=>'William','%b'=>'Robert','%c'=>'Smith')),
  	        '#required'=>TRUE,
  	        '#weight'=>2,
  ),
  	    'family_name'=>array(
  	        '#type'=>'textfield',
  	        '#title'=>t('Family name'),
  	        '#description'=>t('You surname/family name (e.g. %c in the above example).',array('%c'=>'Smith')),
  	        '#required'=>TRUE,
  	        '#weight'=>3,
  ),
  	    'middle_name'=>array(
  	        '#type'=>'textfield',
  	        '#title'=>t('Any middle names'),
            '#description'=>t('Write in full any middle names you may have (e.g. %b in the above example). Use space to delimit multiple ones.',array('%b'=>'Robert')),
  	      	'#weight'=>4,
  ),
    	'suffix'=>array(
  	        '#type'=>'textfield',
  	        '#title'=>t('Suffix'),
  			'#description'=>t('Give any name suffixes.'),
  	        '#weight'=>5,
  ),

  );
  //buttons
  $form['buttons']=array(
  	'#weight' => 15,
  	'submit' => array(
       '#type' => 'submit',
       '#value' => t('Submit data'),
       '#weight' => 1,
  ),
  	'clear' => array(
        '#attributes' => array('title' => t('Clear the form')),
        '#type' => 'submit',
        '#value' => 'Reset this form',
        '#submit' => array('gmod_dbsf_form_clear'),
        '#weight' => 2,
  )
  );
  return $form;

}

function gmod_dbsf_chado_contact_new_form_validate($form,&$form_state){
  //the data needs to be an integer
  $full_name =
  check_plain(
  implode(' ',array(
  $form_state['values']['contact_details']['first_name'],
  $form_state['values']['contact_details']['prefix'],
  $form_state['values']['contact_details']['middle_name'],
  $form_state['values']['contact_details']['family_name'],
  $form_state['values']['contact_details']['suffix']
  )
  )
  );
  $form_state['values']['contact_details']['full_name'] = trim(preg_replace('/\s\s+/', ' ', $full_name));

  if (!is_numeric($form_state['values']['contact_details']['contact_type'])){
    form_set_error('contact_details][contact_type','Invalid method to add data.');
  }
}

function gmod_dbsf_chado_contact_new_form_submit($form,&$form_state){
  //$form_state['values']['contact_details']['full_name']
  //$form_state['values']['contact_details']['contact_type']
  if (!empty($form_state['values']['contact_details']['description'])){
    gmod_dbsf_db_execute('chado_edit','gmod_dbsf_chado_get_add_contact',$form_state['values']['contact_details']['full_name'],$form_state['values']['contact_details']['contact_type'],$form_state['values']['contact_details']['description']);
  }else {
    gmod_dbsf_db_execute('chado_edit','gmod_dbsf_chado_get_add_contact',$form_state['values']['contact_details']['full_name'],$form_state['values']['contact_details']['contact_type']);
  }

}

function gmod_dbsf_chado_contact_edit_form($form_state){
  $form = array();
  $contacts = gmod_dbsf_db_execute('chado','gmod_dbsf_chado_get_contacts');
  if (empty($contacts)){
    drupal_set_message('No contacts have been defined!','warning',FALSE);
    return $form;
  }
  // format arrays into options
  $contacts2 = array();
  foreach ($contacts as $id=>$data){
    $contacts2[$data['type_name']][$id] = $data['name'];
  }
  $contacts = $contacts2; unset($contacts2);

  $form['delete'] = array(
      '#type'=>'select',
      '#title'=>t('Available contacts'),
      '#options' => $contacts,
      '#size' => 1,
      '#weight'=>0,
  );

  //buttons
  $form['buttons']=array(
  	'#weight' => 15,
  	'submit' => array(
       '#type' => 'submit',
       '#value' => t('Delete'),
       '#weight' => 1,
  ));
  return $form;
}

function gmod_dbsf_chado_contact_edit_form_submit($form,&$form_state){
  $delete = check_plain($form_state['values']['delete']);
  if (empty($delete)){return;}
  gmod_dbsf_db_execute('chado_edit','gmod_dbsf_chado_get_add_contact',$delete,NULL,NULL,TRUE);
}

function gmod_dbsf_chado_contact_table(){
  $table_caption    = 'Registered contacts';
  $table_header     = array('Type','Name', 'Description');
  $table_rows       = array();
  $table_attributes = array('class' => 'vertical sortable','width' => '750px');
  $table = '';
  $contacts = gmod_dbsf_db_execute('chado','gmod_dbsf_chado_get_contacts');
  if (empty($contacts)){
    return $table;
  }
  foreach ($contacts as $id=>$data){
    $table_rows[] = array($data['type_name'], $data['name'], $data['description']);
  }
  $table = theme_table($table_header, $table_rows, $table_attributes, $table_caption);
  return $table;
}

function gmod_dbsf_chado_project_contact_table(){
  $table_caption    = 'Contacts linked with projects';
  $table_header     = array('Project name', 'Contacts');
  $table_rows       = array();
  $table_attributes = array('class' => 'vertical sortable','width' => '750px');
  $table = '';
  $linked_data = gmod_dbsf_db_execute('chado','gmod_dbsf_chado_get_project_contacts');
  if (empty($linked_data)){
    return $table;
  }
  $array = array();
  foreach ($linked_data as $id=>$data){
    $array[$data['project_name']][] = $data['contact_name'];
  }
  foreach ($array as $project=>$data){
    $table_rows[] = array($project,implode(', ',$data));
  }
  if (!empty($table_rows)){
    $table = theme_table($table_header, $table_rows, $table_attributes, $table_caption);
  }
  return $table;
}

function gmod_dbsf_chado_project_table(){
  $table_caption    = 'Registered projects';
  $table_header     = array('Name', 'Description');
  $table_rows       = array();
  $table_attributes = array('class' => 'vertical sortable','width' => '750px');
  $table = '';
  $projects = gmod_dbsf_db_execute('chado','gmod_dbsf_chado_get_projects');
  if (empty($projects)){
    return $table;
  }
  foreach ($projects as $id=>$data){
    $table_rows[] = array($data['name'],$data['description']);
  }
  $table = theme_table($table_header, $table_rows, $table_attributes, $table_caption);
  return $table;
}

function gmod_dbsf_chado_link_project_contact_form($form_state){
  $form = array();
  $contacts = gmod_dbsf_db_execute('chado','gmod_dbsf_chado_get_contacts');
  $projects = gmod_dbsf_db_execute('chado','gmod_dbsf_chado_get_projects');

  if (empty($contacts)){
    drupal_set_message('No contacts have been defined!','warning',FALSE);
  }
  if (empty($projects)){
    drupal_set_message('No projects have been defined!','warning',FALSE);
  }
  if (empty($contacts) || empty($projects)){
    return $form;
  }
  // format arrays into options
  $contacts2 = array();
  $projects2 = array();
  foreach ($projects as $id=>$data){
    $projects2[$id] = $data['name'];
  }
  $projects = $projects2; unset($projects2);
  foreach ($contacts as $id=>$data){
    $contacts2[$data['type_name']][$id] = $data['name'];
  }
  $contacts = $contacts2; unset($contacts2);


  $form['link'] = array(
    '#type' => 'fieldset',
    '#tree'=>TRUE,
    '#title' => t('Associating contacts with projects'),
      'contact' => array(
      '#type'=>'select',
      '#title'=>t('Contact'),
      '#options' => $contacts,
      '#required'=>TRUE,
      '#size' => 1,
      '#weight'=>0,
  ),
    'project'=>array(
    	'#type'=>'select',
    	'#title'=>t('Project'),
    	'#options' => $projects,
    	'#required'=>TRUE,
  		'#size' => 1,
    	'#weight'=>1,
  ),
    'action'=>array(
    	'#type'=>'radios',
        '#options'=>array('LINK'=>'Link','UNLINK'=>'Unlink'),
    	'#title'=>t('Action'),
        '#default_value'=>'LINK',
    	'#weight'=>3,
  ),

  );
  //buttons
  $form['buttons']=array(
  	'#weight' => 15,
  	'submit' => array(
       '#type' => 'submit',
       '#value' => t('Submit data'),
       '#weight' => 1,
  ),
  	'clear' => array(
        '#attributes' => array('title' => t('Clear the form')),
        '#type' => 'submit',
        '#value' => 'Reset this form',
        '#submit' => array('gmod_dbsf_form_clear'),
        '#weight' => 2,
  )
  );
  return $form;
}

function gmod_dbsf_chado_link_project_contact_form_submit($form,&$form_state){
  $project = check_plain($form_state['values']['link']['project']);
  $contact = check_plain($form_state['values']['link']['contact']);
  $action = check_plain($form_state['values']['link']['action']);
  if (!empty($project) && !empty($contact)){
    gmod_dbsf_db_execute('chado_edit','gmod_dbsf_chado_link_unlink_project_contact',$project,$contact,$action);
  }
}

function gmod_dbsf_chado_add_project_data_form($form_state){
  $form = array();
  $projects = gmod_dbsf_db_execute('chado','gmod_dbsf_chado_get_projects');
  $project_options=array();
  if (empty($projects)){
    drupal_set_message('No projects have been defined!','warning',FALSE);
    return $form;
  }
  foreach ($projects as $id=>$data){
    $project_options[$id] = $data['name'];
  }
  $projects = $project_options; unset($project_options);
  $file_size_allowed = gmod_dbsf_get_add_var('gmod_dbsf_chado_project_data_upload_size',1);
  $validators_file   = array('file_validate_size' => array($file_size_allowed * 1048576));

  $file_upload_desc = t('Please give a text file, not a MS-Word or other document. It must have the first row with column names. You can upload up to %m Mb.', array('%m' => gmod_dbsf_get_add_var('biosoftware_bench_upload_size')));
    $form = array(
    'project_id'=>array(
      '#type'=>'select',
      '#title'=>t('Available projects'),
      '#options' => $projects,
      '#size' => 1,
      '#weight'=>-5,
  ),
  'data' => array(
  	'#type' => 'fieldset',
    '#title' => 'Add new data to project',
    '#tree' => TRUE,
    'user_file' => array(
        '#type' => 'file',
        '#title' => t('Upload delimited (CSV) file'),
        '#description' => $file_upload_desc,
        '#tree' => FALSE,
  ),
      'file_type' => array(
      	'#type' => 'radios',
      	'#title' => t('How is this file delimited?'),
      	'#options'=>array('comma'=>'Commas','tab'=>'Tabs'),
        '#required'=>TRUE,
        '#default_value'=>'comma',
  ),'quote_type' => array(
      	'#type' => 'radios',
      	'#title' => t('If any fields are quoted, are they quoted with double (") or single quotes (\')?'),
      	'#options'=>array('double'=>'Double quotes','single'=>'Single quotes'),
        '#required'=>TRUE,
        '#default_value'=>'double',
  ),
    'name_column' => array(
        		'#type' => 'textfield',
            	'#size' => 5,
        		'#title' => t('Optionally column number for names'),
      			'#description' => t('You can provide a column number for names (numbering starts from 1). Otherwise row number is used (starts from 1)'),
  ),
  	'long' => array(
        		'#type' => 'textfield',
            	'#size' => 5,
        		'#title' => t('Which column number is the longitude?'),
      			'#description' => t('Counting from 1'),
                '#required' =>TRUE,
  ),
    'lat' => array(
        		'#type' => 'textfield',
            	'#size' => 5,
        		'#title' => t('Which column number is the latitude?'),
            	'#description' => t('Counting from 1'),
  				'#required' =>TRUE,
  ),'habitat_col' => array(
        		'#type' => 'textfield',
            	'#size' => 5,
        		'#title' => t('Which column numbers relate to habitat variables you would like to store?'),
            	'#description' => t('Counting from 1. Use comma or space to separate. Dash (-) to define ranges.'),
  	),
  	'genotype_col' => array(
        		'#type' => 'textfield',
            	'#size' => 5,
        		'#title' => t('Which column numbers relate to genotypes (markers)?'),
            	'#description' => t('Counting from 1. Use comma or space to separate. Dash (-) to define ranges.'),
    )
  )
  );


    $form['buttons'] = array(
    '#weight' => 10,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Submit'),
  )
  );
  //required for files:
  $form['#attributes']['enctype'] = "multipart/form-data";
  return $form;

}

function gmod_dbsf_chado_add_project_data_form_validate($form,&$form_state){
  $form_state['values']['project_id'] = check_plain($form_state['values']['project_id']);
  $form_state['values']['data']['file_type'] = check_plain($form_state['values']['data']['file_type']);
  $form_state['values']['data']['name_column'] = trim(check_plain($form_state['values']['data']['name_column']));
  $form_state['values']['data']['habitat_col'] = trim(check_plain($form_state['values']['data']['habitat_col']));
  $form_state['values']['data']['genotype_col'] = trim(check_plain($form_state['values']['data']['genotype_col']));
  $form_state['values']['data']['long'] = trim(check_plain($form_state['values']['data']['long']));
  $form_state['values']['data']['lat'] = trim(check_plain($form_state['values']['data']['lat']));

  if (!isset($_FILES['files'])){
    form_set_error('user_file', t('Sorry but you must provide a file'));
  }
  if ($form_state['values']['data']['file_type'] != 'comma' && $form_state['values']['data']['file_type'] != 'tab'){
    form_set_error('data][file_type',t('Delimiter must be comma or tabs!'));
  }

  if (!is_numeric($form_state['values']['data']['lat'])){
    form_set_error('data][lat',t('Column number for Latitude must be a number!'));
  }
  if (!is_numeric($form_state['values']['data']['long'])){
    form_set_error('data][long',t('Column number for Longitude must be a number!'));
  }
  if (!is_numeric($form_state['values']['data']['name'])){
    form_set_error('data][name',t('Column number for Name must be a number!'));
  }
}

function gmod_dbsf_chado_add_project_data_form_submit($form,&$form_state){
  $operations = array();
  $file_size_allowed = gmod_dbsf_get_add_var('gmod_dbsf_chado_project_data_upload_size');
  $validators_file   = array('file_validate_size' => array($file_size_allowed * 1048576));
  $tmpdir            = file_directory_temp();
  $timestamp      = time();
  $sessionid      = session_id();
  $par = '';

  if (!empty($form_state['values']['data']['name'])){
    $par .= ' -experiment_uname_column '.$form_state['values']['data']['name'];
  }
  if (!empty($form_state['values']['data']['long'])){
    $par .= ' -experiment_longitude_column '.$form_state['values']['data']['long'];
  }
  if (!empty($form_state['values']['data']['lat'])){
    $par .= ' -experiment_latitude_column '.$form_state['values']['data']['lat'];
  }
  if (!empty($form_state['values']['data']['genotype_col'])){
    $t = str_replace(',', ' ', $form_state['values']['data']['genotype_col']);
    $t = preg_replace('/\s\s+/',' ',$t);
    preg_match_all('/(\d+)\-(\d+)/', $t,$matches,PREG_SET_ORDER );
    $i = 0;
    while (!empty($matches[$i])){
      if (!empty($matches[$i][0]) && !empty($matches[$i][1])){
        for ($k=$matches[$i][0];$k<=$matches[$i][1];$k++){
          $val[]=$k;
        }
      }
      $i++;
    }
    $t = preg_replace('/\s+\d+\-\d+\s+/',' ',$t);
    $val = array_merge((array)$val,explode(' ',$t));
    $par .= ' -experiment_genotype_columns '.implode(' ',$val);
  }
  if (!empty($form_state['values']['data']['habitat_col'])){
    $t = str_replace(',', ' ', $form_state['values']['data']['habitat_col']);
    $t = preg_replace('/\s\s+/',' ',$t);
    preg_match_all('/(\d+)\-(\d+)/', $t,$matches,PREG_SET_ORDER );
    $i = 0;
    while (!empty($matches[$i])){
      if (!empty($matches[$i][0]) && !empty($matches[$i][1])){
        for ($k=$matches[$i][0];$k<=$matches[$i][1];$k++){
          $val[]=$k;
        }
      }
      $i++;
    }
    $t = preg_replace('/\s+\d+\-\d+\s+/',' ',$t);
    $val = array_merge((array)$val,explode(' ',$t));
    $par .= ' -experiment_habitat_columns '.implode(' ',$val);
  }
  $form_uid = gmod_dbsf_create_uid($sessionid, $timestamp, 'gmod_dbsf_chado_project_data');
  if (isset($_FILES['files']) && !empty($_FILES['files']['tmp_name']['user_file'])) {
    $file = file_save_upload('user_file', $validators_file, $tmpdir, FILE_EXISTS_RENAME);
    if (empty($file)) {
      form_set_error('user_file', t('Sorry your file was not saved. Maybe it is too large (>%file_size_allowed Mb)? Otherwise, '. l('contact', 'contact') .' the administrator  (quote %uid).', array('%uid' => $form_uid, '%file_size_allowed' => $file_size_allowed)));
    }
    else {
      file_set_status($file, FILE_STATUS_TEMPORARY);
      $batch_file_data['uid'] = $form_uid;
      $batch_file_data['infile'] = $file->filepath;
      $batch_file_data['par'] = $par;
    }
  }else{
    form_set_error('user_file', t('Sorry but you must provide a file'));
  }
  $operations[] = array('gmod_dbsf_chado_add_project_data_upload_csv', array($batch_file_data));

  $batch        = array(
    'title' => t('Preparing data needed for adding to Chado...'),
    'operations' => $operations,
    'init_message' => t('Starting submission...'),
    'progress_message' => t('@remaining operations remaining...'),
    'error_message' => t('Your submission encountered an error.'),
    'finished' => 'gmod_dbsf_chado_add_project_data_batch_finished',
    'file' => drupal_get_path('module', 'gmod_dbsf_chado') .'/includes/gmod_dbsf_chado_projects.inc',
    'progressive'=> TRUE,
  );

  batch_set($batch);
  // We MUST do a batch_process, otherwise the form api attempts to load the entire
  // $form into $batch (as a variable and as a db) resulting in a 450 Mb variable
  batch_process();
}

function gmod_dbsf_chado_add_project_data_upload_csv($batch_file_data, &$context){
  if (empty($batch_file_data)){return;}
  $script_exec = drupal_get_path('module', 'geogeneticsWizard') .'/scripts/geotree.pl';
  $infile      = $batch_file_data['infile'];
  $par        = $batch_file_data['par'];
  if (empty($infile)) {
      $context['finished'] =1;
      return;
  }
  $par .= " -experiment_data_files '$infile' ";
  $output = array();
  $context['results']['cmd'] = $script_exec . $par;
  $context['results']['uid'] = $batch_file_data['uid'];
  $context['results']['json']= exec($script_exec . $par." 2>&1",$output);
  $context['results']['output'] = $output;
  file_delete($infile);
}

function gmod_dbsf_chado_add_project_data_batch_finished($success, $results, $operations){

}

/**
 * Get project or add a new project
 *
 * Here, we require a description for every project.
 *
 * @param $name
 * A project_id (integer) or unique name (text) in order to get an existing project or
 * a unique name (text) in order to add a new one
 * @param $description
 * A description for a new project (triggers adding a new project or updating an existing one)
 *
 * @return
 * Project ID if you're adding a new project or an array with 'id' and 'name'
 * if retrieving an existing project
 *
 */
function gmod_dbsf_chado_get_add_project($name,$description=NULL,$delete=FALSE){
  if (empty($name)){
    return;
  }
  $sql_select = "SELECT project_id as id,name from {project} WHERE ";
  $sql_insert = "INSERT INTO {project} (name,description) VALUES ('%s','%s')";
  $sql_update_description = "UPDATE {project} set description='%s' WHERE project_id=%d";
  $sql_delete = "DELETE from {project} WHERE project_id=%d";
  if (is_numeric($name)){
    $sql_select .= ' project_id=%d';
  }else{
    $sql_select .= " name='%s'";
  }

  if (empty($description)){
    //retrieve existing project
    $res = db_fetch_array(db_query($sql_select,$name));
    if (empty($res)){
      return FALSE;
    }
    if (empty($delete)){
      return $res;
    }else{
      db_query($sql_delete,$res['id']);
      $res = db_fetch_array(db_query($sql_select,$name));
      if (empty($res)){
        return TRUE;
      }else{
        return FALSE;
      }
    }
  }else{
    //adding new one unless it exists already (in that case it updates the description)
    //check if it exists, if yes get details and trigger update. otherwise if false then add it
    $check = gmod_dbsf_chado_get_add_project($name);
    if (!empty($check)){
      db_query($sql_update_description,$description,$check['id']);
    }else{
      db_query($sql_insert,$name,$description);
      //check if it was added ok
      $check = gmod_dbsf_chado_get_add_project($name);
      if (empty($check)){
        gmod_dbsf_drupal_set_message(t('I failed to add new project %name',array('%name'=>$name)),'error');
      }
    }
  }
}


function gmod_dbsf_chado_get_add_contact($name,$type=NULL,$description=NULL,$delete=FALSE){
  if (empty($name) || (!is_numeric($name) && empty($type)) ){return;}
  // type is required fo insert but not description
  $sql_select = "SELECT contact_id as id,name from {contact} WHERE ";

  $sql_check = $sql_select;
  $sql_insert = "INSERT INTO {contact} (name) VALUES ('%s')";
  $sql_update_term ='';

  $sql_update_description = "UPDATE {contact} set description='%s' WHERE contact_id=%d";
  $sql_delete = "DELETE from {contact} WHERE contact_id=%d";
  $check =array();
  if (is_numeric($name)){
    $sql_select .= " contact_id=%d ";
    $sql_check .= " contact_id=$name ";
  }else{
    $sql_select .= " name='%s' ";
    $sql_check  .= " name='$name' ";
  }

  $sql_recent_add = $sql_check .' AND type_id is NULL';

  if (!empty($type) && is_numeric($type)){
      $sql_update_term = "UPDATE {contact} set type_id=$type WHERE contact_id=%d";
      $sql_select .= ' AND type_id=%d ';
      $sql_check .=  " AND type_id=$type ";

  }elseif (!empty($type)){
    $sql_update_term = "UPDATE {contact} set type_id=(SELECT cvterm_id from {cvterm} where name='$type') WHERE contact_id=%d";
    $sql_select .= " AND type_id=(SELECT cvterm_id from {cvterm} where name='%s') ";
    $sql_check .= " AND type_id=(SELECT cvterm_id from {cvterm} where name='$type') ";
  }
  $check = db_fetch_array(db_query($sql_check));

  if (empty($check)){
    //does not exist
    if (!empty($delete)){
      return TRUE;
    }
    if (empty($type)){
      // we are not adding; report it does not exist
      return FALSE;
    }else{
      // we are adding
      db_query($sql_insert,$name);
      $check = db_fetch_array(db_query($sql_recent_add));
      if (empty($check)){
        gmod_dbsf_drupal_set_message(t('I failed to add new contact %name',array('%name'=>$name)),'error');
        return FALSE;
      }
      if (!empty($type)){
        db_query($sql_update_term,$check['id']);
      }
      if (!empty($description)){
        db_query($sql_update_description,$description,$check['id']);
      }
      gmod_dbsf_drupal_set_message(t('Added new contact %name',array('%name'=>$check['name'])),'warning');
      return $check;
    }
  }else{
    if (!empty($delete)){
      db_query($sql_delete,$check['id']);
      $check = db_fetch_array(db_query($sql_check));
      if (empty($check)){
        gmod_dbsf_drupal_set_message(t('Deleted contact %name',array('%name'=>$check['name'])),'warning');
        return TRUE;
      }else{
        return FALSE;
      }
    }else{
      return $check;
    }
  }
}

function gmod_dbsf_chado_get_contacts(){
  $sql_select = "SELECT contact_id as id,contact.name as contact,description,type_id,cvterm.name as term_name from {contact} JOIN {cvterm} ON type_id=cvterm_id";
  $res = db_query($sql_select);
  $return = array();
  while ($data = db_fetch_array($res)){
    $return[$data['id']]= array('name' => $data['contact'],'description'=>$data['description'],'type_name'=>$data['term_name'],'type_id'=>$data['type_id']);
  }
  return $return;
}

function gmod_dbsf_chado_get_project_contacts(){
  $sql_select = "SELECT project_contact_id as id,pc.project_id,pc.contact_id,p.name as project_name,c.name as contact_name from {project_contact} as pc "
  ." JOIN {project} as p ON pc.project_id=p.project_id JOIN {contact} c ON pc.contact_id=c.contact_id";
  $res = db_query($sql_select);
  $return = array();
  while ($data = db_fetch_array($res)){
    $return[$data['id']]= array('contact_id'=>$data['contact_id'],'project_id'=>$data['project_id'],'contact_name'=>$data['contact_name'],'project_name'=>$data['project_name']);
  }
  return $return;
}

function gmod_dbsf_chado_get_projects(){
  $sql_select = "SELECT project_id as id,name,description from {project}";
  $res = db_query($sql_select);
  $return = array();
  while ($data = db_fetch_array($res)){
    $return[$data['id']]= array('name' => $data['name'],'description'=>$data['description']);
  }
  return $return;
}

function gmod_dbsf_chado_link_unlink_project_contact($project,$contact,$action='LINK'){
 if (empty($project)||empty($contact)){return;}
 $sql_select = "SELECT project_contact_id as id from {project_contact} WHERE project_id=%d AND contact_id=%d ";
 $sql_link = "INSERT INTO {project_contact} (project_id,contact_id) VALUES (%d,%d) ";
 $sql_unlink = "DELETE FROM {project_contact} WHERE project_id=%d AND contact_id=%d";
 if (!is_numeric($contact)){
   $contact = gmod_dbsf_chado_get_add_contact($contact);
 }if (!is_numeric($project)){
   $project = gmod_dbsf_chado_get_add_contact($project);
 }
 $check = db_fetch_array(db_query($sql_select,$project,$contact));
 if (empty($project)||empty($contact)){return FALSE;}
  switch($action){
    case 'LINK':
      if (!empty($check)){
        gmod_dbsf_drupal_set_message(t('Contact and project were already linked'),'warning');
        return TRUE;
      }
      db_query($sql_link,$project,$contact);
      $check = db_fetch_array(db_query($sql_select,$project,$contact));
      if (!empty($check)){
        gmod_dbsf_drupal_set_message(t('Contact and project have been linked'),'warning');
        return TRUE;
      }
      break;
    case 'UNLINK':
      if (empty($check)){
        gmod_dbsf_drupal_set_message(t('Contact and project were not linked'),'warning');
        return TRUE;
      }
      db_query($sql_unlink,$project,$contact);
      $check = db_fetch_array(db_query($sql_select,$project,$contact));
      if (empty($check)){
        gmod_dbsf_drupal_set_message(t('Contact and project have been unlinked'),'warning');
        return TRUE;
      }
      break;
  }
  return FALSE;
}

