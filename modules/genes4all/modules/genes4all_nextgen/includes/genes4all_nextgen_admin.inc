<?php
// $Id$

/**
 * @file
 * Admin functions for genes4all_nextgen
 * Developed by Alexie Papanicolaou
 * @see http://insectacentral.org
 *
 * TODO
 * Actions
 * Allow more than one organism, feature, library?
 *
 */

function genes4all_nextgen_admin_page(){
  $rendered_admin_form = drupal_get_form('genes4all_nextgen_admin_new_form');
  $table_array=gmod_dbsf_db_execute('chado','_genes4all_nextgen_build_alignment_table');
  $rendered_admin_form .= drupal_get_form('genes4all_nextgen_admin_existing_form', $table_array);
  return $rendered_admin_form ;
}

function _genes4all_nextgen_build_alignment_table(){
  $table_array = array();
  $table_array['attributes'] = array('class' => 'padded-table');
  $table_array['caption'] = 'Existing Next Gen data in database.';
  $table_array['header'] = array(
    'Path', 'Description','Organism','Library','Feature','Date','Verified'
    );
    $show_ng_sql = "SELECT genes4all_nextgen_id as id,full_path,description,date(date_added) as date_added,curator_verified from {genes4all_nextgen} where obsolete is NULL OR obsolete = 0";
    $lib_ng_sql = "SELECT name,library_id as id from {library} where library_id IN (SELECT library_id from {genes4all_nextgen_library} where genes4all_nextgen_id=%d)";
    $org_ng_sql = "SELECT genus,species,organism_id as id from {organism} where organism_id IN (SELECT organism_id from {genes4all_nextgen_organism} where genes4all_nextgen_id=%d)";
    $feat_ng_sql = "SELECT name,feature_id as id from {feature} where feature_id IN (SELECT feature_id from {genes4all_nextgen_feature} where genes4all_nextgen_id=%d)";
    $ng_res = db_query($show_ng_sql);
    while ($ng_data = db_fetch_array($ng_res)){
      $organisms=array();
      $libraries=array();
      $features=array();
      $res = db_query($org_ng_sql,$ng_data['id']);
      while ($data = db_fetch_array($res)){
        $organisms[]=$data['genus'].' '.$data['species'];
      }
      $res = db_query($lib_ng_sql,$ng_data['id']);
      while ($data = db_fetch_array($res)){
        $libraries[]=$data['name'];
      }
      $res = db_query($feat_ng_sql,$ng_data['id']);
      while ($data = db_fetch_array($res)){
        $features[]=$data['name'];
      }
      if (empty($organisms)){
        $organisms[]='N/A';
      }if (empty($libraries)){
        $libraries[]='N/A';
      }if (empty($features)){
        $features[]='N/A';
      }
      $verified = empty($ng_data['curator_verified']) ? 'No' : 'Yes';
      $table_array['data'][]= array(
      array('data' => basename($ng_data['full_path']),'id'=>$ng_data['id']),
      array('data' => $ng_data['description']),
      array('data' => (count($organisms)>2) ? count($organisms) : implode(', ',$organisms)),
      array('data' => (count($libraries)>2) ? count($libraries) : implode(', ',$libraries)),
      array('data' => (count($features)>2) ? count($features) : implode(', ',$features)),
      array('data' => $ng_data['date_added']),
      array('data' => $verified),
      );
    }
    return $table_array;
}

function genes4all_nextgen_admin_new_form($form_state){
  $library_names=array(0=>'None');
  $library_data=gmod_dbsf_db_execute('chado','gmod_dbsf_get_add_library','select_options');
  $organism_data = gmod_dbsf_db_execute('chado','gmod_dbsf_getspecies',NULL,1);
  $organism_names = array(0=>'None')+(array)$organism_data['species_array'];
  $form = array(
    '#type' => 'fieldset',
    '#title' => t('Next Generation datasets'),
    '#tree' => TRUE,
    '#description' => t('Here you can add new NGS datasets if they are formatted as !s or BAM files. Each NGS dataset can belong to a feature, an organism, a library or a combination. Ancillary software such as genes4all and JBrowse will be able to utilise these datasets and display the relevant information to the user.',
  array('!s'=>l('SAM','http://samtools.sourceforge.net/',array('attributes'=>array('target'=>'_blank'),'absolute'=>TRUE)))),
    'server_path' => array(
      '#type' => 'textfield',
      '#title' => t('Full path on server'),
      '#description' => t('If the file resides on the server already, please give the full path.'),
  ),
    'upload_file'=>array(
      '#type' => 'file',
      '#title'=> t('Upload file'),
      '#description'=> t('Alternatively, you may upload a file to the server.')
  ),
    'description' => array(
      '#type' => 'textfield',
      '#title' => t('Description'),
      '#description' => t('Please give an informative description to this dataset.'),
      '#required'=> TRUE,
  ),
    'linked_to' => array(
      '#type'=>'fieldset',
      '#title' => t('This NGS will be linked to one or more of these database tables.'),
      '#description'=>t('At least one these table types must be used. NB: Unless you use the SAM header, the file will not be checked to see if the reference sequence exists in the database.'),
      '#tree' => TRUE,
        'grab_refs'=>array(
          '#type'=>'checkbox',
          '#title'=>'Get references from SAM header',
          '#description'=>'Check this to acquire which sequences are linked to this file using the alignment index (accurate but slower)',
          '#default_value'=>1,
  ),
  //TODO allow multiple orgs, libs and features (CSV)
      'organism'=>array(
        '#type' => 'select',
        '#title' => t('Organism'),
        '#description' => t('Choose the organism binomial species name if you wish to link to a species.'),
        '#options'=>$organism_names,
  ),
      'library'=>array(
        '#type' => 'select',
        '#title' => t('Library'),
        '#description' => t('Choose the library name if you wish to link to a library.'),
        '#options'=>$library_names,
  ),
      'feature'=>array(
        '#type' => 'textfield',
        '#title' => t('Feature'),
        '#description' => t('Please give the unique name of a feature if you wish to link to only a single feature.'),
  ),
  ),
  );
  $form['buttons'] = array(
    '#type' => 'submit',
    '#value' => t('Submit'),
  );
  $form['#attributes']['enctype'] = "multipart/form-data";
  return $form;
}

function genes4all_nextgen_admin_new_form_validate($form,&$form_state){
  $save_sam_dir = file_create_path() .'/nextgen';
  $file_size_allowed = 100000000; //gmod_dbsf_get_add_var('genes4all_nextgen_upload_size');
  $validators_file   = array('file_validate_size' => array($file_size_allowed * 1024));
  if (empty($_FILES['files']['tmp_name']['upload_file']) && empty($form_state['values']['server_path'])){
    form_set_error('server_path','You must give a file as a full path or an upload.');
    return FALSE;
  }elseif(!empty($form_state['values']['server_path'])){
    $form_state['values']['server_path']=check_plain($form_state['values']['server_path']);
    if (!file_exists($form_state['values']['server_path']) ){
      form_set_error('server_path',t('File %s was not found on the server.',array('%s'=>$form_state['values']['server_path'])));
      return FALSE;
    }else{
      //file test if it is a sam or bam?
      $file_type=genes4all_nextgen_check_sam_file_type($form_state['values']['server_path']);
      if (empty($file_type)){
        form_set_error('server_path',t('File %s does not seem to be a SAM or BAM file or is not readable by server.',array('%s'=>$form_state['values']['server_path'])));
        return FALSE;
      }else{
        $form_state['values']['file_type']=$file_type;
      }
    }
  }elseif(!empty($_FILES['files']['tmp_name']['upload_file'])){
    if (!isset($_FILES['files']['tmp_name']['upload_file']) ){
      form_set_error('upload_file',t('File %s was not uploaded successfully.',array('%s'=>$form_state['values']['upload_file'])));
      return FALSE;
    }else {
      //file test if it is a sam or bam?
      $file_type=genes4all_nextgen_check_sam_file_type($form_state['values']['server_path']);
      if (empty($file_type)){
        form_set_error('server_path',t('File %s does not seem to be a SAM or BAM file.',array('%s'=>$form_state['values']['server_path'])));
        return FALSE;
      }else{
        $form_state['values']['file_type']=$file_type;
        $file = file_save_upload('upload_file', $validators_file, $save_sam_dir, FILE_EXISTS_RENAME);
        if (empty($file)) {
          form_set_error('upload_file', t('Sorry your file was not saved. Maybe it is too large (>%file_size_allowed Mb)?', array('%file_size_allowed' => $file_size_allowed)));
          return FALSE;
        }else{
          file_set_status($file, FILE_STATUS_PERMANENT);
          $form_state['values']['server_path']=$file->filepath;
        }
      }
    }
  }
  if (empty($form_state['values']['linked_to']['grab_refs'])  && empty($form_state['values']['linked_to']['feature'])
  && empty($form_state['values']['linked_to']['library']) && empty($form_state['values']['linked_to']['organism'])){
    form_set_error('linked_to',t('You must provide at least one table to link with: organism, library or a feature'));
    return FALSE;
  }
  if (!empty($form_state['values']['linked_to']['feature'])){
    $form_state['values']['linked_to']['feature'] = check_plain($form_state['values']['linked_to']['feature']);
    $feature_id_array = gmod_dbsf_db_execute('chado','gmod_dbsf_get_feature_name_id',$form_state['values']['linked_to']['feature']);
    if (empty($feature_id_array)){
      form_set_error('feature',t('Sorry, feature %s was not found.',array('%s'=>$form_state['values']['linked_to']['feature'])));
      return FALSE;
    }
    $form_state['values']['linked_to']['feature'] =$feature_id_array;
  }if (!empty($form_state['values']['linked_to']['organism'])){
    $form_state['values']['linked_to']['organism'] = check_plain($form_state['values']['linked_to']['organism']);
    $organism_id_array = gmod_dbsf_db_execute('chado','gmod_dbsf_decypher_organism_id',$form_state['values']['linked_to']['organism']);
    if (empty($organism_id_array)){
      form_set_error('organism',t('Sorry, organism %s was not found.',array('%s'=>$form_state['values']['linked_to']['organism'])));
      return FALSE;
    }
    $form_state['values']['linked_to']['organism'] = $organism_id_array;
  }if (!empty($form_state['values']['linked_to']['library'])){
    $form_state['values']['linked_to']['library'] = check_plain($form_state['values']['linked_to']['library']);
    $library_id_array = gmod_dbsf_db_execute('chado','gmod_dbsf_decypher_library_id',$form_state['values']['linked_to']['library']);
    if (empty($library_id_array)){
      form_set_error('library',t('Sorry, library %s was not found.',array('%s'=>$form_state['values']['linked_to']['library'])));
      return FALSE;
    }
    $form_state['values']['linked_to']['library'] = $library_id_array;
  }if (!empty($form_state['values']['linked_to']['grab_refs'])){
    $form_state['values']['linked_to']['grab_refs'] = check_plain($form_state['values']['linked_to']['grab_refs']);
  }
  $form_state['values']['description']=check_plain($form_state['values']['description']);
}

function genes4all_nextgen_admin_new_form_submit($form,&$form_state){
  $operations = array();
  $samtools_exec = gmod_dbsf_drupal_root(). drupal_get_path('module', 'genes4all_nextgen') . '/bin/samtools';

  if (!file_exists($samtools_exec) || !is_executable($samtools_exec)){
    drupal_set_message(t('There is a problem with the samtools executable. Please contact your sysadmin.'),'error');
    return FALSE;
  }
  $file_path = $form_state['values']['server_path'];
  $file_type = $form_state['values']['file_type'];
  $description = $form_state['values']['description'];
  // TODO facilitate multiple organisms, features, libraries
  $feature_id_array = $form_state['values']['linked_to']['feature'];
  $library_id = $form_state['values']['linked_to']['library']['library_id'];
  $organism_id = $form_state['values']['linked_to']['organism']['organism_id'];
  $grab_refs = $form_state['values']['linked_to']['grab_refs'];
  if(empty($file_path) || !file_exists($file_path)){
    drupal_set_message('Something went terribly wrong with trying to save the file...','error');
    return FALSE;
  }
  $operations = array();
  $operations[]=array('_genes4all_nextgen_process_alignment_file',array($file_path,$file_type,$description,$feature_id_array,$library_id,$organism_id,$grab_refs));
  $operations[]=array('_genes4all_nextgen_create_bam_alignment',array());
  $operations[]=array('_genes4all_nextgen_process_alignment_file_sort',array());
  $operations[]=array('_genes4all_nextgen_process_alignment_file_index',array());
  batch_set(array(
    'title' => t('Uploading your alignments into the database...'),
    'operations' => $operations,
  	'file' => drupal_get_path('module', 'genes4all_nextgen') .'/includes/genes4all_nextgen_admin.inc',
    'init_message' => t('Uploading your alignments into the database...'),
    'progress_message' => t('@remaining operations remaining...'),
    'error_message' => t('Your submission encountered an error.'),
    'finished' => 'genes4all_nextgen_batch_finished',
  ));
  batch_process();
}

function genes4all_nextgen_batch_finished($success, $results, $operations) {
  $real_success = $results['real_success'];
  if (empty($success) || empty($real_success)){
    $error_operation = reset($operations);
    drupal_set_message(t('An error occurred while processing your request. Please '
    . l('contact', 'contact') ." the administrator."),'error');
    return FALSE;
  }
  $file_path = $results['file_path'];
  $description= $results['description'];
  $feature_id_array=$results['feature_id_array'];
  $library_id=$results['library_id'];
  $organism_id=$results['organism_id'];
  $grab_references =$results['grab_refs'];
  // Now we have an indexed BAM file as $file_path. Time to store it in the database
  $next_gen_id = gmod_dbsf_db_execute('chado_edit','_genes4all_nextgen_store_alignment',$file_path,$description);
  if (empty($next_gen_id)){
    drupal_set_message(t('Failed to link your data with the database.'),'error');
    return FALSE;
  }else{
    drupal_set_message(t('Linked %s with the database.',array('%s'=>$file_path)),'warning');
  }
  // depending on the type of metadata submitted (library_id, organism_id, feature_id), store in those tables
  if (!empty($organism_id)){
    $message = gmod_dbsf_db_execute('chado_edit','_genes4all_nextgen_link_alignment','organism',$next_gen_id,$organism_id);
    drupal_set_message($message,'warning');
  }
  if (!empty($library_id)){
    $message =gmod_dbsf_db_execute('chado_edit','_genes4all_nextgen_link_alignment','library',$next_gen_id,$library_id);
    drupal_set_message($message,'warning');
  }
  if (!empty($feature_id_array)){
    foreach ($feature_id_array as $feature_id=>$feature_name){
      $message = gmod_dbsf_db_execute('chado_edit','_genes4all_nextgen_link_alignment','feature',$next_gen_id,$feature_id);
      drupal_set_message($message,'warning');
    }
  }if (!empty($grab_references)){
    $features = genes4all_nextgen_alignment_references($file_path);
    if (empty($features)){
      drupal_set_message(t('Failed to link your data with the database.'),'error');
      return FALSE;
    }
    foreach ($features as $feature_name){
      $message = gmod_dbsf_db_execute('chado_edit','_genes4all_nextgen_link_alignment','feature',$next_gen_id,$feature_name);
    }
  }
  return TRUE;
}

function _genes4all_nextgen_process_alignment_file_index(&$context){
  $samtools_exec = $context['results']['samtools'];
  $file_path = $context['results']['file_path'];
  if (empty($samtools_exec) || empty($file_path)){return;}
  // Index BAM file
  $exec = exec("$samtools_exec index $file_path");
  if (strpos($exec,'fail')!==FALSE){
    drupal_set_message($exec,'error');
    _batch_finished();
    return FALSE;
  }
  drupal_set_message($exec);
  $context['results']['file_path']=gmod_dbsf_drupal_root().$context['results']['file_path'];
  $context['results']['real_success'] = 1;
}

function _genes4all_nextgen_process_alignment_file_sort(&$context){
  $samtools_exec = $context['results']['samtools'];
  $file_path = $context['results']['file_path'];
  if (empty($samtools_exec) || empty($file_path)){return;}
  // Sort BAM file
  $exec = exec("$samtools_exec sort $file_path $file_path");
  if (strpos($exec,'fault')!==FALSE){
    drupal_set_message($exec,'error');
    _batch_finished();
    return FALSE;
  }
  file_delete($file_path);
  rename("$file_path.bam",$file_path);
  $context['message'] = t('Indexing BAM file.');
}

function _genes4all_nextgen_create_bam_alignment(&$context){
  $samtools_exec = gmod_dbsf_drupal_root(). drupal_get_path('module', 'genes4all_nextgen') . '/bin/samtools';
  $context['results']['samtools']=$samtools_exec ;
  $file_type = $context['results']['file_type'];
  $file_path = $context['results']['file_path'];
  if (empty($samtools_exec) || empty($file_path) || empty($file_type)){return;}
  switch ($file_type){
    case 'SAM':
      $exec = exec("$samtools_exec view -F 0x0004 -b -S -o $file_path.bam $file_path");
      if (strpos($exec,'Abort')!==FALSE){
        drupal_set_message($exec,'error');
        _batch_finished();
        return FALSE;
      }
      file_delete($file_path);
      $file_path = $file_path.'.bam';
      break;
    case 'BAM':
      // remove unaligned
      $exec = exec("$samtools_exec view -F 0x0004 -b -o $file_path.aligned $file_path");
      file_delete($file_path);
      rename("$file_path.aligned",$file_path);
      break;
  }
  $context['results']['file_path']=$file_path ;
  $context['message'] = t('Sorting BAM file.');
}
function _genes4all_nextgen_process_alignment_file($file_path,$file_type,$description=NULL,$feature_id_array=NULL,$library_id=NULL,$organism_id=NULL,$grab_refs=NULL,&$context){
  // this function needs to  be broken up at each 'message' step into separate functions in order to allow processing feedback to the user
  $save_sam_dir = file_create_path() .'/nextgen/';
  switch ($file_type){
    case 'SAM':
      // Convert to BAM
      file_copy($file_path,$save_sam_dir,'FILE_EXISTS_RENAME');
      if (strrpos($file_path,'.gz')===(strlen($file_path)-3)){
        $context['message'] = t('Processing gzipped SAM file to BAM.');
        $file_path = preg_replace('/.gz$/', '', $file_path);
        exec("gunzip -c $file_path.gz > $file_path");
        file_delete("$file_path.gz");
      }elseif (strrpos($file_path,'.bz2')===(strlen($file_path)-4)){
        $context['message'] = t('Processing bzipped SAM file to BAM.');
        $file_path = preg_replace('/.bz2$/', '', $file_path);
        exec("bunzip2 -k $file_path.bz2");
        file_delete("$file_path.bz2");
      }else{
        $context['message'] = t('Processing uncompressed SAM file to BAM.');
      }
      break;
    case 'BAM':
      file_copy($file_path,$save_sam_dir,'FILE_EXISTS_RENAME');
      if (strrpos($file_path,'.gz')===(strlen($file_path)-3)){
        $context['message'] = t('Processing gzipped BAM file.');
        $file_path = preg_replace('/.gz$/', '', $file_path);
        exec("gunzip -c $file_path.gz > $file_path");
        file_delete("$file_path.gz");
      }elseif (strrpos($file_path,'.bz2')===(strlen($file_path)-4)){
        $context['message'] = t('Processing bzipped BAM file.');
        $file_path = preg_replace('/.bz2$/', '', $file_path);
        exec("bunzip2 -k $file_path.bz2");
        file_delete("$file_path.bz2");
      }else{
        $context['message'] = t('Processing BAM file.');
      }
  }
  //sandbox does not work between operations Drupal 6.21-dev
  $context['results']['file_type'] = $file_type;
  $context['results']['file_path']=$file_path;
  $context['results']['description']=$description;
  $context['results']['feature_id_array']=$feature_id_array;
  $context['results']['library_id']=$library_id;
  $context['results']['organism_id']=$organism_id;
  $context['results']['grab_refs'] = $grab_refs;
}

function genes4all_nextgen_admin_existing_form($form_state, $table_array){
  if (empty($table_array['data'])) {
    drupal_set_message('No Next Gen data have been previously linked to your database.','warning');
    return;
  }
  $row_check_data = array();
  $form = array();
  foreach ($table_array['data'] as $row => $value_array) {
    $path = $value_array[0]['data'];
    $path_id = $value_array[0]['id'];
    $row_check_data[$path_id] = wordwrap($path,36,"<br>",TRUE);
  }
  $form             = array();
  $form['data']     = array('#value' => $table_array);
  $form['checkboxes'] = array(
    '#type' => 'checkboxes',
    '#options' => $row_check_data,
  );
  $form['buttons'] = array(
      'delete' => array(
        '#type' => 'submit',
        '#value' => t('Delete'),
        '#weight' => 15,
  ),'curator_validate' => array(
        '#type' => 'submit',
        '#value' => t('Verify'),
        '#weight' => 10,
  ),
  );
  return $form;
}

function genes4all_nextgen_admin_existing_form_submit($form,&$form_state){
  $op = $form_state['values']['op'];
  switch ($op){
    case 'Delete':
      foreach ($form_state['values']['checkboxes'] as $path_id=>$value){
        if (empty($value)){
          continue;
        }
        $res = gmod_dbsf_db_execute('chado_edit','_genes4all_nextgen_delete_alignment',$path_id);
        if ($res[0] ===FALSE){
          drupal_set_message(t('There has been a problem with deleting %s.',array('%s'=>$res[1])),'error');
        }else{
          drupal_set_message(t('Alignment file %s successfully removed from the database.',array('%s'=>$res[1])),'warning');
        }
      }
      break;
    case 'Verify':
      foreach ($form_state['values']['checkboxes'] as $path_id=>$value){
        if (empty($value)){
          continue;
        }
        $res = gmod_dbsf_db_execute('chado_edit','_genes4all_nextgen_validate_alignment',$path_id);
        if ($res[0] ==FALSE){
          drupal_set_message(t('Curator verification flag has now been removed from alignment file %s.',array('%s'=>$res[1])),'warning');
        }else{
          drupal_set_message(t('Alignment %s has now been flagged as verified by a curator.',array('%s'=>$res[1])),'warning');
        }
      }
      break;
  }
}

function genes4all_nextgen_alignment_references($filename=NULL){
  if (empty($filename) || !file_exists($filename)){return FALSE;}
  $reference_sequences = array();
  $samtools_exec = gmod_dbsf_drupal_root(). drupal_get_path('module', 'genes4all_nextgen') . '/bin/samtools';
  $sam_args = escapeshellcmd(" idxstats $filename");
  $exec_value='';
  $ret = exec($samtools_exec.$sam_args."|egrep -o '^\S+'",$reference_sequences,$exec_value);
  if (strpos($ret,' fail '!==FALSE)){
    drupal_set_message(t('There was a problem with getting headers from %f.',array('%f'=>$filename)),'error');
    return FALSE;
  }
  //remove '*'
  array_pop($reference_sequences);
  return $reference_sequences;
}

function _genes4all_nextgen_store_alignment($file_path=NULL,$description=NULL){
  if (empty($file_path)||empty($description)){
    return FALSE;
  }

  $next_gen_id = '';
  $check_sql = "SELECT genes4all_nextgen_id as id from {genes4all_nextgen} where full_path = '$file_path'";
  $insert_sql = "INSERT INTO {genes4all_nextgen} (full_path,description) VALUES ('$file_path','$description')";
  $exists  = db_fetch_array(db_query($check_sql));
  if (empty($exists['id'])){
    db_query($insert_sql);
    $exists  = db_fetch_array(db_query($check_sql));
    if (empty($exists['id'])){
      return;
    }else{
      $next_gen_id = $exists['id'];
    }
  }else{
    $next_gen_id = $exists['id'];
  }
  return $next_gen_id ;
}

function _genes4all_nextgen_link_alignment($type=NULL,$next_gen_id=NULL,$type_id=NULL){
  if (empty($type)||empty($next_gen_id)||empty($type_id)){
    return FALSE;
  }
  if (!is_numeric($type_id)) {
    $sql_select = 'SELECT '.$type.'_id as id from {'.$type."} WHERE uniquename='$type_id'";
    $res = db_fetch_array(db_query($sql_select,$type_id));
    if (empty($res)){
      return t('Failed to link your data with new %a. Couldn\'t find its ID',array('%a'=>$type));
    }
    $type_id = $res['id'];
  }
  $check_sql = 'SELECT genes4all_nextgen_'.$type.'_id as id from {genes4all_nextgen_'.$type."} where genes4all_nextgen_id=$next_gen_id and ".$type.'_id = %d ';
  $insert_sql = 'INSERT INTO {genes4all_nextgen_'.$type.'} ('.$type."_id,genes4all_nextgen_id) VALUES (%d,$next_gen_id)";
  $exists  = db_fetch_array(db_query($check_sql,$type_id));
  if (empty($exists['id'])){
    db_query($insert_sql,$type_id);
    $exists  = db_fetch_array(db_query($check_sql,$type_id));
    if (empty($exists['id'])){
      return t('Failed to link your data with new %a.',array('%a'=>$type));
    }else{
      return t('Linked with new %a.',array('%a'=>$type));
    }
  }
  return t('File was already linked with %a.',array('%a'=>$type));
}

function _genes4all_nextgen_delete_alignment($id){
  if (empty($id)){
    return;
  }
  $check_sql = "SELECT genes4all_nextgen_id as id,full_path as path from {genes4all_nextgen} where genes4all_nextgen_id=$id";
  $delete_sql ="DELETE from {genes4all_nextgen} where genes4all_nextgen_id=$id";
  $data  = db_fetch_array(db_query($check_sql));
  db_query($delete_sql);
  $exists  = db_fetch_array(db_query($check_sql));
  if (empty($exists['id'])){
    file_delete($data['path']);
    file_delete($data['path'].'.bai');
    return array(TRUE,$data['path']);
  }else{
    return array(FALSE,$data['path']);
  }
}

function _genes4all_nextgen_validate_alignment($id){
  if (empty($id)){
    return;
  }
  $validate_sql = "UPDATE {genes4all_nextgen} SET curator_verified = 1 WHERE genes4all_nextgen_id=$id";
  $devalidate_sql = "UPDATE {genes4all_nextgen} SET curator_verified = 0 WHERE genes4all_nextgen_id=$id";
  $check_sql = "SELECT full_path as path,curator_verified as switch from {genes4all_nextgen} where genes4all_nextgen_id=$id";
  $data  = db_fetch_array(db_query($check_sql));
  if (empty($data['switch'])){
    db_query($validate_sql);
    $data  = db_fetch_array(db_query($check_sql));
    return array($data['switch'],$data['path']);
  }else{
    db_query($devalidate_sql);
    $data  = db_fetch_array(db_query($check_sql));
    return array($data['switch'],$data['path']);
  }
}
