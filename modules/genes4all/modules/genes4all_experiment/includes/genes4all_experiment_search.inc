<?php
// $Id$

/**
 * @file
 * include file for searching for particular experiments
 */

function genes4all_experiment_page_search(){
  $return_text = "<p>This page allows you to search for a specific experiment.</p>";
  $return_text .= drupal_get_form('genes4all_experiment_page_search_form');
  return $return_text;
}

function genes4all_experiment_page_search_form($form_state){
  $search_type_array = array('name'=>'Unique name','id'=>'Internal numeric ID','abstract'=>'Abstract has word','species'=>'Target species','last_name'=>'Author last name');
  $search_form = array(
    'field'=>array(
      '#type' => 'fieldset',
      '#title' => 'Search criteria',
    'search_type' => array(
      '#title' => t('Select the type of search you wish to perform'),
      '#weight' => 0,
      '#type' => 'select',
      '#size' => 1,
      '#options' => $search_type_array,
      '#required' => TRUE,
  ),
    'text'=>array(
      '#title' => t('Search keyword/phrase'),
      '#weight' => 1,
      '#type' => 'textfield',
      '#size' => 15,
      '#required' => TRUE,
  ),
    'exact'=>array(
      '#title'=> t('Exact match'),
      '#weight' => 2,
      '#type' => 'checkbox',
      '#default_value' =>TRUE,
      '#description' => t('If checked search will be faster but the experiment must have the exact keyword/phrase (only that) stored for your search criterion. Mandatory for searching with a numerical identifier. Otherwise, any experiment containing the keyword/phrase will be found.'),
  ),
    'buttons' => array(
      '#type'=>'fieldset',
      '#title'=>'Actions',
      '#weight' => 5,
      'submit' => array(
      '#type' => 'submit',
      '#value' => t('Start'),
      '#weight' => 1,
  ),'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'button',
      '#value' => t('Reset'),
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 2,
  ),
  ),
  )
  );
  return $search_form ;
}

function genes4all_experiment_page_search_form_validate($form,&$form_state){
  if ($form_state['values']['search_type'] == 'id'){
    if (!is_numeric($form_state['values']['text'])){
      form_set_error('search_type','Identifiers must be numeric.');
    }elseif(empty($form_state['values']['exact'])){
      form_set_error('text','Searching usign an identifier use "exact".');
    }
  }
  return;
}

function genes4all_experiment_page_search_form_submit($form,&$form_state){
  $term = check_plain($form_state['values']['text']);
  $exact = $form_state['values']['exact'];
  $search_sql='';
  $search_results = array();
  switch ($form_state['values']['search_type']){
    case 'name':
      if (!empty($exact)){
        drupal_goto('genes4all/experiment/view/'.$answer);
      }else{
        $search_sql ="SELECT study_id as id from {gmod_dbsf_study} where uniquename~'$term'";
      }
      break;
    case 'id':
      drupal_goto('genes4all/experiment/view/'.$answer);
      break;
    case 'abstract':
      if (!empty($exact)){
        $search_sql ="SELECT study_id as id from {gmod_dbsf_study} where description='$term'";
      }else{
      $search_sql ="SELECT study_id as id from {gmod_dbsf_study} where description~'$term'";
      }
      break;
    case 'last_name':
      if (!empty($exact)){
        $search_sql ="SELECT distinct study_id as id from {gmod_dbsf_study} where pub_id IN "
        ."( SELECT pub_id from {gmod_dbsf_pub_author} pa JOIN {gmod_dbsf_author} a ON a.author_id=pa.author_id "
        ." WHERE last_names='$term')";
      }else{
      $search_sql ="SELECT distinct study_id as id from {gmod_dbsf_study} where pub_id IN "
       ."( SELECT pub_id from {gmod_dbsf_pub_author} pa JOIN {gmod_dbsf_author} a ON a.author_id=pa.author_id "
       ." WHERE last_names~'$term')";
      }
      break;
    case 'species':
      if (!empty($exact)){
        $search_sql ='SELECT distinct study_id as id from {gmod_dbsf_resource} where resource_id IN '
        .' (SELECT resource_id from {gmod_dbsf_resource} res JOIN {gmod_dbsf_organism} org ON res.organism_id=org.organism_id '
        ." WHERE genus='$term' OR species='$term' OR tax_class='$term' OR tax_order='$term' "
        ." OR tax_family='$term' OR common_name ='$term' OR abbreviation='$term' OR genus||' '||species='$term')";
      }else{
      $search_sql ='SELECT distinct study_id as id from {gmod_dbsf_resource} where resource_id IN '
        .' (SELECT resource_id from {gmod_dbsf_resource} res JOIN {gmod_dbsf_organism} org ON res.organism_id=org.organism_id '
        ." WHERE genus~'$term' OR species~'$term' OR tax_class~'$term' OR tax_order~'$term' OR tax_family~'$term' OR common_name ~'$term' OR abbreviation~'$term')";
      }
      break;
      form_set_error('Invalid search type');
      return;
  }
  $res = db_query($search_sql);
  while ($res_data = db_fetch_array($res)){
    $search_results[]=$res_data['id'];
  }
  if (empty($search_results)){
    form_set_error('',t('No experiments found with these search criteria'));
    return;
  }
  $q = implode($search_results,',');
    drupal_goto('genes4all/experiment/view/'.$q);
}