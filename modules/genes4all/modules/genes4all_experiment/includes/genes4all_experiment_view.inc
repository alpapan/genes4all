<?php
// $Id$

/**
 * @file
 * include file for viewing experiments
 */

/**
 * View all experiments
 *
 * @param $round
 * Integer for page to view. Each page has $how_many experiments
 *
 * @return unknown_type
 */
function genes4all_experiment_page_view_experiment_all($round = 0) {
  $dirpath = drupal_get_path('module', 'genes4all');
  $next_link = l('<img src="/'. $dirpath .'/images/next-button.png">', 'genes4all/experiment/ueberview/'. ($round + 1),
    array(
      'html' => TRUE,
      'alias' => TRUE,
      'attributes' => array('title' => 'Get next experiments'),
    )
  );
  if (!empty($round)) {
    $previous_link = l('<img src="/'. $dirpath .'/images/previous-button.png">', 'genes4all/experiment/ueberview/'. ($round - 1),
      array(
        'html' => TRUE,
        'alias' => TRUE,
        'attributes' => array('title' => 'Get previous experiments'),
      )
    );
  }
  $how_many=5;

  $return_text = '<p>' . 'Overview of all experiments submitted.<br>'. $previous_link;

  $study_select    = "SELECT study_id as id from {gmod_dbsf_study}";
  $study_res       = db_query_range($study_select, $round * $how_many, ($round * $how_many) + $how_many);
  $study_res_ahead = db_query_range($study_select, $round + 1 * $how_many, ($round + 1 * $how_many) + $how_many);
  $study_ahead     = db_fetch_array($study_res_ahead);
  if (!empty($study_ahead)) {
    $return_text .= $next_link;
  }
  $return_text .= '</p>';
  while ($row = db_fetch_array($study_res)) {
    $experiment_id   = $row['id'];
    $pub_select      = "SELECT pub_id as id from {gmod_dbsf_study} where study_id=%d";
    $resource_select = "select sr.resource_id as id from  {gmod_dbsf_study_resource} sr join {gmod_dbsf_resource} r ON " . " sr.resource_id=r.resource_id join {gmod_dbsf_cvterm} cvterm ON type_id=cvterm_id  where cvterm.name='%s' AND  study_id=%d";
    $feature_select  = "select sr.feature_id as id from  {gmod_dbsf_study_feature} sr join {gmod_dbsf_feature} r ON sr.feature_id=r.feature_id" . " join {gmod_dbsf_cvterm} cvterm ON type_id=cvterm_id where cvterm.name='%s' AND study_id=%d";

    $res = db_fetch_array(db_query($pub_select, $experiment_id));
    $pub_id = $res['id'];

    $res = db_fetch_array(db_query($feature_select, 'gene_target', $experiment_id));
    $target_id = $res['id'];

    $res = db_fetch_array(db_query($feature_select, 'rnai_construct', $experiment_id));
    $construct_id = $res['id'];

    $res = db_fetch_array(db_query($resource_select, 'experimental animals', $experiment_id));
    $animals_id = $res['id'];

    $res = db_fetch_array(db_query($resource_select, 'delivery protocol', $experiment_id));
    $delivery_id = $res['id'];

    $res = db_fetch_array(db_query($resource_select, 'assay protocol', $experiment_id));
    $assay_id = $res['id'];

    $table_study     = genes4all_experiment_tabulate_study_data($experiment_id, 'study_id', TRUE);
    $table_pub       = genes4all_experiment_tabulate_study_data($pub_id, 'pub_id', TRUE);
    $table_target    = genes4all_experiment_tabulate_study_data($target_id, 'target_id', TRUE);
    $table_construct = genes4all_experiment_tabulate_study_data($construct_id, 'construct_id', TRUE);
    $table_animal    = genes4all_experiment_tabulate_study_data($animals_id, 'animals_id', TRUE);
    $table_delivery  = genes4all_experiment_tabulate_study_data($delivery_id, 'delivery_id', TRUE);
    $table_assay     = genes4all_experiment_tabulate_study_data($assay_id, 'assay_id', TRUE);
    $return_text .= $table_study . $table_pub . $table_target . $table_construct . $table_animal . $table_delivery . $table_assay;
  }
  return $return_text;
}

/**
 * Helper function to get experiment details
 *
 * @param $experiment_id
 * Can be a numerical identifier or uniquename
 *
 * @return Tabulated HTML
 */
function _genes4all_experiment_experiment_details($experiment_id){
  $return_text = '';
  if (empty($experiment_id)) {
    return;
  }
  if (!is_numeric($experiment_id)) {
    $study_select_sql = "SELECT study_id from {gmod_dbsf_study} where uniquename='%s'";
    $experiment_id    = db_fetch_array(db_query($study_select_sql, $experiment_id));
    $experiment_id    = $experiment_id['study_id'];
  }
  if (!empty($experiment_id) && is_numeric($experiment_id)) {
    $study_select_sql = "SELECT name from {gmod_dbsf_cvterm} where " . "cvterm_id=(SELECT type_id from {gmod_dbsf_study} where study_id=%d)";
    $access           = db_fetch_array(db_query($study_select_sql, $experiment_id));
    $access           = $access['name'];
    if ($access == 'unfinished' && user_access('View all experiments') !== TRUE) {
      drupal_set_message(t('Sorry, but the submission of this study is not finalized and therefore we cannot show it to the public.'), 'error');
      return FALSE;
    }
    $pub_select      = "SELECT pub_id as id from {gmod_dbsf_study} where study_id=%d";
    $resource_select = "select sr.resource_id as id from {gmod_dbsf_study_resource} sr join {gmod_dbsf_resource} r ON " . " sr.resource_id=r.resource_id join {gmod_dbsf_cvterm} cvterm ON type_id=cvterm_id  where cvterm.name='%s' AND study_id=%d";
    $feature_select  = "select sr.feature_id as id from  {gmod_dbsf_study_feature} sr join {gmod_dbsf_feature} r ON sr.feature_id=r.feature_id" . " join {gmod_dbsf_cvterm} cvterm ON type_id=cvterm_id  where cvterm.name='%s' AND study_id=%d";

    $res = db_fetch_array(db_query($pub_select, $experiment_id));
    $pub_id = $res['id'];

    $res = db_fetch_array(db_query($feature_select, 'gene_target', $experiment_id));
    $target_id = $res['id'];

    $res = db_fetch_array(db_query($feature_select, 'rnai_construct', $experiment_id));
    $construct_id = $res['id'];

    $res = db_fetch_array(db_query($resource_select, 'experimental animals', $experiment_id));
    $animals_id = $res['id'];

    $res = db_fetch_array(db_query($resource_select, 'delivery protocol', $experiment_id));
    $delivery_id = $res['id'];

    $res = db_fetch_array(db_query($resource_select, 'assay protocol', $experiment_id));
    $assay_id = $res['id'];

    $table_study     = genes4all_experiment_tabulate_study_data($experiment_id, 'study_id', TRUE);
    if (empty($table_study)){
      return t('Seems you provided an incorrect identifier');
    }
    $table_pub       = genes4all_experiment_tabulate_study_data($pub_id, 'pub_id', TRUE);
    $table_target    = genes4all_experiment_tabulate_study_data($target_id, 'target_id', TRUE);
    $table_construct = genes4all_experiment_tabulate_study_data($construct_id, 'construct_id', TRUE);
    $table_animal    = genes4all_experiment_tabulate_study_data($animals_id, 'animals_id', TRUE);
    $table_delivery  = genes4all_experiment_tabulate_study_data($delivery_id, 'delivery_id', TRUE);
    $table_assay     = genes4all_experiment_tabulate_study_data($assay_id, 'assay_id', TRUE);
    $return_text = $table_study . $table_pub . $table_target . $table_construct . $table_animal . $table_delivery . $table_assay;
    return $return_text;
  }else{
    return t('Seems you provided an incorrect identifier');
  }
}

/**
 * View experiments
 *
 * @param $experiment_id
 * Can be an array or string. Can be a numerical identifier or uniquename
 *
 * @return Tabulated HTML
 */
function genes4all_experiment_page_view_experiment($experiment_id = 0) {
  $return_text = t("<p>This page allows you to view details of specific experiments.</p>");
  if (strpos($experiment_id,',')!==FALSE){
    $experiment_id = explode ($experiment_id,',');
  }
  if (empty($experiment_id)) {
    return $return_text.t('No experiment found with this identifier.');
  }elseif (is_array($experiment_id)){
   foreach ($experiment_id as $id){
    $return_text.= _genes4all_experiment_experiment_details($id);
   }
  }else{
    $return_text.= _genes4all_experiment_experiment_details($experiment_id);
  }
  return $return_text;
}

