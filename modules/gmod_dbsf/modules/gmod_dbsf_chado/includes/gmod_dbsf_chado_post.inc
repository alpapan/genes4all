<?php


/**
 * @file
 *
 * HTTP POST handling for Drupal and Chado
 *
 * This module is to provide Ajax interactivity with JS solutions
 * for example it can
 * a) redirect data from Ext to Drupal, handling security enroute
 * b) Receive a POST request (eg an ajax call), do something within
 *    drupal (e.g. run a batch involving data preprocessing and R )
 *    and then send a reply to the caller. Optionally, the Ajax call
 *    may have redirected us to a Drupal window
 *
 * @see geoGenes.js
 *
 * @author Nathan Dunn <nathan.dunn@UNSWalumni.com>
 * @author Alexie Papanicolaou
 * @see complaints mailto:alexie@butterflybase.org
 *
 */

/*
 * Check if user matches key session and give a response
 * if the user does not have an active session that matches
 * this key, then we give an informative signal for Ajax
 *
 * @TODO AP what is an informative signal?
 *
 */
function gmod_dbsf_chado_post_security($sid=NULL){
  $permission = FALSE;
  // if user is accessing a webservice with the browser without
  // an sid, they can still proceed if they have a recent drupal window logged in
  if (empty($sid)){
    global $user;
    if (!empty($user)){
      $sid = $user->sid;
    }
  }

  $session_uid = "SELECT DISTINCT uid from sessions WHERE sid='%s'";
  $res = db_fetch_array(db_query($session_uid,$sid));
  $uid = $res['uid'];
  // inactive session:
  if (empty($uid)){
    //then check if the permissions of 'Access webservices via Drupal' have been set for anonymous
    // in that case, we allow everyone!
    $perm_array = array();
    $get_perms_sql = "SELECT perm from {permission} where rid IN ((select rid from {role} where name='anonymous user'))";
    $res = db_query($get_perms_sql,$uid);
    while ($data = db_fetch_array($res)){
      $perm_array=explode(',',$data['perm']);
      foreach ($perm_array as $key=>$value){
        $perm_array[$key]= trim($value);
      }
      if (in_array('Access POST (analytics)', $perm_array)){
        $permission = TRUE;
      }
    }
  }else{
    $perm_array = array();
    $get_perms_sql = "SELECT perm from {permission} where rid IN ((select rid from {role} where name='authenticated user'),(SELECT rid from {users_roles} where uid=%d))";
    $res = db_query($get_perms_sql,$uid);
    while ($data = db_fetch_array($res)){
      $perm_array=explode(',',$data['perm']);
      foreach ($perm_array as $key=>$value){
        $perm_array[$key]= trim($value);
      }
      if (in_array('Access POST (analytics)', $perm_array)){
        $permission = TRUE;
      }
    }
  }
  return $permission;
}

/*
 *
 * Store POST requests in SESSION
 *
 * @TODO AP check security implementation
 * @return TRUE on success
 *
*/
function gmod_dbsf_chado_post_receive() {
  if (empty($_POST)){return FALSE;}
  $data = array();
  $data = json_decode($_POST['gmod_dbsf_chado']);

  $authentication = gmod_dbsf_chado_post_security($data->sid);
/*
  if (empty($authentication)){
    $_SESSION['gmod_dbsf_chado'] = 'FAILED AUTHENTICATION';
    return $authentication;
  }else{
    $_SESSION['gmod_dbsf_chado']=json_encode($data);
  }
    return $authentication;*/

}

/**
 *
 * Converts an object to an associative array recursively
 * @see http://codesnippets.joyent.com/posts/show/1641
 * @return array
 */
function gmod_dbsf_chado_object_to_array($data) {
    if(is_array($data) || is_object($data)) {
        $result = array();
        foreach($data as $key => $value) {
            $result[$key] = gmod_dbsf_chado_object_to_array($value);
        }
        return $result;
    }
    return $data;
}

/**
 *
 * Echo the posted data to browser
 *
 * this is a temporary function, to be used as a skeleton
 * for a more interactive one
 *
 * @return html
 */
function gmod_dbsf_chado_post_echo() {
  if (empty($_SESSION)){
    header("HTTP/1.0 400 Bad Request",TRUE,400);
    return ' ';
  }
    $output ='';
    if (empty($_SESSION['gmod_dbsf_chado'])){
      $output.='No posted information stored in session';
    }elseif ($_SESSION['gmod_dbsf_chado'] == 'FAILED AUTHENTICATION') {
      header("HTTP/1.0 401 Unauthorized",TRUE,401);
      return 'FAILED AUTHENTICATION';
    }else{
         $output = $_SESSION['gmod_dbsf_chado'];
         //dpm(json_decode($output));
         json_decode($output);
         //unset($_SESSION['gmod_dbsf_chado']);
    }
    if (empty($output)){
      return ' ';
    }
    return serialize($output);
}

/*
 * Receive a POST request, do something and return some data
 *
 */
function gmod_dbsf_chado_post_request(){
  if (empty($_POST)){return FALSE;}
  $authentication = gmod_dbsf_chado_post_security($_POST['sid']);
  if (empty($authentication)){
    return 'FAILED AUTHENTICATION';
  }
  //do something.
}


