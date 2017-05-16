<?php
// $Id$

/**
 * @file
 * include file for adding new experiments
 */
function genes4all_experiment_page_added() {
  $text = '<p>You have completed submitting your experiment. You may review/edit in the future using your passkey and revisiting the '. l('experiment database', 'genes4all/experiment/add') .'.</p>';
  $overview = genes4all_experiment_tabulate_study_data('all', 'all');
  gmod_dbsf_clear_session();
  gmod_dbsf_clear_session(array('study_id', 'passkey', 'pub_id', 'animals_id', 'target_id', 'construct_id', 'delivery_id', 'assay_id',
      'study_type', 'study_name', 'pub_name', 'animals_name', 'target_name', 'construct_name', 'delivery_name', 'assay_name',
    ));
  return $text . $overview;
}

/**
 * Add new experiments
 *
 * @param $lib_id
 *
 * @return unknown_type
 */
function genes4all_experiment_page_add() {
  $rendered_form = '<p>Welcome,<br>here you can submit results for your experiments. Negative results are particularly appreciated. If you think there some kind information which is important but no option is available, please leave feedback.</p>';
  $rendered_form .= "<p><b>NB: Please note that we have identified issues with Internet Explorer 7 being incapable of rendering CSS pages correctly. This is known to Microsoft and has been corrected in IE8. Please use the latest version, or use another browser (such as Firefox).</b> You can read the ". l('Wikipedia entry', 'http://en.wikipedia.org/wiki/Acid2#Microsoft.27s_response') .".</p>";
  $rendered_form .= "<p>Please read ". l('the tutorial', 'tutorial_experiments', array(alias => TRUE, 'attributes' => array('target' => '_blank'))) ." before proceeding.</p>";
  $study_type_sql = "SELECT name from {gmod_dbsf_cvterm} as cvterm where cvterm_id=(SELECT type_id from {gmod_dbsf_study} as study where study_id=%d)";
  $study_id       = check_plain($_SESSION['study_id']);
  $study_name     = check_plain($_SESSION['study_name']);
  if (empty($study_id)) {
    //clear session
    if (!empty($_SESSION['passkey']) || !empty($_SESSION['pub_id']) || !empty($_SESSION['animals_id']) || !empty($_SESSION['target_id']) || !empty($_SESSION['construct_id']) || !empty($_SESSION['delivery_id']) || !empty($_SESSION['assay_id'])) {
      gmod_dbsf_clear_session(array('passkey', 'pub_id', 'animals_id', 'target_id', 'construct_id', 'delivery_id', 'assay_id',
          'pub_name', 'animals_name', 'target_name', 'construct_name', 'delivery_name', 'assay_name',
        ));
    }
    $rendered_form .= drupal_get_form('genes4all_experiment_study_form');
    return $rendered_form;
  }
  else {
    $study_type = db_fetch_array(db_query($study_type_sql, $study_id));
    $study_type = $study_type['name'];
    //$_SESSION['study_type'] holds what the study will become when it is finalized.
    //$study_type holds the status of the study which ==$_SESSION['study_type'] if it
    // is finalized or 'unfinished' if it is not.
    if ($study_type == 'unfinished') {
      drupal_set_message(t('You are editing %study which is currently <strong>not</strong> finalized.', array('%study' => $study_name)), 'warning', FALSE);
    }
    else {
      drupal_set_message(t('You are editing %study which has been set as finalized.', array('%study' => $study_name)), 'warning', FALSE);
    }
    switch ($_SESSION['study_type']) {
      case 'RNAi':
        return genes4all_experiment_rnai_page();
    }
  }
}

function genes4all_experiment_rnai_page() {
  $form_pub       = '';
  $form_animal    = '';
  $form_target    = '';
  $form_construct = '';
  $form_delivery  = '';
  $form_assay     = '';
  $rendered_form  = '';
  $passkey        = $_SESSION['passkey'];
  // Security
  if (empty($passkey)) {
    gmod_dbsf_clear_session(array('study_id', 'passkey', 'pub_id', 'animals_id', 'target_id', 'construct_id', 'delivery_id', 'assay_id',
        'study_type', 'study_name', 'pub_name', 'animals_name', 'target_name', 'construct_name', 'delivery_name', 'assay_name',
      ));

    $rendered_form .= drupal_get_form('genes4all_experiment_study_form');
    return $rendered_form;
  }

  if (empty($_SESSION['pub_id'])) {
    $form_pub = drupal_get_form('genes4all_experiment_pub_populate_form');
    $form_pub .= drupal_get_form('genes4all_experiment_pub_form');
  }
  else {
    // present a table of the values selected already
    // and an option to reset via gmod_dbsf_clear_session
    $form_pub = drupal_get_form('genes4all_experiment_reset_session_form' . '_pub_id', 'pub_id');
    $form_pub .= genes4all_experiment_tabulate_study_data($_SESSION['pub_id'], 'pub_id');
  }
  if (empty($_SESSION['animals_id'])) {
    $form_animal = drupal_get_form('genes4all_experiment_animals_populate_form');
    $form_animal .= drupal_get_form('genes4all_experiment_animals_form');
  }
  else {
    $form_animal = drupal_get_form('genes4all_experiment_reset_session_form' . '_animals_id', 'animals_id');
    $form_animal .= genes4all_experiment_tabulate_study_data($_SESSION['animals_id'], 'animals_id');
  }

  if (empty($_SESSION['construct_id'])) {
    $form_construct = drupal_get_form('genes4all_experiment_construct_populate_form');
    $form_construct .= drupal_get_form('genes4all_experiment_construct_form');
  }
  else {
    $form_construct = drupal_get_form('genes4all_experiment_reset_session_form' . '_construct_id', 'construct_id');
    $form_construct .= genes4all_experiment_tabulate_study_data($_SESSION['construct_id'], 'construct_id');
  }
  if (empty($_SESSION['target_id'])) {
    $form_target = drupal_get_form('genes4all_experiment_target_populate_form');
    $form_target .= drupal_get_form('genes4all_experiment_target_form');
  }
  else {
    $form_target = drupal_get_form('genes4all_experiment_reset_session_form' . '_target_id', 'target_id');
    $form_target .= genes4all_experiment_tabulate_study_data($_SESSION['target_id'], 'target_id');
  }
  if (empty($_SESSION['delivery_id'])) {
    $form_delivery = drupal_get_form('genes4all_experiment_delivery_populate_form');
    $form_delivery .= drupal_get_form('genes4all_experiment_delivery_form');
  }
  else {
    $form_delivery = drupal_get_form('genes4all_experiment_reset_session_form' . '_delivery_id', 'delivery_id');
    $form_delivery .= genes4all_experiment_tabulate_study_data($_SESSION['delivery_id'], 'delivery_id');
  }
  if (empty($_SESSION['assay_id'])) {
    $form_assay = drupal_get_form('genes4all_experiment_assay_populate_form');
    $form_assay .= drupal_get_form('genes4all_experiment_assay_form');
  }
  else {
    $form_assay = drupal_get_form('genes4all_experiment_reset_session_form' . '_assay_id', 'assay_id');
    $form_assay .= genes4all_experiment_tabulate_study_data($_SESSION['assay_id'], 'assay_id');
  }

  $overview = drupal_get_form('genes4all_experiment_study_finished_form') . genes4all_experiment_tabulate_study_data('all', 'all');

  $page_tabs = array(
    'main' => array(
      '#type' => 'tabset',
      'genes' => array(
        '#type' => 'tabpage',
        '#title' => 'Sequences',
        'genes_set' => array(
          '#type' => 'tabset',
          'construct' => array(
            '#type' => 'tabpage',
            '#title' => 'RNAi construct',
            '#content' => '<h1>Silencing RNAi construct information</h1>'. $form_construct,
            '#weight' => 2,
          ),
          'target' => array(
            '#type' => 'tabpage',
            '#title' => 'Target gene',
            '#content' => '<h1>Target gene information</h1>'. $form_target,
            '#weight' => 1,
          ),
        ),
      ),
      'resources' => array(
        '#type' => 'tabpage',
        '#title' => 'Resources',
        'resources_set' => array(
          '#type' => 'tabset',
          'Animals' => array(
            '#type' => 'tabpage',
            '#title' => 'Experimental animals',
            '#content' => '<h1>Experimental animals</h1>'. $form_animal,
            '#weight' => 3,
          ),
          'delivery' => array(
            '#type' => 'tabpage',
            '#title' => 'Delivery protocol',
            '#content' => '<h1>Delivery protocol</h1>'. $form_delivery,
            '#weight' => 4,
          ),
          'assay' => array(
            '#type' => 'tabpage',
            '#title' => 'Assay protocol',
            '#content' => '<h1>Assay protocol</h1>'. $form_assay,
            '#weight' => 5,
          ),
        ),
      ),
      'finalize' => array(
        '#type' => 'tabpage',
        '#title' => 'Finalize',
        'finalize_set' => array(
          '#type' => 'tabset',
          'pub' => array(
            '#type' => 'tabpage',
            '#title' => 'Publication data ',
            '#content' => '<h1>Publications</h1>'. $form_pub,
            '#weight' => 0,
          ),
          'validate' => array(
            '#type' => 'tabpage',
            '#title' => 'Overview & Validate',
            '#content' => '<h1>Overview</h1>'. $overview,
            '#weight' => 6,
          ),
        ),
      ),
    ),
  );
  $rendered_form = drupal_get_form('genes4all_experiment_reset_session_form');
  $rendered_form .= tabs_render($page_tabs);
  $rendered_form .= drupal_get_form('genes4all_experiment_reset_session_form');

  //dpm($_SESSION);
  return $rendered_form;
}

function genes4all_experiment_study_form($form_state) {
  global $user;
  $user_email = ($user->uid) > 0 ? $user->mail : '';
  $study_types = array(
    'RNAi' => 'RNAi silencing experiment',
  );

  $experiment_type_form = array(
    'study_type' => array(
      '#title' => 'Select the type of experiment you want to submit',
      '#weight' => 0,
      '#type' => 'select',
      '#size' => 1,
      '#options' => $study_types,
      '#required' => TRUE,
    ),
  );
  $choose_study_form = array(
    'study_name' => array(
      '#weight' => 1,
      '#type' => 'textfield',
      '#title' => 'Friendly name for your study',
      '#description' => t('Please give a name for your study that is unique to you. This is used to avoid having to type the details again: if you submitted these details before, type the name you gave.'),
      '#required' => TRUE,
      '#autocomplete_path' => 'genes4all/experiment/study/autocomplete',
    ), 'submitter_email' => array(
      '#weight' => 2,
      '#type' => 'textfield',
      '#title' => 'Your email address',
      '#description' => t('Please provide your (the submitter) email address where we can contact you if there are any issues. This can be the same or different from the "communicating author" responsible for the study you are submitting. If you have an account here, feel free to login first. Unlike the communicating author, these addresses are never made public.'),
      '#required' => TRUE,
      '#default_value' => $user_email,
    ), 'passkey' => array(
      '#type' => 'textfield',
      '#title' => t('Passkey'),
      '#description' => t('Give a simple latin alphanumeric (A-Z,a-z,0-9) passkey so you and your collaborators can access/edit your data in the future but prevent others from accessing/editing them. CAUTION: this is NOT stored in a secured fashion so DO NOT use any of your secure passwords. You may, however, use the same passkey for multiple submissions.'),
      '#required' => TRUE,
      '#weight' => 3,
      '#size' => 15,
      '#maxlength' => 15,
    ), 'study_comments' => array(
      '#title' => 'Abstract / Additional information on study',
      '#type' => 'textarea',
      '#weight' => 4,
      '#description' => t("Feel free to give or update a description or <i>abstract</i> for your study. If you're updating an existing study, you don't need to retype the description. This field will also appear under Finalize->Publication on the next page."),
    ),
  );
  $form['buttons'] = array(
    'submit_upper' => array(
      '#type' => 'submit',
      '#value' => t('Start'),
      '#weight' => -25,
    ),
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Start'),
      '#weight' => 25,
    ),
    'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'button',
      '#value' => t('Reset'),
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 26,
    ),
  );
  return $choose_study_form + $experiment_type_form + $form['buttons'];
}

function genes4all_experiment_study_form_validate($form, &$form_state) {
  $passkey  = $form_state['values']['passkey'];
  $email    = $form_state['values']['submitter_email'];
  $is_valid = valid_email_address($email);
  if (empty($is_valid)) {
    form_set_error('submitter_email', t('Your email address does not look valid.'));
  }
  elseif (strlen($email > 40)) {
    form_set_error('submitter_email', t('Your email address is longer than expected. Are you sure it is valid? If it is, please .') . l('contact', 'contact') . t(' an administrator to help you.')
    );
  }
  if (!empty($passkey)) {
    if (ctype_alnum($passkey) !== TRUE) {
      form_set_error('passkey', t('Passkey is not alphanumeric.'));
    }
    elseif (strlen($passkey) > 15) {
      form_set_error('passkey', t('Please keep the passkey shorter than 16 characters (%f).',array('%f'=>strlen($passkey))));
    }
    else {
      $_SESSION['passkey'] = $passkey;
    }
  }
}

function genes4all_experiment_study_form_submit($form, &$form_state) {
  $check = genes4all_experiment_passkey_check($form_state['values']['study_name'], $form_state['values']['passkey']);
  if ($check === FALSE) {
    drupal_set_message(t('You did not provide the correct passphrase for this study. Are you sure you own this study? Perhaps you should create a new one?'), 'error', TRUE);
    gmod_dbsf_form_clear($form, &$form_state);
    return FALSE;
  }
  $return_array = genes4all_experiment_study_add($form_state['values']);
  //dpm($return_array);
  if (empty($return_array['study_id'])) {
    drupal_set_message(t('I was not able to store the study.'), 'error', TRUE);
    gmod_dbsf_form_clear($form, &$form_state);
    return FALSE;
  }
  $_SESSION['study_type'] = $form_state['values']['study_type'];
  $_SESSION['study_name'] = $form_state['values']['study_name'];
  $_SESSION['study_id'] = $return_array['study_id'];
  // this is not generalized yet.
  // only _id is needed really...
  if (!empty($return_array['study_status'])) {
    $_SESSION['study_status'] = $return_array['study_status'];
  }
  if (!empty($return_array['pub_name'])) {
    $_SESSION['pub_name'] = $return_array['pub_name'];
    $_SESSION['pub_id'] = $return_array['pub_id'];
  }
  if (!empty($return_array['rnai_construct'])) {
    $_SESSION['construct_name'] = $return_array['rnai_construct']['name'];
    $_SESSION['construct_id'] = $return_array['rnai_construct']['id'];
  }
  if (!empty($return_array['gene_target'])) {
    $_SESSION['target_name'] = $return_array['gene_target']['name'];
    $_SESSION['target_id'] = $return_array['gene_target']['id'];
  }
  if (!empty($return_array['assay protocol'])) {
    $_SESSION['assay_name'] = $return_array['assay protocol']['name'];
    $_SESSION['assay_id'] = $return_array['assay protocol']['id'];
  }
  if (!empty($return_array['delivery protocol'])) {
    $_SESSION['delivery_name'] = $return_array['delivery protocol']['name'];
    $_SESSION['delivery_id'] = $return_array['delivery protocol']['id'];
  }
  if (!empty($return_array['experimental animals'])) {
    $_SESSION['animals_name'] = $return_array['experimental animals']['name'];
    $_SESSION['animals_id'] = $return_array['experimental animals']['id'];
  }
  return FALSE;
}

function genes4all_experiment_pub_populate_form($form_state) {
  $data = genes4all_experiment_publication_add();
  if (empty($data)) {
    $data = array('0' => 'None saved by your passkey');
  }
  $form = array(
    'existing' => array(
      '#type' => 'fieldset',
      '#title' => 'Existing publications',
      'pub_id' => array(
        '#type' => 'select',
        '#options' => $data,
      ),
      'populate' => array(
        '#type' => 'submit',
        '#value' => 'Populate',
      ),
    ),
  );
  return $form;
}

function genes4all_experiment_pub_populate_form_submit($form, &$form_state) {
  $_SESSION['pub_id'] = $form_state['values']['pub_id'];
}

function genes4all_experiment_pub_form($form_state) {
  $form = array();
  $study_id          = check_plain($_SESSION['study_id']);
  $check_sql         = "SELECT description from {gmod_dbsf_study} WHERE study_id=%d";
  $study_description = '';
  if (!empty($study_id)) {
    $study_description = db_fetch_array(db_query($check_sql, $study_id));
    $study_description = $study_description['description'];
  }
$form[] = array(
    'basic' => array(
      '#type' => 'fieldset',
      '#weight' => 0,
      '#title' => t('Publication'),
      'author_first1' => array(
        '#type' => 'textfield',
        '#title' => t('Communicating author first/middle name(s)'),
        '#required' => TRUE,
      ),
      'author_last1' => array(
        '#type' => 'textfield',
        '#title' => t('Communicating author family name(s)'),
        '#required' => TRUE,
      ),
      'author_email1' => array(
        '#type' => 'textfield',
        '#title' => t('Communicating author email address'),
        '#required' => TRUE,
      ),
    ),
    'advanced' => array(
      '#type' => 'fieldset',
      '#title' => t('Advanced'),
      '#weight' => 1,
      'pub_db' => array(
        '#type' => 'textfield',
        '#title' => t('Database with publication'),
        '#description' => t('If published, a public database which indexes your publication, e.g. Pubmed, WoK etc.'),
        '#size' => 15,
        '#autocomplete_path' => 'genes4all/experiment/dbname/autocomplete',
      ),
      'pub_db_accession' => array(
        '#type' => 'textfield',
        '#title' => t('Publication id in above database'),
        '#description' => t('Only needed if your publication is published and indexed in the above database.'),
        '#weight' => 2,
      ),
    ),
    'abstract' => array(
      '#type'=> 'fieldset',
      '#title' => t('Extras'),
      '#type' => 'textarea',
      '#title' => 'Abstract',
      '#description' => t('Feel free to give an abstract or a description of your study. This is the same as the "Additional information on study" information requested in the beginning of the submission.'),
      '#default_value' => $study_description,
      '#weight' => 3,
    ),
  );

  $form[] = array(
    '#weight' => -10,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Store'),
      '#weight' => 0,
    ),
    'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'button',
      '#value' => t('Reset'),
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 1,
    ),
  );

  $form[] = array(
    '#weight' => 25,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Store'),
      '#weight' => 0,
    ),
    'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'button',
      '#value' => t('Reset'),
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 1,
    ),
  );

  return $form;
}

function genes4all_experiment_pub_form_validate($form, &$form_state) {
  $email = $form_state['values']['author_email1'];
  $check = valid_email_address($email);
  if ($check !== 1) {
    form_set_error("author_email1", t('Email address (%email) is not valid.', array('%email' => $email)));
  }
}

function genes4all_experiment_pub_form_submit($form, &$form_state) {
  $form_state['values']['pub_uniquename'] = $_SESSION['study_name'] .' publication';
  $description = check_plain($form_state['values']['abstract']);
  $pub_id = genes4all_experiment_publication_add($form_state['values']);
  if (empty($pub_id)) {
    drupal_set_message(t('I was not able to store the publication data.'), 'error', TRUE);
    gmod_dbsf_form_clear($form, &$form_state);
    return FALSE;
  }
  $_SESSION['pub_id'] = $pub_id;
  $_SESSION['pub_name'] = $form_state['values']['pub_uniquename'];

  $success = genes4all_experiment_authors_add($form_state['values'], $pub_id);
  if (!empty($success === FALSE)) {
    drupal_set_message(t('I was not able to store the authorship data.'), 'error', TRUE);
    gmod_dbsf_form_clear($form, &$form_state);
    return FALSE;
  }
  if (!empty($description)) {
    $study_id = check_plain($_SESSION['study_id']);
    $update_sql = "UPDATE {gmod_dbsf_study} set description='%s' where study_id=%d";
    db_query($update_sql, $description, $study_id);
  }

  //dpm($form_state);
}

function genes4all_experiment_target_populate_form($form_state) {
  $data = genes4all_experiment_feature_get('', 'gene_target');

  if (empty($data)) {
    $data = array('0' => 'None saved by your passkey');
  }
  $form = array(
    'existing' => array(
      '#type' => 'fieldset',
      '#title' => 'Existing targets',
      'target_id' => array(
        '#type' => 'select',
        '#options' => $data,
      ),
      'populate' => array(
        '#type' => 'submit',
        '#value' => 'Populate',
      ),
    ),
  );
  return $form;
}

function genes4all_experiment_target_populate_form_submit($form, &$form_state) {
  $_SESSION['target_id'] = $form_state['values']['target_id'];
}

function genes4all_experiment_target_form($form_state) {
  $target_aln_region_array = array(
    '3utr' => "3'UTR",
    '5utr' => "5'UTR",
    'intron' => 'Intron',
    'cds' => 'CDS',
  );

  $form[] = array(
    'target' => array(
      '#weight' => 0,
      '#type' => 'fieldset',
      '#title' => 'Target gene information',
      'target_seq' => array(
        '#weight' => 0,
        '#type' => 'textarea',
        '#title' => 'DNA sequence of gene',
      ),
      'aln_region' => array(
        '#title' => 'Regions of target aligning with construct',
        '#type' => 'checkboxes',
        '#options' => $target_aln_region_array,
        '#required' => TRUE,
        '#weight' => 1,
      ),
      'taxonomy' => array(
        '#type' => 'fieldset',
        '#title' => "Taxonomy",
        '#description' => t('Species from which this data originates'),
        '#collapsible' => TRUE,
        '#weight' => 10,
        'tax_class' => array(
          '#type' => 'textfield',
          '#title' => t('Class'),
          '#required' => TRUE,
          '#weight' => -2,
          '#size' => 20,
          '#default_value' => 'e.g. Insecta',
        ),
        'tax_order' => array(
          '#type' => 'textfield',
          '#title' => t('Order'),
          '#required' => TRUE,
          '#weight' => -1,
          '#size' => 20,
          '#default_value' => 'e.g. Lepidoptera',
        ),
        'tax_family' => array(
          '#type' => 'textfield',
          '#title' => t('Family'),
          '#required' => TRUE,
          '#weight' => 0,
          '#size' => 20,
          '#default_value' => 'e.g. Noctuidae',
        ),
        'genus' => array(
          '#type' => 'textfield',
          '#title' => t('Genus'),
          '#description' => t('the first epithet/attribute of the binomial latin name'),
          '#required' => TRUE,
          '#weight' => 1,
          '#size' => 20,
        ),
        'species' => array(
          '#type' => 'textfield',
          '#title' => t('Species'),
          '#description' => t('the second epithet/attribute of the binomial latin name'),
          '#required' => TRUE,
          '#weight' => 2,
          '#size' => 20,
        ),
        'ncbi_taxid' => array(
          '#type' => 'textfield',
          '#title' => t('NCBI Taxonomy ID'),
          '#description' => t('Obtainable from the NCBI ') . l('Taxonomy database', 'http://www.ncbi.nlm.nih.gov/Taxonomy/', array(
              'attributes' => array('target' => '_blank'),
            )
          ),
          '#weight' => 4,
          '#size' => 8,
        ),
      ),
    ),
    'advanced' => array(
      '#type' => 'fieldset',
      '#title' => 'Advanced',
      '#weight' => 20,
      'target_db' => array(
        '#weight' => 0,
        '#type' => 'textfield',
        '#title' => 'Repository name',
        '#description' => t('The name of the database or repository from which you retrieved the sequence or its orthologue resides. Just leave it blank if none was used.'),
        '#size' => 17,
        '#autocomplete_path' => 'genes4all/experiment/dbname/autocomplete',
      ),
      'target_dbxref' => array(
        '#weight' => 1,
        '#type' => 'textfield',
        '#title' => 'Database gene ID',
        '#description' => t('The ID given to the gene of interest by the above database. Leave it as it is if none was used.'),
        '#size' => 17,
      ),
    ),
  );

  $form[] = array(
    '#weight' => -10,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Store'),
      '#weight' => 0,
    ),
    'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'button',
      '#value' => t('Reset'),
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 1,
    ),
  );
  $form[] = array(
    '#weight' => 25,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Store'),
      '#weight' => 0,
    ),
    'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'button',
      '#value' => t('Reset'),
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 1,
    ),
  );


  return $form;
}

function genes4all_experiment_target_form_validate($form, &$form_state) {
  $ncbi_taxid = $form_state['values']['ncbi_taxid'];
  $term_error = t('A taxonomic term cannot have spaces in it!');
  if (strpos($form_state['values']['tax_class'], ' ') !== FALSE) {
    form_set_error('tax_class', $term_error);
  }
  if (strpos($form_state['values']['tax_order'], ' ') !== FALSE) {
    form_set_error('tax_order', $term_error);
  }
  if (strpos($form_state['values']['tax_family'], ' ') !== FALSE) {
    form_set_error('tax_family', $term_error);
  }
  if (strpos($form_state['values']['genus'], ' ') !== FALSE) {
    form_set_error('genus', $term_error);
  }
  if (strpos($form_state['values']['species'], ' ') !== FALSE) {
    form_set_error('species', $term_error);
  }

  $taxa = array(
    'tax_class' => $form_state['values']['tax_class'],
    'tax_order' => $form_state['values']['tax_order'],
    'tax_family' => $form_state['values']['tax_family'],
    'genus' => $form_state['values']['genus'],
    'species' => $form_state['values']['species'],
  );
  foreach ($taxa as $key => $value) {
    if (!ctype_alpha($value)) {
      form_set_error($key, t('Highlighted taxon does not seem to be composed of only alphabetic characters.'));
    }
  }
  if (!empty($ncbi_taxid)) {
    if (!is_numeric($ncbi_taxid)) {
      form_set_error("ncbi_taxid", t('NCBI Taxonomy ID is not valid.'));
    }
  }
  $dna = $form_state['values']['target_seq'];
  if (!empty($dna)) {
    $dna = gmod_dbsf_validate_seq_dna($dna);
    if ($dna == FALSE || is_numeric($dna)) {
      form_set_error('target_seq', t('Target sequence does not seem to be a valid DNA sequence (%dna errors).', array('%dna' => $dna)));
    }
    else {
      $form_state['values']['target_seq'] = $dna;
    }
  }
}

function genes4all_experiment_target_form_submit($form, &$form_state) {
  // create fake name
  $form_state['values']['target_name'] = $_SESSION['study_name'] .' target';
  $feature_id = genes4all_experiment_feature_add($form_state['values'], $form_state['values']['target_name'], 'gene_target');
  if (!empty($feature_id)) {
    $_SESSION['target_id'] = $feature_id;
    $_SESSION['target_name'] = $form_state['values']['target_name'];
  }
  else {
    drupal_set_message(t('I was not able to store the feature data.'), 'error', FALSE);
    gmod_dbsf_form_clear($form, &$form_state);
    return FALSE;
  }
}

function genes4all_experiment_construct_populate_form($form_state) {
  $data = genes4all_experiment_feature_get('', 'rnai_construct');
  if (empty($data)) {
    $data = array('0' => 'None saved by your passkey');
  }
  $form = array(
    'existing' => array(
      '#type' => 'fieldset',
      '#title' => 'Existing constructs',
      'construct_id' => array(
        '#type' => 'select',
        '#options' => $data,
      ),
      'populate' => array(
        '#type' => 'submit',
        '#value' => 'Populate',
      ),
    ),
  );

  return $form;
}

function genes4all_experiment_construct_populate_form_submit($form, &$form_state) {
  $_SESSION['construct_id'] = $form_state['values']['construct_id'];
}

function genes4all_experiment_construct_form($form_state) {
  //TODO these should come from the CV...
  $rna_probe_type_array = array(
    'siRNA' => 'siRNA',
    'hairpinRNA' => 'hairpinRNA',
    'dsRNA' => 'dsRNA',
  );
  $purification_method_array = array(
    'salt precipitation' => 'Salt Precipitation',
    'ethanol precipitation' => 'Ethanol Precipitation',
    'isopropanol precipitation' => 'Isopropanol Precipitation',
    'column-based' => 'Column-based',
    'phenol_chloroform' => 'Phenol-Chloroform extraction',
    'no_purification' => 'No purification performed',
    'bought_purified' => 'Bought purified and commercial method unknown',
    'unknown' => 'Unknown if/how purified',
    'other' => 'Other',
  );
  $annealing_method_array = array(
    'simultaneously with transcription' => 'Simultaneously with transcription',
    'subsequent to transcription' => 'Subsequent to transcription',
  );

  $form = array(
    'advanced' => array(
      '#type' => 'fieldset',
      '#weight' => 10,
      '#title' => 'Advanced',
      'construct_db' => array(
        '#title' => 'Repository name',
        '#description' => t('The name of the database or repository which holds sequence information for your construct. Just leave it as local feature if none is used.'),
        '#weight' => 0,
        '#type' => 'textfield',
        '#size' => 17,
        '#autocomplete_path' => 'genes4all/experiment/dbname/autocomplete',
      ),
      'construct_dbxref' => array(
        '#weight' => 1,
        '#type' => 'textfield',
        '#title' => 'Database gene ID',
        '#description' => t('The ID given to the construct by the above database.'),
        '#size' => 17,
      ),
    ),
    'basic' => array(
      '#type' => 'fieldset',
      '#weight' => 0,
      '#title' => 'Basic info',
      'rna_probe_type' => array(
        '#weight' => 4,
        '#type' => 'radios',
        '#title' => 'Type of construct',
        '#description' => t('siRNA: double-stranded RNA molecules that are 20-25 nucleotides in length.<br>hairpinRNA: RNA that makes a tight hairpin turn'),
        '#options' => $rna_probe_type_array,
        '#required' => TRUE,
      ),
      'construct_seq' => array(
        '#weight' => 3,
        '#type' => 'textarea',
        '#title' => 'Sequence (as DNA)',
        '#required' => TRUE,
      ),
      'construct_protocol' => array(
        '#weight' => 5,
        '#type' => 'fieldset',
        '#title' => 'Construct preparation protocol',
        // 'construct_protocol_invitro'=>array(
        //   '#weight' => 1,
        //   '#title'=>'Was construct prepared in-vitro?',
        //   '#type'=>'checkbox',
        // ),
        'construct_protocol_kit' => array(
          '#weight' => 2,
          '#title' => t('<i>If construct was prepared in-vitro</i>, what was the name of kit/company used to prepare construct'),
          '#type' => 'textfield',
        ),
        'purification' => array(
          '#weight' => 4,
          '#type' => 'radios',
          '#title' => 'Purification method',
          '#required' => TRUE,
          '#options' => $purification_method_array,
        ),
        'annealing_method' => array(
          '#weight' => 4,
          '#type' => 'radios',
          '#title' => 'Annealing method of dsRNA',
          '#required' => TRUE,
          '#options' => $annealing_method_array,
        ),
        'construct_protocol_details' => array(
          '#weight' => 3,
          '#type' => 'textarea',
          '#title' => 'Detailed protocol',
          '#required' => TRUE,
        ),
      ),
    ),
  );
  $form[] = array(
    '#weight' => -10,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Store'),
      '#weight' => 0,
    ),
    'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'button',
      '#value' => t('Reset'),
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 1,
    ),
  );
  $form[] = array(
    '#weight' => 25,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Store'),
      '#weight' => 0,
    ),
    'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'button',
      '#value' => t('Reset'),
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 1,
    ),
  );
  return $form;
}

function genes4all_experiment_construct_form_validate($form, &$form_state) {
  $dna = $form_state['values']['construct_seq'];
  if (!empty($dna)) {
    $dna = gmod_dbsf_validate_seq_dna($dna);
    if ($dna == FALSE || is_numeric($dna)) {
      form_set_error('construct_seq', t('Construct sequence does not seem to be a valid DNA sequence (%dna errors).', array('%dna' => $dna)));
    }
    else {
      $form_state['values']['construct']['construct_seq'] = $dna;
    }
  }
}

function genes4all_experiment_construct_form_submit($form, &$form_state) {
  // create fake name
  $form_state['values']['construct_name'] = $_SESSION['study_name'] .' RNAi construct';
  $feature_id = genes4all_experiment_feature_add($form_state['values'], $form_state['values']['construct_name'], 'rnai_construct');
  if (!empty($feature_id)) {
    $_SESSION['construct_id'] = $feature_id;
    $_SESSION['construct_name'] = $form_state['values']['construct_name'];
  }
  else {
    drupal_set_message(t('I was not able to store the feature data.'), 'error', FALSE);
    gmod_dbsf_form_clear($form, &$form_state);
    return FALSE;
  }
}

function genes4all_experiment_animals_populate_form($form_state) {
  $data = genes4all_experiment_resource_get('', 'experimental animals');
  if (empty($data)) {
    $data = array('0' => 'None saved by your passkey');
  }
  $form = array(
    'existing' => array(
      '#type' => 'fieldset',
      '#title' => 'Existing experimental animals',
      'animals_id' => array(
        '#type' => 'select',
        '#options' => $data,
      ),
      'populate' => array(
        '#type' => 'submit',
        '#value' => 'Populate',
      ),
    ),
  );

  return $form;
}

function genes4all_experiment_animals_populate_form_submit($form, &$form_state) {
  $_SESSION['animals_id'] = $form_state['values']['animals_id'];
}

function genes4all_experiment_animals_form($form_state) {
  $form = array();
  $stage_array = array(
    'unfertilized egg' => 'Unfertilized egg',
    'embryo' => 'Embryo',
    'larva' => 'Larva',
    'pupa' => 'Pupa',
    'adult' => 'Adult',
  );
  $origin_array = array(
    'lab_colony' => 'Lab-colony',
    'field_collected' => 'Field-collected',
  );
  $infection_array = array(
    'virus' => 'Virus',
    'bacteria' => 'Bacteria',
    'fungi' => 'Fungi',
    'nematode' => 'Nematode',
    'microsporidia' => 'Microsporidia',
    'mite' => 'Mites',
    'unknown' => 'Present but unknown',
  );
  $form[] = array(
    'colony' => array(
      '#weight' => 0,
      '#type' => 'fieldset',
      '#title' => 'Details of lab stock',
      'advanced' => array(
        '#type' => 'fieldset',
        '#weight' => 15,
        '#title' => 'Advanced',
        'resource_db' => array(
          '#weight' => 0,
          '#type' => 'textfield',
          '#title' => 'Stock center name',
          '#description' => t('The name of the stock center or repository which held/holds the material used (e.g. Bayer). Just leave it empty if none is used.'),
          '#size' => 17,
          '#autocomplete_path' => 'genes4all/experiment/dbname/autocomplete',
        ),
        'resource_dbxref' => array(
          '#weight' => 1,
          '#type' => 'textfield',
          '#title' => 'Stock ID',
          '#description' => t('The ID given by the above repository.'),
          '#size' => 17,
        ),
      ),
      'colony_infection' => array(
        '#weight' => 2,
        '#title' => 'Any known colony infection',
        '#type' => 'checkboxes',
        '#options' => $infection_array,
      ),
    ),
    'indiv_bioassay' => array(
      '#type' => 'fieldset',
      '#title' => "Details regarding individuals to which RNAi agent was delivered",
      '#collapsible' => TRUE,
      '#weight' => 4,
      'origin' => array(
        '#type' => 'radios',
        '#options' => $origin_array,
        '#required' => TRUE,
        '#title' => 'Origin of used material',
        '#weight' => 0,
      ),
      'geolocation' => array(
        '#type' => 'textfield',
        '#title' => t('Geolocation if field-collected'),
        '#description' => l(
          t('GPS coordinates in Hours, Minutes, Seconds if known'), 'http://www.fcc.gov/mb/audio/bickel/DDDMMSS-decimal.html',
          //TODO go to advanced_help
          array('attributes' => array(
              'title' => t('How to find GPS coordinates from Google Maps: A) Find your location via Google Maps. B) Centre on the location, Right mouse on the exact spot and choose "Center map here" C) Click on "Link to this Page". D) This will display a URL. Embedded in the URL you will see "ll=". Immediately following that are the coordinates, in this example 51.500795,-0.142264. The first is N (S if negative) the second is E (W if negative). If you want to convert these coordinates to degrees, minutes and seconds, click to go the Federal Communications Commission converter'),
              'target' => '_blank',
            ))
        ),
        '#weight' => 1,
        '#size' => 16,
      ),
      'dev_stage' => array(
        '#title' => 'Developmental stage',
        '#type' => 'radios',
        '#options' => $stage_array,
        '#required' => TRUE,
        '#weight' => 2,
      ),
      'indiv_infection' => array(
        '#weight' => 3,
        '#title' => 'Known infection on studied individuals',
        '#type' => 'checkboxes',
        '#options' => $infection_array,
      ),
      'pretreatment_protocol' => array(
        '#weight' => 5,
        '#title' => 'Pretreatment protocol',
        '#description' => t('If animals where pretreated in any way (anti-infection agent, starved, induced), please give the full protocol.'),
        '#type' => 'textarea',
      ),
    ),
    'taxonomy' => array(
      '#type' => 'fieldset',
      '#title' => "Taxonomy",
      '#description' => t('Species from which this data originates'),
      '#collapsible' => TRUE,
      '#weight' => 10,
      'tax_class' => array(
        '#type' => 'textfield',
        '#title' => t('Class'),
        '#required' => TRUE,
        '#weight' => -2,
        '#size' => 20,
        '#default_value' => 'e.g. Insecta',
      ),
      'tax_order' => array(
        '#type' => 'textfield',
        '#title' => t('Order'),
        '#required' => TRUE,
        '#weight' => -1,
        '#size' => 20,
        '#default_value' => 'e.g. Lepidoptera',
      ),
      'tax_family' => array(
        '#type' => 'textfield',
        '#title' => t('Family'),
        '#required' => TRUE,
        '#weight' => 0,
        '#size' => 20,
        '#default_value' => 'e.g. Noctuidae',
      ),
      'genus' => array(
        '#type' => 'textfield',
        '#title' => t('Genus'),
        '#required' => TRUE,
        '#weight' => 1,
        '#size' => 20,
      ),
      'species' => array(
        '#type' => 'textfield',
        '#title' => t('Species'),
        '#description' => t('the second epithet/attribute of the binomial latin name'),
        '#required' => TRUE,
        '#weight' => 2,
        '#size' => 20,
      ),
      'ncbi_taxid' => array(
        '#type' => 'textfield',
        '#title' => t('NCBI Taxonomy ID'),
        '#description' => t('Obtainable from the NCBI ') . l('Taxonomy database', 'http://www.ncbi.nlm.nih.gov/Taxonomy/', array(
            'attributes' => array('target' => '_blank'),
          )
        ),
        '#weight' => 4,
        '#size' => 8,
      ),
    ),
  );



  $form[] = array(
    '#weight' => -10,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Store'),
      '#weight' => 0,
    ),
    'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'button',
      '#value' => t('Reset'),
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 1,
    ),
  );
  $form[] = array(
    '#weight' => 25,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Store'),
      '#weight' => 0,
    ),
    'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'button',
      '#value' => t('Reset'),
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 1,
    ),
  );

  return $form;
}

function genes4all_experiment_animals_form_validate($form, &$form_state) {
  $ncbi_taxid = $form_state['values']['ncbi_taxid'];
  $term_error = t('A taxonomic term cannot have spaces in it!');
  if (strpos($form_state['values']['tax_class'], ' ') !== FALSE) {
    form_set_error('tax_class', $term_error);
  }
  if (strpos($form_state['values']['tax_order'], ' ') !== FALSE) {
    form_set_error('tax_order', $term_error);
  }
  if (strpos($form_state['values']['tax_family'], ' ') !== FALSE) {
    form_set_error('tax_family', $term_error);
  }
  if (strpos($form_state['values']['genus'], ' ') !== FALSE) {
    form_set_error('genus', $term_error);
  }
  if (strpos($form_state['values']['species'], ' ') !== FALSE) {
    form_set_error('species', $term_error);
  }
  $taxa = array(
    'tax_class' => $form_state['values']['tax_class'],
    'tax_order' => $form_state['values']['tax_order'],
    'tax_family' => $form_state['values']['tax_family'],
    'genus' => $form_state['values']['genus'],
    'species' => $form_state['values']['species'],
  );
  foreach ($taxa as $key => $value) {
    if (!ctype_alpha($value)) {
      form_set_error($key, t('Highlighted taxon does not seem to be composed of only alphabetic characters.'));
    }
  }
  if (!empty($ncbi_taxid)) {
    if (!is_numeric($ncbi_taxid)) {
      form_set_error("ncbi_taxid", t('NCBI Taxonomy ID is not valid.'));
    }
  }
}

function genes4all_experiment_animals_form_submit($form, &$form_state) {
  // create fake name
  $form_state['values']['animals_name'] = $_SESSION['study_name'] .' experimental animals';
  $resource_id = genes4all_experiment_resource_add($form_state['values'], $form_state['values']['animals_name'], 'experimental animals');
  if (!empty($resource_id)) {
    $_SESSION['animals_id'] = $resource_id;
    $_SESSION['animals_name'] = $form_state['values']['animals_name'];
  }
  else {
    drupal_set_message(t('I was not able to store the resource data.'), 'error', FALSE);
    gmod_dbsf_form_clear($form, &$form_state);
    return FALSE;
  }
}

function genes4all_experiment_delivery_populate_form($form_state) {
  $data = genes4all_experiment_resource_get('', 'delivery protocol');
  if (empty($data)) {
    $data = array('0' => 'None saved by your passkey');
  }
  $form = array(
    'existing' => array(
      '#type' => 'fieldset',
      '#title' => 'Existing delivery protocols',
      'delivery_id' => array(
        '#type' => 'select',
        '#options' => $data,
      ),
      'populate' => array(
        '#type' => 'submit',
        '#value' => 'Populate',
      ),
    ),
  );

  return $form;
}

function genes4all_experiment_delivery_populate_form_submit($form, &$form_state) {
  $_SESSION['delivery_id'] = $form_state['values']['delivery_id'];
}

function genes4all_experiment_delivery_form($form_state) {
  $form = array();
  $delivery_method_array = array(
    'injection' => 'Injection',
    'feeding' => 'Feeding',
  );
  $control_method_array = array(
    'non_specific_dsRNA' => 'Non specific dsRNA',
    'buffer' => 'Buffer',
  );

  $form = array(
    'delivery_protocol' => array(
      '#type' => 'fieldset',
      '#title' => 'Delivery protocol',
      '#weight' => 1,
      'construct_concentration' => array(
        '#type' => 'textfield',
        '#title' => t('Construct working concentration (µg/µl)'),
        '#description' => t('This is in order to estimate the amount of the pure RNAi construct. Please give 0 if it is unknown.'),
        '#weight' => 1,
        '#size' => 5,
        '#required' => TRUE,
      ),
      'delivery_organ' => array(
        '#type' => 'textfield',
        '#title' => t('Organ/tissue construct was delivered to'),
        '#description' => t('Provide the organ or tissue (fat-body, hemocytes, muscles, wings, eyes, antenna, etc.) to which the construct was directly delivered. If feeding (or other whole animal approach) was used as the delivery method, then write feeding or other approach. DO NOT give your protocol here.'),
        '#weight' => 2,
        '#required' => TRUE,
      ), 'delivery_detail' => array(
        '#type' => 'textarea',
        '#title' => t('Detailed delivery protocol'),
        '#description' => t('Optionally you may give a detailed protocol used for delivery here.'),
        '#weight' => 7,
      ),
      'delivery_method' => array(
        '#type' => 'checkboxes',
        '#title' => t('Delivery method'),
        '#weight' => 2,
        '#options' => $delivery_method_array,
        '#required' => TRUE,
      ),
      'delivery_buffer' => array(
        '#type' => 'textfield',
        '#title' => t('Delivery buffer'),
        '#weight' => 3,
        '#size' => 15,
      ),
      'adjuvant_name' => array(
        '#type' => 'textfield',
        '#title' => t('Adjuvant/carrier'),
        '#weight' => 4,
        '#size' => 15,
      ),
      'adjuvant_amount' => array(
        '#type' => 'textfield',
        '#title' => t('Estimated delivery amount (µl/mg of insect)'),
        '#description' => t('Delivery amount of with the buffer/adjuvant/carrier. Please give 0 if it is really unknown, but this violates the Minimum Information Criteria for publication.'),
        '#size' => 5,
        '#required' => TRUE,
      ),
      'replicates' => array(
        '#type' => 'textfield',
        '#title' => t('Number of replicates'),
        '#weight' => 6,
        '#size' => 5,
        '#default_value' => 0,
        '#required' => TRUE,
      ),
    ), 'control_set' => array(
      '#type' => 'fieldset',
      '#title' => 'Quality control',
      '#weight' => 2,
      'control_method' => array(
        '#type' => 'checkboxes',
        '#title' => t('Control'),
        '#options' => $control_method_array,
      ),
      'control_method_detail' => array(
        '#type' => 'textarea',
        '#title' => t('Control protocol'),
      ),
    ),
  );

  $form[] = array(
    '#weight' => -10,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Store'),
      '#weight' => 0,
    ),
    'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'button',
      '#value' => t('Reset'),
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 1,
    ),
  );
  $form[] = array(
    '#weight' => 25,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Store'),
      '#weight' => 0,
    ),
    'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'button',
      '#value' => t('Reset'),
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 1,
    ),
  );
  return $form;
}

function genes4all_experiment_delivery_form_validate($form, &$form_state) {
  $constr_conc = $form_state['values']['construct_concentration'];
  $adj_amnt    = $form_state['values']['adjuvant_amount'];
  $replicates  = $form_state['values']['replicates'];
  if (!empty($constr_conc) && !is_numeric($constr_conc)) {
    form_set_error('construct_concentration', t('Construct concentration must be a number.'));
  }
  if (!empty($adj_amnt) && !is_numeric($adj_amnt)) {
    form_set_error('adjuvant_amount', t('Delivered amount must be a number.'));
  }
  if (!empty($replicates) && !is_numeric($replicates)) {
    form_set_error('replicates', t('Replicates must be a number.'));
  }
}

function genes4all_experiment_delivery_form_submit($form, &$form_state) {
  // create fake name
  $form_state['values']['delivery_name'] = $_SESSION['study_name'] .' delivery protocol';
  $resource_id = genes4all_experiment_resource_add($form_state['values'], $form_state['values']['delivery_name'], 'delivery protocol');
  if (!empty($resource_id)) {
    $_SESSION['delivery_id'] = $resource_id;
    $_SESSION['delivery_name'] = $form_state['values']['delivery_name'];
  }
  else {
    drupal_set_message(t('I was not able to store the resource data.'), 'error', FALSE);
    gmod_dbsf_form_clear($form, &$form_state);
    return FALSE;
  }
}

function genes4all_experiment_assay_populate_form($form_state) {
  $data = genes4all_experiment_resource_get('', 'assay protocol');
  if (empty($data)) {
    $data = array('0' => 'None saved by your passkey');
  }
  $form = array(
    'existing' => array(
      '#type' => 'fieldset',
      '#title' => 'Existing assay protocols',
      'assay_id' => array(
        '#type' => 'select',
        '#options' => $data,
      ),
      'populate' => array(
        '#type' => 'submit',
        '#value' => 'Populate',
      ),
    ),
  );
  return $form;
}

function genes4all_experiment_assay_populate_form_submit($form, &$form_state) {
  $_SESSION['assay_id'] = $form_state['values']['assay_id'];
}

function genes4all_experiment_assay_form($form_state) {
  $form = array();
  // TODO: the friendly names should go as definitions into the cvterm table?
  $stage_array = array(
    'unfertilized egg' => 'Unfertilized egg',
    'embryo' => 'Embryo',
    'larva' => 'Larva',
    'pupa' => 'Pupa',
    'adult' => 'Adult',
  );
  $detection_method_array = array(
    'qPCR' => 'Quantitative-PCR',
    'RT-PCR (semi-quantitative)' => 'RT-PCR (semi-quantitative)',
    'N.blot' => 'Northern-blot',
    'W.blot' => 'Western-blot',
    'enzymatic_activity' => 'Enzymatic activity',
    'phenotype' => 'Phenotype',
  );
  $silencing_levels_array = array(
    'none' => 'No silencing',
    'low' => 'Low silencing',
    'high' => 'High silencing',
  );
  $dsrna_processing_array = array(
    'not checked' => 'Not checked',
    'undetected' => 'Not detected',
    'detected' => 'Detected',
  );

  $form = array(
    'detection' => array(
      '#type' => 'fieldset',
      '#title' => 'Detection',
      '#collapsible' => TRUE,
      'detection_method' => array(
        '#type' => 'checkboxes',
        '#title' => t('Detection method'),
        '#options' => $detection_method_array,
        '#required' => TRUE,
      ),
      'detect_time' => array(
        '#type' => 'textfield',
        '#title' => t('Detection time (hours post-delivery)'),
        '#description' => t('You must provide the time you did the detection. If you really do not know, provide 0.'),
        '#weight' => 1,
        '#size' => 5,
        '#required' => TRUE,
      ),
      'dev_stage' => array(
        '#title' => 'Developmental stage at detection',
        '#type' => 'radios',
        '#options' => $stage_array,
        '#required' => TRUE,
      ), 'assay_organ' => array(
        '#type' => 'textfield',
        '#title' => t('Organ/tissue that was assayed for changes in gene expression'),
        '#description' => t('Provide the organ or tissue (e.g. fat-body, hemocytes, muscles, wings, eyes, antenna, etc.) or state it was the whole organism.'),
        '#weight' => 6,
        '#required' => TRUE,
      ),
    ),
    'result_set' => array(
      '#type' => 'fieldset',
      '#title' => 'Study results',
      '#collapsible' => TRUE,
      'silencing_level' => array(
        '#title' => 'Gene silencing levels relative to control',
        '#type' => 'radios',
        '#options' => $silencing_levels_array,
        '#required' => TRUE,
      ),
      'silencing_accurate_level' => array(
        '#title' => 'Percentage of silencing',
        '#type' => 'textfield',
        '#description' => t('If you used an accurate method of estimating silencing level, please provide the percentage of silencing here.'),
        '#size' => 5,
      ),
      'detect_dsrna_processing' => array(
        '#title' => 'Detection of dsRNA processing (siRNA formation)',
        '#type' => 'radios',
        '#options' => $dsrna_processing_array,
        '#required' => TRUE,
      ),
    ),
  );

  $form[] = array(
    '#weight' => -10,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Store'),
      '#weight' => 0,
    ),
    'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'button',
      '#value' => t('Reset'),
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 1,
    ),
  );
  $form[] = array(
    '#weight' => 25,
    'submit' => array(
      '#type' => 'submit',
      '#value' => t('Store'),
      '#weight' => 0,
    ),
    'clear' => array(
      '#attributes' => array('title' => t('Clear the form')),
      '#type' => 'button',
      '#value' => t('Reset'),
      '#validate' => array('gmod_dbsf_form_clear'),
      '#weight' => 1,
    ),
  );

  return $form;
}

function genes4all_experiment_assay_form_validate($form, &$form_state) {
  $time = $form_state['values']['detect_time'];
  $level = $form_state['values']['silencing_accurate_level'];

  if (!empty($time) && !is_numeric($time)) {
    form_set_error('detect_time', t('Detection time must be a number (hours).'));
  }
  if (!empty($level) && !is_numeric($level)) {
    form_set_error('silencing_accurate_level', t('Silencing level must be just a number (e.g. no percentage sign).'));
  }
}

function genes4all_experiment_assay_form_submit($form, &$form_state) {
  // create fake name
  $form_state['values']['assay_name'] = $_SESSION['study_name'] .' assay protocol';
  $resource_id = genes4all_experiment_resource_add($form_state['values'], $form_state['values']['assay_name'], 'assay protocol');
  if (!empty($resource_id)) {
    $_SESSION['assay_id'] = $resource_id;
    $_SESSION['assay_name'] = $form_state['values']['assay_name'];
  }
  else {
    drupal_set_message(t('I was not able to store the resource data.'), 'error', FALSE);
    gmod_dbsf_form_clear($form, &$form_state);
    return FALSE;
  }
}

function genes4all_experiment_study_finished_form($form_state) {
  $form = array();
  if (empty($_SESSION['study_id']) || empty($_SESSION['study_type']) || empty($_SESSION['passkey']) || empty($_SESSION['pub_id']
    ) || empty($_SESSION['target_id']) || empty($_SESSION['construct_id']) || empty($_SESSION['animals_id'])
    || empty($_SESSION['delivery_id']) || empty($_SESSION['assay_id'])
  ) {

    $form[] = array(
      '#weight' => -10,
      'submit' => array(
        '#type' => 'submit',
        '#value' => t('Not ready yet'),
        '#disabled' => TRUE,
        '#weight' => 0,
      ),
    );
  }
  else {
    $form[] = array(
      '#weight' => -10,
      'submit' => array(
        '#type' => 'submit',
        '#value' => t('Finalize'),
        '#weight' => 0,
      ),
    );
  }
  return $form;
}

function genes4all_experiment_study_finished_form_validate($form, &$form_state) {
  if (empty($_SESSION['study_id']) || empty($_SESSION['study_type']) || empty($_SESSION['passkey']) || empty($_SESSION['pub_id']
    ) || empty($_SESSION['target_id']) || empty($_SESSION['construct_id']) || empty($_SESSION['animals_id'])
    || empty($_SESSION['delivery_id']) || empty($_SESSION['assay_id'])
  ) {
    form_set_error('', t('You have not filled out all the required forms.'));
  }
}

function genes4all_experiment_study_finished_form_submit($form, &$form_state) {
  $study_id  = $_SESSION['study_id'];
  $type_name = $_SESSION['study_type'];
  $passkey   = $_SESSION['passkey'];
  $data      = array(
    'pub_id' => $_SESSION['pub_id'],
    'target_id' => $_SESSION['target_id'],
    'construct_id' => $_SESSION['construct_id'],
    'animals_id' => $_SESSION['animals_id'],
    'delivery_id' => $_SESSION['delivery_id'],
    'assay_id' => $_SESSION['assay_id'],
  );

  $success = genes4all_experiment_study_finished($study_id, $type_name, $passkey, $data);
  if (!empty($success === TRUE)) {
    drupal_goto('genes4all/experiment/added');
  }
  else {
    drupal_set_message(t('I was not able to store your experiment.'), 'error', FALSE);
    gmod_dbsf_form_clear($form, &$form_state);
    return FALSE;
  }
}

function genes4all_experiment_study_finished($study_id, $type_name, $passkey, $data = NULL) {
  if (empty($study_id) || empty($type_name) || empty($passkey)) {
    return FALSE;
  }
  $approved_types = array('RNAi');
  if (!in_array($type_name, $approved_types)) {
    drupal_set_message(t('%type_name is not an approved type.', array('%type_name' => $type_name)), 'error', FALSE);
    return FALSE;
  }
  $type_id = gmod_dbsf_get_add_cv_withnames('study_type', $type_name);

  //security: check if id is real and get existing type
  $check_sql      = "SELECT study_id,type_id from {gmod_dbsf_study} WHERE study_id=$study_id and passkey='$passkey'";
  $res_check      = db_fetch_array(db_query($check_sql));
  $study_id       = $res_check['study_id'];
  $type_id_stored = $res_check['type_id'];
  $update_sql     = "UPDATE {gmod_dbsf_study} set %s='%s' where study_id=$study_id";
  if (empty($study_id)) {
    return FALSE;
  }
  if ($type_id != $type_id_stored) {
    db_query($update_sql, 'type_id', $type_id);
  }
  // now add the various data.
  $pub_id = $data['pub_id'];
  $target_id = $data['target_id'];
  $construct_id = $data['construct_id'];
  $experimental_animal_id = $data['animals_id'];
  $delivery_protocol_id = $data['delivery_id'];
  $assay_protocol_id = $data['assay_id'];
  if (empty($target_id) || empty($pub_id) || empty($construct_id) || empty($experimental_animal_id)
    || empty($delivery_protocol_id) || empty($assay_protocol_id)
  ) {
    return FALSE;
  }
  // pub
  $update_sql = "UPDATE {gmod_dbsf_study} set pub_id=$pub_id where study_id=$study_id";
  db_query($update_sql);
  //features
  $select_feature_sql = "SELECT feature_id from {gmod_dbsf_study_feature} where study_id=$study_id AND feature_id IN " . "(SELECT feature_id from {gmod_dbsf_feature} as feature JOIN {gmod_dbsf_cvterm} as cvterm ON cvterm.cvterm_id=feature.type_id " . " where cvterm.name='%s')";
  $delete_feature_sql = "DELETE FROM {gmod_dbsf_study_feature} where study_id=$study_id and feature_id IN " . "(SELECT feature_id from {gmod_dbsf_feature} as feature JOIN {gmod_dbsf_cvterm} as cvterm ON cvterm.cvterm_id=feature.type_id " . " where cvterm.name='%s')";
  $insert_feature_sql = "INSERT INTO {gmod_dbsf_study_feature} (study_id,feature_id) VALUES ($study_id,%d)";
  $res                = db_fetch_array(db_query($select_feature_sql, 'gene_target'));
  if (empty($res['feature_id']) || $res['feature_id'] != $target_id) {
    db_query($delete_feature_sql, 'gene_target');
    db_query($insert_feature_sql, $target_id);
  }
  unset($res);
  $res = db_fetch_array(db_query($select_feature_sql, 'rnai_construct'));
  if (empty($res['feature_id']) || $res['feature_id'] != $construct_id) {
    db_query($delete_feature_sql, 'rnai_construct');
    db_query($insert_feature_sql, $construct_id);
  }
  unset($res);


  $select_resource_sql = "SELECT resource_id from {gmod_dbsf_study_resource} where study_id=$study_id AND resource_id IN " . "(SELECT resource_id from {gmod_dbsf_resource} as resource JOIN {gmod_dbsf_cvterm} as cvterm ON cvterm.cvterm_id=resource.type_id " . " where cvterm.name='%s')";
  $delete_resource_sql = "DELETE FROM {gmod_dbsf_study_resource} where study_id=$study_id and resource_id IN " . "(SELECT resource_id from {gmod_dbsf_resource} as resource JOIN {gmod_dbsf_cvterm} as cvterm ON cvterm.cvterm_id=resource.type_id " . " where cvterm.name='%s')";
  $insert_resource_sql = "INSERT INTO {gmod_dbsf_study_resource} (study_id,resource_id) VALUES ($study_id,%d)";

  $res = db_fetch_array(db_query($select_resource_sql, 'experimental animals'));
  if (empty($res['resource_id']) || $res['resource_id'] != $experimental_animal_id) {
    db_query($delete_resource_sql, 'experimental animals');
    db_query($insert_resource_sql, $experimental_animal_id);
  }
  unset($res);
  $res = db_fetch_array(db_query($select_resource_sql, 'delivery protocol'));
  if (empty($res['resource_id']) || $res['resource_id'] != $delivery_protocol_id) {
    db_query($delete_resource_sql, 'delivery protocol');
    db_query($insert_resource_sql, $delivery_protocol_id);
  }
  unset($res);
  $res = db_fetch_array(db_query($select_resource_sql, 'assay protocol'));
  if (empty($res['resource_id']) || $res['resource_id'] != $assay_protocol_id) {
    db_query($delete_resource_sql, 'assay protocol');
    db_query($insert_resource_sql, $assay_protocol_id);
  }
  unset($res);
  //$res= db_fetch_array(db_query($select_resource_sql,'result set'));
  //if (empty($res['resource_id']) || $res['resource_id'] !=$result_set_id){
  //  db_query($delete_resource_sql,'result set');
  //  db_query($insert_resource_sql,$result_set_id);
  //}unset($res);

  return TRUE;
}



/*
 * return array with study_id pub_name, name for each resource & feature
 */
function genes4all_experiment_study_add($data) {
  if (empty($data) || empty($data['study_name'])) {
    return FALSE;
  }
  $return_array    = array();
  $submitter_email = check_plain(trim($data['submitter_email']));
  $uname           = check_plain(trim($data['study_name']));
  $db              = 'local study';
  $description     = check_plain(trim($data['study_comments']));
  $accession       = $uname;
  $db              = strtolower(check_plain(trim($db)));
  $passkey         = check_plain(trim($data['passkey']));
  $dbxref_id       = gmod_dbsf_get_add_dbxref_withnames($db, $accession);

  $check_sql = "SELECT study_id,description,type_id,submitter_email from {gmod_dbsf_study} WHERE uniquename='%s'";
  $res_check = db_fetch_array(db_query($check_sql, $uname));
  $study_id = $res_check['study_id'];
  $type_id = $res_check['type_id'];
  $description_stored = $res_check['description'];
  $submitter_email_stored = $res_check['submitter_email'];
  $update_sql = "UPDATE {gmod_dbsf_study} set %s='%s' where study_id=%d";

  if (empty($study_id)) {
    $insert_sql = "INSERT into {gmod_dbsf_study} (submitter_email,uniquename,dbxref_id,passkey,type_id,pub_id) VALUES ('%s','%s','%s','%s',"
    . '(SELECT cvterm_id from {gmod_dbsf_cvterm} cvterm JOIN {gmod_dbsf_cv} cv ON cv.cv_id=cvterm.cv_id '
    . "where cvterm.name='unfinished' AND cv.name='study_type'),1)";
    db_query($insert_sql,$submitter_email, $uname, $dbxref_id, $passkey);
    $res_check = db_fetch_array(db_query($check_sql, $uname));
    $study_id = $res_check['study_id'];
    if (empty($study_id)) {
      return FALSE;
    }
    $return_array['study_id'] = $study_id;
  }
  else {
    $check_pass = "SELECT passkey,pub_id,type_id from {gmod_dbsf_study} WHERE study_id=$study_id";
    $check_res = db_fetch_array(db_query($check_pass));
    if ($check_res['passkey'] !== $passkey) {
      return FALSE;
    }
    $pub_id = $check_res['pub_id'];
    $return_array['study_id'] = $study_id;
    $return_array['study_status'] = $type_id;
    if ($pub_id > 1) {
      $study_pub_select = "SELECT uniquename from {gmod_dbsf_pub} where pub_id=$pub_id";
      $res = db_fetch_array(db_query($study_pub_select));
      $return_array['pub_name'] = $res['uniquename'];
      $return_array['pub_id'] = $pub_id;
      $type_select = "SELECT name from {gmod_dbsf_cvterm} where cvterm_id=$type_id";
      $type_res = db_fetch_array(db_query($type_select));
      $return_array['finished'] = $type_res['name'];
      $study_feature_select = 'SELECT feature.feature_id,feature.uniquename,cvterm.name as type_name from {gmod_dbsf_feature} as feature JOIN {gmod_dbsf_cvterm} as cvterm on feature.type_id=cvterm.cvterm_id '. " where feature.feature_id IN (SELECT feature_id from {gmod_dbsf_study_feature} as study_feature where study_feature.study_id=$study_id)";
      $study_resource_select = 'SELECT resource.resource_id,resource.uniquename,cvterm.name as type_name from {gmod_dbsf_resource} as resource JOIN {gmod_dbsf_cvterm} as cvterm on resource.type_id=cvterm.cvterm_id '. " where resource.resource_id IN (SELECT resource_id from {gmod_dbsf_study_resource} as study_resource where study_resource.study_id=$study_id)";
      $feat_res = db_query($study_feature_select);
      while ($feat_data = db_fetch_array($feat_res)) {
        $return_array[$feat_data['type_name']]['name'] = $feat_data['uniquename'];
        $return_array[$feat_data['type_name']]['id'] = $feat_data['feature_id'];
      }
      $res_res = db_query($study_resource_select);
      while ($res_data = db_fetch_array($res_res)) {
        $return_array[$res_data['type_name']]['name'] = $res_data['uniquename'];
        $return_array[$res_data['type_name']]['id'] = $res_data['resource_id'];
      }
    }
  }

  if (!empty($description) && $description != $description_stored) {
    db_query($update_sql, 'description', $description, $study_id);
  }
  elseif (!empty($submitter_email) && $submitter_email != $submitter_email_stored) {
    db_query($update_sql, 'submitter_email', $submitter_email, $study_id);
  }

  return ($return_array);
}

/**
 *
 * Add data to pub table(s) if one does not exist already.
 *
 * @return $pub_id string integer
 */
function genes4all_experiment_publication_add($data = NULL) {
  //dpm($_SESSION);
  //dpm($data);
  $passkey = check_plain($_SESSION['passkey']);
  $array = array();

  if (empty($data['pub_uniquename']) && empty($data['pub_id'])) {
    // if empty return array of all pubs.
    $select = "SELECT uniquename,passkey,pub_id from {gmod_dbsf_pub}";
    $res = db_query($select);
    while ($row = db_fetch_array($res)) {
      if ($row['passkey'] == $passkey) {
        $array[$row['pub_id']] = $row['uniquename'];
      }
    }
    return $array;
  }
  elseif (!empty($data['pub_id']) && is_numeric($data['pub_id'])) {
    return $data['pub_id'];
  }
  elseif (!empty($data['pub_uniquename'])) {
    // or add publication
    $uname     = check_plain(trim($data['pub_uniquename']));
    $db        = $data['pub_db'] ? $data['pub_db'] : 'local publication';
    $accession = $data['pub_db_accession'] ? $data['pub_db_accession'] : 'Unpublished';
    $accession = check_plain(trim($accession));
    $db        = strtolower(check_plain(trim($db)));

    $dbxref_id_user = gmod_dbsf_get_add_dbxref_withnames($db, $accession);
    $check_pub      = "SELECT pub_id,dbxref_id from {gmod_dbsf_pub} WHERE uniquename='$uname'";
    $res_pub        = db_fetch_array(db_query($check_pub));
    $pub_id         = $res_pub['pub_id'];
    $dbxref_id      = $res_pub['dbxref_id'];

    if (empty($pub_id)) {
      $insert_pub = "INSERT into {gmod_dbsf_pub} (uniquename,dbxref_id,passkey) VALUES ('$uname',$dbxref_id_user,'$passkey')";
      db_query($insert_pub);
      $res_pub = db_fetch_array(db_query($check_pub));
      $pub_id = $res_pub['pub_id'];
    }
    else {
      if (!empty($dbxref_id_user) && $dbxref_id_user != $dbxref_id) {
        $update_pub = "UPDATE {gmod_dbsf_pub} SET dbxref_id=$dbxref_id_user WHERE pub_id=$pub_id";
        db_query($update_pub);
      }
    }
    // we now have pub_id.
    return $pub_id;
  }
}

function genes4all_experiment_resource_add($data, $uname, $type) {
  $passkey = check_plain($_SESSION['passkey']);
  if (empty($data) || empty($uname) || empty($type)) {
    return FALSE;
  }
  $approved_types = array('experimental animals', 'result set', 'assay protocol', 'delivery protocol');
  //dpm($data);
  if (!in_array($type, $approved_types)) {
    drupal_set_message(t('%type is not an approved type.', array('%type' => $type)), 'error', FALSE);
    return FALSE;
  }

  $cvterm_check = "SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . "WHERE cv.name='%s' AND cvterm.name='%s'";
  $cvterm_insert = "INSERT INTO {gmod_dbsf_cvterm} (name,cv_id) VALUES ('%s',(SELECT cv_id from {gmod_dbsf_cv} where name='%s'))";
  $resource_cvterm_insert = "INSERT INTO {gmod_dbsf_resource_cvterm} (resource_id,cvterm_id) VALUES (%d," . "(SELECT cvterm_id FROM {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " WHERE cvterm.name='%s' AND cv.name='%s'))";
  $resource_cvterm_update = "UPDATE {gmod_dbsf_resource_cvterm} set cvterm_id=%d WHERE resource_id=%d";
  $resource_cvterm_delete = "DELETE from {gmod_dbsf_resource_cvterm} WHERE resource_id=%d " . "AND cvterm_id=(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN " . " {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id where cvterm.name='%s' AND cv.name='%s')";
  $resource_cvterm_check = "SELECT cvterm.name,cv.name from {gmod_dbsf_resource_cvterm} as resource_cvterm " . " JOIN {gmod_dbsf_cvterm} cvterm ON resource_cvterm.cvterm_id=cvterm.cvterm_id JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " where resource_id=%d";
  $resource_update = "UPDATE {gmod_dbsf_resource} as resource set %s='%s' where resource_id=%d and passkey='$passkey'";
  $resourceprop_check = "SELECT resourceprop_id from {gmod_dbsf_resourceprop} as resourceprop WHERE resource_id=%d and type_id=%d";
  $resourceprop_update = "UPDATE {gmod_dbsf_resourceprop} SET value='%s' WHERE resource_id=%d and type_id=%d";
  $resourceprop_insert = "INSERT INTO {gmod_dbsf_resourceprop} (resource_id,type_id,value) VALUES (%d,%d,'%s')";

  $uname = check_plain(trim($uname));
  // STORED data
  $resource_data = genes4all_experiment_resource_get($uname);
  $resource_id   = $resource_data['general']['resource_id'];
  $cvterms       = $resource_data['cvterm'];

  // FROM user form
  $db = $data['resource_db'] ? $data['resource_db'] : 'local resource';
  $accession = $data['resource_dbxref'] ? $data['resource_dbxref'] : $uname;

  $db         = strtolower(check_plain(trim($db)));
  $accession  = check_plain(trim($accession));
  $tax_order  = '';
  $genus      = '';
  $tax_family = '';
  $tax_class  = '';
  $ncbi_taxid = '';
  if (!empty($data['tax_order'])) {
    $tax_class          = check_plain(trim($data['tax_class']));
    $data['tax_class']  = $tax_class;
    $tax_order          = check_plain(trim($data['tax_order']));
    $data['tax_order']  = $tax_order;
    $tax_family         = check_plain(trim($data['tax_family']));
    $data['tax_family'] = $tax_family;
    $genus              = check_plain(trim($data['genus']));
    $data['genus']      = $genus;
    $species            = check_plain(trim($data['species']));
    $data['species']    = $species;
    $ncbi_taxid         = check_plain(trim($data['ncbi_taxid']));
    $data['ncbi_taxid'] = $ncbi_taxid;
  }
  $organism_id = genes4all_experiment_get_add_taxonomy($data);
  if (empty($organism_id)) {
    $organism_id = 1;
    //general one
  }
  $dbxref_id = gmod_dbsf_get_add_dbxref_withnames($db, $accession);

  if (!empty($resource_id)) {
    // update if needed
    if ($dbxref_id != $resource_data['general']['dbxref_id']) {
      db_query($resource_update, 'dbxref_id', $dbxref_id, $resource_id);
    }
    if ($organism_id != $resource_data['general']['organism_id']) {
      db_query($resource_update, 'organism_id', $organism_id, $resource_id);
    }
    if ($type == 'experimental animals') {
      genes4all_experiment_resource_add_animal($resource_id, $data, $cvterms);
    }
    elseif ($type == 'delivery protocol') {
      genes4all_experiment_resource_add_delivery_protocol($resource_id, $data, $cvterms);
    }
    elseif ($type == 'assay protocol') {
      genes4all_experiment_resource_add_assay_protocol($resource_id, $data, $cvterms);
    }
    elseif ($type == 'result set') {
      genes4all_experiment_resource_add_result_set($resource_id, $data, $cvterms);
    }
  }
  else {
    // create it
    $db_res          = db_fetch_array(db_query($cvterm_check, 'resource_type', $type));
    $type_id         = $db_res['cvterm_id'];
    $resource_insert = "INSERT INTO {gmod_dbsf_resource} (uniquename,dbxref_id,type_id,organism_id,passkey) ". " VALUES ('$uname',$dbxref_id,$type_id,$organism_id,'$passkey')";
    db_query($resource_insert);
    $resource_check = "SELECT resource_id from {gmod_dbsf_resource} where uniquename='$uname'";
    $db_res         = db_fetch_array(db_query($resource_check));
    $resource_id    = $db_res['resource_id'];
    if (empty($resource_id)) {
      drupal_set_message(t('I was not able to store the resource data at this time.'), 'error', FALSE);
      return FALSE;
    }
    // now add the seq/cvterms
    if ($type == 'experimental animals') {
      genes4all_experiment_resource_add_animal($resource_id, $data, $cvterms);
    }
    elseif ($type == 'delivery protocol') {
      genes4all_experiment_resource_add_delivery_protocol($resource_id, $data, $cvterms);
    }
    elseif ($type == 'assay protocol') {
      genes4all_experiment_resource_add_assay_protocol($resource_id, $data, $cvterms);
    }
    elseif ($type == 'result set') {
      genes4all_experiment_resource_add_result_set($resource_id, $data, $cvterms);
    }
  }
  return $resource_id;
}

function genes4all_experiment_feature_add($data, $uname, $type) {
  $passkey = check_plain($_SESSION['passkey']);
  if (empty($data) || empty($uname) || empty($type)) {
    return FALSE;
  }
  $approved_types = array('gene_target', 'rnai_construct');
  if (!in_array($type, $approved_types)) {
    drupal_set_message(t('%type is not an approved type.', array('%type' => $type)), 'error', FALSE);
    return FALSE;
  }

  $cvterm_check = "SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id
  " . "WHERE cv.name='%s' AND cvterm.name='%s'";
  $cvterm_insert = "INSERT INTO {gmod_dbsf_cvterm} (name,cv_id) VALUES ('%s',(SELECT cv_id from {gmod_dbsf_cv} where name='%s'))";
  $feature_cvterm_insert = "INSERT INTO {gmod_dbsf_feature_cvterm} (feature_id,cvterm_id) VALUES (%d," . "(SELECT cvterm_id FROM {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " WHERE cvterm.name='%s' AND cv.name='%s'))";
  $feature_cvterm_update = "UPDATE {gmod_dbsf_feature_cvterm} set cvterm_id=%d WHERE feature_id=%d";
  $feature_cvterm_delete = "DELETE from {gmod_dbsf_feature_cvterm} WHERE feature_id=%d " . "AND cvterm_id=(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN " . " {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id where cvterm.name='%s' AND cv.name='%s')";
  $feature_cvterm_check = "SELECT cvterm.name,cv.name from {gmod_dbsf_feature_cvterm} as feature_cvterm " . " JOIN {gmod_dbsf_cvterm} cvterm ON feature_cvterm.cvterm_id=cvterm.cvterm_id JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " where feature_id=%d";
  $feature_update = "UPDATE {gmod_dbsf_feature} as feature set %s='%s' where feature_id=%d and passkey='$passkey'";
  $featureprop_check = "SELECT featureprop_id from {gmod_dbsf_featureprop} as featureprop WHERE feature_id=%d and type_id=%d";
  $featureprop_update = "UPDATE {gmod_dbsf_featureprop} SET value='%s' WHERE feature_id=%d and type_id=%d";
  $featureprop_insert = "INSERT INTO {gmod_dbsf_featureprop} (feature_id,type_id,value) VALUES (%d,%d,'%s')";

  $uname = check_plain(trim($uname));
  // STORED data
  $feat_data  = genes4all_experiment_feature_get($uname);
  $feature_id = $feat_data['general']['feature_id'];
  $cvterms    = $feat_data['cvterm'];

  // FROM user form
  $db        = '';
  $accession = '';
  $seq       = '';

  if ($type == 'gene_target') {
    $db        = $data['target_db'] ? $data['target_db'] : 'local feature';
    $accession = $data['target_dbxref'] ? $data['target_dbxref'] : $uname;
    $seq       = check_plain(trim($data['target_seq']));
  }
  elseif ($type == 'rnai_construct') {
    $db        = $data['construct_db'] ? $data['construct_db'] : 'local feature';
    $accession = $data['construct_dbxref'] ? $data['construct_dbxref'] : $uname;
    $seq       = check_plain(trim($data['construct_seq']));
  }
  $db                 = strtolower(check_plain(trim($db)));
  $tax_class          = check_plain(trim($data['tax_class']));
  $data['tax_class']  = $tax_class;
  $accession          = check_plain(trim($accession));
  $tax_order          = check_plain(trim($data['tax_order']));
  $data['tax_order']  = $tax_order;
  $tax_family         = check_plain(trim($data['tax_family']));
  $data['tax_family'] = $tax_family;
  $genus              = check_plain(trim($data['genus']));
  $data['genus']      = $genus;
  $species            = check_plain(trim($data['species']));
  $data['species']    = $species;
  $ncbi_taxid         = check_plain(trim($data['ncbi_taxid']));
  $data['ncbi_taxid'] = $ncbi_taxid;

  $organism_id = genes4all_experiment_get_add_taxonomy($data);

  if (empty($organism_id)) {
    $organism_id = 1;
    //general one
  }
  $dbxref_id = gmod_dbsf_get_add_dbxref_withnames($db, $accession);

  if (!empty($feature_id)) {
    // update if needed
    if (!empty($seq) && $seq != $feat_data['general']['sequence']) {
      db_query($feature_update, 'residues', $seq, $feature_id);
    }
    if ($dbxref_id != $feat_data['general']['dbxref_id']) {
      db_query($feature_update, 'dbxref_id', $dbxref_id, $feature_id);
    }
    if ($organism_id != $feat_data['general']['organism_id']) {
      db_query($feature_update, 'organism_id', $organism_id, $feature_id);
    }
    if ($type == 'gene_target') {
      genes4all_experiment_feature_add_gene($feature_id, $data, $cvterms);
    }
    elseif ($type == 'rnai_construct') {
      genes4all_experiment_feature_add_construct($feature_id, $data, $cvterms);
    }
  }
  else {
    // create it
    $db_res         = db_fetch_array(db_query($cvterm_check, 'feature_type', $type));
    $type_id        = $db_res['cvterm_id'];
    $feature_insert = "INSERT INTO {gmod_dbsf_feature} (uniquename,dbxref_id,type_id,organism_id,passkey) ". " VALUES ('$uname',$dbxref_id,$type_id,$organism_id,'$passkey')";
    db_query($feature_insert);
    $feature_check = "SELECT feature_id from {gmod_dbsf_feature} where uniquename='$uname'";
    $db_res        = db_fetch_array(db_query($feature_check));
    $feature_id    = $db_res['feature_id'];
    if (empty($feature_id)) {
      drupal_set_message(t('I was not able to store the feature data at this time.'), 'error', FALSE);
      return FALSE;
    }
    // now add the seq/cvterms
    if (!empty($seq)) {
      db_query($feature_update, 'residues', $seq, $feature_id);
    }
    if ($type == 'gene_target') {
      genes4all_experiment_feature_add_gene($feature_id, $data, $cvterms);
    }
    elseif ($type == 'rnai_construct') {
      genes4all_experiment_feature_add_construct($feature_id, $data, $cvterms);
    }
  }
  return $feature_id;
}

function genes4all_experiment_featureprop_edit_text($feature_id, $cvname, $cvterm_name, $value_new = NULL) {
  if (empty($feature_id) || empty($cvname) || empty($cvterm_name)) {
    return FALSE;
  }

  $featureprop_check  = "SELECT value from {gmod_dbsf_featureprop} as featureprop WHERE feature_id=%d " . "AND type_id=(SELECT cvterm_id FROM {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " WHERE cvterm.name='%s' AND cv.name='%s')";
  $featureprop_update = "UPDATE {gmod_dbsf_featureprop} SET value='%s' WHERE feature_id=%d and type_id=" . "(SELECT cvterm_id FROM {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " WHERE cvterm.name='%s' AND cv.name='%s')";
  $featureprop_insert = "INSERT INTO {gmod_dbsf_featureprop} (value,feature_id,type_id) VALUES ('%s',%d," . "(SELECT cvterm_id FROM {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " WHERE cvterm.name='%s' AND cv.name='%s'))";
  $featureprop_delete = "DELETE FROM {gmod_dbsf_featureprop} WHERE feature_id=%d and type_id IN" . "(SELECT cvterm_id FROM {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " WHERE cv.name='%s')";

  $db_res = db_fetch_array(db_query($featureprop_check, $feature_id, $cvterm_name, $cvname));
  $check = $db_res['value'];
  if ($check == $value_new) {
    return FALSE;
  }
  elseif (!empty($check) && $check !== $value_new) {
    db_query($featureprop_delete, $feature_id, $cvname);
    db_query($featureprop_insert, $value_new, $feature_id, $cvterm_name, $cvname);
  }
  elseif (empty($check) && !empty($value_new)) {
    db_query($featureprop_insert, $value_new, $feature_id, $cvterm_name, $cvname);
  }
}

function genes4all_experiment_featurecvterm_edit_radio($feature_id, $cvname, $cvterm_name_existing = NULL, $cvterm_name_new = NULL) {
  if (empty($feature_id) || empty($cvname)) {
    return FALSE;
  }
  $feature_cvterm_insert = "INSERT INTO {gmod_dbsf_feature_cvterm} (feature_id,cvterm_id) VALUES (%d," . "(SELECT cvterm_id FROM {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " WHERE cvterm.name='%s' AND cv.name='%s'))";
  $feature_cvterm_update = "UPDATE {gmod_dbsf_feature_cvterm} set cvterm_id=%d WHERE feature_id=%d";
  $feature_cvterm_delete = "DELETE from {gmod_dbsf_feature_cvterm} WHERE feature_id=%d " . "AND cvterm_id IN (SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN " . " {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id where cv.name='%s')";
  //
  //dpm($cvterm_name_existing);
  if (!empty($cvterm_name_existing)) {
    if (!in_array($cvterm_name_new, $cvterm_name_existing)) {
      db_query($feature_cvterm_delete, $feature_id, $cvname);
      db_query($feature_cvterm_insert, $feature_id, $cvterm_name_new, $cvname);
    }
  }
  else {
    db_query($feature_cvterm_insert, $feature_id, $cvterm_name_new, $cvname);
  }
}

/**
 *
 * @param $feature_id feature_id
 * @param $data data from user
 * @param $cvterms data from database
 */
function genes4all_experiment_feature_add_gene($feature_id, $data, $cvterms = NULL) {
  if (empty($feature_id) || empty($data)) {
    return FALSE;
  }

  //aln region
  $array = gmod_dbsf_checkboxes_results2array($data['aln_region']);
  genes4all_experiment_cvterm_edit_checkboxes($array, $cvterms['aln_region'], 'aln_region', 'feature', $feature_id);
}

function genes4all_experiment_cvterm_edit_checkboxes($array, $cvterm_array = NULL, $cv_name, $type, $type_id) {
  if (empty($cv_name) || empty($type) || empty($type_id)) {
    // error was with array-> it can be null when nothing has been selected?
    //dpm(array('array'=>$array,'cv_name'=>$cv_name,'type'=>$type,'type_id'=>$type_id,'cvterm_array'=>$cvterm_array));
    drupal_set_message(t('Something went wrong with genes4all_experiment_cvterm_edit_checkboxes. Please let the administrator know via the feedback button.'), 'error', FALSE);
    return FALSE;
  }
  $feature_cvterm_insert = "INSERT INTO {gmod_dbsf_feature_cvterm} (feature_id,cvterm_id) VALUES (%d," . "(SELECT cvterm_id FROM {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " WHERE cvterm.name='%s' AND cv.name='%s'))";
  $feature_cvterm_delete = "DELETE from {gmod_dbsf_feature_cvterm} WHERE feature_id=%d " . "AND cvterm_id=(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN " . " {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id where cvterm.name='%s' AND cv.name='%s')";
  $resource_cvterm_insert = "INSERT INTO {gmod_dbsf_resource_cvterm} (resource_id,cvterm_id) VALUES (%d," . "(SELECT cvterm_id FROM {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " WHERE cvterm.name='%s' AND cv.name='%s'))";
  $resource_cvterm_delete = "DELETE from {gmod_dbsf_resource_cvterm} WHERE resource_id=%d " . "AND cvterm_id=(SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN " . " {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id where cvterm.name='%s' AND cv.name='%s')";

  $sql_insert = '';
  $sql_delete = '';
  if ($type == 'resource') {
    $sql_insert = $resource_cvterm_insert;
    $sql_delete = $resource_cvterm_delete;
  }
  elseif ($type == 'feature') {
    $sql_insert = $feature_cvterm_insert;
    $sql_delete = $feature_cvterm_delete;
  }
  else {
    drupal_set_message(t('Illegal type requested (%type).', array('%type' => $type)), 'error', FALSE);
    return FALSE;
  }

  // if array is empty, then we delete all data from there.
  if (empty($array)) {
    if (!empty($cvterm_array)) {
      foreach ($cvterm_array as $element) {
        db_query($sql_delete, $type_id, $element, $cv_name);
      }
    }
  }
  else {
    foreach ($array as $element) {
      $element = check_plain(trim($element));
      if (empty($cvterm_array) || !in_array($element, $cvterm_array)) {
        db_query($sql_insert, $type_id, $element, $cv_name);
      }
    }
    if (!empty($cvterm_array)) {
      foreach ($cvterm_array as $element) {
        if (!in_array($element, $array)) {
          db_query($sql_delete, $type_id, $element, $cv_name);
        }
      }
    }
  }
}

/**
 *
 * @param $feature_id feature_id
 * @param $data data from user
 * @param $cvterms data from database (if any)
 */
function genes4all_experiment_feature_add_construct($feature_id, $data, $cvterms = NULL) {
  $rna_probe_type = $data['rna_probe_type'];
  $purification = $data['purification'];
  $annealing_method = $data['annealing_method'];
  $construct_protocol_kit = $data['construct_protocol_kit'];
  $construct_protocol_details = $data['construct_protocol_details'];
  //$construct_protocol_invitro=$data['construct_protocol_invitro'];

  //if (!empty($construct_protocol_invitro)){$construct_protocol_invitro='TRUE';}
  // feature_cvterm: rna_probe_type, ['construct_protocol'][purification,annealing_method]
  // radios: check any cvterm in cv. if it the same, leave it. if it is not: delete and insert
  genes4all_experiment_featurecvterm_edit_radio($feature_id, 'rna_probe_type', $cvterms['rna_probe_type'], $rna_probe_type);
  genes4all_experiment_featurecvterm_edit_radio($feature_id, 'purification', $cvterms['purification'], $purification);
  genes4all_experiment_featurecvterm_edit_radio($feature_id, 'annealing_method', $cvterms['annealing_method'], $annealing_method);
  // featureprop: construct_protocol[construct_protocol_kit,$construct_protocol_invitro,construct_protocol_details]
  genes4all_experiment_featureprop_edit_text($feature_id, 'construct_protocol', 'construct_protocol_kit', $construct_protocol_kit);
  genes4all_experiment_featureprop_edit_text($feature_id, 'construct_protocol', 'construct_protocol_details', $construct_protocol_details);
  //genes4all_experiment_featureprop_edit_text($feature_id,'construct_protocol','construct_protocol_invitro',$construct_protocol_invitro);
}

function genes4all_experiment_resource_add_animal($resource_id, $data, $cvterms = NULL) {
  if (empty($resource_id) || empty($data)) {
    return FALSE;
  }
  $geolocation = $data['geolocation'];
  $dev_stage = $data['dev_stage'];
  $indiv_infection = $data['indiv_infection'];
  $pretreatment_protocol = $data['pretreatment_protocol'];
  $colony_infection = $data['colony_infection'];
  $origin = $data['origin'];


  //colony: colony infection[checkboxes]
  if (!empty($colony_infection)) {
    $array = gmod_dbsf_checkboxes_results2array($colony_infection);
    genes4all_experiment_cvterm_edit_checkboxes($array, $cvterms['colony_infection'], 'colony_infection', 'resource', $resource_id);
  }
  if (!empty($origin)) {
    genes4all_experiment_resourcecvterm_edit_radio($resource_id, 'origin', $cvterms['origin'], $origin);
  }
  if (!empty($indiv_infection)) {
    $array = gmod_dbsf_checkboxes_results2array($indiv_infection);
    genes4all_experiment_cvterm_edit_checkboxes($array, $cvterms['indiv_infection'], 'indiv_infection', 'resource', $resource_id);
  }
  if (!empty($geolocation)) {
    genes4all_experiment_resourceprop_edit_text($resource_id, 'animal_properties', 'geolocation', $geolocation);
  }
  if (!empty($pretreatment_protocol)) {
    genes4all_experiment_resourceprop_edit_text($resource_id, 'animal_properties', 'pretreatment_protocol', $pretreatment_protocol);
  }
  if (!empty($dev_stage)) {
    genes4all_experiment_resourcecvterm_edit_radio($resource_id, 'dev_stage', $cvterms['dev_stage'], $dev_stage);
  }
}

function genes4all_experiment_resource_add_assay_protocol($resource_id, $data, $cvterms = NULL) {
  if (empty($resource_id) || empty($data)) {
    return FALSE;
  }

  // Cs
  $detection_method = $data['detection_method'];
  // T
  $detect_time = $data['detect_time'];
  // Rs
  $dev_stage = $data['dev_stage'];
  // Rs
  $silencing_level = $data['silencing_level'];
  // T
  $silencing_accurate_level = $data['silencing_accurate_level'];
  // Rs
  $detect_dsrna_processing = $data['detect_dsrna_processing'];
  // T
  $assay_organ = $data['assay_organ'];

  if (!empty($silencing_level)) {
    genes4all_experiment_resourcecvterm_edit_radio($resource_id, 'silencing_level', $cvterms['silencing_level'], $silencing_level);
  }
  if (!empty($silencing_accurate_level)) {
    genes4all_experiment_resourceprop_edit_text($resource_id, 'result_set_properties', 'silencing_accurate_level', $silencing_accurate_level);
  }
  if (!empty($detect_dsrna_processing)) {
    genes4all_experiment_resourcecvterm_edit_radio($resource_id, 'detect_dsrna_processing', $cvterms['detect_dsrna_processing'], $detect_dsrna_processing);
  }
  if (!empty($detection_method)) {
    $array = gmod_dbsf_checkboxes_results2array($detection_method);
    genes4all_experiment_cvterm_edit_checkboxes($array, $cvterms['detection_method'], 'detection_method', 'resource', $resource_id);
  }
  if (!empty($detect_time)) {
    genes4all_experiment_resourceprop_edit_text($resource_id, 'assay_protocol_properties', 'detect_time', $detect_time);
  }
  if (!empty($dev_stage)) {
    genes4all_experiment_resourcecvterm_edit_radio($resource_id, 'dev_stage', $cvterms['dev_stage'], $dev_stage);
  }
  if (!empty($assay_organ)) {
    genes4all_experiment_resourceprop_edit_text($resource_id, 'assay_protocol_properties', 'assay_organ', $assay_organ);
  }
}

function genes4all_experiment_resource_add_delivery_protocol($resource_id, $data, $cvterms = NULL) {
  if (empty($resource_id) || empty($data)) {
    return FALSE;
  }
  // Cs
  $control_method = $data['control_method'];
  // T
  $delivery_detail = $data['delivery_detail'];
  // T
  $control_method_detail = $data['control_method_detail'];
  // T
  $construct_concentration = $data['construct_concentration'];
  // Cs
  $delivery_method = $data['delivery_method'];
  // T
  $delivery_buffer = $data['delivery_buffer'];
  // T
  $adjuvant_name = $data['adjuvant_name'];
  // T
  $adjuvant_amount = $data['adjuvant_amount'];
  // T
  $replicates = $data['replicates'];
  // T
  $delivery_organ = $data['delivery_organ'];

  if (!empty($control_method)) {
    $array = gmod_dbsf_checkboxes_results2array($control_method);
    genes4all_experiment_cvterm_edit_checkboxes($array, $cvterms['control_method'], 'control_method', 'resource', $resource_id);
  }
  if (!empty($delivery_detail)) {
    genes4all_experiment_resourceprop_edit_text($resource_id, 'delivery_protocol_properties', 'delivery_detail', $delivery_detail);
  }
  if (!empty($control_method_detail)) {
    genes4all_experiment_resourceprop_edit_text($resource_id, 'delivery_protocol_properties', 'control_method_detail', $control_method_detail);
  }
  if (!empty($delivery_method)) {
    $array = gmod_dbsf_checkboxes_results2array($delivery_method);
    genes4all_experiment_cvterm_edit_checkboxes($array, $cvterms['delivery_method'], 'delivery_method', 'resource', $resource_id);
  }
  if (!empty($construct_concentration)) {
    genes4all_experiment_resourceprop_edit_text($resource_id, 'delivery_protocol_properties', 'construct_concentration', $construct_concentration);
  }
  if (!empty($delivery_buffer)) {
    genes4all_experiment_resourceprop_edit_text($resource_id, 'delivery_protocol_properties', 'delivery_buffer', $delivery_buffer);
  }
  if (!empty($adjuvant_name)) {
    genes4all_experiment_resourceprop_edit_text($resource_id, 'delivery_protocol_properties', 'adjuvant_name', $adjuvant_name);
  }
  if (!empty($adjuvant_amount)) {
    genes4all_experiment_resourceprop_edit_text($resource_id, 'delivery_protocol_properties', 'adjuvant_amount', $adjuvant_amount);
  }
  if (!empty($replicates)) {
    genes4all_experiment_resourceprop_edit_text($resource_id, 'delivery_protocol_properties', 'replicates', $replicates);
  }
  if (!empty($delivery_organ)) {
    genes4all_experiment_resourceprop_edit_text($resource_id, 'delivery_protocol_properties', 'delivery_organ', $delivery_organ);
  }
}

function genes4all_experiment_resourceprop_edit_text($resource_id, $cvname, $cvterm_name, $value_new = NULL) {
  if (empty($resource_id) || empty($cvname) || empty($cvterm_name)) {
    return FALSE;
  }

  $resourceprop_check  = "SELECT value from {gmod_dbsf_resourceprop} as resourceprop WHERE resource_id=%d " . "AND type_id=(SELECT cvterm_id FROM {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " WHERE cvterm.name='%s' AND cv.name='%s')";
  $resourceprop_update = "UPDATE {gmod_dbsf_resourceprop} SET value='%s' WHERE resource_id=%d and type_id=" . "(SELECT cvterm_id FROM {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " WHERE cvterm.name='%s' AND cv.name='%s')";
  $resourceprop_insert = "INSERT INTO {gmod_dbsf_resourceprop} (value,resource_id,type_id) VALUES ('%s',%d," . "(SELECT cvterm_id FROM {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " WHERE cvterm.name='%s' AND cv.name='%s'))";
  $resourceprop_delete = "DELETE FROM {gmod_dbsf_resourceprop} WHERE resource_id=%d and type_id IN " . "(SELECT cvterm_id FROM {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " WHERE cv.name='%s')";
  $value_new           = check_plain(trim($value_new));
  // this can be avoided, as we already know the value when this function is called.
  $db_res = db_fetch_array(db_query($resourceprop_check, $resource_id, $cvterm_name, $cvname));
  $check = $db_res['value'];
  if ($check == $value_new) {
    return FALSE;
  }
  elseif (!empty($check) && $check !== $value_new) {
    db_query($resourceprop_delete, $resource_id, $cvname);
    db_query($resourceprop_insert, $value_new, $resource_id, $cvterm_name, $cvname);
  }
  elseif (empty($check) && !empty($value_new)) {
    db_query($resourceprop_insert, $value_new, $resource_id, $cvterm_name, $cvname);
  }
}

function genes4all_experiment_resourcecvterm_edit_radio($resource_id, $cvname, $cvterm_name_existing = NULL, $cvterm_name_new = NULL) {
  if (empty($resource_id) || empty($cvname)) {
    return FALSE;
  }
  $resource_cvterm_insert = "INSERT INTO {gmod_dbsf_resource_cvterm} (resource_id,cvterm_id) VALUES (%d," . "(SELECT cvterm_id FROM {gmod_dbsf_cvterm} as cvterm JOIN {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id " . " WHERE cvterm.name='%s' AND cv.name='%s'))";
  $resource_cvterm_update = "UPDATE {gmod_dbsf_resource_cvterm} set cvterm_id=%d WHERE resource_id=%d";
  $resource_cvterm_delete = "DELETE from {gmod_dbsf_resource_cvterm} WHERE resource_id=%d " . "AND cvterm_id IN (SELECT cvterm_id from {gmod_dbsf_cvterm} as cvterm JOIN " . " {gmod_dbsf_cv} as cv ON cv.cv_id=cvterm.cv_id where cv.name='%s')";
  $cvterm_name_new = check_plain(trim($cvterm_name_new));
  if (!empty($cvterm_name_existing)) {
    if (!in_array($cvterm_name_new, $cvterm_name_existing)) {
      db_query($resource_cvterm_delete, $resource_id, $cvname);
      db_query($resource_cvterm_insert, $resource_id, $cvterm_name_new, $cvname);
    }
  }
  else {
    db_query($resource_cvterm_insert, $resource_id, $cvterm_name_new, $cvname);
  }
}

function genes4all_experiment_get_add_taxonomy($data) {
  if (empty($data)) {
    return FALSE;
  }
  $tax_class  = $data['tax_class'];
  $tax_order  = $data['tax_order'];
  $tax_family = $data['tax_family'];
  $genus      = $data['genus'];
  $species    = $data['species'];
  if (empty($tax_order) || empty($tax_family) || empty($genus) || empty($species)) {
    return FALSE;
  }

  $update_org = "UPDATE {gmod_dbsf_organism} set ncbi_taxid=". $data['ncbi_taxid'] ." WHERE organism_id=%d";
  $select_org = "SELECT organism_id,ncbi_taxid from {gmod_dbsf_organism} where tax_class='$tax_class' AND tax_order='$tax_order'". " AND tax_family='$tax_family' AND genus='$genus' AND species='$species' ";

  $res        = db_fetch_array(db_query($select_org));
  $org_id     = $res['organism_id'];
  $ncbi_taxid = $res['ncbi_taxid'];
  if (!empty($org_id)) {
    //ncbi_taxid is optional so org may exist without and then someone may want to add it:
    // update ncbi_tax id if is not set or is different
    if (!empty($data['ncbi_taxid']) && $data['ncbi_taxid'] !== $ncbi_taxid) {
      db_query($update_org, $org_id);
    }
    return ($org_id);
  }
  else {
    //add it
    $insert_org = "INSERT INTO {gmod_dbsf_organism} (tax_class,tax_order,tax_family,genus,species)". " VALUES ('$tax_class','$tax_order','$tax_family','$genus','$species')";
    db_query($insert_org);
    $res        = db_fetch_array(db_query($select_org));
    $org_id     = $res['organism_id'];
    $ncbi_taxid = $res['ncbi_taxid'];
    if (empty($org_id)) {
      drupal_set_message(t('I was not able to store organism data.'), 'error', FALSE);
      return FALSE;
    }
    if (!empty($data['ncbi_taxid']) && $data['ncbi_taxid'] !== $ncbi_taxid) {
      db_query($update_org, $org_id);
    }
    return $org_id;
  }
}

function genes4all_experiment_authors_add($data, $pub_id) {
  if (empty($data) || empty($pub_id)) {
    return FALSE;
  }
  $success = 0;
  $total = 0;
  //$authors_address=$data['addresses'];
  $authors = array();
  foreach ($data as $key => $value) {
    if (preg_match('/author_first(\d+)/', $key, $match)) {
      $authors[$match[1]]['first'] = check_plain(trim($value));
    }
    elseif (preg_match('/author_last(\d+)/', $key, $match)) {
      $authors[$match[1]]['last'] = check_plain(trim($value));
    }
    elseif (preg_match('/author_email(\d+)/', $key, $match)) {
      $authors[$match[1]]['email'] = check_plain(trim($value));
    }
  }
  /*foreach ($authors_address as $key=>$value){
   if (preg_match('/address(\d+)/',$key,$match)){
   $authors[$match[1]]['address']=check_plain(trim($value));
   }
   }*/



  if (empty($authors)) {
    return FALSE;
  }

  $delete_pub_author = "DELETE from {gmod_dbsf_pub_author} where pub_id=$pub_id";
  db_query($delete_pub_author);

  foreach ($authors as $order => $details) {
    $total++;
    $first = "'". $details['first'] ."'";
    $last  = "'". $details['last'] ."'";
    $email = "'". $details['email'] ."'";
    // $address="'".$details['address']."'";

    $check_author = 'SELECT author_id from {gmod_dbsf_author} where'. " first_names=$first AND last_names=$last AND email=$email";
    $res_author   = db_fetch_array(db_query($check_author));
    $author_id    = $res_author['author_id'];
    if (empty($author_id)) {
      $insert_author = 'INSERT INTO {gmod_dbsf_author} (first_names,last_names,email)'. " VALUES ($first,$last,$email)";
      db_query($insert_author);
      $res_author = db_fetch_array(db_query($check_author));
      $author_id = $res_author['author_id'];
    }
    else {
      $update_author = "UPDATE {gmod_dbsf_author} SET %s='%s' where author_id=$author_id";
      $check1        = db_query($update_author, 'first_names', $details['first']);
      $check2        = db_query($update_author, 'last_names', $details['last']);
      $check3        = db_query($update_author, 'email', $details['email']);
      //$check4=db_query($update_author,'address',$details['address']);
      if ($check1 === FALSE || $check2 === FALSE || $check3 === FALSE) {
        drupal_set_message(t('Could not update author %first %last (%email).', array('%first' => $first, '%last' => $last, '%email' => $email)), 'error', TRUE);
      }
    }
    $insert_pub_author = 'INSERT INTO {gmod_dbsf_pub_author} (pub_id,author_id,rank) VALUES'. " ($pub_id,$author_id,$order)";
    db_query($insert_pub_author);
    $check_pub_author = "SELECT pub_id from {gmod_dbsf_pub_author} where author_id=$author_id AND pub_id=$pub_id";
    $res              = db_fetch_array(db_query($check_pub_author));
    $check            = $res['pub_id'];
    if (!empty($check)) {
      $success++;
    }
  }
  if (!empty($success === $total)) {
    return TRUE;
  }
  else {
    return FALSE;
  }
}

