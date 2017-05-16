<?php
// $Id$

/**
 * @file
 * Autocomplete functions use in genes4all and submodules
 */

/**
 * Autocomplete features
 *
 * @see user_autocomplete
 *
 */
function genes4all_feature_autocomplete($string = '') {
  $matches = array();
  gmod_dbsf_db_execute('chado');
  if (!empty($string)) {
    $result = db_query_range(
      "SELECT feature.uniquename FROM feature feature join cvterm as cvterm ON " . " cvterm.cvterm_id=feature.type_id " . " WHERE LOWER(feature.uniquename) LIKE LOWER('%s%%') AND cvterm.name IN ('ORF','polypeptide','contig','gene','mRNA','CDS')", $string, 0, 10
    );
    while ($res = db_fetch_object($result)) {
      $matches[$res->uniquename] = check_plain($res->uniquename);
    }
    if (empty($matches)) {
      $result = db_query_range(
        "SELECT feature.uniquename FROM feature feature join cvterm as cvterm ON " . " cvterm.cvterm_id=feature.type_id " . " WHERE cvterm.name IN ('ORF','polypeptide','contig','gene','mRNA','CDS')", 0, 10
      );
      while ($res = db_fetch_object($result)) {
        $matches[$res->uniquename] = $res->uniquename;
      }
    }
  }
  gmod_dbsf_db_execute();
  drupal_json($matches);
}

/**
 * Autocomplete curation terms
 *
 * @see user_autocomplete
 *
 */
function genes4all_curate_autocomplete_free_term($string = NULL) {
  $matches = array();
  if (!empty($string)) {
    gmod_dbsf_db_execute('chado');
    //Chado specific SQL
    $result = db_query_range(
      "SELECT name FROM cvterm WHERE cv_id IN "
    ."(SELECT cv_id from cv where name IN ('free_curation_terms','sequence','biological_process','molecular_function','cellular_component','EC','KEGG_PATHWAY'))"
    ." AND LOWER(name) LIKE LOWER('%s%%')", $string, 0, 10
    );
    while ($res = db_fetch_object($result)) {
      $name = $res->name;
      $matches[$name] = $name;
    }
    if (empty($matches)) {
      //Chado specific SQL
      $result = db_query_range(
        "SELECT name FROM cvterm WHERE cv_id IN "
      ."(SELECT cv_id from cv where name IN ('free_curation_terms','sequence','biological_process','molecular_function','cellular_component','EC','KEGG_PATHWAY'))"
      ." order by name", $string, 0, 10
      );
      while ($res = db_fetch_object($result)) {
        $name = check_plain($res->name);
        $matches[$name] = $name;
      }
    }
    gmod_dbsf_db_execute();
  }
  drupal_json($matches);
}

/**
 * Autocomplete library
 *
 * @see user_autocomplete
 *
 */
function genes4all_curate_autocomplete_library_cvterm($string = NULL) {
  $matches = array();
  if (!empty($string)) {
    gmod_dbsf_db_execute('chado');
    //Chado specific SQL
    $result = db_query_range(
      "select name from cvterm where cvterm_id IN (select distinct cvterm_id from library_cvterm) AND name ~LOWER('%s') order by name", $string, 0, 20
    );
    while ($res = db_fetch_object($result)) {
      $name = $res->name;
      $matches[$name] = $name;
    }
    if (empty($matches)) {
      //Chado specific SQL
      $result = db_query_range(
        "select name from cvterm where cvterm_id IN (select distinct cvterm_id from library_cvterm) order by name", $string, 0, 10
      );
      while ($res = db_fetch_object($result)) {
        $name = $res->name;
        $matches[$name] = $name;
      }
    }
    gmod_dbsf_db_execute();
  }
  drupal_json($matches);
}

/**
 * Autocomplete library dbname
 *
 * @see user_autocomplete
 *
 */
function genes4all_curate_autocomplete_library_dbname($string = NULL) {
  $matches = array();
  if (!empty($string)) {
    gmod_dbsf_db_execute('chado');
    //Chado specific SQL
    $result = db_query_range(
      "SELECT name from db where db_id IN (select distinct db_id from library_dbxref as ld join dbxref as d on d.dbxref_id=ld.dbxref_id) AND name LIKE LOWER('%s%%') order by name", $string, 0, 10
    );
    while ($res = db_fetch_object($result)) {
      $name = $res->name;
      $matches[$name] = $name;
    }
    if (empty($matches)) {
      //Chado specific SQL
      $result = db_query_range(
        "SELECT name from db where db_id IN (select distinct db_id from library_dbxref  as ld join dbxref as d on d.dbxref_id=ld.dbxref_id) order by name", $string, 0, 10
      );
      while ($res = db_fetch_object($result)) {
        $name = check_plain($res->name);
        $matches[$name] = $name;
      }
    }
    gmod_dbsf_db_execute();
  }
  drupal_json($matches);
}

/**
 * Autocomplete genus of organism
 *
 * @see user_autocomplete
 *
 */
function genes4all_curate_autocomplete_genus($string = NULL) {
  $matches = array();
  if (!empty($string)) {
    gmod_dbsf_db_execute('chado');
    //Chado specific SQL
    $result = db_query_range(
      "SELECT genus as name FROM organism WHERE lower(genus) LIKE LOWER('%s%%') order by genus", $string, 0, 10
    );
    while ($res = db_fetch_object($result)) {
      $name = $res->name;
      $matches[$name] = $name;
    }
    if (empty($matches)) {
      //Chado specific SQL
      $result = db_query_range(
        "SELECT genus as name FROM organism ORDER BY genus", $string, 0, 10
      );
      while ($res = db_fetch_object($result)) {
        $name = check_plain($res->name);
        $matches[$name] = $name;
      }
    }
    gmod_dbsf_db_execute();
  }
  drupal_json($matches);
}

/**
 * Autocomplete species of organism
 *
 * @see user_autocomplete
 *
 */
function genes4all_curate_autocomplete_species($string = NULL) {
  $matches = array();
  if (!empty($string)) {
    gmod_dbsf_db_execute('chado');
    //Chado specific SQL
    $result = db_query_range(
      "SELECT species as name FROM organism WHERE lower(species) LIKE LOWER('%s%%') order by species", $string, 0, 10
    );
    while ($res = db_fetch_object($result)) {
      $name = $res->name;
      $matches[$name] = $name;
    }
    if (empty($matches)) {
      //Chado specific SQL
      $result = db_query_range(
        "SELECT species as name FROM organism ORDER BY species", $string, 0, 10
      );
      while ($res = db_fetch_object($result)) {
        $name = check_plain($res->name);
        $matches[$name] = $name;
      }
    }
    gmod_dbsf_db_execute();
  }
  drupal_json($matches);
}

/**
 * Autocomplete publication of experiment
 *
 * @see user_autocomplete
 *
 */
function genes4all_experiment_publication_autocomplete($string = '') {
  $passkey = check_plain($_SESSION['passkey']);
  $matches = array();
  if ($string) {
    $result = db_query_range("SELECT uniquename FROM {gmod_dbsf_pub} WHERE LOWER(uniquename) LIKE LOWER('%s%%') AND passkey='$passkey'", $string, 0, 10);
    while ($res = db_fetch_object($result)) {
      $matches[$res->uniquename] = $res->uniquename;
    }
    if (empty($matches)) {
      $result = db_query_range("SELECT uniquename FROM {gmod_dbsf_pub} order by uniquename", $string, 0, 10);
      while ($res = db_fetch_object($result)) {
        $matches[$res->uniquename] = $res->uniquename;
      }
    }
  }
  drupal_json($matches);
}

/**
 * Autocomplete dbname of experiment
 *
 * @see user_autocomplete
 *
 */
function genes4all_experiment_dbname_autocomplete($string = NULL) {
  $matches = array();
  if (!empty($string)) {
    $result = db_query_range("SELECT name FROM {gmod_dbsf_db} WHERE LOWER(name) LIKE LOWER('%s%%')", $string, 0, 10);
    while ($res = db_fetch_object($result)) {
      $matches[$res->name] = check_plain($res->name);
    }
    if (empty($matches)) {
      $result = db_query_range("SELECT name FROM {gmod_dbsf_db} order by name", $string, 0, 10);
      while ($res = db_fetch_object($result)) {
        $matches[$res->name] = check_plain($res->name);
      }
    }
  }
  drupal_json($matches);
}

/**
 * Autocomplete target database
 *
 * @see user_autocomplete
 *
 */
function genes4all_experiment_targetdb_autocomplete($string = '') {
  $passkey = check_plain($_SESSION['passkey']);
  $matches = array();
  if (!empty($string) && !empty($passkey)) {
    $result = db_query_range(
      "SELECT uniquename FROM {gmod_dbsf_feature} as experiment_feature join {gmod_dbsf_cvterm} as experiment_cvterm ON " . " experiment_cvterm.cvterm_id=experiment_feature.type_id ". "  WHERE LOWER(experiment_feature.uniquename) LIKE LOWER('%s%%') AND experiment_cvterm.name='gene_target' AND passkey='$passkey'", $string, 0, 10
    );
    while ($res = db_fetch_object($result)) {
      $matches[$res->uniquename] = check_plain($res->uniquename);
    }
    if (empty($matches)) {
      $result = db_query_range(
        "SELECT uniquename FROM {gmod_dbsf_feature} as experiment_feature join {gmod_dbsf_cvterm} as experiment_cvterm ON " . " experiment_cvterm.cvterm_id=experiment_feature.type_id ". "  WHERE experiment_cvterm.name='gene_target' AND passkey='$passkey' order by uniquename", $string, 0, 10
      );
      while ($res = db_fetch_object($result)) {
        $matches[$res->uniquename] = $res->uniquename;
      }
    }
  }
  drupal_json($matches);
}

/**
 * Autocomplete construct database
 *
 * @see user_autocomplete
 *
 */
function genes4all_experiment_constructdb_autocomplete($string = '') {
  $passkey = check_plain($_SESSION['passkey']);
  $matches = array();
  if (!empty($string) && !empty($passkey)) {
    $result = db_query_range(
      "SELECT uniquename FROM {gmod_dbsf_feature} as experiment_feature join {gmod_dbsf_cvterm} as experiment_cvterm ON " . " experiment_cvterm.cvterm_id=experiment_feature.type_id ". "  WHERE LOWER(experiment_feature.uniquename) LIKE LOWER('%s%%') AND experiment_cvterm.name='rnai_construct' AND passkey='$passkey'", $string, 0, 10
    );
    while ($res = db_fetch_object($result)) {
      $matches[$res->uniquename] = check_plain($res->uniquename);
    }
    if (empty($matches)) {
      $result = db_query_range(
        "SELECT uniquename FROM {gmod_dbsf_feature} as experiment_feature join {gmod_dbsf_cvterm} as experiment_cvterm ON " . " experiment_cvterm.cvterm_id=experiment_feature.type_id ". "  WHERE experiment_cvterm.name='rnai_construct' AND passkey='$passkey' order by uniquename", $string, 0, 10
      );
      while ($res = db_fetch_object($result)) {
        $matches[$res->uniquename] = $res->uniquename;
      }
    }
  }
  drupal_json($matches);
}

/**
 * Autocomplete resource database
 *
 * @see user_autocomplete
 *
 */
function genes4all_experiment_resourcedb_autocomplete($string = '') {
  $passkey = check_plain($_SESSION['passkey']);
  $matches = array();
  if (!empty($string) && !empty($passkey)) {
    $result = db_query_range(
      "SELECT uniquename FROM {gmod_dbsf_resource} ". " WHERE LOWER(uniquename) LIKE LOWER('%s%%') AND passkey='$passkey'", $string, 0, 10
    );
    while ($res = db_fetch_object($result)) {
      $matches[$res->uniquename] = check_plain($res->uniquename);
    }
    if (empty($matches)) {
      $result = db_query_range(
        "SELECT uniquename FROM {gmod_dbsf_resource} ". "  WHERE passkey='$passkey' order by uniquename", $string, 0, 10
      );
      while ($res = db_fetch_object($result)) {
        $matches[$res->uniquename] = $res->uniquename;
      }
    }
  }
  drupal_json($matches);
}

/**
 * Autocomplete study database
 *
 * @see user_autocomplete
 *
 */
function genes4all_experiment_studydb_autocomplete($string = '') {
  $matches = array();
  if (!empty($string)) {
    $result = db_query_range(
      "SELECT uniquename FROM {gmod_dbsf_study} " . "  WHERE LOWER(uniquename) LIKE LOWER('%s%%')", $string, 0, 10
    );
    while ($res = db_fetch_object($result)) {
      $matches[$res->uniquename] = check_plain($res->uniquename);
    }
    if (empty($matches)) {
      $result = db_query_range(
        "SELECT uniquename FROM {gmod_dbsf_study} " . "  order by uniquename", $string, 0, 10
      );
      while ($res = db_fetch_object($result)) {
        $matches[$res->uniquename] = $res->uniquename;
      }
    }
  }
  drupal_json($matches);
}

