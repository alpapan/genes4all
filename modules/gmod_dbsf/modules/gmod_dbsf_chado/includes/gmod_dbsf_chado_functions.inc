<?php

/**
 * Convert a Chado uniquename to an ID
 *
 * Used by for variety of Chado tables
 * @todo expand this function for all tables and also make it
 * work with Drupal too, then move it to gmod_dbsf.module
 *
 * @param $name
 * uniquename to fetch
 * @param $table_type
 * resource, feature, organism,pub or library
 * @param $reverse
 * assume that $name is an internal ID and fetch the uniquename
 * @param $reconnect
 * whether to reconnect to Drupal or leave the connection in Chado
 *
 * @return unknown_type
 * the internal ID or the uniquename if $reverse was used
 */
function gmod_dbsf_chado_name2id($name, $table_type, $reverse = FALSE, $reconnect = TRUE) {
  if (empty($name) || empty($table_type)) {
    return;
  }
  $key        = '';
  $table_pkey = '';
  $table_name = '';
  if ($table_type == 'resource') {
    $key = 'resource_id';
    $table_name = 'resource';
  }
  elseif ($table_type == 'feature') {
    $key = 'feature_id';
    $table_name = 'feature';
  }
  elseif ($table_type == 'organism') {
    $key = 'organism_id';
    $table_name = 'organism';
  }
  elseif ($table_type == 'library') {
    $key = 'library_id';
    $table_name = 'library';
  }
  elseif ($table_type == 'pub') {
    $key = 'pub_id';
    $table_name = 'pub';
  }
  else {
    drupal_set_message(t('Currently only resource, organism, library, pub and feature table types are allowed'), 'error', FALSE);
    return FALSE;
  }
  $sql_select = "SELECT $key as id FROM $table_name WHERE uniquename='$name'";
  if (!empty($reverse)) {
    $sql_select = "SELECT uniquename as id FROM $table_name WHERE $key='$name'";
  }
  $res = gmod_dbsf_query_execute('chado',$sql_select,NULL,$reconnect);
  return $res[0]['id'];
}

/**
 * Get NCBI taxonomy id for an organism
 *
 * @param $organism_id
 * The organism_id FROM organism table
 *
 */
function gmod_dbsf_chado_get_ncbi_taxid($organism_id) {
  if (empty($organism_id)) {
    return FALSE;
  }
  $sql_select = "SELECT ncbi_taxid FROM organism_ncbi_taxa WHERE organism_id=$organism_id";
  $res = db_fetch_array(db_query($sql_select));
  return $res['ncbi_taxid'];
}

/**
 * Get latest est2assembly accession ID
 *
 * @param $organism_id
 * The organism_id from organism table
 *
 * @return
 * Associative array with 'latest_assembly' (e.g. Ab) and
 * 'next_serial' (int) the next free serial
 */
function gmod_dbsf_chado_latest_assembly_serial($variables) {
  $organism_id = $variables['org'];
  $prefix      = $variables['prefix'];
  $db_id       = $variables['db'];
  // drupal has a temporary table, and therefore it will be faster to not query by type_id
  // also it will not have type contig (or indeed, it is not mandatory to have any other type)
  $drupal = $variables['drupal'];

  if (empty($organism_id) || !is_numeric($organism_id)) {
    return FALSE;
  }
  $sql_serial = "SELECT ic_accession_serial(uniquename) as id FROM feature WHERE ". " uniquename~'^$db_id' AND uniquename~'$prefix' AND organism_id=$organism_id " . " AND type_id = (SELECT cvterm_id FROM cvterm WHERE name='contig' AND cv_id=" . "(SELECT cv_id FROM cv WHERE name='sequence')) order by ic_accession_serial(uniquename) DESC limit 1";
  $sql_assembly = "SELECT ic_accession_assembly(uniquename) as id FROM feature WHERE ". " uniquename~'^$db_id' AND organism_id=$organism_id " . " AND type_id = (SELECT cvterm_id FROM cvterm WHERE name='contig' AND cv_id=" . "(SELECT cv_id FROM cv WHERE name='sequence')) order by ic_accession_assembly(uniquename) DESC limit 1";

  if (!empty($drupal)) {
    $sql_serial = "SELECT ic_accession_serial(uniquename) as id FROM feature WHERE ". " uniquename~'^$db_id' AND uniquename~'$prefix' AND organism_id=$organism_id " . " order by ic_accession_serial(uniquename) DESC limit 1";
    $sql_assembly = "SELECT ic_accession_assembly(uniquename) as id FROM feature WHERE ". " uniquename~'^$db_id' AND organism_id=$organism_id " . " order by ic_accession_assembly(uniquename) DESC limit 1";
  }
  $assembly_res = db_fetch_array(db_query($sql_assembly));
  $serial_res = db_fetch_array(db_query($sql_serial));
  if (!empty($assembly_res)) {
    return array('latest_assembly' => $assembly_res['id'], 'next_serial' => $serial_res['id'] + 1);
  }
  return FALSE;
}

/**
 * Add a library database cross-reference
 *
 * @param $lib_id
 * The internal library_id
 * @param $dbxref_id
 * The internal dbxref_id
 *
 * @return integer; the library_dbxref_id
 */
function gmod_dbsf_chado_add_library_dbxref($lib_id, $dbxref_id) {
  if (empty($lib_id) || empty($dbxref_id)) {
    return;
  }

  $check_sql  = "SELECT library_dbxref_id as id FROM library_dbxref WHERE dbxref_id=%d AND library_id=%d";
  $insert_sql = "INSERT INTO library_dbxref (library_id,dbxref_id) VALUES (%d,%d)";
  $check      = db_fetch_array(db_query($check_sql, $dbxref_id, $lib_id));
  if (empty($check['id'])) {
    db_query($insert_sql, $lib_id, $dbxref_id);
    $check = db_fetch_array(db_query($check_sql, $dbxref_id, $lib_id));
    if (empty($check['id'])) {
      drupal_set_message(t('Failed to insert date into your library. Perhaps a permissions error?'), 'error');
      return FALSE;
    }
  }
  return $check['id'];
}

/**
 * Add a library controlled vocabulary term
 *
 * @param $lib_id
 * The internal library_id
 * @param $cvterm_id
 * The CV term ID
 * @param $pub_id
 * The publication ID for storing this CVterm
 *
 * @return integer; the library_cvterm_id
 */
function gmod_dbsf_chado_add_library_cvterm($lib_id, $cvterm_id, $pub_id) {
  if (empty($lib_id) || empty($cvterm_id) || empty($pub_id)) {
    return;
  }

  $check_sql  = "SELECT library_cvterm_id as id FROM library_cvterm WHERE cvterm_id=%d AND library_id=%d";
  $insert_sql = "INSERT INTO library_cvterm (library_id,cvterm_id,pub_id) VALUES (%d,%d,%d)";
  $check      = db_fetch_array(db_query($check_sql, $cvterm_id, $lib_id));
  if (empty($check['id'])) {
    db_query($insert_sql, $lib_id, $cvterm_id, $pub_id);
    $check = db_fetch_array(db_query($check_sql, $cvterm_id, $lib_id));
    if (empty($check['id'])) {
      drupal_set_message(t('Failed to insert some data into your library. Perhaps a permissions error?'), 'error');
      return FALSE;
    }
  }
  return $check['id'];
}

/**
 * @param metadata : array('description'=>'','url_prefix'=>'','url'=>'')
 *
 * @see http://lsrn.org/lsrn/registry.html
 */
function gmod_dbsf_chado_load_dbxref_from_file($db,$db_metadata=NULL,$filename,$delimiter=",",$quotes='"',$noupdate=FALSE){
  if (empty($db)||empty($filename)){return;}
  if (!file_exists($filename)||!is_readable($filename)){
    gmod_dbsf_drupal_set_message('Cannot find or read file...','error');
    return FALSE;
  }
  if (empty($quotes)){$quotes='"';}
  if (empty($delimiter)){$delimiter=',';}
  ini_set('auto_detect_line_endings', TRUE);
  $dbxref_descr=array();
  $dbxref_names=array();
  $dbxref_properties=array();
  if (($handle = fopen($filename, 'r')) !== FALSE) {
    while (($data = fgetcsv($handle, 0, "$delimiter","$quotes")) !== FALSE) {
      if (empty($data[0])||empty($data[1])){continue;}
      $dbxref_names[]=check_plain($data[0]);
      $dbxref_descr[]=check_plain($data[1]);
      $slice = array_slice($data,2);
      if (!empty($slice)){
        $dbxref_properties[]=$slice;
      }
    }
    fclose($handle);
    ini_set('auto_detect_line_endings', FALSE);
  }else{
    gmod_dbsf_drupal_set_message('Cannot read file...','error');
    return FALSE;
  }
  if (!empty($dbxref_descr)&&!empty($dbxref_names)){
    $db_id=gmod_dbsf_add_db($db,$dbxref_names,$db_metadata,$dbxref_descr,$dbxref_properties,$noupdate);
    if (!empty($db_id)){
      return $db_id;
    }
  }
  return FALSE;
}

function gmod_dbsf_chado_load_ontology_from_file($cv,$cvdef=NULL,$filename,$delimiter=",",$quotes='"'){
  if (empty($cv)||empty($filename)){return;}
  if (!file_exists($filename)||!is_readable($filename)){
    gmod_dbsf_drupal_set_message('Cannot find or read file...','error');
    return FALSE;
  }
  if (empty($quotes)){$quotes='"';}
  if (empty($delimiter)){$delimiter=',';}
  ini_set('auto_detect_line_endings', TRUE);
  $cvterm_defs=array();
  $cvterm_names=array();
  $cvterm_properties=array();
  if (($handle = fopen($filename, 'r')) !== FALSE) {
    while (($data = fgetcsv($handle, 0, "$delimiter","$quotes")) !== FALSE) {
      //$data[0] is cvterm.name and $data[1] is definition
      if (empty($data[0])||empty($data[1])){continue;}
      $cvterm_names[]=check_plain($data[0]);
      $cvterm_defs[]=check_plain($data[1]);
      $slice = array_slice($data,2);
      if (!empty($slice)){
        $cvterm_properties[]=$slice;
      }
    }
    fclose($handle);
    ini_set('auto_detect_line_endings', FALSE);
  }else{
    gmod_dbsf_drupal_set_message('Cannot read file...','error');
    return FALSE;
  }
  if (!empty($cvterm_defs)&&!empty($cvterm_names)){
    $cv_id=gmod_dbsf_add_cv($cv,$cvterm_names,$cvdef,$cvterm_defs,$cvterm_properties);
    if (!empty($cv_id)){
      return $cv_id;
    }
  }
  return FALSE;
}

function gmod_dbsf_chado_count_library_feature_term($id,$term_name='EST',$parent=NULL){
  $library_ests = "SELECT count(*) from {library_feature} as lf JOIN {feature} as f ON lf.feature_id=f.feature_id WHERE f.type_id = %d AND library_id=%d";
  if (!empty($parent)){
    $parent = gmod_dbsf_get_feature_name_id($parent);
    $parent = $parent['id'];
    if (!empty($parent)){
      // -1 to link READ and EST... a hack.
      $library_ests .= " AND f.feature_id IN (select feature_id-1 from {featureloc} where srcfeature_id = $parent)";
    }
  }
  $feature_cvterm_id = gmod_dbsf_decypher_cvterm($term_name,'sequence');
  $feature_count = 0;
  if (!empty($est_cvterm_id)){
    $res = db_query($library_ests,$feature_cvterm_id,$id);
    $feature_count = $res['count'];
  }
  return $feature_count;
}

/**
 * Decide if an accession is a chado library id or a library name
 * TODO: decide how to support dbest (a prefix -> dbEST:\d+)?
 *
 * @param $id
 *
 * @return array
 */
function gmod_dbsf_chado_determine_library_accession($id) {
  $lib_id = '';
  $lib_name = '';
  // is cDNA BAC etc
  $lib_type = '';
  // decide if it is a lib name or a lib id.
  if (preg_match('/^\d+$/', $id)) {
    $lib_id = $id;
  }
  else {
    $lib_name = $id;
  }
  //TODO
  // get lib_id/lib_name
  if (!empty($lib_id)) {
    $lib_info_sql     = 'SELECT library_id as lib_id,library.name as name,cvterm.name as type_name,organism_id,uniquename as accession ' . ' FROM library join cvterm on type_id=cvterm_id'. " WHERE library_id= $lib_id";
    $lib_result       = db_query($lib_info_sql);
    $db_lib_data      = db_fetch_array($lib_result);
    $lib_name         = $db_lib_data['name'];
    $lib_type         = $db_lib_data['type_name'];
    $lib_org_id       = $db_lib_data['organism_id'];
    $lib_accession_id = $db_lib_data['accession'];
  }
  elseif (!empty($lib_name)) {
    $lib_info_sql     = 'SELECT library_id as lib_id,library.name as name,cvterm.name as type_name,organism_id,uniquename as accession ' . ' FROM library join cvterm on type_id=cvterm_id'. " WHERE library.name='$lib_name'";
    $lib_result       = db_query($lib_info_sql);
    $db_lib_data      = db_fetch_array($lib_result);
    $lib_id           = $db_lib_data['lib_id'];
    $lib_type         = $db_lib_data['type_name'];
    $lib_org_id       = $db_lib_data['organism_id'];
    $lib_accession_id = $db_lib_data['accession'];
  }
  //return empty if we failed to find it.
  if (empty($lib_id) || empty($lib_name)) {
    return;
  }
  else {
    return array('id' => $lib_id, 'name' => $lib_name, 'type' => $lib_type, 'organism' => $lib_org_id, 'accession' => $lib_accession_id);
  }
}

/**
 * Get properties of a library
 *
 * @param $lib_id
 * the library internal ID
 *
 * @return associative array
 * with the type name as key
 */
function gmod_dbsf_chado_getlibrary_properties($lib_id) {
  if (empty($lib_id)) {
    return;
  }
  $return_array = array();
  $sql_select   = "SELECT cvterm.name as type_name,lp.value FROM libraryprop as lp JOIN ". "cvterm as cvterm ON lp.type_id=cvterm.cvterm_id WHERE library_id=$lib_id";
  $db_res       = db_query($sql_select);
  while ($row = db_fetch_array($db_res)) {
    $return_array[$row['type_name']][] = $row['value'];
  }
  return $return_array;
}

/**
 * Get library database cross-references as a URL
 *
 * @param $lib_id
 * the library internal ID
 *
 * @return associative array
 * with the name and url as consecutive keys. the value is the
 * url prefix + the accession
 */
function gmod_dbsf_chado_getlibrary_dbxrefs($lib_id) {
  if (empty($lib_id)) {
    return;
  }
  $return_array = array();
  $sql_select   = "SELECT dbxref.accession,db.name,db.urlprefix FROM library_dbxref as ld JOIN dbxref as dbxref on dbxref.dbxref_id=ld.dbxref_id join db as db on db.db_id=dbxref.db_id  WHERE library_id=$lib_id";
  $db_res       = db_query($sql_select);
  while ($row = db_fetch_array($db_res)) {
    $return_array[$row['name']]['accessions'][] = $row['accession'];
    if (!empty($row['urlprefix'])) {
      $return_array[$row['name']]['url'][] = $row['urlprefix'] . $row['accession'];
    }
  }
  return $return_array;
}

/**
 *
 * @param $variables
 *
 * @return unknown_type
 */
function gmod_dbsf_chado_getlibraries_tags($variables) {
  $organism_id_array = !empty($variables['organism_id']) ? $variables['organism_id'] : '';
  $lib_id_array      = !empty($variables['lib_id']) ? $variables['lib_id'] : '';
  $tag_array         = array();
  // SQL result to harvest later
  $result    = '';
  $lib_str   = gmod_dbsf_array2csv($lib_id_array);
  $tag_array = array();
  // get all organisms;
  $org_data_sql = "SELECT genus,species,common_name FROM {organism} WHERE organism_id=%d";
  $lib_name_sql = "SELECT library_id as id,name as lib_name, organism_id from {library} WHERE 1=1 ";
  $organism_id_array = gmod_dbsf_array2csv($organism_id_array);
  if (!empty($lib_str)) {
    $lib_name_sql  .= " AND library_id IN ($lib_str) ";
  }if (empty($organism_id_array)) {
    // only libraries with features
    //$lib_name_sql .= ' AND organism_id  IN (SELECT organism_id FROM {type_feature_count_org_mat}) ';
    //$lib_name_sql .= ' AND library_id  IN (SELECT library_id FROM {library_feature}) ';
  }else{
    $lib_name_sql .= " AND organism_id IN ($organism_id_array)";
  }
  $lib_name_sql .=' ORDER BY name';
  $get_overview_tags = 'SELECT library_cvterm.cvterm_id as term_id,'
  . ' cv.name as tag,cvterm.name as value FROM library_cvterm '
  . ' join cvterm on library_cvterm.cvterm_id=cvterm.cvterm_id join cv on cvterm.cv_id=cv.cv_id '
  . ' WHERE library_cvterm.library_id=%d';
  $lib_res = db_query($lib_name_sql);
  while ($lib_data = db_fetch_array($lib_res)) {
    $access = gmod_dbsf_chadoapi_check_restriction('select', 'library', $lib_data['id'], FALSE);
    if ($access === FALSE) {
      continue;
    }
    $tag_array[$lib_data['lib_name']]['library_name'][$lib_data['id']]=$lib_data['lib_name'];
    $org_data = db_fetch_array(db_query($org_data_sql,$lib_data['organism_id']));
    $tag_array[$lib_data['lib_name']]['species'] = array(
      'genus' => $org_data['genus'],
      'species' => $org_data['species'],
    );
    if (!empty($org_data['common_name'])){
      $tag_array[$lib_data['lib_name']]['species']['common_name'] = $org_data['common_name'];
    }
    $cv_result = db_query($get_overview_tags,$lib_data['id']);
    while ($cv_data = db_fetch_array($cv_result)) {
      if (empty($tag_array[$lib_data['lib_name']][$cv_data['tag']]) && !empty($cv_data['term_id'])&& !empty($cv_data['value'])) {
        $tag_array[$lib_data['lib_name']][$cv_data['tag']][$cv_data['term_id']] = $cv_data['value'];
      }
    }
  }
  return $tag_array;
}

/**
 * Convert a feature to GFF
 *
 * This function is slower than the ones using materialized views
 * but  it does not depend on them.
 *
 * @param $feature_id
 * the feature internal ID or uniquename
 * @param $advanced
 * UNUSED. A placeholder for future ability to fetch data in a
 * more advanced and complete way
 *
 * @return associative array
 * with 'gff' and 'fsa' as keys. Each value is the relevant data stored as a string
 */
function gmod_dbsf_chado_feature2gff($feature_id, $advanced = FALSE) {
  if (empty($feature_id)) {
    return;
  }
  $feature_name = '';
  $feature_seq  = '';
  $gff_outdata  = '';
  $fsa_outdata  = '';
  $feature_type = '';
  $feature_org  = '';
  if (!is_numeric($feature_id)) {
    $feature_name  = $feature_id;
    $sql           = "SELECT f.feature_id,f.name as alias,f.residues,cvterm.name as term_name,organism.abbreviation as org_abbr " . " FROM feature f JOIN cvterm cvterm on cvterm_id=type_id JOIN organism organism ON organism.organism_id=f.organism_id" . " WHERE uniquename='%s'";
    $res           = db_fetch_array(db_query($sql, $feature_name));
    $feature_id    = $res['feature_id'];
    $feature_alias = $res['alias'];
    $feature_seq   = $res['residues'];
    $feature_type  = $res['term_name'];
    $feature_org   = $res['org_abbr'];
    if (empty($feature_id)) {
      return;
    }
  }
  else {
    $sql           = "SELECT uniquename,f.name as alias,residues,cvterm.name as term_name,organism.abbreviation as org_abbr " . " FROM feature f JOIN cvterm cvterm on cvterm_id=type_id JOIN organism organism ON organism.organism_id=f.organism_id" . " WHERE feature_id=%d";
    $res           = db_fetch_array(db_query($sql, $feature_id));
    $feature_name  = $res['uniquename'];
    $feature_alias = $res['alias'];
    $feature_seq   = $res['residues'];
    $feature_type  = $res['term_name'];
    $feature_org   = $res['org_abbr'];
  }
  if (empty($feature_name) || empty($feature_type)) {
    return;
  }
  $gff_feature_alias = '';
  if (!empty($feature_alias)) {
    // so that it is included in GFF string only if it exists.
    $gff_feature_alias = ';Alias='. $feature_alias;
    $feature_alias = ' '. $feature_alias;
  }

  if (!empty($advanced)) {
    /*
     *TODO
     $sql_all_select =
     "SELECT sf.name AS ref, dbx.accession AS source, cv.name AS type, fl.fmin + 1 AS fstart, fl.fmax AS fend, "
     ." af.significance AS score, fl.strand, fl.phase, f.seqlen, f.name, f.uniquename as uname, f.residues as seq, org.abbreviation as abbr "
     ." FROM feature f "
     ." LEFT JOIN featureloc fl ON f.feature_id = fl.feature_id "
     ." LEFT JOIN feature sf ON fl.srcfeature_id = sf.feature_id "
     ." LEFT JOIN feature_dbxref fd ON f.feature_id = fd.feature_id "
     ." LEFT JOIN dbxref dbx ON dbx.dbxref_id = fd.dbxref_id "
     ." AND (dbx.db_id IN ( SELECT db.db_id "
     ." FROM db "
     ." WHERE db.name::text = 'GFF_source'::text)) "
     ." LEFT JOIN cvterm cv ON f.type_id = cv.cvterm_id "
     ." LEFT JOIN analysisfeature af ON f.feature_id = af.feature_id "
     ." LEFT JOIN organism org ON f.organism_id=org.organism_id "
     ." WHERE sf.name='%s'";
     //dpm($sql_all_select);
     $res =db_query($sql_all_select,$feature_name);
     if (!empty($res)){
     if (!file_exists($outfile)){
     $gff_outdata .= "##gff-version 3\n";
     }  if (!file_exists($outfile.'.fsa')){
     $fsa_outdata .= "##gff-version 3\n##FASTA";
     }
     while ($row=db_fetch_array($res)){

     $ref=!empty($row['ref']) ? $row['ref']: $feature_name;
     $source = !empty($row['source']) ? $row['source']: 'source';
     // never empty
     $type = $row['type'];
     $seqlen=$row['seqlen'];
     $seq   = $row['seq'];
     $start = !empty($row['fstart']) ? $row['fstart']: 1;
     $end   = !empty($row['fend']) ? $row['fend']: $seqlen;
     if (empty($end) && !empty($seq)){
     //self-ref
     $end=strlen($seq);
     }
     $score  = !empty($row['score']) ? $row['score']: '.';
     $strand = !empty($row['strand']) ? $row['strand']: '.';
     $phase  = !empty($row['phase']) ? $row['phase']: '.';

     // never empty
     $uname =$row['uname'];
     $name = !empty($row['name']) ? $row['name']: $uname;
     $org_abbr = !empty($row['abbr']) ? $row['abbr']: 'unknown';
     $gff_outdata .="$ref\t$source\t$type\t$start\t$end\t$score\t$strand\t$phase\tID=$uname;Name=$name;organism=$org_abbr;score=$score";
     if (!empty($seq)){
     $seq=wordwrap($seq,80,"\n",TRUE);
     $fsa_outdata .=">$name\n$seq\n";
     }
     }
     unset($res);
     }
     */
  }
  else {
    $seqlength   = strlen($feature_seq);
    $gff_outdata = "$feature_name\tgmod_dbsf_curation\t$feature_type\t1\t$seqlength\t.\t+\t.\tID=$feature_name;Name=$feature_name". "$gff_feature_alias;organism=$feature_org\n";
    $fsa_outdata = ">". $feature_name . $feature_alias ."\n$feature_seq\n";
  }
  // the data will be appended to the GFF and FASTA outfiles.
  return array('gff' => $gff_outdata, 'fsa' => $fsa_outdata);
}

/**
 * Upload GFF3 annotation files to Chado
 *
 * This function is not new sequence data, the reference must already exist in Chado.
 * It prepares the commandline program ic_chado_loadcv.pl (part of est2assembly) so
 * the data must be compatible to that one (BLAST, EC, KEGG, InterProScan)
 * This function is used as part of the BATCH API
 * @see batch_set
 * @see genes4all_curate_load_approved in genes4all_curate_feature.inc
 *
 *
 * @param $dsn
 * The database DSN to connect to
 * @param $gff_files
 * One or more GFF files as 2D-associative array with the first key being the ID
 * from the batch relevant functions.
 * @param $context
 * The batch_set $context variable is automatically parsed when this function is
 * called via a batch-implementing function
 *
 * @return
 * updates the $context['results'] variable which is eventually parsed into $results
 * in the batch API. The associative has the 'chado']['par' and 'chado']['out elements
 * set to the cmdline arguments and output file respectively
 */
function gmod_dbsf_chado_gff2chado($dsn, $gff_files, &$context) {
  $execdir = './'. drupal_get_path('module', 'gmod_dbsf') .'/scripts';
  $chado_load_exec = $execdir .'/ic_chado_loadcv.pl';
  if (!file_exists($chado_load_exec)) {
    drupal_set_message(t('Failed to find executable %chado_load_exec.', array('%chado_load_exec' => $chado_load_exec)), 'error');
    return FALSE;
  }
  elseif (!is_executable($chado_load_exec)) {
    drupal_set_message(t('Cannot execute the %f program.', array('%f' => $chado_load_exec)), 'error');
    return FALSE;
  }
  if (empty($gff_files) || empty($dsn)) {
    return;
  }
  $gff_files_str = '';
  foreach ($gff_files as $id) {
    foreach ($id as $file) {
      $gff_files_str .= ' '. implode(' ', $file);
    }
  }

  $chado_load_exec = 'nice '. $chado_load_exec;
  $output          = array();
  $tmpdir          = file_directory_temp();
  $par_chado       = '';
  $par_chado       = " -dsn $dsn -gff $gff_files_str 2>&1";
  // load annotation data only to chado: first non-blast
  exec($chado_load_exec . $par_chado, $output);
  $context['results']['chado']['par'] = $par_chado;
  $context['results']['chado']['out'] = $output;
}

/**
 * Import a GFF3 to a seqfeature database
 *
 * @see gmod_dbsf_chado_gff2chado
 *
 * @param $db_array
 * @param $gff_files
 * @param $context
 *
 * @return unknown_type
 */
function gmod_dbsf_chado_gff2seqfeature($db_array, $gff_files, &$context) {
  $execdir = './'. drupal_get_path('module', 'gmod_dbsf') .'/scripts';
  $seqfeat_load_exec = $execdir .'/bp_seqfeature_load.pl';
  if (!file_exists($seqfeat_load_exec)) {
    drupal_set_message(t('Failed to find executable %seqfeat_load_exec.', array('%seqfeat_load_exec' => $seqfeat_load_exec)), 'error');
    return;
  }
  elseif (!is_executable($seqfeat_load_exec)) {
    drupal_set_message(t('Cannot execute the %f program.', array('%f' => $seqfeat_load_exec)), 'error');
    return FALSE;
  }
  if (empty($gff_files) || empty($db_array)) {
    return;
  }
  $gff_files_str = '';
  if (is_array($gff_files)) {
    $gff_files_str = implode(' ', $gff_files);
  }
  $seqfeat_load_exec = 'nice '. $seqfeat_load_exec;
  $output            = array();
  $tmpdir            = file_directory_temp();
  foreach ($db_array['dbnames'] as $dbname) {
    $par = '';
    $test_connection = pg_connect("dbname=". $dbname ." host=". $db_array['host'] ." port=". $db_array['port'] ." user=". $db_array['user'] .' password='. $db_array['password']);
    if ($test_connection === FALSE) {
      drupal_set_message(t('Database %dbname either does not exists or cannot connect to it.', array('%dbname' => $dbname)), 'error');
    }
    else {
      pg_close($test_connection);
      $dsn = "dbname=". $dbname .";host=". $db_array['host'] .";port=". $db_array['port'];
      $par = " -a DBI::Pg -T $tmpdir -d '$dsn' " . ' -u '. $db_array['user'] ." -p '". $db_array['password'] ."' -f $gff_files_str  2>&1";
      //load all data to seqfeature
      exec($seqfeat_load_exec . $par, $output);
      //dpm($seqfeat_load_exec.$par);
    }
    $context['results']['seqfeature']['par'][] = $par;
  }
  //for BATCH API
  $context['results']['seqfeature']['out'] = $output;
}

/**
 * Get assembly in a format suitable for GFF
 *
 * @see
 *
 * @param $features (unique)names of features.
 *
 * @return unknown_type
 */
function gmod_dbsf_chado_get_assembly($variables) {
  if (empty($variables)) {
    return;
  }
  if (!is_array($variables)) {
    $variables = array($variables);
  }

  $features        = !empty($variables['features']) ? $variables['features'] : $variables;
  if (empty($features)) {
    return ('No features requested');
  }
  if (!empty($variables['outfile'])) {
    ini_set('auto_detect_line_endings', TRUE);
    if (!$outhandle = fopen($variables['outfile'], 'at')) {
      drupal_set_message(t('Could not create %outfile.', array('%outfile' => $variables['outfile'])), 'error');
      return FALSE;
    }
    fwrite($outhandle, "##gff-version 3\n");
  }
  $assembly_array  = array();
  $fasta_array     = array();
  $sql_fasta       = "SELECT residues FROM feature WHERE uniquename = '%s'";
  $sql_prot2contig = 'SELECT name FROM feature WHERE type_id = ' . "(SELECT cvterm_id FROM cvterm WHERE name = 'contig' AND cv_id=(SELECT cv_id FROM cv WHERE name='sequence'))" . ' AND feature_id IN (SELECT srcfeature_id FROM featureloc WHERE feature_id = ' . " (SELECT feature_id FROM feature WHERE uniquename = '%s'))";

  $sql_assembly = 'SELECT ref,source,type,fstart,fend,score,strand,phase,name,seqlen FROM gff3view_mat ' . " WHERE ref = '%s' OR name = '%s' ORDER by type";
  if (!empty($variables['no_constituents'])) {
    $sql_assembly = 'SELECT ref,source,type,fstart,fend,score,strand,phase,name,seqlen FROM gff3view_mat ' . " WHERE name = '%s'";
  }
  foreach ($features as $feature) {
    $feature_query = $feature;
    // check to see if it is a polypeptide, if yes, get contig name.
    // don't forget changing to Aorf adding _0 and closing ' and 2 brackets
    if (strstr($feature, 'Apep') || strstr($feature, 'Aorf')) {
      $orf_id = $feature;
      if (strstr($orf_id, 'Apep')) {
        $orf_id = str_replace('Apep', 'Aorf', $orf_id);
      }
      $orf_id .= '_0';
      if (!empty($orf_id)) {
        $result        = db_query($sql_prot2contig, $orf_id);
        $res           = db_fetch_array($result);
        $feature_query = $res['name'];
      }
    }
    // get reads and contigs. if a row has no ref, then it is a contig (also type==contig).
    // then foreach name, $name=str_replace('_assembly','',$name);
    // then foreach name get sequence as above.
    $res = db_query($sql_assembly, $feature_query, $feature_query);
    while ($dbdata = db_fetch_array($res)) {
      $name_strip = str_replace('_assembly', '', $dbdata['name']);
      $name       = $dbdata['name'];
      $strand     = "+";
      if ($dbdata['strand'] == -1 || $dbdata['strand'] == '-') {
        $strand = '-';
      }
      $assembly_array[$name] = array(
        'ref' => $dbdata['ref'],
        'source' => $dbdata['source'],
        'type' => $dbdata['type'],
        'fstart' => $dbdata['fstart'],
        'fend' => $dbdata['fend'],
        'score' => $dbdata['score'],
        'strand' => $strand,
        'phase' => $dbdata['phase'],
        'tags' => array(),
      );

      if ($dbdata['type'] == 'contig') {
        $assembly_array[$name]['ref'] = $name;
        $assembly_array[$name]['fstart'] = '1';
        $assembly_array[$name]['fend'] = $dbdata['seqlen'];
        $assembly_array[$name]['tags']['ID'] = $name;
      }
      elseif ($dbdata['type'] == 'read') {
        $assembly_array[$name]['tags']['ID'] = $name_strip;
        $sql_read_target = 'SELECT fstart,fend FROM gff3view_mat WHERE ref='. "(SELECT regexp_replace ('$name' , '_assembly','')) AND name='$name'";
        $result_read = db_query($sql_read_target);
        $dbdata_read = db_fetch_array($result_read);
        $assembly_array[$name]['tags']['Target'] = $name_strip .' '. $dbdata_read['fstart'] .' '. $dbdata_read['fend'];
        $assembly_array[$name]['tags']['Parent'] = $dbdata['ref'];
        $assembly_array[$name]['tags']['#target']['id'] = $name_strip;
        $assembly_array[$name]['tags']['#target']['value'] = ' '. $dbdata_read['fstart'] .' '. $dbdata_read['fend'];
      }
      unset($dbdata);
      if (empty($fasta_array[$name_strip])) {
        $dbdata_fasta = db_fetch_array(db_query($sql_fasta, $name_strip));
        $fasta_array[$name_strip] = str_replace('X', 'N', $dbdata_fasta['residues']);
      }
      if (!empty($variables['outfile'])) {
        foreach ($assembly_array as $name => $value_array) {
          $data = '';
          foreach ($value_array as $column => $value) {
            if (empty($value)) {
              $value = '.';
            }
            if (is_array($value)) {
              foreach ($value as $key => $r) {
                // if it starts with a hash, ignore it.
                if (strpos($key, '#') === FALSE) {
                  $data .= $key .'='. $r .';';
                }
              }
            }
            else {
              $data .= $value ."\t";
            }
          }
          $data = trim($data);
          fwrite($outhandle, $data ."\n");
        }
        $assembly_array = array();
      }
    }
  }
  //end each feature
  if (!empty($variables['outfile'])) {
    fwrite($outhandle, "##FASTA\n");
    foreach ($fasta_array as $feature => $residues) {
      fwrite($outhandle, ">$feature\n$residues\n");
    }
    fclose($outhandle);
    return $variables['outfile'];
  }
  else {
    return (array($assembly_array, $fasta_array));
  }
}

/**
 * Get an assembly as GFF
 *
 * helper function used by exploration module
 * @see gmod_dbsf_chado_get_assembly
 * @see gmod_dbsf_download_bulk_get
 *
 * @param string or array of reference name(s)
 *
 * @return HTML
 */
function gmod_dbsf_chado_assembly2gff($ref_names) {

  if (empty($ref_names)) {
    return;
  }
  $assembly       = '';
  $return_array   = gmod_dbsf_chado_get_assembly($ref_names);
  $assembly_array = $return_array[0];
  $fasta_array    = $return_array[1];
  if (!empty($assembly_array)) {
    $assembly = '<pre><font size=-3><br>##gff-version3<br>';
    foreach ($assembly_array as $id => $value_array) {
      foreach ($value_array as $key => $value) {
        if (is_array($value)) {
          //Then key is probably tags
          foreach ($value as $tag => $data) {
            if (!is_array($data)) {
              $assembly .= $tag .'='. $data .';';
            }
          }
        }
        elseif (empty($value) && ($key !== 'residues' || $key !== 'target')) {
          $assembly .= "\t" . '.' . "\t";
        }
        elseif ($key !== 'residues' && $key !== 'target') {
          $assembly .= $value ."\t";
        }
      }
      $assembly = rtrim($assembly, ';');
      $assembly = trim($assembly);
      $assembly .= '<br>';
      if (!empty($assembly_array[$id]['tags']['#target'])) {
        $assembly .= $assembly_array[$id]['tags']['#target']['id'] ."\t" . 'basecaller' . "\t" . 'EST' . "\t" . '1' . "\t". strlen($fasta_array[$assembly_array[$id]['tags']['#target']['id']]) ."\t" . '.' . "\t" . '+' . "\t" . '.' . "\t" . 'ID='. $assembly_array[$id]['tags']['#target']['id'] .'<br>';
      }
    }
    // now FASTA
    $fasta_output = '';
    if (!empty($fasta_array)) {
      $assembly .= '##FASTA<br>';
      foreach ($fasta_array as $id => $sequence) {
        $sequence_wrap = wordwrap($sequence, 80, '<br/>', TRUE);
        $fasta_output .= "<p>>$id<br>$sequence_wrap</p>";
      }
    }
    $assembly .= $fasta_output .'</font></pre>';
  }
  return $assembly;
}

/**
 *
 * @param $features
 *
 * @return unknown_type
 */
function gmod_dbsf_chado_get_fasta($variables) {
  if (empty($variables)) {
    return;
  }
  if (!is_array($variables)) {
    $variables = array($variables);
  }
  $features = !empty($variables['features']) ? $variables['features'] : $variables;
  if (empty($features)) {
    return ('No features requested');
  }
  if (!empty($variables['outfile'])) {
    ini_set('auto_detect_line_endings', TRUE);
    if (!$outhandle = fopen($variables['outfile'], 'at')) {
      drupal_set_message(t('Could not create %outfile.', array('%outfile' => $variables['outfile'])), 'error');
      return FALSE;
    }
  }
  $fasta_array = array();
  $feature_query = '';
  if (count($features) > 100) {
    $sql = "SELECT name,residues FROM feature WHERE uniquename ='%s'";
    foreach ($features as $feature) {
      $result = db_query($sql, $feature);
      while ($dbdata = db_fetch_array($result)) {
        if (empty($dbdata['residues'])) {
          continue;
        }
        $fasta_array[$dbdata['name']] = $dbdata['residues'];
        if (!empty($variables['outfile'])) {
          fwrite($outhandle, ">". $dbdata['name'] ."\n". wordwrap($dbdata['residues'], 80, "\n", TRUE) ."\n");
          $fasta_array = array();
        }
      }
    }
  }
  else {
    foreach ($features as $feature) {
      $feature_query .= "'$feature',";
    }
    $feature_query = rtrim($feature_query, ",");
    $sql           = 'SELECT name,residues FROM feature WHERE uniquename IN ('. $feature_query .')';
    $result        = db_query($sql);
    while ($dbdata = db_fetch_array($result)) {
      if (empty($dbdata['residues'])) {
        continue;
      }
      $fasta_array[$dbdata['name']] = $dbdata['residues'];
      if (!empty($variables['outfile'])) {
        fwrite($outhandle, ">". $dbdata['name'] ."\n". wordwrap($dbdata['residues'], 80, "\n", TRUE) ."\n");
        $fasta_array = array();
      }
    }
  }
  if (!empty($variables['outfile'])) {
    fclose($outhandle);
    return $variables['outfile'];
  }
  return ($fasta_array);
}


/**
 * Get BLAST databases that were used for searching
 *
 * @TODO utilize gmod_dbsf_get_feature_cvterm_and_prop()
 *
 * @param $feature_id
 * @param $limit
 *
 * @return unknown_type
 */
function gmod_dbsf_chado_get_feature_blastdbs($feature_id, $limit = NULL) {
  if (empty($feature_id)) {
    return;
  }
  $blastdb_pro_against_array = array();
  $blastdb_nuc_against_array = array();
  //Chado-specific SQLs so no {} fields
  /*
  $blastdb_against_sql = 'SELECT value from featureprop where type_id ='
  . " (select cvterm_id from cvterm where name ='searched_against') "
  . " and feature_id= $feature_id ORDER by lower(value)";
  */
  $blastdb_pro_against_sql =
   "SELECT DISTINCT ON (lower(db.name)) db.name,db.description from db JOIN dbxref ON db.db_id=dbxref.db_id where db.db_id != (SELECT db_id from db where name='GFF_source')"
   ." AND dbxref_id IN (SELECT  dbxref_id from feature_dbxref where feature_id IN (select feature_id from feature where type_id IN (SELECT cvterm_id from cvterm "
   . "where name='protein_match' and cv_id=(SELECT cv_id from cv where name='sequence'))"
   ." AND feature_id IN (SELECT feature_id from featureloc where srcfeature_id=$feature_id))) ORDER by lower(db.name)";
   $blastdb_nuc_against_sql =
   "SELECT DISTINCT ON (lower(db.name)) db.name,db.description from db JOIN dbxref ON db.db_id=dbxref.db_id where db.db_id != (SELECT db_id from db where name='GFF_source')"
   ." AND dbxref_id IN (SELECT  dbxref_id from feature_dbxref where feature_id IN (select feature_id from feature where type_id IN (SELECT cvterm_id from cvterm "
   ." where name='nucleotide_match' and cv_id=(SELECT cv_id from cv where name='sequence'))"
   ." AND feature_id IN (SELECT feature_id from featureloc where srcfeature_id=$feature_id))) ORDER by lower(db.name)";

   if (!empty($limit)) {
     $blastdb_pro_against_sql .= " LIMIT $limit ";
     $blastdb_nuc_against_sql .= " LIMIT $limit ";
   }
   //Protein matches
   $result = db_query($blastdb_pro_against_sql);
   while ($db_blastdb = db_fetch_array($result)) {
     $blastdb_pro_against_array[] = array($db_blastdb['name'],$db_blastdb['description']);
   }
   unset($result);
   //nucleotide matches
   $result = db_query($blastdb_nuc_against_sql);
   while ($db_blastdb = db_fetch_array($result)) {
     $blastdb_nuc_against_array[] = array($db_blastdb['name'],$db_blastdb['description']);
   }

   return (array($blastdb_pro_against_array,$blastdb_nuc_against_array));
}

/**
 * Get BLAST databases that were used for searching
 *
 * @TODO utilize gmod_dbsf_get_feature_cvterm_and_prop()
 *
 * @param $feature_id
 * @param $limit
 *
 * @return unknown_type
 */
function gmod_dbsf_chado_get_featuregff_blastdbs($feature_id) {
  if (empty($feature_id)) {
    return;
  }
  $blastdb_pro_against_array = array();
  $blastdb_nuc_against_array = array();
  $blastdb_pro_against_sql =
   "SELECT DISTINCT dbxref.accession from feature_dbxref fd JOIN dbxref ON dbxref.dbxref_id=fd.dbxref_id "
   ." WHERE dbxref.db_id = (SELECT db_id from db where name='GFF_source') "
   .' AND fd.feature_id IN (select feature_id from feature where type_id IN (SELECT cvterm_id from cvterm '
   ." WHERE name='protein_match' and cv_id=(SELECT cv_id from cv where name='sequence'))"
   ." AND feature_id IN (SELECT feature_id from featureloc where srcfeature_id=$feature_id)) ORDER by accession";
   $blastdb_nuc_against_sql =
   "SELECT DISTINCT dbxref.accession from feature_dbxref fd JOIN dbxref ON dbxref.dbxref_id=fd.dbxref_id "
   ." WHERE dbxref.db_id = (SELECT db_id from db where name='GFF_source') "
   .' AND fd.feature_id IN (select feature_id from feature where type_id IN (SELECT cvterm_id from cvterm '
   ." WHERE name='nucleotide_match' and cv_id=(SELECT cv_id from cv where name='sequence'))"
   ." AND feature_id IN (SELECT feature_id from featureloc where srcfeature_id=$feature_id)) ORDER by accession";
   //Protein matches
   $result = db_query($blastdb_pro_against_sql);
   while ($db_blastdb = db_fetch_array($result)) {
     $blastdb_pro_against_array[] = $db_blastdb['accession'];
   }
   unset($result);
   //nucleotide matches
   $result = db_query($blastdb_nuc_against_sql);
   while ($db_blastdb = db_fetch_array($result)) {
     $blastdb_nuc_against_array[] = $db_blastdb['accession'];
   }

   return (array($blastdb_pro_against_array,$blastdb_nuc_against_array));
}


/**
 * Get keywords that appeared on BLAST.
 *
 * @TODO utilize gmod_dbsf_get_feature_cvterm_and_prop()
 *
 * @param $feature_id
 * @param $limit
 *
 *
 * @return unknown_type
 */
function gmod_dbsf_chado_get_feature_blast_values($feature_id, $limit = NULL) {
  if (empty($feature_id)) {
    return;
  }
  $blastdb_values_array = array();
  //Chado-specific SQLs so no {} fields
  $blastdb_values_sql = 'SELECT value from featureprop where type_id ='
  . "(select cvterm_id from cvterm where name ='inferred from electronic annotation' "
  . " and cv_id=(select cv_id from cv where name='Evidence Codes')) "
  . ' and feature_id ='. $feature_id;
  if (!empty($limit)) {
    $blastdb_values_sql .= " LIMIT $limit";
  }
  $result = db_query($blastdb_values_sql);
  while ($db_blastdb = db_fetch_array($result)) {
    $blastdb_values_array[] = $db_blastdb['value'];
  }
  return ($blastdb_values_array);
}

