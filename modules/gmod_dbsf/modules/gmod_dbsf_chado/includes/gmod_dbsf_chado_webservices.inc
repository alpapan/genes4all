<?php

/**
 * @file
 * Web services for gmod_dbsf_chado Geospatial Genetics
 *
 * This module implements simple webservices for Drupal
 * i) Chado, focusing on the natural diversity module.
 * ii) The Atlas of Living Australia, @see http://ala.org
 * iii) TRIN @see trin.org.au
 * iv) others?
 *
 * @author Nathan Dunn <nathan.dunn@UNSWalumni.com>
 * @author Alexie Papanicolaou
 * @see complaints mailto:alexie@butterflybase.org
 *
 * @TODO ND: server side import of external DB/LSID data - ie sequences from GenBank
 * @TODO AP: move chado webservices as a gmod_dbsf_chado submodule
 * @TODO AP: move some functions to gmod_dbsf
 * @TODO AP: cleanup and ensure it is Drupal-compatible (esp. sql)
 */


/*
 * Check if a user is permitted to access the webservices by providing a SID
 *
 * @param SID key, if none, then the current SID for the active user is used
 *
 * @return FALSE on failure, TRUE on success
*/
function gmod_dbsf_chado_sid_perm($sid=NULL) {
    $permission = FALSE;
    // if user is accessing a webservice with the browser without
    // an sid, they can still proceed if they have a recent drupal window logged in
    if (empty($sid)) {
        global $user;
        if (!empty($user)) {
            $sid = $user->sid;
        }
    }

    $session_uid = "SELECT DISTINCT uid from sessions WHERE sid='%s'";
    $res = db_fetch_array(db_query($session_uid,$sid));
    $uid = $res['uid'];
    // inactive session:
    if (empty($uid)) {
        //then check if the permissions of 'Access webservices via Drupal' have been set for anonymous
        // in that case, we allow everyone!
        $perm_array = array();
        $get_perms_sql = "SELECT perm from {permission} where rid IN ((select rid from {role} where name='anonymous user'))";
        $res = db_query($get_perms_sql,$uid);
        while ($data = db_fetch_array($res)) {
            $perm_array=explode(',',$data['perm']);
            foreach ($perm_array as $key=>$value) {
                $perm_array[$key]= trim($value);
            }
            if (in_array('Access webservices via Drupal', $perm_array)) {
                $permission = TRUE;
            }
        }
    }else {
        $perm_array = array();
        $get_perms_sql = "SELECT perm from {permission} where rid IN ((select rid from {role} where name='authenticated user'),(SELECT rid from {users_roles} where uid=%d))";
        $res = db_query($get_perms_sql,$uid);
        while ($data = db_fetch_array($res)) {
            $perm_array=explode(',',$data['perm']);
            foreach ($perm_array as $key=>$value) {
                $perm_array[$key]= trim($value);
            }
            if (in_array('Access webservices via AJAX', $perm_array) ||
                    in_array('Access webservices via Drupal', $perm_array)
            ) {
                $permission = TRUE;
            }
        }
    }
    return $permission;
}



/***************
 * URL parser for chado web services queries using drupal filepath translation
 *
 *	full_url/database/table/id/outputFormat/session-key
 *
 * @param table string the name of the table in the chado database you wish to base your query upon
 * @param id integer the primary key of the row you are interested in
 *        if a string is entered all available results will be shown
 * @param outputFormat string the format you wish the information to be returned in.
 *        Currently supports:
 *          - json (default)
 *          - raw (print_r of the database query for debugging, although json is
 *                 arguably easier for this if you have JSONview installed)
 *          - geoJSON (for selected tables that are traversible to the nd_geolocation table)
 *          - csv
 *          - fasta (for selected tables that are traversible to the feature table)
 *
 *        Would be nice to support:
 *          - xml
 *          - chado.xml
 *
 * @TODO ND: Access Control by Project / User access tables (see db_lookup_perm)
 * @TODO ND: return values only from within selected projects or experiments, accept text value searches within the $primary_key parameter
 *       Note: this may be difficult as not all table attributes are consistent within the chado schema
 *         (except for <table/type>_id. and to a lesser extent table_prop /tableprop)
 * @TODO ND: more flexible search queries / parameters (possibly by accepting HTTP POST/GET requests
 *       There is only so much that can be done through drupal URL routing.
 * @TODO ND: Cascading feature queries to support a wider range of tables
 * @TODO ND: parse outputformat by lastarg.format
 *  .json
 *  .wms
 *  .etc etc
 *
 * @TODO ND: Implement WSDL descriptors
 *
 * @return void This function writes raw output to the browser in order to act as a web service.
 *
 */
function gmod_dbsf_chado_ws_read($db=NULL,$table=NULL, $primary_key=NULL, $outputFormat='json',$sid=NULL) {
    $output = '';

    //debugging
    echo $db.$table.$primary_key=NULL.$outputFormat.$sid;
    //header("HTTP/1.0 403 Forbidden",TRUE,403);
    //header("HTTP/1.0 400 Bad Request",TRUE,400);
    //  support variable passing via path, POST and GET

    //ND Note:URL Router should comply with REST architecture

    /* Breaks REST conventions, removed.
    if (empty($db) && !empty($_POST['db'])){
      $db = $_POST['db'];
    }
    if  (empty($db) && !empty($_GET['db'])){
      $db = $_GET['db'];
    }
    if (empty($table) && !empty($_POST['table'])){
      $table = $_POST['table'];
    }
    if (empty($table) && !empty($_GET['table'])){
      $table = $_GET['table'];
    }
    if (empty($primary_key) && !empty($_POST['id'])){
      $primary_key = $_POST['id'];
    }
    if (empty($primary_key) && !empty($_GET['id'])){
      $primary_key = $_GET['id'];
    }
    if (empty($outputFormat) && !empty($_POST['outputFormat'])){
      $outputFormat = $_POST['outputFormat'];
    }
    if (empty($outputFormat) && !empty($_GET['outputFormat'])){
      $outputFormat = $_GET['outputFormat'];
    }
    if (empty($sid) && !empty($_POST['sid'])){
      $sid = $_POST['sid'];
    }
    if (empty($sid) && !empty($_GET['sid'])){
      $sid = $_GET['sid'];
    }*/


    //force lowercase parameters
    $db             =strtolower($db);
    $table          =strtolower($table);
    $primary_key             =strtolower($primary_key);
    $outputFormat   =strtolower($outputFormat);
    $sid            =strtolower($sid);


    //allows alternative names for tables
    $table          =gmod_dbsf_chado_table_alias($table);

    //Checks
    if (empty($db) || empty($table) || empty($primary_key)) {
        $output = 'Insufficient parameters';
        if ($outputFormat == 'json' || $outputFormat == 'geojson') {
            $output = json_encode($output);
        }
        header("HTTP/1.0 400 Bad Request",TRUE,400);
        echo $output;
        exit;
    }

    if (!is_numeric($primary_key) && $primary_key !=='*') {
        $output ='Unknown parameter';
        if ($outputFormat == 'json' || $outputFormat == 'geojson') {
            $output = json_encode($output);
        }
        header("HTTP/1.0 400 Bad Request",TRUE,400);
        echo $output;
        exit;
    }


    $permitted = gmod_dbsf_chado_sid_perm($sid);
    if (empty($permitted)) {
        $output = 'Failed authentication';
        if ($outputFormat == 'json' || $outputFormat == 'geojson') {
            $output = json_encode($output);
        }
        header("HTTP/1.0 401 Unauthorized",TRUE,401);
        echo $output;
        exit;
    }


    //Limit options to databases that the user is allowed to view
    //database names are configured in drupal settings
    $allowedDBs = array('chado');
    if (!in_array($db,$allowedDBs)) {
        $output = 'Invalid or Unauthorised Database';
        if ($outputFormat == 'json' || $outputFormat == 'geojson') {
            $output = json_encode($output);
        }
        header("HTTP/1.0 403 Forbidden",TRUE,403);
        echo $output;
        exit;
    }
    //end checks


    //OUTPUT FORMATS
    //GeoJSON queries are handled in a different manner than single table db queries
    //as the output must contain Long/Lat (X/Y)
    //and therefore involve a traversal to nd_geolocation
    //$dbresult = array();

    if ($outputFormat === "geojson") {
        //needs a specific query to traverse to nd_geolocation from each table
        //traverses additional tables for more information (optional)
        //$dbresult = gmod_dbsf_chado_chado_lookup_geo_cascade($db,$table, $primary_key);
        $dbresult = gmod_dbsf_chado_chado_lookup_geo($db,$table, $primary_key);
    }
    else {
        $dbresult = gmod_dbsf_chado_db_lookup($db,$table, $primary_key);
    }
    
    exit;
}

/* 
 * Formats output according to parameters
 */
function gmod_dbsf_chado_ws_output_format($output,$format){
    gmod_dbsf_chado_array_remove_nulls($dbresult);

    if (empty($dbresult)) {
        header("HTTP/1.0 400 Bad Request",TRUE,400);
        $output = 'No results';
        if ($outputFormat == 'json' || $outputFormat == 'geojson') {
            $output = json_encode($output);
        }
    }elseif ($outputFormat === 'geojson') {
        header("Content-type: application/json; subtype=geojson");
        header("Pragma: no-cache");
        header("Expires: 0");
        $output = gmod_dbsf_chado_chado_to_GeoJSON($dbresult, $table);
    }elseif ($outputFormat === "json") {
        header("Content-type: application/json");
        header("Pragma: no-cache");
        header("Expires: 0");
        $output = drupal_to_js($dbresult);
    }elseif ($outputFormat === "csv") {
        header("Content-type: application/csv");
        //output as a file
        header("Content-Disposition: attachment; filename=GeospatialGenetics_".$table.$primary_key.".csv");
        header("Pragma: no-cache");
        header("Expires: 0");
        //Do we need to do any cascading here? at present we are just showing a single table with no joins
        $output = gmod_dbsf_chado_array_to_csv($dbresult);
    }elseif ($outputFormat === "fasta") {
        header("Content-type: application/x-fasta");
        header("Content-Disposition: attachment; filename=GeospatialGenetics_".$table.$primary_key.".fasta");
        header("Pragma: no-cache");
        header("Expires: 0");
        //the web service only supports direct requests with a feature ID at present.
        //ie a number or *
        $output = gmod_dbsf_chado_array_to_fasta($dbresult);
    }elseif ($outputFormat === "raw") {
        $output = print_r($dbresult,TRUE);
        //echo "ERROR:html output is implemented within the JavaScript Interface and NOT as a web service.";
    }elseif ($outputFormat === "xml") {
        header("HTTP/1.0 400 Bad Request",TRUE,400);
        $output = "ERROR:xml not yet implemented";
        //@see http://gmod.org/wiki/Chado_XML
        //@see http://doc.bioperl.org/bioperl-dev/ide/emacs/bioperl-mode/test/test-path-1/Bio/SeqIO/chadoxml.html
        //@see http://health2api.com/?p=9543
    }elseif ($outputFormat === "html") {
        header("HTTP/1.0 400 Bad Request",TRUE,400);
        $output = "ERROR:html output is implemented within the JavaScript Interface and NOT as a web service.";
    }elseif ($outputFormat === "table") {
        // use drupal to construct a themed table:
        header("HTTP/1.0 400 Bad Request",TRUE,400);
        $output = "ERROR:html table not yet implemented";
    }else {
        header("HTTP/1.0 400 Bad Request",TRUE,400);
        $output = "ERROR:INVALID FORMAT!";
    }
    echo $output;
}

/**
 * Generatting CSV formatted string from an array.
 * By Sergey Gurevich.
 * @see http://www.codehive.net/PHP-Array-to-CSV-1.html
 */
function gmod_dbsf_chado_array_to_csv($array, $header_row = true, $col_sep = ",", $row_sep = "\n", $qut = '"') {
    $output='';
    if (!is_array($array) or !is_array($array[0])) return false;

    //Header row.
    if ($header_row) {
        foreach ($array[0] as $key => $val) {
            //Escaping quotes.
            $key = str_replace($qut, "$qut$qut", $key);
            $output .= "$col_sep$qut$key$qut";
        }
        $output = substr($output, 1)."\n";
    }
    //Data rows.
    foreach ($array as $key => $val) {
        $tmp = '';
        foreach ($val as $cell_key => $cell_val) {
            //Escaping quotes.
            $cell_val = str_replace($qut, "$qut$qut", $cell_val);
            $tmp .= "$col_sep$qut$cell_val$qut";
        }
        $output .= substr($tmp, 1).$row_sep;
    }
    return $output;
}


/*
 * @TODO AP: implement within gmod_dbsf
 *
 * @TODO ND: Lookup against dbxref/ genbank ID
 GenBank                           gi|gi-number|gb|accession|locus
 EMBL Data Library                 gi|gi-number|emb|accession|locus
 DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
 NBRF PIR                          pir||entry
 Protein Research Foundation       prf||name
 SWISS-PROT                        sp|accession|name
 Brookhaven Protein Data Bank (1)  pdb|entry|chain
 Brookhaven Protein Data Bank (2)  entry:chain|PDBID|CHAIN|SEQUENCE
 Patents                           pat|country|number
 GenInfo Backbone Id               bbs|number
 General database identifier       gnl|database|identifier
 NCBI Reference Sequence           ref|accession|locus
 Local Sequence identifier         lcl|identifier


 *  @see http://code.google.com/p/soframework/source/browse/trunk/main/_modules/biophp/inc/fasta.php?r=53
 *  @see http://en.wikipedia.org/wiki/FASTA_format
 *
 * @return fasta formatted string
*/
function gmod_dbsf_chado_array_to_fasta($array) {
    $output="";

    $rows = count($array); //faster than running the count function call on each repetition in the below for loop
    //Creates a new element in the $feature array for each row returned from the $sql
    if ($rows < 1)
        return null;
    for ($row = 0; $row < $rows; $row++) {
        if (!empty($array[$row]["residues"])) {
            $output.= '>gnl|biosecurity|'.$array[$row]["feature_id"]." ".$array[$row]["uniquename"]."\n";
            $output.= wordwrap($array[$row]["residues"],80,"\n",TRUE);
            $output.= "\n";
        }

    }
    return $output;
}

/**
 * Converts an array of database results containing latitude, longitude and properties into GeoJSON for map display
 * @param type $result
 *   An array containing latitude, longitude and some kind of identifier + metadata The result of your database query
 * @param type $table
 * @return string Returns a GeoJSON formatted string
 *
 * @TODO ND: inject source, sourceid and provider into results
 * @TODO ND: standardise naming and ids
 * @TODO ND: Error handling when lat/long not found in result
 * @TODO ND: Remove GG_Title - no longer needed and too specific to geospatial_genetics
 */
function gmod_dbsf_chado_chado_to_GeoJSON($result, $table) {
    if (empty($result)) {
        return FALSE;
    }
    //The organisation's LSID authority identifier
    $authority = 'biosecurity.cbf.csiro.au';
    //could put db name here?
    $namespace_pre = 'geospatial_genetics_chado_';
    $rows = count($result); //faster than running the count function call on each repetition in the below for loop
    //Creates a new element in the $feature array for each row returned from the $sql
    if ($rows>0) {
        for ($row = 0; $row < $rows; $row++) {
            //$rowval = array_slice($result,$row,1);Might be faster when dealing with many rows/columns
            //$primary_key = $result[$row]["nd_experiment_id"];
            //$nd_geolocation_id = $result[$row]["nd_geolocation_id"];

            $GG_title = 'Undefined';
            switch ($table) {
                case 'nd_geolocation':
                    $GG_title = "Geolocation ID:" . $result[$row]["nd_geolocation_id"];
                    break;
                case 'stock':
                    if (array_key_exists("stock_name",$result[$row])) {
                        $GG_title = "Stock:" . $result[$row]["stock_name"];
                    }elseif (array_key_exists("name",$result[$row])) {
                        $GG_title = "Stock:" . $result[$row]["name"];
                    }
                    break;
                case 'nd_experiment':
                    $GG_title = "Experiment ID:" . $result[$row]["nd_experiment_id"];
                    break;
                case 'project':
                    $GG_title = "Project ID:" . $result[$row]["project_id"];
                    break;
            }

            //add new feature in row $row
            $features[$row] = array(
                    "geometry" => array(
                            "type" => "Point",
                            "coordinates" => array(
                                    $result[$row]["longitude"],
                                    $result[$row]["latitude"]
                            ),
                    ),
                    "type" => "Feature",
                    "properties" => $result[$row],
                    //IDs are LSIDs - URN:LSID:<Authority>:<Namespace>:<ObjectID>[:<Version>]
                    "id" => 'urn:lsid:'.$authority.':'.$namespace_pre.$table.':'.$result[$row][$table."_id"]
            );
            $features[$row]["properties"]["GG_title"] = $GG_title; //better solution would be to have $GG_title as the first element in the array?
            $features[$row]["properties"]["sourcetable"]=$table;

            //Add in a feature id if it can be linked
            if (array_key_exists("stock_id",$result[$row])) {
                $sql="  SELECT
                          feature.feature_id
                        FROM
                          public.feature,
                          public.feature_genotype,
                          public.genotype,
                          public.stock_genotype,
                          public.stock
                        WHERE
                          feature.feature_id = feature_genotype.feature_id AND
                          feature_genotype.genotype_id = genotype.genotype_id AND
                          genotype.genotype_id = stock_genotype.genotype_id AND
                          stock.stock_id = stock_genotype.stock_id AND
                          stock.stock_id =". $result[$row]["stock_id"];
                //If the stock ID can be matched to a feature then return the matching id in the output
                $result2 = gmod_dbsf_query_execute('chado', $sql, null);
                if (!empty($result2[0]["feature_id"])) {
                    $features[$row]["properties"]["feature_id"]=$result2[0]["feature_id"];
                }
            }
            //Stack sort
            $features[$row]["properties"]=gmod_dbsf_chado_array_key_to_top('latitude',$features[$row]["properties"]);
            $features[$row]["properties"]=gmod_dbsf_chado_array_key_to_top('longitude', $features[$row]["properties"]);
            $features[$row]["properties"]=gmod_dbsf_chado_array_key_to_top('GG_title', $features[$row]["properties"]);

        }
    }
    //wrap features as a FeatureCollection
    $output = array(
            "type" => "FeatureCollection",
            "features" => $features
    );
    header('Content-Type: application/json');
    echo drupal_to_js($output);
    exit;
}

/**
 * Returns chado data associated with a geolocation
 *
 * If no ID is specified all results are listed
 *
 * @param $table = the table you would like to search for an id
 *      Currently supported tables are :
 *      - nd_geolocation
 *      - nd_experiment
 *      - stock
 * Defaults to nd_experiment if no other valid table is specified.
 *
 * @param $primary_key = the primary key you wish to look up in the specifed table.
 * (ie. "1" would look for nd_experiment.nd_experiment ID = 1).
 * "*" Can be used to list all results.
 *
 * @TODO ND: test to see if string wrapping operations as in chadoToGeoJSON.php are faster/use less memory
 * than calling all values into nested arrays
 * @TODO ND: support for more tables
 * @TODO ND: pull in props
 * @TODO ND: Genericise similar to gmod_dbsf_chado_db_lookup
 * @TODO ND: remove null columns
 *
 * @return array returns an array of form long(X)/lat(Y)/ID/DATUM/[METADATA] of 'features' of the format GeoJSON to display nd_geolocations specified (or all geolocations if not specified)
 *
 */
function gmod_dbsf_chado_chado_lookup_geo($db='chado',$table=NULL, $primary_key=NULL, $alldata=TRUE) {
    $accessCheck=gmod_dbsf_chado_db_lookup_perm($db,$table,$primary_key);
    if (!($accessCheck==='ACCESS_GRANTED')) {
        return $accessCheck;
    }

    //Tables currently supported for geospatial queries (must be able to traverse to nd_geolocation)
    $supportedTables = array(
            "nd_geolocation",
            "nd_experiment",
            'project',
            "stock"
    );
    if (!in_array($table, $supportedTables)) {
        return FALSE;
    }

    //SQL QUERIES
    // AP: rewrite these to be more STND & efficient
    // and return only what is needed
    $sql ='';
    switch ($table) {
        case 'stock':
            $sql = "SELECT *
            FROM
              public.stock as s,
              public.nd_experiment_stock as es,
              public.nd_experiment as e,
              public.nd_geolocation as g
            WHERE
	          g.nd_geolocation_id = e.nd_geolocation_id AND
    	      e.nd_experiment_id = es.nd_experiment_id AND
              es.stock_id = s.stock_id
                ";
            break;
        case 'nd_experiment':
            $sql = "SELECT *
            FROM
              public.stock as s,
              public.nd_experiment_stock as es,
              public.nd_experiment as e,
              public.nd_geolocation as g
            WHERE
	          g.nd_geolocation_id = e.nd_geolocation_id AND
    	      e.nd_experiment_id = es.nd_experiment_id AND
              es.stock_id = s.stock_id
            ";
            break;
        case 'project':
            $sql = 'SELECT *
        	FROM
              public.stock as s,
              public.nd_experiment_stock as es,
        	  public.nd_experiment_project as ep,
              public.nd_experiment as e,
        	  public.nd_geolocation as g
            WHERE
              g.nd_geolocation_id = e.nd_geolocation_id AND
    	      e.nd_experiment_id = es.nd_experiment_id AND
              e.nd_experiment_id = ep.nd_experiment_id AND
    	      es.stock_id = s.stock_id
              ';
            break;
            //default to nd_geolocation
            $sql = "SELECT * FROM public.nd_geolocation";
    };

    $result =array();
    if (is_numeric($primary_key)) {
        //if a numeric ID is specified then only return the record for that specific ID
        switch ($table) {
            case 'stock':
            //Change SQL to select stocks with geolocations and retrieve all NOT NULL data
                $sql .= " AND s.stock_id=%d "; // parameterized . $primary_key;
                break;
            case 'nd_experiment':
                $sql .= " AND e.nd_experiment_id=%d ";
                break;
            case 'project':
                $sql .= " AND ep.project_id=%d ";
                break;
                //default to nd_geolocation
                $sql .= " AND g.nd_geolocation_id=%d ";

        }
        $result = gmod_dbsf_query_execute($db, $sql, $primary_key);
    }
    elseif ($alldata) {
        //if wildcards are allowed then return all results
        $result = gmod_dbsf_query_execute($db, $sql);
    }
    //if empty result then show error
    if (empty($result)) {
        return FALSE;
    }
    return $result;
}


/**
 * Returns chado data associated with a geolocation, cacading information from other tables
 *
 * UNUSED
 *
 * If no ID is specified all results are listed
 *
 * @param $table = the table you would like to search for an id
 *      Currently supported tables are :
 *      - nd_geolocation
 *      - nd_experiment
 *      - stock
 * Defaults to nd_experiment if no other valid table is specified.
 *
 * @param $primary_key = the primary key you wish to look up in the specifed table.
 * (ie. "1" would look for nd_experiment.nd_experiment ID = 1).
 * "*" Can be used to list all results.
 *
 * @return array returns an array of 'features' of the format GeoJSON to display nd_geolocations specified (or all geolocations if not specified)
 *
 * @TODO ND: test to see if string wrapping operations as in chadoToGeoJSON.php are faster/use less memory
 * than calling all values into nested arrays
 * @TODO ND: support for more tables
 * @TODO ND: pull in props
 * @TODO ND: Genericise similar to gmod_dbsf_chado_chado
 *
 *
 */
function gmod_dbsf_chado_chado_lookup_geo_cascade($db='chado',$table=NULL, $primary_key=NULL, $alldata=TRUE) {
    $accessCheck=gmod_dbsf_chado_db_lookup_perm($db,$table,$primary_key);
    if (!($accessCheck==='ACCESS_GRANTED')) {
        return $accessCheck;
    }

    //Tables supported for geospatial queries (must be able to traverse to nd_geolocation)
    $supportedTables = array(
            "nd_geolocation",
            "nd_experiment",
            "stock"
    );
    if (!in_array($table, $supportedTables)) {
        return FALSE;
    }

    //**SQL QUERIES
    if ($table === "stock") {
        $sql = "SELECT
              stock.uniquename AS stock_uniquename,
              nd_geolocation.longitude,
              nd_geolocation.latitude,
              nd_geolocation.nd_geolocation_id,
              nd_experiment.nd_experiment_id,
              nd_experiment.type_id nd_experiment_type_id,
              nd_experiment_stock.type_id AS experiment_stock_type_id,
              stock.stock_id,
              stock.name AS stock_name,
              stock.description AS stock_description,
              stock.type_id AS stock_type_id,
              organism.abbreviation AS organism_abbreviation,
              organism.genus AS organism_genus,
              organism.species AS organism_species,
              organism.common_name AS organism_common_name,
              organism.comment AS organism_comment,
              organism.organism_id,
              nd_geolocation.description AS nd_geolocation_description,
              nd_geolocation.geodetic_datum
            FROM
              public.nd_experiment,
              public.nd_geolocation,
              public.nd_experiment_stock,
              public.stock,
              public.organism
            WHERE
              nd_experiment.nd_experiment_id = nd_experiment_stock.nd_experiment_id AND
              nd_geolocation.nd_geolocation_id = nd_experiment.nd_geolocation_id AND
              nd_experiment_stock.stock_id = stock.stock_id AND
              organism.organism_id = stock.organism_id";
    }
    else {
        $sql = "SELECT
              nd_geolocation.nd_geolocation_id,
              nd_geolocation.longitude,
              nd_geolocation.latitude,
              nd_geolocation.description AS geolocation_description,
              nd_experiment.nd_experiment_id,
              nd_experiment.type_id AS cvterm_experiment
            FROM
              public.nd_experiment,
              public.nd_geolocation,
              public.cvterm
            WHERE
              nd_experiment.nd_geolocation_id = nd_geolocation.nd_geolocation_id AND
              cvterm.cvterm_id = nd_experiment.type_id";
    }


    //return $q_nd_geolocation_id if specified or all geolocations if not numeric (ie null/invalid);
    //text based searches by cvterm?
    if (is_numeric($primary_key)) {
        if ($table === "nd_experiment") {
            $sql .= " AND nd_experiment.nd_experiment_id=%d "; // parameterized . $primary_key;
        }
        elseif ($table === "nd_geolocation") {
            $sql .= " AND nd_geolocation.nd_geolocation_id=%d "; // parameterized . $primary_key;
        }
        elseif ($table === "stock") {
            //Change SQL to select stocks with geolocations and retrieve all NOT NULL data
            $sql .= " AND stock.stock_id=%d "; // parameterized . $primary_key;
        } //else default to experiment
        else {
            $sql .= " AND nd_experiment.nd_experiment_id=%d ";
        }
        $result = gmod_dbsf_query_execute($db, $sql, $primary_key);
    }
    elseif ($alldata) {
        //if wildcards are allowed then return all results
        $result = gmod_dbsf_query_execute($db, $sql);
    }
    else {
        $result = NULL;
    }
    if (empty($result)) {
        return array('Invalid'=>"primary key");
    }
    return $result;

}

/* *****************************************************
 *
 * Look up any value from a chado database table
 * @param table string
 *   The table you wish to work with
 * @param id integer
 *   The primary key serial of the row you are looking for
 * @param boolean allowWildcard
 *   Allow the whole table to be returned, defaults to true
 * @return array (dbresult) as [row][columnname][value] or <error message>
 *
 * @TODO ND: add authority:namespace:id or some such (could be provider, source, id too) to output
 *
 * @TODO ND: re implement access check as a callback to strip protected columns and null columns from the output
 *       or alternatively implement db_lookup_strip
 *
 * @TODO AP: move to gmod_dbsf_chado
*/
function gmod_dbsf_chado_db_lookup($db='chado',$table=NULL, $primary_key=NULL, $alldata=TRUE) {
    $accessCheck=gmod_dbsf_chado_db_lookup_perm($db,$table,$primary_key);
    if (!($accessCheck==='ACCESS_GRANTED')) {
        return $accessCheck;
    }
    /*
    column level access control
    if ($table==='feature') {
        //do not return residues if the table 'feature' is selected. These can be accessed via fasta
        $select = "feature_id, dbxref_id, organism_id, name, uniquename, seqlen, md5checksum, type_id, is_analysis, is_obsolete, timeaccessioned, timelastmodified";
    }
    */
    $select='*';
    if ($primary_key == 'count') {
        $sql = "SELECT COUNT(*) FROM " . $table;
        $result = gmod_dbsf_query_execute($db, $sql,$table);
    }elseif (is_numeric($primary_key)) {
        //if the id is a number, query the DB for it.
        $sql = "SELECT ".$select." FROM {" . $table . "} WHERE " .$table . "_id=%d ";
        $result = gmod_dbsf_query_execute($db, $sql, $primary_key);
    }elseif ($alldata) {
        $sql = "SELECT ".$select." FROM {" . $table.'}';
        $result = gmod_dbsf_query_execute($db, $sql);
    }else {
        $result = FALSE;
    }
    return $result;
}
/***********************************************
 * Access control for database queries
 *
 * @return array listing any access control errors
 * @TODO ND: change error messages to DPMS
 * @TODO ND: Access control
 *   - Global wildcard access or wildcard for certain dbs only.
 *   - [DONE] Web services level access control in gmod_dbsf_chado_menu() parser
 *              Drupal Usergroup (ie. auth_user, admin etc)
 *   - [DONE] Table level access control
 *              Hardcoded in gmod_dbsf_chado_chado_lookup
 *     [TODO] ND: access control by user group
 *   - [TODO] ND: Column level access control ie. dissalow access to certain columns
 *            Lookup tables, rows within the admin panel, ie. (if in column in table, disallowed columns delete from array prior to output
 *               NOTE: will need to encrypt DB communications for this method to be secure.
 *               Alternatively you can build in select <column> from table but a lot of work.
 *  -  [TODO] ND: Access by member of chado project
 *             Can configure by project membership but will need lookup tables for user/project and project/table/id
 *             Ie. if a member of a project submits an item as part of a project it will go into project/table/id
 *                 Alternatively you can partition groups into their own databases
 *  -  [TODO] ND: Access by drupal UID
 *             Can configure by project membership but will need lookup tables - probably easier to match drupal UID to project and cascade
 * -   [TODO] ND: Serials Only
 *
 * -   [TODO] ND: Callback to run a query function ie this(dbQueryFunction, args) to handle access control and output sanity
 *
 * @TODO ND: refactor to minimise database calls and call $available tables into a global variable on init.
 * @TODO ND: cron to clean up database serials so that they run from 0 - maxnum
 *
 * @returns array showing an error or 'ACCESS_GRANTED"
 *    @TODO ND: Return array of ALLOW/DENY columns within table
 *          OR specify function as callback and strip out illegal results
 *    @TODO ND: Standardise output
 */
function gmod_dbsf_chado_db_lookup_perm($db=NULL,$table=NULL, $primary_key=NULL) {
    //not implemented yet - if false it needs to check $primary_key for non numeric if false and deny or rewrite char to NULL.
    //May not want to implement as it may intefere with cascades.
    $alldata=TRUE;

    //disallow access to the following restricted tables
    $restrictTables=FALSE;
    $restrictedTables = array(
            "feature",
    );
    if (in_array($table, $restrictedTables) && $restrictTables===TRUE) {
        header('Content-Type: application/json');
        echo drupal_to_js(array($table => 'Table restricted'));
        exit;

    }
    //does not stop caller for some reason

    //Support only the following tables if specifyAllowedTables=TRUE:
    $specifyAllowedTables=TRUE;
    $allowedTables = array(
            'nd_geolocation',
            'nd_experiment',
            'nd_experiment_stock',
            'nd_experiment_project',
            'stock',
            'feature',
            'project',
            'organism',
            'cvterm',
            'cv','db','dbxref'
    );
    if (!in_array($table, $allowedTables) && $specifyAllowedTables===TRUE) {
        header('Content-Type: application/json');
        echo drupal_to_js(array($table => 'Table not allowed'));
        exit;
    }

    //check if table appears in database
    $availableTables = gmod_dbsf_query_execute($db,
            "SELECT table_name FROM information_schema.tables where information_schema.tables.table_schema='public'",
            $table);
    $result = gmod_dbsf_query_execute(
            $db, "SELECT table_name FROM information_schema.tables where information_schema.tables.table_schema='public' AND table_name='%s'",
            $table);
    if (empty($result)) {
        return 'Table not found';
    }
    return "ACCESS_GRANTED";
}

/**
 * Returns the name of a table if aliased
 *
 * @param string $tableName
 *
 * @return string The correct name of the table
 */
function gmod_dbsf_chado_table_alias($tableName) {
    //ALIASES
    if ($tableName === "nd_experiment" || $tableName === "experiment" || $tableName === "experiments" || $tableName === "exp") {
        return "nd_experiment";
    }
    if ($tableName === "nd_geolocation" || $tableName === "geo") {
        return "nd_geolocation";
    }
    if (stristr($tableName,"type")) {
        return "cvterm";
    }
    elseif ($tableName === "stock" || $tableName === "sample" || $tableName === "specimen") {
        return "stock";
    }
    else {
        return $tableName;
    }
}

/* Moves an array key => value pair to the top of an array if found
 * ie. pull the row from the array and then stack.push
 *
 * @param $key string OR integer the array key we are looking for
 *
 * @param $array array the input array
 *
 * @return an array with the new values
*/
function gmod_dbsf_chado_array_key_to_top($key, $array) {
    if (isset($array[$key])) {
        //set first element to the row corresponding to $key
        $newarray = array($key => $array[$key]);
        //remove the old row from the array
        unset($array[$key]);
        //append the rest of the array after the row we have moved to the first element of the new array
        return array_merge($newarray,$array);
    }
    //dpm($key.' not found in '.$array);
    return $array; //return old array if key not found
}

/*
 * Steps recursively through an array and removes null values
*/
function gmod_dbsf_chado_array_remove_nulls(&$arr) {
    if (empty($arr)) {
        return FALSE;
    }
    elseif (!is_array($arr)) {
        return $arr;
    }
    foreach($arr as $key => &$item) {
        if(empty($arr[$key])) {
            unset($arr[$key]);
        }elseif(is_array($arr[$key])) {
            gmod_dbsf_chado_array_remove_nulls($arr[$key]);
        }
    }
}




