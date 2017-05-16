<?php
// $Id$

/**
 * @file
 * include file for viewing features and creating genes pages
 */

function genes4all_explore_lookup_menu(){
  $items = array();
  $include_path = drupal_get_path('module','genes4all').'/modules/genes4all_explore/includes/';
  $items['genes4all/info/feature-ajax/assembly'] = array(
    'title' => 'Feature AJAX assembly',
    'page arguments' => array(4),
    'page callback' => 'genes4all_explore_page_feature_info_assembly',
    'access arguments' => array('Access genes4all_explore'),
    'file' => 'genes4all_explore_lookup.inc',
    'file path' => $include_path,
    'type' => MENU_CALLBACK,
  );
  $items['genes4all/info/feature-ajax/annotation'] = array(
    'title' => 'Feature AJAX annotation',
    'page arguments' => array(4),
    'page callback' => 'genes4all_explore_page_feature_info_annotation',
    'access arguments' => array('Access genes4all_explore'),
    'file' => 'genes4all_explore_lookup.inc',
    'file path' => $include_path,
    'type' => MENU_CALLBACK,
  );
  $items['genes4all/info/feature-ajax/sequence'] = array(
    'title' => 'Feature AJAX sequence',
    'page arguments' => array(4),
    'page callback' => 'genes4all_explore_page_feature_info_sequence',
    'access arguments' => array('Access genes4all_explore'),
    'file' => 'genes4all_explore_lookup.inc',
    'file path' => $include_path,
    'type' => MENU_CALLBACK,
  );
  return $items;
}


function genes4all_explore_page_checkout_feature($accession = NULL) {
  if (empty($accession)) {
    return '';
  }
  $return_text = '';
  $feature_data = array();
  $accession = check_plain(trim($accession));
  if (is_numeric($accession)) {
    drupal_goto('http://www.ncbi.nlm.nih.gov/sites/gquery', "term=$accession");
  }
  //est2assembly hacks
  if (strpos($accession,':ORF') !==FALSE){
    $accession=substr($accession,0,strpos($accession,':ORF'));
  }elseif(strpos($accession,'_0') !==FALSE){
    $accession=substr($accession,0,strpos($accession,'_0'));
  }elseif (strpos($accession,'_assembly') !==FALSE){
    $accession=substr($accession,0,strpos($accession,'_assembly'));
  }
  //BOOKMARK DEVEL; comment out following to enable development
  $feature_data =cache_get($accession,'cache_genes4all_explore_feature');
  if (empty($feature_data)) {
    $feature_data = gmod_dbsf_db_execute('chado', 'genes4all_explore_feature_data', $accession);
    cache_set($accession, $feature_data , 'cache_genes4all_explore_feature', time()+7*3600*24);
  }else{
    $feature_data=$feature_data->data;
  }
  if (empty($feature_data)) {
    $ncbi_link = l("the NCBI Entrez database", 'http://www.ncbi.nlm.nih.gov/sites/gquery', array('query' => array('term' => $accession)));
    $uniprot_link = l("at the UniProt database", 'http://www.uniprot.org/uniprot', array('query' => array('sort'=>'score','query' => $accession)));
    $return_text = t("<p>I'm sorry, but there are no data stored for feature %s</p>",array('%s'=>$accession));
    $return_text .= t("<p>You can try to find the accession at:");
    $return_text .="<br>".$ncbi_link;
    $return_text .="<br>".$uniprot_link;
    $return_text .="</p>";
    return t($return_text);
  }
  $access = gmod_dbsf_db_execute('chado', 'gmod_dbsf_chadoapi_check_restriction', 'select', 'feature', $feature_data['id']);

  if ($access != TRUE) {
    $return_text =t("Sorry, but this is part of a private dataset and you are not authorized to view this feature. You can try to !contact_link us to request permission and we will inform the owner of the data of your request.",array('!contact_link'=>l('contact','contact')));
    return $return_text;
  }
  //method to get all properties and CVs: $feature_properties = gmod_dbsf_db_execute('chado','gmod_dbsf_get_feature_cvterm_and_prop',$feature_id);
  return (array('accession'=>$accession,'text'=>$return_text,'data'=>$feature_data));
}

/*
 * The information page. Give it a feature id and it tell you all about it
 */
function genes4all_explore_page_feature_info($accession = NULL) {
  if (empty($accession)){
    $accession = $_GET['feature_id'];
    if (empty($accession)) {
      return t("I'm sorry, you did not specify a feature");
    }
  }
  $accession_data = genes4all_explore_page_checkout_feature($accession);
  if (!is_array($accession_data)){
    return t($accession_data);
  }
  $accession = $accession_data['accession'];
  $return_text = $accession_data['text'];
  $feature_data = $accession_data['data'];
  if (empty($feature_data)) {
    drupal_set_message(t('Warning, no data could be retrieved for this feature.'), 'error');
    return $return_text;
  }
  elseif (!empty($feature_data['obsolete'])) {
    drupal_set_message(t('Warning, this feature is marked as obsolete.'), 'warning');
  }
  $return_basic      = genes4all_explore_page_feature_info_basic($feature_data);
  $return_graphs     = genes4all_explore_page_feature_info_plots($accession);
  $return_annotation = genes4all_explore_page_feature_info_annotation($accession);
  $page_tabs = array('results' => array('#type' => 'tabset'));
  $page_tabs['results'][] = array(
    '#type' => 'tabpage',
    '#title' => 'Basic',
    '#content' => $return_basic,
    '#weight' => -1,
    '#selected' => TRUE,
  );
  if (!empty($return_annotation)) {
    $page_tabs['results'][] = array(
      '#type' => 'tabpage',
      '#title' => 'Annotation',
      '#content' => $return_annotation,
//not work with table      '#ajax_url' => url('genes4all/info/feature-ajax/annotation/'.$feature_data['accession']),
      '#weight' => 1,
    );
  }
  /* I've decided to remove the assembly as no-one seems to be using it (also it no longer makes sense
   * with RNA-seq data...
   if ($feature_data['type'] == 'contig') {
   $page_tabs['results'][] = array(
   '#type' => 'tabpage',
   '#title' => 'Assembly',
   '#ajax_url' => url('genes4all/info/feature-ajax/assembly/'.$feature_data['accession']),
   '#weight' => 15,
   );
   }
   */
  if (!empty($feature_data['residues'])) {
    // even though we already have residues? This function recalls it in order to format
    // it and PostGres should have cached the result?
    $page_tabs['results'][] = array(
      '#type' => 'tabpage',
      '#title' => 'Sequence',
      '#ajax_url' => url('genes4all/info/feature-ajax/sequence/'.$feature_data['accession']),
      '#weight' => 10,
    );
  }
  if (!empty($return_graphs)) {
    $page_tabs['results'][] = array(
      '#type' => 'tabpage',
      '#title' => 'Graphs',
      '#description'=> 'Graphs along the sequence are dynamically generated for coverage (optionally per library), presence of BLAST hits and if provided, Next-Gen data. Hover over the graph to highlight the value.',
      '#content' => $return_graphs, // ajax does not work
      '#weight' => 5,
    );
  }
  $return_text .= tabs_render($page_tabs);
  return $return_text;
}

function genes4all_explore_page_feature_info_assembly($accession=NULL){
  if (empty($accession)){return;}
  $feature_data =cache_get($accession,'cache_genes4all_explore_feature');
  $feature_data=$feature_data->data;
  if (empty($feature_data)){return;}
  echo $feature_data['assembly'];
}

function genes4all_explore_page_feature_info_related($feature_data){
  $related_features = "";
  if (empty($feature_data['links'])){return;}
  foreach ($feature_data['links'] as $type => $id_array) {
    if (empty($id_array)) {
      continue;
    }
    if (!is_array($id_array)) {
      $related_features .= "<br>A $type is available at ". l($id_array, 'genes4all/info/feature', array('query' => array('feature_id' => $id_array)));
    }
    else {
      foreach ($id_array as $id) {
        if (empty($id)) {
          continue;
        }
        $related_features .= "<br>A $type is available at ". l($id, 'genes4all/info/feature', array('query' => array('feature_id' => $id)));
      }
    }
  }
  return $related_features;
}


function genes4all_explore_page_feature_info_sequence($accession=NULL){
  if (empty($accession)){return;}
  $feature_data =cache_get($accession,'cache_genes4all_explore_feature');
  $feature_data=$feature_data->data;
  if (empty($feature_data)){return;}
  $sequence = '';
  if (!empty($feature_data['residues'])) {
    $sequence = wordwrap($feature_data['residues'], 80, '<br/>', TRUE);
  }
  if (!empty($sequence)){
    echo "<h2>Sequence</h2><pre>>gnl|Insectacentral|".$feature_data['accession']."<br/>$sequence</pre>";;
  }
}

/*
 * TODO Check if it is all ok
 */
function genes4all_explore_page_feature_info_curated($feature_name){
  $submitted_text = '';
  if (db_table_exists('genes4all_approval_queue')) {
    //Drupal-specific SQL
    $submitted_sql    = "SELECT u.uid,aq.submitter_name as who, aq.date_added as when , aq.date_curated as curated from {genes4all_approval_queue} as aq" . " JOIN {users} u ON u.name=aq.submitter_name where aq.uniquename='%s'";
    $subs_res         = db_fetch_array(db_query($submitted_sql, $feature_name));
    $subs_res['when'] = str_replace(' ', ' at ', preg_replace('/:\d\d\.\d.+$/', '', $subs_res['when']));
    if (!empty($subs_res['curated'])) {
      $subs_res['curated'] = str_replace(' ', ' at ', preg_replace('/:\d\d\.\d.+$/', '', $subs_res['curated']));
    }
    $submitted_text .= !empty($subs_res['who']) ? ' It was submitted by '. l($subs_res['who'], 'user', array('query' => $subs_res['uid'])) .' on '. $subs_res['when'] .'.' : ' It was batch submitted.';
    if (!empty($subs_res['who'])) {
      $submitted_text .= !empty($subs_res['curated']) ? ' It was approved on '. $subs_res['curated'] .'.' : ' It has not been approved yet.';
    }
  }
  return $submitted_text;
}
function genes4all_explore_page_feature_info_basic($feature_data=NULL){
  if (empty($feature_data)){return;}
  //LINKS to Gbrowse and contig/orf/polypetide
  $related_features =genes4all_explore_page_feature_info_related($feature_data);
  $browser_links=genes4all_explore_page_feature_info_browse($feature_data);
  //Curation info
  $submitted_text = genes4all_explore_page_feature_info_curated($feature_data['uniquename']);
  // Species info
  $species_info = $feature_data['genus'] .' '. $feature_data['species'];
  if ($feature_data['org_taxid']){
    $ncbitax_url = l($species_info,
    'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi',
    array(
      'attributes' => array('target' => '_blank'), 'external' => TRUE,
      'query' => array('id' => $feature_data['org_taxid']),
    )
    );
  }else{
    $ncbitax_url = $species_info;
  }
  $lastmodified      = preg_replace('/:\d+\.\d+$/', '', $feature_data['lastmodified']);
  $lastmodified      = str_replace(' ', ' at ', $lastmodified);

  $return_basic = "<p>Feature ".$feature_data['accession']." is a ".$feature_data['type'];
  if (!empty($feature_data['seqlength_text'])) {
    $return_basic .= " which is ".$feature_data['seqlength_text']." long";
  }
  if (!empty($feature_data['est_count'])) {
    $return_basic .= " and composed of ".$feature_data['est_count']." reads";
  }
  if (!empty($feature_data['prediction_method'])) {
    $return_basic .= " and was predicted using the ".$feature_data['prediction_method']." software";
  }
  $return_basic .= '.';
  if (!empty($feature_data['snp_info'])) {
    $return_basic .= $feature_data['snp_info'];
  }

  if (!empty($feature_data['species_common'])) {
    $species_common = '('. $feature_data['species_common'] .')';
  }

  $return_basic .= "<br> It was derived from $ncbitax_url";
  if (!empty($species_common)){
    $return_basic .=" $species_common";
  }
  $return_basic .='.';
  $return_basic .= "<br> $submitted_text It was last modified on $lastmodified.<br>";
  // curation link
  $return_basic .=l(
  theme_image(drupal_get_path('module', 'genes4all') .'/images/curate_this.png', 'Curate this feature', 'Curate this feature')
  ,'genes4all/feature/curate/'.$feature_data['accession'],array('html'=>TRUE,'alias'=>TRUE)
  );
  $return_basic .= '</p>';

  $themed_lib_table = genes4all_explore_page_feature_info_library($feature_data['id']);
  if (!empty($themed_lib_table)) {
    $return_basic .= $themed_lib_table;
  }
  if (!empty($browser_links) || !empty($related_features)){
    $return_basic .= "<h2>Links</h2><p>";
    if (!empty($browser_links['gbrowse'])) {
      $return_basic .= "You can view its ".$browser_links['gbrowse'].".";
    }elseif(!empty($browser_links['jbrowse'])){
      $return_basic .= "You can view its ".$browser_links['jbrowse'].".";
    }
    $return_basic .= $related_features."</p>";
  }
  return $return_basic ;
}

/**
 * Prepare plots for feature graphing with the Dojo Toolkit
 *
 * @param $accession
 *
 * @return HTML of graph
 *
 */
function genes4all_explore_page_feature_info_plots($accession=NULL){
  $return_graphs='';
  if (empty($accession)){return;}
  $javascript = '';
  $feature_data =cache_get($accession,'cache_genes4all_explore_feature');
  if (!empty($feature_data)){
    $feature_data=$feature_data->data;
  }
  // Feature data must have been in the cache already
  if (empty($feature_data)){
    drupal_set_message(t('Error EL1. This should not have happened. Please contact a system administrator.'),'error');
    return;
  }
  $module_path = drupal_get_path('module', 'genes4all');
  if (!empty($feature_data['coverage_array']) || !empty($feature_data['blast_coverage_array'])
  || !empty($feature_data['nextgen_coverage_array'][$feature_data['uniquename']])
   || !empty($feature_data['lib_coverage_array'])) {
    drupal_add_css($module_path.'lib/dojo/resources/dojo.css', 'module', 'all');
    drupal_add_css($module_path.'/lib/dijit/themes/tundra/tundra.css', 'module');
    drupal_add_css($module_path.'/lib/dojox/form/resources/RangeSlider.css','module');
    drupal_add_js($module_path.'/lib/dojo/dojo.js');

    $cov_graph='';
    $cov_javascript = '';
    if (!empty($feature_data['coverage_array'])) {
      $cov_javascript = 'chart1.addPlot("coverage", {type: "Lines",markers: true, tension:"X",hAxis: "x", vAxis: "coverage y"});';
      $max = $feature_data['coverage_array_max'];
      $max += ceil($max * 0.10);
      $cov_graph.='
    chart1.addAxis("coverage y", {vertical: true, includeZero: true, fixUpper: "major", max: '.$max.'});
    chart1.addSeries("Coverage", '.drupal_to_js($feature_data['coverage_array']).',{plot: "coverage", stroke: {color: "blue", width: 2},fill: "lightblue"});
    var animcova = new dc.action2d.Magnify(chart1,"coverage",{scale: 4});
    var animcovb = new dc.action2d.Tooltip(chart1,"coverage");
    ';
    }
    $blast_javascript = '';
    $blast_graph = '';
    if (!empty($feature_data['blast_coverage_array'])) {
      $blast_javascript = 'chart1.addPlot("blast", {type: "Lines",markers: true, tension:"X",hAxis: "x", vAxis: "blast y"});';
      $max = $feature_data['blast_coverage_array_max'];
      $max += ceil($max * 0.10);
      $blast_graph.='
    chart1.addAxis("blast y", {vertical: true, includeZero: true,leftBottom: false, fixUpper: "major", max: '.$max.'});
    chart1.addSeries("BLAST results on 2ndary axis", '.drupal_to_js($feature_data['blast_coverage_array']).',{plot: "blast", stroke: {color: "red", width: 2},fill: "pink"});
    var animblasta = new dc.action2d.Magnify(chart1,"blast",{scale: 4});
    var animblastb = new dc.action2d.Tooltip(chart1,"blast");
    ';
    }
    $lib_javascript = '';
    $lib_graph='';
    if (!empty($feature_data['lib_coverage_array'])) {
      $k=1;
      foreach ($feature_data['lib_coverage_array'] as $lib_name => $lib_data){
        if (empty($lib_name)){continue;}
        $lib_javascript .='
   		chart1.addPlot("'.$k.'_lib_coverage", {type: "Lines",markers: true, tension:"X",hAxis: "x", vAxis: "coverage y"});
   		';
        $lib_graph.='
        chart1.addSeries("Coverage of '.$lib_name.'", '
        .drupal_to_js($lib_data).',{plot: "'.$k.'_lib_coverage", stroke: {width: 2}});
        var animlib'.$k.'a = new dc.action2d.Magnify(chart1,"'.$k.'_lib_coverage",{scale: 4});
        var animlib'.$k.'b = new dc.action2d.Tooltip(chart1,"'.$k.'_lib_coverage");
       ';
        $k++;
      }
    }

    // this is not the proper way to pass variables, but it works; feel free to fix
    $javascript = '
  dojo.require("dojox.charting.Chart2D");
  dojo.require("dojox.charting.themes.Julie");
  dojo.require("dojox.charting.action2d.Magnify");
  dojo.require("dojox.charting.action2d.Tooltip");
  dojo.require("dojox.charting.widget.Legend");
  var dc = dojox.charting;

  makeCharts = function() {
        var chart1 = new dojox.charting.Chart2D("coverage_graphs");
        '
        .$cov_javascript
        .$lib_javascript
        .'
		chart1.addPlot("grid", {hAxis: "x", vAxis: "coverage y",type: "Grid",hMajorLines: false,hMinorLines: false,vMajorLines: true,vMinorLines: false});
        '.$blast_javascript
        .'
        chart1.addAxis("x",{includeZero: true});
        '
        .$cov_graph
        .$lib_graph
        .$blast_graph
        .'
        chart1.setTheme(dojox.charting.themes.Julie);
        chart1.render();
        var legend1 = new dojox.charting.widget.Legend({chart: chart1, horizontal: false}, "legend1");
	';
  //TODO nextgen
    if (!empty($feature_data['nextgen_coverage_array']) ){
      $nextgen_javascript ='';
      $nextgen_graph ='';
      $max =$feature_data['nextgen_coverage_array']['all']['max'];
      $max += ceil($max * 0.10);
      $nextgen_javascript = '
      chart2.addAxis("nextgen y", {vertical: true, includeZero: true, max: '.$max.'});
      ';
      $k=1;
      foreach ($feature_data['nextgen_coverage_array'] as $aln_description=>$nextgen_data){
        if (empty($nextgen_data['values'])){continue;}
        $nextgen_javascript.='
        chart2.addPlot("'.$k.'_nextgen", {type: "Lines",markers: true, tension:"X",hAxis: "x", vAxis: "nextgen y"});
   		';
        $nextgen_graph .='
        chart2.addSeries("Next-Gen coverage of '.$aln_description.'", '
        .drupal_to_js($nextgen_data['values']).',{plot: "'.$k.'_nextgen", stroke: {width: 2}});
        var animngs'.$k.'a = new dc.action2d.Magnify(chart2,"'.$k.'_nextgen",{scale: 4});
        var animngs'.$k.'b = new dc.action2d.Tooltip(chart2,"'.$k.'_nextgen");
       ';
        $k++;
      }
      $javascript.='var chart2 = new dojox.charting.Chart2D("ngs_graphs");'
        .$nextgen_javascript
        .'
		chart2.addPlot("grid", {hAxis: "x", vAxis: "nextgen y",type: "Grid",hMajorLines: false,hMinorLines: false,vMajorLines: true,vMinorLines: false});
        chart2.addAxis("x",{includeZero: true});
        '
        .$nextgen_graph
        .'
        chart2.setTheme(dojox.charting.themes.Julie);
        chart2.render();
        var legend2 = new dojox.charting.widget.Legend({chart: chart2, horizontal: false}, "legend2");
        ';
    }

    $javascript.='}; dojo.addOnLoad(makeCharts);';
    drupal_add_js($javascript,'inline');
    $return_graphs = '
    <div id="coverage_graphs" style="width: 500px; height: 300px; margin: 5px auto 0px auto;"></div>
  	<div id="legend1" style="margin: 5px auto 0px auto;"></div>';
    if (!empty($feature_data['nextgen_coverage_array']) ){
      $return_graphs .='<br/>
    	<div id="ngs_graphs" style="width: 500px; height: 300px; margin: 5px auto 0px auto;"></div>
  		<div id="legend2" style="margin: 5px auto 0px auto;"></div>
  		';
    }
   }
  return $return_graphs;
}

function genes4all_explore_page_feature_info_browse($feature_data){
  $return_array=array();
  // We will add JBrowse links only to specific data types (e.g. not SNP)
  $jbrowse_link ='';
  if (module_exists('jbrowse')){
    if (in_array($feature_data['type'],array('contig','ORF','polypeptide'))){
      $jbrowse_link = l(t('JBrowse page'),'jbrowse/'.$feature_data['accession'],array(
      'attributes' => array('target' => '_blank'), 'alias' => TRUE));
    }
    $return_array['jbrowse']=$jbrowse_link;
  }
  //GBrowse if present, then no jbrowse
  $gbrowse_link = array();
  if (!empty($feature_data['gbrowse'])) {
    $gbrowse_link = l(t('Gbrowse page'), $feature_data['gbrowse'],
    array(
        'attributes' => array('target' => '_blank'), 'external' => TRUE,
        'query' => array('name' => $feature_data['accession']),
    )
    );
    $return_array['gbrowse']=$gbrowse_link;
  }
  return $return_array;
}

function genes4all_explore_page_feature_info_library($feature_id){
  //Page 1 - Basic
  //Library
  drupal_add_js(drupal_get_path('module', 'gmod_dbsf') .'/includes/gmod_dbsf_sorttable.js');
  $themed_library_data  = cache_get($feature_id.'_feature_summary','cache_genes4all_explore_library');
  if (empty($themed_library_data)){
    $library_data = gmod_dbsf_db_execute('chado', 'genes4all_explore_getlibraries_details', $feature_id);
    if (empty($library_data)){
      cache_set($feature_id.'_feature_summary', 'none', 'cache_genes4all_explore_library',time()+7*3600*24);
      return '';
    }else{
      $themed_library_data = theme_table($library_data['header'], $library_data['rows'], $library_data['attributes'], $library_data['caption']);
      cache_set($feature_id.'_feature_summary', $themed_library_data, 'cache_genes4all_explore_library',time()+7*3600*24);
      return $themed_library_data;
    }
  }elseif($themed_library_data->data == 'none'){
    return '';
  }else{
    return $themed_library_data->data;
  }
}

function genes4all_explore_page_feature_info_annotation($accession=NULL){
  if (empty($accession)){return;}

  $feature_data =cache_get($accession,'cache_genes4all_explore_feature');
  if (empty($feature_data)){return;}
  $feature_data=$feature_data->data;
  //TODO limit to top 10??; phylogenet blast hits
  $annotation = '';
  //more can be : blastdbs,blast_values,blast_taxa, feature_dbxref,feature_cvterm
  $more_url = 'genes4all/info/showmore';
  if (!empty($feature_data['blastdbs'])) {
    $text = "<h3>This feature has been searched with BLAST versus the following databases</h3><ul class='item-list'>";
    $i = (int)0;
    $mored = FALSE;
    foreach ($feature_data['blastdbs'] as $blastdb) {
      $i++;
      if ($i <= 14) {
        $text .= '<li><span>'. $blastdb[0] .'</span><span class="hovertip">'.$blastdb[1].'</span></li>';
      }
      elseif ($i == 15) {
        $mored = TRUE;
        $more_link = l('more...', $more_url,
        array('html' => TRUE,
            'attributes' => array('target' => '_blank'),
            'query' => array('feature_id' => $feature_data['id'], 'info' => 'blastdbs'),
        )
        );
        $text .= $more_link;
      }
    }
    if (empty($mored)){
      $details = l('Details', $more_url,
      array('html' => TRUE,
            'attributes' => array('target' => '_blank'),
            'query' => array('feature_id' => $feature_data['id'], 'info' => 'blastdbs'),
      )
      );
      $text .= $details;
    }
    $text .= '</ul>';
    $annotation .= $text;
  }
  if (!empty($feature_data['blast_values'])) {
    $text = "<h3>These keywords have been identified from the BLAST description lines:</h3><ul class='item-list'>";
    $i = (int)0;
    foreach ($feature_data['blast_values'] as $blast) {
      $i++;
      if ($i <= 10) {
        $text .= '<li >'. $blast .'</li>';
      }
      elseif ($i == 11) {
        $more_link = l('more...', $more_url,
        array('html' => TRUE, 'attributes' => array('target' => '_blank'),
            'query' => array('feature_id' => $feature_data['id'], 'info' => 'blast_values'),
        )
        );
        $text .= $more_link;
      }
    }
    $text .= '</ul>';
    $annotation .= $text;
  }
  if (!empty($feature_data['blast_orgs'])) {
    $text = "<h3>These organisms have a protein which is significantly similar to this feature:</h3><ul class='item-list'>";
    $i = (int)0;
    foreach ($feature_data['blast_orgs'] as $blast) {
      $i++;
      if ($i <= 10) {
        //unfortunately, l() seems not to work with  li???
        $blast_org_link = l($blast['name'] .'<br>',
          'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi',
        array(
            'attributes' => array('target' => '_blank'),
            'query' => array('id' => $blast['ncbi']),
            'absolute' => TRUE, 'html' => TRUE, 'external' => TRUE,
        )
        );
        $text .= $blast_org_link;
      }
      elseif ($i == 11) {

        $more_link = l('more...', $more_url, array('html' => TRUE, 'attributes' => array('target' => '_blank'), 'query' => array('feature_id' => $feature_data['id'], 'info' => 'blast_taxa')));
        $text .= $more_link;
      }
    }
    $text .= '</ul>';
    $annotation .= $text;
  }

  if (!empty($feature_data['dbxref'])) {
    $text = "Further annotations:<br><ul class='item-list'>";
    $i = (int)0;
    foreach ($feature_data['dbxref'] as $dbxref) {
      $i++;
      if ($i <= 14) {
        $text .= '<li>';
        if (!empty($dbxref['db_name'])) {
          $text .= $dbxref['db_name'] .' found that this feature';
        }
        if (!empty($dbxref['dbxref_name'])) {
          $text .= ' has a motif belonging to the '. $dbxref['dbxref_name'];
        }

        if (!empty($dbxref['dbxref_property'])) {
          $dbxref['dbxref_property'] = strtolower($dbxref['dbxref_property']);
          $text .= ' '. $dbxref['dbxref_property'];
        }
        $text .= '</li>';
      }
      elseif ($i == 15) {
        $more_link = l('more...', $more_url, array('html' => TRUE, 'attributes' => array('target' => '_blank'), 'query' => array('feature_id' => $feature_data['id'], 'info' => 'feature_dbxref')));
        $text .= $more_link;
      }
    }
    $text .= "</ul>";
    $annotation .= $text;
  }

  if (!empty($feature_data['cvterm'])) {
    $text = "Further annotations:<br><ul class='item-list'>";
    $i = (int)0;
    foreach ($feature_data['cvterm'] as $cvterm) {
      $i++;
      if ($i <= 14) {
        $text .= '<li>';
        if (!empty($cvterm['cv_name'])) {
          $cvterm['cv_name'] = str_replace('_', ' ', $cvterm['cv_name']);
          $text .= ' from the '. $cvterm['cv_name'] .' ontology: ';
        }
        if (!empty($cvterm['cvterm_name'])) {
          $text .= 'A similarity to a known protein which has an annotation of '. $cvterm['cvterm_name'];
        }
        $text .= '</li>';
      }
      elseif ($i == 15) {
        $more_link = l('more...', $more_url, array('html' => TRUE, 'attributes' => array('target' => '_blank'), 'query' => array('feature_id' => $feature_data['id'], 'info' => 'feature_cvterm')));
        $text .= $more_link;
      }
    }
    $text .= "</ul>";
    $annotation .= $text;
  }
  if (module_exists('genes4all_curate')){
      $existing_annotation_table = genes4all_explore_page_feature_annotations_table($feature_data['uniquename']);
      $annotation .= theme_table($existing_annotation_table['header'],$existing_annotation_table['data'],$existing_annotation_table['attributes'],$existing_annotation_table['caption']);
  }
  if (!empty($annotation)){
    return "<h2>Annotation</h2><ul class='item-list'>$annotation</ul>";
    //if it was ajax, it would use echo and not return
  }
}


function genes4all_explore_page_feature_annotations_table($feature_name){
  if (empty($feature_name)) {return;}
  $feature_id = gmod_dbsf_chado_name2id($feature_name, 'feature');
  //$feat_properties  = gmod_dbsf_db_execute('chado', 'gmod_dbsf_get_feature_cvterm_and_prop', $feature_name);
  $local_annotations =genes4all_curation_get_all_curation_actions('feature',$feature_id);
  $table_rows       = array();
  $table_header     = array('Ontology', 'Term','1st submitter (date)','Agreements','Disagreements');
  if (!empty($local_annotations['cvterm'])) {
    foreach ($local_annotations['cvterm'] as $type => $data_array) {
      foreach ($data_array as $id=>$term) {
        $text =$type;
        $term_text = $term;
        if ($type ==='free_curation_terms'){$text1 = 'Free';}
        if (!empty($local_annotations['cvterm_def'][$type][$id][1])){
          $term_text= '<span>'.$term.'</span><span class="hovertip">'
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
        array('data' => $text),
        array('data' => $term_text),
        array('data' => $first_submitter.' ('.$curator_data['date'][0].')'),
        array('data' => $agreements),
        array('data' => $disagreements)
        );
      }
    }
  }if (!empty($local_annotations['dbxref'])) {
    foreach ($local_annotations['dbxref'] as $type => $data_array) {
      foreach ($data_array as $id=>$term) {
        $text =$type;
        $term_text = $term;
        if ($type ==='free_curation_terms'){$text1 = 'Free';}
        if (!empty($local_annotations['dbxref_def'][$type][$id][1])){
          $term_text= '<span>'.$term.'</span><span class="hovertip">'
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
        array('data' => $text),
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
  ,'caption'=>t('Overview of annotations from human curators')
  ,'header'=>$table_header);
}


/**
 *
 * @param $lib_accession
 *
 * @return unknown_type
 */
function genes4all_explore_page_library_info($lib_accession = NULL) {
  if (empty($lib_accession)) {
    $lib_accession = $_GET['library_id'];
    if (empty($lib_accession)) {
      return t("I'm sorry, you did not specify a library accession");
    }
  }
  $lib_accession = check_plain(trim($lib_accession));

  $array = gmod_dbsf_db_execute('chado', 'gmod_dbsf_chado_determine_library_accession', $lib_accession);
  if (empty($array)) {
    drupal_set_message(t('Cannot find any data for library %lib_accession. Are you sure it is supposed to exist? Please leave a feedback.', array('%lib_accession' => $lib_accession)), 'error');
    return '';
  }
  drupal_add_js(drupal_get_path('module', 'gmod_dbsf') .'/includes/gmod_dbsf_sorttable.js');
  $lib_id        = $array['id'];
  $lib_name      = $array['name'];
  $lib_type      = $array['type'];
  $lib_org_id    = $array['organism'];
  $lib_accession = $array['accession'];

  $return_text = '';
  $lib_table = cache_get($lib_id.'_info','cache_genes4all_explore_library');
  if (empty($lib_table)){
    $lib_table = gmod_dbsf_db_execute('chado', 'genes4all_get_lib_info', $array);
    cache_set($lib_id.'_info', $lib_table, 'cache_genes4all_explore_library',time()+7*3600*24);
  }else{
    $lib_table =$lib_table->data;
  }

  //$return_text.='theme_table(array('Tag','Value'),$rows,$attributes,$caption).'</div>';
  $graph_data =cache_get($lib_id.'_graphs','cache_genes4all_explore_library');
  if (empty($graph_data)) {
    $graph_data=gmod_dbsf_db_execute('chado', 'genes4all_explore_library_get_go', array('lib_id' => $lib_id));
    cache_set($lib_id.'_graphs', $graph_data, 'cache_genes4all_explore_library'); //cache permanently
  }else{
    $graph_data=$graph_data->data;
  }
  $graphs =genes4all_explore_library_get_go_graphs($graph_data);

  if (!empty($lib_table)) {
    $return_text .= (!empty($graphs))
    ? '<div style="align:right; position: relative;float: right;">' . call_user_func_array('theme_table',$lib_table)  .'</div>'
    : '<p>'.call_user_func_array('theme_table',$lib_table).'</p>';
  }

  if (!empty($graphs)) {
    $return_text.='<div style="align: left;position: relative;float: left;">';
    foreach ($graphs as $graph) {
      $return_text .= $graph;
    }
    $return_text.='</div>';
  }
  if (empty($return_text)) {
    $return_text .= t("I'm sorry but there is no annotation stored in the database for this library.");
  }
  return $return_text;
}

/*
 * A page where a lot of text can be displayed on its own
 */
function genes4all_explore_page_moreinfo($feature_id = NULL) {
  if (empty($feature_id)) {
    $feature_id = check_plain($_GET['feature_id']);
  }
  $info = check_plain($_GET['info']);

  if (!empty($feature_id) && !empty($info)) {
    //info can be : blastdbs,blast_values,blast_taxa, feature_dbxref,feature_cvterm
    $return_text = '';
    switch ($info) {
      case 'blastdbs':
        //two indexes for blastdbs: [0] is protein_matches and [1] is nucleotide matches
        $return_array = gmod_dbsf_db_execute('chado', 'gmod_dbsf_chado_get_feature_blastdbs', $feature_id);
        $return_text .= t('<h3>This feature had BLAST hits with the following database</h3>')
        ."Searched in protein space:<ul class='item-list'>";
        foreach ($return_array[0] as $blastdb) {
          $return_text .= '<li><span>'. $blastdb[0] .'</span><span class="hovertip">'.$blastdb[1].'</span></li>';
        }
        $return_text .= "</ul>Searched in nucleotide space:<ul class='item-list'>";
        foreach ($return_array[1] as $blastdb) {
          $return_text .= '<li><span>'. $blastdb[0] .'</span><span class="hovertip">'.$blastdb[1].'</span></li>';
        }
        break;

      case 'blast_values':
        $return_array = gmod_dbsf_db_execute('chado', 'gmod_dbsf_chado_get_feature_blast_values', $feature_id);
        $return_text = t('These keywords have been identified from the BLAST description lines') .":<br/><ul class='item-list'>";
        foreach ($return_array as $blast) {
          $return_text .= '<li>'. $blast .'</li>';
        }
        break;

      case 'blast_taxa':
        $return_array = gmod_dbsf_db_execute('chado', 'genes4all_explore_get_feature_blast_taxa', $feature_id);
        $return_text = t('These organisms have a protein which is significantly similar to this feature') .":<br/><ul class='item-list'>";
        foreach ($return_array as $blast) {
          //unfortunately, l() seems not to work with  li???
          $blast_org_link = l($blast['name'] .'<br>',
            'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi',
          array(
              'attributes' => array('target' => '_blank'),
              'query' => array('id' => $blast['ncbi']),
              'absolute' => TRUE, 'html' => TRUE,
          )
          );
          $return_text .= $blast_org_link;
        }
        break;

      case 'feature_dbxref':
        $return_array = gmod_dbsf_db_execute('chado', 'genes4all_explore_get_feature_dbxref_data', $feature_id);
        $return_text = t('Further annotations') .":<br><ul class='item-list'>";
        foreach ($return_array as $dbxref) {
          $return_text .= '<li>';
          if (!empty($dbxref['db_name'])) {
            $return_text .= $dbxref['db_name'] .' found that this feature';
          }
          if (!empty($dbxref['dbxref_name'])) {
            $return_text .= t(' has a motif belonging to the ') . $dbxref['dbxref_name'];
          }
          if (!empty($dbxref['dbxref_property'])) {
            $dbxref['dbxref_property'] = strtolower($dbxref['dbxref_property']);
            $return_text .= ' '. $dbxref['dbxref_property'];
          }
          $return_text .= '</li>';
        }
        break;

      case 'feature_cvterm':
        $return_array = gmod_dbsf_db_execute('chado', 'genes4all_explore_get_feature_cvterm_data', $feature_id);
        $return_text = t('Further annotations') .":<br><ul class='item-list'>";
        foreach ($return_array as $cvterm) {
          $return_text .= '<li>';
          if (!empty($cvterm['cv_name'])) {
            $cvterm['cv_name'] = str_replace('_', ' ', $cvterm['cv_name']);
            $return_text .= t(' from the ') . $cvterm['cv_name'] .' ontology: ';
          }
          if (!empty($cvterm['cvterm_name'])) {
            $return_text .= t('A similarity to a known protein which has an annotation of ') . $cvterm['cvterm_name'];
          }
          $return_text .= '</li>';
        }
        break;
    }
    $return_text .= '</ul>';
    if (!empty($return_text)) {
      return $return_text;
    }
    else {
      return t('No more info to show');
    }
  }
  else {
    return t('No more info to show');
  }
}


/*
 * Fetch all the information for a particular feature and return it an as assoc array
 * it can accept either the feature id or the feature name
 * should this be for more than one feature id ?
 */
function genes4all_explore_feature_data($accession) {
  if (empty($accession)) {
    return FALSE;
  }

  $array = genes4all_explore_determine_feature_accession($accession);
  if (empty($array)) {
    return FALSE;
  }
  $feature_id = $array['id'];
  $feature_name = $array['name'];
  unset($array);
  // Chado specific function so Chado SQLs
  $feature_info_sql = 'SELECT feature.organism_id,genus,species,abbreviation,common_name,' . 'feature.name,feature.uniquename,residues,cvterm.name as type_name,feature.is_obsolete,timelastmodified ' . ' from feature JOIN organism ON feature.organism_id=organism.organism_id JOIN cvterm ON feature.type_id=cvterm.cvterm_id' . ' where feature_id='. $feature_id;
  $feature_data_auxilary = db_fetch_array(db_query($feature_info_sql));
  if (empty($feature_data_auxilary['uniquename'])) {
    return FALSE;
  }

  //NEXT gen; currently only for 'contig'
  $nextgen_data = array();


  $organism_sql = 'SELECT dbxref.accession from dbxref where dbxref_id = ' . ' (select dbxref_id from organism_dbxref where organism_id ='. $feature_data_auxilary['organism_id'] .')';
  $org_result   = db_query($organism_sql);
  $db_org_data  = db_fetch_array($org_result);

  //TODO: for the following 2, i made separate functions which should replace the calls below:
  //genes4all_explore_get_feature_dbxref_data,genes4all_explore_get_feature_cvterm_data

  $dbxref_data   = array();
  $dbxref_sql    = 'SELECT dbxref_name,db_name,dbxref_property from dbxref_organism_mat where feature_id='. $feature_id .' LIMIT 15';
  $dbxref_result = db_query($dbxref_sql);
  while ($db_dbxref_data = db_fetch_array($dbxref_result)) {
    $dbxref_data[] = $db_dbxref_data;
  }

  $cvterm_data   = array();
  $cvterm_sql    = 'SELECT cvterm_name,cv_name from cvterm_organism_mat where feature_id='. $feature_id .' LIMIT 15';
  $cvterm_result = db_query($cvterm_sql);
  while ($db_cvterm_data = db_fetch_array($cvterm_result)) {
    $cvterm_data[] = $db_cvterm_data;
  }

  //Links to contig, polypeptide, orf.
  $link_array = array();
  //prot to contig via ORF ref
  $sql_prot2contig = 'SELECT feature_id,uniquename from feature where type_id IN ' . " (select cvterm_id from cvterm where name ~ 'contig|mRNA' and cv_id=(select cv_id from cv where name='sequence'))" . ' and feature_id IN (select srcfeature_id from featureloc where feature_id = ' . " (select feature_id from feature where uniquename ='";
  //contig to orf reference/placeholder, we have the feature_id
  // to get peptide strip _0 replace /Aorf/Apep/
  //don't forget the closing bracket
  $sql_contig2orf_ref = "SELECT uniquename from feature where type_id =(select cvterm_id from cvterm where name = 'ORF'" . " AND cv_id=(select cv_id from cv where name='sequence'))"
  . " AND feature_id IN (select feature_id from featureloc where srcfeature_id =";

  //number of ESTs
  $est_count = '';
  $est_count_sql = 'select count(feature_id) as est_count from feature,cvterm where '
  . " feature.type_id=cvterm.cvterm_id and cvterm.name='read'"
  . ' and feature_id IN (select feature_id from featureloc where srcfeature_id=';

  //prediction method
  $prediction_method = '';
  $sql_prot_pred_method = 'select name as prediction_method from analysis where analysis_id =(select analysis_id from '
  . ' analysisfeature where feature_id=';

  // if a EST get parent via 'read'
  $sql_est2contig = "SELECT srcfeature_id as contig_feature_id,feature.uniquename from featureloc join feature on featureloc.srcfeature_id=feature.feature_id where " . " feature.type_id=(select cvterm_id from cvterm join cv on cvterm.cv_id=cv.cv_id where cvterm.name='contig' and cv.name='sequence')" . " AND featureloc.feature_id IN (select featureloc.feature_id from featureloc join feature on featureloc.feature_id=feature.feature_id where ". " srcfeature_id= $feature_id AND feature.type_id=("
  . "select cvterm_id from cvterm join cv on cvterm.cv_id=cv.cv_id where cvterm.name='read' and cv.name='sequence'" . "))";

  //if a SNP get parent via
  $sql_snp2contig = "SELECT srcfeature_id as contig_feature_id,feature.uniquename from featureloc join feature on featureloc.srcfeature_id=feature.feature_id where " . " feature.type_id=(select cvterm_id from cvterm join cv on cvterm.cv_id=cv.cv_id where cvterm.name='contig' and cv.name='sequence')"
  . " AND featureloc.feature_id= $feature_id";
  switch ($feature_data_auxilary['type_name']) {
    case 'contig':
      $nextgen_data = genes4all_nextgen_coverage_data($feature_name,TRUE);
      // Allow more than one protein prediction
      $dbres = db_query($sql_contig2orf_ref . $feature_id .')');
      while ($res = db_fetch_array($dbres)) {
        $orf_ref = $res['uniquename'];
        $orf_id = str_replace("_0", "", $orf_ref);
        if (strpos($orf_id, 'Corf')) {
          $pept_id = str_replace("Corf", "Cpep", $orf_id);
        }
        else {
          $pept_id = str_replace("Aorf", "Apep", $orf_id);
        }
        $pept_id = gmod_dbsf_get_feature_name_id($pept_id, 'STRING');
        $orf_id = gmod_dbsf_get_feature_name_id($orf_id, 'STRING');
        $link_array['ORF'][] = $orf_id;
        $link_array['polypeptide'][] = $pept_id;
      }
      $res = db_fetch_array(db_query($est_count_sql . $feature_id .')'));
      $est_count = $res['est_count'];
      break;

    case 'mRNA':
      $dbres = db_query($sql_contig2orf_ref . $feature_id .')');
      while ($res = db_fetch_array($dbres)) {
        $orf_ref = $res['uniquename'];
        $orf_id = str_replace("_0", "", $orf_ref);
        if (strpos($orf_id, 'Corf')) {
          $pept_id = str_replace("Corf", "Cpep", $orf_id);
        }
        else {
          $pept_id = str_replace("Aorf", "Apep", $orf_id);
        }
        $pept_id = gmod_dbsf_get_feature_name_id($pept_id, 'STRING');
        $orf_id = gmod_dbsf_get_feature_name_id($orf_id, 'STRING');
        $link_array['ORF'][] = $orf_id;
        $link_array['polypeptide'][] = $pept_id;
      }
      $res = db_fetch_array(db_query($est_count_sql . $feature_id .')'));
      $est_count = $res['est_count'];
      break;

    case 'polypeptide':
      if (strpos($feature_name, 'Cpep')) {
        $orf_id = str_replace("Cpep", "Corf", $feature_name);
      }
      else {
        $orf_id = str_replace("Apep", "Aorf", $feature_name);
      }
      $res = db_fetch_array(db_query($sql_prot2contig . $orf_id ."_0'))"));
      $contig_id = gmod_dbsf_get_feature_name_id($res['uniquename'], 'STRING');
      $orf_id = gmod_dbsf_get_feature_name_id($orf_id, 'STRING');
      $link_array['contig'] = $contig_id;
      $link_array['ORF'] = $orf_id;
      $res2 = db_fetch_array(db_query($sql_prot_pred_method . $feature_id .')'));
      $prediction_method = $res2['prediction_method'];
      break;

    case 'ORF':
      if (strpos($feature_name, 'Corf')) {
        $pept_id = str_replace("Corf", "Cpep", $feature_name);
      }
      else {
        $pept_id = str_replace("Aorf", "Apep", $feature_name);
      }
      $pept_id = gmod_dbsf_get_feature_name_id($pept_id, 'STRING');
      $res = db_fetch_array(db_query($sql_prot2contig . $feature_name ."_0'))"));
      $contig_id = gmod_dbsf_get_feature_name_id($res['uniquename'], 'STRING');
      $link_array['contig'] = $contig_id;
      $link_array['polypeptide'] = $pept_id;
      break;

    case 'EST':
      $res = db_fetch_array(db_query($sql_est2contig));
      $contig_id = gmod_dbsf_get_feature_name_id($res['uniquename'], 'STRING');
      $link_array['contig'] = $contig_id;
      $contig_feature_id = $res['contig_feature_id'];
      unset($res);
      $dbres = db_query($sql_contig2orf_ref . $contig_feature_id .')');
      while ($res = db_fetch_array($dbres)) {
        $orf_ref = $res['uniquename'];
        $orf_id = str_replace("_0", "", $orf_ref);
        if (strpos($orf_id, 'Corf')) {
          $pept_id = str_replace("Corf", "Cpep", $orf_id);
        }
        else {
          $pept_id = str_replace("Aorf", "Apep", $orf_id);
        }
        $pept_id = gmod_dbsf_get_feature_name_id($pept_id, 'STRING');
        $orf_id = gmod_dbsf_get_feature_name_id($orf_id, 'STRING');
        $link_array['ORF'][] = $orf_id;
        $link_array['polypeptide'][] = $pept_id;
      }
      break;

    case 'SNP':
      $res = db_fetch_array(db_query($sql_snp2contig));
      $contig_id = gmod_dbsf_get_feature_name_id($res['uniquename'], 'STRING');
      $link_array['contig'] = $contig_id;
      $contig_feature_id = $res['contig_feature_id'];
      unset($res);
      $dbres = db_query($sql_contig2orf_ref . $contig_feature_id .')');
      while ($res = db_fetch_array($dbres)) {
        $orf_ref = $res['uniquename'];
        $orf_id = str_replace("_0", "", $orf_ref);
        if (strpos($orf_id, 'Corf')) {
          $pept_id = str_replace("Corf", "Cpep", $orf_id);
        }
        else {
          $pept_id = str_replace("Aorf", "Apep", $orf_id);
        }
        $pept_id = gmod_dbsf_get_feature_name_id($pept_id, 'STRING');
        $orf_id = gmod_dbsf_get_feature_name_id($orf_id, 'STRING');
        $link_array['ORF'][] = $orf_id;
        $link_array['polypeptide'][] = $pept_id;
      }
      break;
  }

  // BLASTS
  //databases blasted against. This is an array with two indexes: [0] is protein_matches and [1] is nucleotides
  $blastdb_against_array = gmod_dbsf_chado_get_feature_blastdbs($feature_id, '15');
  //blast description array
  $blastdb_values_array = gmod_dbsf_chado_get_feature_blast_values($feature_id, '11');
  //blast organisms
  $blastdb_orgs_array = genes4all_explore_get_feature_blast_taxa($feature_id, '11');
  // if multiple features are being searched, then the following array will have to increase by one level?
  $feature_data = array(
    'blast_orgs' => $blastdb_orgs_array,
    'blast_values' => $blastdb_values_array,
    'blastdbs' => $blastdb_against_array[0]+$blastdb_against_array[1],
    'prediction_method' => $prediction_method,
    'est_count' => $est_count,
    'links' => $link_array,
    'dbxref' => $dbxref_data,
    'cvterm' => $cvterm_data,
    'org_id' => $feature_data_auxilary['organism_id'],
    'org_taxid' => $db_org_data['accession'],
    'id' => $feature_id,
    'accession' => $feature_name,
    'type' => $feature_data_auxilary['type_name'],
    'uniquename' => $feature_data_auxilary['uniquename'],
    'lastmodified' => $feature_data_auxilary['timelastmodified'],
    'genus' => $feature_data_auxilary['genus'],
    'species' => $feature_data_auxilary['species'],
    'species_abbr' => $feature_data_auxilary['abbreviation'],
    'species_common' => $feature_data_auxilary['common_name'],
  );
  /* decided to remove the assembly option
   $feature_data['assembly'] = gmod_dbsf_chado_assembly2gff($feature_name);
   */
  if (!empty($feature_data_auxilary['residues'])) {
    $feature_data_auxilary['residues'] = str_replace('-', '', $feature_data_auxilary['residues']);
    $feature_data['seqlength'] = strlen($feature_data_auxilary['residues']);
    $feature_data['residues'] = $feature_data_auxilary['residues'];
    if ($feature_data_auxilary['type_name'] == 'polypeptide') {
      $feature_data['seqlength_text'] = $feature_data['seqlength'] .' amino acids';
    }
    else{
      $feature_data['seqlength_text'] = $feature_data['seqlength'] .' base pairs';
    }
  }

  if ($feature_data_auxilary['type_name'] == 'SNP') {
    $feature_data['snp_info'] = genes4all_explore_feature_snp_data($feature_id);
  }

  if ($feature_data_auxilary['is_obsolete'] == 't') {
    $feature_data['obsolete'] = 'true';
  }
  //GRAPHS
  if (!empty($nextgen_data)){
    $feature_data['nextgen_coverage_array'] = $nextgen_data[$feature_name];
  }
  if (!empty($feature_data_auxilary['residues'])) {
    $res = genes4all_explore_plot_blast_coverage($feature_data_auxilary['type_name'],$feature_id,$feature_data['seqlength'],$feature_name);
    $feature_data['blast_coverage_array']=$res[0];
    $feature_data['blast_coverage_array_max'] = $res[1];
    unset($res);
    $res = genes4all_explore_plot_coverage($feature_data_auxilary['type_name'],$feature_name,$feature_data['seqlength'],$feature_id);
    $feature_data['coverage_array'] = $res[0];
    $feature_data['coverage_array_max'] = $res[1];
    $feature_data['lib_coverage_array'] = $res[2];
    unset($res);
  }
  //GBROWSE LINK
  if (!module_exists('jbrowse')){
    $feature_data['gbrowse'] =  genes4all_explore_gbrowse_link($feature_data_auxilary['type_name'],$db_org_data['accession']);
  }
  return $feature_data;
}

/**
 *
 */
function genes4all_explore_gbrowse_link($type_name=NULL,$db_org_accession=NULL){
  if (empty($type_name) || empty($db_org_accession)){return;}
  $data_object = '';
  $gbrowse_url = gmod_dbsf_get_add_var('genes4all_explore_gbrowse_link');
  if ($type_name == 'ORF') {
    $data_object = '_orf';
  }
  elseif ($type_name == 'polypeptide') {
    $data_object = '_prot';
  }
  elseif ($type_name == 'contig' || $feature_data_auxilary['type_name'] == 'mRNA') {
    $data_object = '_caf';
  }
  elseif ($type_name == 'EST') {
    $data_object = '_caf';
  }
  elseif ($type_name == 'SNP') {
    $data_object = '_caf';
  }
  if (!empty($data_object)) {
    $gbrowse_url .=  $db_org_accession . $data_object .'/';
  }
  return $gbrowse_url;
}

/**
 *
 */
function genes4all_explore_plot_blast_coverage($feature_type=NULL,$feature_id=NULL,$seqlength=NULL,$feature_name=NULL){
  if (empty($feature_type)||empty($feature_id)||empty($seqlength)){return;}
  $blast_coverage_array = array();
  $graph_window =15;
  $blast_coverage_sql='';
  if ($feature_type == 'contig' || $feature_type == 'ORF' ) {
    //decrease this to make it smoother.
    $graph_window = 20;
    if ($seqlength > 2000) {
      $graph_window = 5 * $graph_window;
    }
    elseif ($seqlength > 1000) {
      $graph_window = 2.5 * $graph_window;
    }
  }
  // for ORF we use the polypeptide BLAST results ;-)
  if ($feature_type == 'ORF'){
    if (empty($feature_name)){return;}
    $feature_name=str_replace('Aorf','Apep',$feature_name);
    $seqlength/=3;
    $blast_coverage_sql = "SELECT count(fl.feature_id) FROM featureloc fl JOIN feature f ON f.feature_id=fl.feature_id join cvterm ON f.type_id=cvterm.cvterm_id WHERE fl.srcfeature_id=(SELECT feature_id from feature where uniquename='$feature_name') AND (cvterm.name='nucleotide_match' OR cvterm.name='protein_match') AND fl.fmin<=%d and fl.fmax>=%d";
  }else{
    $blast_coverage_sql = "SELECT count(fl.feature_id) FROM featureloc fl JOIN feature f ON f.feature_id=fl.feature_id join cvterm ON f.type_id=cvterm.cvterm_id WHERE fl.srcfeature_id='$feature_id' AND (cvterm.name='nucleotide_match' OR cvterm.name='protein_match') AND fl.fmin<=%d and fl.fmax>=%d";
  }
  $max = (int)0;
  for ($i = 0; $i < $seqlength; $i += $graph_window) {
    $result = db_query($blast_coverage_sql, $i + $graph_window, $i);
    $db_result = db_fetch_array($result);

    if ($feature_type == 'ORF'){
      $blast_coverage_array[] = array('y' => $db_result['count'], 'x' => $i*3);
    }else{
      $blast_coverage_array[] = array('y' => $db_result['count'], 'x' => $i);
    }
    if ($db_result['count']>$max){
      $max = $db_result['count'];
    }
  }

  if (!empty($blast_coverage_array) && $max!==0) {
    return array($blast_coverage_array,$max);
  }else{
    return FALSE;
  }
}

/**
 * Create data for read coverage
 */
function genes4all_explore_plot_coverage($feature_type=NULL,$feature_name=NULL,$seqlength=NULL,$feature_id=NULL){
  if (empty($feature_type)||empty($feature_name)||empty($seqlength)||empty($feature_id)){return;}
  $coverage_array = array();
  $lib_coverage_array = array(); // library specific coverage
  if ($feature_type == 'contig') {
    $graph_window = 20;
    if ($seqlength > 2000) {
      $graph_window = 5 * $graph_window;
    }
    elseif ($seqlength > 1000) {
      $graph_window = 2.5 * $graph_window;
    }
    $max = (int)0;
    //$coverage_sql = "SELECT count(fl.feature_id) FROM featureloc fl JOIN feature f ON f.feature_id=fl.feature_id join cvterm ON f.type_id=cvterm.cvterm_id WHERE fl.srcfeature_id=(select feature_id from feature where uniquename='$feature_name') AND cvterm.name='read' AND fl.fmin<=%d and fl.fmax>=%d";
    //gff3view is FASTER. Count removed to facilitate per-library queries
    $coverage_sql = "SELECT feature_id FROM gff3view_mat WHERE ref='$feature_name' and type='read' and fstart<=%d and fend>=%d";
    // the library specific coverage is a bit more difficult: in the GFF the reads are present but for library we use the EST
    // the relationship between EST and read is very simple (read is EST + _assembly) especially with feature_id:
    // the EST is feature.id(read)+1 so read is feature_id(EST)-1. From there, we can get the library_id from library_feature.
    // then we need to make a coverage graph for each library.name
    // The most efficient way to get the coverage is from the gff3view. For each query, how many of them belong to which library?
    $feature2library_sql = "SELECT library.name as libname from {library} JOIN library_feature on library.library_id = library_feature.library_id"
    ." WHERE library_feature.feature_id=%d-1";
    $k = (int)0;
    for ($i = 0; $i < $seqlength; $i += $graph_window) {
      $result       = db_query($coverage_sql, $i + $graph_window, $i);
      $count        = pg_num_rows($result);
      $coverage_array[] = array('y' => $count, 'x' => $i);
      while ($db_result   = db_fetch_array($result)){
        $lib_result       = db_query($feature2library_sql,$db_result['feature_id']);
        $lib_db_result    = db_fetch_array($lib_result);
        if (empty($lib_db_result['libname'])){
          continue;
        }
        if (!isset($lib_coverage_array[$lib_db_result['libname']])){
          $lib_coverage_array[$lib_db_result['libname']]=array();
        }
        if (!isset($lib_coverage_array[$lib_db_result['libname']][$k])){
          $lib_coverage_array[$lib_db_result['libname']][$k] = array('y'=>0,'x'=>$i);
        }
        $lib_coverage_array[$lib_db_result['libname']][$k]['y']++;
      }
      if ($count > $max){
        $max = $count;
      }
      $k++;
    }
  }
  // ensure lib_coverage array is fully populated for each library
  // also annoyingly, the graph needs an indexed array...
  if (!empty($lib_coverage_array)){
    foreach ($lib_coverage_array as $libname=>$data){
      $k = (int)0;
      for ($i = 0; $i < $seqlength; $i += $graph_window) {
        if (!isset( $lib_coverage_array[$libname][$k])){
          $lib_coverage_array[$libname][$k] = array('y'=>0,'x'=>$i);
        }
        $k++;
      }
      ksort($lib_coverage_array[$libname]);
    }
  }
  if (!empty($coverage_array)) {
    return array($coverage_array,$max,$lib_coverage_array);
  }
  else{
    return FALSE;
  }
}

/**
 *
 */
function genes4all_explore_feature_snp_data($feature_id=NULL){
  if (!empty($feauture_id)){return;}
  $data = gmod_dbsf_get_feature_cvterm_and_prop($feature_id);
  $headers=array();
  $columns =array();
  foreach ($data['prop'] as $key=>$value){
    switch ($key){
      case 'orf ID':
        $columns[]=l($value[0],'genes4all/info/feature',array(
          'alias'=>TRUE,'query'=>array('feature_id'=>$value[0])));
        $headers[]='Belongs to ORF';
        break;
      case 'codon_change':
        $columns[]=$value[0];
        $headers[]='Affects codon';
        break;
      case 'substitution_type':
        $columns[]=$value[0];
        $headers[]='Type of substitution';
        break;
      case 'minor_frequency':
        $columns[]=$value[0];
        $headers[]='Minor allele frequency';
        break;
      case 'alleles':
        $columns[]=implode('/',$value);
        $headers[]='Alleles';
        break;
      case 'total_reads_position':
        $columns[]=$value[0];
        $headers[]='Reads supporting position';
        break;
      case 'position in codon':
        $columns[]=$value[0];
        $headers[]='Codon site';
        break;
      case 'quality_SNP':
        $columns[]=$value[0];
        $headers[]='Consensus base quality';
        break;
    }
  }
  if (!empty($columns)) {
    return theme_table($headers,array(0=>$columns),array('class' => 'genes4all_explore-library-table sortable'),'<h2>Marker Information</h2>');;
  }
}


/**
 *
 * @param $feature_id
 * @param $limit
 *
 * @return unknown_type
 */
function genes4all_explore_get_feature_blast_taxa($feature_id, $limit = NULL) {
  if (empty($feature_id)) {
    return;
  }
  $blastdb_orgs_array = array();
  //Chado-specific SQLs so no {} fields
  $blastdb_orgs_sql = 'SELECT accession,description from feature_dbxref join dbxref ON dbxref.dbxref_id=feature_dbxref.dbxref_id where '. " db_id =(select db_id from db where name='NCBI_TAXONOMY') and feature_id= $feature_id";
  if (!empty($limit)) {
    $blastdb_orgs_sql .= " LIMIT $limit ";
  }
  $result = db_query($blastdb_orgs_sql);
  while ($db_blastdb = db_fetch_array($result)) {
    $blastdb_orgs_array[] = array('ncbi' => $db_blastdb['accession'], 'name' => $db_blastdb['description']);
  }
  return ($blastdb_orgs_array);
}

/**
 *
 * @param $feature_id
 * @param $limit
 *
 * @return unknown_type
 */
function genes4all_explore_get_feature_dbxref_data($feature_id, $limit = NULL) {
  if (empty($feature_id)) {
    return;
  }
  $dbxref_data = array();
  //Chado-specific SQLs so no {} fields
  $dbxref_sql = 'SELECT dbxref_name,db_name,dbxref_property from dbxref_organism_mat where feature_id='. $feature_id;
  if (!empty($limit)) {
    $dbxref_sql .= " LIMIT $limit ";
  }
  $dbxref_result = db_query($dbxref_sql);
  while ($db_dbxref_data = db_fetch_array($dbxref_result)) {
    $dbxref_data[] = $db_dbxref_data;
  }
  return ($dbxref_data);
}

/**
 *
 * @param $feature_id
 * @param $limit
 *
 * @return unknown_type
 */
function genes4all_explore_get_feature_cvterm_data($feature_id, $limit = NULL) {
  if (empty($feature_id)) {
    return;
  }
  $cvterm_data = array();
  //Chado-specific SQLs so no {} fields
  $cvterm_sql = 'SELECT cvterm_name,cv_name from cvterm_organism_mat where feature_id='. $feature_id;
  if (!empty($limit)) {
    $cvterm_sql .= " LIMIT $limit ";
  }
  $cvterm_result = db_query($cvterm_sql);
  while ($db_cvterm_data = db_fetch_array($cvterm_result)) {
    $cvterm_data[] = $db_cvterm_data;
  }
  return ($cvterm_data);
}

/*
 * Get library details
 *
 * @param $accession
 *  Feature accession that belongs to library
 *  is either an integer (internal chado id) or varchar (IC accession)
 * @param $organism_id_array
 *  is one or more organisms which in lack of feature_id can be used to
 * generate a library tablee for that species
 * @param $library_id_array
 *  likewise as above but for libraries
 *
 * @return data for table.
 */
function genes4all_explore_getlibraries_details($accession = NULL, $organism_id_array = NULL, $library_id_array = NULL) {
  // contains basic lib data
  $lib_array = array();
  // once above is populated, get additional information. get tags to serve as table headers in all libraries
  $extra_data = array();
  // populates extra_data
  $get_tags_sql = "SELECT cv.name as tag,cvterm.name as value from library_cvterm join "
  . " cvterm on library_cvterm.cvterm_id=cvterm.cvterm_id"
  . " join cv on cvterm.cv_id=cv.cv_id where "
  . " library_id =%d order by cv.name";
  //if $accession is an array, then it has been called via db_execute
  if (!empty($accession) && is_array($accession)) {
    $library_id_array  = (array)$accession[2];
    $organism_id_array = (array)$accession[1];
    $accession         = $accession[0];
  }

  //if nothing is requested, then we return a table of all libraries
  if (empty($accession) && empty($organism_id_array) && empty($library_id_array)) {

    $result = db_query('select library_id,organism_id,name from library order by name');
    while ($dbdata = db_fetch_array($result)) {
      $friendly = $dbdata['name'];
      $access = gmod_dbsf_chadoapi_check_restriction('select', 'library', $dbdata['library_id'], FALSE);
      if ($access == TRUE) {
        $lib_array[$dbdata['library_id']] = array('org_id' => $dbdata['organism_id'], 'name' => $dbdata['name']);
        $extra_info_result = db_query($get_tags_sql, $dbdata['library_id']);
        while ($extra = db_fetch_array($extra_info_result)) {
          $extra_data[$dbdata['library_id']][$extra['tag']] = $extra['value'];
        }
      }
    }
  }


  //TODO for each case create a new type of table
  elseif (!empty($accession)) {
    $array = genes4all_explore_determine_feature_accession($accession);
    if (empty($array)) {
      return;
    }
    $feature_id   = $array['id'];
    $feature_name = $array['name'];
    $feature_type = $array['type'];
    unset($array);

    //we have a single feature we want to get libraries for.
    // this feature can be an est, in which case the library_feature table holds the link
    // or it can be a contig which contains multiple ests in which case we get first the ests
    // or it can be a peptide or orf in which case we must first get the contig!

    //First get all library ids
    switch ($feature_type) {
      case 'EST':
        $lib_sql = "SELECT library_id,organism_id,name from library where library_id IN ("
        . " SELECT library_id from library_feature where feature_id= $feature_id)"
        . " order by name";
        $lib_result = db_query($lib_sql);
        // only one library id per feature_id.
        $basic = db_fetch_array($lib_result);
        $lib_array[$basic['library_id']] = array(
          'org_id' => $basic['organism_id'],
          'name' => $basic['name'],
        );
        break;

      case 'contig':
        $lib_sql = "SELECT library_id,organism_id,name from library where library_id IN ("
        . "select library_id from library_feature where "
        . " feature_id IN (select srcfeature_id from featureloc join feature on featureloc.srcfeature_id=feature.feature_id where "
        . " feature.type_id=(select cvterm_id from cvterm join cv on cvterm.cv_id=cv.cv_id where "
        . " cvterm.name='EST' and cv.name='sequence') AND featureloc.feature_id IN ("
        . " select featureloc.feature_id from featureloc join feature on featureloc.feature_id=feature.feature_id where "
        . " srcfeature_id= $feature_id AND feature.type_id=("
        . " select cvterm_id from cvterm join cv on cvterm.cv_id=cv.cv_id where "
        . " cvterm.name='read' and cv.name='sequence'))))"
        . " order by name";
        $lib_result = db_query($lib_sql);
        while ($basic = db_fetch_array($lib_result)) {
          $lib_array[$basic['library_id']] = array(
            'org_id' => $basic['organism_id'],
            'name' => $basic['name'],
          );
        }
        break;

      case 'ORF':
        // seems remo's implementation is not good enough to database orf with contig so we hack it using the name
        // add a _0 to get the reference and grab the srcfeature
        $lib_sql = "SELECT library_id,organism_id,name from library where library_id IN (" . "SELECT library_id from library_feature where feature_id IN (" . " SELECT srcfeature_id from featureloc join feature on featureloc.srcfeature_id=feature.feature_id where " . " feature.type_id=(select cvterm_id from cvterm join cv on cvterm.cv_id=cv.cv_id where cvterm.name='EST' and cv.name='sequence') AND featureloc.feature_id IN (" . " select featureloc.feature_id from featureloc join feature on featureloc.feature_id=feature.feature_id where" . " srcfeature_id=(select distinct srcfeature_id from featureloc join feature on featureloc.srcfeature_id=feature.feature_id where featureloc.feature_id=(". " select feature_id from feature where name='$feature_name'||'_0') and feature.type_id=(" . " select cvterm_id from cvterm join cv on cvterm.cv_id=cv.cv_id where cvterm.name='contig' and cv.name='sequence')) and feature.type_id=(" . " select cvterm_id from cvterm join cv on cvterm.cv_id=cv.cv_id where cvterm.name='read' and cv.name='sequence'))))" . " order by name";
        $lib_result = db_query($lib_sql);
        while ($basic = db_fetch_array($lib_result)) {
          $lib_array[$basic['library_id']] = array(
            'org_id' => $basic['organism_id'],
            'name' => $basic['name'],
          );
        }
        break;

      case 'polypeptide':
        // another hack, we will directly get the Aorf_0 by substituting Apep
        $feature_orf = str_replace('Apep', 'Aorf', $feature_name);
        $lib_sql     = "SELECT library_id,organism_id,name from library where library_id IN (" . "SELECT library_id from library_feature where feature_id IN (" . " SELECT srcfeature_id from featureloc join feature on featureloc.srcfeature_id=feature.feature_id where " . " feature.type_id=(select cvterm_id from cvterm join cv on cvterm.cv_id=cv.cv_id where cvterm.name='EST' and cv.name='sequence') AND featureloc.feature_id IN (" . " select featureloc.feature_id from featureloc join feature on featureloc.feature_id=feature.feature_id where" . " srcfeature_id=(select distinct srcfeature_id from featureloc join feature on featureloc.srcfeature_id=feature.feature_id where featureloc.feature_id=(". " select feature_id from feature where name='$feature_orf'||'_0') and feature.type_id=(" . " select cvterm_id from cvterm join cv on cvterm.cv_id=cv.cv_id where cvterm.name='contig' and cv.name='sequence')) and feature.type_id=(" . " select cvterm_id from cvterm join cv on cvterm.cv_id=cv.cv_id where cvterm.name='read' and cv.name='sequence'))))" . " order by name";
        $lib_result  = db_query($lib_sql);
        while ($basic = db_fetch_array($lib_result)) {
          $lib_array[$basic['library_id']] = array(
            'org_id' => $basic['organism_id'],
            'name' => $basic['name'],
          );
        }
        break;
    }
    foreach ($lib_array as $lib_id => $basic) {
      $friendly = $basic['name'];
      $extra_info_result = db_query($get_tags_sql, $lib_id);
      while ($extra = db_fetch_array($extra_info_result)) {
        $extra_data[$lib_id][$extra['tag']] = $extra['value'];
      }
    }
  }
  elseif (!empty($library_id_array)) {
    $lib_str = '';
    foreach ($library_id_array as $lib) {
      $lib_str .= $lib .',';
    }
    $lib_str = rtrim($lib_str, ',');
    $result = db_query(
      "SELECT library_id,organism_id,name from library ". " where library_id IN ($lib_str)" . "order by name"
      );
      while ($dbdata = db_fetch_array($result)) {
        $friendly = $dbdata['name'];

        // get tags to serve as table headers in all libraries
        // if it is not empty, then there are restrictions
        $access = gmod_dbsf_chadoapi_check_restriction('select', 'library', $dbdata['library_id'], FALSE);
        if ($access == TRUE) {
          $extra_info_result = db_query($get_tags_sql, $dbdata['library_id']);
          $lib_array[$dbdata['library_id']] = array(
          'org_id' => $dbdata['organism_id'], 'name' => $dbdata['name'],
          );
          while ($extra = db_fetch_array($extra_info_result)) {
            $extra_data[$dbdata['library_id']][$extra['tag']] = $extra['value'];
          }
        }
      }
  }
  elseif (!empty($organism_id_array)) {
    $org_str = '';
    foreach ($organism_id_array as $org) {
      $org_str .= $org .',';
    }
    $org_str = rtrim($org_str, ',');
    $result = db_query(
      "SELECT library_id,organism_id,name from library ". " where organism_id IN ($org_str)" . "order by name"
      );
      while ($dbdata = db_fetch_array($result)) {
        $friendly = $dbdata['name'];

        // get tags to serve as table headers in all libraries
        $access = gmod_dbsf_chadoapi_check_restriction('select', 'library', $dbdata['library_id'], FALSE);
        if ($access == TRUE) {
          $extra_info_result = db_query($get_tags_sql, $dbdata['library_id']);
          $lib_array[$dbdata['library_id']] = array('org_id' => $dbdata['organism_id'], 'name' => $dbdata['name']);
          while ($extra = db_fetch_array($extra_info_result)) {
            $extra_data[$dbdata['library_id']][$extra['tag']] = $extra['value'];
          }
        }
      }
  }
  if (empty($lib_array)) {
    return;
  }
  $headers_unfriendly = array('name', 'dev_stage', 'sex', 'tissue_type');
  $headers = array('Library name', 'Developmental stage', 'Sex', 'Tissue');

  $rows = array();
  foreach ($lib_array as $lib_id => $library_data) {
    foreach ($headers_unfriendly as $tag) {
      if ($tag == 'name') {
        $name = l($library_data['name'], 'genes4all/info/library',
        //alias needed
        array('alias' => TRUE, 'query' => array('library_id' => $lib_id))
        );
        $data = 'no additional data';
        if (!empty($extra_data)) {
          foreach ($extra_data[$lib_id] as $tag1 => $value) {
            if ($tag1 != 'dev_stage' && $tag1 != 'sex' && $tag1 != 'tissue_type' && $tag1 != 'name') {
              $data .= $tag1 .' => '. $value .', ';
            }
          }
          $data = rtrim($data, ', ');
        }
        $name = '<span>'. $name .'</span>' . '<span class="hovertip">'. $data .'</span>';
        $rows[$lib_id][] = $name;
      }
      elseif (!empty($library_data[$tag])) {
        $rows[$lib_id][] = $library_data[$tag];
      }
      elseif (!empty($extra_data[$lib_id][$tag])) {
        $rows[$lib_id][] = $extra_data[$lib_id][$tag];
      }
      else {
        $rows[$lib_id][] = 'unknown';
      }
    }
  }
  return array(
    'header' => $headers,
    'rows' => $rows,
    'attributes' => array('class' => 'genes4all_explore-library-table sortable'),
    'caption' => '<h2>Library information</h2>',
  );
}

/**
 *
 * @param $lib_id array parsed later by gmod_dbsf_array2csv
 * @param $limit integer limit graph data to top $limit.
 * @param in future implementation of the y-axis scaling function
 *
 * @return array of graphs
 */
function genes4all_explore_library_get_go($lib_id, $limit = 100, $base = 'log10') {
  if (empty($lib_id)) {
    return;
  }
  $lib_id = gmod_dbsf_array2csv($lib_id);
  // counts of each cvterm for this library,build through SQL. associated info attached
  $cvterms = array();
  // counts of each dbxref for this library,build through SQL. associated info attached
  $dbxrefs = array();
  $snps_array = array();
  // returned
  $graphs = array();
  // get ests from library_feature -> get contigs -> get peptide -> get cvs
  //-> get cvterms from feature_cvterm for GO EC KEGG and feature_dbxref for iprscan
  // Chado-specific SQLs so no {} fields
  $get_cvterm_sql = "select cv.name as tag,cvterm.name as name, count(cvterm.name) as value_count, dbxref.accession from " . " library_feature_cvterm_mat join cvterm on library_feature_cvterm_mat.cvterm_id=cvterm.cvterm_id join cv on cvterm.cv_id=cv.cv_id " . " join dbxref on dbxref.dbxref_id=cvterm.dbxref_id ". " where library_id IN ($lib_id) group by cv.name,cvterm.name,dbxref.accession order by value_count DESC";
  $get_dbxref_sql = "select db.name as tag,dbxref.accession,dbxref.description, count(dbxref.accession) as value_count from " . " library_feature_dbxref_mat join dbxref on library_feature_dbxref_mat.dbxref_id=dbxref.dbxref_id join db on dbxref.db_id=db.db_id ". " where library_id IN ($lib_id) group by db.name,dbxref.accession,dbxref.description order by value_count DESC";
  $get_snp_sql    = "SELECT cvterm.name as tag,value as name,count(value) as value_count from library_snpprop_mat join cvterm on type_id=cvterm_id where cvterm.name IN ('codon_change','position in codon','substitution_type','non-coding')". " AND library_id IN ($lib_id) group by cvterm.name,value order by cvterm.name,value ";
  if (!empty($limit)) {
    $get_dbxref_sql .= " LIMIT $limit";
  }

  // GO EC KEGG
  $cvterm_dbresult = db_query($get_cvterm_sql);
  while ($cvterm_data = db_fetch_array($cvterm_dbresult)) {
    if (empty($limit) || !empty($limit) && $limit > count($cvterms[$cvterm_data['tag']])) {
      $cvterms[$cvterm_data['tag']][$cvterm_data['name']] = array(
        'count' => $cvterm_data['value_count'],
        'accession' => $cvterm_data['accession'],
        'label' => $cvterm_data['name'],
      );
    }
  }
  /*
   * For GO we want the GO slim
   */
  $cvterms = gmod_dbsf_db_execute('go', 'genes4all_explore_go_slims', $cvterms);
  gmod_dbsf_db_execute('chado');
  $graph_cvterms = array();
  if (!empty($cvterms)) {
    foreach ($cvterms as $cv => $term) {
      foreach ($term as $id => $data) {
        // skip the top level GO
        if ($id === $cv) {
          continue;
        }
        $label = $id.': '.$data['count'];
        if ($id !== $data['label']) {
          $data['label'] = rtrim($data['label'], ', ');
          $label = wordwrap($id.': '.$data['count']."\n(". $data['label'] .")");
        }
        $graph_cvterms[$cv][] = array(
          'y' => log10($data['count']),
          'tooltip' => $label,
        );
      }
    }
    unset($cvterms);
  }
  //IPR
  $ipr_count = 0;
  $ipr_x_legend = 'Domains';
  if (!empty($limit)) {
    $ipr_x_legend = "Top $limit domains";
  }
  $dbxref_dbresult = db_query($get_dbxref_sql);
  while ($dbxref_data = db_fetch_array($dbxref_dbresult)) {
    $ipr_count++;
    $dbxrefs[$dbxref_data['tag']][] = array(
      'y' => log10($dbxref_data['value_count']),
      'tooltip' => $dbxref_data['description'].': '.$dbxref_data['value_count']." (". $dbxref_data['accession'] .")",
    );
  }
  if (!empty($dbxrefs)) {
    if ($ipr_count > 10) {
      $ipr_max = 0;
      foreach ($dbxrefs as $db => $term) {
        foreach ($term as $data) {
          if ($data['y'] > $ipr_max) {
            $ipr_max = $data['y'];
          }
        }
      }
      if (!empty($ipr_max)) {
        $ipr_x_legend .= " with $base > ". sprintf('%.2f', $ipr_max / 3);
        foreach ($dbxrefs as $db => $term) {
          foreach ($term as $key => $data) {
            if ($data['y'] < ($ipr_max / 3)) {
              unset($dbxrefs[$db][$key]);
            }
          }
        }
      }
    }
  }
  //SNPs
  // we want to make the following graphs: boxplot of frequencies (difficult?) and
  // bar graph of number of syn/non-syn/non-coding/transition/transversions
  // // non-coding(TRUE), substitution_type(transition/transversion), position in codon (1,2,3) codon_change (non-synonymous,synonymous,unknown)
  $snp_width = array();
  $snp_dbresult = db_query($get_snp_sql);
  while ($snp_data = db_fetch_array($snp_dbresult)) {
    $tag = $snp_data['tag'];
    $name = $snp_data['name'];
    if ($name === 'unknown') {
      continue;
    }
    if (empty($snp_width[$tag])) {
      $snp_width[$tag] = 250;
    }
    if ($tag === 'non-coding') {
      $tag             = 'Effect on codon';
      $name            = 'Non coding';
      $snp_width[$tag] = 400;
    }
    elseif ($tag === 'codon_change') {
      $tag = 'Effect on codon';
    }
    elseif ($tag === 'substitution_type') {
      $tag = 'Type of substitution';
    }
    elseif ($tag === 'position in codon') {
      $tag = 'Codon position';
      if ($name == 1) {
        $name = '1st';
      }
      elseif ($name == 2) {
        $name = '2nd';
      }
      elseif ($name == 3) {
        $name = '3rd';
      }
    }
    $snps_array[$tag][] = array('y' => $snp_data['value_count'], 'tooltip' => $name.': '.$snp_data['value_count']);
  }
  return array(
    'graph_cvterms' => $graph_cvterms,
    'snp_width' => $snp_width,
    'snp_data' => $snps_array,
    'dbxrefs' => $dbxrefs,
  );
}

/*
 * Create graphs for genes4all/info/library using js (dojo)
 *
 * @see genes4all_explore_library_get_go()
 * @param $graph_data
 * return value from genes4all_explore_library_get_go()
 *
 * @return array
 * of HTML code (one per graph)
 */
function genes4all_explore_library_get_go_graphs($graph_data, $limit = 100, $base = 'log10') {
  if (empty($graph_data)){return;}
  $module_path = drupal_get_path('module', 'genes4all');
  drupal_add_css($module_path.'lib/dojo/resources/dojo.css', 'module', 'all');
  drupal_add_css($module_path.'/lib/dijit/themes/tundra/tundra.css', 'module');
  drupal_add_js($module_path.'/lib/dojo/dojo.js');

  $html_code=array();
  $js_code='
  dojo.require("dojox.charting.Chart2D");
  dojo.require("dojox.charting.themes.Distinctive");
  var dc = dojox.charting;
  dojo.require("dojox.charting.action2d.Tooltip");
  dojo.require("dojox.charting.action2d.Highlight");
  dojo.require("dojox.charting.widget.Legend");
  dojo.require("dojox.charting.action2d.Magnify");
  makeCharts = function() {
  ';
  if (!empty($graph_data['graph_cvterms']['biological_process'])) {
    $js_code.='
        var chart_bio = new dc.Chart2D("bio_graph");
        chart_bio.addPlot("default", {type: "Columns",gap: 5});
        chart_bio.addAxis("y",{includeZero: true, vertical: true});
        chart_bio.addSeries("Biological process", '.drupal_to_js($graph_data['graph_cvterms']['biological_process']).');
        chart_bio.setTheme(dojox.charting.themes.Distinctive);
        var anim_bio_tooltip = new dc.action2d.Tooltip(chart_bio, "default");
        var anim_bio_2 = new dc.action2d.Highlight(chart_bio, "default");
        chart_bio.render();
    ';
    $html_code['go'].='<h2>Biological process</h2><div id="bio_graph" style="width: 500px; height: 300px; margin: 5px auto 0px auto;"></div>';
  }
  if (!empty($graph_data['graph_cvterms']['molecular_function'])) {
    $js_code.='
        var chart_mf = new dc.Chart2D("mf_graph");
        chart_mf.addPlot("default", {type: "Columns",gap: 5});
        chart_mf.addAxis("y",{includeZero: true, vertical: true});
        chart_mf.addSeries("Molecular Function", '.drupal_to_js($graph_data['graph_cvterms']['molecular_function']).');
        chart_mf.setTheme(dojox.charting.themes.Distinctive);
        var anim_mf_tooltip = new dc.action2d.Tooltip(chart_mf, "default");
        var anim_mf_2 = new dc.action2d.Highlight(chart_mf, "default");
        chart_mf.render();
    ';
    $html_code['go'].='<h2>Molecular Function</h2><div id="mf_graph" style="width: 500px; height: 300px; margin: 5px auto 0px auto;"></div>';
  }
  if (!empty($graph_data['graph_cvterms']['cellular_component'])) {
    $js_code.='
        var chart_cc = new dc.Chart2D("cc_graph");
        chart_cc.addPlot("default", {type: "Columns",gap: 5});
        chart_cc.addAxis("y",{includeZero: true, vertical: true});
        chart_cc.addSeries("Cellular Component", '.drupal_to_js($graph_data['graph_cvterms']['cellular_component']).');
        chart_cc.setTheme(dojox.charting.themes.Distinctive);
        var anim_cc_tooltip = new dc.action2d.Tooltip(chart_cc, "default");
        var anim_cc_2 = new dc.action2d.Highlight(chart_cc, "default");
        chart_cc.render();
    ';
    $html_code['go'].='<h2>Cellular Component</h2><div id="cc_graph" style="width: 500px; height: 300px; margin: 5px auto 0px auto;"></div>';
  }
  if (!empty($graph_data['graph_cvterms']['EC'])) {
    $js_code.='
        var chart_ec = new dc.Chart2D("ec_graph");
        chart_ec.addPlot("default", {type: "Columns",gap: 1});
        chart_ec.addAxis("y",{includeZero: true, vertical: true});
        chart_ec.addSeries("Enzyme Classification", '.drupal_to_js($graph_data['graph_cvterms']['EC']).');
        chart_ec.setTheme(dojox.charting.themes.Distinctive);
        var anim_ec_tooltip = new dc.action2d.Tooltip(chart_ec, "default");
        var anim_ec_2 = new dc.action2d.Highlight(chart_ec, "default", {highlight: "gold"});
        chart_ec.render();
    ';
    $html_code['ec']='<h2>Enzyme Classifications</h2><div id="ec_graph" style="width: 700px; height: 300px; margin: 5px auto 0px auto;"></div>';
  }
  if (!empty($graph_data['graph_cvterms']['KEGG_PATHWAY'])) {
    $js_code.='
        var chart_kegg = new dc.Chart2D("kegg_graph");
        chart_kegg.addPlot("default", {type: "Columns",gap: 1});
        chart_kegg.addAxis("y",{includeZero: true, vertical: true});
        chart_kegg.addSeries("KEGG Pathway", '.drupal_to_js($graph_data['graph_cvterms']['KEGG_PATHWAY']).');
        chart_kegg.setTheme(dojox.charting.themes.Distinctive);
        var anim_kegg_tooltip = new dc.action2d.Tooltip(chart_kegg, "default");
        var anim_kegg_2 = new dc.action2d.Highlight(chart_kegg, "default", {highlight: "gold"});
        chart_kegg.render();
    ';
    $html_code['kegg']='<h2>KEGG Pathway</h2><div id="kegg_graph" style="width: 700px; height: 300px; margin: 5px auto 0px auto;"></div>';
  }
  if (!empty($graph_data['dbxrefs']['InterPro'])) {
    $js_code.='
        var chart_ip = new dc.Chart2D("ip_graph");
        chart_ip.addPlot("default", {type: "Columns",gap: 1});
        chart_ip.addAxis("y",{includeZero: true, vertical: true});
        chart_ip.addSeries("InterPro Domains", '.drupal_to_js($graph_data['dbxrefs']['InterPro']).');
        chart_ip.setTheme(dojox.charting.themes.Distinctive);
        var anim_ip_tooltip = new dc.action2d.Tooltip(chart_ip, "default");
        var anim_ip_2 = new dc.action2d.Highlight(chart_ip, "default", {highlight: "gold"});
        chart_ip.render();
    ';
    $html_code['ip']='<h2>InterPro Domains</h2><div id="ip_graph" style="width: 700px; height: 300px; margin: 5px auto 0px auto;"></div>';
  }
  if (!empty($graph_data['snp_data'])) {
    $html_code['snp'].='<h2>SNP markers</h2>';
    foreach ($graph_data['snp_data'] as $series=>$points){
      $i++;
      $js_code.='
        var chart_snp_'.$i.' = new dc.Chart2D("snp_graph_'.$i.'");
        chart_snp_'.$i.'.addAxis("y",{includeZero: true, vertical: true});
        chart_snp_'.$i.'.addPlot("default", {type: "Columns",gap: 10});
        chart_snp_'.$i.'.addSeries("SNP markers", '.drupal_to_js($points).');
        chart_snp_'.$i.'.setTheme(dojox.charting.themes.Distinctive);
        var anim_snp_'.$i.'_tooltip = new dc.action2d.Tooltip(chart_snp_'.$i.', "default");
        var anim_snp_'.$i.'_2 = new dc.action2d.Highlight(chart_snp_'.$i.', "default");
        chart_snp_'.$i.'.render();
      ';
      $html_code['snp'].='<h3>'.$series.'</h3><div id="snp_graph_'.$i.'" style="width: 300px; height: 200px; margin: 5px auto 0px auto;"></div>';
    }
  }

  $js_code .='};dojo.addOnLoad(makeCharts);';
  drupal_add_js($js_code,'inline');
  return $html_code;
}

/**
 *
 * @param $cvterms
 *
 * @return unknown_type
 */
function genes4all_explore_go_slims($cvterms) {
  if (empty($cvterms)) {
    return;
  }

  /*
   * so i converted the GO's mysql schema to postgres and loaded the term database
   * this query will get the immediate parent so we will loop through till we get the top level.
   */



  //GO database specific SQLs so no {} fields
  $go_sql = "select object.name as object from graph_path " . " join term as subject on subject.id=graph_path.term2_id join term as object on object.id=graph_path.term1_id where " . " subject.acc='%s' and distance=(select max(distance)-2 from graph_path " . " join term as subject on subject.id=graph_path.term2_id where subject.acc='%s')";

  foreach ($cvterms as $cv => $term) {
    if ($cv == 'biological_process' || $cv == 'molecular_function' || $cv == 'cellular_component') {
      foreach ($term as $id => $data) {
        $data['accession'] = "GO:". $data['accession'];
        $result            = db_query($go_sql, $data['accession'], $data['accession']);
        $is_child          = 0;
        while ($dbdata = db_fetch_array($result)) {
          $slim_id = $dbdata['object'];
          if ($slim_id === 'all') {
            continue;
          }
          $cvterms[$cv][$slim_id]['count'] += $data['count'];
          $cvterms[$cv][$slim_id]['accession'] .= " ". $data['accession'];
          $cvterms[$cv][$slim_id]['label'] .= $id .", ";
          $is_child = 1;
        }
        if ($is_child === 1) {
          unset($cvterms[$cv][$id]);
        }
      }
    }
  }
  return ($cvterms);
}


/*
 * to be called through gmod_dbsf_db_execute
 * Take data array prepared
 * from submitted form and return the SQL to fetch features according to offset and limit
 */
function genes4all_explore_results_prepare_sql($data) {
  //Chado specific function
  //make sure the db is reset to drupal if redirecting out of this function
  // the use of return is safe however.
  if (!isset($data) || !isset($data['jointype'])) {
    return "No data offered";
  }
  if ($data['jointype'] !== 'intersect' && $data['jointype'] !== 'union') {
    return "Wrong join type";
  }
  if (
  !isset($data['dbxref_array']) && !isset($data['cvterm_array'])
  && !isset($data['blast_string'])
  && !isset($data['library']) && !isset($data['markers'])
  ) {
    return "Nothing to do...";
  }

  // selects which organisms to look for
  $organism_sql = '';
  // The following are selecting for a feature id :
  // queries the dbxref_organism_mat
  $dbxref_sql = '';
  // queries the cvterm_organism_mat
  $cvterm_sql = '';
  // complex query for the cvterm to get blasts. TODO simply by using keywords. DONOT load blast.gff into chado therefore.
  $blast_sql = '';
  $blast_string = '';
  // library
  $lib_sql = '';
  // markers
  $marker_sql = '';
  $boolean_join = '';
  if ($data['jointype'] !== 'intersect') {
    $boolean_join = ' AND ';
  }
  elseif ($data['jointype'] !== 'union') {
    $boolean_join = ' OR ';
  }
  //join different queries together
  $jointype = ' '. $data['jointype'] .' ';

  if (!empty($data['blast_string'])) {
    $blast_string = check_plain($data['blast_string']);
    $blast_string = trim($blast_string);
  }
  if (!empty($data['organism_array'])) {
    $species_query = '';
    foreach ($data['organism_array'] as $species) {
      $species_query .= "$species,";
    }
    $species_query = rtrim($species_query, ",");
    if (!empty($species_query)) {
      $organism_sql = " AND organism_id IN ($species_query) ";
    }
  }

  if (!empty($data['library'])) {

    //rest is cvterm.name and cv.name
    $lib_query_cv = '';
    $lib_query_name = '';
    //TODO
    $lib_sql = "SELECT feature_id from library_feature_contigs_mat where  ";

    //  ."(select cvterm_id from cvterm join cv on cvterm.cv_id=cv.cv_id where cvterm.name='contig' and cv.name='sequence')"
    //  ." AND feature_id IN (SELECT srcfeature_id from featureloc where feature_id "
    //  ."IN (SELECT feature_id from featureloc where srcfeature_id IN "
    //  ."(SELECT feature_id from library_feature where ";

    foreach ($data['library'] as $tag => $type) {
      if (empty($type)) {
        continue;
      }
      if ($tag == 'library_name') {
        foreach ($type as $value) {
          $lib_query_name .= "'$value',";
        }
      }
      else {
        foreach ($type as $key => $value) {
          if (empty($value)) {
            continue;
          }
          $lib_query_cv .= $value .',';
        }
      }
    }

    if (!empty($lib_query_cv)) {
      $lib_query_cv = rtrim($lib_query_cv, ",");
      $lib_sql .= ' library_id IN (select library_id from library_cvterm '. "  where cvterm_id IN ($lib_query_cv) )";
      // if both
      if (!empty($lib_query_name)) {
        $lib_query_name = rtrim($lib_query_name, ",");
        $lib_sql .= $boolean_join;
        $lib_sql .= " library_id IN (SELECT library_id from library where name IN ". " ($lib_query_name))";
      }
    }
    elseif (!empty($lib_query_name)) {
      $lib_query_name = rtrim($lib_query_name, ",");
      $lib_sql .= " library_id IN (SELECT library_id from library where name IN ". " ($lib_query_name))";
    }
  }
  //MARKERS;
  if (!empty($data['markers'])) {
    if (!empty($data['markers']['SNP'])) {
      //TODO CONSIDER using library_feature_contigs_mat and library_snpprop_mat
      $marker_sql = 'SELECT srcfeature_id as feature_id from featureloc where feature_id IN (' . ' select feature_id from feature where type_id =' . " (select cvterm_id from cvterm where name='SNP' and " . " cv_id = (select cv_id from cv where name='sequence'))";

      //TODO now add the various acquired limits
      // which are stored in featureprop. the cv of the type id is feature_property
      /* main snp:
      total_reads_position  int
      codon_change          synonymous/non-synonymous
      substitution_type   transition/transversion
      organism            abbrev of organism
      minor_frequency     float (2 dec)
      alleles             char(1)
      alleles             char(1)
      position in codon   int 1-3
      quality_SNP         int
      orf ID                orf id which ref for orf snp
      --- orf snp:
      alleles              char(1)
      alleles              char(1)
      total_reads_position  int
      codon_change          synonymous/non-synonymous
      position in codon   int (1-3)
      quality_SNP          int
      substitution_type   transition/transversion
      organism            abbrev of organism

      $marker_sql='';
      foreach ($data['markers'] as $key=>$value)
      {// key is SNP
      foreach ($value as $key2=>$value2){
      $marker_sql.="$key2=>$value2 ;";
      }
      }

      'force' => t('Any'),
      'syn' => t('Is coding & causes a synonymous change'),
      'non_syn' => t('Is coding & a non-synonymous change'),
      'transition' => t('Is a transition'),
      'transversion' => t('Is a transversion'),
      'coding' => t('Is part of a open-reading-frame'),
      'non-coding' => t('Is not part of a open-reading-frame'),
      [SNP] = array(
      force=>force,synonymous=>synonymous,non-synonymous=>non-synonymous,
      transition=>transition,transversion=>transversion,
      coding=>coding,non-coding=>non-coding
      )
      [SSR] = array()force=>force
      */



      $intercept = ' INTERSECT ';
      $intercept_flag = '';

      foreach ($data['markers']['SNP'] as $tag) {
        switch ($tag) {
          case 'force':
            break2;
          case 'synonymous':
            if (empty($intercept_flag)) {
              $marker_sql .= $intercept;
              $intercept_flag = 1;
            }
            $marker_sql .= " SELECT feature_id from featureprop where " . "type_id=(select cvterm_id from cvterm where cv_id = (select cv_id from cv where name='feature_property')" . " and name='codon_change' ) AND value='synonymous' ";
            if (
            !empty($data['markers']['SNP']['non-synonymous'])
            || !empty($data['markers']['SNP']['transition'])
            || !empty($data['markers']['SNP']['transversion'])
            || !empty($data['markers']['SNP']['coding'])
            || !empty($data['markers']['SNP']['non-coding'])
            ) {
              $marker_sql .= $jointype;
              $jointype_flag = '';
            }
            break;

          case 'non-synonymous':
            if (empty($intercept_flag)) {
              $marker_sql .= $intercept;
              $intercept_flag = 1;
            }
            $marker_sql .= " SELECT feature_id from featureprop where " . "type_id=(select cvterm_id from cvterm where cv_id = (select cv_id from cv where name='feature_property')" . " AND name='codon_change' ) AND value='non-synonymous' ";
            if (
            !empty($data['markers']['SNP']['transition'])
            || !empty($data['markers']['SNP']['transversion'])
            || !empty($data['markers']['SNP']['coding'])
            || !empty($data['markers']['SNP']['non-coding'])
            ) {
              $marker_sql .= $jointype;
              $jointype_flag = '';
            }
            break;

          case 'transition':
            if (empty($intercept_flag)) {
              $marker_sql .= $intercept;
              $intercept_flag = 1;
            }
            $marker_sql .= " SELECT feature_id from featureprop where type_id=(select cvterm_id from cvterm where cv_id = (select cv_id from cv where name='feature_property')" . " AND name='substitution_type' ) AND value='transition' ";
            if (
            !empty($data['markers']['SNP']['transversion'])
            || !empty($data['markers']['SNP']['coding'])
            || !empty($data['markers']['SNP']['non-coding'])
            ) {
              $marker_sql .= $jointype;
              $jointype_flag = '';
            }
            break;

          case 'transversion':
            if (empty($intercept_flag)) {
              $marker_sql .= $intercept;
              $intercept_flag = 1;
            }
            $marker_sql .= " SELECT feature_id from featureprop where type_id=(select cvterm_id from cvterm where cv_id = (select cv_id from cv where name='feature_property')" . " AND name='substitution_type' ) AND value='transversion' ";
            if (
            !empty($data['markers']['SNP']['coding'])
            || !empty($data['markers']['SNP']['non-coding'])
            ) {
              $marker_sql .= $jointype;
              $jointype_flag = '';
            }
            break;

          case 'coding':
            if (empty($intercept_flag)) {
              $marker_sql .= $intercept;
              $intercept_flag = 1;
            }
            $marker_sql .= " SELECT feature_id from featureprop where " . "type_id=(select cvterm_id from cvterm where cv_id = (select cv_id from cv where name='feature_property') " . " AND name='orf ID' ) ";
            if (
            !empty($data['markers']['SNP']['non-coding'])
            ) {
              $marker_sql .= $jointype;
              $jointype_flag = '';
            }
            break;

          case 'non-coding':
            if (empty($intercept_flag)) {
              $marker_sql .= $intercept;
              $intercept_flag = 1;
            }
            $marker_sql .= " SELECT feature_id from featureprop where " . "type_id=(select cvterm_id from cvterm where cv_id = (select cv_id from cv where name='feature_property') " . " AND name='non-coding' ) ";
            break;
        }
        // end options
      }
      //END tag


      $marker_sql .= ')';
    }
    // END SNP
  }
  //END markers
  //DBXREFs (e.g IPRSCAN)
  if (!empty($data['dbxref_array'])) {
    $dbxref_query = '';
    // UNION
    if ($data['jointype'] == 'union') {
      foreach ($data['dbxref_array'] as $dbxref) {
        if (!empty($dbxref)) {
          $dbxref_query .= "$dbxref,";
        }
      }
      $dbxref_query = rtrim($dbxref_query, ",");
      $dbxref_sql = "SELECT feature_id FROM feature_dbxref WHERE dbxref_id IN ($dbxref_query) ";
    }

    //INTERSECT
    // rather than making a single sql query
    // the most cost-efficient approach i can make is to call the subquery for each dbxref id and put values in hash
    // the value of the hash element increments for each hit. then select only those hits that have maximum number of hit.
    elseif ($data['jointype'] == 'intersect') {
      $temp_array = array();
      // table is well-indexed so no need to limit by organism at this stage.
      $temp_sql = 'select feature_id FROM feature_dbxref WHERE dbxref_id=%d';
      $requested_hits = 0;
      foreach ($data['dbxref_array'] as $dbxref) {
        if (!empty($dbxref)) {
          $requested_hits++;
          $temp_result = db_query($temp_sql, $dbxref);
          while ($temp_data = db_fetch_array($temp_result)) {
            $temp_array[$temp_data['feature_id']]++;
          }
        }
      }
      // purge thost that don't match all terms
      // trimmed array of the above
      $temp_feat_array = array();
      foreach ($temp_array as $key => $value) {
        if ($value === $requested_hits) {
          $temp_feat_array[$key] = $value;
        }
      }
      if (!empty($temp_feat_array)) {
        $temp_str = gmod_dbsf_array2csv($temp_feat_array, NULL, ',', TRUE);
        $dbxref_sql = "SELECT feature_id from feature WHERE feature_id IN ($temp_str) ";
      }
    }
  }
  //CVTERMS (EC,GO,KEGG)
  if (!empty($data['cvterm_array'])) {
    if ($data['jointype'] == 'union') {
      $cvterm_query = gmod_dbsf_array2csv($data['cvterm_array']);
      $cvterm_sql = "SELECT feature_id from feature_cvterm where cvterm_id IN ($cvterm_query) ";
    }
    elseif ($data['jointype'] == 'intersect') {
      // SEE dbxref above for reasoning of method:
      $temp_array     = array();
      $requested_hits = 0;
      $temp_sql       = 'select feature_id FROM feature_cvterm WHERE cvterm_id=%d';
      foreach ($data['cvterm_array'] as $cvterm) {
        if (!empty($cvterm)) {
          $requested_hits++;
          $temp_result = db_query($temp_sql, $cvterm);
          while ($temp_data = db_fetch_array($temp_result)) {
            $temp_array[$temp_data['feature_id']]++;
          }
        }
      }
      $temp_feat_array = array();
      foreach ($temp_array as $key => $value) {
        if ($value === $requested_hits) {
          $temp_feat_array[$key] = $value;
        }
      }
      if (!empty($temp_feat_array)) {
        $temp_str = gmod_dbsf_array2csv($temp_feat_array, NULL, ',', TRUE);
        $cvterm_sql = "SELECT feature_id from feature WHERE feature_id IN ($temp_str) ";
      }
    }
  }

  if (isset($blast_string) && !empty($blast_string)) {
    //TODO improve blast search method
    // blast_string becomes lower case
    // blast_string to explode to array using spaces as delimited
    // then build a new sql based string such 'term1','term2' etc
    // fix sql command for blast. we are interested 'select feature_id from featureprop where value'
    // first we do a AND AND search and then UNION the OR OR search
    // one problem with the UNION approach is the sorting as the data are always resorted....
    $blast_string = strtolower($blast_string);
    //remove double spaces and trim
    $blast_string    = preg_replace('/\s\s+/', ' ', $blast_string);
    $blast_string    = trim($blast_string);
    $blast_values    = explode(" ", $blast_string);
    $blast_value_str = '';
    foreach ($blast_values as $value) {
      $blast_value_str .= "'$value',";
    }
    $blast_value_str = rtrim($blast_value_str, ",");

    // the following is an OR
    $blast_sql = "SELECT feature_id from blast_organism_mat where value IN ($blast_value_str) ";

    /**
     //   'SELECT feature_id from featureprop where type_id='
     //   ."(select cvterm_id from cvterm where name='inferred from electronic annotation' and is_obsolete=0 limit 1) "
     //  ."and value IN ($blast_value_str) "
     //  ;
     //  $blast_sql=
     //' select feature_id from feature where feature_id IN '
     //.'(select srcfeature_id from featureloc where feature_id IN '
     //.' (select feature_id from featureprop where type_id='
     //."(select cvterm_id from cvterm where name='Note' and cv_id= "
     //."(select cv_id from cv where name='feature_property')) and value~*'"
     //.$blast_string
     //."')) and type_id IN"
     //." (select cvterm_id from cvterm where name IN ('polypeptide','contig') and cv_id=(select cv_id from cv where name='sequence')) ";
     */
  }

  $sql = '';
  if (!empty($dbxref_sql)) {
    $sql = $dbxref_sql;
    if (!empty($cvterm_sql)) {
      $sql .= $jointype . $cvterm_sql;
    }
    if (!empty($blast_sql)) {
      $sql .= $jointype . $blast_sql;
    }
    if (!empty($lib_sql)) {
      $sql .= $jointype . $lib_sql;
    }
    if (!empty($marker_sql)) {
      $sql .= $jointype . $marker_sql;
    }
  }
  elseif (!empty($cvterm_sql)) {
    $sql = $cvterm_sql;
    if (!empty($blast_sql)) {
      $sql .= $jointype . $blast_sql;
    }
    if (!empty($lib_sql)) {
      $sql .= $jointype . $lib_sql;
    }
    if (!empty($marker_sql)) {
      $sql .= $jointype . $marker_sql;
    }
  }
  elseif (!empty($blast_sql)) {
    $sql = $blast_sql;
    if (!empty($lib_sql)) {
      $sql .= $jointype . $lib_sql;
    }
    if (!empty($marker_sql)) {
      $sql .= $jointype . $marker_sql;
    }
  }
  elseif (!empty($lib_sql)) {
    $sql = $lib_sql;
    if (!empty($marker_sql)) {
      $sql .= $jointype . $marker_sql;
    }
  }
  elseif (!empty($marker_sql)) {
    $sql = $marker_sql;
  }

  // if nothing has been asked, exit.
  else {
    db_set_active('default');
    drupal_goto('genes4all/noresults');
    return;
  }
  //TODO Further TESTING
  $sql = "select feature_id from feature where feature_id IN (". $sql .")";
  if (!empty($organism_sql)) {
    $sql .= $organism_sql;
  }
  $sql .= "order by type_id,length(feature.name),feature.name";
  if (!empty($sql)) {
    return $sql;
  }
}

