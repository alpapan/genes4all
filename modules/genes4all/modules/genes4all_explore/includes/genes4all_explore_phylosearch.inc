<?php
// $Id$

/**
 * @file
 * include file for searching for features using a phylogeny
 * and annotation
 */

/**
 * Render the main page and form
 *
 * @return unknown_type
 */
function genes4all_explore_page_phylosearch() {
  $species_array = array();
  if (!empty($_GET['species'])){
    $species_array =  explode(',',check_plain($_GET['species']));
  }
  $table_array = '';
  $encoded     = '';
  $data        = array();
  $level = !empty($_GET['level']) ? check_plain($_GET['level']) :(int)1;

  if (!empty($_GET['bookmark'])) {
    $encoded = $_GET['bookmark'];
    $data = gmod_dbsf_decode($encoded);
    if (!empty($data)) {
      // Build SQL argument without limits
      $sql = gmod_dbsf_db_execute('chado', 'genes4all_explore_results_prepare_sql', $data);
      // Run SQL (with limits), fetch Data. If no data, then this function stops as there is a redirect to genes4all_explore_noresults()
      $table_array = gmod_dbsf_db_execute('chado', 'genes4all_explore_results_prepare_data', array('sql' => $sql, 'encoded' => $encoded, 'level' => $level,'limit'=>50));
      if (!empty($table_array)) {
        $rendered_form = drupal_get_form('genes4all_explore_result_form', $table_array);
        $species_list = $_SESSION['species'];
        unset($_SESSION['species']);
        return '<p>Perform a new '. l('advanced search', 'genes4all/search/phylogeny')
        .' or another query the '
        .l('same species','genes4all/search/phylogeny',array('query'=>array('species'=>$species_list)))
        .'.</p>'. $rendered_form;
      }
    }
  }
  //Otherwise present form to select species or to search within a species selection
  if (!empty($species_array)){
    return drupal_get_form('genes4all_explore_search_species_form',$species_array);
  }
  //return drupal_get_form('genes4all_explore_filters_form');
  return drupal_get_form('genes4all_explore_select_species_form');
}

function genes4all_explore_species_friendly($species_array){
  if (!empty($species_array)) {
    $species_array_friendly = array();
    $limit = implode(',',$species_array);
    if (!empty($limit)) {
      $return_array = gmod_dbsf_db_execute('chado', 'gmod_dbsf_getspecies', $limit);
      $array = $return_array['species_array'];
      if (empty($array)) {
        return;
      }
      foreach ($array as $name) {
        $species_array_friendly[] = $name;
      }
    }
    else {
      $species_array_friendly[] = "None selected.";
    }

    $form['filters']['species'] = array(
        '#type' => 'fieldset',
        '#title' => t('Limit by species'),
        '#collapsible' => FALSE,
        '#weight' => -1,
    );

    $form['filters']['species']['species_show'] = array(
        '#title' => 'Species',
        '#description' => t('You have selected these species.'),
        '#type' => 'item',
        '#value' => theme('item_list', $species_array_friendly, '', 'ul', array('class' => 'italics')),
        '#weight' => 1,
    );
    $form['filters']['species']['species_select'] = array(
        '#title' => 'Species',
        '#type' => 'value',
        '#value' => $species_array,
    );
  }
  return $form;
}

function genes4all_explore_cv_friendly($species_array){
  $return_array_cv = gmod_dbsf_db_execute('chado', 'genes4all_explore_getcvs', $species_array);
  $dbxref_array = array();
  $iea_array    = array();
  $ipr_array    = array(
      'active_site' => array('None available'),
      'binding_site' => array('None available'),
      'conserved_site' => array('None available'),
      'domain' => array('None available'),
      'family' => array('None available'),
      'ptm' => array('None available'),
      'region' => array('None available'),
      'repeat' => array('None available'),
  );
  $go_array = array(
      'process' => array('None available'),
      'location' => array('None available'),
      'function' => array('None available'),
  );
  $ec_array      = array('None available');
  $kegg_array    = array('None available');
  $library_array = array('None available');

  //CACHING library user specific?
  global $user;
  $user_id=$user->uid;
  $md5=md5(implode(',',(array)$species_array));
  $library_array = cache_get('LIB_'.$user_id.'_'.$md5,'cache_genes4all_explore');
  if (empty($library_array)){
    $library_array = gmod_dbsf_db_execute('chado', 'gmod_dbsf_chado_getlibraries_tags', array('organism_id' => $species_array));
    cache_set('LIB_'.$user_id.'_'.$md5, $library_array , 'cache_genes4all_explore');
  }else{
    $library_array = $library_array ->data;
  }

  // we are calling this before storing the species array on purpose
  $return_array_db = gmod_dbsf_db_execute('chado', 'genes4all_explore_getdbxrefs', $species_array);
  // ALL dbxref (incl species, iprscan)
  if (!empty($return_array_db['dbxref_array'])) {
    $dbxref_array = $return_array_db['dbxref_array'];
  }
  // ALL CVs (incl species, go, ec, kegg)
  if (!empty($return_array_cv['cv_array'])) {
    $iea_array = $return_array_cv['cv_array'];
  }

  //KEGG
  if (!empty($iea_array['KEGG_PATHWAY'])) {
    $kegg_array = $iea_array['KEGG_PATHWAY'];
    asort($kegg_array);
  }
  //GO
  $go_array = array();
  if (!empty($iea_array['biological_process'])) {
    $go_array['process'] = $iea_array['biological_process'];
    asort($go_array['process']);
  }
  else {
    $go_array['process'] = array('None available');
  }
  if (!empty($iea_array['cellular_component'])) {
    $go_array['location'] = $iea_array['cellular_component'];
    asort($go_array['location']);
  }
  else {
    $go_array['location'] = array('None available');
  }
  if (!empty($iea_array['molecular_function'])) {
    $go_array['function'] = $iea_array['molecular_function'];
    asort($go_array['function']);
  }
  else {
    $go_array['function'] = array('None available');
  }

  //EC
  if (!empty($iea_array['EC'])) {
    $ec_array = $iea_array['EC'];
    asort($ec_array);
  }

  //IPR
  $ipr_array = array();
  if (!empty($dbxref_array['InterPro']['Active_site'])) {
    $ipr_array['active_site'] = $dbxref_array['InterPro']['Active_site'];
    asort($ipr_array['active_site']);
  }
  else {
    $ipr_array['active_site'] = array('None available');
  }
  if (!empty($dbxref_array['InterPro']['Binding_site'])) {
    $ipr_array['binding_site'] = $dbxref_array['InterPro']['Binding_site'];
    asort($ipr_array['binding_site']);
  }
  else {
    $ipr_array['binding_site'] = array('None available');
  }
  if (!empty($dbxref_array['InterPro']['Conserved_site'])) {
    $ipr_array['conserved_site'] = $dbxref_array['InterPro']['Conserved_site'];
    asort($ipr_array['conserved_site']);
  }
  else {
    $ipr_array['conserved_site'] = array('None available');
  }
  if (!empty($dbxref_array['InterPro']['Domain'])) {
    $ipr_array['domain'] = $dbxref_array['InterPro']['Domain'];
    asort($ipr_array['domain']);
  }
  else {
    $ipr_array['domain'] = array('None available');
  }
  if (!empty($dbxref_array['InterPro']['Family'])) {
    $ipr_array['family'] = $dbxref_array['InterPro']['Family'];
    asort($ipr_array['family']);
  }
  else {
    $ipr_array['family'] = array('None available');
  }
  if (!empty($dbxref_array['InterPro']['PTM'])) {
    $ipr_array['ptm'] = $dbxref_array['InterPro']['PTM'];
    asort($ipr_array['ptm']);
  }
  else {
    $ipr_array['ptm'] = array('None available');
  }
  if (!empty($dbxref_array['InterPro']['Region'])) {
    $ipr_array['region'] = $dbxref_array['InterPro']['Region'];
    asort($ipr_array['region']);
  }
  else {
    $ipr_array['region'] = array('None available');
  }
  if (!empty($dbxref_array['InterPro']['Repeat'])) {
    $ipr_array['repeat'] = $dbxref_array['InterPro']['Repeat'];
    asort($ipr_array['repeat']);
  }
  else {
    $ipr_array['repeat'] = array('None available');
  }
  return array('ipr'=>$ipr_array,'ec'=>$ec_array,'go'=>$go_array,'kegg'=>$kegg_array,'library'=>$library_array);
}

/**
 * Implementation of hook_form to query data after a species has been selected
 */
function genes4all_explore_search_species_form($form_state,$species_array){
  global $user;
  $user_id=$user->uid;
  if (empty($species_array)){
    drupal_set_message('Erroneous selection of species.','error');
    return '';
  }
  $arrayed_data = genes4all_explore_cv_friendly($species_array);

  // Now we have all the data, we can build the form
  $form = array();
  $form['filters'] = array(
        '#type' => 'fieldset',
      '#title' => t('Fetch'),
      '#collapsible' => FALSE,
      '#description'=>t('A feature is a genomics object in the database. Please select a method to query the database. Hover your mouse over the buttons to see a detailed description.'),
      '#weight' => -1,
  )
  + genes4all_explore_species_friendly($species_array);
  $form['filters']['match']= array(
      	'#weight' => 25,
        '#type'=>'fieldset',
        '#title' =>t(' features matching '),
        '#collapsible' => FALSE,
      	'grab' => array(
        	'union' => array(
            '#attributes' => array('title' => t('Get features which match any one of the above selected terms.')),
            '#type' => 'submit',
            '#value' => t('ANY of your selected criteria (OR)'),
            '#weight' => 1,
      ),
          'intersect' => array(
            '#attributes' => array('title' => t('Get features which match all terms you selected above.')),
            '#type' => 'submit',
          	'#value' => t('ALL of your selected criteria (AND)'),
            '#weight' => 2,
      ),
      ),
      );

  $form['filters']['IEA'] = array(
      '#type' => 'fieldset',
      '#title' => t('Inference by Electronic Annotation'),
      '#description' => t('Inference by Electronic Annotation is assigned by similarity to known proteins (UniProt) ' . 'which have an ontology term assigned to them by a curator.'),
      '#collapsible' => TRUE,
      '#collapsed' => TRUE,
      '#weight' => 5,
      'GO' => array(
        '#type' => 'fieldset',
        '#title' => t('Gene Ontology'),
        '#description' => t('The GO project has developed three structured controlled vocabularies (ontologies)' . ' that describe gene products in terms of their associated biological processes, cellular components ' . 'and molecular functions in a species-independent manner. The use of GO terms by collaborating databases ' . 'facilitates uniform queries across them. The controlled vocabularies are structured so that they can ' . 'be queried at different levels: for example, you can use GO to find all the gene products in the mouse ' . 'genome that are involved in signal transduction, or you can zoom in on all the receptor tyrosine kinases.'),
        '#collapsible' => TRUE,
        '#collapsed' => TRUE,
        '#weight' => 1,
        'P' => array(
          '#type' => 'fieldset',
          '#title' => t('Biological Process'),
          '#description' => t('A biological process is series of events accomplished by one or more ordered' . ' assemblies of molecular functions. Examples of broad biological process terms are cellular physiological ' . 'process or signal transduction. Examples of more specific terms are pyrimidine metabolic process ' . 'or alpha-glucoside transport. It can be difficult to distinguish between a biological process ' . 'and a molecular function, but the general rule is that a process must have more than one distinct steps.'),
          '#collapsible' => TRUE,
          '#collapsed' => TRUE,
          '#weight' => 1,
          'go_process_select' => array(
            '#description' => t('Please select one or more Biological Process GO terms (using the Shift/Ctl/Option key). Use the Ctl/Option key to deselect.'),
            '#multiple' => '1',
            '#type' => 'select',
            '#options' => $arrayed_data['go']['process'],
            '#weight' => 6,
  ),
  ),
        'F' => array(
          '#type' => 'fieldset',
          '#title' => t('Molecular Function'),
          '#description' => t('Molecular function describes activities, such as catalytic or binding activities,' . ' that occur at the molecular level. GO molecular function terms represent activities rather ' . 'than the entities (molecules or complexes) that perform the actions, and do not specify where or when, ' . 'or in what context, the action takes place. Molecular functions generally correspond to activities that' . ' can be performed by individual gene products, but some activities are performed by assembled complexes' . ' of gene products. Examples of broad functional terms are catalytic activity, transporter activity, ' . 'or binding; examples of narrower functional terms are adenylate cyclase activity or Toll receptor binding.'),
          '#collapsible' => TRUE,
          '#collapsed' => TRUE,
          '#weight' => 5,
          'go_function_select' => array(
            '#description' => t('Please select one or more Molecular Function GO terms (using the Shift/Ctl/Option key). Use the Ctl/Option key to deselect.'),
            '#multiple' => '1',
            '#type' => 'select',
            '#options' => $arrayed_data['go']['function'],
            '#weight' => 3,
  ),
  ),
        'L' => array(
          '#type' => 'fieldset',
          '#title' => t('Cellular localization'),
          '#description' => t('Localization to a component of a cell, but with the proviso ' . 'that it is part of some larger object; this may be an anatomical structure ' . '(e.g. rough endoplasmic reticulum or nucleus) or a gene product group ' . '(e.g. ribosome, proteasome or a protein dimer). ' . 'See the !url on the cellular component ontology for more details.',
  array('!url' => l('documentation', 'http://www.geneontology.org/GO.component.guidelines.shtml',
  array('absolute' => TRUE, 'attributes' => array('target' => '_blank'))
  ))
  ),
          '#collapsible' => TRUE,
          '#collapsed' => TRUE,
          '#weight' => 1,
          'go_location_select' => array(
            '#description' => t('Please select one or more Cellular localization GO terms (using the Shift/Ctl/Option key). Use the Ctl/Option key to deselect.'),
            '#multiple' => '1',
            '#type' => 'select',
            '#options' => $arrayed_data['go']['location'],
            '#weight' => 1,
  ),
  ),
  ),
      'KEGG' => array(
        '#type' => 'fieldset',
        '#title' => t('Pathways from the Kyoto Encyclopedia'),
        '#description' => t('Collection of manually drawn pathway maps representing our knowledge ' . 'on the molecular interaction and reaction networks. We assign KEGG terms by IEA ' . '- similarity to annotated proteins. KEGG has a biomedical focus and therefore many ' . 'pathways will make little sense to insects.'
        ),
        '#collapsible' => TRUE,
        '#collapsed' => TRUE,
        '#weight' => 10,
        'kegg_select' => array(
          '#description' => t('Please select one or more KEGG pathway terms (using the Shift/Ctl/Option key). Use the Ctl/Option key to deselect.'),
          '#multiple' => '1',
          '#type' => 'select',
          '#options' => $arrayed_data['kegg'],
          '#weight' => 1,
        ),
        ),
      'EC' => array(
        '#type' => 'fieldset',
        '#title' => t('Enzyme Commission'),
        '#description' => t('EC numbers are assigned to enzymes according to the reaction they catalyse.'),
        '#collapsible' => TRUE,
        '#collapsed' => TRUE,
        '#weight' => 5,
        'ec_select' => array(
          '#description' => t('Please select one or more Enzyme Classification terms (using the Shift/Ctl/Option key). Use the Ctl/Option key to deselect.'),
          '#multiple' => '1',
          '#type' => 'select',
          '#title' => 'EC',
          '#options' => $arrayed_data['ec'],
          '#weight' => 8,
        ),
        ),
        );
        $form['filters']['motif'] = array(
      '#type' => 'fieldset',
      '#title' => t('Protein motif signatures'),
      '#description' => t('Another method of annotation is to detect conserved signatures of certain motifs. Please note that due to the high computational power required, this is not provided for all species.'),
      '#collapsible' => TRUE,
      '#collapsed' => TRUE,
      '#weight' => 2,
      'IPR' => array(
        '#type' => 'fieldset',
        '#title' => t('InterProScan'),
        '#description' => t('InterPro is a database of protein families, domains, regions, repeats and sites in which ' . 'identifiable features found in known proteins can be applied to new protein sequences. InterProScan is a tool ' . 'that combines different protein signature recognition methods.'),
        '#collapsible' => TRUE,
        '#collapsed' => TRUE,
        '#weight' => 0,
        'Active_site' => array(
          '#type' => 'fieldset',
          '#title' => t('Active site'),
          '#description' => t('The active site of an enzyme contains the catalytic and binding sites. The structure and chemical ' . 'properties of the active site allow the recognition and binding of the substrate. The active site is usually a big ' . 'pocket or cleft surrounded by amino acid- and other side chains at the surface of the enzyme that contains residues ' . 'responsible for the substrate specificity (charge, hydrophobicity, steric hindrance) and catalytic residues which ' . 'often act as proton donors or acceptors or are responsible for binding a cofactor.'
          ),
          '#collapsible' => TRUE,
          '#collapsed' => TRUE,
          '#weight' => 8,
          'ipr_active_site_select' => array(
            '#description' => t('Please select one or more InterProScan terms (using the Shift/Ctl/Option key). Use the Ctl/Option key to deselect.'),
            '#multiple' => '1',
            '#type' => 'select',
            '#options' => $arrayed_data['ipr']['active_site'],
            '#weight' => 1,
          ),
          ),
        'Binding_site' => array(
          '#type' => 'fieldset',
          '#title' => t('Binding site'),
          '#description' => t('A binding site is a region on a protein, DNA, or RNA to which ligands, specific ' . 'molecules and ions form a chemical bond.'
          ),
          '#collapsible' => TRUE,
          '#collapsed' => TRUE,
          '#weight' => 6,
          'ipr_binding_site_select' => array(
            '#description' => t('Please select one or more InterProScan terms (using the Shift/Ctl/Option key). Use the Ctl/Option key to deselect.'),
            '#multiple' => '1',
            '#type' => 'select',
            '#options' => $arrayed_data['ipr']['binding_site'],
            '#weight' => 1,
          ),
          ),
        'Conserved_site' => array(
          '#type' => 'fieldset',
          '#title' => t('Conserved site'),
          '#description' => t('Conserved sites are created from PROSITE patterns which are not a post-translational ' . 'modification or do not have a binding or catalytic activity but are conserved across members of a protein family.'),
          '#collapsible' => TRUE,
          '#collapsed' => TRUE,
          '#weight' => 9,
          'ipr_conserved_site_select' => array(
            '#description' => t('Please select one or more InterProScan terms (using the Shift/Ctl/Option key). Use the Ctl/Option key to deselect.'),
            '#multiple' => '1',
            '#type' => 'select',
            '#options' => $arrayed_data['ipr']['conserved_site'],
            '#weight' => 1,
          ),
          ),
        'Domain' => array(
          '#type' => 'fieldset',
          '#title' => t('Protein domain'),
          '#description' => t('A protein domain is a part of protein sequence and structure that can evolve, function, and exist ' . 'independently of the rest of the protein chain.'
          ),
          '#collapsible' => TRUE,
          '#collapsed' => TRUE,
          '#weight' => 3,
          'ipr_domain_select' => array(
            '#description' => t('Please select one or more InterProScan terms (using the Shift/Ctl/Option key). Use the Ctl/Option key to deselect.'),
            '#multiple' => '1',
            '#type' => 'select',
            '#options' => $arrayed_data['ipr']['domain'],
            '#weight' => 1,
          ),
          ),
        'Family' => array(
          '#type' => 'fieldset',
          '#title' => t('Protein Family'),
          '#description' => t('A protein family is a group of evolutionarily related proteins, and is often nearly synonymous with gene family.'),
          '#collapsible' => TRUE,
          '#collapsed' => TRUE,
          '#weight' => 1,
          'ipr_family_select' => array(
            '#description' => t('Please select one or more InterProScan terms (using the Shift/Ctl/Option key). Use the Ctl/Option key to deselect.'),
            '#multiple' => '1',
            '#type' => 'select',
            '#options' => $arrayed_data['ipr']['family'],
            '#weight' => 1,
          ),
          ),
        'PTM' => array(
          '#type' => 'fieldset',
          '#title' => t('Other post-translational modification'),
          '#description' => t('Posttranslational modification (PTM) is the chemical modification of a protein after its translation.'),
          '#collapsible' => TRUE,
          '#collapsed' => TRUE,
          '#weight' => 10,
          'ipr_ptm_select' => array(
            '#description' => t('Please select one or more InterProScan terms (using the Shift/Ctl/Option key). Use the Ctl/Option key to deselect.'),
            '#multiple' => '1',
            '#type' => 'select',
            '#options' => $arrayed_data['ipr']['ptm'],
            '#weight' => 1,
          ),
          ),
        'Region' => array(
          '#type' => 'fieldset',
          '#title' => t('Region'),
          '#description' => t('A protein region contains certain defining characteristics but unlike a domain, it cannot exist independently.'),
          '#collapsible' => TRUE,
          '#collapsed' => TRUE,
          '#weight' => 4,
          'ipr_region_select' => array(
            '#description' => t('Please select one or more InterProScan terms (using the Shift/Ctl/Option key). Use the Ctl/Option key to deselect.'),
            '#multiple' => '1',
            '#type' => 'select',
            '#options' => $arrayed_data['ipr']['region'],
            '#weight' => 1,
          ),
          ),
        'Repeat' => array(
          '#type' => 'fieldset',
          '#title' => t('Repeat'),
          '#description' => t('Conserved protein regions which form repeats with specific functions, such as the ankyrin repeat.'),
          '#collapsible' => TRUE,
          '#collapsed' => TRUE,
          '#weight' => 5,
          'ipr_repeat_select' => array(
            '#description' => t('Please select one or more InterProScan terms (using the Shift/Ctl/Option key). Use the Ctl/Option key to deselect.'),
            '#multiple' => '1',
            '#type' => 'select',
            '#options' => $arrayed_data['ipr']['repeat'],
            '#weight' => 1,
          ),
          ),
          ),
          );
          $form['filters']['blast'] = array(
      '#type' => 'fieldset',
      '#title' => t('Search for pre-annotation'),
      '#collapsible' => TRUE,
      '#collapsed' => TRUE,
      '#weight' => 1,
      'blast_search' => array(
        '#description' => t('Use this text box to search the BLAST descriptions.'),
        '#type' => 'textfield',
        '#title' => 'BLAST annotation',
        '#weight' => 1,
          ),
          );

          //add markers
          $form['filters']['Markers'] = array(
      '#type' => 'fieldset',
      '#title' => t('Markers'),
      '#description' => t('Markers are predicted using est2assembly and given unique permanent IDs.' . ' If you want, you may require one or more marker types in order to fetch' . ' contigs for a specific utility (e.g. neutral marker design).'
      ),
      '#collapsible' => TRUE,
      '#collapsed' => TRUE,
      '#weight' => 10,
      'SNP' => array(
        '#type' => 'fieldset',
        '#description' => t('Note that searching for the databased markers may be time consuming. Please be patient.'),
        '#collapsible' => TRUE,
        '#collapsed' => TRUE,
        '#title' => t('Single Nucleotide Polymorphisms'),
        'SNP_check' => array(
          '#type' => 'checkboxes',
      //'#title' => t('Single Nucleotide Polymorphisms'),
          '#description' => t('You may select one or more marker types to refine your query. You may select mutually exclusive terms (obviously, then you must use the "ANY selected criteria" retrieval method).'),
          '#options' => array(
            'force' => t('Any'),
            'synonymous' => t('Is coding & causes a synonymous change'),
            'non-synonymous' => t('Is coding & a non-synonymous change'),
            'transition' => t('Is a transition'),
            'transversion' => t('Is a transversion'),
            'coding' => t('Is part of a open-reading-frame'),
            'non-coding' => t('Is not part of a open-reading-frame'),
      ),
      ),
        'SNP_mapped' => array(
          '#type' => 'checkbox',
          '#title' => t('Mapped to genome'),
          '#disabled' => TRUE,
      ),
      ),
      'SSR' => array(
        '#type' => 'fieldset',
        '#collapsible' => TRUE,
        '#collapsed' => TRUE,
        '#title' => t('Simple (tandem) sequence repeats'),
        'SSR_check' => array(
          '#disabled' => TRUE,
          '#type' => 'checkboxes',
      //   '#title' => t('Simple (tandem) sequence repeats'),
          '#description' => t('You may select one or more marker types to refine your query. You may select mutually exclusive terms (obviously, then you must use the "ANY selected criteria" retrieval method).'),
          '#options' => array(
            'force' => t('Any'),
      ),
      ),
        'SSR_mapped' => array(
          '#type' => 'checkbox',
          '#title' => t('Mapped to genome'),
          '#disabled' => TRUE,
      ),
      ),
      );
      //add library data
      $form['filters']['Library'] = array(
      '#type' => 'fieldset',
      '#title' => t('Library data'),
      '#description' => t('If available, library data is extracted from dbEST and the Short Read Archive, trusting the user provided data and therefore not curated.' . ' If the submitter has provided useful information, such as developmental stage or tissue,' . ' you may use it to build more complex queries or refine your search.'
      ),
      '#collapsible' => TRUE,
      '#collapsed' => TRUE,
      '#weight' => 1,
      );
      if (!empty($arrayed_data['library'])) {
        foreach ($arrayed_data['library'] as $lib_name=>$tag_data){
          foreach ($tag_data as $tag => $values) {
            if (empty($tag) || $tag == 'species') {
              continue;
            }
            $tag_title = $tag;
            $lib_weight = 10;
            //add weight, skip and rename
            if ($tag == 'dev_stage') {
              $lib_weight = 0;
              $tag_title = 'Developmental stage';
            }
            elseif ($tag == 'sex') {
              $lib_weight = 1;
            }
            elseif ($tag == 'tissue_type') {
              $lib_weight = 2;
            }
            elseif ($tag == 'Organ') {
              $lib_weight = 3;
            }
            elseif ($tag == 'Site_1') {
              $tag_title = 'Restriction enzyme 1';
              $lib_weight = 19;
            }
            elseif ($tag == 'Site_2') {
              $tag_title = 'Restriction enzyme 2';
              $lib_weight = 20;
            }
            elseif ($tag == 'library name') {
              $lib_weight = -10;
            }
            elseif ($tag == 'map') {
              $tag_title = 'Cytogenetic map';
            }

            if ($tag != 'WARNING') {
              $tag = strtolower($tag);
              $tag = str_replace('_', ' ', $tag);
              $tag_title = ucwords($tag_title);
              $form['filters']['Library'][$tag] = array(
            '#type' => 'fieldset',
            '#title' => t($tag_title),
            '#collapsible' => TRUE,
            '#collapsed' => TRUE,
            '#weight' => $lib_weight,
            'lib_'. $tag .'_select' => array(
              '#tree' => TRUE,
              '#description' => t('You may select one or more terms (using the Shift/Ctl/Option key). Use the Ctl/Option key to deselect.'),
              '#multiple' => '1',
              '#type' => 'select',
              '#options' => $values,
              ),
              );
            }
          }}
      }
      $form['buttons'] = array(
        '#type'=>'fieldset',
      'all'=>array(
        '#value' => t(' or ALL features of type: '),
        '#attributes' => array('title' => t('Get all features (limited to the types you have selected) for this phylogeny selection. Other selections are ignored and this may result in a slow request.')),
        '#type' => 'submit',
      	'#weight' => 1,
      ),
      'data_types'=>array(
	    '#type' => 'checkboxes',
        '#description'=>t('Select one or more types to limit your search. If you select none, then all feature types will be searched.'),
      	'#default_value'=>TRUE,
	    '#options' => array(
          'contig' => t('Contig'),
          'polypeptide' => t('Peptide'),
          'ORF' => t('ORF'),
          'read' => t('EST, mRNA / non-NGS read'),
        	'SNP' => t('SNP'),
        ),
      	'#weight' => 2,
      ),
      );
      $form['clear'] = array(
        '#attributes' => array('title' => t('Clear the form')),
        '#type' => 'submit',
        '#value' => 'Reset this form',
        '#submit' => array('gmod_dbsf_form_clear'),
        '#weight' => 15,

      );
      $form['restart'] = array(
        '#attributes' => array('title' => t('Start from the beginning')),
        '#type' => 'submit',
        '#value' => 'Restart species selection',
        '#submit' => array('gmod_dbsf_form_clear_norebuild'),
        '#weight' => 15,
      );
      $form_state['values']['species_select']=$species_array;
      return $form;
}

/**
 * Get terms from result array, get feature_id linked to them, get anchillary information
 * return array which will tabulate results to genes4all_explore_results
 *
 * @param $form
 * @param $form_state
 *
 * @return unknown_type
 */
function genes4all_explore_search_species_form_submit($form, &$form_state) {
  //take $form_state['values'] and build a data array(cvterm,dbxref) which doesn't distinguish on their origin
  // passed to execute
  $data = array();
  // returned from execute
  $return_array = array();
  $jointype = '';
  if (!empty($form_state['values'])) {
    $data = array();
    $data['organism_array']=$form_state['values']['species_select'];
    if ($form_state['clicked_button']['#id'] == 'edit-intersect') {
      $data['jointype'] = 'intersect';
    }
    elseif ($form_state['clicked_button']['#id'] == 'edit-union') {
      $data['jointype']  = 'union';
    }
    // because no jointype means get all features:
    if (!empty($data['jointype'])){
    $data['dbxref_array'] = array_merge(
    $form_state['values']['ipr_active_site_select'],
    $form_state['values']['ipr_binding_site_select'],
    $form_state['values']['ipr_conserved_site_select'],
    $form_state['values']['ipr_domain_select'],
    $form_state['values']['ipr_family_select'],
    $form_state['values']['ipr_ptm_select'],
    $form_state['values']['ipr_region_select'],
    $form_state['values']['ipr_repeat_select']
    );
    $data['cvterm_array'] = array_merge(
    $form_state['values']['go_process_select'],
    $form_state['values']['go_function_select'],
    $form_state['values']['go_location_select'],
    $form_state['values']['kegg_select'],
    $form_state['values']['ec_select']
    );
    //sanitize blastsearch user input
    $data['blast_string'] = check_plain($form_state['values']['blast_search']);
    //database the library and markers
    //HACK HACK
    foreach ($form['#post'] as $key => $value) {
      if (!empty($value) && strpos($key, 'lib_') !== FALSE) {
        $newkey = str_replace('lib_', '', $key);
        $newkey = str_replace('_select', '', $newkey);
        $data['library'][$newkey] = $value;
      }
    }
    foreach ($form_state['values'] as $key => $value) {
      if (!empty($value) && ($key === 'SNP_check' || $key === 'SSR_check')) {
        $newkey = str_replace('_check', '', $key);
        //only include non-0 values
        foreach ($value as $type => $snp_result) {
          if ($snp_result !== 0) {
            $data['markers'][$newkey][$type] = $snp_result;
          }
        }
      }
    }
    }else{
      //only when jointype is empty
      $data['data_types']=$form_state['values']['data_types'];
    }
    if (!empty($data)){
      // we add $_SESSION['species'] to allow to pass links on which species is being select when the table
      // is presented. remember to unset when it is no longer needed.
      $_SESSION['species'] = implode(',',$form_state['values']['species_select']);
      //function url() doesnot work with $form_state['redirect'] (and neither giving the url directly)
      drupal_goto('genes4all/search/phylogeny',array('level'=>1,'bookmark'=>gmod_dbsf_encode($data)));
    }
  }
  else {
    drupal_set_message('No results found.','error');
    return;
  }

}

/**
 * Implementation of hook_form
 * This is the search form which allows to search by species IEA or blast annotation
 *
 * @param $form_state
 */
function genes4all_explore_select_species_form($form_state) {
  $phylogeny_array = array();
  $array = gmod_dbsf_db_execute('chado', 'gmod_dbsf_getspecies');
  if (empty($array)) {
    drupal_set_message('Not species found in the database!','error');
    return '';
  }
  $species_array   = ($array['species_array']);
  $phylogeny_array = ($array['phylogeny_array']);
  $button = array(
    'description'=>t('Limit your query to these organisms.'),
    'title'=>t('Add to filter & continue'),
    'validate'=>TRUE,
  );
  $form = gmod_dbsf_phylogeny_tree($phylogeny_array, $button);
  return $form;
}

/**
 * Implementation of hook_form_submit
 *
 * Redirects to genes4all_explore_page_phylosearch with a species GET array
 * composed of organism_id values of selected species
 */
function genes4all_explore_select_species_form_submit($form, &$form_state) {
  $data = array();
  if (!empty($form_state['values'])) {
    foreach ($form_state['values'] as $key => $value) {
      if (!empty($value) && !is_array($value) && preg_match('/^species-(\d+)/', $value, $matches)) {
        $data[$matches[1]] = $matches[1];
        //unset($form_state['values'][$key]);
      }
    }
    if (!empty($data)){
      drupal_goto('genes4all/search/phylogeny',array('species'=>implode(',',$data)));
    }else{
      drupal_set_message('No species selected!','error');
    }
  }else{
    drupal_set_message('No species selected!','error');
  }
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
  if (!isset($data)) {
    return FALSE;
  }
  if (empty($data['jointype'])){
    //we are getting all the features for this organism
    $sql = "SELECT feature_id from {feature} WHERE 1=1 ";
    $sql .= !empty($data['organism_array'])
      ? ' AND organism_id IN ('.implode(',', $data['organism_array']).') '
      :'';
      if (!empty($data['data_types'])){
        $data_to_add = array();
        foreach ($data['data_types'] as $name=>$is_set){
          if (!empty($is_set)){
            $data_to_add[] =$name;
          }
        }
        if (!empty($data_to_add)){
          $sql .= " AND type_id IN (SELECT cvterm_id from cvterm where cv_id=(SELECT cv_id from cv where name='sequence') AND name IN ("
          . gmod_dbsf_array2csv($data_to_add,TRUE)
          . '))';
        }
      }
    return $sql;
  }
  if ($data['jointype'] !== 'intersect' && $data['jointype'] !== 'union') {
    return FALSE;
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
  $species_query = '';
  if (!empty($data['organism_array'])) {
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
    //TODO TODO
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
    $message = t('It seems you haven\'t made any queries. Click back or ');
    if (!empty($_SESSION['species'])){
      $message .=l('try again','genes4all/search/phylogeny',array('query'=>array('species'=>$_SESSION['species']))).'.';
    }else{
      $message .=l('try again','genes4all/search/phylogeny').'.';
    }
    drupal_set_message($message,'error');
    drupal_goto('genes4all/noresults' . '/phylogeny');
  }
  //TODO Further TESTING
  $sql = "select feature_id from feature where feature_id IN (". $sql .")";
  if (!empty($organism_sql)) {
    $sql .= $organism_sql;
  }
  $sql .= ' order by type_id,length(feature.name),feature.name';
  if (!empty($sql)) {
    return $sql;
  }
}

/**
 * Run the SQL prepared from genes4all_explore_results_prepare_sql and calculate:
 * feature_ids to be tabulated
 * data is an array with keys:
 * dbxref_array
 * cvterm_array
 * blast_string
 * organism_array -> if empty then not used
 * jointype which expects intersect or union (lc)
 *  and causes a fetch feature_id for UNION (OR) or INTERSECT (AND) of the above (except organism which is always intersect).
 *  and whatever other data...
 *
 * @param $data is array(an sql,encoded of sql)
 *
 * @return unknown_type
 * @see http://php.net/manual/en/security.database.sql-injection.php
 */
function genes4all_explore_results_prepare_data($variables) {
  $is_chado=gmod_dbsf_is_chado();
  $offset = 0;
  // This is used for the links.
  $level = !empty($variables['level']) ? $variables['level'] : 0;
  // In all reality, the time required for the $sql query is not that much
  // influenced by the size of $limit, a limit of 10 and 50 require
  // almost the same amount of time.
  $limit = !empty($variables['limit']) ? $variables['limit'] : 10;

  if (!empty($level)) {
    $offset = ($level * $limit) - $limit;
  }
  // Grab one more from limit in order to decide if we are building a pager.
  $sql = $variables['sql'] ." OFFSET $offset LIMIT ". ($limit + 1);
  //$total_hits
  $hits_sql = str_replace('SELECT feature_id ','SELECT count(feature_id) ',$variables['sql']);
  $count = db_fetch_array(db_query($hits_sql));
  $count = $count['count'];
  $md5=md5($sql);
  if ($is_chado==TRUE){db_set_active();}
  $table_array = cache_get('SEARCH_'.$md5,'cache_genes4all_explore');
  if ($is_chado==TRUE){db_set_active('chado');}
  if (!empty($table_array)){
    return $table_array->data;
  }
  // Chado specific function, so Chado SQLs (and no {} fields)
  $result = db_query($sql);
  $feature_array = array();
  $org_array     = array();
  $link_array    = array();
  $feature_query = '';

  while ($dbdata = db_fetch_array($result)) {
    if (!isset($result_array[$dbdata['feature_id']]['exist'])) {
      $feature_query .= $dbdata['feature_id'] .",";
    }
  }
  unset($result);
  $feature_query = rtrim($feature_query, ",");
  if (empty($feature_query)) {
    //Nothing has been found, do not continue. It expects a table array for genes4all_explore_result_form so let's give a standard error page.
    //we redirect to genes4all/noresults (genes4all_explore_noresults)
    db_set_active('default');
    drupal_set_message('Sorry, no results where found. Try a new search with '
    .l('the same species','genes4all/search/phylogeny/',array('query'=>array('species'=>$_SESSION['species'])))
    .' or '
    .l('start anew','genes4all/search/phylogeny'),'error').'.';
    unset($_SESSION['species']);
    drupal_goto('genes4all/noresults' . '/phylogeny');
  }

  // continue with populating feature array
  // name, genus, species,common name
  // type_id, sequence
  //if (!isset($feature_query)||empty($feature_query)){$feature_query=int(0);}
  $sql_data1 = 'select feature.name as feature_name,organism_id,feature_id,feature.type_id,cvterm.name as type_name' . ' from feature,cvterm' . ' where feature.type_id=cvterm.cvterm_id and'. " feature_id IN ($feature_query) ";
  $sql_data1 .= ' order by type_id,length(feature.name),feature.name';

  // not the most clean implementation, but good enough
  // don't forget closing bracket
  $sql_data2_cont = 'select count(feature_id) as est_count from feature,cvterm where ' . "feature.type_id=cvterm.cvterm_id and cvterm.name='read'" . ' and feature_id IN (select feature_id from featureloc where srcfeature_id=';

  // don't forget closing bracket
  $sql_data2_pept = 'select name as prediction_method from analysis where analysis_id =(select analysis_id from ' . 'analysisfeature where feature_id=';

  // don't forget changing to Aorf adding _0 and closing ' and 2 brackets
  //TODO rename this sql variable to $sql_prot2contig
  // TODO create/use function to do this.
  $sql_data3_pept = 'select feature_id,name from feature where type_id = ' . "(select cvterm_id from cvterm where name = 'contig' and cv_id=(select cv_id from cv where name='sequence'))" . ' and feature_id IN (select srcfeature_id from featureloc where feature_id = ' . " (select feature_id from feature where uniquename ='";

  // this can be expanded with ncbi information which is why we are not using the existing function.
  // here, we are calling it once for all species so not calling it foreach feature
  $sql_data4      = 'select organism_id,genus,species,common_name from organism where organism_id IN (';
  $organism_query = ' ';
  $org_temp_array = array();
  // this can be expanded, we are calling it once per cvterm.type
  $sql_data5         = 'select cvterm_id,name from cvterm where cvterm_id IN (';
  $cvterm_query      = ' ';
  $cvterm_temp_array = array();
  $result            = db_query($sql_data1);
  $column_type ='';
  while ($dbdata = db_fetch_array($result)) {
    $type_specific = '';
    //TODO fix type data by parsing array (also in other arrays where similar implem.
    if ($dbdata['type_name'] == 'contig') {
      if (empty($column_type)){$column_type .=' No. ESTs';}
      $result2       = db_query($sql_data2_cont . $dbdata['feature_id'] .')');
      $res           = db_fetch_array($result2);
      $type_specific = $res['est_count'];
    }
    elseif ($dbdata['type_name'] == 'polypeptide') {
      if (empty($column_type)){$column_type .=' Prediction method';}
      $result2       = db_query($sql_data2_pept . $dbdata['feature_id'] .')');
      $res           = db_fetch_array($result2);
      $type_specific = $res['prediction_method'];
      $orf_id        = str_replace("Apep", "Aorf", $dbdata['feature_name']);
      $orf_id .= '_0';
      $result3 = db_query($sql_data3_pept . $orf_id ."'))");
      $res = db_fetch_array($result3);
      $link_array[$dbdata['feature_name']]['contig_name'] = $res['name'];
      $link_array[$dbdata['feature_name']]['contig_internalid'] = $res['feature_id'];
    }

    $org_temp_array[$dbdata['organism_id']] = 'exist';
    $cvterm_temp_array[$dbdata['type_id']] = 'exist';

    //$feature_array[$dbdata['feature_id']] = array(
    $feature_array[] = array(
      'feature.name' => $dbdata['feature_name'],
      'organism_id' => $dbdata['organism_id'],
      'type_id' => $dbdata['type_id'],
      'type_name' => $dbdata['type_name'],
      'type_specific' => $type_specific,
    );
  }
  if (empty($column_type)){$column_type .=' No. ESTs';}
  foreach ($org_temp_array as $key => $value) {
    $organism_query .= "$key,";
  }
  foreach ($cvterm_temp_array as $key => $value) {
    $cvterm_query .= "$key,";
  }
  $organism_query = rtrim($organism_query, ",");
  $cvterm_query = rtrim($cvterm_query, ",");
  //if (!isset($cvterm_query)||empty($cvterm_query)){$cvterm_query=int(0);}
  //if (!isset($organism_query)||empty($organism_query)){$organism_query=int(0);}
  $result = db_query($sql_data4 . $organism_query .')');
  while ($dbdata = db_fetch_array($result)) {
    $org_array[$dbdata['organism_id']] = array(
      'genus' => $dbdata['genus'],
      'species' => $dbdata['species'],
      'common_name' => $dbdata['common_name'],
    );
  }
  $result = db_query($sql_data5 . $cvterm_query .')');
  while ($dbdata = db_fetch_array($result)) {
    $cvterm_array[$dbdata['cvterm_id']]['name'] = $dbdata['name'];
  }


  $table_attributes = array('class' => 'padded-table');
  $table_header = array(
    'Name', 'Type', $column_type, 'Organism', 'Links',
  );


  //TODO option to order by other columns? currently only contig id (from links)
  $sorts      = array();
  $sorts1     = array();
  $sorts2     = array();
  $table_data = array();

  // Result count is the number of rows fetched with in the offset/limit.
  // The maximum we can have is actually limit+1, which if it exists then a pager is built.
  $result_count = 0;
  foreach ($feature_array as $feature_data) {
    $result_count++;
    $contig_id = !empty($link_array[$feature_data['feature.name']]['contig_name']) ? $link_array[$feature_data['feature.name']]['contig_name'] : $feature_data['feature.name'];
    $type_sorter = 999999999;
    $name_sorter = 9;
    $feature_data['organism'] = '<i>'. $org_array[$feature_data['organism_id']]['genus'] .' '. $org_array[$feature_data['organism_id']]['species'] .'</i><br>'. $org_array[$feature_data['organism_id']]['common_name'];
    $gbrowse_url = gmod_dbsf_get_add_var('genes4all_explore_gbrowse_link');
    if (!empty($gbrowse_url)) {
      preg_match('/^([A-Z]+)(\d+)/', $feature_data['feature.name'], $gbrowse_links);
      $feature_data['links'] = l(
      t('feature composition'), $gbrowse_url . $gbrowse_links[2] ."_caf/?name=". $contig_id,
      array(
          'absolute' => TRUE,
          'attributes' => array(
            'target' => '_blank',
            'title' => "See composition of $contig_id using GBrowse.",
      ),
      )
      );
    }
    elseif(module_exists('jbrowse')){
      if (in_array($feature_data['type_name'],array('contig','ORF','polypeptide'))){
        $feature_data['links'] =
        l(t('feature composition'),'jbrowse/'.$feature_data['feature.name'],array(
      'attributes' => array(
      	'target' => '_blank',
        'title' => "See composition of $contig_id using JBrowse.",
        ),
      'alias' => TRUE));
      }
    }
    else {
      $feature_data['links'] = '<em>Feature not activated</em>';
    }
    if ($feature_data['type_name'] == 'contig') {
      $type_sorter = $feature_data['type_specific'];
      $name_sorter = '1';
    }
    //TODO the following is very server memory intensive. what shall we do??
    $table_data[] = array(
    array('data' => $feature_data['feature.name']),
    array('data' => $feature_data['type_name']),
    array('data' => $feature_data['type_specific']),
    array('data' => $feature_data['organism']),
    array('data' => $feature_data['links']),
    );
    //$sorts[] = $contig_idnumber[0];
    //$sorts1[] = $gbrowse_links[1];
    //$sorts2[] = $gbrowse_links[2];
  }


  // Resetting to DRUPAL
  db_set_active('default');
  $dirname = drupal_get_path('module', 'genes4all');
  $next_link = array();
  if ($result_count > $limit) {
    // because we fetched $limit+1 records, we delete the last record.
    $result_count--;array_pop($table_data);
    $next_link =l(
    theme_image($dirname .'/images/next-button.png', 'Get next IDs', 'Get next IDs'), 'genes4all/search/phylogeny',
    array(
          'html' => TRUE,
          'alias' => TRUE,
          'attributes' => array('class' => 'genes4all_explore-center'),
          'query' => array('level' => $level + 1, 'bookmark' => $variables['encoded']),
    )
    );
  }
  // reset to string for display
  else {
    $next_link = '';
  }
  $bookmark = l(
  theme_image($dirname .'/images/bookmark_page.png', 'Bookmark', 'Bookmark'), 'genes4all/search/phylogeny',
  array(
      'html' => TRUE,
      'alias' => TRUE,
      'attributes' => array('class' => 'genes4all_explore-center'),
      'query' => array('bookmark' => $variables['encoded'], 'level' => $level),
  )
  );
  // only populated if limit is greater than one
  $previous_link = array();
  if ($level > 1) {
    $previous_link = l(
    theme_image($dirname .'/images/previous-button.png', 'Get previous IDs', 'Get previous IDs'), 'genes4all/search/phylogeny',
    array(
        'html' => TRUE,
        'alias' => TRUE,
        'attributes' => array('class' => 'genes4all_explore-center'),
        'query' => array('level' => $level - 1, 'bookmark' => $variables['encoded']),
    )
    );
  }
  // reset to string for display.
  else {
    $previous_link = '';
  }

  $table_caption = "<p>Showing $result_count out of $count results; sorted by Database ID";
  $table_caption .= ".</p>";
  $table_array = array(
    '#links' => array('previous' => $previous_link, 'next' => $next_link, 'bookmark' => $bookmark),
    'header' => $table_header,
    'caption' => $table_caption,
    'attributes' => $table_attributes,
    'data' => $table_data,
  );
  cache_set('SEARCH_'.$md5,$table_array , 'cache_genes4all_explore', time()+30*3600*24); //month
  return $table_array;
}

/**
 * Returns IEA data. If no species is provided, then it returns all species as well
 *
 * @param $species_array optional
 *
 * @return assoc. array with variable names as keys
 */
function genes4all_explore_getcvs($species_array=NULL) {
  global $user;
  $user_id=$user->uid;
  $cv_array = array();
  $return_array = array();
  $is_chado=gmod_dbsf_is_chado();
  $md5=md5(implode(',',(array)$species_array));
  if ($is_chado==TRUE){db_set_active();}
  $return_array = cache_get('CV_'.$user_id.'_'.$md5,'cache_genes4all_explore');
  if ($is_chado==TRUE){db_set_active('chado');}
  if (empty($return_array)){
    if (!empty($species_array)) {
      $species_query = gmod_dbsf_array2csv($species_array,TRUE);
      //if (!isset($species_query)||empty($species_query)){$species_query=int(0);}
      $sql = "select cvterm_name,cv_name,cvterm_id from cvterm_organism_mat where". " organism_id IN ($species_query)";
      $result = db_query($sql);
      while ($dbdata = db_fetch_array($result)) {
        // we dopn't want to use distinct:
        if (!isset($cv_array[$dbdata['cv_name']][$dbdata['cvterm_id']])) {
          $cv_array[$dbdata['cv_name']][$dbdata['cvterm_id']] = $dbdata['cvterm_name'];
        }
      }
      $return_array= array('cv_array' => $cv_array);
    }
    // if no limitation by species
    else {
      $sql             = "select cvterm_name,cv_name,cvterm_id from cvterm_organism_mat";
      $result          = db_query($sql);
      while ($dbdata = db_fetch_array($result)) {
        // we dopn't want to use distinct:
        if (!isset($cv_array[$dbdata['cv_name']][$dbdata['cvterm_id']])) {
          $cv_array[$dbdata['cv_name']][$dbdata['cvterm_id']] = $dbdata['cvterm_name'];
        }
      }
      $return_array=array('cv_array' => $cv_array);
    }
    if ($is_chado==TRUE){db_set_active();}
    cache_set('CV_'.$user_id.'_'.$md5, $return_array , 'cache_genes4all_explore');
    if ($is_chado==TRUE){db_set_active('chado');}
  }else{
    $return_array=$return_array->data;
  }
  return $return_array;
}

/**
 *
 * @param $species_array
 *
 * @return unknown_type
 */
function genes4all_explore_getdbxrefs($species_array) {
  $dbxref_array = array();
  $return_array = array();
  $is_chado=gmod_dbsf_is_chado();
  global $user;
  $user_id=$user->uid;
  $md5=md5(implode(',',(array)$species_array));
  if ($is_chado==TRUE){db_set_active();}
  $return_array = cache_get('DBXREF_'.$user_id.'_'.$md5,'cache_genes4all_explore');
  if ($is_chado==TRUE){db_set_active('chado');}
  if (empty($return_array)){

    if (!empty($species_array)) {
      $species_query = gmod_dbsf_array2csv($species_array,TRUE);
      //if (!isset($species_query)||empty($species_query)){$species_query=int(0);}
      $sql = "select dbxref_name,db_name,dbxref_id,dbxref_property from dbxref_organism_mat where". " organism_id IN ($species_query)";

      $result = db_query($sql);
      while ($dbdata = db_fetch_array($result)) {
        // we don't want to use distinct:
        if (!isset($dbxref_array[$dbdata['db_name']][$dbdata['dbxref_id']])) {
          $dbxref_array[$dbdata['db_name']][$dbdata['dbxref_property']][$dbdata['dbxref_id']] = $dbdata['dbxref_name'];
        }
      }
    }
    // if no limitation by species
    else {
      $sql = "select dbxref_name,db_name,dbxref_id,dbxref_property from dbxref_organism_mat";
      $result = db_query($sql);
      while ($dbdata = db_fetch_array($result)) {
        // we dopn't want to use distinct:
        if (!isset($dbxref_array[$dbdata['db_name']][$dbdata['dbxref_id']])) {
          $dbxref_array[$dbdata['db_name']][$dbdata['dbxref_property']][$dbdata['dbxref_id']] = $dbdata['dbxref_name'];
        }
      }
    }
    $return_array=array('dbxref_array' => $dbxref_array);
    if ($is_chado==TRUE){db_set_active();}
    cache_set('DBXREF_'.$user_id.'_'.$md5, $return_array , 'cache_genes4all_explore');
    if ($is_chado==TRUE){db_set_active('chado');}
  }else{
    $return_array=$return_array->data;
  }
  return $return_array;
}
