<?php
// $Id$

/**
 * @file
 * Installs blat
 */



/**
 * Helper function to add blat software
 *
 * Use it as a template for adding plugins
 *
 */
function biosoftware_bench_add_software_blat() {
  $software   = 'blat';
  $dbxref_id  = gmod_dbsf_get_add_dbxref_withnames('jim_kent@pacbell.net', $software);
  $check_sql  = "SELECT software_id as id from {gmod_dbsf_software} where uniquename='%s'";
  $cvterm_id  = gmod_dbsf_get_add_cv_withnames('resource_type', 'software option');
  $insert_sql = "INSERT INTO {gmod_dbsf_software} (uniquename,description) VALUES ('%s','%s')";
  $link_sql   = 'INSERT INTO {gmod_dbsf_software_resource} (software_id,resource_id) VALUES (' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='$software')" . ",(SELECT resource_id from {gmod_dbsf_resource} where uniquename='%s')" . ')';
  if (empty($dbxref_id) || empty($cvterm_id)) {
    drupal_set_message(t('There has been a problem getting the dbxref_id and/or cvterm_id for %software / "software option" from the database (!dbxref_id,!cvterm_id).', array('!dbxref_id' => $dbxref_id, '!cvterm_id' => $cvterm_id,'%software'=> $software)));
    return FALSE;
  }
  /*
-minMatch=N
-minScore=N
-minIdentity=N
-maxGap=N
-repMatch=N
-mask=type
-qMask=type
-repeats=type
-minRepDivergence=NN
-trimT
-noTrimA
-trimHardA
-fastMap
 -fine
-maxIntron=N
-extendThroughN
   mask types are
   lower - mask out lower cased sequence
   upper - mask out upper cased sequence
   file.out - mask database according to RepeatMasker file.out
   */
  $options_array = array(
    'blat_t' => array(
    	'description'=>t('Database type.'),
  		'options'=>array('DNA','Protein','DNA and translate in 6 frames.'),
  ),
    'blat_q' => array(
    'description'=>t('Query type.'),
  	'options'=>array('DNA','RNA','Protein','DNA and translate in 6 frames.','RNA and translate in 3 frames'),
  ),
    'blat_tileSize'=> array(
    'description'=>t('Search tile defines the size of the match that triggers an alignment. Usually between 8 and 12. If empty, default is 11 for DNA and 5 for protein.'),
  ),
  	'blat_stepSize'=> array(
'description'=>t('Spacing between tiles. If empty, default is equal to size of search tile, i.e. non-overlapping matches.'),
  ),
  	'blat_oneOff'=> array(
'description'=>t('Allow one mismatched residue in tile to still trigger an alignment.'),
  ),
  'blat_minMatch'=> array(
'description'=>t('Minimum number of tile matches to report alignment. Usually set from 2 to 4. If empty, default is 2 for nucleotide, 1 for protein.'),
  ),
  'blat_minScore'=> array(
'description'=>t('Minimum score to report alignment. Score is equal to matches minus the mismatches minus some sort of gap penalty. If empty, default is 30.'),
  ),
  'blat_minIdentity'=> array(
'description'=>t('Minimum % sequence identity. If empty, default is 90 for nucleotide searches, 25 for protein or translated protein searches.'),
  ),
  'blat_maxGap'=> array(
'description'=>t('Maximum gap size between tiles in an alignment. Usually between 0 to 3. If empty, default is 2.'),
  ),
  'blat_repMatch'=> array(
'description'=>t('Number of a tile may be repeated before it is marked as overused, i.e. a potential repeat. Typically this is 256 for tiles of size 12, 1024 for tile size 11, 4096 for tile size 10. If empty, default is 1024. Also affected by stepSize. When stepSize is halved repMatch is doubled to compensate.'),
  ),
  'blat_mask'=> array(
'description'=>t('Mask out repeats in reference sequence. Alignments will not be started in masked region but may extend through it in nucleotide searches. Masked areas are ignored entirely in protein or translated searches.'),
'options'=>array('lower-case','UPPER-case','RepeatMasked reference'),
  ),
  'blat_qMask'=> array(
'description'=>t('Mask out repeats in query sequence. Alignments will not be started in masked region but may extend through it in nucleotide searches. Masked areas are ignored entirely in protein or translated searches.'),
'options'=>array('lower-case','UPPER-case','RepeatMasked reference'),
  ),
  'blat_repeats'=> array(
'description'=>t('Repeat bases will not be masked in any way, but matches in repeat areas will be reported separately from matches in other areas in output.'),
'options'=>array('lower-case','UPPER-case','RepeatMasked reference'),
  ),
  // only with RepeatMasker
  'blat_minRepDivergence'=> array(
'description'=>t('Minimum percent divergence of repeats to allow them to be unmasked. If empty, default is 15. Only relevant for masking using RepeatMasker.'),
  ),
  // cDNA stuff:
  'blat_trimT'=> array(
'description'=>t('Trim leading poly-Ts.'),
  'input_type'=>'checkbox',
  ),
  // default is on and option is actually -noTrimA
  'blat_TrimA'=> array(
'description'=>t('Trim trailing poly-As.'),
  'input_type'=>'checkbox',
  ),
  'blat_trimHardA'=> array(
'description'=>t('Remove poly-A tail from qSize as well as output.'),
  'input_type'=>'checkbox',
  ),
  'blat_fastMap'=> array(
'description'=>t('Run for fast DNA/DNA remapping i.e. not allowing introns, requiring high percentage identity'),
  'input_type'=>'checkbox',
  ),
  'blat_fine' => array(
'description'=>t('For high quality mRNAs look harder for small initial and terminal exons. Not recommended for ESTs.'),
  'input_type'=>'checkbox',
  ),
  'blat_maxIntron'=> array(
'description'=>t('Sets maximum intron size. If empty, default is 750000.'),
  ),
  'blat_extendThroughN'=>array(
'description'=>t('Allows extension of alignment through large blocks of N\'s.'),
  ),
  );
  foreach ($options_array as $option => $data) {
    $name = !empty($data['name']) ? $data['name']: $option;
    biosoftware_bench_add_software_option($name, $dbxref_id, $cvterm_id, $data['description']);
    $input_type = 'textfield';
    if (!empty($data['options'])){
      gmod_dbsf_populate_prop('resource', $name, $cvterm_id,$data['options'] );
      $input_type = 'select';
    }else{
      $input_type = !empty($data['input_type']) ? $data['input_type'] : 'textfield';
    }
    gmod_dbsf_populate_cvterms('resource', $name, 'input type', array($input_type));
    //Now link them to the software
    db_query($link_sql, $name);
  }
  //Store allowed dataset suffixes. Rank 0 is the switch specifying if ANY or ALL of
  //the suffixes must be present for a dataset to be allowed
  //remember to include the "." if there is one before the suffix
  $suffix_cvterm_id = gmod_dbsf_get_add_cv_withnames('software_setting', 'dataset_suffix');
  $dataset_suffixes = gmod_dbsf_populate_prop('software', 'blat', $suffix_cvterm_id, array('ALL', '.repeatmasker.out'));
  return TRUE;
}

