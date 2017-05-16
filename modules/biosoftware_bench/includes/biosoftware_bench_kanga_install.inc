<?php
// $Id$

/**
 * @file
 * Installs kanga
 */



/**
 * Helper function to add kanga software
 *
 * Use it as a template for adding plugins
 *
 */
function biosoftware_bench_add_software_kanga() {
  $software   = 'kanga';
  $dbxref_id  = gmod_dbsf_get_add_dbxref_withnames('stuart.stephen@csiro.au', $software);
  $check_sql  = "SELECT software_id as id from {gmod_dbsf_software} where uniquename='%s'";
  $cvterm_id  = gmod_dbsf_get_add_cv_withnames('resource_type', 'software option');
  $insert_sql = "INSERT INTO {gmod_dbsf_software} (uniquename,description) VALUES ('%s','%s')";
  $link_sql   = 'INSERT INTO {gmod_dbsf_software_resource} (software_id,resource_id) VALUES (' . "(SELECT software_id from {gmod_dbsf_software} where uniquename='$software')" . ",(SELECT resource_id from {gmod_dbsf_resource} where uniquename='%s')" . ')';
  if (empty($dbxref_id) || empty($cvterm_id)) {
    drupal_set_message(t('There has been a problem getting the dbxref_id and/or cvterm_id for %software / "software option" from the database (!dbxref_id,!cvterm_id).', array('!dbxref_id' => $dbxref_id, '!cvterm_id' => $cvterm_id,'%software'=> $software)));
    return FALSE;
  }
  /*
   from version 1.2
   kanga:
   * -m, --mode=<int>          alignment processing mode: 0 - standard sensitivity, 1 - more sensitive (slower), 2 - ultra sensitive (slowest), 3 - less sensitive (quicker)
   *  -a, --microindellen=<int> accept microInDels inclusive of this length: 0 to 20 (default = 0 or no microIndels)
   *  -A, --splicejunctlen=<int> aligning RNA-seq, force flank trim, accept splice junctions separated by at most this distance: 25 to 200000 (default = 0 for DNA non-spliced aligning)
   *  -C, --colorspace          process for colorspace (SOLiD)
   *  -k, --nomultireads        do not accept multiple reads aligning to same loci (default is to accept all multiple reads)
   *  -q, --quality=<int>       fastq quality scoring - 0 - Sanger, 1 = Illumina 1.3+, 2 = Solexa < 1.3, 3 = Ignore quality (default = 3)
   *  -r, --mlmode=<int>        processing mode for reads mapped to multiple loci: 0 slough, 1 stats only, 2 rand, 3 cluster with uniques only, 4 cluster with uniques + other multi, 5 report all match loci up to '-R<limit>' (default is 0)
   *  -y, --trim5=<int>         trim this number of bases from 5' end of reads when loading raw reads (default is 0)
   *  -Y, --trim3=<int>         trim this number of bases from 3' end of reads when loading raw reads (default is 0)
   *  -R, --maxmulti=<int>      allow any read to match at most this many genome loci then process acording to mlmode (default is 10)
   *  -b, --bisulfite           process for bisulfite methylation patterning
   *  -e, --editdelta=<int>     accepted matches must be at least this Hamming edit distance from the next best match (default is 1, max 2)
   *  -s, --substitutions=<int> accept upto this number of aligner induced substitutions in alignments (default is 1/10th of individual read lengths)
   *  -n, --maxns=<int>         maximum number of indeterminate 'N's in reads before treating read as unalignable (default is 1, max 5)
   *  -x, --minflankexacts=<int> trim matching reads on 5' and 3' flanks until at least this number of exactly matching bases in flanks (default is 0 or no trimming)
   *  -p, --snpreadsmin=<int>   minimum read coverage at loci before processing for SNP determination (default is 0 or no SNP processing)
   *  -P, --snpreadsprob=<dbl>  only call a SNP if probability of false positive binominal(n,k,0.05) <= this (default is 0.001, range is 0-0.1)
   *  -t, --title=<string>      track title
   *  -U, --pemode=<int>        paired end processing mode: 0 - none, 1 - paired end with recover orphan ends, 2 - paired end no orphan recovery (default: 0)
   *  -d, --pairminlen=<int>    accept paired end alignments with apparent length of at least this (default = 100)
   *  -D, --pairmaxlen=<int>    accept paired end alignments with apparent length of at most this (default = 300)
   *  -u, --pair=<file>         if raw paired end processing then input read pairs from this raw paired end file
   -f, --FileLogLevel=<int>  Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug
   -F, --log=<file>          diagnostics log file
   -M, --format=<int>        output format: 0 - CSV loci only, 1 - CSV loci + match sequence, 2 - CSV loci + read sequence, 3 - CSV loci + read + match sequence, 4 - UCSC BED, 5 - SAM format with accepted aligned, 6 - SAM format with all reads (default: 0)
   -i, --in=<file>           input from these raw sequencer read files, wildcards allowed if single ended
   -I, --sfx=<file>          align against this suffix array (kangax generated) file
   -o, --out=<file>          output alignments to this file
   -j, --nonealign=<file>    output unalignable reads to this file (requires '-M3' output mode specified)
   -J, --multialign=<file>   output multialigned reads to this file (requires '-M3' output mode specified)
   -O, --stats=<file>        output aligner induced substitution distribution stats (requires '-M3' output mode specified) or paired end length distributions to this file
   -H, --hamming=<file>      Optional pre-generated Hamming edit distance file for targeted genome
   -Z, --chromexclude=<string> high priority - regular expressions defining chromosomes to exclude
   -z, --chromeinclude=<string> low priority - regular expressions defining chromosomes to include
   -T threads

   kangar:
   *   -m, --mode=<int>          processing mode: 0 - single end create, 1 - paired end create, 2 - output statistics, 3 - dump as fasta (default = 0)
   *  -k, --removedups          remove duplicate reads retaining only one (default is not to remove any duplicated reads)
   *  -q, --quality=<int>       fastq quality scoring - 0 - Sanger, 1 = Illumina 1.3+, 2 = Solexa < 1.3, 3 = Ignore quality (default = 3)
   *  -t, --trim5=<int>         trim this many bases off leading sequence 5' end (default = 0)
   *  -T, --trim3=<int>         trim this many bases off trailing sequence 3' end (default = 0)
   -i, --in=<file>           input from these raw sequencer read files, wildcards allowed if single ended
   -u, --pair=<file>         if paired end processing then input read pairs from this paired end file
   -o, --out=<file>          output accepted processed reads to this file
   -n, --numreadslimit=<int> limit number of reads (or dumps) in each input file to this many - 0 (default) if no limit
   -f, --FileLogLevel=<int>  Level of diagnostics written to screen and logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug
   -F, --log=<file>          diagnostics log file

   kangax:
   -f, --FileLogLevel=<int>  Level of diagnostics written to logfile 0=fatal,1=errors,2=info,3=diagnostics,4=debug
   -F, --log=<file>          diagnostics log file
   *  -r, --ref=<string>        reference species
   *  -d, --descr=<string>      full description
   *  -t, --title=<string>      short title
   *  -m, --mode=<int>          Processing mode, 0=standard, 1=bisulphite index
   *  -C, --colorspace          Generate for colorspace (SOLiD)
   -s, --blockseqlen=<int>   Generated suffix blocks to hold at most this length (MB) concatenated sequences (default = 3300)
   -i <file>                 input from wildcarded kangas or fasta files
   -o <file>                 output suffix array file
   -D, --TargDep             Generate target file only if missing or older than any of the independent source files
   */
  $options_array = array(
'kanga_m' => array(
'description'=>t('Alignment sensitivity.'),
'options'=>array(''),
  ),
'kanga_a' => array(
'description'=>t('Accept microInDels inclusive of this length.'),
  ),
'kanga_A' => array(
'description'=>t('Used when aligning RNA-seq: force flank trim & accept splice junctions separated by at most this distance: (default is 0; 25 - 200000).'),
  ),
'kanga_C' => array(
'description'=>t('Process for colorspace.'),
'input_type'=>'checkbox',
  ),
'kanga_k' => array(
'description'=>t('Do not accept multiple reads aligning to same loci.'),
'input_type'=>'checkbox',
  ),
'kanga_q' => array(
'description'=>t('Query fastq quality scoring.'),
'options'=>array(''),
  ),
'kanga_r' => array(
'description'=>t('Processing mode for reads mapped to multiple loci.'),
'options'=>array(''),
  ),
'kanga_y' => array(
'description'=>t('Trim this number of bases from 5\' end of reads when loading raw reads.'),
  ),
'kanga_Y' => array(
'description'=>t('Trim this number of bases from 3\' end of reads when loading raw reads.'),
  ),
'kanga_R' => array(

'description'=>t('Allow any read to match at most this many genome loci before processing reads matching at multiple locations.'),
  ),
'kanga_b' => array(
'description'=>t('Process for bisulfite methylation patterning.'),
'input_type'=>'checkbox',
  ),
'kanga_e' => array(
'description'=>t('Accepted matches must be at least this Hamming edit distance from the next best match.'),
'options'=>array(''),
  ),
'kanga_s' => array(
'description'=>t('Accept upto this number of aligner induced substitutions in alignments (default is 1/10th of individual read lengths).'),
  ),
'kanga_n' => array(
'description'=>t('Maximum number of indeterminate \'N\'s in reads before treating read as unalignable (default is 1, max 5).'),
  ),
'kanga_x' => array(
'description'=>t('Trim matching reads on 5\' and 3\' flanks until at least this number of exactly matching bases in flanks (default is 0 i.e. no trimming).'),
  ),
'kanga_p' => array(
'description'=>t('Minimum read coverage at loci before processing for SNP determination (default is 0 i.e. no SNP processing).'),
  ),
'kanga_P' => array(
'description'=>t('Only call a SNP if probability of false positive binominal => (n,k,0.05) (default is 0.001, range is 0-0.1).'),
  ),
'kanga_t' => array(
'description'=>t('Title.'),
  ),
'kanga_U' => array(
'description'=>t('Paired end processing mode.'),
'options'=>array(''),
  ),
'kanga_d' => array(
'description'=>t('Minimum Paired-End distance (default = 100).'),
  ),
'kanga_D' => array(
'description'=>t('Maximum Paired-End distance (default = 300).'),
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
    db_query($link_sql, $option);
  }
  $suffix_cvterm_id = gmod_dbsf_get_add_cv_withnames('software_setting', 'dataset_suffix');
  $dataset_suffixes = gmod_dbsf_populate_prop('software', 'kanga', $suffix_cvterm_id, array('ALL', '.kangax','.color.kangax'));
  return TRUE;
}

