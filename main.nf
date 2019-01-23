#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/ddamsproteomics
========================================================================================
 nf-core/ddamsproteomics Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/ddamsproteomics
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     nf-core/ddamsproteomics v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/ddamsproteomics --mzmls '*.mzML' --tdb swissprot_20181011.fa --mods assets/tmtmods.txt -profile standard,docker

    Mandatory arguments:
      --mzmls                       Path to mzML files
      --mzmldef                     Alternative to --mzml: path to file containing list of mzMLs 
                                    with sample set and fractionation annotation (see docs)
      --tdb                         Path to target FASTA protein database
      --mods                        Path to MSGF+ modification file (two examples in assets folder)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test

    Options:
      --isobaric VALUE              In case of isobaric, specify: tmt10plex, tmt6plex, itraq8plex, itraq4plex
      --activation VALUE            Specify activation protocol: hcd (DEFAULT), cid, etd for isobaric 
                                    quantification. Not necessary for other functionality.
      --normalize                   Normalize isobaric values by median centering on channels of protein table
      --genes                       Produce gene table (i.e. ENSG or gene names from Swissprot)
      --symbols                     Produce gene symbols table (i.e. gene names when using ENSEMBL DB)
      --martmap FILE                Necessary when using ENSEMBL FASTA database, tab-separated file 
                                    with information from Biomart. An example can be found at
                                    https://github.com/nf-core/test-datasets/raw/ddamsproteomics/testdata/
      --fractions                   Fractionated samples, 
      --hirief                      IEF fractionated samples, implies --fractions, allows delta pI calculation
      --pipep FILE                  File containing peptide sequences and their isoelectric points. Example
                                    can be found in https://github.com/nf-core/test-datasets/raw/ddamsproteomics/
      --onlypeptides                Do not produce protein or gene level data
      --noquant                     Do not produce isobaric or MS1 quantification data
      --quantlookup FILE            Use previously generated SQLite lookup database containing spectra 
                                    quantification data when e.g. re-running. Need to match exactly to the
                                    mzML files of the current run
      --fastadelim VALUE            FASTA header delimiter in case non-standard FASTA is used, to be used with
                                    --genefield
      --genefield VALUE             Number to determine in which field of the FASTA header (split 
                                    by --fastadelim) the gene name can be found.


    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.email = false
params.plaintext_email = false

params.martmap = false
params.isobaric = false
params.instrument = 'qe' // Default instrument is Q-Exactive
params.activation = 'hcd' // Only for isobaric quantification
params.outdir = 'results'
params.normalize = false
params.genes = false
params.symbols = false
params.fastadelim = false
params.genefield = false
params.quantlookup = false
params.fractions = false
params.hirief = false
params.pipep = false
params.onlypeptides = false
params.noquant = false

// Validate and set file inputs
fractionation = (params.hirief || params.fractions)
mods = file(params.mods)
if( !mods.exists() ) exit 1, "Modification file not found: ${params.mods}"
tdb = file(params.tdb)
if( !tdb.exists() ) exit 1, "Target fasta DB file not found: ${params.tdb}"

// Files which are not standard can be checked here
if (params.martmap) {
  martmap = file(params.martmap)
  if( !martmap.exists() ) exit 1, "Biomart ENSEMBL mapping file not found: ${params.martmap}"
}
if (params.pipep) {
  trainingpep = file(params.pipep)
  if( !trainingpep.exists() ) exit 1, "Peptide pI data file not found: ${params.pipep}"
}
output_docs = file("$baseDir/docs/output.md")

// set constant variables
accolmap = [peptides: 12, proteins: 14, genes: 17, assoc: 18]


// parse inputs that combine to form values or are otherwise more complex.
setdenoms = [:]
if (!(params.noquant) && params.isobaric) {
  params.denoms.tokenize(' ').each{ it -> x=it.tokenize(':'); setdenoms.put(x[0], x[1..-1])}
}
plextype = params.isobaric ? params.isobaric.replaceFirst(/[0-9]+plex/, "") : 'false'
normalize = (!params.noquant && params.normalize && params.isobaric)


// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Check workDir/outdir paths to be S3 buckets if running on AWSBatch
// related: https://github.com/nextflow-io/nextflow/issues/813
if( workflow.profile == 'awsbatch') {
    if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}


// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/ddamsproteomics v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/ddamsproteomics'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['mzMLs']        = params.mzmls
summary['Target DB']    = params.tdb
summary['Modifications'] = params.mods
summary['Instrument'] = params.instrument
summary['Isobaric tags'] = params.isobaric
summary['Isobaric activation'] = params.activation
summary['Isobaric median normalization'] = params.normalize
summary['Output genes'] = params.genes
summary['Output symbols'] = params.symbols
summary['Custom FASTA delimiter'] = params.fastadelim 
summary['Custom FASTA gene field'] = params.genefield
summary['Premade quant data SQLite'] = params.quantlookup
summary['Fractionated sample'] = fractionation
summary['HiRIEF'] = params.hirief 
summary['peptide pI data'] = params.pipep
summary['Only output peptides'] = params.onlypeptides
summary['Do not quantify'] = params.noquant
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-ddamsproteomics-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/ddamsproteomics Workflow Summary'
    section_href: 'https://github.com/nf-core/ddamsproteomics'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    msgf_plus | head -n1 > v_msgf.txt
    hardklor | head -n1 > v_hk.txt || true
    kronik | head -n2 > v_kr.txt
    percolator -h |& head -n1 > v_perco.txt || true
    msspsmtable --version > v_mss.txt
    source activate openms-2.4.0
    IsobaricAnalyzer |& grep Version > v_openms.txt || true
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

if (params.mzmlPaths) {
  // Profile 'test' delivers mzmlPaths
  Channel
    .from(params.mzmlPaths)
    .set { mzml_in }
}
else if (!params.mzmldef) {
  Channel
    .fromPath(params.mzmls)
    .map { it -> [it, 'NA'] }
    .set { mzml_in }
} else {
  Channel
    .from(file("${params.mzmldef}").readLines())
    .map { it -> it.tokenize('\t') }
    .set { mzml_in }
}


// Parse mzML input to get files and sample names etc
// get setname, sample name (baseName), input mzML file. 
// Set platename to samplename if not specified. 
// Set fraction name to NA if not specified
mzml_in
  .map { it -> [it[1], file(it[0]).baseName.replaceFirst(/.*\/(\S+)\.mzML/, "\$1"), file(it[0]), it[2] ? it[2] : it[1], it[3] ? it[3] : 'NA' ]}
  .tap{ sets; strips; mzmlfiles; mzml_isobaric; mzml_quant; mzml_msgf }
  .count()
  .set{ amount_mzml }

// Set names are first item in input lists, collect them for PSM tables and QC purposes
sets
  .map{ it -> it[0] }
  .unique()
  .tap { setnames_psm } 
  .collect()
  .map { it -> [it] }
  .into { setnames_featqc; setnames_psmqc }

// Strip names for HiRIEF fractionation are third item, 
strips
  .map { it -> it[3] }
  .unique()
  .toList()
  .set { strips_for_deltapi }


/*
* Step 1: Extract quant data from peptide spectra
*/

process quantifySpectra {
  when: !params.quantlookup && !params.noquant

  input:
  set val(setname), val(sample), file(infile), val(platename), val(fraction) from mzml_quant
  file(hkconf) from Channel.fromPath("$baseDir/assets/hardklor.conf").first()

  output:
  set val(sample), file("${sample}.kr"), file(infile) into kronik_out
  set val(sample), file("${infile}.consensusXML") optional true into isobaricxml

  script:
  activationtype = [hcd:'High-energy collision-induced dissociation', cid:'Collision-induced dissociation', etd:'Electron transfer dissociation'][params.activation]
  massshift = [tmt:0.0013, itraq:0.00125, false:0][plextype]
  """
  # Run hardklor on config file with added line for in/out files
  # then run kronik on hardklor and quant isobaric labels if necessary
  hardklor <(cat $hkconf <(echo "$infile" hardklor.out))
  kronik -c 5 -d 3 -g 1 -m 8000 -n 600 -p 10 hardklor.out ${sample}.kr
  source activate openms-2.4.0
  ${params.isobaric ? "IsobaricAnalyzer  -type $params.isobaric -in $infile -out \"${infile}.consensusXML\" -extraction:select_activation \"$activationtype\" -extraction:reporter_mass_shift $massshift -extraction:min_precursor_intensity 1.0 -extraction:keep_unannotated_precursor true -quantification:isotope_correction true" : ''}
  """
}


// Collect all mzMLs into single item to pass to lookup builder and spectra counter
mzmlfiles
  .buffer(size: amount_mzml.value)
  .map { it.sort( {a, b -> a[1] <=> b[1]}) } // sort on sample for consistent .sh script in -resume
  .map { it -> [it.collect() { it[0] }, it.collect() { it[2] }, it.collect() { it[3] } ] } // lists: [sets], [mzmlfiles], [plates]
  .into { mzmlfiles_all; mzmlfiles_all_count }


process createSpectraLookup {

  when: !params.quantlookup

  input:
  set val(setnames), file(mzmlfiles), val(platenames) from mzmlfiles_all

  output:
  file 'mslookup_db.sqlite' into newspeclookup 

  script:
  """
  msslookup spectra -i ${mzmlfiles.join(' ')} --setnames ${setnames.join(' ')}
  """
}


// Collect all isobaric quant XML output for quant lookup building process
isobaricxml
  .ifEmpty(['NA', 'NA', 'NA'])
  .buffer(size: params.isobaric ? amount_mzml.value : 1)
  .map { it.sort({a, b -> a[0] <=> b[0]}) }
  .map { it -> [it.collect() { it[0] }, it.collect() { it[1] }] } // samples, isoxml
  .set { isofiles_sets }

// Collect all MS1 kronik output for quant lookup building process
kronik_out
  .ifEmpty(['NA', 'NA'])
  .buffer(size: amount_mzml.value)
  .map { it.sort({a, b -> a[0] <=> b[0]}) }
  .map { it -> [it.collect() { it[0] }, it.collect() { it[1] }, it.collect() { it[2] }] } // samples, kronikout, mzml
  .set { krfiles_sets }


// Need to populate channels depending on if a pre-made quant lookup has been passed
// even if not needing quant (--noquant) this is necessary or NF will error
if (params.noquant && !params.quantlookup) {
  newspeclookup
    .into { quant_lookup; spec_lookup; countlookup }
} else if (!params.quantlookup) {
  newspeclookup
    .into { spec_lookup; countlookup }
} else {
  Channel
    .fromPath(params.quantlookup)
    .into { quant_lookup; countlookup }
} 


process quantLookup {

  when: !params.quantlookup && !params.noquant

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {it == 'db.sqlite' ? 'quant_lookup.sql' : null }

  input:
  file lookup from spec_lookup
  set val(isosamples), file(isofns) from isofiles_sets
  set val(krsamples), file(krfns), file(mzmls) from krfiles_sets

  output:
  file 'db.sqlite' into newquantlookup

  script:
  if (params.isobaric)
  """
  # SQLite lookup needs copying to not modify the input file which would mess up a rerun with -resume
  cat $lookup > db.sqlite
  msslookup ms1quant --dbfile db.sqlite -i ${krfns.join(' ')} --spectra ${mzmls.join(' ')} --quanttype kronik --mztol 20.0 --mztoltype ppm --rttol 5.0 
  msslookup isoquant --dbfile db.sqlite -i ${isofns.join(' ')} --spectra ${isosamples.collect{ x -> x + '.mzML' }.join(' ')}
  """
  else
  """
  # SQLite lookup needs copying to not modify the input file which would mess up a rerun with -resume
  cat $lookup > db.sqlite
  msslookup ms1quant --dbfile db.sqlite -i ${krfns.join(' ')} --spectra ${mzmls.join(' ')} --quanttype kronik --mztol 20.0 --mztoltype ppm --rttol 5.0 
  """
}


if (!params.quantlookup && !params.noquant) {
  newquantlookup
    .set { quant_lookup }
} 

mzmlfiles_all_count
  .merge(countlookup)
  .set { specfilein }


process countMS2perFile {

  input:
  set val(setnames), file(mzmlfiles), val(platenames), file(speclookup) from specfilein

  output:
  set val(setnames), file(mzmlfiles), val(platenames), file('amount_spectra_files') into specfilems2

  script:
  """
  sqlite3 $speclookup "SELECT mzmlfilename, COUNT(*) FROM mzml JOIN mzmlfiles USING(mzmlfile_id) JOIN biosets USING(set_id) GROUP BY mzmlfilename" > amount_spectra_files
  """
}


if (fractionation) { 
  specfilems2.set { scans_platecount }
} else {
  specfilems2
    .map { it -> [it[3], ['noplates']] }
    .into { scans_platecount; scans_result }
}


process countMS2sPerPlate {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true 
  when: fractionation

  input:
  set val(setnames), file(mzmlfiles), val(platenames), file('nr_spec_per_file') from scans_platecount

  output:
  set file('scans_per_plate'), val(splates) into scans_perplate

  script:
  splates = [setnames, platenames].transpose().collect() { "${it[0]}_${it[1]}" }
  """
  #!/usr/bin/env python
  platesets = [\"${splates.join('", "')}\"]
  platescans = {p: 0 for p in platesets}
  fileplates = {fn: p for fn, p in zip([\"${mzmlfiles.join('", "')}\"], platesets)}
  with open('nr_spec_per_file') as fp:
      for line in fp:
          fn, scans = line.strip('\\n').split('|')
          platescans[fileplates[fn]] += int(scans)
  with open('scans_per_plate', 'w') as fp:
      for plate, scans in platescans.items():
          fp.write('{}\\t{}\\n'.format(plate, scans))
  """
}

if (fractionation) {
  scans_perplate.set { scans_result }
}

/*
* Step 2: Identify peptides
*/

process createTargetDecoyFasta {
 
  input:
  file(tdb)

  output:
  file('db.fa') into concatdb
  set file(tdb), file("decoy_${tdb}") into searchdbs 

  script:
  """
  tryprev.py $tdb
  cat $tdb decoy_${tdb} > db.fa
  """
}


process msgfPlus {

  input:
  set val(setname), val(sample), file(x), val(platename), val(fraction) from mzml_msgf
  file(db) from concatdb
  file mods

  output:
  set val(setname), val(sample), file("${sample}.mzid") into mzids
  set val(setname), file("${sample}.mzid"), file('out.mzid.tsv'), val(platename), val(fraction) into mzidtsvs
  
  script:
  msgfprotocol = [tmt:4, itraq:2, false:0][plextype]
  msgfinstrument = [velos:1, qe:3, false:0][params.instrument]
  """
  msgf_plus -Xmx16G -d $db -s $x -o "${sample}.mzid" -thread 12 -mod $mods -tda 0 -t 10.0ppm -ti -1,2 -m 0 -inst ${msgfinstrument} -e 1 -protocol ${msgfprotocol} -ntt 2 -minLength 7 -maxLength 50 -minCharge 2 -maxCharge 6 -n 1 -addFeatures 1
  msgf_plus -Xmx3500M edu.ucsd.msjava.ui.MzIDToTsv -i "${sample}.mzid" -o out.mzid.tsv
  rm ${db.baseName.replaceFirst(/\.fasta/, "")}.c*
  """
}

mzids
  .groupTuple()
  .set { mzids_2pin }


process percolator {

  input:
  set val(setname), val(samples), file('mzid?') from mzids_2pin

  output:
  set val(setname), file('perco.xml') into percolated

  """
  echo $samples
  mkdir mzids
  count=1;for sam in ${samples.join(' ')}; do ln -s `pwd`/mzid\$count mzids/\${sam}.mzid; echo mzids/\${sam}.mzid >> metafile; ((count++));done
  msgf2pin -o percoin.xml -e trypsin -P "decoy_" metafile
  percolator -j percoin.xml -X perco.xml -N 500000 --decoy-xml-output -y
  """
}


mzidtsvs
  .groupTuple()
  .join(percolated)
  .set { mzperco }


process svmToTSV {

  input:
  set val(setname), file('mzident????'), file('mzidtsv????'), val(platenames), val(fractions), file(perco) from mzperco 

  output:
  set val(setname), val('target'), file('tmzidperco') into tmzidtsv_perco
  set val(setname), val('decoy'), file('dmzidperco') into dmzidtsv_perco

  script:
  """
  perco_to_tsv.py -p $perco --plates ${platenames.join(' ')} --fractions ${fractions.join(' ')}
  """
}

// Collect percolator data of target/decoy and feed into PSM table creation
tmzidtsv_perco
  .concat(dmzidtsv_perco)
  .groupTuple(by: 1)
  .combine(quant_lookup)
  .set { prepsm }

// Set strips to false if not running hirief
if (params.hirief) {
  strips_for_deltapi
    .map { it -> [it, trainingpep] }
    .set { stripannot }
} else {
  strips_for_deltapi
    .map { it -> [it, false, false]}
    .set { stripannot }
}


/*
* Step 3: Post-process peptide identification data
*/

process createPSMTable {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {["target_psmlookup.sql", "target_psmtable.txt", "decoy_psmtable.txt"].contains(it) ? it : null}

  input:
  set val(setnames), val(td), file('psms?'), file('lookup') from prepsm
  set file(tdb), file(ddb) from searchdbs
  set val(allstrips), file(trainingpep) from stripannot

  output:
  set val(td), file("${outpsms}") into psm_result
  set val(td), file({setnames.collect() { it + '.tsv' }}) into setpsmtables
  set val(td), file("${psmlookup}") into psmlookup

  script:
  psmlookup = "${td}_psmlookup.sql"
  outpsms = "${td}_psmtable.txt"
  """
  msspsmtable merge -i psms* -o psms.txt
  msspsmtable conffilt -i psms.txt -o filtpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PSM q-value'
  msspsmtable conffilt -i filtpsm -o filtpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'peptide q-value'
  # SQLite lookup needs copying to not modify the input file which would mess up a rerun with -resume
  cat lookup > $psmlookup
  msslookup psms -i filtpep --dbfile $psmlookup ${params.onlypeptides ? '' : "--fasta ${td == 'target' ? tdb : "${ddb} --decoy"}"} ${params.martmap ? "--map ${martmap}" : ''}
  msspsmtable specdata -i filtpep --dbfile $psmlookup -o prepsms.txt
  ${!params.noquant ? "msspsmtable quant -i prepsms.txt -o qpsms.txt --dbfile $psmlookup --precursor ${params.isobaric && td=='target' ? '--isobaric' : ''}" : 'mv prepsms.txt qpsms.txt'}
  sed 's/\\#SpecFile/SpectraFile/' -i qpsms.txt
  ${!params.onlypeptides ? "msspsmtable genes -i qpsms.txt -o gpsms --dbfile $psmlookup" : ''}
  ${!params.onlypeptides ? "msslookup proteingroup -i qpsms.txt --dbfile $psmlookup" : ''}
  ${!params.onlypeptides ? "msspsmtable proteingroup -i gpsms -o ${params.hirief ? "pgpsms" : "$outpsms"} --dbfile $psmlookup" : 'mv qpsms.txt pgpsms'}
  ${params.hirief ? "peptide_pi_annotator.py -i $trainingpep -p pgpsms --o $outpsms --stripcolpattern Strip --pepcolpattern Peptide --fraccolpattern Fraction --strippatterns ${allstrips.join(' ')} --intercepts ${allstrips.collect() { params.strips[it].intercept}.join(' ')} --widths ${allstrips.collect() { params.strips[it].fr_width}.join(' ')} --ignoremods \'*\'" : ''}
  msspsmtable split -i ${outpsms} --bioset
  """
}

// Collect setnames and merge with PSM tables for peptide table creation
setnames_psm
  .toList()
  .map { it -> [it.sort()]}
  .set { setlist_psm }
setpsmtables
  .map { it -> [it[0], it[1] instanceof java.util.List ? it[1] : [it[1]] ] }
  .map{ it -> [it[0], it[1].sort { a, b -> a.baseName.tokenize('.')[0] <=> b.baseName.tokenize('.')[0] }] } // names are setnames, sort on them then merge with sorted setnames
  .merge(setlist_psm)
  .transpose()
  .set { psm_pep }


process psm2Peptides {

  input:
  set val(td), file('psms'), val(setname) from psm_pep
  
  output:
  set val(setname), val(td), file("${setname}_linmod"), file('proteinratios') into pepslinmod
  set val(setname), val('peptides'), val(td), file("${setname}_linmod") into peptides_out
  set val(setname), file('normratiosused') optional true into normratios
  set val(setname), val(td), file(psms), file('proteins'), val('proteins') into proteins
  set val(setname), val(td), file(psms), file('genes'), val('genes') into genes
  set val(setname), val(td), file(psms), file('symbols'), val('assoc') into symbols
  set val(setname), file('proteinratios') optional true into proteinratios

  script:
  col = accolmap.peptides + 1  // psm2pep adds a column
  isoquant = !params.noquant && params.isobaric && td == 'target'
  """
  # Create peptide table from PSM table, picking best scoring unique peptides
  msspeptable psm2pep -i psms -o peptides --scorecolpattern svm --spectracol 1 ${!params.noquant && params.isobaric && td == 'target' ? "--isobquantcolpattern plex" : "" } ${!params.noquant ? "--ms1quantcolpattern area" : ""}
  # Move peptide sequence to first column
  paste <( cut -f ${col} peptides) <( cut -f 1-${col-1},${col+1}-500 peptides) > peptide_table.txt
  # Create empty protein/gene/gene-symbol tables with only the identified accessions, will be filled later
  echo Protein accession |tee proteins genes symbols
  tail -n+2 psms|cut -f ${accolmap.proteins}|grep -v '\\;'|grep -v "^\$"|sort|uniq >> proteins
  tail -n+2 psms|cut -f ${accolmap.genes}|grep -v '\\;'|grep -v "^\$"|sort|uniq >> genes
  tail -n+2 psms|cut -f ${accolmap.assoc}|grep -v '\\;'|grep -v "^\$"|sort|uniq >> symbols
  # Do isobaric quantification if necessary
  ${normalize && td == 'target' ? "msspsmtable isoratio -i psms -o proteinratios --protcol ${accolmap.proteins} --targettable proteins --isobquantcolpattern plex --minint 0.1 --denompatterns ${setdenoms[setname].join(' ')}" : 'touch proteinratios'}
  ${isoquant ? "msspsmtable isoratio -i psms -o pepisoquant --targettable peptide_table.txt --protcol ${accolmap.peptides} --isobquantcolpattern plex --minint 0.1 --denompatterns ${setdenoms[setname].join(' ')} ${normalize ? '--normalize median --norm-ratios proteinratios' : ''} > normratiosused" : ''}
  ${isoquant ? "mv pepisoquant peptide_table.txt" : ''}
  # Create linear modeled q-values of peptides (modeled svm scores vs q-values) for more protein-FDR precision.
  msspeptable modelqvals -i peptide_table.txt -o ${setname}_linmod --scorecolpattern svm --fdrcolpattern '^q-value'
  """
}


// Different amount of processes depending on genes and gene symbols are desired
// Input for proteins, genes and symbols is identical at this stage so tap and concat
// onto itself.
if (params.genes && params.symbols) { 
  pepslinmod
    .tap { pepsg; pepss }
    .concat(pepsg, pepss)
    .set { pepslinmod_prot }
  proteins
    .concat(genes, symbols)
    .join(pepslinmod_prot, by: [0,1])
    .set { prepgs_in }
} else if (params.genes) { 
  pepslinmod
    .tap { pepsg }
    .concat(pepsg)
    .set { pepslinmod_prot }
  proteins
    .concat(genes)
    .join(pepslinmod_prot, by: [0,1])
    .set { prepgs_in }
} else { 
  proteins
    .join(pepslinmod, by: [0,1])
    .set { prepgs_in }
}


/*
* Step 4: Infer and quantify proteins and genes
*/

process prepProteinGeneSymbolTable {

  when: !params.onlypeptides

  input:
  set val(setname), val(td), file('psms'), file('proteins'), val(acctype), file('peplinmod'), file('pratios') from prepgs_in

  output:
  set val(setname), val(acctype), val(td), file('bestpeptides') into bestpep

  script:
  if (!params.noquant && params.isobaric && td == 'target')
  """
  mssprottable ms1quant -i proteins -o protms1 --psmtable psms --protcol ${accolmap[acctype]}
  msspsmtable isoratio -i psms -o proteintable --protcol ${accolmap[acctype]} --targettable protms1 --isobquantcolpattern plex --minint 0.1 --denompatterns ${setdenoms[setname].join(' ')} ${normalize && td == 'target' ? '--norm-ratios pratios --normalize median': ''}
  mssprottable bestpeptide -i proteintable -o bestpeptides --peptable peplinmod --scorecolpattern ${acctype == 'proteins' ? '\'^q-value\'' : '\'linear model\''} --logscore --protcol ${accolmap[acctype] + 1}
  """
  else
  """
  ${td == 'target' && !params.noquant ? "mssprottable ms1quant -i proteins -o proteintable --psmtable psms --protcol ${accolmap[acctype]}" : 'mv proteins proteintable'}
  mssprottable bestpeptide -i proteintable -o bestpeptides --peptable peplinmod --scorecolpattern ${acctype == 'proteins' ? '\'^q-value\'' : '\'linear model\''} --logscore --protcol ${accolmap[acctype] + 1}
  """
}

tbestpep = Channel.create()
dbestpep = Channel.create()
bestpep
  .groupTuple(by: [0,1])
  .transpose()
  .choice(tbestpep, dbestpep) { it[2] == 'target' ? 0 : 1 }


process proteinFDR {
  
  when: !params.onlypeptides
  input:
  set val(setname), val(acctype), val(td), file('tbestpep') from tbestpep
  set val(setname), val(acctype), val(td), file('dbestpep') from dbestpep
  set file(tfasta), file(dfasta) from searchdbs

  output:
  set val(setname), val(acctype), file("${setname}_protfdr") into protfdrout
  script:
  if (acctype == 'genes')
  """
  mssprottable pickedfdr --picktype fasta --targetfasta $tfasta --decoyfasta $dfasta ${params.fastadelim ? "--fastadelim \'${params.fastadelim}\' --genefield ${params.genefield}" : ''} -i tbestpep --decoyfn dbestpep -o ${setname}_protfdr
  """
  else
  """
  mssprottable ${acctype == 'proteins' ? 'protfdr' : 'pickedfdr --picktype result'} -i tbestpep --decoyfn dbestpep -o ${setname}_protfdr
  """
}

// Prepare for set-tables merging into a side-by-side table
peptides_out
  .filter { it[2] == 'target' }
  .map { it -> [it[0], it[1], it[3]] }
  .set { peptides_to_merge }

if (!params.onlypeptides) {
  peptides_to_merge
    .concat(protfdrout)
    .groupTuple(by: 1)
    .set { ptables_to_merge }
} else {
  peptides_to_merge
    .groupTuple(by: 1)
    .set { ptables_to_merge }
}

psmlookup
  .filter { it[0] == 'target' }
  .collect()
  .map { it[1] }
  .set { tlookup }

/*
* Step 5: Create reports
*/

process proteinPeptideSetMerge {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "proteintable" ? "${outname}_table.txt": null}

  input:
  set val(setnames), val(acctype), file(tables) from ptables_to_merge
  file(lookup) from tlookup
  
  output:
  set val(acctype), file('proteintable') into featuretables
  set val(acctype), file('proteintable') into featqc_getpeptable

  script:
  outname = (acctype == 'assoc') ? 'symbols' : acctype
  """
  # SQLite lookup needs copying to not modify the input file which would mess up a rerun with -resume
  cat $lookup > db.sqlite
  msslookup ${acctype == 'peptides' ? 'peptides --fdrcolpattern \'^q-value\' --peptidecol' : 'proteins --fdrcolpattern \'q-value\' --protcol'} 1 --dbfile db.sqlite -i ${tables.join(' ')} --setnames ${setnames.join(' ')} ${!params.noquant ? "--ms1quantcolpattern area" : ""}  ${!params.noquant && params.isobaric ? '--psmnrcolpattern quanted --isobquantcolpattern plex' : ''} ${acctype in ['genes', 'assoc'] ? "--genecentric ${acctype}" : ''}
  ${acctype == 'peptides' ? 'msspeptable build' : 'mssprottable build --mergecutoff 0.01'} --dbfile db.sqlite -o proteintable ${!params.noquant && params.isobaric ? '--isobaric' : ''} ${!params.noquant ? "--precursor": ""} --fdr ${acctype in ['genes', 'assoc'] ? "--genecentric ${acctype}" : ''} ${params.onlypeptides ? "--noncentric" : ''}
  sed -i 's/\\#/Amount/g' proteintable
  """
}

psm_result
  .filter { it[0] == 'target' }
  .merge(scans_result)
  .map { it -> [it[0], it[1], it[2], it[3].unique()] }
  .set { targetpsm_result }


process psmQC {
  input:
  set val(td), file('psms'), file('scans'), val(plates) from targetpsm_result
  val(setnames) from setnames_psmqc
  output:
  set val('psms'), file('psmqc.html') into psmqccollect
  val(plates) into qcplates
  // TODO no proteins == no coverage for pep centric
  script:
  """
  qc_psms.R ${setnames[0].size()} ${fractionation ? 'TRUE' : 'FALSE'} ${plates.join(' ')}
  echo "<html><body>" > psmqc.html
  for graph in psm-scans missing-tmt miscleav
    do
    [[ -e \$graph ]] && paste -d \\\\0  <(echo "<div class=\\"chunk\\" id=\\"\${graph}\\"><img src=\\"data:image/png;base64,") <(base64 -w 0 \$graph) <(echo '"></div>') >> psmqc.html
    done 
  for graph in retentiontime precerror fryield msgfscore
    do
    for plateid in ${plates.join(' ')}
      do
      plate="PLATE___\${plateid}___\${graph}"
      [[ -e \$plate ]] && paste -d \\\\0  <(echo "<div class=\\"chunk \$plateid\\" id=\\"\${graph}\\"><img src=\\"data:image/png;base64,") <(base64 -w 0 \$plate) <(echo '"></div>') >> psmqc.html
      done 
    done
  echo "</body></html>" >> psmqc.html
  """
}

featqc_getpeptable
  .filter { it[0] == 'peptides' }
  .map { it -> it[1] }
  .set { featqc_peptides }
featuretables
  .merge(setnames_featqc)
  .combine(featqc_peptides)
  .view()
  .set { featqcinput }

normratios
  .toList()
  .map { it -> [it.collect() { it[0] }.sort(), it.collect() { it[1] }.sort()] }
  .set{ allsetnormratios }


process featQC {

  input:
  set val(acctype), file('feats'), val(setnames), file(peptable) from featqcinput
  set val(setnames), file('norm?') from allsetnormratios
  output:
  set val(acctype), file('featqc.html') into qccollect

  script:
  """
  ${normalize ? "count=1;for setn in ${setnames.join(' ')}; do echo '' >> norm\${count} ; tail -n+2 norm\${count} | sed \$'s/ - /\t'\${setn}\$'\t/'; ((count++)); done >> normtable" : ''}
  qc_protein.R ${setnames.size()} ${acctype} $peptable ${normalize ? 'normtable' : ''}
  echo "<html><body>" > featqc.html
  for graph in featyield precursorarea coverage isobaric nrpsms nrpsmsoverlapping percentage_onepsm normfac ms1nrpeps;
    do
    [ -e \$graph ] && paste -d \\\\0  <(echo "<div class=\\"chunk\\" id=\\"\${graph}\\"><img src=\\"data:image/png;base64,") <(base64 -w 0 \$graph) <(echo '"></div>') >> featqc.html
    done 
  echo "</body></html>" >> featqc.html
  """
}

qccollect
  .concat(psmqccollect)
  .toList()
  .map { it -> [it.collect() { it[0] }, it.collect() { it[1] }] }
  .set { collected_feats_qc }


process collectQC {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(acctypes), file('feat?') from collected_feats_qc
  val(plates) from qcplates

  output:
  file('qc.html')

  script:
  """
  count=1; for ac in ${acctypes.join(' ')}; do mv feat\$count \$ac.html; ((count++)); done
  qc_collect.py $params.name ${params.hirief ? "hirief" : "nofrac"} ${plates.join(' ')}
  """
}


/* 
 * STEP 3 - Output Description HTML
*/
process output_documentation {
    tag "$prefix"

    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}



/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/ddamsproteomics] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/ddamsproteomics] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/ddamsproteomics] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/ddamsproteomics] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/ddamsproteomics] Pipeline Complete"

}
