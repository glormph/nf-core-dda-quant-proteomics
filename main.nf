#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/lehtio-quant-proteomics
========================================================================================
 nf-core/lehtio-quant-proteomics Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/lehtio-quant-proteomics
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     nf-core/lehtio-quant-proteomics v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/lehtio-quant-proteomics --reads '*_R{1,2}.fastq.gz' -profile standard,docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --genome                      Name of iGenomes reference
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test

    Options:
      --singleEnd                   Specifies that the input is single end reads

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference

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

params.isobaric = false
params.activation = 'hcd' // Only for isobaric quantification
params.outdir = 'results'
params.normalize = false
params.genes = false
params.symbols = false
params.fastadelim = false
params.genefield = false
params.speclookup = false
params.quantlookup = false
params.hirief = false
params.onlypeptides = false
params.noquant = false

// Validate and set file inputs
mods = file(params.mods)
if( !mods.exists() ) exit 1, "Modification file not found: ${params.mods}"
tdb = file(params.tdb)
if( !tdb.exists() ) exit 1, "Target fasta DB file not found: ${params.tdb}"

// file templates
hkconf = file('assets/hardklor.conf')
qcknitrplatepsms = file('assets/knitr_psms_perplate.Rhtml')
qcknitrnofrpsms = file('assets/knitr_psms_nofr.Rhtml')
qcknitrpsms = file('assets/knitr_psms.Rhtml')
qcknitrprot = file('assets/knitr_prot.Rhtml')
qcknitrnormfac = file('assets/knitr_iso_norm.Rhtml')

if (params.martmap) {
  martmap = file(params.martmap)
  if( !martmap.exists() ) exit 1, "Biomart ENSEMBL mapping file not found: ${params.martmap}"
}
if (params.pipep) {
  trainingpep = file(params.pipep)
  if( !trainingpep.exists() ) exit 1, "Peptide pI data file not found: ${params.pipep}"
}
output_docs = file("$baseDir/docs/output.md")

// Parse input variables values
accolmap = [peptides: 12, proteins: 14, genes: 17, assoc: 18]
setdenoms = [:]
if (!(params.noquant) && params.isobaric) {
  params.denoms.tokenize(' ').each{ it -> x=it.tokenize(':'); setdenoms.put(x[0], x[1..-1])}
}
normalize = (params.normalize && params.isobaric) ? true: false
activations = [hcd:'High-energy collision-induced dissociation', cid:'Collision-induced dissociation', etd:'Electron transfer dissociation']
activationtype = activations[params.activation]
plextype = params.isobaric ? params.isobaric.replaceFirst(/[0-9]+plex/, "") : 'false'
massshift = [tmt:0.0013, itraq:0.00125, false:0][plextype]
msgfprotocol = [tmt:4, itraq:2, false:0][plextype]
instrument = params.instrument ? params.instrument : 'qe'
msgfinstrument = [velos:1, qe:3, false:0][instrument]


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

nf-core/lehtio-quant-proteomics v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/lehtio-quant-proteomics'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['mzMLs']        = params.reads
summary['Target DB']    = params.tdb
summary['Modifications']= params.mods
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
    id: 'nf-core-lehtio-quant-proteomics-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/lehtio-quant-proteomics Workflow Summary'
    section_href: 'https://github.com/nf-core/lehtio-quant-proteomics'
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
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

if (params.mzmlPaths) {
  // Profile 'test' delivers mzmlPaths
  Channel
    .from(params.mzmlPaths)
    .view()
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

mzml_in
  .tap { sets }
  .map { it -> [file(it[0]), it[1], it[2] ? it[2] : it[1], it[3] ? it[3] : 'NA' ]} // create file, set plate to setname, and fraction to NA if there is none
    .view()
  .tap { strips }
  .map { it -> [it[1], it[0].baseName.replaceFirst(/.*\/(\S+)\.mzML/, "\$1"), it[0], it[2], it[3]] }
  .tap{ mzmlfiles; mzml_isobaric; mzml_hklor; mzml_msgf }
  .count()
  .set{ amount_mzml }

sets
  .map{ it -> it[1] }
  .unique()
  .tap { setnames_psm } 
  .collect()
  .map { it -> [it] }
  .into { setnames_featqc; setnames_psmqc }

strips
  .map { it -> it[2] }
  .unique()
  .toList()
  .set { strips_for_deltapi }

process IsobaricQuant {

  when: !params.noquant && params.isobaric && !params.quantlookup

  input:
  set val(setname), val(sample), file(infile), val(platename), val(fraction)from mzml_isobaric

  output:
  set val(sample), file("${infile}.consensusXML") into isobaricxml

  """
  source activate openms-2.4.0
  IsobaricAnalyzer  -type $params.isobaric -in $infile -out "${infile}.consensusXML" -extraction:select_activation "$activationtype" -extraction:reporter_mass_shift $massshift -extraction:min_precursor_intensity 1.0 -extraction:keep_unannotated_precursor true -quantification:isotope_correction true 
  """
}


process hardklor {
  when: !params.quantlookup && !params.noquant

  input:
  set val(setname), val(sample), file(infile), val(platename), val(fraction) from mzml_hklor
  file hkconf

  output:
  set val(sample), file('hardklor.out'), file(infile) into hk_out

  """
  cp $hkconf config
  echo "$infile" hardklor.out >> config
  hardklor config
  """
}


process kronik {
  when: !params.quantlookup && !params.noquant

  input:
  set val(sample), file('hardklor.out'), file(mzml) from hk_out 

  output:
  set val(sample), file("${sample}.kr"), file(mzml) into kronik_out

  """
  kronik -c 5 -d 3 -g 1 -m 8000 -n 600 -p 10 hardklor.out ${sample}.kr
  """
}


mzmlfiles
  .buffer(size: amount_mzml.value)
  .map { it.sort( {a, b -> a[1] <=> b[1]}) } // sort on sample for consistent .sh script in -resume
  .map { it -> [it.collect() { it[0] }, it.collect() { it[2] }, it.collect() { it[3] } ] } // lists: [sets], [mzmlfiles], [plates]
  .into { mzmlfiles_all; mzmlfiles_all_count }


if (params.speclookup && !params.quantlookup) {
  Channel
    .fromPath(params.speclookup)
    .into{ spec_lookup; countlookup }
} 
if (!params.speclookup && params.quantlookup) {
  Channel
    .fromPath(params.quantlookup)
    .into { spec_lookup; quant_lookup; countlookup }
} 


process createSpectraLookup {

  when: !(params.speclookup || params.quantlookup)

  input:
  set val(setnames), file(mzmlfiles), val(platenames) from mzmlfiles_all

  output:
  file 'mslookup_db.sqlite' into newspeclookup 

  script:
  """
  msslookup spectra -i ${mzmlfiles.join(' ')} --setnames ${setnames.join(' ')}
  """
}


isoquant_amount = params.isobaric ? amount_mzml.value : 1
isobaricxml
  .ifEmpty(['NA', 'NA', 'NA'])
  .buffer(size: isoquant_amount)
  .map { it.sort({a, b -> a[0] <=> b[0]}) }
  .map { it -> [it.collect() { it[0] }, it.collect() { it[1] }] } // samples, isoxml
  .set { isofiles_sets }

kronik_out
  .ifEmpty(['NA', 'NA'])
  .buffer(size: amount_mzml.value)
  .map { it.sort({a, b -> a[0] <=> b[0]}) }
  .map { it -> [it.collect() { it[0] }, it.collect() { it[1] }, it.collect() { it[2] }] } // samples, kronikout, mzml
  .set { krfiles_sets }


if (params.noquant && !(params.speclookup || params.quantlookup)) {
  newspeclookup
    .into { quant_lookup; spec_lookup; countlookup }
} else if (!(params.speclookup || params.quantlookup)) {
  newspeclookup
    .into { spec_lookup; countlookup }
}

process quantLookup {

  when: !params.quantlookup && !params.noquant

  input:
  file lookup from spec_lookup
  set val(isosamples), file(isofns) from isofiles_sets
  set val(krsamples), file(krfns), file(mzmls) from krfiles_sets

  output:
  file 'db.sqlite' into newquantlookup

  script:
  if (params.isobaric)
  """
  cp $lookup db.sqlite
  msslookup ms1quant --dbfile db.sqlite -i ${krfns.join(' ')} --spectra ${mzmls.join(' ')} --quanttype kronik --mztol 20.0 --mztoltype ppm --rttol 5.0 
  msslookup isoquant --dbfile db.sqlite -i ${isofns.join(' ')} --spectra ${isosamples.collect{ x -> x + '.mzML' }.join(' ')}
  """
  else
  """
  cp $lookup db.sqlite
  msslookup ms1quant --dbfile db.sqlite -i ${krfns.join(' ')} --spectra ${mzmls.join(' ')} --quanttype kronik --mztol 20.0 --mztoltype ppm --rttol 5.0 
  """
}


if (!params.quantlookup && !params.noquant) {
  newquantlookup
    .into { quant_lookup }
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


if (params.hirief) {
  specfilems2.set { scans_platecount }
} else {
  specfilems2
    .map { it -> [it[3], ['noplates']] }
    .into { scans_platecount; scans_result }
}


process countMS2sPerPlate {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true 
  when: params.hirief

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

if (params.hirief) {
  scans_perplate.set { scans_result }
}

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

tmzidtsv_perco
  .concat(dmzidtsv_perco)
  .groupTuple(by: 1)
  .combine(quant_lookup)
  .set { prepsm }

if (params.hirief) {
  strips_for_deltapi
    .map { it -> [it, trainingpep] }
    .set { stripannot }
} else {
  strips_for_deltapi
    .map { it -> [it, false, false]}
    .set { stripannot }
}

process createPSMTable {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: {["target_psmlookup.sql", "target_psmtable.txt", "decoy_psmtable.txt"].contains(it) ? it : null}

  input:
  set val(setnames), val(td), file('psms?'), file('lookup') from prepsm
  set file(tdb), file(ddb) from searchdbs
  set val(allstrips), file(trainingpep) from stripannot

  output:
  set val(td), file("${td}_psmtable.txt") into psm_result
  set val(td), file({setnames.collect() { it + '.tsv' }}) into setpsmtables
  set val(td), file("${td}_psmlookup.sql") into psmlookup

  script:
  """
  msspsmtable merge -i psms* -o psms.txt
  msspsmtable conffilt -i psms.txt -o filtpsm --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'PSM q-value'
  msspsmtable conffilt -i filtpsm -o filtpep --confidence-better lower --confidence-lvl 0.01 --confcolpattern 'peptide q-value'
  cp lookup psmlookup
  msslookup psms -i filtpep --dbfile psmlookup ${params.onlypeptides ? '' : "--fasta ${td == 'target' ? tdb : "${ddb} --decoy"}"} ${params.martmap ? "--map ${martmap}" : ''}
  msspsmtable specdata -i filtpep --dbfile psmlookup -o prepsms.txt
  ${!params.noquant ? "msspsmtable quant -i prepsms.txt -o qpsms.txt --dbfile psmlookup --precursor ${params.isobaric && td=='target' ? '--isobaric' : ''}" : 'mv prepsms.txt qpsms.txt'}
  sed 's/\\#SpecFile/SpectraFile/' -i qpsms.txt
  ${!params.onlypeptides ? "msspsmtable genes -i qpsms.txt -o gpsms --dbfile psmlookup" : ''}
  ${!params.onlypeptides ? "msslookup proteingroup -i qpsms.txt --dbfile psmlookup" : ''}
  ${!params.onlypeptides ? "msspsmtable proteingroup -i gpsms -o pgpsms --dbfile psmlookup" : 'mv qpsms.txt pgpsms'}
  ${params.hirief ? "peptide_pi_annotator.py -i $trainingpep -p pgpsms --o dppsms --stripcolpattern Strip --pepcolpattern Peptide --fraccolpattern Fraction --strippatterns ${allstrips.join(' ')} --intercepts ${allstrips.collect() { params.strips[it].intercept}.join(' ')} --widths ${allstrips.collect() { params.strips[it].fr_width}.join(' ')} --ignoremods \'*\'" : ''}
  msspsmtable split -i ${params.hirief ? 'dppsms' : 'pgpsms'} --bioset
  mv ${params.hirief ? 'dppsms' : 'pgpsms'} ${td}_psmtable.txt
  mv psmlookup ${td}_psmlookup.sql
  """
}

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
  set val(td), val(setname), file('psms'), file('peptides') into prepep
  """
  msspeptable psm2pep -i psms -o peptides --scorecolpattern svm --spectracol 1 ${!params.noquant && params.isobaric && td == 'target' ? "--isobquantcolpattern plex" : "" } ${!params.noquant ? "--ms1quantcolpattern area" : ""}
  """
}


process shuffleHeaderPeptidesMakeProttables {
 
  input:
  set val(td), val(setname), file(psms), file('peptides') from prepep

  output:
  set val(setname), val(td), file(psms), file('peptide_table.txt') into peptides
  set val(setname), val(td), file(psms), file('proteins'), val('proteins') into proteins
  set val(setname), val(td), file(psms), file('genes'), val('genes') into genes
  set val(setname), val(td), file(psms), file('symbols'), val('assoc') into symbols

  // protein, gene and symbol table are outputted regardless of necessity
  script:
  col = accolmap.peptides + 1  // psm2pep adds a column
  """
  paste <( cut -f ${col} peptides) <( cut -f 1-${col-1},${col+1}-500 peptides) > peptide_table.txt
  echo Protein accession |tee proteins genes symbols
  tail -n+2 psms|cut -f ${accolmap.proteins}|grep -v '\\;'|grep -v "^\$"|sort|uniq >> proteins
  tail -n+2 psms|cut -f ${accolmap.genes}|grep -v '\\;'|grep -v "^\$"|sort|uniq >> genes
  tail -n+2 psms|cut -f ${accolmap.assoc}|grep -v '\\;'|grep -v "^\$"|sort|uniq >> symbols
  """
}

proteins
  .tap { proteins_pre }
  .filter { it[1] == 'target' }
  .set { tprot_norm }

process ratioNormalizeProteins {

  when: normalize

  input:
  set val(setname), val(td), file('psms'), file('proteins'), val(acctype) from tprot_norm
  output:
  set val(setname), file('proteinratios') into proteinratios
  """
  msspsmtable isoratio -i psms -o proteinratios --protcol ${accolmap.proteins} --targettable proteins --isobquantcolpattern plex --minint 0.1 --denompatterns ${setdenoms[setname].join(' ')}
  """
}

tpep = Channel.create()
dpep = Channel.create()
peptides.choice(tpep, dpep) { it -> it[1] == 'target' ? 0 : 1}

if (normalize) {
  tpep
    .join(proteinratios)
    .concat(dpep)
    .set { peptable_quant }
} else {
  tpep
    .concat(dpep)
    .set { peptable_quant }
}

process postprodPeptideTable {

  input:
  set val(setname), val(td), file('psms'), file('peptides'), file(pratios) from peptable_quant

  output:
  set val(setname), val(td), file("${setname}_linmod"), file(pratios) into pepslinmod
  set val(setname), val('peptides'), val(td), file("${setname}_linmod") into peptides_out
  set val(setname), file('normratiosused') optional true into normratios

  script:
  if (!params.noquant && params.isobaric && td=='target')
  """
  msspsmtable isoratio -i psms -o pepisoquant --targettable peptides --protcol ${accolmap.peptides} --isobquantcolpattern plex --minint 0.1 --denompatterns ${setdenoms[setname].join(' ')} ${normalize ? "--normalize median --norm-ratios $pratios" : ''} > normratiosused
  mv pepisoquant peptide_table.txt
  msspeptable modelqvals -i peptide_table.txt -o ${setname}_linmod --scorecolpattern svm --fdrcolpattern '^q-value'
  """
  else
  """
  mv peptides peptide_table.txt
  msspeptable modelqvals -i peptide_table.txt -o ${setname}_linmod --scorecolpattern svm --fdrcolpattern '^q-value'
  """
}


if (params.genes && params.symbols) { 
  pepslinmod.tap { pepsg; pepss }.concat(pepsg, pepss).set { pepslinmod_prot }
  proteins_pre.concat(genes).concat(symbols).set { pgstables } 
} else if (params.genes) { 
  pepslinmod.tap { pepsg }.concat(pepsg).set { pepslinmod_prot }
  proteins_pre.concat(genes).set { pgstables } 
} else { 
  pepslinmod.set { pepslinmod_prot }
  proteins_pre.set { pgstables }
}


pgstables
  .join(pepslinmod_prot, by: [0,1])
  .set { prepgs_in }

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
  msspsmtable isoratio -i psms -o proteintable --protcol ${accolmap[acctype]} --targettable protms1 --isobquantcolpattern plex --minint 0.1 --denompatterns ${setdenoms[setname].join(' ')} ${normalize ? '--norm-ratios pratios --normalize median': ''}
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

process proteinPeptideSetMerge {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true, saveAs: { it == "proteintable" ? "${outname}_table.txt": null}

  input:
  set val(setnames), val(acctype), file(tables) from ptables_to_merge
  file(lookup) from tlookup
  
  output:
  set val(acctype), file('proteintable') into featuretables

  script:
  outname = (acctype == 'assoc') ? 'symbols' : acctype
  """
  cp $lookup db.sqlite
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
  set val('psms'), file('knitr.html') into psmqccollect
  file('*_psms.html') optional true into platepsmscoll
  // TODO no proteins == no coverage for pep centric
  script:
  if (params.hirief)
  """
  #!/usr/bin/env Rscript
  library(ggplot2)
  library(reshape2)
  library(knitr)
  feats = read.table("psms", header=T, sep="\\t", comment.char = "", quote = "")
  feats\$plateID = paste(feats\$Biological.set, feats\$Strip, sep='_')
  nrsets=${setnames[0].size()}
  amount_ms2 = read.table("scans")
  knitr::knit2html("$qcknitrpsms", output="knitr.html")
  for (plateid in c(${plates.collect() {"\"${it}\"" }.join(',') })) {
    knitr::knit2html("$qcknitrplatepsms", output=paste(plateid, "psms.html", sep="_"))
  }
  file.remove('knitr_psms.html')
  """
  else
  """
  #!/usr/bin/env Rscript
  library(ggplot2)
  library(reshape2)
  library(knitr)
  nrsets=${setnames[0].size()}
  feats = read.table("psms", header=T, sep="\\t", comment.char = "", quote = "")
  amount_ms2 = read.table("scans", sep="|", header=F)
  knitr::knit2html("$qcknitrnofrpsms", output="knitr.html")
  file.remove('knitr_psms.html')
  """
}

featuretables
  .merge(setnames_featqc)
  .set { featqcinput }

normratios
  .toList()
  .map { it -> [it.collect() { it[0] }.sort(), it.collect() { it[1] }.sort()] }
  .set{ allsetnormratios }


process normRatioParse {

  input:
  set val(setnames), file('norm?') from allsetnormratios

  output:
  file('normtable_sets') into normtable
  """
  count=1;for setn in ${setnames.join(' ')}; do echo "" >> norm"\$count" ; tail -n+2 norm"\$count" | sed \$'s/ - /\t'"\$setn"\$'\t/'; ((count++)); done >> normtable_sets
  """
}


process featQC {

  input:
  set val(acctype), file('feats'), val(setnames) from featqcinput
  file('normtable') from normtable
  output:
  set val(acctype), file('knitr.html') into qccollect
  file('normalizefactors.html') optional true into normhtml

  script:
  """
  #!/usr/bin/env Rscript
  library(ggplot2)
  library(forcats)
  library(reshape2)
  library(knitr)
  nrsets=${setnames.size()}
  feats = read.table("feats", header=T, sep="\\t", comment.char = "", quote = "")
  feattype="$acctype"
  knitr::knit2html("$qcknitrprot", output="knitr.html")
  ${normalize ? 'normtable="normtable"' : ''}
  ${normalize ? "knitr::knit2html(\"$qcknitrnormfac\", output=\"normalizefactors.html\")": ''}
  """
}

qccollect
  .concat(psmqccollect)
  .toList()
  .map { it -> [it.collect() { it[0] }, it.collect() { it[1] }] }
  .set { collected_feats_qc }

if (!params.hirief) {
  Channel.from([1]).set { platepsmscoll }
}
process collectQC {

  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  set val(acctypes), file('feat?') from collected_feats_qc
  file('norm.html') from normhtml
  file(ppsms) from platepsmscoll

  output:
  file('qc.html')

  script:
  if (params.hirief)
  """
  count=1; for ac in ${acctypes.join(' ')}; do mv feat\$count \$ac.html; ((count++)); done
  qc_collect.py $params.name hirief ${ppsms.join(' ')}
  """
  else
  """
  count=1; for ac in ${acctypes.join(' ')}; do mv feat\$count \$ac.html; ((count++)); done
  qc_collect.py $params.name nofrac
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
    def subject = "[nf-core/lehtio-quant-proteomics] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/lehtio-quant-proteomics] FAILED: $workflow.runName"
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
          log.info "[nf-core/lehtio-quant-proteomics] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/lehtio-quant-proteomics] Sent summary e-mail to $params.email (mail)"
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

    log.info "[nf-core/lehtio-quant-proteomics] Pipeline Complete"

}
