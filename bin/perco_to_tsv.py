#!/usr/bin/env python3
import re
import argparse
from glob import glob


def count_missed_cleavage(full_pepseq, count=0):
    '''Regex .*[KR][^P] matches until the end and checks if there is a final
    charachter so this will not match the tryptic residue'''
    pepseq = re.sub('[\+\-]\d*.\d*', '', full_pepseq)
    match = re.match('.*[KR][^P]', pepseq)
    if match:
        count += 1
        return count_missed_cleavage(match.group()[:-1], count)
    else:
        return count

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--perco', dest='perco')
parser.add_argument('--fractions', dest='fractions', nargs='+')
parser.add_argument('--plates', dest='plates', nargs='+')
args = parser.parse_args()

mzidtsvfns = sorted(glob('mzidtsv*'))
mzidfns = sorted(glob('mzident*'))
fractions = args.fractions
perco = args.perco
plates = args.plates
from app.readers import pycolator, xml, tsv, mzidplus
import os
ns = xml.get_namespace_from_top(perco, None) 
psms = {p.attrib['{%s}psm_id' % ns['xmlns']]: p for p in pycolator.generate_psms(perco, ns)}
decoys = {True: 0, False: 0}
for psm in sorted([(pid, float(p.find('{%s}svm_score' % ns['xmlns']).text), p) for pid, p in psms.items()], reverse=True, key=lambda x:x[1]):
    pdecoy = psm[2].attrib['{%s}decoy' % ns['xmlns']] == 'true'
    decoys[pdecoy] += 1
    try:
        psms[psm[0]] = {'decoy': pdecoy, 'svm': psm[1], 'qval': decoys[True]/decoys[False]}  # T-TDC
    except ZeroDivisionError:
        psms[psm[0]] = {'decoy': pdecoy, 'svm': psm[1], 'qval': 1}  # T-TDC
decoys = {'true': 0, 'false': 0}
for svm, pep in sorted([(float(x.find('{%s}svm_score' % ns['xmlns']).text), x) for x in pycolator.generate_peptides(perco, ns)], reverse=True, key=lambda x:x[0]):
    decoys[pep.attrib['{%s}decoy' % ns['xmlns']]] += 1
    try:
        [psms[pid.text].update({'pepqval': decoys['true']/decoys['false']}) for pid in pep.find('{%s}psm_ids' % ns['xmlns'])]
    except ZeroDivisionError:
        [psms[pid.text].update({'pepqval': 1}) for pid in pep.find('{%s}psm_ids' % ns['xmlns'])]
oldheader = tsv.get_tsv_header(mzidtsvfns[0])
header = oldheader + ['percolator svm-score', 'PSM q-value', 'peptide q-value', 'Strip', 'Fraction', 'missed_cleavage']
with open('tmzidperco', 'w') as tfp, open('dmzidperco', 'w') as dfp:
    tfp.write('\t'.join(header))
    dfp.write('\t'.join(header))
    for fnix, mzidfn in enumerate(mzidfns):
        mzns = mzidplus.get_mzid_namespace(mzidfn)
        inpsms = tsv.generate_tsv_psms(mzidtsvfns[fnix], oldheader)
        for specidr in mzidplus.mzid_spec_result_generator(mzidfn, mzns):
            for specidi in specidr.findall('{%s}SpectrumIdentificationItem' % mzns['xmlns']):
                psm = next(inpsms)
                # percolator psm ID is: samplename_SII_scanindex_rank_scannr_charge_rank
                scanindex, rank = specidi.attrib['id'].replace('SII_', '').split('_')
                scan = {x.split('=')[0]: x.split('=')[1] for x in specidr.attrib['spectrumID'].split(' ')}['scan']
                outpsm = {k: v for k,v in psm.items()}
                spfile = os.path.splitext(psm['#SpecFile'])[0]
                try:
                    percopsm = psms['{fn}_SII_{ix}_{rk}_{sc}_{ch}_{rk}'.format(fn=spfile, ix=scanindex, sc=scan, rk=rank, ch=psm['Charge'])]
                except KeyError:
                    continue
                outpsm.update({'percolator svm-score': percopsm['svm'], 'PSM q-value': percopsm['qval'], 'peptide q-value': percopsm['pepqval'], 'Strip': plates[fnix], 'Fraction': fractions[fnix], 'missed_cleavage': count_missed_cleavage(outpsm['Peptide'])})
                if percopsm['decoy']:
                    dfp.write('\n')
                    dfp.write('\t'.join([str(outpsm[k]) for k in header]))
                else:
                    outpsm['Protein'] = ';'.join([x for x in outpsm['Protein'].split(';') if 'decoy_' not in x])
                    tfp.write('\n')
                    tfp.write('\t'.join([str(outpsm[k]) for k in header]))

