#!/usr/bin/env python3

from jinja2 import Template
from lxml.html import parse, tostring
import sys
from collections import OrderedDict
import os


ppsms = {}
template = sys.argv[1]
searchname = sys.argv[2]
frac = sys.argv[3]
plateids = sys.argv[4:] 

templatetype = os.path.splitext(os.path.basename(template))[0]
with open(template) as fp: 
    main = Template(fp.read())
with open('psms.html') as fp:
    psmsel = parse(fp).find('body').findall('div')

psms = {x.attrib['id']: tostring(x, encoding='unicode') for x in psmsel if x.attrib['class'] == 'chunk'}
if frac == 'hirief':
    fryield = 'Fraction yield'
    for plateid in plateids:
        ppsms[plateid] = {x.attrib['id']: tostring(x, encoding='unicode') for x in psmsel if x.attrib['class'] == 'chunk {}'.format(plateid)}
else:
    fryield = 'Yield'
    ppsms['No plate'] = {x.attrib['id']: tostring(x, encoding='unicode') for x in psmsel if x.attrib['class'] == 'chunk noplates'}

titles = {'psm-scans': '# PSMs and scans', 'miscleav': 'Missed cleavages',
          'missing-tmt': 'Isobaric missing values', 'fryield': fryield,
          'retentiontime': 'Retention time', 'precerror': 'Precursor error',
          'msgfscore': 'MSGF Score',
          'featyield': 'Identifications', 'isobaric': 'Isobaric intensities',
          'precursorarea': 'Precursor area intensity',
          'nrpsms': '# PSMs used for isobaric quantitation per identification',
          'nrpsmsoverlapping': '# PSMs used for isobaric quantitation per identification for only complete overlapping set',
          'percentage_onepsm': 'Percentage of identifications with >1 quantifying PSM in the complete overlapping set',
          'ms1nrpeps': '# peptides with MS1 quant per protein (top 3 used)',
          'coverage': 'Overall protein coverage',
}
featnames = {
        'qc_light': {'peptides': 'Peptides', 'proteins': 'Proteins', 'genes': 'Proteins (genecentric)'},
        'qc_full': {'assoc': 'Gene symbols', 'peptides': 'Peptides', 'proteins': 'Proteins', 'genes': 'Genes'},
        }
# FIXME  use this for ppsms!
#PSMs/protein for quant (median)

tablefields = {
        'nr_sets': 'IDed in # overlapping sets', 
        'Set': 'Experiment set', 
        'no_pep_proteins': 'Peptides/protein (unique, median)',
        'no_pep_genes': 'Peptides/protein (genecentric, unique, median)',
        'nr_proteins': 'Proteins ID and quant. (1%FDR)',
        'nr_genes': 'Proteins ID and quant. (gene centric, 1%FDR)',
        'nr_assoc': 'Proteins ID and quant. (symbol centric, 1%FDR)',
        'Non-shared (unique)': 'Peptides (unique, 1%FDR)',
        'psmcount': 'PSMs (total)',
        }

graphs = OrderedDict()
feattypes = {
    'qc_light': ['peptides', 'proteins', 'genes'],
    'qc_full': ['peptides', 'proteins', 'genes', 'assoc'],
    }

for feat in feattypes[templatetype]:
    try:
        with open('{}.html'.format(feat)) as fp:
            graphs[feat] = {x.attrib['id']: tostring(x, encoding='unicode') for x in parse(fp).find('body').findall('div') if x.attrib['class'] == 'chunk'}
    except IOError as e:
        print(feat, e)

def parse_table(fn):
    table = {'_rows': {}}
    with open(fn) as fp:
        header = next(fp).strip('\n').split('\t')
        table['_fields'] = header
        for line in fp:
            line = line.strip('\n').split('\t')
            line = {header[x]: line[x] for x in range(0,len(line))}
            table['_rows'][line[header[0]]] = line
    return table

summaries = {'qc_light': 'summary_light', 'qc_full': 'summary'}
sumtable = parse_table(summaries[templatetype])
overlaptables = {}
for feat in feattypes[templatetype]:
    try:
        overlap = parse_table('{}_overlap'.format(feat))
    except IOError:
        pass
    else:
        overlaptables[feat] = overlap
if templatetype == 'qc_light' and 'genes' in overlaptables:
    overlaptables.pop('proteins')
    
with open('{}.html'.format(templatetype), 'w') as fp:
    fp.write(main.render(sumtable=sumtable, overlap=overlaptables, tablefields=tablefields, hirief=frac, searchname=searchname, titles=titles, featnames=featnames[templatetype], psms=psms, firstplate=sorted(ppsms.keys())[0], ppsms=ppsms, features=graphs))
