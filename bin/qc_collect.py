#!/usr/bin/env python3

from jinja2 import Template
from lxml.html import parse, tostring
import sys
from collections import OrderedDict
import os

main = Template("""<!DOCTYPE html>
<html lang="en">
<head>
    <title>Lehtio proteomics QC report</title>
<link rel="stylesheet" type="text/css" href="https://cdnjs.cloudflare.com/ajax/libs/bulma/0.6.2/css/bulma.min.css">
</head>
<body>
<div class="container">
  {{ logo }} 
  <hr>
  <h3 class="title is-3">Protein/peptide level QC</h3>
{% for graphtype in ["featyield", "precursorarea", "isobaric", "nrpsms", "nrpsmsoverlapping", "percentage_onepsm"] %}
  {% if graphtype in features[features.keys()[0]] %}
  <h4 class="title is-4">{{ titles[graphtype] }}</h4>
  <div class="columns">
    {% for feat in features %}
    <div class="column">
      <h5 class="title is-5">{{ featnames[feat] }}</h5>
      {{ features[feat][graphtype] }}
    </div>
    {% endfor %}
    {% if graphtype == "isobaric" %}
    <div class="column">
      <h5 class="title is-5">Median centering</h5>
        {{ norm }}
    </div>
    {% endif %}
  </div>
<hr>
{% endif %}
{% endfor %}
</div>
{% if 'proteins' in features %}
<div class="container">
  <h4 class="title is-4">Overall protein coverage</h3>
    {{ features.proteins.coverage }}
</div>
<hr>
{% endif %}
<div class="container">
  <h3 class="title is-3">PSM level QC</h3>
  {% if hirief == 'hirief' %}
  <div class="columns">
    {% for graph in psms %}
    <div class="column">
      <h5 class="title is-5">{{ titles[graph] }}</h3>
      {{ psms[graph] }}
    </div>
    {% endfor %}
  </div>
  {% endif %}
</div>
{% for graph in ppsms[firstplate] %}
<div class="container">
  <h4 class="title is-4">{{ titles[graph] }}</h4>
{% for plate, graphs in ppsms|dictsort %}
<div class="container">
  <h5 class="title is-5">Plate: {{ plate }}</h5>
  {{ ppsms[plate][graph] }}
</div>
{% endfor %}
</div>
{% endfor %}
</body>
</html>
""")


# FIXME
# PSMs
# coverage if protein
# venn diagrams
ppsms = {}
searchname = sys.argv[1]
frac = sys.argv[2]
if frac == 'hirief':
    fryield = 'Fraction yield'
    for ppsm in sys.argv[3:]:
        plate = os.path.basename(ppsm).replace('_psms.html', '')
        with open(ppsm) as fp:
            ppsms[plate] = {x.attrib['id']: tostring(x, encoding='unicode') for x in parse(fp).find('body').findall('div') if x.attrib['class'] == 'chunk'}
    
    with open('psms.html') as fp:
        psms = {x.attrib['id']: tostring(x, encoding='unicode') for x in parse(fp).find('body').findall('div') if x.attrib['class'] == 'chunk'}
else:
    fryield = 'Yield'
    psms = {}
    with open('psms.html') as fp:
        ppsms['No plate'] = {x.attrib['id']: tostring(x, encoding='unicode') for x in parse(fp).find('body').findall('div') if x.attrib['class'] == 'chunk'}
    
titles = {'psm-scans': '# PSMs and scans', 'miscleav': 'Missed cleavages',
          'missing-tmt': 'Isobaric missing values', 'fr-yield': fryield,
          'retentiontime': 'Retention time', 'prec-error': 'Precursor error',
          'featyield': 'Identifications', 'isobaric': 'Isobaric intensities',
          'precursorarea': 'Precursor area intensity',
          'nrpsms': '# PSMs used for isobaric quantitation per identification',
          'nrpsmsoverlapping': '# PSMs used for isobaric quantitation per identification for only complete overlapping set',
          'percentage_onepsm': 'Percentage of identifications with >1 quantifying PSM in the complete overlapping set',
          'coverage': 'Overall protein coverage',
}
featnames = {'assoc': 'Gene symbols', 'peptides': 'Peptides', 'proteins': 'Proteins', 'genes': 'Genes'}

graphs = OrderedDict()
for feat in ['peptides', 'proteins', 'genes', 'assoc']:
    try:
        with open('{}.html'.format(feat)) as fp:
            graphs[feat] = {x.attrib['id']: tostring(x, encoding='unicode') for x in parse(fp).find('body').findall('div') if x.attrib['class'] == 'chunk'}
    except IOError as e:
        print(feat, e)
try:
    with open('norm.html') as fp:
        normgraph = [tostring(x, encoding='unicode') for x in parse(fp).find('body').findall('div') if x.attrib['class'] == 'chunk'][0]
except IOError:
    print('No normalization file')
    normgraph = False
except AssertionError:
    print('Normalization file not XML, assuming no normalizing done')
    normgraph = False

with open('qc.html', 'w') as fp:
    fp.write(main.render(hirief=frac, searchname=searchname, titles=titles, featnames=featnames, psms=psms, firstplate=sorted(ppsms.keys())[0], ppsms=ppsms, features=graphs, norm=normgraph))
