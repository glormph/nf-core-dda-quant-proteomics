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
{% for graphtype in ["featyield", "precursorarea", "isobaric", "nrpsms", "nrpsmsoverlapping", 
                     "percentage_onepsm", "ms1nrpeps"] %}
  {% if graphtype in features['peptides'] or ('proteins' in features and graphtype in features['proteins']) %}
  <h4 class="title is-4">{{ titles[graphtype] }}</h4>
  <div class="columns">
    {% for feat in features %}
    <div class="column">
      <h5 class="title is-5">{{ featnames[feat] }}</h5>
      {{ features[feat][graphtype] }}
    </div>
    {% endfor %}
    {% if graphtype == "isobaric" and 'normfac' in features['proteins'] %}
    <div class="column">
      <h5 class="title is-5">Median centering</h5>
        {{ features['proteins']['normfac']}}
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
  {% else %}
    {% for graph in psms %}
      <h5 class="title is-5">{{ titles[graph] }}</h5>
      {{ psms[graph] }}
    {% endfor %}
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
plateids = sys.argv[3:] # FIXME  use this for ppsms!

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
featnames = {'assoc': 'Gene symbols', 'peptides': 'Peptides', 'proteins': 'Proteins', 'genes': 'Genes'}

graphs = OrderedDict()
for feat in ['peptides', 'proteins', 'genes', 'assoc']:
    try:
        with open('{}.html'.format(feat)) as fp:
            graphs[feat] = {x.attrib['id']: tostring(x, encoding='unicode') for x in parse(fp).find('body').findall('div') if x.attrib['class'] == 'chunk'}
    except IOError as e:
        print(feat, e)

with open('qc.html', 'w') as fp:
    fp.write(main.render(hirief=frac, searchname=searchname, titles=titles, featnames=featnames, psms=psms, firstplate=sorted(ppsms.keys())[0], ppsms=ppsms, features=graphs))
