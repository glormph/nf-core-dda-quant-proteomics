#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'nf-core/ddamsproteomics': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'MSGF+': ['v_msgf.txt', r"([0-9\.]+)"],
    'Hardklor': ['v_hk.txt', r"([0-9\.]+)"],
    'Kronik': ['v_kr.txt', r"([0-9\.]+)"],
    'Percolator': ['v_perco.txt', r"([0-9\.]+)"],
    'msstitch': ['v_mss.txt', r"(\S+)"],
    'OpenMS': ['v_openms.txt', r"Version: ([0-9A-Z\-\.]+)"],
}

results = OrderedDict()
results['nf-core/ddamsproteomics'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'nf-core/ddamsproteomics-software-versions'
section_name: 'nf-core/ddamsproteomics Software Versions'
section_href: 'https://github.com/nf-core/ddamsproteomics'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")
