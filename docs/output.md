# nf-core/ddamsproteomics: Output

This document describes the output produced by the pipeline. 


## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [MSGF+](#msgf) - Peptide identification search engine
* [Percolator](#percolator) - Target-decoy scoring
* [OpenMS](#openms) - Quantification of isobaric tags
* [Hardklor/Kronik](#hardklor) - Quantification of precursor peptides
* [Msstitch](#msstitch) - Post processing, protein inference

## MSGF+
[MSGF+](https://omics.pnl.gov/software/ms-gf) (aka MSGF+ or MSGFPlus) performs peptide identification by scoring MS/MS spectra against peptides derived from a protein sequence database.


## Percolator
[Percolator](http://percolator.ms/) is a semi-supervised machine learning program to better separate target vs decoy peptide scoring.


## OpenMS
[OpenMS](http://www.openms.de/) is a library that contains a large amount of tools for MS data analysis. This workflow uses its isobaric quantification program to extract peak intenstities for isobaric multiplex samples.


## Hardklor/Kronik
[Hardklor](https://proteome.gs.washington.edu/software/hardklor/) identifies peptide features in MS1 data by reducing isotope distributions to a single monoisotopic mass. It in tandem with its [Kronik](https://github.com/mhoopmann/kronik) utility (which summarizes the Hardklor results from LC/MS data) can be used to quantify MS1 peptide features.


## Msstitch
[Msstitch](https://github.com/glormph/msstitch) is a package to merge identification and quantification PSM data, reporting PSM, peptide, protein and gene tables, adding q-values, quantitfications, protein groups, etc.
