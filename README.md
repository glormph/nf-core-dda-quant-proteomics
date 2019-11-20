# lehtiolab/ddamsproteomics
**A Quantitative MS proteomics analysis pipeline**

[![Build Status](https://travis-ci.org/lehtiolab/ddamsproteomics.svg?branch=master)](https://travis-ci.org/lehtiolab/ddamsproteomics)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.04.1-brightgreen.svg)](https://www.nextflow.io/)
[![DOI](https://zenodo.org/badge/219955514.svg)](https://zenodo.org/badge/latestdoi/219955514)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/lehtiolab/ddamsproteomics.svg)](https://hub.docker.com/r/lehtiolab/ddamsproteomics)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
This workflow identifies peptides in mzML input data using [MSGF+](https://github.com/MSGFPlus/msgfplus), and [Percolator](https://github.com/percolator/percolator/), quantifies isobarically labeled samples with [OpenMS](https://github.com/openms/openms), and precursor peptides with [Hardklor]()/[Kronik](), and processes that output to formatted peptide and protein/gene tables using [Msstitch](https://github.com/glormph/msstitch). A considerable l

  
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


### Documentation
The lehtiolab/ddamsproteomics pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)
