# nf-core/ddamsproteomics: Troubleshooting

## Input files not found

If no file is picked up then something is wrong with your input file declaration

The path must be enclosed in quotes (`'` or `"`)


## Data organization
The pipeline can't take a list of multiple input files - it takes a glob expression. If your input files are scattered in different paths then we recommend that you generate a directory with symlinked files. 

## Data organization
The pipeline can't take a list of multiple input files - it takes a glob expression. If your input files are scattered in different paths then we recommend that you generate a directory with symlinked files. 

## Extra resources and getting help
If you still have an issue with running the pipeline then feel free to contact us.
Have a look at the [pipeline website](https://github.com/nf-core/ddamsproteomics) to find out how.

If you have problems that are related to Nextflow and not our pipeline then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).
