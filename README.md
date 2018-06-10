# gelSeq
Reconstruction of T cell receptor sequences from single-cell RNA-seq data.

## Contents ##
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Setup](#setup)
4. [Usage](#using-gelSeq)
	- [*plate*](#split fastq files to cells and run gelSeq on each cell )
    - [*cell*](#run IgBlast and output statistics about the most abundant VDJ sequence)


## Introduction
This tool was build to process single-cell RNA-seq data of TCRB gene and return an output with statistics about the most abundant VDJ sequence of each cell. 

The code is based on TraCer (https://github.com/Teichlab/tracer) architecture with adjustments to another version of input.


## Installation
gelSeq is written in Python and so can just be downloaded, made executable (with `chmod u+x tracer`) and run or run with `python gelSeq`.
Download the latest version and accompanying files from https://github.com/DiklaGel/VDJ_pipeline. 
gelSeq relies on several additional tools and Python modules that you should install.

### Pre-requisites

#### Software 
1. [IgBLAST](http://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone) - required for analysis of assembled contigs. (ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/).
2. [makeblastdb](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ ) - **optional** but required if you want to use TraCeR's `build` mode to make your own references.

##### Installing IgBlast 
Downloading the executable files from `ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/<version_number>` is not sufficient for a working IgBlast installation. You must also download the `internal_data` directory (ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data) and put it into the same directory as the igblast executable. This is also described in the igblast README file.

You should also ensure to set the `$IGDATA` environment variable to point to the location of the IgBlast executable. For example run `export IGDATA=/<path_to_igblast>/igblast/1.4.0/bin`.

gelSeq uses a configuration file to point it to the locations of files that it needs and a couple of other options.
An example configuration file is included in the repository - `gelSeq.conf`.

TraCeR looks for the configuration file, in descending order of priority, from the following sources:
1. The `-c` option used at run time for any of the gelSeq modes
2. The default global location of `~/.gelSeqrc`
3. If all else fails, the provided example `tracer.conf` in the directory is used

**Important:** If you  specify relative paths in the config file these will be used as relative to the main installation directory. For example, `resources/Mmus/igblast_dbs` will resolve to `/<wherever you installed tracer>/tracer/resources/Mmus/igblast_dbs`.

### External tool locations 
gelSeq will look in your system's `PATH` for external tools. You can override this behaviour by editing your `gelSeq.conf`.
Edit `gelSeq.conf` (or a copy) so that the paths within the `[tool_locations]` section point to the executables for all of the required tools.

	[tool_locations]
	#paths to tools used by TraCeR for alignment, quantitation, etc
	igblast_path = /path/to/igblastn

#### IgBLAST options 
##### Receptor type 
    igblast_seqtype = TCR

Type of sequence to be analysed. 

## Using gelSeq 
gelSeq has two modes: *plate* and *cell*