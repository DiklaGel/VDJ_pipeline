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
gelSeq is written in Python and so can just be downloaded, and run with `python3.5 gelSeq.py`.
Download the latest version and accompanying files from https://github.com/DiklaGel/VDJ_pipeline. 
gelSeq relies on several additional tools and Python modules that you should install.

### Pre-requisites

#### Software 
1. [IgBLAST](http://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone) - required for analysis of assembled contigs. (ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/).
2. [makeblastdb](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ ) 

##### Installing IgBlast 
Downloading the executable files from `ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/<version_number>` is not sufficient for a working IgBlast installation. You must also download the `internal_data` directory (ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/internal_data) and put it into the same directory as the igblast executable. This is also described in the igblast README file.

You should also ensure to set the `$IGDATA` environment variable to point to the location of the IgBlast executable. For example run `export IGDATA=/<path_to_igblast>/igblast/1.4.0/bin`.

gelSeq uses a configuration file to point it to the locations of files that it needs and a couple of other options.
An example configuration file is included in the repository - `gelSeq.conf`.

gelSeq looks for the configuration file, in descending order of priority, from the following sources:
1. The `-c` option used at run time for any of the gelSeq modes
2. The default global location of `~/.gelSeqrc`
3. If all else fails, the provided example `gelSeq.conf` in the directory is used

**Important:** If you  specify relative paths in the config file these will be used as relative to the main installation directory. For example, `resources/Mmus/igblast_dbs` will resolve to `/<wherever you installed gelSeq>/gelSeq/resources/Mmus/igblast_dbs`.

### External tool locations 
gelSeq will look in your system's `PATH` for external tools. You can override this behaviour by editing your `gelSeq.conf`.
Edit `gelSeq.conf` (or a copy) so that the paths within the `[tool_locations]` section point to the executables for all of the required tools.

	[tool_locations]
	#paths to tools used by gelSeq for alignment, quantitation, etc
	igblast_path = /path/to/igblastn

#### IgBLAST options 
##### Receptor type 
    igblast_seqtype = TCR

Type of sequence to be analysed. 

## Using gelSeq 
gelSeq has two modes: *plate* and *cell*

#### *plate*: Process fastq files from single plate, split the reads by cell basrcodes run gelSeq on each cell

##### Usage:
    gelSeq.py plate [-h] [--ncores <CORES>] [--config_file <CONFIG_FILE>]
                 [--resume_with_existing_files] [--species {Mmus,Hsap}]
                 [--receptor_name RECEPTOR_NAME] [--loci [LOCI [LOCI ...]]]
                 [--full] [--filter FILTER]
                 <FASTQ1> <FASTQ2> <PLATE_NAME> <OUTPUT_DIR>


##### Positional arguments:
    <FASTQ1>              first fastq file - read1
    <FASTQ2>              second fastq file - read2
    <PLATE_NAME>          name of plate for file labels
    <OUTPUT_DIR>          directory for output as <output_dir>/<plate_name>

##### Optional arguments:
    -h, --help            show this help message and exit
    --ncores <CORES>, -p <CORES> number of processor cores to use (default: 1)
    --config_file <CONFIG_FILE>, -c <CONFIG_FILE> config file to use (default: ~/.gelseqrc)
    --resume_with_existing_files, -r look for existing intermediate files and use those instead of starting from scratch (default: False)
    --species {Mmus,Hsap}, -s {Mmus,Hsap} Species to use for reconstruction (default: Hsap)
    --receptor_name RECEPTOR_NAME Name of receptor to reconstruct (default: TCR)
    --loci [LOCI [LOCI ...]]
                            Space-separated list of loci to reconstruct for
                            receptor (default: ['A', 'B'])
    --full                Continue the full process - after splitting to cells,
                            create new job for each cell (default: False)
    --filter FILTER       umis with more than filter reads (with respect to quantile) will be saved (default: 0.96)
 
#### *cell*: Reconstruct TCR sequences from RNAseq reads for a single cell
  
##### Usage:
    gelSeq.py cell [-h] [--ncores <CORES>] [--config_file <CONFIG_FILE>]
                 [--resume_with_existing_files] [--species {Mmus,Hsap}]
                 [--receptor_name RECEPTOR_NAME] [--loci LOCI [LOCI ...]]
                 <FASTA> <CELL_NAME> <OUTPUT_DIR>

##### Positional arguments:
    <FASTA>               fasta file
    <CELL_NAME>           name of cell for file labels
    <OUTPUT_DIR>          directory for output as <output_dir>/<cell_name>

##### Optional arguments:
    -h, --help            show this help message and exit
    --ncores <CORES>, -p <CORES> 
        number of processor cores to use (default: 1)
    --config_file <CONFIG_FILE>, -c <CONFIG_FILE>
                        config file to use (default: ~/.gelseqrc)
    --resume_with_existing_files, -r
                        look for existing intermediate files and use those instead of starting from scratch (default: False)
    --species {Mmus,Hsap}, -s {Mmus,Hsap}
                        Species to use for reconstruction (default: Hsap)
    --receptor_name RECEPTOR_NAME
                        Name of receptor to reconstruct (default: TCR)
    --loci LOCI [LOCI ...]
                        Space-separated list of loci to reconstruct for receptor (default: ['A', 'B'])


