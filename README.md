# Sister-Pore-C

This repository contains an adaption of the [Pore-C preprocessing pipeline](https://github.com/nanoporetech/Pore-C-Snakemake) to be able to deal with sister chromatid specific contacts. Specifically, data suitable for this pipeline is labelled with the synthetic nucleotide analogue BrdU on one strand of each sister chromatid.

## Installation

To install the sister-pore-c pipeline, you only need to clone this repository using:

```
git clone https://github.com/nanoporetech/Pore-C-Snakemake.git
```

To create the conda environment for running the pipeline, change into the directory of it and run:

```
conda env create
```

To activate it, run

```
conda activate pore_c_snakemake
```


## Usage

### Config files

To specify the input files, genome assemblies for alignment and output folder, adjustments to a set of config files are needed:

#### **basecalls.tsv**

| **run_id**  | **enzyme**       | **refgenome_ids**           | **biospecimen** | **fastq_path**                 | **fast5_directory**     | **sequencing_summary_path**     |
|-------------|------------------|-----------------------------|-----------------|--------------------------------|-------------------------|---------------------------------|
| Sample name | Digestion enzyme | Id for the reference genome | -               | Path to basecalled fastq files | Path to fast5 directory | Path to sequencing summary path |

This files contains paths to input files. To run the pipeline, the following files are needed:

-  The path to the basecall fastq files
- The path to the fast5 directory
- The path to the sequencing summary path (mapping of read_ids to fast5 files)

In addition, each row should contain a `run_id`, which specifies the name of a given sample (**Note**: `run_id` cannot contain the characters `_` as this is needed for path construction and parsing ). The enzyme field should contain the name of a restriction enzyme that is recognized by biopython.

#### **config.yaml**

This file contains numerous parameters for parts of the pipeline, most of which do not need to be changed. The only parameter that needs changing is the `output_dir` parameter that specifies the output directory.

#### **references.tsv**

This file contains the paths and ids for the reference genomes to be used. In most cases, this fill will contain a single (uncommented) line of the form:

```
refgenome_ids   PATH
```

Note that the id of the reference genome needs to correspond to the id used in `basecalls.tsv`.

### Pipeline submission


To run the pipeline, active the `pore_c_snakemake` environment (see above) and run the submission script (`submit_snakemake_job.sh`):

```
conda activate pore_c_snakemake
./submit_snakemake_job.sh
```

This will start the pipeline in the current terminal process. If you want to be able to run the pipeline in the background and come back to inspecting its outputs, I recommend using a terminal multiplexer such as [tmux](https://github.com/tmux/tmux/wiki). This allows you to switch quickly between different processes and detach from them. The snakemake process that orchestrates the pipeline is not resource heavy and can (usually) be run on the login node of HPC clusters. To start a tmux session, run

```
tmux new -s {SESSION_NAME}
```
If you want to detach again, hits CTRL-b d

and to reattach run

```
tmux a -t {SESSION_NAME}
```

## High-level analysis steps

## Specific task dependencies

## Troubleshooting

