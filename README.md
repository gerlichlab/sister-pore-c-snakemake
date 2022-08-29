# Sister-Pore-C

This repository contains an adaption of the [Pore-C preprocessing pipeline](https://github.com/nanoporetech/Pore-C-Snakemake) to be able to deal with sister chromatid specific contacts. Specifically, data suitable for this pipeline is labeled with the synthetic nucleotide analog BrdU on one strand of each sister chromatid.

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


## Config files

To specify the input files, genome assemblies for alignment, and output folder, adjustments to a set of config files are needed:

### **basecalls.tsv**

| **run_id**  | **enzyme**       | **refgenome_ids**           | **biospecimen** | **fastq_path**                 | **fast5_directory**     | **sequencing_summary_path**     |
|-------------|------------------|-----------------------------|-----------------|--------------------------------|-------------------------|---------------------------------|
| Sample name | Digestion enzyme | Id for the reference genome | -               | Path to basecalled fastq files | Path to fast5 directory | Path to sequencing summary path |

This file contains paths to input files. To run the pipeline, the following files are needed:

-  The path to the basecall fastq files
- The path to the fast5 directory
- The path to the sequencing summary path (mapping of read_ids to fast5 files)

In addition, each row should contain a `run_id`, which specifies the name of a given sample (**Note**: `run_id` cannot contain the characters `_` as this is needed for path construction and parsing ). The enzyme field should contain the name of a restriction enzyme that is recognized by biopython.

### **config.yaml**

This file contains numerous parameters for parts of the pipeline, most of which need not be changed. The only parameter that needs changing is the `output_dir` parameter, that specifies the output directory.

### **references.tsv**

This file contains the paths and ids for the reference genomes to be used. In most cases, this fill will contain a single (uncommented) line of the form:

```
refgenome_ids   PATH
```

Note that the id of the reference genome needs to correspond to the id used in `basecalls.tsv`.

## Pipeline submission


To run the pipeline, active the `pore_c_snakemake` environment (see above) and run the submission script (`submit_snakemake_job.sh`):

```
conda activate pore_c_snakemake
./submit_snakemake_job.sh
```

This will start the pipeline in the current terminal process. If you want to be able to run the pipeline in the background and come back to inspect its outputs, I recommend using a terminal multiplexer such as [tmux](https://github.com/tmux/tmux/wiki). This allows you to switch quickly between different processes and detach from them. The snakemake process that orchestrates the pipeline is not resource-heavy and can (usually) be run on the login node of HPC clusters. To start a tmux session, run

```
tmux new -s {SESSION_NAME}
```
If you want to detach again, hits CTRL-b d

and to reattach run

```
tmux a -t {SESSION_NAME}
```
**Attention: Currently, there seems to be an issue when entering into tmux with an activated conda environment. Snakemake will have issues switching environments for the individual steps. Deactivate conda before entering the tmux session. Once inside, reactivate it.**

## Analysis steps

### Pairs-based analysis

This pipeline extends the pore-c snakemake pipeline and builds analysis steps on top of its output. Specifically, the following steps are performed:

- BrdU calls are generated by using [DNAscent v2](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07736-6). To this end, first a BrdU index is created based on the sequencing summary file and the fast5 input. Then, DNAscent is run using the generated index, the aligned bam files and the reference genome fasta. The environment to run this analysis is in a docker container that is publicly available at `docker://gerlichlab/dnascent-docker:latest`. The output of this step is a `*.detect` file for each bam input that contains BrdU-probabilities for each base of each read. Note that at this step, the analysis is run on a batch-level to parallelize efficiently. Each sample is split into batches of 100k reads (see the original pipeline for a more detailed explanation).
- The next step is to generate a "label library". A label library in this context is a mapping of read-IDs to their respective (binary) labeling state. This step takes as input all the detect files of a given sample and generates a pickle file containing the label library as a python dictionary. On a technical note, the key of this dictionary is not the actual read id since a given read id may contain multiple aligned fragments. Instead, the key is composed of the read_id, alignment chromosome, alignment start and alignment end. This is assumed to be a unique identifier for aligned fragments.
- After generation of a label library, pairs are assigned their respective labeling state using the pairs output form the conventional porec-pipeline and the generated label library. The output of this step is a file with assigned pairs
- Then, pairs are split into cis, trans and combined contacts based on different rules. The `labelled_only` rule assigns sister-identity only based on labeled reads, whereas the `all_reads` rule assigns sister-identity based on all reads. The `labelled_only` rule is more conservative but results in a lower yield.
- After splitting of pairs, coolers are generated.

### Higher-order contact-based analysis

Pore-C allows analysis of higher-order contacts and this pipeline implements the generation of contact triplets with assigned sister identities. A contact triplet are three aligned fragments that have been sequenced on the same read. Generation of triplets is done analogously to generating duplets (or pairs) by taking all possible triplet combinations for each given read. The following pipeline steps are concerned with this task:

- First, triplets are expanded and saved as a triplet_contact parquet file
- Then, triplets for each basecall batch are merged
- Finally, sister identity is annotated using the merged triplet contacts and the label library mentioned above. The final output is a triplet parquet file.

## Troubleshooting
If you get an error like this when trying to run the pipeline:
 ```
MissingInputException in line 45 of /groups/gerlich/shared/fede-pore-c/sister-pore-c-snakemake/Snakefile:
Missing input files for rule all:
...
 ```
 You most likely have an error in one of the paths of the config/basecalls.tsv, or you are using an editor that replaces tabs with spaces. This will mess with the structure of the .tsv format.

If you run into ou of memory or timeout issues increase the values for your step in `config\cluster_config.yml` be careful not to increase them over the limit of your cluster. Since if slurm can not find the computational resources (e.g. a node with more RAM then currently exists) the job will not fail but slurm will be stuck waiting for those to appear magically.

## Step-by-step instructions for running the pipeline in the gerlich lab
Download the data from the facility into the ngs folder using the
wget command on the cluster:
```
cd /groups/ngs
wget -c --no-check-certificate --auth-no-challenge --user NAME.SURNAME --ask-password LINK_THAT_YOU_FIND_IN_FORSKALLE
```
Create in your experiment folder the following folders:
```
mkdir sequencing_data
mkdir sequencing_results
```
Pull this repository into your experiment folder:
```
git clone https://github.com/nanoporetech/Pore-C-Snakemake.git
```
Your experiment folder should now look like this:
```
├──EXPERIMENT_ID
    ├──sequencing_data
    ├──sequencing_results
    └──sister-pore-c-snakemake 
```
Unzip and copy the data from the ngs folder into the experiment folder, preferably not do this on the login node:
```
tar -xzvf /groups/gerlich/ngs/???.tgz -C /groups/gerlich/experiments/Experiments_???/???/sequencing_data/
```
Update `config\basecalls.tsv` to fit your experiment an example of a multiplexed PromethION flowcell can be found in experiment 005548.

Update the output folder in `config\config.yaml` if you do not want the results in the upstream `sequencing_results` folder we just created.

Go into tmux (you can choose any SESSION_NAME):
```
tmux new -s SESSION_NAME
```
Activate the environment (you will find above how to create it) and run the script:
```
conda activate pore_c_snakemake
./submit_snakemake_job.sh
```

Detach from tmux with CTRL-b d

Now wait a few hours, and then you reattach to tmux and check if it has finished successfully with: 

```
tmux a -t SESSION_NAME
```
