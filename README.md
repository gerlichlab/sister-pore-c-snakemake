# Sister-Pore-C

This repository contains an adaption of the [Pore-C preprocessing pipeline](https://github.com/nanoporetech/Pore-C-Snakemake) to be able to deal with sister chromatid specific contacts. Specifically, data suitable for this pipeline is labelled with the synthetic nucleotide analogue BrdU on one strand of each sister chromatid.

## Installation

To install the sister-pore-c pipeline, you only need to clone this repository using:

```
git clone https://github.com/gerlichlab/sister-pore-c-snakemake.git
```

To create the conda environment for running the pipeline, change into the directory of it and run:

```
conda env create
```

To activate it, run

```
conda activate pore-c-snakemake
```


## Config files

To specify the input files, genome assemblies for alignment, and output folder, adjustments to a set of config files are needed:

### **basecalls.tsv**

| **run_id**  | **enzyme**       | **refgenome_ids**           | **biospecimen** | **fastq_path**                 | **fast5_directory**     | **sequencing_summary_path**     |
|-------------|------------------|-----------------------------|-----------------|--------------------------------|-------------------------|---------------------------------|
| Sample name | Digestion enzyme | Id for the reference genome | -               | Path to basecalled fastq files | Path to fast5 directory | Path to sequencing summary path |

This file contains paths to input files. To run the pipeline, the following files are needed:

- The path to the basecall fastq files
- The path to the fast5 directory
- The path to the sequencing summary path (mapping of read_ids to fast5 files)

In addition, each row should contain a `run_id`, which specifies the name of a given sample (**Note**: `run_id` cannot contain the characters `_` as this is needed for path construction and parsing ). The enzyme field should contain the name of a restriction enzyme that is recognized by biopython.

### **config.yaml**

This file contains numerous parameters for parts of the pipeline, most of which need not be changed. It contains the batch size as a parameter, per default set to 10k reads; setting it to low might stall the snakemake solver. The `output_dir` parameter is set per default relative to the current folder at: `../pipeline_results`

### **references.tsv**

This file contains the paths and ids for the reference genomes to be used. In most cases, you will use the hg19 snip corrected HELA reference genome created by Michael Mitter:

```
refgenome_ids   PATH
hg19    /groups/gerlich/members/MichaelMitter/Reference_genomes/Fasta/hg19_SNPs/hg19_SNPs.fa
```

Note that the id of the reference genome needs to correspond to the id used in `basecalls.tsv`.

### **cluster_config.yml**
Contains the time, cores and memory requested for the different jobs. It is set only to request the minimal amount necessary. In some instances, memory and time need to be increased as slurm (the cluster schedular) will terminate a job if the requested time or memory is exceeded. 

## Pipeline submission


To run the pipeline, active the `pore_c_snakemake` environment (see above) and run the submission script (`submit_snakemake_job.sh`):

This will start snakemake with the `nohub` command (this allows you to log out and keep snakemake running in the background) and pipes the outputs to `sisterporec_$DATE.log`. It will keep displaying the progress of the log file. You can press Ctrl-c to stop the displaying, this will not stop snakemake, it will keep running in the background. To stop snakemake you will need to use htop. If you want to resume the live viewing, you can use the command:
`tail -n 1000 -f MOST_RECENT_LOGFILE`

## Analysis steps

### Pairs-based analysis

This pipeline extends the pore-c snakemake pipeline and builds analysis steps on top of its output. Specifically, the following steps are performed:

- BrdU calls are generated by using [DNAscent v2](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07736-6). To this end, first, a BrdU index is created based on the sequencing summary file and the fast5 input. Then, DNAscent is run using the generated index, the aligned bam files and the reference genome fasta. The environment to run this analysis is in a docker container that is publicly available at `docker://gerlichlab/dnascent-docker:latest`. The output of this step is a `*.detect` file for each bam input that contains BrdU-probabilities for each base of each read. Note that at this step, the analysis is run on a batch-level to parallelize efficiently. Each sample is split into batches of 100k reads (see the original pipeline for a more detailed explanation).
- The next step is to generate a "label library". A label library in this context is a mapping of read-IDs to their respective (binary) labelling state. This step takes as input all the detect files of a given sample and generates a pickle file containing the label library as a python dictionary. On a technical note, the key of this dictionary is not the actual read id since a given read id may contain multiple aligned fragments. Instead, the key is composed of the read_id, alignment chromosome, alignment start and alignment end. This is assumed to be a unique identifier for aligned fragments.
- After the generation of a label library, pairs are assigned their respective labelling state using the pairs output from the conventional porec-pipeline and the generated label library. The output of this step is a file with assigned pairs
- Then, pairs are split into cis, trans and combined contacts based on different rules. The `labelled_only` rule assigns sister-identity only based on labelled reads, whereas the `all_reads` rule assigns sister-identity based on all reads. The `labelled_only` rule is more conservative but results in a lower yield.
- After splitting of pairs, coolers are generated.

### Higher-order contact-based analysis

Pore-C allows analysis of higher-order contacts and this pipeline implements the generation of contact triplets with assigned sister identities. A contact triplet are three aligned fragments that have been sequenced on the same read. Generation of triplets is done analogously to generating duplets (or pairs) by taking all possible triplet combinations for each given read. The following pipeline steps are concerned with this task:

- First, triplets are expanded and saved as a triplet_contact parquet file
- Then, triplets for each basecall batch are merged
- Finally, sister identity is annotated using the merged triplet contacts and the label library mentioned above. The final output is a triplet parquet file.

## Troubleshooting
- If you get an error like this when trying to run the pipeline:
 ```
InputFunctionException in rule import_basecalls in file /groups/gerlich/sequencing_data/155/sister-pore-c-snakemake/rules/reads.smk, line 1:
Error:
  KeyError: '319165'
Wildcards:
  enzyme=NlaIII
  run_id=319165
 ```
 Please use a run_id that does not contain a pure number! (i.e. add a letter at the beginning a319165)
 
 - If you get an error like this when trying to run the pipeline:
    ```
    MissingInputException in line 45 of /groups/gerlich/shared/fede-pore-c/sister-pore-c-snakemake/Snakefile:
    Missing input files for rule all:
    ...
    ```
    You most likely have an error in one of the paths of the config/basecalls.tsv, or you are using an editor that replaces tabs with spaces. This will mess with the structure of the .tsv format.
- If you run into out-of-memory or timeout issues, increase the values for your step in `config\cluster_config.yml`. Be careful not to increase them over the limit of your cluster. Since if slurm can not find the computational resources (e.g. a node with more RAM than currently exists), the job will not fail, but slurm will be stuck waiting for those to appear magically.
- If you see in your `sister-pore-c-snakemake/sisterporec_????.log` something like
  ```
  Finished job 68.
  2 of 6 steps (33%) done
  ```
  and you do not have jobs queued anymore for your user (check with `squeue -u FIRSTNAME.LASTNAME`).
  One of your slurmjobs failed due to using too many resources and was killed without snakemake noticing.
  Kill the running snakemake via htop using `SIGKILL`.
  
  Update the `config\cluster_config.yml` as in the paragraph above.
  
  Then, unlock the `sister-pore-c-snakemake` folder with `snakemake --unlock` before you restart the pipeline.
  
  Tip: If you can not find which jobs failed, you can simply rerun the pipeline again without changing `config\cluster_config.yml` and the new run should start and fail again at the job that failed before. Now look at the new `sister-pore-c-snakemake/sisterporec_????.log` for the last submitted jobs e.g. `Submitted batch job 63213725` run the command `jobinfo 63213725` it will tell you what rule that job was and why the job failed. Now increase `config\cluster_config.yml` accordingly.
  
## Step-by-step instructions for running the pipeline in the gerlich lab
Download the data from the facility into the ngs folder using the
wget command on the cluster:
```
cd /groups/ngs
wget -c --no-check-certificate --auth-no-challenge --user NAME.SURNAME --ask-password LINK_THAT_YOU_FIND_IN_FORSKALLE
```
or preferable ask ngs facility for the full command with checksum, it should look like this:
```
wget -c --no-check-certificate --auth-no-challenge --user NAME.SURNAME --ask-password https://ngs.vbcf.ac.at/filemanager/byurl/a9c5298_20220802_1117_2A_PAK48337_ee3198ee.tar
md5sum -c <(echo "e81a8822fdf8803d63af88ad5cdd7ee1  a9c5298_20220802_1117_2A_PAK48337_ee3198ee.tar") 
```

Create in your experiment folder the following folders:
```
mkdir sequencing_data
mkdir pipeline_results
```
Pull this repository into your experiment folder:
```
git clone https://github.com/gerlichlab/sister-pore-c-snakemake.git
```
To unzip and copy the data from the ngs folder into the experiment folder, preferably not do this on the login node
create a file called to_unzip.sh that contains a line like this:
```
tar -xvf /groups/gerlich/ngs/???.tgz -C /groups/gerlich/sequencing_data/???/sequencing_data/
```
This will allow us to delete the unzipped data after the pipeline is finished but help us unzip the data if we need to rerun the pipeline again with the following command:
```
bash to_unzip.sh
```
Create a new `basecall.tsv` to fit your experiment; an example of a multiplexed PromethION flowcell can be found in `config\basecalls.tsv`.
Copy your basecall.tsv and overwrite the sample one in the config folder of the pipeline by:
```
cp basecalls.tsv sister-pore-c-snakemake/config/.
```
Last but not least, link the Experiment that either produced this data or were the follow-up analysis is taking place with the following command:
```
ln -s ../../experiments/Experiments_00??00/00????/ exp????
```
Your experiment folder should now look like this:
```
├──EXPERIMENT_ID
    ├──basecalls.tsv
    ├──exp????
    ├──pipeline_results
    ├──sequencing_data
    ├──sister-pore-c-snakemake
    └──to_unzip.sh
```
Now run the pipeline; running more than two at a time will not make it faster since you will be limited by the slurm quotas at our institute.
```
conda activate pore-c-snakemake
./submit_snakemake_job.sh
```
This will start snakemake with the `nohub` command (this allows you to log out and keep snakemake running in the background) and pipes the outputs to `sisterporec_$DATE.log`. It will keep displaying the progress of the log file. You can press Ctrl-c to stop the displaying, this will not stop snakemake, it will keep running in the background. To stop snakemake you will need to use htop. If you want to resume the live viewing, you can use the command:
`tail -n 1000 -f MOST_RECENT_LOGFILE`

If you see in the `sister-pore-c-snakemake/sisterporec_????.log` the two final lines:
```
? of ? steps (100%) done
Complete log: .snakemake/log/???.snakemake.log
```
you can now delete the content of `sequencing_data`; otherwise, check above for what to do if one of the slurmjob runs out of resources.

