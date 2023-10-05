#!/usr/bin/env bash
DATE=$(date +"%Y-%m-%d-%H:%M:%S")
echo "========================================================
Log of the pipeline run is written to 'sisterporec_$DATE.log' and can be observed in this terminal.
It will contain a message when finished. This terminal can be closed.
========================================================" > sisterporec_$DATE.log

nohup snakemake --keep-going --scheduler greedy --use-conda --use-singularity --rerun-incomplete --cluster "sbatch -t {cluster.time} -c {cluster.nodes} --mem {cluster.memory} -o slurm-out/%j.out"\
          --cluster-config config/cluster_config.yml --jobs 300 --output-wait 60 >> sisterporec_$DATE.log 2>&1 &

tail -n 1000 -f sisterporec_$DATE.log
