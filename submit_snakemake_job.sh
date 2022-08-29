snakemake --use-conda --use-singularity --rerun-incomplete --cluster "sbatch -t {cluster.time} -c {cluster.nodes} --mem {cluster.memory}"\
          --cluster-config config/cluster_config.yml --jobs 300 --output-wait 60
