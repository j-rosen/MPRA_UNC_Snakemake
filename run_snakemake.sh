#!/bin/bash

# Command line submission for MPRAsnakeflow pipeline via Longleaf cluster
# Note the latency time of 60 seconds. Test runs at default (5 seconds)
# did not find necessary files in time.
module add python
module add r


# Directory named logs required for sbatch processing
if [[ ! -e logs ]]; then
  mkdir logs
fi


# SLURM command for execution
sbatch -t 8:00:00 -o MPRA_UNC_Snakemake.log --wrap "snakemake --latency-wait 60 --configfile config/config.yaml --cluster 'sbatch --nodes=1 --ntasks={cluster.threads} --mem={cluster.mem} -t {cluster.time} -p {cluster.queue} -o {cluster.output}' --jobs 45 --cluster-config config/sbatch.yaml"
