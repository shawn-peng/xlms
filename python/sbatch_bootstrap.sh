#!/bin/bash

#sbatch cpu_batch_run.sh gaussian

bootstrap_start=$1
bootstrap_end=$2

for i in $(seq $bootstrap_start $bootstrap_end);
do
	echo $i;
	sbatch cpu_batch_run_bootstrap.sh gaussian $i;
done

