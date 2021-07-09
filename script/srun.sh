#!/bin/sh

# runs benchmark using slurm's job system. No huge automation, just some
# defaults.

# stderr & stdout are redirected by default to caller's stderr & stdout.
srun \
	--partition training \
	--gres=gpu \
	--mem-per-cpu=4gb \
	--nodes=1 \
	--time=00:01:00 \
	--ntasks=1 \
	-N 1 \
	./srun.sh
