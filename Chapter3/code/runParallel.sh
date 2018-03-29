#!/bin/bash
### script to run an mpi job using 12-core or less (using only one 12-core node)
### Set the job name
#PBS -N testingNicheModels

### Specify the PI group for this job
### List of PI groups available to each user can be found with "va" command
#PBS -W group_list=neomartinez

### Set the queue for this job as windfall
#PBS -q standard

### jobtype parameter determines initial queue placement
###PBS -l jobtype=serial
### Set the number of nodes, cores and memory that will be used for this job
#PBS -l select=1:ncpus=28:mem=168gb:pcmem=6gb

### Specify "wallclock time" required for this job, hhh:mm:ss
#PBS -l walltime=24:00:00
### Specify total cpu time required for this job, hhh:mm:ss
### total cputime = walltime * ncpus
#PBS -l cput=672:00:00

#PBS -m bea
#PBS -M guj6y@email.arizona.edu

### Load required modules/libraries if needed (openmpi example)
### Use "module avail" command to list all available modules
### NOTE: /usr/share/Modules/init/csh -CAPITAL M in Modules
###     source /usr/share/Modules/init/csh


cd ~/Dissertation/Chapter3/code

module load matlab/r2016b

### run your executable program with begin and end date and time output
date
matlab -nodisplay -nodesktop -nosplash < NicheModelTests.m > ./run.log
date

