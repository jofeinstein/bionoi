#!/bin/bash 
#PBS -l nodes=1:ppn=20
#PBS -l walltime=6:00:00
#PBS -A loni_bionoi01
#PBS -q workq 
#PBS -N bionoi
#PBS -j oe 

export DAT_PATH=/work/$USER/bionoi_project

cd $DAT_PATH/

cat $PBS_NODEFILE | grep qb | sort | uniq > pssh-bionoi.lst

export PYTHONPATH=/project/michal/apps/pssh-2.3.1

/project/michal/apps/pssh-2.3.1/bin/pssh -h pssh-bionoi.lst -t 7200 "rm -rf /var/scratch/$USER ; mkdir -p /var/scratch/$USER"
/project/michal/apps/pssh-2.3.1/bin/pssh -h pssh-bionoi.lst -t 7200 "tar -xzf /work/$USER/bionoi_project/bae-data-mol2.tar.gz -C /var/scratch/$USER/"
/project/michal/apps/pssh-2.3.1/bin/pssh -h pssh-bionoi.lst -t 7200 "tar -xzf /work/$USER/bionoi_project/bae-pops.tar.gz -C /var/scratch/$USER/"
/project/michal/apps/pssh-2.3.1/bin/pssh -h pssh-bionoi.lst -t 7200 "tar -xzf /work/$USER/bionoi_project/bae-profile.tar.gz -C /var/scratch/$USER/"
/project/michal/apps/pssh-2.3.1/bin/pssh -h pssh-bionoi.lst -t 7200 "tar -xzf /work/$USER/bionoi_project/bionoi.tar.gz -C /var/scratch/$USER/"

rm pssh-bionoi.lst

export JOBS_PER_NODE=5

# parallel command options

module load gnuparallel/20170122

PARALLEL="parallel -j $JOBS_PER_NODE --slf $PBS_NODEFILE --wd /work/$USER/bionoi_project --joblog runtask.log"

# gnu-parallel launch serial tasks

$PARALLEL -a /work/$USER/bionoi_project/bionoi.lst sh /work/$USER/bionoi_project/bionoi_worker.sh {}
