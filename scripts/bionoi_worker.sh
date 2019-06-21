#!/bin/bash

INP_DAT=$1

TMP_DIR=/var/scratch/$USER/tmp-$(echo $INP_DAT | cut -b 1-7)-$RANDOM

mkdir -p $TMP_DIR

cd $TMP_DIR/

for PDB_ID in `echo $INP_DAT | sed s/\:/\ /g`
do
 
 mkdir $PDB_ID
 
 /work/$USER/anaconda3/bin/python /var/scratch/$USER/bionoi/main.py -mol /var/scratch/$USER/bae-data-mol2/$PDB_ID.mol2 -pop /var/scratch/$USER/bae-pops/$(echo $PDB_ID | cut -b 1-5).out -profile /var/scratch/$USER/bae-profile/$(echo $PDB_ID | cut -b 1-5).profile -out $PDB_ID -colorby properties
 
 tar -cf $PDB_ID.tar $PDB_ID
 
 gzip -9 $PDB_ID.tar
 
 rm -rf $PDB_ID
 
 mv $PDB_ID.tar.gz /work/$USER/bionoi_project/bae-images/
 
done

cd $TMP_DIR/../

rm -rf $TMP_DIR

exit 0
