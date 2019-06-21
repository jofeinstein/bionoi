#!/bin/bash

PDB_ID=$1

export DAT_PATH=/home/jfeinst/Projects/bionoi

TMP_DIR=/var/scratch/jfeinst/tmp-$PDB_ID-$RANDOM

mkdir -p $TMP_DIR
mkdir -p $DAT_PATH/clusteroutput/$PDB_ID

cd $TMP_DIR/

M=`echo $PDB_ID | cut -b -5`

cp $DAT_PATH/protein_files/bae-pops/$M.out $TMP_DIR/
cp $DAT_PATH/protein_files/bae-profile/$M.profile $TMP_DIR/
cp $DAT_PATH/protein_files/bae-data-mol2/$PDB_ID.mol2 $TMP_DIR/


python3 $DAT_PATH/main.py -mol $PDB_ID.mol2 -pop $M.out -profile $M.profile -out $DAT_PATH/clusteroutput/$PDB_ID -colorby properties


cd $TMP_DIR/../

rm -rf $TMP_DIR

exit 0
