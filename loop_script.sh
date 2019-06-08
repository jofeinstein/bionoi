#!/bin/sh

ls protein_files/bae-pops |awk -F. '{print $1}'> name.lst

for x in `cat name.lst`
do
	profile=$x.profile
	popsa=$x.out
	mol=$(ls ./protein_files/bae-data-mol2/|grep  "$x"|awk '{print $1}')
	mol2=$(echo $mol|awk '{print $1}')
	mkdir -p output/$x
	python3 main.py -mol ./protein_files/bae-data-mol2/$mol2 -pop ./protein_files/bae-pops/$popsa -profile ./protein_files/bae-profile/$profile -out ./output/$x -alpha 0.8 -size 256 -dpi 256 -colorby properties
done

rm name.lst
