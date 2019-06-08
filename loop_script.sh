#!/bin/sh

ls testfiles/popsa |awk -F. '{print $1}'> name.lst

for x in `cat name.lst`
do
	profile=$x.profile
	popsa=$x.out
	mol=$(ls mol2/|grep  "$x"|awk '{print $1}')
	mol2=$(echo $mol|awk '{print $1}')
	mkdir -p output/$x
        python3 main.py -profile ./testfiles/profile/$profile -mol ./testfiles/mol2/$mol2 -popsa ./testfiles/popsa/$popsa -out ./output/$x -colorby properties -alpha 0.8 -size 256 -dpi 256
done

rm name.lst
