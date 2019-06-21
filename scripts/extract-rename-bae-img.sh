find bae-images | shuf -n 5|cut -b 12-18 > bionoi_data.lst
 
for x in `cat bionoi_data.lst`
do
	mkdir -p /home/jfeinst/Projects/bionoi_files/tmp
	tar -xzvf /home/jfeinst/Projects/bionoi_files/bae-images/$x.tar.gz -C /home/jfeinst/Projects/bionoi_files/
	cd /home/jfeinst/Projects/bionoi_files/$x
	for f in * ;  do mv "$f" "$(basename "$(pwd)")"_"$f" ;  done
	mv /home/jfeinst/Projects/bionoi_files/$x/* /home/jfeinst/Projects/bionoi_files/test
	cd ..
	rm -rf $x
done

rm bionoi_data.lst
