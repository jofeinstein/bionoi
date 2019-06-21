find bae-data-mol2 -type f | shuf -n 1000 | cut -b 15-26 > bionoi_data.lst

for x in `cat bionoi_data.lst`
do
	python3 main.py -mol ./bae-data-mol2/$x -direction 1 -rot_angle 1 -flip 1 -out ./dataset/
	mv ./dataset/XOY+_r0_OO.jpg ./dataset/$(echo $x|cut -b 1-7).jpg
done
