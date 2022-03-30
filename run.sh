#!/bin/bash

#requires imagemagick
function conversion(){
	for file in *.ppm
	do
		convert $file $file.png
		rm $file
	done
}

set -e

./compile.sh
rm img_*.ppm || :
if [[ "$1" == "quality" ]]
then
	./build/renderer -s scenes/bunny_kitten.json -w 1280 -h 960 -sw 1664 -sh 1248 -a 2
elif [[ "$1" == "fast" ]]
then
	./build/renderer -s scenes/bunny_kitten.json -w 300 -h 300 -d 0.5 -a 2
else
	./build/renderer -s scenes/bunny_kitten.json -w 150 -h 150 -d 0.5 -a 1
fi
conversion
