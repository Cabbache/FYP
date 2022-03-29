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
if [[ "$1" == "fast" ]]
then
	./build/renderer -s scene.json -w 150 -h 150 -d 0.5 -a 1
else
	./build/renderer -s scene.json -d 0.5 -w 1664 -h 1248 -sw 1664 -sh 1248
fi
conversion
