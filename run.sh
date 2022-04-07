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
#SCENE="bunny_kitten.json"
SCENE="dragons.json"
if [[ "$1" == "quality" ]]
then
	./build/renderer -s "scenes/$SCENE" -w 1280 -h 960 -sw 1664 -sh 1248 -a 2
elif [[ "$1" == "fast" ]]
then
	./build/renderer -s "scenes/$SCENE" -w 300 -h 300 -sw 1664 -sh 1248 -d 0.5 -a 2 -sg 1.8
else
	./build/renderer -s "scenes/$SCENE" -w 150 -h 150 -d 0.5 -a 1
fi
conversion
