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
SPATH="cubes_path.json"
#SCENE="dragon_simple.json"
#SCENE="bunnies.json"
#SCENE="bunnies_simple.json"
#SCENE="cornell.json"
SCENE="cube.json"
#SCENE="planes.json"
#SCENE="sibenik.json"
#SCENE="bunny_kitten.json"
#SCENE="dragons.json"
if [[ "$1" == "quality" ]]
then
	./build/renderer -s "scenes/$SCENE" -w 1280 -h 960 -sw 1664 -sh 1248 -a 2 -ppm
elif [[ "$1" == "fast" ]]
then
	./build/renderer -s "scenes/$SCENE" -w 300 -h 300 -sw 1664 -sh 1248 -a 5 -sg 1.8
elif [[ "$1" == "path" ]]
then
	#./build/renderer -s "scenes/$SCENE" -w 1280 -h 960 -a 2 -sg 1.8 -p "$SPATH" -f 120
	./build/renderer -s "scenes/$SCENE" -w 300 -h 300 -a 1 -sg 1.8 -p "$SPATH" -f 30
	#./build/renderer -s "scenes/$SCENE" -w 30 -h 30 -a 1 -sg 1.8 -p "$SPATH" -f 30
else
	./build/renderer -s "scenes/$SCENE" -w 150 -h 150 -a 1 -sg 1.3
fi
conversion
