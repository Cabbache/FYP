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
./build/renderer -s scene.json
conversion
