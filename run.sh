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

#compile
./compile.sh

#allow bigger stack size, otherwise there will be segfault
ulimit -s unlimited

./main

conversion
