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
g++ main.cpp vec3.h tiny_obj_loader.h -march=native -fopenmp -Ofast -o main

#allow bigger stack size, otherwise there will be segfault
ulimit -s unlimited

./main

conversion
