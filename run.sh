#!/bin/bash

#requires imagemagick
function conversion(){
	convert img_$1.ppm img_$1.png
	rm img_$1.ppm
}

set -e

#compile
g++ main.cpp vec3.h tiny_obj_loader.h -march=native -fopenmp -Ofast -o main

#allow bigger stack size, otherwise there will be segfault
ulimit -s unlimited

for X in `seq 0 3 360`
do
	./main 3 $X > img_$X.ppm
	conversion $X &
	echo $X
done
