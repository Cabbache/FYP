#!/bin/bash
for file in *.ppm;do ./ppmtopng.sh $file;done
cat $(ls *.ppm.png | sort -V) | ffmpeg -framerate 24 -i - -c:v libx264 animation.mp4
