#!/bin/bash
cat $(ls *.ppm.png | sort -V) | ffmpeg -framerate 24 -i - -c:v libx264 -pix_fmt yuv420p animation.mp4
