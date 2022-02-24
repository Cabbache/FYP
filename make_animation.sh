#!/bin/bash
cat $(ls *.ppm.png | sort -V) | ffmpeg -framerate 24 -i - -c:v libx264 animation.mp4
