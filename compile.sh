#!/bin/bash
g++ main.cpp vec3.h tiny_obj_loader.h -march=native -fopenmp -Ofast -o main
ulimit -s unlimited
echo use command 'ulimit -s unlimited' before running main
