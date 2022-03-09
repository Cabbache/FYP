#!/bin/bash
g++ main.cpp mesh.cpp mesh.h KDTreeCPU.cpp KDTreeCPU.h KDTreeStructs.h KDTreeStructs.cpp vec3.h tiny_obj_loader.h -march=native -fopenmp -Ofast -o main
#g++ main.cpp vec3.h tiny_obj_loader.h -g -fopenmp -o main
