#ifndef MESH_H
#define MESH_H

#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "vec3.h"
#include "structures.h"

class mesh{

public:
	mesh();
	mesh(std::string fileName);
	std::vector<std::string> tokenizeString(std::string str);
	std::istream& safeGetline(std::istream& is, std::string& t);
	void translate(vec3);
	void scale(float);

	int numTris, numVerts;
	vec3* tris;
	vec3* verts;
	Volume bb;
};

#endif
