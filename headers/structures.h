#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "vec3.h"
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

typedef struct hitInfo{
	vec3 ray;
	vec3 origin;

	bool hit;
	double t;

	double beta;
	double gamma;
} hitInfo;

typedef struct boxHitInfo{
	vec3 ray;
	vec3 origin;

	bool hit;
	double tmin;
	double tmax;
} boxHitInfo;

typedef struct Triangle{
	vec3 p[3];
	bool reflective;

	bool operator==(const Triangle& triangle) const{
		return p[0] == triangle.p[0] && p[1] == triangle.p[1] && p[2] == triangle.p[2];
	}

	struct HashFunction{
		size_t operator()(const Triangle& triangle) const
    {
			hashFuncVec hfv;
			return hfv(triangle.p[0]) ^ hfv(triangle.p[1]) ^ hfv(triangle.p[2]);
    }
	};

} Triangle;

typedef unordered_map<vec3_int, unordered_set<Triangle, Triangle::HashFunction>, hashFuncVec, equalsFunc> GridMap;

typedef struct Volume{
	vec3 min;
	vec3 max;
} Volume;

typedef struct SDF{
	vec3 origin;
	vec3 corner;
	vec3 dimensions;
	double resolution;
	double ***values;
} SDF;

typedef struct Obj{
	Volume bounds;
	SDF sdf;
	struct {
		double resolution;
		GridMap map;
	} grid;
} Obj;

//used for the function that returns sdf value given world coordinate
//if world coordinate is outside sdf domain, inside = false
typedef struct SDFResult{
	bool inside;
	double value;
} SDFResult;

typedef pair<const Obj*, double> Ohd;

#endif
