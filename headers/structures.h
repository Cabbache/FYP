#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "vec3.h"
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

typedef struct boxHitInfo{
	vec3 ray;
	vec3 origin;

	bool hit;
	float tmin;
	float tmax;
} boxHitInfo;

typedef struct Triangle{
	vec3 p[3];

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

typedef struct hitInfo{
	vec3 ray;
	vec3 origin;

	bool hit;
	float t;
	Triangle tri;

	float beta;
	float gamma;
} hitInfo;

typedef unordered_map<vec3_int, unordered_set<Triangle, Triangle::HashFunction>, hashFuncVec, equalsFunc> GridMap;

typedef struct Volume{
	vec3 min;
	vec3 max;
} Volume;

typedef struct SDF{
	vec3 origin;
	vec3 corner;
	vec3 dimensions;
	float resolution;
	float ***values;
} SDF;

////////////////////////////////////////////////////
// const.
////////////////////////////////////////////////////

const float KD_TREE_EPSILON = 0.00001f;

////////////////////////////////////////////////////
// enums.
////////////////////////////////////////////////////

enum SplitAxis {
	X_AXIS = 0,
	Y_AXIS = 1,
	Z_AXIS = 2
};

enum AABBFace {
	LEFT = 0,
	FRONT = 1,
	RIGHT = 2,
	BACK = 3,
	TOP = 4,
	BOTTOM = 5
};

////////////////////////////////////////////////////
// structs.
////////////////////////////////////////////////////

struct Ray
{
	vec3 origin;
	vec3 dir;
};


////////////////////////////////////////////////////
// classes.
////////////////////////////////////////////////////

class KDTreeNode
{
public:
	KDTreeNode( void );
	~KDTreeNode( void );

	Volume bbox;
	KDTreeNode *left;
	KDTreeNode *right;
	int num_tris;
	int *tri_indices;

	SplitAxis split_plane_axis;
	float split_plane_value;

	bool is_leaf_node;

	// One rope for each face of the AABB encompassing the triangles in a node.
	KDTreeNode *ropes[6];

	int id;

	bool isPointToLeftOfSplittingPlane( const vec3 &p ) const;
	KDTreeNode* getNeighboringNode( vec3 p );

	// Debug method.
	void prettyPrint( void );
};

//#include "KDTreeCPU.h"
class KDTreeCPU;

typedef struct Material{
	bool isLight;
	float specularity; //0 to 1
	float absorption; //0 to 1
	vec3 color;
} Material;

typedef struct Obj{
	Volume bounds;
	KDTreeCPU *kdtree;
	SDF sdf;
	Material material;
	struct {
		float resolution;
		GridMap map;
	} grid;
} Obj;


//used for the function that returns sdf value given world coordinate
//if world coordinate is outside sdf domain, inside = false
typedef struct SDFResult{
	bool inside;
	float value;
} SDFResult;

typedef pair<const Obj*, float> Ohd;
#endif
