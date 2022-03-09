#ifndef KD_TREE_STRUCTS_H
#define KD_TREE_STRUCTS_H

#include "vec3.h"

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

class KDTreeNodeGPU
{
public:
	KDTreeNodeGPU( void );

	Volume bbox;
	int left_child_index;
	int right_child_index;
	int first_tri_index;
	int num_tris;

	int neighbor_node_indices[6];

	SplitAxis split_plane_axis;
	float split_plane_value;

	bool is_leaf_node;

	bool isPointToLeftOfSplittingPlane( const vec3 &p ) const;
	int getNeighboringNodeIndex( vec3 p );

	// Debug method.
	void prettyPrint( void );
};

#endif
