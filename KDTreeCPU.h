#ifndef KD_TREE_CPU_H
#define KD_TREE_CPU_H

#include <limits>
#include "KDTreeStructs.h"
#include "intersections.h"

////////////////////////////////////////////////////
// Constants.
////////////////////////////////////////////////////

//const int NUM_TRIS_PER_NODE = 2000;
const int NUM_TRIS_PER_NODE = 10;
const bool USE_TIGHT_FITTING_BOUNDING_BOXES = false;
//const float INFINITY = std::numeric_limits<float>::max();


////////////////////////////////////////////////////
// KDTreeCPU.
////////////////////////////////////////////////////
class KDTreeCPU
{
public:
	KDTreeCPU( int num_tris, vec3 *tris, int num_verts, vec3 *verts );
	~KDTreeCPU( void );

	// Public traversal method that begins recursive search.
	bool intersect(hitInfo &hitinfo, vec3 &hit_point) const;
	bool singleRayStacklessIntersect( const vec3 &ray_o, const vec3 &ray_dir, float &t, vec3 &hit_point) const;

	void buildRopeStructure( void );

	// kd-tree getters.
	KDTreeNode* getRootNode( void ) const;
	int getNumLevels( void ) const;
	int getNumLeaves( void ) const;
	int getNumNodes( void ) const;

	// Input mesh getters.
	int getMeshNumVerts( void ) const;
	int getMeshNumTris( void ) const;
	vec3* getMeshVerts( void ) const;
	vec3* getMeshTris( void ) const;

	// Debug methods.
	void printNumTrianglesInEachNode( KDTreeNode *curr_node, int curr_depth=1 );
	void printNodeIdsAndBounds( KDTreeNode *curr_node );

private:
	// kd-tree variables.
	KDTreeNode *root;
	int num_levels, num_leaves, num_nodes;

	// Input mesh variables.
	int num_verts, num_tris;
	vec3 *verts, *tris;

	KDTreeNode* constructTreeMedianSpaceSplit( int num_tris, int *tri_indices, Volume bounds, int curr_depth );

	// Private recursive traversal method.
	bool intersect( KDTreeNode *curr_node, hitInfo &hitinfo) const;
	bool singleRayStacklessIntersect( KDTreeNode *curr_node, const vec3 &ray_o, const vec3 &ray_dir, float &t_entry, float &t_exit) const;

	// Rope construction.
	void buildRopeStructure( KDTreeNode *curr_node, KDTreeNode *ropes[], bool is_single_ray_case=false );
	void optimizeRopes( KDTreeNode *ropes[], Volume bbox );

	// Bounding box getters.
	SplitAxis getLongestBoundingBoxSide( vec3 min, vec3 max );
	Volume computeTightFittingBoundingBox( int num_verts, vec3 *verts );
	Volume computeTightFittingBoundingBox( int num_tris, int *tri_indices );

	// Triangle getters.
	float getMinTriValue( int tri_index, SplitAxis axis );
	float getMaxTriValue( int tri_index, SplitAxis axis );
};

#endif
