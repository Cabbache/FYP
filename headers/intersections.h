#ifndef INTERSECTIONS_H
#define INTERSECTIONS_H

#include "vec3.h"
#include "structures.h"

class Intersections{
public:
	static bool aabbIntersect(Volume bbox, vec3 ray_o, vec3 ray_dir, float &t_near, float &t_far )
	{
		vec3 dirfrac( 1.0f / ray_dir.x(), 1.0f / ray_dir.y(), 1.0f / ray_dir.z() );

		float t1 = ( bbox.min.x() - ray_o.x() ) * dirfrac.x();
		float t2 = ( bbox.max.x() - ray_o.x() ) * dirfrac.x();
		float t3 = ( bbox.min.y() - ray_o.y() ) * dirfrac.y();
		float t4 = ( bbox.max.y() - ray_o.y() ) * dirfrac.y();
		float t5 = ( bbox.min.z() - ray_o.z() ) * dirfrac.z();
		float t6 = ( bbox.max.z() - ray_o.z() ) * dirfrac.z();

		float tmin = std::max( std::max( std::min( t1, t2 ), std::min( t3, t4 ) ), std::min( t5, t6 ) );
		float tmax = std::min( std::min( std::max( t1, t2 ), std::max( t3, t4 ) ), std::max( t5, t6 ) );

		// If tmax < 0, ray intersects AABB, but entire AABB is behind ray, so reject.
		if ( tmax < 0.0f ) {
			return false;
		}

		// If tmin > tmax, ray does not intersect AABB.
		if ( tmin > tmax ) {
			return false;
		}

		t_near = tmin;
		t_far = tmax;
		return true;
	}

	static bool boxHit(hitInfo &info, const Volume &box){
		float tmin = (box.min.x() - info.origin.x()) / info.ray.x(); 
		float tmax = (box.max.x() - info.origin.x()) / info.ray.x(); 

		if (tmin > tmax) {
			float tmp = tmin;
			tmin = tmax;
			tmax = tmp;
		}

		float tymin = (box.min.y() - info.origin.y()) / info.ray.y(); 
		float tymax = (box.max.y() - info.origin.y()) / info.ray.y(); 

		if (tymin > tymax) {
			float tmp = tymin;
			tymin = tymax;
			tymax = tmp;
		}

		if ((tmin > tymax) || (tymin > tmax)) {
			info.hit = false;
			return false;
		}

		if (tymin > tmin) 
			tmin = tymin; 

		if (tymax < tmax) 
			tmax = tymax; 

		float tzmin = (box.min.z() - info.origin.z()) / info.ray.z(); 
		float tzmax = (box.max.z() - info.origin.z()) / info.ray.z();

		if (tzmin > tzmax){
			float tmp = tzmin;
			tzmin = tzmax;
			tzmax = tmp;
		}

		if ((tmin > tzmax) || (tzmin > tmax)) {
			info.hit = false;
			return false;
		}

		if (tzmin > tmin) 
			tmin = tzmin; 

		if (tzmax < tmax) 
			tmax = tzmax; 

		info.hit = true;
		info.t = tmin;
		return true;
	}

	//pages 36 - 38, triangle - ray intersection test
	static bool triHit(const Triangle &tri, hitInfo &info){
		double a = tri.p[0].x() - tri.p[1].x();
		double b = tri.p[0].y() - tri.p[1].y();
		double c = tri.p[0].z() - tri.p[1].z();

		double d = tri.p[0].x() - tri.p[2].x();
		double e = tri.p[0].y() - tri.p[2].y();
		double f = tri.p[0].z() - tri.p[2].z();

		double g = info.ray.x();
		double h = info.ray.y();
		double i = info.ray.z();

		vec3 diff = tri.p[0] - info.origin;
		double j = diff.x();
		double k = diff.y();
		double l = diff.z();

		double M = a*(e*i - h*f) + b*(g*f - d*i) + c*(d*h - e*g);

		double beta = j*(e*i - h*f) + k*(g*f - d*i) + l*(d*h - e*g);
		double gamma = i*(a*k - j*b) + h*(j*c - a*l) + g*(b*l - k*c);
		double t = f*(a*k - j*b) + e*(j*c - a*l) + d*(b*l - k*c);

		beta /= M;
		gamma /= M;
		t /= -M;

		if (!(info.hit = beta + gamma < 1 && beta > 0 && gamma > 0 && t > 0.0001))
			return false;

		info.t = t;
		info.beta = beta;
		info.gamma = gamma;
		info.tri = tri;
		return true;
	}
};

#endif
