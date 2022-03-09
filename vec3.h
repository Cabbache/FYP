#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <set>
#include <limits>
#include <unordered_set>

using namespace std;

class vec3 {
	public:
		vec3() : e{0,0,0} {}
		vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}

		double x() const { return e[0]; }
		double y() const { return e[1]; }
		double z() const { return e[2]; }

		void setX(double x){ e[0] = x; }
		void setY(double y){ e[1] = y; }
		void setZ(double z){ e[2] = z; }

		vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
		double operator[](int i) const { return e[i]; }
		double& operator[](int i) { return e[i]; }

		bool operator==(const vec3 &v) const{
			return e[0] == v.e[0] && e[1] == v.e[1] && e[2] == v.e[2];
		}

		vec3& operator+=(const vec3 &v) {
			e[0] += v.e[0];
			e[1] += v.e[1];
			e[2] += v.e[2];
			return *this;
		}

		vec3& operator*=(const double t) {
			e[0] *= t;
			e[1] *= t;
			e[2] *= t;
			return *this;
		}

		vec3& operator/=(const double t) {
			return *this *= 1/t;
		}

		vec3& max(const vec3 &v){
			for (int i = 0;i < 3;i++)
				e[i] = std::max(e[i], v[i]);
			return *this;
		}

		vec3& min(const vec3 &v){
			for (int i = 0;i < 3;i++)
				e[i] = std::min(e[i], v[i]);
			return *this;
		}

		double length() const {
			return sqrt(length_squared());
		}

		double length_squared() const {
			return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
		}

	public:
		double e[3];
};

class vec3_int{
	public:
		int x,y,z;

		vec3_int(int x, int y, int z): x(x), y(y), z(z){};
		vec3_int(vec3 point){
			this->x = (int)round(point[0]);
			this->y = (int)round(point[1]);
			this->z = (int)round(point[2]);
		}

		bool operator< (const vec3_int &point) const{
			return this->x < point.x;
		}
};

struct hashFuncVec{
	size_t operator()(const vec3_int &k) const{
		size_t h1 = std::hash<int>()(k.x);
		size_t h2 = std::hash<int>()(k.y);
		size_t h3 = std::hash<int>()(k.z);
		return (h1 ^ (h2 << 1)) ^ h3;
	}
};

struct equalsFunc{
	bool operator()( const vec3_int& lhs, const vec3_int& rhs ) const{
		return (lhs.x == rhs.x) && (lhs.y == rhs.y) && (lhs.z == rhs.z);
	}
};

inline std::ostream& operator<<(std::ostream &out, const vec3_int &v) {
	return out << v.x << ' ' << v.y << ' ' << v.z;
}

// Type aliases for vec3
using point3 = vec3;   // 3D point
using color = vec3;	// RGB color

// vec3 Utility Functions
inline std::ostream& operator<<(std::ostream &out, const vec3 &v) {
	return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline vec3 operator+(const vec3 &u, const vec3 &v) {
	return vec3(u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vec3 operator-(const vec3 &u, const vec3 &v) {
	return vec3(u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vec3 operator*(const vec3 &u, const vec3 &v) {
	return vec3(u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vec3 operator*(double t, const vec3 &v) {
	return vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline vec3 operator*(const vec3 &v, double t) {
	return t * v;
}

inline vec3 operator/(vec3 v, double t) {
	return (1/t) * v;
}

inline double dot(const vec3 &u, const vec3 &v) {
	return u.e[0] * v.e[0]
		 + u.e[1] * v.e[1]
		 + u.e[2] * v.e[2];
}

inline vec3 cross(const vec3 &u, const vec3 &v) {
	return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
				u.e[2] * v.e[0] - u.e[0] * v.e[2],
				u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vec3 unit_vector(vec3 v) {
	return v / v.length();
}

inline vec3 rotateY(vec3 v, double degrees){
	double radians = 3.1415926535 * degrees / 180;
	return vec3(
		v.e[0] * cos(radians) + v.e[2] * sin(radians),
		v.e[1],
		-v.e[0] * sin(radians) + v.e[2] * cos(radians)
	);
}

typedef struct Volume{
	vec3 min;
	vec3 max;
} Volume;

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

typedef struct hitInfo{
	vec3 ray;
	vec3 origin;

	bool hit;
	double t;

	double beta;
	double gamma;

	Triangle tri;
} hitInfo;

typedef unordered_map<vec3_int, unordered_set<Triangle, Triangle::HashFunction>, hashFuncVec, equalsFunc> GridMap;

typedef struct SDF{
	vec3 origin;
	vec3 corner;
	vec3 dimensions;
	double resolution;
	double ***values;
} SDF;

typedef struct Obj{
	vec3 *verts;
	vec3 *faces;
	int numVerts;
	int numFaces;
	Volume bounds;
	SDF sdf;
} Obj;

//used for the function that returns sdf value given world coordinate
//if world coordinate is outside sdf domain, inside = false
typedef struct SDFResult{
	bool inside;
	double value;
} SDFResult;

#endif
