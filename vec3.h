#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

using std::sqrt;

class vec3 {
	public:
		vec3() : e{0,0,0} {}
		vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}

		double x() const { return e[0]; }
		double y() const { return e[1]; }
		double z() const { return e[2]; }

		vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
		double operator[](int i) const { return e[i]; }
		double& operator[](int i) { return e[i]; }

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
		vec3_int(int x, int y, int z): x(x), y(y), z(z){};
		vec3_int(vec3 point){
			this->x = (int)round(point[0]);
			this->y = (int)round(point[1]);
			this->z = (int)round(point[2]);
		}
		int x,y,z;
};

struct hashFunc{
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

#endif
