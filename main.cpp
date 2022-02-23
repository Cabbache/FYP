#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>
#include <set>
#include <limits>

#include "vec3.h"
#include "tiny_obj_loader.h"

using namespace std;

typedef struct Triangle{
	vec3 p[3];
	bool reflective;
} Triangle;

typedef struct Volume{
	vec3 min;
	vec3 max;
} Volume;

//triangles only
typedef struct hitInfo{
	vec3 ray;
	vec3 origin;

	bool hit;
	double t;

	double beta;
	double gamma;
} hitInfo;

typedef struct SDF{
	vec3 origin;
	vec3 corner;
	vec3 dimensions;
	double resolution;
	double ***values;
} SDF;

//used for the function that returns sdf value given world coordinate
//if world coordinate is outside sdf domain, inside = false
typedef struct SDFResult{
	bool inside;
	double value;
} SDFResult;

static int gcount = 0;

//test with 1000 triangles 200x150:
//38 seconds - no oct
//25 seconds - oct leafel = 25
//27 seconds - oct leafel = 5
//23 seconds - oct leafel = 50
class OCT{
	public:
		OCT(const vector<Triangle> &triangles, int leafElements=50){
			box = {
				vec3(0,0,0),
				vec3(0,0,0)
			};
			for (int i = 0;i < triangles.size();i++){
				Triangle tri = triangles.at(i);
				for (int j = 0;j < 3;j++){
					box.min.min(tri.p[j]);
					box.max.max(tri.p[j]);
				}
			}

			if (isLeaf = triangles.size() <= leafElements){
				leafTri = triangles;
				return;	
			}

			sep = vec3(
				(box.max[0] + box.min[0]) / 2,
				(box.max[1] + box.min[1]) / 2,
				(box.max[2] + box.min[2]) / 2
			);
			
			vector<Triangle> octs[8];
			for (int i = 0;i < triangles.size();i++){
				Triangle tri = triangles.at(i);
				std::set<int> adds;
				for (int j = 0;j < 1;j++){ //TODO investigate if j < 1 or j < 3
					adds.insert(
						(tri.p[j].x() > sep.x()) << 2 |
						(tri.p[j].y() > sep.y()) << 1 |
						(tri.p[j].z() > sep.z()) << 0
					);
				}
				for (auto it = adds.begin();it != adds.end();it++)
					octs[*it].push_back(tri);
			}
			
			for (int i = 0;i < 8;i++){
				if (octs[i].size() == triangles.size() && (isLeaf = true)){
					leafTri = triangles;
					return;
				}
			}

			for (int i = 0;i < 8;i++)
				children[i] = new OCT(octs[i]);
		}

		void getGroup(vec3 &ray, vec3 &origin, vector<Triangle> &tris) {
			if (isLeaf){
				tris = leafTri;
				return;
			}

			for (int i = 0;i < 8;i++){
				hitInfo info;
				children[i]->boxHit(ray, origin, info);
				if (!info.hit)
					continue;
				vector<Triangle> ctris;
				children[i]->getGroup(ray, origin, ctris);
				tris.insert(tris.end(), ctris.begin(), ctris.end());
			}
		}

	private:
		Volume getBox(){
			return box;
		}

		vec3 getSep(){
			return sep;
		}

		//https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
		void boxHit(const vec3 &origin, const vec3 &ray, hitInfo &info){
			float tmin = (box.min.x() - origin.x()) / ray.x(); 
			float tmax = (box.max.x() - origin.x()) / ray.x(); 
 
			if (tmin > tmax) {
				float tmp = tmin;
				tmin = tmax;
				tmax = tmp;
			}
	 
			float tymin = (box.min.y() - origin.y()) / ray.y(); 
			float tymax = (box.max.y() - origin.y()) / ray.y(); 
	 
			if (tymin > tymax) {
				float tmp = tymin;
				tymin = tymax;
				tymax = tmp;
			}
	 
			if ((tmin > tymax) || (tymin > tmax)) {
				info.hit = false;
				return;
			}
	 
			if (tymin > tmin) 
					tmin = tymin; 
	 
			if (tymax < tmax) 
					tmax = tymax; 
	 
			float tzmin = (box.min.z() - origin.z()) / ray.z(); 
			float tzmax = (box.max.z() - origin.z()) / ray.z();
	 
			if (tzmin > tzmax){
				float tmp = tzmin;
				tzmin = tzmax;
				tzmax = tmp;
			}
	 
			if ((tmin > tzmax) || (tzmin > tmax)) {
				info.hit = false;
				return;
			}
	 
			if (tzmin > tmin) 
					tmin = tzmin; 
	 
			if (tzmax < tmax) 
					tmax = tzmax; 
	 
	 		info.hit = true;
			info.t = tmin;
		}

		OCT *children[8];
		vector<Triangle> leafTri;
		Volume box;
		vec3 sep;
		bool isLeaf;
};

//maps world coordinates to sdf coordinates and returns value of sdf
bool getValue(const SDF &sdf, const vec3 &world, SDFResult &res){
	if(
		world.x() < sdf.origin.x() || world.x() > sdf.corner.x() ||
		world.y() < sdf.origin.y() || world.y() > sdf.corner.y() ||
		world.z() < sdf.origin.z() || world.z() > sdf.corner.z()
	){
		return res.inside = false;
	}
	res.inside = true;
	vec3 distance = world - sdf.origin;
	distance /= sdf.resolution;
	res.value = sdf.values[(int)distance[0]][(int)distance[1]][(int)distance[2]];
	return true;
}

//loads sdf from sdf file
void loadSDF(SDF &sdf, string filename){
	ifstream sfile(filename);
	if (!sfile.is_open()){
		cerr << "Failed to load SDF file" << endl;
		return;
	}
	string line;
	string value;

	getline(sfile, line);
	stringstream ss(line);
	getline(ss, value, ' ');
	sdf.dimensions[0] = stoi(value);
	getline(ss, value, ' ');
	sdf.dimensions[1] = stoi(value);
	getline(ss, value, ' ');
	sdf.dimensions[2] = stoi(value);

	sdf.values = new double**[(int)sdf.dimensions[0]];
	for (int a = 0;a < sdf.dimensions[0];a++){
		sdf.values[a] = new double*[(int)sdf.dimensions[1]];
		for (int b = 0;b < sdf.dimensions[1];b++)
			sdf.values[a][b] = new double[(int)sdf.dimensions[2]];
	}

	getline(sfile, line);
	ss.str(line);
	ss.clear();
	getline(ss, value, ' ');
	sdf.origin[0] = stod(value);
	getline(ss, value, ' ');
	sdf.origin[1] = stod(value);
	getline(ss, value, ' ');
	sdf.origin[2] = stod(value);

	getline(sfile, line);
	sdf.resolution = stod(line);

	sdf.corner = sdf.origin + sdf.resolution*sdf.dimensions;

	//line number technically not 0
	for (unsigned long line_number = 0;getline(sfile, line);++line_number){
		sdf.values[
			line_number % (int)sdf.dimensions[0]
		][
			(int)(line_number/sdf.dimensions[0]) % (int)sdf.dimensions[1]
		][
			(int)(line_number/(sdf.dimensions[0]*sdf.dimensions[1])) % (int)sdf.dimensions[2]
		] = stod(line);
	}
}

//loads triangles from obj file
vector<Triangle> loadTriangles(string inputfile){
	vector<Triangle> triangles;

	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;

	std::string warn;
	std::string err;

	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str());

	if (!warn.empty()) {
		std::cerr << warn << std::endl;
	}

	if (!err.empty()) {
		std::cerr << err << std::endl;
	}

	if (!ret) {
		exit(1);
	}

	// Loop over shapes
	for (size_t s = 0; s < shapes.size(); s++) {
		// Loop over faces(polygon)
		size_t index_offset = 0;
		for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
			size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);

			// Loop over vertices in the face.
			//cout << fv << endl;
			Triangle tr;
			tr.reflective = false;
			for (size_t v = 0; v < fv; v++) {
				// access to vertex
				tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];

				tinyobj::real_t vx = attrib.vertices[3*size_t(idx.vertex_index)+0];
				tinyobj::real_t vy = attrib.vertices[3*size_t(idx.vertex_index)+1];
				tinyobj::real_t vz = attrib.vertices[3*size_t(idx.vertex_index)+2];

				tr.p[v] = vec3(vx, vy, vz);

				// Check if `normal_index` is zero or positive. negative = no normal data
				if (idx.normal_index >= 0) {
					tinyobj::real_t nx = attrib.normals[3*size_t(idx.normal_index)+0];
					tinyobj::real_t ny = attrib.normals[3*size_t(idx.normal_index)+1];
					tinyobj::real_t nz = attrib.normals[3*size_t(idx.normal_index)+2];
				}

				// Check if `texcoord_index` is zero or positive. negative = no texcoord data
				if (idx.texcoord_index >= 0) {
					tinyobj::real_t tx = attrib.texcoords[2*size_t(idx.texcoord_index)+0];
					tinyobj::real_t ty = attrib.texcoords[2*size_t(idx.texcoord_index)+1];
				}
				// Optional: vertex colors
				// tinyobj::real_t red   = attrib.colors[3*size_t(idx.vertex_index)+0];
				// tinyobj::real_t green = attrib.colors[3*size_t(idx.vertex_index)+1];
				// tinyobj::real_t blue  = attrib.colors[3*size_t(idx.vertex_index)+2];
			}

			triangles.push_back(tr);
			index_offset += fv;

			// per-face material
			shapes[s].mesh.material_ids[f];
		}
	}

	return triangles;
}

double drand(){
	return (double)rand()/(double)RAND_MAX;
}

//pages 36 - 38, triangle - ray intersection test
void triHit(const Triangle &tri, hitInfo &info){
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
		return;
	info.t = t;
	info.beta = beta;
	info.gamma = gamma;
}

//ray marching, returns vector of sdf points on ray
vector<double> rayValues(const SDF &sdf, const vec3 &origin, vec3 ray){
	ray = unit_vector(ray) * sdf.resolution;
	SDFResult sres;
	vector<double> values;
	vec3 p;
	for (p = origin;getValue(sdf, p, sres);p += ray)
		values.push_back(sres.value);
	int exp = (p-origin).length()/sdf.resolution;
	int actual = values.size();
	return values;
}

vec3 get_color(vec3 origin, vec3 ray, const vector<Triangle> &triangles, int depth=0){

	//sometimes rays get stuck inside the 3d model
	if (depth > 40)
		return vec3(0,0,0);
		
	Triangle closestTri;
	hitInfo closest = {
		ray,
		origin,
		false,
		std::numeric_limits<float>::max()
	};

	for (Triangle tri : triangles){
		hitInfo check = {
			ray,
			origin,
			false,
			0.0f
		};

		triHit(tri, check);
		if (!check.hit)
			continue;

		if (closest.t < check.t)
			continue;
		closest = check;
 		closestTri = tri;
	}

	//if no triangles hit color the background
	if (!closest.hit){
		vec3 pp = 127*(unit_vector(ray)+vec3(1,1,1));
		return vec3((int)pp[1], (int)pp[2], (int)pp[0]);
	}

	vec3 hitpoint = origin + (ray * closest.t);
	
	//if triangle is mirror
	if (closestTri.reflective){
		vec3 normal = unit_vector(
			cross(
				closestTri.p[0] - closestTri.p[1],
				closestTri.p[0] - closestTri.p[2]
			)
		);
		vec3 reflected = unit_vector(ray - 2*dot(ray, normal)*normal);

		//fuzzy mirror
		//vec3 deflect = unit_vector(vec3(drand()-0.5, drand()-0.5, drand()-0.5));
		//reflected += deflect / 20;

		return get_color(hitpoint, reflected, triangles, ++depth);
	}

	//if hit border of non mirror triangle, make border black
	if (closest.beta < 0.01 || closest.gamma < 0.01 || (1-closest.beta-closest.gamma) < 0.01)
		return vec3(0,0,0);
	
	//color the non mirror triangle
	vec3 color = unit_vector(
		vec3(
			(hitpoint - closestTri.p[0]).length(),
			(hitpoint - closestTri.p[1]).length(),
			(hitpoint - closestTri.p[2]).length()
		)
	)*255;
	return color;
}

int main(int argc, char **argv){
	const unsigned int image_width = 200;
	const unsigned int image_height = 150;

	const unsigned int aliasing_iters = 2;
	const double angle = 3.0;
	//const double cam_distance = 0.35;
	const double cam_distance = 0.3;
	const double aspect = (double)image_width / image_height;

	cerr << "Loading obj" << endl;
	vector<Triangle> triangles = loadTriangles("bunny_1k.obj");
	cerr << "Loading complete" << endl;

	auto global_start = std::chrono::system_clock::now();
	for (double angle_loop = 0;angle_loop < 360;angle_loop += angle){
		ofstream ppm("img_"+to_string(int(angle_loop))+".ppm");

		//ppm image headers
		ppm << "P3" << endl
		<< image_width << " " << image_height << endl
		<< "255" << endl;

		//frame buffer
		double frame_width = 1;
		double frame_height = frame_width / aspect;
		double eye_frame_distance = 1; //or focal length?

		vec3 camera_origin = vec3(0,0.2,-cam_distance);
		camera_origin = rotateY(camera_origin, angle_loop);

		vec3 frame_topleft = vec3(-frame_width/2, frame_height/2, eye_frame_distance);
		vec3 image[image_height*image_width];
		
		cerr << "Starting timer" << endl;
		auto start = std::chrono::system_clock::now();
		#pragma omp parallel for
		for (int y = 0;y < image_height;y++){
			for (int x = 0;x < image_width;x++){
				vec3 average(0,0,0);

				//anti aliasing loop
				for (int a = 0;a < aliasing_iters;a++){
					vec3 ray = frame_topleft;
					ray += vec3(
						frame_width * ((x+drand()) / (double)image_width),
						-frame_height * ((y+drand()) / (double)image_height),
						0
					);
					ray = rotateY(ray, angle_loop);
					average += get_color(camera_origin, ray, triangles);
				}
				average /= aliasing_iters;
				image[y*image_width + x] = average;
			}
		}
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
		cerr << "duration: " << elapsed_seconds.count() << endl;

		cerr << "Writing image to file" << endl;

		for (int i = 0;i < image_height*image_width;i++){
			vec3 color = image[i];
			ppm << (int)color[0] << " " << (int)color[1] << " " << (int)color[2] << endl;
		}
		ppm.close();
	}
	auto global_end = std::chrono::system_clock::now();
	std::chrono::duration<double> global_elapsed = global_end - global_start;
	cerr << "total duration: " << global_elapsed.count() << endl;
	return 0;
}
