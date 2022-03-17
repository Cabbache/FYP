#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>
#include <set>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "vec3.h"
#include "tiny_obj_loader.h"

using namespace std;

double total_marchtime = 0;
const double sres_to_gres = 2;

uint32_t total_cells = 0;
uint32_t total_hits = 0;

typedef struct hitInfo{
	vec3 ray;
	vec3 origin;

	bool hit;
	double t;

	double beta;
	double gamma;
} hitInfo;

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
	vector<Triangle> triangles;
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

//maps world coordinates to sdf coordinates and returns value of sdf
bool getValue(const SDF &sdf, const vec3 &world, SDFResult &res){
	if(
		world.x() < sdf.origin.x() || world.x() >= sdf.corner.x() - sdf.resolution/2 ||
		world.y() < sdf.origin.y() || world.y() >= sdf.corner.y() - sdf.resolution/2 ||
		world.z() < sdf.origin.z() || world.z() >= sdf.corner.z() - sdf.resolution/2
	){
		return res.inside = false;
	}
	res.inside = true;
	vec3 distance = world - sdf.origin;
	distance /= sdf.resolution;
	res.value = sdf.values[(int)round(distance[0])][(int)round(distance[1])][(int)round(distance[2])];
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

//loads triangles from obj file into triangles
void loadTriangles(string inputfile, vector<Triangle> &triangles){
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

void marchSDF(const SDF &sdf, hitInfo &info){
	//assumes info.ray is unit vector
	SDFResult sres;
	vec3 p;
	info.hit = false;
	for (p = info.origin;getValue(sdf, p, sres);p += info.ray*abs(sres.value)){
		if (sres.value <= 0){
			info.hit = true;
			info.t = (p-info.origin).length();
			return;
		}
	}
}

void boxHit(hitInfo &info, const Volume &box){
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
		return;
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
		return;
	}

	if (tzmin > tmin) 
		tmin = tzmin; 

	if (tzmax < tmax) 
		tmax = tzmax; 

	info.hit = true;
	info.t = tmin;
}

vec3 get_color(vec3 origin, vec3 ray, const vector<Obj> &world, int depth=0){
	//sometimes rays get stuck inside the 3d model
	if (depth > 40)
		return vec3(0,0,0);
	ray = unit_vector(ray);

	Triangle closestTri;
	hitInfo closest = {
		ray,
		origin,
		false,
		std::numeric_limits<float>::max(),
	};

	vector<Ohd> boxes_hit;
	
	for (auto it = world.begin(); it != world.end(); ++it){
		hitInfo boxintersect = {
			ray,
			origin,
			false,
			0,
		};
		boxHit(boxintersect, it->bounds);
		if (boxintersect.hit)
			boxes_hit.push_back(make_pair(&(*it), boxintersect.t));
	}

	//order boxes by their distance from ray origin
	sort(
		boxes_hit.begin(),
		boxes_hit.end(),
		[](const Ohd &o1, const Ohd &o2) -> bool{
			return o1.second > o2.second;
		}
	);

	for (Ohd ohd: boxes_hit){
		hitInfo march = {
			ray,
			origin + (ohd.second * ray), //start marching from ray box intersection
			false,
			0,
		};

		auto start = std::chrono::system_clock::now();
		marchSDF(ohd.first->sdf, march);
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
		total_marchtime += elapsed_seconds.count();

		if (!march.hit)
			continue;

		vec3 hitpoint = (march.origin + (march.ray * march.t)) / ohd.first->grid.resolution;

		set<vec3_int> points;
		points.insert(hitpoint); //to start checking from this cell
		for (int a = -1;a <= 1;a++)
		for (int b = -1;b <= 1;b++)
		for (int c = -1;c <= 1;c++){
			if (abs(a) + abs(b) + abs(c) != 1)
				continue;
			points.insert(
				vec3_int(
					hitpoint + vec3(
						a*ohd.first->sdf.resolution*60,
						b*ohd.first->sdf.resolution*60,
						c*ohd.first->sdf.resolution*60
					)
				)
			);
		}

		total_hits++;

		for (set<vec3_int>::iterator it = points.begin();it != points.end();++it){
			if (ohd.first->grid.map.count(*it) != 1)
				continue;
			total_cells++;
			for (Triangle tri : ohd.first->grid.map.at(*it)){
				hitInfo check = {
					ray,
					origin,
					false,
					0.0f
				};
				triHit(tri, check);
				if (!check.hit)
					continue;
				
				closest = check;
				closestTri = tri;

				goto end_iter;
			}
		}
		end_iter:;
	}

	//if no triangles hit, color the background
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

		return get_color(hitpoint, reflected, world, ++depth);
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
	//vec3 color(127 * ((unit_vector(hitpoint - (obj.bounds.min + obj.bounds.max) / 2)) + vec3(1,1,1)));
	//return vec3(color.z(), color.x(), color.y());
}

Volume getBoundingVolume(const vector<Triangle> &triangles){
	vec3 min(
		std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max(),
		std::numeric_limits<double>::max()
	);
	vec3 max(
		std::numeric_limits<double>::min(),
		std::numeric_limits<double>::min(),
		std::numeric_limits<double>::min()
	);
	for (Triangle tri : triangles){
		for (int i = 0;i < 3;i++)
		for (int j = 0;j < 3;j++){
			min.e[j] = std::min(min.e[j], tri.p[i].e[j]);
			max.e[j] = std::max(max.e[j], tri.p[i].e[j]);
		}
	}
	return Volume{min, max};
}

//This does not translate the gridmap
void translateObj(Obj &object, vec3 translation){
	//shfit triangles
	for (auto it = object.triangles.begin();it != object.triangles.end();++it){
		it->p[0] += translation;
		it->p[1] += translation;
		it->p[2] += translation;
	}
	object.bounds.min += translation;
	object.bounds.max += translation;
	
	object.sdf.origin += translation;
	object.sdf.corner += translation;
}

int main(int argc, char **argv){
	const unsigned int image_width = 1536;
	const unsigned int image_height = 1152;

	//const unsigned int image_width = 160;
	//const unsigned int image_height = 120;
	const unsigned int aliasing_iters = 2;
	const double angle = 1.0;
	const double aspect = (double)image_width / image_height;
	const string objFiles = "loads.txt";

	vector<Obj> world;
	cerr << "loading objs" << endl;

	ifstream obsFile(objFiles);
	if (!obsFile.is_open()){
		cerr << "Failed to open objs file" << endl;
		return 1;
	}

	vec3 translations[] = {
		vec3(0,0,0),
		vec3(1,0,0),
		vec3(0,3,0),
		vec3(0,0,1),
		vec3(-1,0,0),
		vec3(0,-1,0),
		vec3(0,0,-1)
	};

	string line;
	int tt = 0;
	while (getline(obsFile, line)){
		cerr << "Loading " << line << endl;
		Obj object;
		loadTriangles(line + ".obj", object.triangles);
		cerr << "loading sdf" << endl;
		loadSDF(object.sdf, line + ".sdf");

		object.bounds = getBoundingVolume(object.triangles);
		translateObj(object, translations[tt++]);

		cout << object.bounds.min << ", " << object.bounds.max << endl;

		cerr << "generating hashmap" << endl;
		object.grid.resolution = object.sdf.resolution * sres_to_gres;
		//walk around on barycentric coordinates of each triangle

		int last = 0;
		for (int i = 0;i < object.triangles.size();i++){
			int percentage = (100 * (float)i / object.triangles.size());
			if (last != percentage && percentage && percentage % 5 == 0){
				cerr << percentage << "%" << endl;
				last = percentage;
			}
			Triangle tri = object.triangles.at(i);
			set<vec3_int> points;
			double a_res = 0.3 * object.grid.resolution / tri.p[0].length();
			double b_res = 0.3 * object.grid.resolution / tri.p[1].length();
			for (double a = 0;a <= 1.0 + a_res;a+=a_res){
				for (double b = 0;b <= 1.0 - min(1.0,a) + b_res;b+=b_res){

					//clipping
					double a_c = min(1.0, a);
					double b_c = min(1.0-a_c, b);
					double c = 1.0 - a_c - b_c;

					vec3_int point(
							(
								a_c*tri.p[0] +
								b_c*tri.p[1] +
								c*tri.p[2]
							) / object.grid.resolution
						);
					points.insert(point);
				}
			}
			for (set<vec3_int>::iterator it = points.begin();it != points.end();++it)
				object.grid.map[*it].insert(tri);
		}

		cerr << "counting triangles" << endl;
		int avg = 0;
		int count = 0;
		for (GridMap::iterator iter = object.grid.map.begin(); iter != object.grid.map.end(); ++iter){
			avg += iter->second.size();
			count++;
		}
		cerr << (avg / (float)count) << " triangles / cell (avg)"<< endl;
		world.push_back(object);
	}

	Obj object = world.at(0);
	const double cam_distance = max(
		max(
			abs(object.bounds.min.x()),
			abs(object.bounds.min.z())
		),
		max(
			abs(object.bounds.max.x()),
			abs(object.bounds.max.z())
		)
	) * 2;

	double total_duration = 0;
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

		vec3 camera_origin = vec3(0,0.5*(object.bounds.max.y() + object.bounds.min.y()),-cam_distance);
		camera_origin = rotateY(camera_origin, angle_loop);

		vec3 frame_topleft = vec3(-frame_width/2, frame_height/2, eye_frame_distance);
		vec3 image[image_height*image_width];
		
		cerr << "Starting timer" << endl;
		auto start = std::chrono::system_clock::now();
		#pragma omp parallel for num_threads(8)
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
					average += get_color(camera_origin, ray, world);
				}
				average /= aliasing_iters;
				image[y*image_width + x] = average;
			}
		}
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
		cerr << "duration: " << elapsed_seconds.count() << endl;
		total_duration += elapsed_seconds.count();
		cerr << "march time: " << (100.0 * total_marchtime / elapsed_seconds.count()) << "%" << endl;
		cerr << total_cells / (double)total_hits << " blocks / hit " << endl;
		total_marchtime = 0.0;
		total_cells = 0;
		total_hits = 0;

		cerr << "Writing image to file" << endl;

		for (int i = 0;i < image_height*image_width;i++){
			vec3 color = image[i];
			ppm << (int)color[0] << " " << (int)color[1] << " " << (int)color[2] << endl;
		}
		ppm.close();
	}
	cerr << "total duration: " << total_duration << endl;
	return 0;
}
