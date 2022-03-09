#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc

#include "vec3.h"
#include "mesh.h"
#include "intersections.h"
#include "KDTreeCPU.h"

using namespace std;

double total_marchtime = 0;
const double sres_to_gres = 2.3;

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
void loadTriangles(string inputfile, Obj &object){
	mesh msh(inputfile);
	object.faces = msh.tris;
	object.verts = msh.verts;
	object.numFaces = msh.numTris;
	object.numVerts = msh.numVerts;
	object.bounds = msh.bb;
}

double drand(){
	return (double)rand()/(double)RAND_MAX;
}

void marchSDF(const SDF &sdf, hitInfo &info){
	//assumes info.ray is unit vector
	SDFResult sres;
	vec3 p;
	info.hit = false;
	for (p = info.origin;getValue(sdf, p, sres);p += info.ray*abs(sres.value)){
		if (sres.value <= sdf.resolution/2.0){ // <= sdf.resolution or <= 0 ?
			info.hit = true;
			info.t = (p-info.origin).length();
			return;
		}
	}
}

vec3 get_color(vec3 origin, vec3 ray, KDTreeCPU *&kdtree, int depth=0){
	//sometimes rays get stuck inside the 3d model
	if (depth > 40)
		return vec3(0,0,0);
		
	ray = unit_vector(ray);

	hitInfo hitinfo{
		ray,
		origin,
		false
	};
	vec3 hitpoint;
	bool hit = kdtree->intersect(hitinfo, hitpoint);

	//if no triangles hit, color the background
	if (!hit){
		vec3 pp = 127*(ray+vec3(1,1,1));
		return vec3((int)pp[1], (int)pp[2], (int)pp[0]);
	}

	//if hit border of non mirror triangle, make border black
	if (hitinfo.beta < 0.01 || hitinfo.gamma < 0.01 || (1-hitinfo.beta-hitinfo.gamma) < 0.01)
		return vec3(0,0,0);
	
	//color the triangle
	vec3 color = unit_vector(
		vec3(
			(hitpoint - hitinfo.tri.p[0]).length(),
			(hitpoint - hitinfo.tri.p[1]).length(),
			(hitpoint - hitinfo.tri.p[2]).length()
		)
	)*255;
	return color;
}

int main(int argc, char **argv){
	const unsigned int image_width = 1280;
	const unsigned int image_height = 960;

	const unsigned int aliasing_iters = 2;
	const double angle = 1.0;
	//const double cam_distance = 0.5;
	const double aspect = (double)image_width / image_height;

	cerr << "loading obj" << endl;
	Obj object;
	loadTriangles("object.obj", object);
	cout << object.bounds.min << ", " << object.bounds.max << endl;
	const double cam_distance = max(
		max(
			object.bounds.min.x(),
			object.bounds.min.z()
		),
		max(
			object.bounds.max.x(),
			object.bounds.max.z()
		)
	) * 4;

	cerr << "loading sdf" << endl;
	loadSDF(object.sdf, "object.sdf");

	cerr << "Generating KDTree" << endl;
	KDTreeCPU *kdtree = new KDTreeCPU(
		object.numFaces,
		object.faces,
		object.numVerts,
		object.verts
	);
	cerr << "Starting render" << endl;

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
		#pragma omp parallel for num_threads(4)
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
					average += get_color(camera_origin, ray, kdtree);
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
		total_marchtime = 0.0;

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
