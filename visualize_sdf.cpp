#include <iostream>
#include <fstream>
#include <sstream>
#include "headers/vec3.h"

using namespace std;

typedef struct SDF{
	vec3 origin;
	vec3 corner;
	vec3 dimensions;
	float resolution;
	float ***values;
} SDF;

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

	sdf.values = new float**[(int)sdf.dimensions[0]];
	for (int a = 0;a < sdf.dimensions[0];a++){
		sdf.values[a] = new float*[(int)sdf.dimensions[1]];
		for (int b = 0;b < sdf.dimensions[1];b++)
			sdf.values[a][b] = new float[(int)sdf.dimensions[2]];
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

int main(int argc, char **argv){
	if (argc != 3){
		cerr << "Usage: " << argv[0] << " [file.sdf] [num_slices (0 means all)]" << endl;
		return 0;
	}

	int num_slices = atoi(argv[2]);

	cerr << "Loading sdf" << endl;
	SDF sdf;
	loadSDF(sdf, argv[1]);

	if (num_slices <= 0)
		num_slices = (int)sdf.dimensions.x();

	cerr << "Writing " << num_slices << " slices" << endl;

	for (int sliceX = 0;sliceX < (int)sdf.dimensions.x();sliceX += (int)sdf.dimensions.x() / num_slices){
		cerr << sliceX << endl;
		ofstream ppm("sdfslice_"+to_string(sliceX)+".ppm");
		ppm << "P3" << endl
		<< sdf.dimensions.y() << " " << sdf.dimensions.z() << endl
		<< "255" << endl;
		
		int slice[(int)sdf.dimensions.y()][(int)sdf.dimensions.z()];

		//figure min and max
		float minvalue_in, maxvalue_in, minvalue_out, maxvalue_out;
		minvalue_in = maxvalue_in = minvalue_out = maxvalue_out = sdf.values[sliceX][0][0];
		for (int y = 0;y < sdf.dimensions.y();y++){
			for (int z = 0;z < sdf.dimensions.z();z++){
				float value = sdf.values[sliceX][y][z];
				if (value > 0){
					if (value > maxvalue_out)
						maxvalue_out = value;
					if (value < minvalue_out)
						minvalue_out = value;
				} else {
					if (value > maxvalue_in)
						maxvalue_in = value;
					if (value < minvalue_in)
						minvalue_in = value;
				}
			}
		}

		//loop again but adjust from 0 to 255
		for (int z = 0;z < sdf.dimensions.z();z++){
			for (int y = 0;y < sdf.dimensions.y();y++){
				float value = sdf.values[sliceX][y][z];
				int adjusted;
				if (value >= 0)
					adjusted = (int)((value - minvalue_out) * (255.0 / (maxvalue_out - minvalue_out)));
				else 
					adjusted = (int)((value - minvalue_in) * (255.0 / (maxvalue_in - minvalue_in)));
				if (value >= 0)
					ppm << adjusted << " " << adjusted << " " << adjusted << endl;
				else
					ppm << "0 " << (int)(0.5*adjusted) << " " << adjusted << endl;
			}
		}

		ppm.close();
	}
	return 0;
}
