#include "mesh.h"
#include <iterator>

using namespace std;

mesh::mesh(){
	numTris = numVerts = 0;
}

std::vector<std::string> mesh::tokenizeString(std::string str){
    std::stringstream strstr(str);
    std::istream_iterator<std::string> it(strstr);
    std::istream_iterator<std::string> end;
    std::vector<std::string> results(it, end);
    return results;
}

std::istream& mesh::safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
            sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
            is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}

mesh::mesh(string fileName){
	ifstream ifile;
	string line;
	ifile.open(fileName.c_str());

	vector<vec3> verts, faces;

	while (safeGetline(ifile, line)) {
		vector<string> tokens = tokenizeString(line);

		if (tokens.size()>0 && strcmp(tokens[0].c_str(),"v")==0){
			verts.push_back(vec3(atof(tokens[1].c_str()),atof(tokens[2].c_str()),atof(tokens[3].c_str())));
		}
		else if (tokens.size()>0 && strcmp(tokens[0].c_str(),"f")==0){
			char* findex1 = strtok (const_cast<char*>(tokens[1].c_str()),"/");
			char* findex2 = strtok (const_cast<char*>(tokens[2].c_str()),"/");
			char* findex3 = strtok (const_cast<char*>(tokens[3].c_str()),"/");
			faces.push_back(vec3(atof(findex1)-1,atof(findex2)-1,atof(findex3)-1));
		}
	}

	vec3* vertData = new vec3[verts.size()];
	vec3* faceData = new vec3[faces.size()];

	vec3 max = vec3(-10000,-10000,-10000);
	vec3 min = vec3( 10000, 10000, 10000);

	for (int i=0; i<verts.size(); i++){
		vertData[i] = verts[i];
		if (verts[i].x()<min.x()){
			min.setX(verts[i].x());
		}
		if (verts[i].y()<min.y()){
			min.setY(verts[i].y());
		}
		if (verts[i].z()<min.z()){
			min.setZ(verts[i].z());
		}
		if (verts[i].x()>max.x()){
			max.setX(verts[i].x());
		}
		if (verts[i].y()>max.y()){
			max.setY(verts[i].y());
		}
		if (verts[i].z()>max.z()){
			max.setZ(verts[i].z());
		}
	}
	for (int i=0; i<faces.size(); i++){
		faceData[i] = faces[i];
	}

	this->tris = faceData;
	this->verts = vertData;
	this->numTris = faces.size();
	this->numVerts = verts.size();

	bb.min = min;
	bb.max = max;
}
