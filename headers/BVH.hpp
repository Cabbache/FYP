#ifndef BVHCLASS
#define BVHCLASS

#include <vector>
#include "vec3.h"

using namespace std;

class BVH;

class BVH{
	public:
		BVH(const vector<Obj> &world){
			this->bounds.min = vec3(
					std::numeric_limits<float>::max(),
					std::numeric_limits<float>::max(),
					std::numeric_limits<float>::max()
				);
			this->bounds.max = vec3(
					std::numeric_limits<float>::min(),
					std::numeric_limits<float>::min(),
					std::numeric_limits<float>::min()
				);
			this->isLeaf = false;
			build(world);
		}

		~BVH(){
			delete half1;
			delete half2;
		}

		static void boxHit(boxHitInfo &info, const Volume &box){
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

			if (tmin < 0 && tmax < 0){
				info.hit = false;
				return;
			}

			info.hit = true;
			info.tmax = tmax;
			info.tmin = tmin > 0 ? tmin:0;
		}

		bool traverse(boxHitInfo &bhi, Obj &obj) const{
			cerr << "traversing" << endl;
			if (this->isLeaf){
				boxHit(bhi, this->leafObject.bounds);
				if (!bhi.hit){
					return false;
				}
				obj = this->leafObject;
				return true;
			}
			boxHitInfo box1 = bhi;
			boxHitInfo box2 = bhi;
			boxHit(box1, this->half1->bounds);
			boxHit(box2, this->half2->bounds);
			if (!box1.hit && !box2.hit){
				return false;
			}
			else if (box1.hit != box2.hit){
				BVH *hitBVH = box1.hit ? this->half1:this->half2;
				return hitBVH->traverse(bhi, obj);
			} else if (box1.hit && box2.hit){
				BVH *hitBVH = box1.tmin  < box2.tmin ? this->half1:this->half2;
				return hitBVH->traverse(bhi, obj);
			}
			return false;
		}

	const Volume& getBounds(){
		return this->bounds;
	}

	private:
		BVH *half1;
		BVH *half2;
		Volume bounds;
		Obj leafObject;
		bool isLeaf;

		void build(const vector<Obj> &world){
			cerr << world.size() << endl;
			if (world.size() <= 1){
				if (world.size() == 0)
					cerr << "Impossible" << endl;
				this->isLeaf = true;
				this->leafObject = world.at(0);
				cerr << this->leafObject.bounds.min << endl;
				cerr << this->leafObject.bounds.max << endl;
				cerr << "-" << endl;
				return;
			}

			for (auto &object : world){
				this->bounds.min.min(object.bounds.min);
				this->bounds.max.max(object.bounds.max);
			}
			
			vec3 boxSize = this->bounds.max - this->bounds.min;
			pair<vector<Obj>, vector<Obj>> splits[3];
			for (int i = 0;i < 3;i++){
				for (auto &object : world){
					if ((object.bounds.max - this->bounds.min)[i] > boxSize[i]/2.){
						splits[i].first.push_back(object);
					} if ((object.bounds.min - this->bounds.min)[i] < boxSize[i]/2.){
						splits[i].second.push_back(object);
					}	
				}
			}
			
			//find partition that is closest to a 50 50 split
			int bestIndex = 0;
			float bestValue = abs(splits[0].first.size() - world.size()/2.);
			for (int i = 1;i < 3;i++){
				float value = abs(splits[i].first.size() - world.size()/2.);
				if (value >= bestValue) continue;
				bestValue = value;
				bestIndex = i;
			}

			half1 = new BVH(splits[bestIndex].first);
			half2 = new BVH(splits[bestIndex].second);
		}
};

#endif
