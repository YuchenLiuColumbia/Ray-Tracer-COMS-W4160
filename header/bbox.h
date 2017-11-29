#ifndef BBOX_H_
#define BBOX_H_

#include "basicmath.h"
#include "ray.h"

class bbox
{
public:

	void buildbox(mypoint p1, mypoint p2, mypoint p3);

	void buildbox(mypoint pos, double r);

	bool hitbox(ray &thisray, double &int_near_r, double &int_far_r);

	void combinebox(bbox box1, bbox box2);

	bool insidebox(const ray thisray);

	//heart of box
	double h[3];


private:

	// bbox_min[0] - xmin, bbox_min[1] - ymin, bbox_min[2] - zmin
	// bbox_max[0] - xmax, bbox_max[1] - ymax, bbox_max[2] - zmax
	double bbox_min[3];
	double bbox_max[3];

};







#endif
