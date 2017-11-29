#include "bbox.h"

double max(const double x, const double y);
double min(const double x, const double y);

void bbox::buildbox(mypoint p1, mypoint p2, mypoint p3)
{
	//set bounding box
	if ((p1[0] > p2[0]) && (p1[0] > p3[0]))
		bbox_max[0] = p1[0];
	else if (p2[0] > p3[0])
		bbox_max[0] = p2[0];
	else
		bbox_max[0] = p3[0];

	if ((p1[1] > p2[1]) && (p1[1] > p3[1]))
		bbox_max[1] = p1[1];
	else if (p2[1] > p3[1])
		bbox_max[1] = p2[1];
	else
		bbox_max[1] = p3[1];

	if ((p1[2] > p2[2]) && (p1[2] > p3[2]))
		bbox_max[2] = p1[2];
	else if (p2[2] > p3[2])
		bbox_max[2] = p2[2];
	else
		bbox_max[2] = p3[2];

	if ((p1[0] < p2[0]) && (p1[0] < p3[0]))
		bbox_min[0] = p1[0];
	else if (p2[0] < p3[0])
		bbox_min[0] = p2[0];
	else
		bbox_min[0] = p3[0];

	if ((p1[1] < p2[1]) && (p1[1] < p3[1]))
		bbox_min[1] = p1[1];
	else if (p2[1] < p3[1])
		bbox_min[1] = p2[1];
	else
		bbox_min[1] = p3[1];

	if ((p1[2] < p2[2]) && (p1[2] < p3[2]))
		bbox_min[2] = p1[2];
	else if (p2[2] < p3[2])
		bbox_min[2] = p2[2];
	else
		bbox_min[2] = p3[2];

	h[0] = (bbox_max[0] + bbox_min[0]) / 2.0;
	h[1] = (bbox_max[1] + bbox_min[1]) / 2.0;
	h[2] = (bbox_max[2] + bbox_min[2]) / 2.0;

}

void bbox::buildbox(mypoint pos, double r)
{
	//set bounding box
	bbox_min[0] = pos[0] - r;
	bbox_max[0] = pos[0] + r;
	bbox_min[1] = pos[1] - r;
	bbox_max[1] = pos[1] + r;
	bbox_min[2] = pos[2] - r;
	bbox_max[2] = pos[2] + r;

	h[0] = (bbox_max[0] + bbox_min[0]) / 2.0;
	h[1] = (bbox_max[1] + bbox_min[1]) / 2.0;
	h[2] = (bbox_max[2] + bbox_min[2]) / 2.0;
}

bool bbox::hitbox(ray &thisray, double &int_near_r, double &int_far_r)
{
	// if it's a hit, then the int_near & int_far would be changed;
	// otherwise, return the original value

	myvector dd = thisray.dir;
	mypoint pp = thisray.pos;
	double int_near = int_near_r;
	double int_far = int_far_r;
	
	//check bounding box
	if (int_near > int_far)
		return false;
	double tempmax, tempmin;

	for (int i = 0; i < 3; i++)
	{
		if (dd[i] > 0)
		{
			tempmin = (bbox_min[i] - pp[i]) / dd[i];
			if (tempmin > int_near)
				int_near = tempmin;
			tempmax = (bbox_max[i] - pp[i]) / dd[i];
			if (tempmax < int_far)
				int_far = tempmax;
		}
		else if (dd[i] < 0)
		{
			tempmin = (bbox_max[i] - pp[i]) / dd[i];
			if (tempmin > int_near)
				int_near = tempmin;
			tempmax = (bbox_min[i] - pp[i]) / dd[i];
			if (tempmax < int_far)
				int_far = tempmax;
		}
		else
		{
			if ((bbox_max[i] - pp[i]) * (bbox_min[i] - pp[i]) >= 0)
				return false;
		}

		if (int_near > int_far)
			return false;
	}
	int_near_r = int_near;
	int_far_r = int_far;
	return true;
}

void bbox::combinebox(bbox box1, bbox box2)
{

	bbox_max[0] = max(box1.bbox_max[0], box2.bbox_max[0]) + 0.0001;
	bbox_max[1] = max(box1.bbox_max[1], box2.bbox_max[1]) + 0.0001;
	bbox_max[2] = max(box1.bbox_max[2], box2.bbox_max[2]) + 0.0001;
	bbox_min[0] = min(box1.bbox_min[0], box2.bbox_min[0]) - 0.0001;
	bbox_min[1] = min(box1.bbox_min[1], box2.bbox_min[1]) - 0.0001;
	bbox_min[2] = min(box1.bbox_min[2], box2.bbox_min[2]) - 0.0001;

	h[0] = (bbox_max[0] + bbox_min[0]) / 2.0;
	h[1] = (bbox_max[1] + bbox_min[1]) / 2.0;
	h[2] = (bbox_max[2] + bbox_min[2]) / 2.0;
	
}

//not just inside box; also it's goint OUTSIDE OF the box
bool bbox::insidebox(const ray thisray)
{
	mypoint pos = thisray.pos;
	myvector dir = thisray.dir;

	if ((bbox_min[0] < pos[0]) && (bbox_max[0] > pos[0]))
		if ((bbox_min[1] < pos[1]) && (bbox_max[1] > pos[1]))
			if ((bbox_min[2] < pos[2]) && (bbox_max[2] > pos[2]))
			{
				myvector dir2;
				dir2.set(pos[0] - h[0], pos[1] - h[1], pos[2] - h[2]);
				if (dir2.dotProduct(dir) > 0) //means that the ray is heading outside of box
					return true;
			}
	
	return false;
}

double min(const double x, const double y)
{
	if (x >= y)
		return y;
	else
		return x;
}

double max(const double x, const double y)
{
	if (x <= y)
		return y;
	else
		return x;
}