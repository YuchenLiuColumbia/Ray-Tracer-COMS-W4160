#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>

#include "basicmath.h"
#include "material.h"
#include "ray.h"
#include "bbox.h"


class Triangle
{

	public:
		
		Triangle() {};

		mypoint p1, p2, p3;
		int vert1, vert2, vert3;
		Material material;

		bool isobjtri;

		void init(mypoint p1_, mypoint p2_, mypoint p3_, Material thismaterial);

		myvector objinit(mypoint p1_, mypoint p2_, mypoint p3_, int vert1_, int vert2_, int vert3_, Material thismaterial);

		bool intersect(double &t, double &tmin, ray &thisray, mypoint &thispoint, myvector &normal_out);
		bool intersect(double &t, double &tmin, ray &thisray, mypoint &thispoint, myvector &normal_out, std::vector<myvector *> &asoftnormals);
/*----------------------FROM-HERE-IS-FOR-BVH-TREE----------------------------------*/



		bbox thisbox;

		void createtree(std::vector<Triangle *> &atriangle, int axis);

		bool hit(ray &thisray, double &int_near, double &int_far, std::vector<Triangle*> &record);

	private:
		bool ifleaf = false;
		Triangle *left = nullptr;
		Triangle *right = nullptr;

};


#endif