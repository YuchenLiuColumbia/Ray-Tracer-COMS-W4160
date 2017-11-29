#ifndef SPHERE_H_
#define SPHERE_H_


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include "basicmath.h"
#include "material.h"
#include "ray.h"
#include "bbox.h"

class Sphere
{

	public:

		Sphere() {};

		mypoint pos;
		double r;
		Material material;
		
		void init(mypoint position, double radius, Material thismaterial);
		
		bool intersect(double &t, ray &thisray, mypoint &ispoint, myvector &normal);

/*----------------------FROM-HERE-IS-FOR-BVH-TREE----------------------------------*/


		bbox thisbox;

		void createtree(std::vector<Sphere *> &asphere, int axis);

		bool hit(ray &thisray, double &int_near, double &int_far, std::vector<Sphere *> &record);

	private:
		bool ifleaf = false;

		Sphere *left = nullptr;
		Sphere *right = nullptr;


};

#endif