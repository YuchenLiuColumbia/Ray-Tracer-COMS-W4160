#ifndef CAMERA_H_
#define CAMERA_H_


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "basicmath.h"
#include "ray.h"
#include "material.h"
#include "triangle.h"
#include "sphere.h"

#include <ImfRgbaFile.h>
#include <ImfStringAttribute.h>
#include <ImfMatrixAttribute.h>
#include <ImfArray.h>


class Camera 
{

	public:

		mypoint eye;

		double d;

		myvector u;
		myvector v;
		myvector w;

		double width;
		double height;
		int pw;
		int ph;

		double nx;
		double ny;

		mypoint ImgVerticalPoint;
		mypoint leftup;

		//function: initial data readin
		void init(mypoint position, 
			      myvector direction, 
				  double focallength, 
				  double imagew, 
				  double imageh, 
				  int pixelw, 
				  int pixelh
				  );

		//for spread ray

		void renderscene(std::vector<Sphere *> &asphere,
						 std::vector<Triangle *> &atriangle,
						 std::vector<pLight *> &plight,
						 std::vector<areaLight *> &arealight,
						 std::vector<myvector *> &asoftnormals,
						 const char BVH[],
						 int samples,
						 int areasamples
						 );
		
		//for initial ray
		ray createray(int ix, int iy);

		//for stratified rays
		ray createray(int ix, int iy, int p, int q, int n);

		//for true reflection
		ray createray(mypoint startpoint, myvector direction);


		void setPixel(int px, int py, myvector rgb_this) {
			Imf::Rgba &pixel = lightimg[px][py];
			pixel.r += rgb_this[0];
			pixel.g += rgb_this[1];
			pixel.b += rgb_this[2];
			pixel.a = 1.0;
		}


		void writescene(aLight alight, const char filename[], int y, int x);

		myvector calcRGB(ray &thisray,
			int ray_type,
			int certainlight,
			int spread_time,
			std::vector<Sphere *> &asphere,
			std::vector<Triangle *> &atriangle,
			std::vector<pLight *> &plight,
			std::vector<areaLight *> &arealight,
			double mint,
			double maxt);

		myvector FcalcRGB(ray &thisray,
			int ray_type,
			int light_type,
			int certainlight,
			int spread_time,
			Sphere *asphere,
			Triangle *atriangle,
			bool ifBVHsphere,
			bool ifBVHtriangle,
			std::vector<pLight *> &plight,
			std::vector<areaLight *> &arealight,
			std::vector<myvector *> &asoftnormals,
			double mint,
			double maxt,
			int areasamples);


	private:
		myvector stdh = myvector(0, 1, 0);
		Imf::Array2D<Imf::Rgba> lightimg;

};

#endif
