// Yuchen Liu, yl3733, l.yuchen@columbia.edu

// Picture Randerer 1.0 - With functions of Setting up Camera, Sphere and Rays
// Most of basicmath.h and parser.cc are used - except some small changes. I will eventually write my own one.

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <ImfRgbaFile.h>
#include <ImfStringAttribute.h>
#include <ImfMatrixAttribute.h>
#include <ImfArray.h>

#include "parser.h"
#include "basicmath.h"
#include "sphere.h"
#include "triangle.h"
#include "ray.h"
#include "light.h"

using namespace std;
using namespace Imf;
using namespace Imath;


int main(int argc, char *argv[])
{

	try
	{

		//initial necessary variables of classes;
		Camera cam;
		vector<Sphere *> asphere;
		vector<Triangle *> atriangle;
		vector<myvector *> asoftnormals;
		vector<pLight *> plight;
		vector<areaLight *> arealight;
		aLight alight;
		Parser p1;
		char *file = argv[1];
		char *BVH;
//		if (argc == 4)
//			BVH = argv[3];
//		else
			BVH = (char*)"1";

		//statement
		cout << file << endl;
		//readin camera & sphere
		p1.parse(file, asphere, atriangle, plight, arealight, alight, asoftnormals, cam);


		int x, y; //x stands to height; y stands to width
		x = cam.ph;
		y = cam.pw;

		int sample, areasample;
		if (argc == 4)
			sample = atoi(argv[3]);

		if (argc == 5)
		{
			sample = atoi(argv[3]);
			areasample = atoi(argv[4]);
		}

		cout << "Uses " << sample << " samples and " << areasample << " area samples" << endl;

		cam.renderscene(asphere, atriangle, plight, arealight, asoftnormals, BVH, sample, areasample);

		cout << "writing output image" << endl;
		cam.writescene(alight, argv[2], y, x);
		cout << "done." << endl;

	}

	catch (const std::exception &exc)
	{
		std::cerr << exc.what() << std::endl;
		return 1;
	}

	return 0;

}



