#ifndef PARSE_H
#define PARSE_H

#include <iostream>
#include <vector>

#include "sphere.h"
#include "camera.h"
#include "light.h"
#include "triangle.h"


class Parser 
{
	public:

		virtual void parse(const char *file,
			std::vector<Sphere *> &asphere,
			std::vector<Triangle *> &atriangle,
			std::vector<pLight *> &plight,
			std::vector<areaLight *> &arealight,
			aLight &alight,
			std::vector<myvector *> &asoftnormals,
			Camera &cam
			);

		void read_wavefront_file(const char *file,
								 std::vector< int > &tris,
								 std::vector< float > &verts,
								 std::vector<myvector *> &asoftnormals
								 );

};

#endif
