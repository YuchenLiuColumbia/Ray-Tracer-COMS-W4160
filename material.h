#ifndef MATERIAL_H_
#define MATERIAL_H_


#include <iostream>
#include <cmath>
#include <vector>

#include "basicmath.h"
#include "ray.h"
#include "light.h"

class Material
{
	public:
		Material() {};

		double dr, dg, db;

		double sr, sg, sb;

		double ir, ig, ib;

		double r;

		void init(double dr_, double dg_, double db_, double sr_, double sg_, double sb_, double r_, double ir_, double ig_, double ib_)
		{
			dr = dr_;
			dg = dg_;
			db = db_;
			sr = sr_;
			sg = sg_;
			sb = sb_;
			ir = ir_;
			ig = ig_;
			ib = ib_;
			r = r_;
		}


		myvector RGB(const myvector &ptl,
			const myvector &ptv,
			const myvector &norm,
			const myvector &shadow_v);

		myvector BackRGB(const myvector &ptl,
			const myvector &ptv,
			const myvector &norm,
			const myvector &shadow_v);

};


#endif