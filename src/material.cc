
#include "material.h"

double max(const double x, const double y);

myvector Material::RGB(const myvector &ptl,
	const myvector &ptv,
	const myvector &norm,
	const myvector &shadow_v)
{

	myvector rgb_this;
	rgb_this.set(0, 0, 0);

	myvector h_nor; //h normal

	h_nor.set(ptl[0] + ptv[0], ptl[1] + ptv[1], ptl[2] + ptv[2]);

    //calculate diffusion light
	double backsidejudge = norm.dotProduct(ptl);
	if (backsidejudge > 0)
	{
		double localr, localg, localb;
		localr = shadow_v[0] * dr * backsidejudge;
		localg = shadow_v[1] * dg * backsidejudge;
		localb = shadow_v[2] * db * backsidejudge;
		rgb_this.set(rgb_this[0] + localr, rgb_this[1] + localg, rgb_this[2] + localb);
	}
	else
		return rgb_this;

	//calculate specular light
	if (h_nor[0] != 0 || h_nor[1] != 0 || h_nor[2] != 0)
	{
		h_nor.normalize();
		double x1 = norm.dotProduct(h_nor);
		if (x1 > 0)
		{
			double phongexp;
			double y1 = r;
			phongexp = pow(x1, y1);
			double localr, localg, localb;
			localr = sr * shadow_v[0] * phongexp;
			localg = sg * shadow_v[1] * phongexp;
			localb = sb * shadow_v[2] * phongexp;
			rgb_this.set(rgb_this[0] + localr, rgb_this[1] + localg, rgb_this[2] + localb);
		}
	}

	return rgb_this;
}

myvector Material::BackRGB(const myvector &ptl,
	const myvector &ptv,
	const myvector &norm,
	const myvector &shadow_v)
{
	myvector rgb_this;
	rgb_this.set(0, 0, 0);

	double backsidejudge = norm.dotProduct(ptl);
	if (backsidejudge > 0)
	{
		double localr, localg, localb;
		localr = shadow_v[0] * 1.0 * backsidejudge;
		localg = shadow_v[1] * 1.0 * backsidejudge;
		localb = shadow_v[2] * 0.0 * backsidejudge;
		rgb_this.set(rgb_this[0] + localr, rgb_this[1] + localg, rgb_this[2] + localb);
	}

	return rgb_this;
}
