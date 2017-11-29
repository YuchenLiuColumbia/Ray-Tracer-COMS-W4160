#ifndef LIGHT_H_
#define LIGHT_H_

#include <iostream>
#include <cmath>

#include "basicmath.h"

class pLight
{
	public:

		pLight() {};

		mypoint loc;
		double r, g, b;

		void init(mypoint loc_, double r_, double g_, double b_)
		{
			loc = loc_;
			r = r_;
			g = g_;
			b = b_;
		}

};

class areaLight
{

	public:
		areaLight() {};

		mypoint center, leftdown;
		myvector dir, u_dir, v_dir;
		double len;
		double r, g, b;

		void init(mypoint center_, myvector dir_, myvector u_dir_, double len_, double r_, double g_, double b_)
		{
			center = center_;
			u_dir = u_dir_;
			len = len_;
			dir = dir_;
			r = r_;
			g = g_;
			b = b_;
			v_dir.crossProduct(dir, u_dir);
			
			u_dir.normalize();
			v_dir.normalize();
			dir.normalize();
			
			leftdown.set(center[0] - 0.5*len*u_dir[0] - 0.5*len*v_dir[0], center[1] - 0.5*len*u_dir[1] - 0.5*len*v_dir[1], center[2] - 0.5*len*u_dir[2] - 0.5*len*v_dir[2]);

		}

};

class aLight
{

public:
	aLight() {};

	double r, g, b;

	void init(double r_, double g_, double b_)
	{
		r = r_;
		g = g_;
		b = b_;
	}

};


#endif