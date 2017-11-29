#include <vector>
#include <cmath>


#include "basicmath.h"
#include "camera.h"
#include "ray.h"
#include "sphere.h"
#include "triangle.h"
#include "light.h"


void Camera::init(mypoint position, 
				  myvector direction, 
				  double focallength, 
				  double imagew, 
				  double imageh, 
				  int pixelw, 
				  int pixelh
				  )
		{
			//readin
			eye = position;
			d = focallength;
			width = imagew;
			height = imageh;
			pw = pixelw;
			ph = pixelh;
			ny = width / pw;
			nx = height / ph;
			w = direction.invert();
			w.normalize();
			u.crossProduct(direction, stdh);
			v.crossProduct(u, direction);

			//normalize
			u.normalize();
			v.normalize();
			w.normalize();


			//calculate vertical point
			ImgVerticalPoint.set(eye[0] - w[0] * d, eye[1] + direction[1] - w[1] * d, eye[2] - w[2] * d);

			//calculate leftup point
			if (ph % 2 == 0)
				leftup.set(ImgVerticalPoint[0] + nx*(ph*0.5 - 0.5)*v[0], ImgVerticalPoint[1] + nx*(ph*0.5 - 0.5)*v[1], ImgVerticalPoint[2] + nx*(ph*0.5 - 0.5)*v[2]);
			else
				leftup.set(ImgVerticalPoint[0] + nx*(ph*0.5 - 1)*v[0], ImgVerticalPoint[1] + nx*(ph*0.5 - 1)*v[1], ImgVerticalPoint[2] + nx*(ph*0.5 - 1)*v[2]);

			if (pw % 2 == 0)
				leftup.set(leftup[0] - ny*(pw*0.5 - 0.5)*u[0], leftup[1] - ny*(pw*0.5 - 0.5)*u[1], leftup[2] - ny*(pw*0.5 - 0.5)*u[2]);
			else
				leftup.set(leftup[0] - ny*(pw*0.5 - 1)*u[0], leftup[1] - ny*(pw*0.5 - 1)*u[1], leftup[2] - ny*(pw*0.5 - 1)*u[2]);

			//create scene
			lightimg.resizeErase(ph, pw); //allocate memory and set to all 0
			int i, j;
			for (i = 0; i < ph; i++)
				for (j = 0; j < pw; j++)
				{
					Imf::Rgba &px = lightimg[i][j];
					{
						px.r = 0;
						px.g = 0;
						px.b = 0;
						px.a = 1;
					}
				}

		}



void Camera::renderscene(std::vector<Sphere *> &asphere,
						 std::vector<Triangle *> &atriangle,
					   	 std::vector<pLight *> &plight,
						 std::vector<areaLight *> &arealight,
						 std::vector<myvector *> &asoftnormals,
						 const char BVH[],
						 int samples,
						 int areasamples
						 )
{
	int x, y; //x stands to height; y stands to width
	int i, j;
	int nn = samples * samples;
	x = ph;
	y = pw;
	int spread_time = 3;
	bool useBVH = true;
	bool sphereBVH = false;
	bool triangleBVH = false;
	Sphere *thissphere = new Sphere;
	Triangle *thistriangle = new Triangle;

	if (BVH[0] == '1')
	{
		if (asphere.size() > 0)
		{
			thissphere->createtree(asphere, 0);
			sphereBVH = true;
		}
		if (atriangle.size() > 0)
		{
			thistriangle->createtree(atriangle, 0);
			triangleBVH = true;
		}
	}
	else
		useBVH = false;

	if(nn == 0)
	for (i = 0; i < x; i++)
	{
		for (j = 0; j < y; j++)
		{
			//create ray of certain pixel
			ray thisray = createray(i, j);

			myvector rgb_this;
			rgb_this.set(0, 0, 0);
			if (!useBVH)
				rgb_this = calcRGB(thisray,0, 0, spread_time, asphere, atriangle, plight, arealight, 0, std::numeric_limits<double>::max());
			else
				rgb_this = FcalcRGB(thisray, 0, 0, 0, spread_time, thissphere, thistriangle, sphereBVH, triangleBVH, plight, arealight, asoftnormals, 0, std::numeric_limits<double>::max(), areasamples);
			setPixel(i, j, rgb_this);

//			std::cout << "row" << i << "column" << j << "finished" << std::endl;

			if (i == 272 && j == 458)
				continue;


		}

		if (i % 50 == 0)
			std::cout << "line" << i << "finished" << std::endl;
	}
	
	else
		for (i = 0; i < x; i++)
		{
			for (j = 0; j < y; j++)
			{
				myvector rgb_this;
				rgb_this.set(0, 0, 0);

				for (int p = 0; p < samples; p++)
					for (int q = 0; q < samples; q++)
					{
						ray thisray = createray(i, j, p, q, samples);

						myvector rgb_temp;
						rgb_temp = FcalcRGB(thisray, 0, 0, 0, spread_time, thissphere, thistriangle, sphereBVH, triangleBVH, plight, arealight, asoftnormals, 0, std::numeric_limits<double>::max(), areasamples);
						rgb_this.set(rgb_this[0] + rgb_temp[0], rgb_this[1] + rgb_temp[1], rgb_this[2] + rgb_temp[2]);
					}

				//now we have got the total rgb_this for all stratifies
				rgb_this.set(rgb_this[0] / nn, rgb_this[1] / nn, rgb_this[2] / nn);
				setPixel(i, j, rgb_this);

				if (i == 272 && j == 458)
					continue;

			}
			
			if (i % 50 == 0)
			std::cout << "row" << i << "finished" << std::endl;
		}
}

void Camera::writescene(aLight alight, const char filename[], int y, int x)
{
	int i, j;
	for (i = 0; i < x; i++)
		for (j = 0; j < y; j++)
		{
			Imf::Rgba &px = lightimg[i][j];
			{
//				px.r += alight.r;
//				px.g += alight.g;
//				px.b += alight.b;
				px.a = 1;


			}
		}
	Imf::Rgba &pixels = lightimg[0][0];
	Imf::RgbaOutputFile file(filename, y, x, Imf::WRITE_RGBA);
	file.setFrameBuffer(&pixels, 1, y);
	file.writePixels(x);
}

ray Camera::createray(int ix, int iy)
{

	mypoint locnow;
	ray theray;

	locnow = leftup;

	//set local screen point
	locnow.set(locnow[0] + iy*ny*u[0], locnow[1] + iy*ny*u[1], locnow[2] + iy*ny*u[2]);
	locnow.set(locnow[0] - ix*nx*v[0], locnow[1] - ix*nx*v[1], locnow[2] - ix*nx*v[2]);

	theray.pos = eye;
	theray.dir = locnow - eye;
	theray.dir.normalize();

	return theray;

}

ray Camera::createray(int ix, int iy, int p, int q, int n)
{
	mypoint locnow;
	ray theray;
	double randnum1 = rand()%10;
	double randnum2 = rand()%10;
	randnum1 /= 10;
	randnum2 /= 10;

	locnow = leftup;

	//set local screen point; this time, left up of each pixel
	locnow.set(locnow[0] + (iy - 0.5 + (p + randnum1)/n)*ny*u[0], locnow[1] + (iy - 0.5 + (p + randnum1)/n)*ny*u[1], locnow[2] + (iy - 0.5 + (p + randnum1)/n)*ny*u[2]);
	locnow.set(locnow[0] - (ix - 0.5 + (q + randnum2)/n)*nx*v[0], locnow[1] - (ix - 0.5 + (q + randnum2)/n)*nx*v[1], locnow[2] - (ix - 0.5 + (q + randnum2)/n)*nx*v[2]);
	
	theray.pos = eye;
	theray.dir = locnow - eye;
	theray.dir.normalize();

	return theray;
}



ray Camera::createray(mypoint startpoint, myvector direction)
{
	ray theray;
	theray.pos = startpoint;
	theray.dir = direction;
	theray.dir.normalize();

	return theray;
}

//default method for ray-tracing.
//tear of time - never been used while implementing soft shadows
myvector Camera::calcRGB(ray &thisray,
	int ray_type,
	int certainlight,
	int spread_time,
	std::vector<Sphere *> &asphere,
	std::vector<Triangle *> &atriangle,
	std::vector<pLight *> &plight,
	std::vector<areaLight *> &arealight,
	double mint,
	double maxt)
{
	myvector rgb_this;
	rgb_this.set(0, 0, 0);

	if (spread_time == 0)
		return rgb_this;

	int aspherelength = asphere.size();
	int atrilength = atriangle.size();

	if (ray_type == 2)
	{
		for (int m = 0; m < aspherelength; m++)
		{
			Sphere *thissphere = asphere[m];
			double t;
			mypoint pointtemp;
			myvector normaltemp;
			if (thissphere->intersect(t, thisray, pointtemp, normaltemp) && (t < maxt && t > mint))
				return rgb_this;
		}

		double tmax = std::numeric_limits<double>::max();

		for (int m = 0; m < atrilength; m++)
		{
			Triangle *thistriangle = atriangle[m];
			double t;
			mypoint pointtemp;
			myvector normaltemp;
			if (thistriangle->intersect(t, tmax, thisray, pointtemp, normaltemp) && (t < maxt && t > mint))
				return rgb_this;
		}

		rgb_this.set(plight[certainlight]->r, plight[certainlight]->g, plight[certainlight]->b);
		return rgb_this;
	}


	double tsmin = std::numeric_limits<double>::max(); //the largest number
	mypoint intersectpoint;
	myvector n_normal;
	Material thismaterial;

	//judge surface intersection
	int frontsurface = -1;
	for (int m = 0; m < aspherelength; m++)
	{
		Sphere *thissphere = asphere[m];
		double t;
		mypoint pointtemp;
		myvector normaltemp;

		if (thissphere->intersect(t, thisray, pointtemp, normaltemp))
			if (t < tsmin)
			{
				frontsurface = m;
				tsmin = t;
				intersectpoint = pointtemp;
				n_normal = normaltemp;
				thismaterial = thissphere->material;
			}
	}

	int fronttri = -1;
	for (int m = 0; m < atrilength; m++)
	{
		Triangle *thistriangle = atriangle[m];
		double t;
		mypoint pointtemp;
		myvector normaltemp;

		if (thistriangle->intersect(t, tsmin, thisray, pointtemp, normaltemp))
			if (t < tsmin)
			{
				frontsurface = m;
				tsmin = t;
				intersectpoint = pointtemp;
				n_normal = normaltemp;
				thismaterial = thistriangle->material;
			}
	}


	if (frontsurface == -1 && fronttri == -1)
		return rgb_this;


	//shadow calculation
	bool ifbackside = false;
	if (thisray.dir.dotProduct(n_normal) > 0)
		ifbackside = true;
	myvector backsidenorm = n_normal.invert();
	myvector ptv = thisray.dir;
	ptv = ptv.invert();

	int plightlength = plight.size();
	for (int n = 0; n < plightlength; n++)
	{
		//set up a vector from the light to the intersectpoint

		pLight *thislight = plight[n];
		myvector ptl = thislight->loc - intersectpoint;

		double thistsq = ptl[0] * ptl[0] + ptl[1] * ptl[1] + ptl[2] * ptl[2];
		double thist = sqrt(thistsq);

		ptl.normalize();

		ray shadowRay = createray(intersectpoint, ptl);
		myvector shadow_v = calcRGB(shadowRay, 2, n, 1, asphere, atriangle, plight, arealight, 0.001, thist);

		if (shadow_v[0] != 0 || shadow_v[1] != 0 || shadow_v[2] != 0)
		{
			myvector rgb_temp;
			rgb_temp.set(0, 0, 0);

			if (ifbackside)
				rgb_temp = thismaterial.BackRGB(ptl, ptv, backsidenorm, shadow_v);
			else
				rgb_temp = thismaterial.RGB(ptl, ptv, n_normal, shadow_v);

			rgb_this.set(rgb_this[0] + rgb_temp[0] / thistsq, rgb_this[1] + rgb_temp[1] / thistsq, rgb_this[2] + rgb_temp[2] / thistsq);

		}
	}


	//true_specular light
	if (thismaterial.ir == 0 && thismaterial.ig == 0 && thismaterial.ib == 0)
		return rgb_this;
	else if (ifbackside)
		return rgb_this;
	else
	{
		double d_N = 2 * thisray.dir.dotProduct(n_normal);
		myvector direction;
		myvector rgb_temp;
		mypoint iipoint;
		rgb_temp.set(0, 0, 0);

		direction.set(thisray.dir[0] - d_N * n_normal[0], thisray.dir[1] - d_N * n_normal[1], thisray.dir[2] - d_N * n_normal[2]);
		direction.normalize();
		iipoint.set(intersectpoint[0] + direction[0] * 0.001, intersectpoint[1] + direction[1] * 0.001, intersectpoint[2] + direction[2] * 0.001);
		ray spreadray = createray(iipoint, direction);
		rgb_temp = calcRGB(spreadray, 1, 0, spread_time - 1, asphere, atriangle, plight, arealight, 0.001, std::numeric_limits<double>::max());
		rgb_this.set(rgb_this[0] + rgb_temp[0] * thismaterial.ir, rgb_this[1] + rgb_temp[1] * thismaterial.ig, rgb_this[2] + rgb_temp[2] * thismaterial.ib);
		return rgb_this;
	}

}


//Bounding box, BVH tree, Monte-carlo method
//Soft shadows, smooth normals
myvector Camera::FcalcRGB(ray &thisray,
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
	int areasamples)
{
	myvector rgb_this;
	rgb_this.set(0, 0, 0);

	if (spread_time == 0)
		return rgb_this;


	if (ray_type == 2)
	{

		std::vector<Sphere *> thissphere;
		double t_n = 0, t_f = std::numeric_limits<double>::max();
		if (ifBVHsphere)
			asphere->hit(thisray, t_n, t_f, thissphere);
		if (thissphere.size() > 0)
		{

			int lengtha = thissphere.size();
			for (int m = 0; m < lengtha; m++)
			{
				double t;
				mypoint pointtemp;
				myvector normaltemp;
				Sphere *tsphere = thissphere[m];

				if (tsphere->intersect(t, thisray, pointtemp, normaltemp) && (t < maxt && t > mint))
					return rgb_this;
			}
		}

		double tmax = std::numeric_limits<double>::max();
		std::vector<Triangle *> thistriangle;
		t_n = 0, t_f = std::numeric_limits<double>::max();
		if (ifBVHtriangle)
			atriangle->hit(thisray, t_n, t_f, thistriangle);
		if (thistriangle.size() > 0)
		{

			int lengtha = thistriangle.size();
			for (int m = 0; m < lengtha; m++)
			{
				double t;
				mypoint pointtemp;
				myvector normaltemp;
				Triangle *ttri = thistriangle[m];

				if (ttri->intersect(t, tmax, thisray, pointtemp, normaltemp) && (t < maxt && t > mint))
					return rgb_this;
			}
		}
		if (light_type == 1)
			rgb_this.set(plight[certainlight]->r, plight[certainlight]->g, plight[certainlight]->b);
		else
			rgb_this.set(arealight[certainlight]->r, arealight[certainlight]->g, arealight[certainlight]->b);
		
		return rgb_this;
	}


	double tsmin = std::numeric_limits<double>::max(); //the largest number
	mypoint intersectpoint;
	myvector n_normal;
	Material thismaterial;

	//judge surface intersection
	bool frontsphere = false;
	std::vector<Sphere *> thissphere;
	double t_n = 0, t_f = std::numeric_limits<double>::max();
	if (ifBVHsphere)
		asphere->hit(thisray, t_n, t_f, thissphere);
	if (thissphere.size() > 0)
	{
		int lengtha = thissphere.size();
		for (int m = 0; m < lengtha; m++)
		{
			double t;
			mypoint pointtemp;
			myvector normaltemp;
			Sphere *tsphere = thissphere[m];
			if (tsphere->intersect(t, thisray, pointtemp, normaltemp) && (t < tsmin))
			{
				frontsphere = true;
				tsmin = t;
				intersectpoint = pointtemp;
				n_normal = normaltemp;
				thismaterial = tsphere->material;
			}
		}
	}

	bool fronttri = false;
	std::vector<Triangle *> thistriangle;
	t_n = 0, t_f = std::numeric_limits<double>::max();
	if (ifBVHtriangle)
		atriangle->hit(thisray, t_n, t_f, thistriangle);
	if (thistriangle.size() > 0)
	{
		int lengtha = thistriangle.size();
		for (int m = 0; m < lengtha; m++)
		{
			double t;
			mypoint pointtemp;
			myvector normaltemp;
			Triangle *ttri = thistriangle[m];
			if (ttri->intersect(t, tsmin, thisray, pointtemp, normaltemp, asoftnormals) && (t < tsmin))
			{
				fronttri = true;
				tsmin = t;
				intersectpoint = pointtemp;
				n_normal = normaltemp;
				thismaterial = ttri->material;
			}
		}
	}


	if ((!frontsphere) && (!fronttri))
		return rgb_this;


	//shadow calculation
	bool ifbackside = false;
	if (thisray.dir.dotProduct(n_normal) > 0)
		ifbackside = true;
	myvector backsidenorm = n_normal.invert();
	myvector ptv = thisray.dir;
	ptv = ptv.invert();

	int plightlength = plight.size();
	for (int n = 0; n < plightlength; n++)
	{
		//set up a vector from the intersectpoint to the light

		pLight *thislight = plight[n];
		myvector ptl = thislight->loc - intersectpoint;

		double thistsq = ptl[0] * ptl[0] + ptl[1] * ptl[1] + ptl[2] * ptl[2];
		double thist = sqrt(thistsq);

		ptl.normalize();

		ray shadowRay = createray(intersectpoint, ptl);
		myvector shadow_v = FcalcRGB(shadowRay, 2, 1, n, 1, asphere, atriangle, ifBVHsphere,ifBVHtriangle, plight, arealight, asoftnormals, 0.001, thist, areasamples);

		if (shadow_v[0]!=0 || shadow_v[1]!=0 || shadow_v[2]!=0)
		{
			myvector rgb_temp;
			rgb_temp.set(0, 0, 0);

			if (ifbackside)
				rgb_temp = thismaterial.BackRGB(ptl, ptv, backsidenorm, shadow_v);
			else
				rgb_temp = thismaterial.RGB(ptl, ptv, n_normal, shadow_v);
			
			rgb_this.set(rgb_this[0] + rgb_temp[0] / thistsq, rgb_this[1] + rgb_temp[1] / thistsq, rgb_this[2] + rgb_temp[2] / thistsq);

		}
	}

	int arealightlength = arealight.size();
	int areasum = areasamples * areasamples;
	for (int n = 0; n < arealightlength; n++)
	{
		//if we are not calculating soft shadow;
		//set up a vector from the intersectpoint to the light's center
		
		if (areasamples == 0)
		{			
			areaLight *thislight = arealight[n];
			myvector ptl = thislight->center - intersectpoint;

			double thistsq = ptl[0] * ptl[0] + ptl[1] * ptl[1] + ptl[2] * ptl[2];
			double thist = sqrt(thistsq);

			ptl.normalize();
			
			// calculate if the ptl vector goes to the facing side of area light
			// if does, then calculate; if doesn't, set to 0
			double facetest = ptl.dotProduct(thislight->dir);
			if (facetest > 0)
				continue;
			
			ray shadowRay = createray(intersectpoint, ptl);
			myvector shadow_v = FcalcRGB(shadowRay, 2, 2, n, 1, asphere, atriangle, ifBVHsphere,ifBVHtriangle, plight, arealight, asoftnormals, 0.1, thist, areasamples);

			if (shadow_v[0]!=0 || shadow_v[1]!=0 || shadow_v[2]!=0)
			{
				myvector rgb_temp;
				rgb_temp.set(0, 0, 0);

				if (ifbackside)
					rgb_temp = thismaterial.BackRGB(ptl, ptv, backsidenorm, shadow_v);
				else
					rgb_temp = thismaterial.RGB(ptl, ptv, n_normal, shadow_v);
			
				rgb_this.set(rgb_this[0] + rgb_temp[0] / thistsq, rgb_this[1] + rgb_temp[1] / thistsq, rgb_this[2] + rgb_temp[2] / thistsq);

			}
		}
		
		else
		{			
			areaLight *thislight = arealight[n];
							
			// for this number might be big, so it would be much easier
			// to calculate if its shined by area light
			myvector ptl_test = thislight->center - intersectpoint;
			double face_test = ptl_test.dotProduct(thislight->dir);
				if (face_test > 0)
					continue;
			
			// now we know it could be shined by this area light
			for (int pp = 0; pp < areasamples; pp++)
				for(int qq = 0; qq < areasamples; qq++)
				{
					mypoint lightfacepoint;
					double randnum1 = rand()%10;
					double randnum2 = rand()%10;
					randnum1 /= 10;
					randnum2 /= 10;
					double percent1 = (pp + randnum1) / areasamples;
					double percent2 = (qq + randnum2) / areasamples;
					lightfacepoint.set(thislight->leftdown[0] + percent1*thislight->len*thislight->u_dir[0] + percent2*thislight->len*thislight->v_dir[0],
									   thislight->leftdown[1] + percent1*thislight->len*thislight->u_dir[1] + percent2*thislight->len*thislight->v_dir[1],
									   thislight->leftdown[2] + percent1*thislight->len*thislight->u_dir[2] + percent2*thislight->len*thislight->v_dir[2]);
					myvector ptl = lightfacepoint - intersectpoint;

					double thistsq = ptl[0] * ptl[0] + ptl[1] * ptl[1] + ptl[2] * ptl[2];
					double thist = sqrt(thistsq);

					ptl.normalize();
			

			
					ray shadowRay = createray(intersectpoint, ptl);
					myvector shadow_v = FcalcRGB(shadowRay, 2, 2, n, 1, asphere, atriangle, ifBVHsphere,ifBVHtriangle, plight, arealight, asoftnormals, 0.1, thist, areasamples);

					if (shadow_v[0]!=0 || shadow_v[1]!=0 || shadow_v[2]!=0)
					{
						myvector rgb_temp;
						rgb_temp.set(0, 0, 0);

						if (ifbackside)
							rgb_temp = thismaterial.BackRGB(ptl, ptv, backsidenorm, shadow_v);
						else
							rgb_temp = thismaterial.RGB(ptl, ptv, n_normal, shadow_v);
			
						rgb_this.set(rgb_this[0] + rgb_temp[0] / (thistsq * areasum), 
									 rgb_this[1] + rgb_temp[1] / (thistsq * areasum), 
									 rgb_this[2] + rgb_temp[2] / (thistsq * areasum));

					}
				}
		}
	}
	
	
	
	

	//true_specular light
	if (thismaterial.ir == 0 && thismaterial.ig == 0 && thismaterial.ib == 0)
		return rgb_this;
	else if (ifbackside)
		return rgb_this;
	else
	{
		double d_N = 2 * thisray.dir.dotProduct(n_normal);
		myvector direction;
		myvector rgb_temp;
		mypoint iipoint;
		rgb_temp.set(0, 0, 0);

		direction.set(thisray.dir[0] - d_N * n_normal[0], thisray.dir[1] - d_N * n_normal[1], thisray.dir[2] - d_N * n_normal[2]);
		direction.normalize();
		iipoint.set(intersectpoint[0] + direction[0] * 0.001, intersectpoint[1] + direction[1] * 0.001, intersectpoint[2] + direction[2] * 0.001);
		ray spreadray = createray(iipoint, direction);
		rgb_temp = FcalcRGB(spreadray, 1, 0, 0, spread_time - 1, asphere, atriangle, ifBVHsphere, ifBVHtriangle, plight, arealight, asoftnormals, 0.001, std::numeric_limits<double>::max(), areasamples);
		rgb_this.set(rgb_this[0] + rgb_temp[0] * thismaterial.ir, rgb_this[1] + rgb_temp[1] * thismaterial.ig, rgb_this[2] + rgb_temp[2] * thismaterial.ib);
		return rgb_this;
	}

}
