

#include "sphere.h"

using namespace std;

bool comparex(const Sphere *s1, const Sphere *s2);
bool comparey(const Sphere *s1, const Sphere *s2);
bool comparez(const Sphere *s1, const Sphere *s2);

void Sphere::init(mypoint position, double radius, Material thismaterial)
{
	pos = position;
	r = radius;
	material = thismaterial;
	thisbox.buildbox(position, radius);
	left = nullptr;
	right = nullptr;
	ifleaf = true;

	return;
}

bool Sphere::intersect(double &t, ray &thisray, mypoint &ispoint, myvector &normal)
{

	myvector dd = thisray.dir;
	mypoint ee = thisray.pos;
	mypoint cc = pos;
	myvector eminusc = ee - cc;
	double r_n2 = r * r;

	// (lf +- sqrt(rf1 - rf2 * rf3)) / denom
	//calculate denominator
	double denom = dd.dotProduct(dd);

	//calculate left fraction
	double lf = dd.dotProduct(eminusc);
	lf = -lf;

	//calculate right fraction
	double rf1 = dd.dotProduct(eminusc);
	rf1 = rf1 * rf1;

	double rf2 = dd.dotProduct(dd);
	double rf3 = eminusc.dotProduct(eminusc);
	rf3 = rf3 - r_n2;

	double rf = rf1 - rf2*rf3;
	if (rf < 0)
		return false;

	rf = sqrt(rf);

	double t1 = (lf + rf) / denom;
	if (t1 < 0) return false;

	double t2 = (lf - rf) / denom;
	if (t2 < 0) return false;

	double tempt = min(t1, t2);
	t = tempt;

	ispoint.set(thisray.pos[0] + thisray.dir[0] * t, thisray.pos[1] + thisray.dir[1] * t, thisray.pos[2] + thisray.dir[2] * t);
	normal = ispoint - pos;
	normal.normalize();

	return true;
}

void Sphere::createtree(std::vector<Sphere *> &asphere, int axis)
{
	int aspherelength = asphere.size();
	left = new Sphere;
	right = new Sphere;

	if (aspherelength == 1)
	{
		left = asphere[0];
		right = nullptr;
		thisbox = left->thisbox;
	}
	else if (aspherelength == 2)
	{
		left = asphere[0];
		right = asphere[1];
		thisbox.combinebox(left->thisbox, right->thisbox);
	}
	else
	{
		//sort the sphere class by axis
		if (axis == 0)
			std::nth_element(asphere.begin(), asphere.begin() + asphere.size() / 2, asphere.end(), comparex);
		else if (axis == 1)
			std::nth_element(asphere.begin(), asphere.begin() + asphere.size() / 2, asphere.end(), comparey);
		else 
			std::nth_element(asphere.begin(), asphere.begin() + asphere.size() / 2, asphere.end(), comparez);
		
		vector <Sphere *> leftsphere;
		vector <Sphere *> rightsphere;
		int mid = aspherelength / 2;
		for (int i = 0; i < mid; i++)
		{
			Sphere *thissphere = asphere[i];
			leftsphere.push_back(thissphere);
		}
		for (int i = mid; i < aspherelength; i++)
		{
			Sphere *thissphere = asphere[i];
			rightsphere.push_back(thissphere);
		}
		int newaxis = (axis + 1) % 3;
		left->createtree(leftsphere, newaxis);
		right->createtree(rightsphere, newaxis);
		thisbox.combinebox(left->thisbox, right->thisbox);

	}
}

bool Sphere::hit(ray &thisray, double &int_near, double &int_far, std::vector<Sphere*> &record)
{
	bool left_hit = false, right_hit = false;

	// in this step, if it's a hit, then the int_near would be changed, fitting to the bounding box
	// if it's not hit, then the value would be change
	// if it hits but the final answer is 'not hit' (targeting gap), 
	//the return value would be false and doesn't influence the result
	if (thisbox.hitbox(thisray, int_near, int_far))
	{
		if ((ifleaf) && left == nullptr && right == nullptr)
		{
			if (thisbox.insidebox(thisray))
				return false;
			else
				return true;
		}
		double nn1 = int_near, nn2 = int_near;
		double ff1 = int_far, ff2 = int_far;
		left_hit = ((left != nullptr) && (left->hit(thisray, nn1, ff1, record)));
		if (left_hit)
			record.push_back(left);

		right_hit = ((right != nullptr) && (right->hit(thisray, nn2, ff2, record)));
		if (right_hit)
			record.push_back(right);

/*		if (left_hit && right_hit)
		{
			if (left->ifleaf)
				lrecord.push_back(left);
			if (right->ifleaf)
				rrecord.push_back(right);

			int llength = lrecord.size();
			int rlength = rrecord.size();

			for (int i = 0; i < llength; i++)
			{
				Sphere *tosphere = lrecord[i];
				record.push_back(tosphere);
			}
			for (int i = 0; i < rlength; i++)
			{
				Sphere *tosphere = rrecord[i];
				record.push_back(tosphere);
			}
			int_near = nn2;
			int_far = ff2;
			return true;
		}
		else if (left_hit)
		{
			if (left->ifleaf)
				lrecord.push_back(left);

			int llength = lrecord.size();

			for (int i = 0; i < llength; i++)
			{
				Sphere *tosphere = lrecord[i];
				record.push_back(tosphere);
			}
			int_near = nn1;
			int_far = ff1;
			return true;
		}
		else if (right_hit)
		{
			if (right->ifleaf)
				rrecord.push_back(right);

			int rlength = rrecord.size();

			for (int i = 0; i < rlength; i++)
			{
				Sphere *tosphere = rrecord[i];
				record.push_back(tosphere);
			}
			int_near = nn2;
			int_far = ff2;
			return true;
		}
*/	
	}
	return false;

}




bool comparex(const Sphere *s1, const Sphere *s2)
{
	if (s1->thisbox.h[0] < s2->thisbox.h[0])
		return true;
	else
		return false;
}

bool comparey(const Sphere *s1, const Sphere *s2)
{
	if (s1->thisbox.h[1] < s2->thisbox.h[1])
		return true;
	else
		return false;
}

bool comparez(const Sphere *s1, const Sphere *s2)
{
	if (s1->thisbox.h[2] < s2->thisbox.h[2])
		return true;
	else
		return false;
}