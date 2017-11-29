
#include "triangle.h"

using namespace std;

bool comparex(const Triangle *s1, const Triangle *s2);
bool comparey(const Triangle *s1, const Triangle *s2);
bool comparez(const Triangle *s1, const Triangle *s2);

void Triangle::init(mypoint p1_, mypoint p2_, mypoint p3_, Material thismaterial)
	{
		p1 = p1_;
		p2 = p2_;
		p3 = p3_;
		material = thismaterial;
		ifleaf = true;
		isobjtri = false;

		thisbox.buildbox(p1_, p2_, p3_);

		return;
	}

myvector Triangle::objinit(mypoint p1_, mypoint p2_, mypoint p3_, int vert1_, int vert2_, int vert3_, Material thismaterial)
{
	p1 = p1_;
	p2 = p2_;
	p3 = p3_;
	vert1 = vert1_;
	vert2 = vert2_;
	vert3 = vert3_;
	material = thismaterial;
	ifleaf = true;
	isobjtri = true;

	thisbox.buildbox(p1_, p2_, p3_);

	myvector ab = p2 - p1;
	myvector ac = p3 - p1;
	myvector normal;
	normal.crossProduct(ab, ac);
	normal.normalize();

	return normal;
}


bool Triangle::intersect(double &t, double &tmin, ray &thisray, mypoint &thispoint, myvector &normal_out)
{

	myvector dd = thisray.dir;
	mypoint pp = thisray.pos;
	double ttemp, gamma, beta;
	
	double a = p1[0] - p2[0];
	double b = p1[1] - p2[1];
	double c = p1[2] - p2[2];

	double d = p1[0] - p3[0];
	double e = p1[1] - p3[1];
	double f = p1[2] - p3[2];

	myvector ab = p2 - p1;
	myvector ac = p3 - p1;
	myvector normal;
	normal.crossProduct(ab, ac);
	normal.normalize();
	myvector normal_inv;
	normal_inv.crossProduct(ac, ab);
	normal.normalize();
	double g = dd[0];
	double h = dd[1];
	double i = dd[2];
	double j = p1[0] - pp[0];
	double k = p1[1] - pp[1];
	double l = p1[2] - pp[2];
	//next we calculate B1 = ei-hf, B2 = gf-di, 
	//B3 = dh-eg, G1 = ak-jb, G2 = jc-al, G3 = bl-kc.
	double B1 = e * i - h * f;
	double B2 = g * f - d * i;
	double B3 = d * h - e * g;
	double G1 = a * k - j * b;
	double G2 = j * c - a * l;
	double G3 = b * l - k * c;

	//calculate M
	double M = a * B1 + b * B2 + c * B3;

	//calculate t
	ttemp = -(f * G1 + e * G2 + d * G3) / M;
	if (ttemp < 0 || ttemp > tmin)
		return false;

	//calculate gamma
	gamma = (i * G1 + h * G2 + g * G3) / M;
	if (gamma < 0 || gamma > 1)
		return false;

	//calculate beta
	beta = (j * B1 + k * B2 + l * B3) / M;
	if (beta < 0 || beta + gamma > 1)
		return false;

	t = ttemp;
	thispoint.set(pp[0] + dd[0] * t, pp[1] + dd[1] * t, pp[2] + dd[2] * t);


		normal_out = normal;
	
	return true;
}

bool Triangle::intersect(double &t, double &tmin, ray &thisray, mypoint &thispoint, myvector &normal_out, std::vector<myvector *> &asoftnormals)
{

	myvector dd = thisray.dir;
	mypoint pp = thisray.pos;
	double ttemp, gamma, beta;

	double a = p1[0] - p2[0];
	double b = p1[1] - p2[1];
	double c = p1[2] - p2[2];

	double d = p1[0] - p3[0];
	double e = p1[1] - p3[1];
	double f = p1[2] - p3[2];

	myvector normal;
	myvector normal_inv;

	myvector ab = p2 - p1;
	myvector ac = p3 - p1;
	normal.crossProduct(ab, ac);
	normal.normalize();
	normal_inv.crossProduct(ac, ab);
	normal.normalize();

	double g = dd[0];
	double h = dd[1];
	double i = dd[2];
	double j = p1[0] - pp[0];
	double k = p1[1] - pp[1];
	double l = p1[2] - pp[2];
	//next we calculate B1 = ei-hf, B2 = gf-di, 
	//B3 = dh-eg, G1 = ak-jb, G2 = jc-al, G3 = bl-kc.
	double B1 = e * i - h * f;
	double B2 = g * f - d * i;
	double B3 = d * h - e * g;
	double G1 = a * k - j * b;
	double G2 = j * c - a * l;
	double G3 = b * l - k * c;

	//calculate M
	double M = a * B1 + b * B2 + c * B3;

	//calculate t
	ttemp = -(f * G1 + e * G2 + d * G3) / M;
	if (ttemp < 0 || ttemp > tmin)
		return false;

	//calculate gamma
	gamma = (i * G1 + h * G2 + g * G3) / M;
	if (gamma < 0 || gamma > 1)
		return false;

	//calculate beta
	beta = (j * B1 + k * B2 + l * B3) / M;
	if (beta < 0 || beta + gamma > 1)
		return false;

	t = ttemp;
	thispoint.set(pp[0] + dd[0] * t, pp[1] + dd[1] * t, pp[2] + dd[2] * t);

	if (isobjtri)
	{
		myvector *N_1 = asoftnormals[vert1];
		myvector *N_2 = asoftnormals[vert2];
		myvector *N_3 = asoftnormals[vert3];

		myvector normaltemp;
		myvector AB, AC, PB, PC, PA;
		AB.set(p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]);
		AC.set(p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]);
		PB.set(p2[0] - thispoint[0], p2[1] - thispoint[1], p2[2] - thispoint[2]);
		PC.set(p3[0] - thispoint[0], p3[1] - thispoint[1], p3[2] - thispoint[2]);
		PA.set(p1[0] - thispoint[0], p1[1] - thispoint[1], p1[2] - thispoint[2]);

		//square ABC
		normaltemp.crossProduct(AB, AC);
		double ABC = normal.dotProduct(normaltemp);

		//square PBC
		normaltemp.crossProduct(PB, PC);
		double PBC = normal.dotProduct(normaltemp);

		//square PCA
		normaltemp.crossProduct(PC, PA);
		double PCA = normal.dotProduct(normaltemp);

		double ALPHA = PBC / ABC;
		double BETA = PCA / ABC;
		double GAMMA = 1 - ALPHA - BETA;

		normal.set(ALPHA * N_1->n0() + BETA * N_2->n0() + GAMMA * N_3->n0(),
		       	   ALPHA * N_1->n1() + BETA * N_2->n1() + GAMMA * N_3->n1(),
			       ALPHA * N_1->n2() + BETA * N_2->n2() + GAMMA * N_3->n2());

		normal.normalize();

	}


	normal_out = normal;

	return true;
}


void Triangle::createtree(std::vector<Triangle *> &atriangle, int axis)
{
	int atrilength = atriangle.size();
	left = new Triangle;
	right = new Triangle;

	if (atrilength == 1)
	{
		left = atriangle[0];
		right = nullptr;
		thisbox = left->thisbox;
	}
	else if (atrilength == 2)
	{
		left = atriangle[0];
		right = atriangle[1];
		thisbox.combinebox(left->thisbox, right->thisbox);
	}
	else
	{
		//sort the sphere class by axis
		if (axis == 0)
			std::nth_element(atriangle.begin(), atriangle.begin() + atriangle.size() / 2, atriangle.end(), comparex);
		else if (axis == 1)
			std::nth_element(atriangle.begin(), atriangle.begin() + atriangle.size() / 2, atriangle.end(), comparey);
		else
			std::nth_element(atriangle.begin(), atriangle.begin() + atriangle.size() / 2, atriangle.end(), comparez);

		vector <Triangle *> lefttriangle;
		vector <Triangle *> righttriangle;
		int mid = atrilength / 2;
		for (int i = 0; i < mid; i++)
		{
			Triangle *thistriangle = atriangle[i];
			lefttriangle.push_back(thistriangle);
		}
		for (int i = mid; i < atrilength; i++)
		{
			Triangle *thistriangle = atriangle[i];
			righttriangle.push_back(thistriangle);
		}

		left->createtree(lefttriangle, (axis + 1) % 3);
		right->createtree(righttriangle, (axis + 1) % 3);
		thisbox.combinebox(left->thisbox, right->thisbox);

	}
}

bool Triangle::hit(ray &thisray, double &int_near, double &int_far, std::vector<Triangle* > &record)
{
	bool left_hit = false, right_hit = false;

	// in this step, if it's a hit, then the int_near would be changed, fitting to the bounding box
	// if it's not hit, then the value would be change
	// if it hits but the final answer is 'not hit' (targeting gap), 
	//the return value would be false and doesn't influence the result
	if (thisbox.hitbox(thisray, int_near, int_far))
	{
		//measure if it's from inside the box
		// if it is, then it should be a reflected ray or shadow ray, then it definitely should ignore itself
		if ((ifleaf) && left == nullptr && right == nullptr)
		{
			if (thisbox.insidebox(thisray))
				return false;
			else
				return true;
		}

		//		vector <Triangle*> lrecord;
		//		vector <Triangle*> rrecord;
		double nn1 = int_near, nn2 = int_near;
		double ff1 = int_far, ff2 = int_far;
		left_hit = ((left != nullptr) && (left->hit(thisray, nn1, ff1, record)));
		if (left_hit)
			record.push_back(left);
		//		if (left_hit)
		//			right_hit = ((right != nullptr) && (right->hit(thisray, nn2, nn1, rrecord)));
		//		else
		right_hit = ((right != nullptr) && (right->hit(thisray, nn2, ff2, record)));
		if (right_hit)
			record.push_back(right);

		/*		if (left_hit && right_hit)
				{
				if (left->ifleaf)
				record.push_back(left);
				if (right->ifleaf)
				rrecord.push_back(right);

				int llength = lrecord.size();
				int rlength = rrecord.size();

				for (int i = 0; i < llength; i++)
				{
				Triangle *totri = lrecord[i];
				record.push_back(totri);
				}
				for (int i = 0; i < rlength; i++)
				{
				Triangle *totri = rrecord[i];
				record.push_back(totri);
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
				Triangle *totri = lrecord[i];
				record.push_back(totri);
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
				Triangle *totri = rrecord[i];
				record.push_back(totri);
				}
				int_near = nn2;
				int_far = ff2;
				return true;
				}
		*/

	}
	return false;

}

bool comparex(const Triangle *s1, const Triangle *s2)
{
	if (s1->thisbox.h[0] < s2->thisbox.h[0])
		return true;
	else
		return false;
}

bool comparey(const Triangle *s1, const Triangle *s2)
{
	if (s1->thisbox.h[1] < s2->thisbox.h[1])
		return true;
	else
		return false;
}

bool comparez(const Triangle *s1, const Triangle *s2)
{
	if (s1->thisbox.h[2] < s2->thisbox.h[2])
		return true;
	else
		return false;
}