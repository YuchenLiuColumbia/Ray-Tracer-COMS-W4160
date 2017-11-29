//basicmath.h is offered by Prof. Michael Reed.
//I did some small changes.

#ifndef BASICMATH_H_
#define BASICMATH_H_

#include <cassert>
#include <iostream>
#include <cmath>

using std::istream;
using std::ostream;


//
//
//
class myvector {
    
public:
    
    myvector () {}
    
    myvector (double i, double j, double k) {
        ijk_[0] = i; ijk_[1] = j; ijk_[2] = k;
    }
    
    void set (double i, double j, double k) {
        ijk_[0] = i; ijk_[1] = j; ijk_[2] = k;
    }

	double n0() {
		return ijk_[0];
	}

	double n1() {
		return ijk_[1];
	}

	double n2() {
		return ijk_[2];
	}


    double operator[] (const int i) const;
    
    void normalize ();
    
    // return scaled version of this vector:
    myvector operator* (const double scale) const {
        return myvector (scale * ijk_[0],scale * ijk_[1],scale * ijk_[2]); 
    }
    
    myvector invert () const { 
        // call the above scaling operator to invert:
        return operator* (-1.0); 
    }
    
    // allow vector addition via +=
    myvector &operator+= (const myvector &v) {
        ijk_[0] = ijk_[0] + v.ijk_[0];
        ijk_[1] = ijk_[1] + v.ijk_[1];
        ijk_[2] = ijk_[2] + v.ijk_[2];
        
        return *this;
    }
   
    // returns the scalar value that is the dotproduct
    // between this and v:
    double dotProduct (const myvector v) const {
        return (ijk_[0] * v.ijk_[0] + ijk_[1] * v.ijk_[1] + ijk_[2] * v.ijk_[2]);
    }
    
    
    // sets this to be the cross product of vectors va and vb:
    // i.e. va X vb
    void crossProduct (const myvector va, const myvector vb) {
        ijk_[0] = va.ijk_[1] * vb.ijk_[2] - va.ijk_[2] * vb.ijk_[1];
        ijk_[1] = va.ijk_[2] * vb.ijk_[0] - va.ijk_[0] * vb.ijk_[2];
        ijk_[2] = va.ijk_[0] * vb.ijk_[1] - va.ijk_[1] * vb.ijk_[0];
    }
    
    // these two functions are used to allow us to read in
    // and print out (or save to a file). they are used by
    // the parser code.
    friend istream &operator>>(istream &is, myvector &v) {
        return is>>v.ijk_[0]>>v.ijk_[1]>>v.ijk_[2];
    }
    
    friend ostream &operator<<(ostream &os, myvector &v) {
        return os<<"<"<<v.ijk_[0]<<","<<v.ijk_[1]<<","<<v.ijk_[2]<<">";
    }
    
    
private:
    
    double ijk_[3];
};




//
//
//
class mypoint {
    
public:
    
    mypoint () {} // default constructor does nothing, not even inits.
    mypoint (double x, double y, double z) {
        xyz_[0] = x; xyz_[1] = y; xyz_[2] = z;
    }

    void set (double x, double y, double z) {
        xyz_[0] = x; xyz_[1] = y; xyz_[2] = z;
    }

	double p0() {
		return xyz_[0];
	}

	double p1() {
		return xyz_[1];
	}

	double p2() {
		return xyz_[2];
	}

    double operator[] (const int i) const;
    
    // subtracting one point from another returns a vector
    myvector operator- (const mypoint &p) const {
        return myvector (
                         xyz_[0] - p.xyz_[0], 
                         xyz_[1] - p.xyz_[1], 
                         xyz_[2] - p.xyz_[2]);
    }
    
    // these two functions are used to allow us to read in
    // and print out (or save to a file). they are used by
    // the parser code.
    friend  istream &operator>>(istream &is, mypoint &p) {
        return is >> p.xyz_[0] >> p.xyz_[1] >> p.xyz_[2];
    }
    
    friend ostream &operator<<(ostream &os, mypoint &p) {
        return os<<"<"<<p.xyz_[0]<<","<<p.xyz_[1]<<","<<p.xyz_[2]<<">";
    }
    
private:
    
    double xyz_[3];
    
};





inline double myvector::operator[]
  (const int i) const
{
    assert(i >= 0 && i < 3);
    return ijk_[i];
}



inline void myvector::normalize ()
{

    assert (! (ijk_[0] == 0.0 && ijk_[1] == 0.0 && ijk_[2] == 0.0));
    
    double len = ijk_[0] * ijk_[0] + 
                 ijk_[1] * ijk_[1] +
                 ijk_[2] * ijk_[2];
    
    len = sqrt (len);
    
    ijk_[0] /= len;
    ijk_[1] /= len;
    ijk_[2] /= len;
    
}



inline double mypoint::operator[]
(const int i) const
{
    assert(i >= 0 && i < 3);
    return xyz_[i];
}

#endif
