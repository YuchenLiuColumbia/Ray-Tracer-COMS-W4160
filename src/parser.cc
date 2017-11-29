
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "parser.h"
#include "basicmath.h"
#include "sphere.h"
#include "triangle.h"
#include "material.h"
#include "light.h"

using namespace std;

//most codes are used offered by Prof. Micheal Reed.
//I did some small changes.

void Parser::parse(
    const char file[],
    std::vector<Sphere *> &asphere,
	std::vector<Triangle *> &atriangle,
	std::vector<pLight *> &plight,
	std::vector<areaLight *> &arealight,
    aLight &alight,
	std::vector<myvector *> &asoftnormals,
	Camera &cam
	)
{
    
    ifstream in(file);
    char buffer[1025];
    string cmd;
    
    int num_cams = 0; // keep track of how many cameras we read in
	Material thismaterial; //keep track of the materials' using
    
    for (int line=1; in.good(); line++) {
        in.getline(buffer,1024);
        buffer[in.gcount()]=0;
        

        cmd="";
        
        istringstream iss(buffer);
        
        iss >> cmd;
        
        if (cmd[0]=='/' || cmd.empty()) {
            // ignore comments or blank lines
            continue;
        } 
        else if (cmd=="s") {
            // Sphere:
            mypoint pos; 
            double r;
            iss >> pos >> r;
			Sphere *thissphere = new Sphere;
			thissphere->init(pos, r, thismaterial);
			asphere.push_back(thissphere);
        } 
        else if (cmd=="t") {
            // Triangle:
			mypoint p1, p2, p3;
			iss >> p1 >> p2 >> p3;
			Triangle *thistriangle = new Triangle;
			thistriangle->init(p1, p2, p3, thismaterial);
			atriangle.push_back(thistriangle);
        }
        else if (cmd=="p") {
			continue;
            // Plane
        }
        else if (cmd=="c") {
            // Camera:
            num_cams++; // keep track of how many we read in

            mypoint pos; myvector dir; 
            double d,iw,ih; 
            int pw,ph;
            iss >> pos >> dir >> d >> iw >> ih >> pw >> ph;
            
			if (dir[0] == 0 && dir[2] == 0)
			{
				std::cerr << "camera error: camera initial direction unable to build coordinates." << endl;
				continue;
			}
            cam.init (pos,dir,d,iw,ih,pw,ph);
        } 
        else if (cmd=="l") {
			string lightcmd;
			iss >> lightcmd;

			if (lightcmd == "p")
			{
				mypoint lp;
				double r, g, b;
				iss >> lp >> r >> g >> b;
				pLight *thisplight = new pLight;
				thisplight->init(lp, r, g, b);
				plight.push_back(thisplight);
			}
			else if (lightcmd == "s")
			{
				mypoint center;
				myvector dir, u_dir;
				double len;
				double r, g, b;
				iss >> center >> dir >> u_dir >> len >> r >> g >> b;
				areaLight *thisarealight = new areaLight;
				thisarealight->init(center, dir, u_dir, len, r, g, b);
				arealight.push_back(thisarealight);
			}
			else if (lightcmd == "a")
			{
				double r, g, b;
				iss >> r >> g >> b;
				alight.init(r, g, b);
			}
			else continue;
        }
        else if (cmd=="m") {
			double dr, dg, db;
			double sr, sg, sb;
			double ir, ig, ib;
			double r;
			iss >> dr >> dg >> db >> sr >> sg >> sb >> r >> ir >> ig >> ib;
			thismaterial.init(dr, dg, db, sr, sg, sb, r, ir, ig, ib);
         }
		else if (cmd == "w")
		{
			std::vector< int > tris;
			std::vector< float > verts;
			char objfile[200];
			iss >> objfile;

			//we set the asoftnormals as a vector to store normals of all vertices
			//so it's important for us to store them in queue
			//everytime we begin to readin a new obj file
			//we begin coding the vertice number from the last readin one
			//使用一个数字来标记一共读入了多少个顶点
			//每当开启一个新的obj文件，就从上一次的累计顶点数开始往下继续存储
			int current_ver = asoftnormals.size();

			read_wavefront_file(objfile, tris, verts, asoftnormals);

			//separate them into triangles;
			int i;
			int numtri = int(tris.size() / 3);
			for (i = 0; i < numtri; i++)
			{

				mypoint p1, p2, p3;
				int vertical1 = tris[i * 3];
				int vertical2 = tris[i * 3 + 1];
				int vertical3 = tris[i * 3 + 2];
				p1.set(verts[vertical1 * 3], verts[vertical1 * 3 + 1], verts[vertical1 * 3 + 2]);
				p2.set(verts[vertical2 * 3], verts[vertical2 * 3 + 1], verts[vertical2 * 3 + 2]);
				p3.set(verts[vertical3 * 3], verts[vertical3 * 3 + 1], verts[vertical3 * 3 + 2]);
				Triangle *thistriangle = new Triangle;

				int newvert1 = vertical1 + current_ver;
				int newvert2 = vertical2 + current_ver;
				int newvert3 = vertical3 + current_ver;
				myvector thisnormal = thistriangle->objinit(p1, p2, p3, newvert1, newvert2, newvert3, thismaterial);
				atriangle.push_back(thistriangle);


				//as we have calculated the normal of this triangle, we add this normal to all vertices
				myvector *new1 = asoftnormals[newvert1];
				myvector *new2 = asoftnormals[newvert2];
				myvector *new3 = asoftnormals[newvert3];
				new1->set(new1->n0() + thisnormal[0],
						  new1->n1() + thisnormal[1],
						  new1->n2() + thisnormal[2]);
				new2->set(new2->n0() + thisnormal[0],
						  new2->n1() + thisnormal[1],
						  new2->n2() + thisnormal[2]);
				new3->set(new3->n0() + thisnormal[0],
						  new3->n1() + thisnormal[1],
						  new3->n2() + thisnormal[2]);

				if (i == 745)
					continue;

			}
		}
        else {
            std::cerr << "Parser error: invalid command at line " << line << endl;
        }

    }
    
    in.close();

	//normalize all N_i of vertices
	if (asoftnormals.size() > 0)
		for (unsigned int n = 0; n < asoftnormals.size(); n++)
		{
			myvector *thisvector = asoftnormals[n];
			if (thisvector->n0() != 0 || thisvector->n1() != 0 || thisvector->n2() != 0)
				thisvector->normalize();
		}
    
    // make sure we read in 1 camera, no more no less 8).
    if (num_cams != 1) {
        std::cerr << "scene file error: exactly ONE camera must be defined." << endl;
    }
}


void Parser::read_wavefront_file(const char *file,
	std::vector< int > &tris,
	std::vector< float > &verts,
	std::vector<myvector *> &asoftnormals
	)
{
	// clear out the tris and verts vectors:
	tris.clear();
	verts.clear();

	ifstream in(file);
	char buffer[1025];
	string cmd;


	for (int line = 1; in.good(); line++) {
		in.getline(buffer, 1024);
		buffer[in.gcount()] = 0;

		cmd = "";

		istringstream iss(buffer);

		iss >> cmd;

		if (cmd[0] == '#' || cmd.empty()) {
			// ignore comments or blank lines
			continue;
		}
		else if (cmd == "v") {
			// got a vertex:

			// read in the parameters:
			double pa, pb, pc;
			iss >> pa >> pb >> pc;

			verts.push_back(pa);
			verts.push_back(pb);
			verts.push_back(pc);

			//push in a new normal for local vectice
			//begin with 0
			myvector *thisvernormal = new myvector;
			thisvernormal->set(0, 0, 0);
			asoftnormals.push_back(thisvernormal);
		}
		else if (cmd == "f") {
			// got a face (triangle)

			// read in the parameters:
			int i, j, k;
			iss >> i >> j >> k;

			// vertex numbers in OBJ files start with 1, but in C++ array
			// indices start with 0, so we're shifting everything down by
			// 1
			tris.push_back(i - 1);
			tris.push_back(j - 1);
			tris.push_back(k - 1);
		}
		else {
			std::cerr << "Parser error: invalid command at line " << line << std::endl;
		}

	}

	std::cout << "found this many tris, verts: " << tris.size() / 3.0 << "  " << verts.size() / 3.0 << std::endl;
	in.close();

}
