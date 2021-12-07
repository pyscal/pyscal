#include "modsystem.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include "string.h"
#include <chrono>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>,
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <map>
#include <string>
#include <any>


double get_abs_distance(const vector<double>& pos1, 
    const vector<double>& pos2, 
	const int& triclinic, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
	const vector<double>& box,
    double& diffx,
    double& diffy,
    double& diffz){

    double abs, ax, ay, az;
    diffx = pos1[0]-pos2[0];
    diffy = pos1[1]-pos2[1];
    diffz = pos1[2]-pos2[2];


    if (triclinic == 1){

        //convert to the triclinic system
        ax = rotinv[0][0]*diffx + rotinv[0][1]*diffy + rotinv[0][2]*diffz;
        ay = rotinv[1][0]*diffx + rotinv[1][1]*diffy + rotinv[1][2]*diffz;
        az = rotinv[2][0]*diffx + rotinv[2][1]*diffy + rotinv[2][2]*diffz;

        //scale to match the triclinic box size
        diffx = ax*box[0];
        diffy = ay*box[1];
        diffz = az*box[2];

        //now check pbc
        //nearest image
        if (diffx> box[0]/2.0) {diffx-=box[0];};
        if (diffx<-box[0]/2.0) {diffx+=box[0];};
        if (diffy> box[1]/2.0) {diffy-=box[1];};
        if (diffy<-box[1]/2.0) {diffy+=box[1];};
        if (diffz> box[2]/2.0) {diffz-=box[2];};
        if (diffz<-box[2]/2.0) {diffz+=box[2];};

        //now divide by box vals - scale down the size
        diffx = diffx/box[0];
        diffy = diffy/box[1];
        diffz = diffz/box[2];

        //now transform back to normal system
        ax = rot[0][0]*diffx + rot[0][1]*diffy + rot[0][2]*diffz;
        ay = rot[1][0]*diffx + rot[1][1]*diffy + rot[1][2]*diffz;
        az = rot[2][0]*diffx + rot[2][1]*diffy + rot[2][2]*diffz;

        //now assign to diffs and calculate distnace
        diffx = ax;
        diffy = ay;
        diffz = az;

        //finally distance
        abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

    }
    else{
        //nearest image
        if (diffx> box[0]/2.0) {diffx-=box[0];};
        if (diffx<-box[0]/2.0) {diffx+=box[0];};
        if (diffy> box[1]/2.0) {diffy-=box[1];};
        if (diffy<-box[1]/2.0) {diffy+=box[1];};
        if (diffz> box[2]/2.0) {diffz-=box[2];};
        if (diffz<-box[2]/2.0) {diffz+=box[2];};
        abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
    }
    return abs;
}

void reset_all_neighbors(vector<py::dict>& atoms){   
    vector<int> tempint;
    vector<int> tempdouble;
    for(unsigned int ti=0; ti<atoms.size(); ti++){
        atoms[ti][py::str("neighbors")] = tempint;
        atoms[ti][py::str("neighbordist")] = tempdouble;
        atoms[ti][py::str("neighborweight")] = tempdouble;
        atoms[ti][py::str("diff")] = tempdouble;
        atoms[ti][py::str("r")] = tempdouble;
        atoms[ti][py::str("phi")] = tempdouble;
        atoms[ti][py::str("theta")] = tempdouble;
        atoms[ti][py::str("cutoff")] = 0.0;
    }
}

void convert_to_spherical_coordinates(double x, 
    double y, 
    double z, 
    double &r, 
    double &phi, 
    double &theta){
    r = sqrt(x*x+y*y+z*z);
    theta = acos(z/r);
    phi = atan2(y,x);
}

void get_all_neighbors_normal(vector<py::dict>& atoms, 
    const int& triclinic,
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box, 
    const double neighbordistance,
    const int filter){

    double d;
    double diffx, diffy, diffz;
    double r,theta, phi;
    int nop = atoms.size();
    vector<double> pos1, pos2, diffi, diffj;
    int type1, type2;

    for (int ti=0; ti<nop; ti++){
        for (int tj=ti+1; tj<nop; tj++){

            //access positions
            pos1 = atoms[ti][py::str("pos")].cast<vector<double>>();
            pos2 = atoms[tj][py::str("pos")].cast<vector<double>>();
            type1 = atoms[ti][py::str("type")].cast<py::int_>();
            type2 = atoms[tj][py::str("type")].cast<py::int_>();

            d = get_abs_distance(pos1, pos2, triclinic, rot, rotinv, box, diffx, diffy, diffz);
            
            if (d < neighbordistance){
                if ((filter == 1) && (type1 != type2)){
                    continue;
                }
                else if ((filter == 2) && (type1 == type2)){
                    continue;
                }
                
                atoms[ti][py::str("neighbors")].cast<py::list>().append(tj);
                atoms[tj][py::str("neighbors")].cast<py::list>().append(ti);
                
                atoms[ti][py::str("neighbordist")].cast<py::list>().append(d);
                atoms[tj][py::str("neighbordist")].cast<py::list>().append(d);

                atoms[ti][py::str("neighborweight")].cast<py::list>().append(1.00);
                atoms[tj][py::str("neighborweight")].cast<py::list>().append(1.00);

                diffi.clear();
                diffi.emplace_back(diffx);
                diffi.emplace_back(diffy);
                diffi.emplace_back(diffz);

                diffj.clear();
                diffj.emplace_back(-diffx);
                diffj.emplace_back(-diffy);
                diffj.emplace_back(-diffz);

                atoms[ti][py::str("diff")].cast<py::list>().append(diffi);
                atoms[tj][py::str("diff")].cast<py::list>().append(diffj);

                convert_to_spherical_coordinates(diffx, diffy, diffz, r, phi, theta);
                atoms[ti][py::str("r")].cast<py::list>().append(r);
                atoms[ti][py::str("phi")].cast<py::list>().append(phi);
                atoms[ti][py::str("theta")].cast<py::list>().append(theta);

                convert_to_spherical_coordinates(-diffx, -diffy, -diffz, r, phi, theta);
                atoms[tj][py::str("r")].cast<py::list>().append(r);
                atoms[tj][py::str("phi")].cast<py::list>().append(phi);
                atoms[tj][py::str("theta")].cast<py::list>().append(theta);

                atoms[ti][py::str("cutoff")] = neighbordistance;
                atoms[tj][py::str("cutoff")] = neighbordistance;
            }
        }
    }
}
