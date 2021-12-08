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

double get_abs_distance(vector<double> pos1, vector<double> pos2, 
	const int& triclinic, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
	const vector<double>& box,
    double& diffx,
    double& diffy,
    double& diffz){

    double abs, ax, ay, az;
    diffx = pos1[0] - pos2[0];
    diffy = pos1[1] - pos2[1];
    diffz = pos1[2] - pos2[2];


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

void reset_all_neighbors(py::dict& atoms){   
    vector<vector<int>> temp2int;
    vector<double> tempdouble;
    vector<vector<double>> temp2double;
    vector<vector<vector<double>>> temp3double;

    for(unsigned int ti=0; ti<atoms.size(); ti++){
        atoms[py::str("neighbors")] = temp2int;
        atoms[py::str("neighbordist")] = temp2double;
        atoms[py::str("neighborweight")] = temp2double;
        atoms[py::str("diff")] = temp3double;
        atoms[py::str("r")] = temp2double;
        atoms[py::str("phi")] = temp2double;
        atoms[py::str("theta")] = temp2double;
        atoms[py::str("cutoff")] = tempdouble;
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

void perturb_atom(vector<Atom>& atoms){
    for(int i=0; i<atoms.size(); i++)
        atoms[i].neighbors.emplace_back(1);
}

void get_all_neighbors_normal(py::dict& atoms,
    const double& neighbordistance,
    const int& triclinic,
    const int& filter, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box)
    {
    
    double d;
    double diffx,diffy,diffz;
    double tempr,temptheta,tempphi;
    vector<double> diffi, diffj;

    //access positions and put it in an array
    vector<vector<double>> positions = atoms[py::str("positions")].cast<vector<vector<double>>>();

    int nop = positions.size();
    vector<vector<int>> neighbors(nop);
    vector<vector<double>> neighbordist(nop);
    vector<vector<double>> neighborweight(nop);
    vector<vector<vector<double>>> diff(nop);
    vector<vector<double>> r(nop);
    vector<vector<double>> phi(nop);
    vector<vector<double>> theta(nop);
    vector<double> cutoff(nop); 

    //now loop and calculate things
    for (int ti=0; ti<nop; ti++){
        for (int tj=ti+1; tj<nop; tj++){
            d = get_abs_distance(positions[ti], positions[tj],
                triclinic, rot, rotinv, box, 
                diffx, diffy, diffz);
            if (d < neighbordistance){
                //if ((filter == 1) && (atoms[ti].type != atoms[tj].type)){
                //    continue;
                //}
                //else if ((filter == 2) && (atoms[ti].type == atoms[tj].type)){
                 //   continue;
                //}
                
                neighbors[ti].emplace_back(tj);
                neighbors[tj].emplace_back(ti);

                neighbordist[ti].emplace_back(d);
                neighbordist[tj].emplace_back(d);

                neighborweight[ti].emplace_back(1.00);
                neighborweight[tj].emplace_back(1.00);

                diffi.clear();
                diffi.emplace_back(diffx);
                diffi.emplace_back(diffy);
                diffi.emplace_back(diffz);

                diffj.clear();
                diffj.emplace_back(-diffx);
                diffj.emplace_back(-diffy);
                diffj.emplace_back(-diffz);

                diff[ti].emplace_back(diffi);
                diff[tj].emplace_back(diffj);
                
                convert_to_spherical_coordinates(diffx, diffy, diffz, tempr, tempphi, temptheta);
                
                r[ti].emplace_back(tempr);
                phi[ti].emplace_back(tempphi);
                theta[ti].emplace_back(temptheta);
                //n_neighbors += 1;
                cutoff[ti] = neighbordistance;

                convert_to_spherical_coordinates(-diffx, -diffy, -diffz, tempr, tempphi, temptheta);

                r[tj].emplace_back(tempr);
                phi[tj].emplace_back(tempphi);
                theta[tj].emplace_back(temptheta);
                //atoms[tj].n_neighbors += 1;
                cutoff[tj] = neighbordistance;

            }
        }
    }

    //calculation over lets assign
    atoms[py::str("neighbors")] = neighbors;
    atoms[py::str("neighbordist")] = neighbordist;
    atoms[py::str("neighborweight")] = neighborweight;
    atoms[py::str("diff")] = diff;
    atoms[py::str("r")] = r;
    atoms[py::str("theta")] = theta;
    atoms[py::str("phi")] = phi;
    atoms[py::str("cutoff")] = cutoff;
}

/*
void get_all_neighbors_normal(vector<Atom>& atoms,
    const double& neighbordistance,
    const int& triclinic,
    const int& filter, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box)
    {
    
    double d;
    double diffx,diffy,diffz;
    double r,theta,phi;
    vector<double> diffi, diffj;

    int nop = atoms.size();

    //clear all arrays
    for (int ti=0; ti<nop; ti++){
        atoms[ti].neighbors.clear();
        atoms[ti].neighbordist.clear();
        atoms[ti].neighborweight.clear();
        atoms[ti].diff.clear();
        atoms[ti].r.clear();
        atoms[ti].phi.clear();
        atoms[ti].theta.clear();
    }    

    for (int ti=0; ti<nop; ti++){
        for (int tj=ti+1; tj<nop; tj++){
            d = get_abs_distance(atoms[ti], atoms[tj],
                triclinic, rot, rotinv, box, 
                diffx, diffy, diffz);
            if (d < neighbordistance){
                if ((filter == 1) && (atoms[ti].type != atoms[tj].type)){
                    continue;
                }
                else if ((filter == 2) && (atoms[ti].type == atoms[tj].type)){
                    continue;
                }
                
                atoms[ti].neighbors.emplace_back(tj);
                atoms[tj].neighbors.emplace_back(ti);

                atoms[ti].neighbordist.emplace_back(d);
                atoms[tj].neighbordist.emplace_back(d);

                atoms[ti].neighborweight.emplace_back(1.00);
                atoms[tj].neighborweight.emplace_back(1.00);

                diffi.clear();
                diffi.emplace_back(diffx);
                diffi.emplace_back(diffy);
                diffi.emplace_back(diffz);

                diffj.clear();
                diffj.emplace_back(-diffx);
                diffj.emplace_back(-diffy);
                diffj.emplace_back(-diffz);

                atoms[ti].diff.emplace_back(diffi);
                atoms[tj].diff.emplace_back(diffj);
                
                convert_to_spherical_coordinates(diffx, diffy, diffz, r, phi, theta);
                
                atoms[ti].r.emplace_back(r);
                atoms[ti].phi.emplace_back(phi);
                atoms[ti].theta.emplace_back(theta);
                atoms[ti].n_neighbors += 1;
                atoms[ti].cutoff = neighbordistance;

                convert_to_spherical_coordinates(-diffx, -diffy, -diffz, r, phi, theta);

                atoms[tj].r.emplace_back(r);
                atoms[tj].phi.emplace_back(phi);
                atoms[tj].theta.emplace_back(theta);
                atoms[tj].n_neighbors += 1;
                atoms[tj].cutoff = neighbordistance;

            }
        }
    }
}

*/