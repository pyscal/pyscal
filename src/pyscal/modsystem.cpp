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


double get_abs_distance(const vector<double>& pos1, const vector<double>& pos2, 
	const int& triclinic, const vector<vector<double>>& rot, const vector<vector<double>>& rotinv,
	const vector<double>& box){

    double abs, ax, ay, az;
    double diffx = pos1[0]-pos2[0];
    double diffy = pos1[1]-pos2[1];
    double diffz = pos1[2]-pos2[2];


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