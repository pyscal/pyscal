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
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <map>
#include <string>
#include <any>
#include "voro++.hh"

using namespace voro;

void get_all_neighbors_voronoi(py::dict& atoms,
    const double neighbordistance,
    const int triclinic,
    const vector<vector<double>> rot, 
    const vector<vector<double>> rotinv,
    const vector<double> box)
    {

    double d;
    double diffx,diffy,diffz;
    double tempr,temptheta,tempphi;
    vector<double> diffi, diffj;

    vector<vector<double>> positions = atoms[py::str("positions")].cast<vector<vector<double>>>();
    vector<bool> mask_1 = atoms[py::str("mask_1")].cast<vector<bool>>();
    vector<bool> mask_2 = atoms[py::str("mask_2")].cast<vector<bool>>();

    int nop = positions.size();
    vector<vector<int>> neighbors(nop);
    vector<vector<double>> neighbordist(nop);
    vector<vector<double>> neighborweight(nop);
    vector<vector<vector<double>>> diff(nop);
    vector<vector<double>> r(nop);
    vector<vector<double>> phi(nop);
    vector<vector<double>> theta(nop);
    vector<double> cutoff(nop); 

    pre_container pcon(0.00, box[0], 0.00, box[1], 0.0, box[2], true, true, true);
    for(int i=0; i<nop; i++){
        pos = atoms[i].gx();
        pos = remap_atom(pos);
        pcon.put(i, pos[0], pos[1], pos[2]);
    }
    pcon.guess_optimal(tnx,tny,tnz);
    //container con(boxdims[0][0],boxdims[1][1],boxdims[1][0],boxdims[1][1],boxdims[2][0],boxdims[2][1],tnx,tny,tnz,true,true,true, nop);
    container con(0.00, boxx, 0.00, boxy, 0.0, boxz, tnx, tny, tnz, true, true, true, nop);
    pcon.setup(con);

} 	