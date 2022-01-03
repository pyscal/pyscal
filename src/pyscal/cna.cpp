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

void get_cna_neighbors_cn12(py::dict& atoms,
	double lattice_constant,
	int style){

    double d;
    double diffx,diffy,diffz;
    double tempr,temptheta,tempphi;
    vector<double> diffi, diffj;

    vector<vector<double>> positions = atoms[py::str("positions")].cast<vector<vector<double>>>();
    vector<vector<int>> atom_temp_neighbors = atoms[py::str("atom_temp_neighbors")].cast<vector<vector<int>>>();;
    vector<vector<double>> atom_temp_neighbordist = atoms[py::str("atom_temp_neighbordist")].cast<vector<vector<double>>>();;

    int nop = positions.size();
    vector<vector<int>> neighbors(nop);
    vector<vector<double>> neighbordist(nop);
    vector<vector<double>> neighborweight(nop);
    vector<vector<vector<double>>> diff(nop);
    vector<vector<double>> r(nop);
    vector<vector<double>> phi(nop);
    vector<vector<double>> theta(nop);
    vector<double> cutoff(nop);     

    
    double factor = 0.854;
    int ncount = 12;
    vector<double> cutoff;

    for (int ti=0; ti<nop; ti++){
        atoms[ti].cutoff = factor*lattice_constant;
        for(int i=0 ; i<ncount; i++){
            int tj = atoms[ti].temp_neighbors[i].index;
            process_neighbor(ti, tj);
        }
    }
}

void get_cna_neighbors_cn14(py::dict& atoms,
	double lattice_constant){
    
    double factor = 1.207;
    int ncount = 14;
    vector<double> cutoff;

    for (int ti=0; ti<nop; ti++){
        atoms[ti].cutoff = factor*lattice_constant;
        for(int i=0 ; i<ncount; i++){
            int tj = atoms[ti].temp_neighbors[i].index;
            process_neighbor(ti, tj);
        }
    }
}

void get_acna_neighbors(int style){
    /*
    A new neighbor algorithm that finds a specified number of 
    neighbors for each atom.
    There are two styles available:
    
    Style 1: For FCC like structures (HCP/ICO)
    Style 2: For BCC structure
    */

    double dist;

    //reset neighbors
    reset_main_neighbors();

    if (style == 1){ 
        for (int ti=0; ti<nop; ti++){
            if (atoms[ti].temp_neighbors.size() > 11){
                double ssum = 0;
                for(int i=0 ; i<12; i++){
                    ssum += atoms[ti].temp_neighbors[i].dist;
                }
                //process sum
                atoms[ti].cutoff = 1.207*ssum/12.00;
                //now assign neighbors based on this
                for(int i=0 ; i<12; i++){
                    int tj = atoms[ti].temp_neighbors[i].index;
                    dist = atoms[ti].temp_neighbors[i].dist;
                    //if (dist <= atoms[ti].cutoff)
                    process_neighbor(ti, tj);
                }                                 
            }
        }
    }
    else if (style == 2){
        for (int ti=0; ti<nop; ti++){
            if (atoms[ti].temp_neighbors.size() > 13){
                double ssum = 0;
                for(int i=0 ; i<8; i++){
                    ssum += 1.1547*atoms[ti].temp_neighbors[i].dist;
                }
                for(int i=8 ; i<14; i++){
                    ssum += atoms[ti].temp_neighbors[i].dist;
                }
                atoms[ti].cutoff = 1.207*ssum/14.00;
                //now assign neighbors based on this
                for(int i=0 ; i<14; i++){
                    int tj = atoms[ti].temp_neighbors[i].index;
                    dist = atoms[ti].temp_neighbors[i].dist;
                    //if (dist <= atoms[ti].cutoff)
                    process_neighbor(ti, tj);
                }                                 
            }
        }
    }
}