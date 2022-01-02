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

double get_number_from_bond(const int lm,
	vector<double>& real_qi,
	vector<double>& imag_qi,
	vector<double>& real_qj,
	vector<double>& imag_qj){

    double sum2ti,sum2tj;
    double realdotproduct,imgdotproduct;
    double connection;

    sum2ti = 0.0;
    sum2tj = 0.0;
    realdotproduct = 0.0;
    imgdotproduct = 0.0;

    for (int mi = 0; mi<2*lm+1 ; mi++){
        sum2ti += real_qi[mi]*real_qi[mi] + imag_qi[mi]*imag_qi[mi];
        sum2tj += real_qj[mi]*real_qj[mi] + imag_qj[mi]*imag_qj[mi];
        realdotproduct += real_qi[mi]*real_qj[mi];
        imgdotproduct  += imag_qi[mi]*imag_qj[mi];
    }

    connection = (realdotproduct+imgdotproduct)/(sqrt(sum2tj)*sqrt(sum2ti));
    return connection;
}

void calculate_bonds(py::dict& atoms,
	const int lm,
	const double threshold,
	const int comparecriteria){
    
    int frenkelcons;
    double scalar;

    string key1, key2;
	key1 = "q"+to_string(lm)+"_real";
	key2 = "q"+to_string(lm)+"_imag";

    vector<vector<double>> q_real = atoms[py::str(key1)].cast<vector<vector<double>>>();
    vector<vector<double>> q_imag = atoms[py::str(key2)].cast<vector<vector<double>>>();
    vector<vector<int>> neighbors = atoms[py::str("neighbors")].cast<vector<vector<int>>>();
    vector<double> bonds;
    int nop = neighbors.size();

    vector<vector<double>> sij(nop);
    vector<double> avg_sij;

    double tempsij;

    for (int ti= 0;ti<nop;ti++){

        frenkelcons = 0;
        tempsij = 0.0;
        
        for (int c = 0; c<neighbors[ti].size(); c++){

            scalar = get_number_from_bond(lm, q_real[ti], q_imag[ti], q_real[neighbors[ti][c]], q_imag[neighbors[ti][c]]);
            sij[ti].emplace_back(scalar);
            
            if (comparecriteria == 0)
                if (scalar > threshold) frenkelcons += 1;
            else
                if (scalar < threshold) frenkelcons += 1;
            
            tempsij += scalar;
        }

        bonds.emplace_back(frenkelcons);
        tempsij = tempsij/double(neighbors[ti].size());
        avg_sij.emplace_back(tempsij);
    }

    atoms[py::str("bonds")] = bonds;
    atoms[py::str("sij")] = sij;
    atoms[py::str("avg_sij")] = avg_sij;
}

