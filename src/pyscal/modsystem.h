#include <iostream>
#include <iostream>
#include <exception>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <vector>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <string>
#include <any>

const double PI = 3.141592653589793;

namespace py = pybind11;
using namespace std;

/*-----------------------------------------------------
    Some utility objects
-----------------------------------------------------*/

struct cell{
  vector<int> members;
  vector<int> neighbor_cells;
};

struct datom{
    double dist;
    int  index;
};

//create another for the sorting algorithm
struct by_dist{
    bool operator()(datom const &datom1, datom const &datom2){
        return (datom1.dist < datom2.dist);
    }
};
/*-----------------------------------------------------
    Neighbor Methods
-----------------------------------------------------*/

double get_abs_distance(vector<double>, vector<double>,
    const int&, 
    const vector<vector<double>>&, 
    const vector<vector<double>>&, 
    const vector<double>&, 
    double&, double&, double&);

void reset_all_neighbors(py::dict&);

void convert_to_spherical_coordinates(double, double, double, 
    double&, double&, double&);

void get_all_neighbors_normal(py::dict& atoms,
    const double neighbordistance,
    const int triclinic,
    const vector<vector<double>> rot,
    const vector<vector<double>> rotinv,
    const vector<double> box);

int cell_index(int, int, int, int, int, int);
vector<int> cell_periodic(int, int, int, int, int, int);

vector<cell> set_up_cells(const vector<vector<double>>&,
    const vector<double>&,
    const double);

void get_all_neighbors_cells(py::dict&,
    const double&,
    const int&,
    const vector<vector<double>>&, 
    const vector<vector<double>>&,
    const vector<double>&);

void get_temp_neighbors_brute(const vector<vector<double>>& positions,
    const vector<bool>& mask_1,
    const vector<bool>& mask_2,
    vector<vector<datom>>& temp_neighbors,
    const int& triclinic,
    const double neighbordistance, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box);

void get_temp_neighbors_cells(const vector<vector<double>>& positions,
    const vector<bool>& mask_1,
    const vector<bool>& mask_2,
    vector<vector<datom>>& temp_neighbors,
    const int& triclinic,
    const double neighbordistance, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box);

int get_all_neighbors_bynumber(py::dict& atoms,
    double& neighbordistance,
    const int& triclinic, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box,
    double prefactor,
    int nns, 
    int usecells,
    int assign);

int get_all_neighbors_sann(py::dict& atoms,
    double& neighbordistance,
    const int& triclinic,
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box,
    double prefactor,
    int usecells);

int get_all_neighbors_adaptive(py::dict& atoms,
    double& neighbordistance,
    const int& triclinic,
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box,
    double prefactor,
    int nlimit,
    double padding, 
    int usecells);


/*-----------------------------------------------------
    Steinhardt Methods
-----------------------------------------------------*/

void calculate_factors(const int lm, 
    vector<vector<double>>& alm,
    vector<vector<double>>& blm, 
    vector<vector<double>>& clm,
    vector<double>& dl,
    vector<double>& el);

vector<vector<double>> calculate_plm(const int lm,
    const double costheta, 
    const double sintheta);

double dfactorial(int l,
    int m);

vector<vector<vector<double>>> calculate_ylm(const int lm,
    const double costheta,
    const double sintheta,
    const double cosphi,
    const double sinphi);

vector<vector<vector<vector<double>>>> calculate_q_atom(const int lm,
    const vector<double>& theta,
    const vector<double>& phi);

void calculate_q(py::dict& atoms,
    const int lm);

double plm(const int l,
    const int m,
    const double theta);

double sph_legendre(const int l,
    const int m,
    const double theta);

void calculate_qlm(const int l, 
    const int m, 
    const double theta, 
    const double phi, 
    double &ylm_real, 
    double &ylm_imag);

void calculate_q_single(py::dict& atoms,
    const int lm);

void calculate_aq_single(py::dict& atoms,
    const int lm);

void calculate_disorder(py::dict& atoms,
    const int lm);

void calculate_bonds(py::dict& atoms,
    const int lm,
    const double threshold,
    const double avgthreshold,
    const double minbonds,
    const int comparecriteria,
    const int criteria);

void extract_cluster(int ti,
    int clusterindex,
    vector<bool>& condition,
    vector<bool>& ghost,
    vector<vector<int>>& neighbors,
    vector<vector<double>>& neighbordist,
    vector<double>& cutoff,
    vector<int>& cluster);

void find_clusters(py::dict& atoms,
    double clustercutoff);