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
#include "atom.h"

namespace py = pybind11;
using namespace std;

/*-----------------------------------------------------
    Some utility objects
-----------------------------------------------------*/

struct cell{
  vector<int> members;
  vector<int> neighbor_cells;
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

void get_all_neighbors_normal(py::dict&,
	const double&,
	const int&,
	const int&,
	const vector<vector<double>>&,
	const vector<vector<double>>&,
	const vector<double>&);

int cell_index(int, int, int, int, int, int);
vector<int> cell_periodic(int, int, int, int, int, int);

vector<cell> set_up_cells(const vector<vector<double>>&,
    const vector<double>&,
    const double);

void get_all_neighbors_cells(py::dict&,
    const double&,
    const int&,
    const int&, 
    const vector<vector<double>>&, 
    const vector<vector<double>>&,
    const vector<double>&);
