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

double get_abs_distance(const vector<double>&, const vector<double>&, 
	const int&, const vector<vector<double>>&, const vector<vector<double>>&,
	const vector<double>&);


