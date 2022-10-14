#include "atom.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include "string.h"
#include <chrono>
#include <pybind11/stl.h>
#include <complex>

//-------------------------------------------------------
// Constructor, Destructor
//-------------------------------------------------------
Atom::Atom( vector<double> poss, int idd, int typ){

    pos = poss;
    id = idd;
    type = typ;
    ghost = 0;

}

Atom::~Atom(){ }
