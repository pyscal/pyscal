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
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <complex>

namespace py = pybind11;
using namespace std;


const double PI = 3.141592653589793;
const int MAXNUMBEROFNEIGHBORS = 300;
const int NILVALUE = 333333;

//create a structure for sorting
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


class Atom{
    /*
    Class to hold the details of an atom. This is a list of all
    members of the class.
    Attributes
    ----------
    Basic atom properties
    ---------------------
    posx : float
        x coordinate of atom
    posy : float
        y coordinate of atom
    posz : float
        z coordinate of atom
    id : int
        id of the atom
    loc : int
        location of the atom in the array of all atoms in
        the system.
    Neighbor related properties
    ---------------------------
    n_neighbors : int
        number of neighbors of the atom
    Cluster related properties
    --------------------------
    frenkelnumber : int
        frenkelnumber of the atom.
    issolid : int
        0 or 1. 1 if atom is solid.
    structure : int
        structure of the atom.
    belongsto : int
        id of the cluster to which atom belongs to.
     */
    public:

        //-------------------------------------------------------
        // Constructor, Destructor
        //-------------------------------------------------------
        Atom(vector<double>, int, int);
        virtual ~Atom();

        //-------------------------------------------------------
        // Basic Atom properties
        //-------------------------------------------------------
        vector<double> pos;
        int id;
        int condition;
        int loc;
        int type;
        py::dict custom;
        int ghost;
        bool mask;

  
        //-------------------------------------------------------
        // Neighbor related properties
        //-------------------------------------------------------
        vector<int> neighbors;
        vector<int> masks;
        vector<double> neighbordist;
        vector<double> neighborweight;
        vector<vector<double>> diff;
        vector<double> r;
        vector<double> phi;
        vector<double> theta;

        vector<datom> temp_neighbors;
        double cutoff;
        int n_neighbors;
        int isneighborset;

};
