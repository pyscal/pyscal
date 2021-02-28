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
        int id;
        int loc;
        double posx,posy,posz;
        vector<double> gx();
        void sx(vector<double>);

        py::dict custom;
        int type;
        int condition;

        //mask for calculations to exclude atoms
        bool mask;

        //ghost properties
        int ghost;

  
        //-------------------------------------------------------
        // Neighbor related properties
        //-------------------------------------------------------
        int neighbors[MAXNUMBEROFNEIGHBORS];
        int masks[MAXNUMBEROFNEIGHBORS];
        double neighbordist[MAXNUMBEROFNEIGHBORS];
        double neighborweight[MAXNUMBEROFNEIGHBORS];
        double n_diffx[MAXNUMBEROFNEIGHBORS];
        double n_diffy[MAXNUMBEROFNEIGHBORS];
        double n_diffz[MAXNUMBEROFNEIGHBORS];
        double n_r[MAXNUMBEROFNEIGHBORS];
        double n_phi[MAXNUMBEROFNEIGHBORS];
        double n_theta[MAXNUMBEROFNEIGHBORS];

        vector<datom> temp_neighbors;
        double cutoff;
        int n_neighbors;
        int isneighborset;

        //function to set neighbors
        void sneighdist(vector<double>);
        vector<double> gneighdist();
        void sneighbors(vector<int> nns);
        vector<int> gneighbors();
        void sneighborweights(vector<double> nns);
        vector<double> gneighborweights();
        void sdistvecs(vector<vector<double>>);
        vector<vector<double>> gdistvecs();
        void slocalangles(vector<vector<double>>);
        vector<vector<double>> glocalangles();
        void find_filtered_neighbors(double);
        vector<vector<int>> next_neighbors;
        vector<vector<int>> next_neighbor_masks;
        vector<vector<double>> next_neighbor_distances;
        vector<int> next_neighbor_counts;
        int head;
 
        //-------------------------------------------------------
        // Q parameter properties
        //-------------------------------------------------------
        double sij[MAXNUMBEROFNEIGHBORS];
        double realQ6[13],imgQ6[13];
        void ssij(vector<double>);
        vector<double> gsij();
        double q[11];
        double aq[11];
        double realq[11][25];
        double imgq[11][25];
        double arealq[11][25];
        double aimgq[11][25];
        vector<complex<double>> get_qcomps(int, bool);
        vector<double> gallq();
        vector<double> gallaq();
        void sallq(vector<double>);
        void sallaq(vector<double>);
        double gq(int);
        void sq(int, double);
        double gq_big(int, bool);
        void sq_big( int, double, bool);
        vector<double> gq_big(vector<int>, bool);
        void sq_big(vector<int>, vector<double>, bool);
        double gaq(int);
        void saq(int, double);
        double sii;
        double disorder;
        double avgdisorder;

        //-------------------------------------------------------
        // Solid related properties
        //-------------------------------------------------------
        int frenkelnumber;
        double avq6q6;
        int belongsto;
        bool lcluster;
        bool issurface;
        bool issolid;
        int structure;


        //-------------------------------------------------------
        // Voronoi related properties
        //-------------------------------------------------------
        int facevertices[MAXNUMBEROFNEIGHBORS];
        int faceverticenumbers[MAXNUMBEROFNEIGHBORS];
        double faceperimeters[MAXNUMBEROFNEIGHBORS];
        int n3, n4, n5, n6;
        vector<vector<double>> edgelengths;
        vector<vector<double>> vertex_positions;
        vector<double> vertex_vectors;
        vector<int> vertex_numbers;
        double volume;
        double avgvolume;
        vector<int> gfacevertices();
        void sfacevertices(vector<int>);
        vector<double> gfaceperimeters();
        void sfaceperimeters(vector<double>);
        vector<vector<double>> gedgelengths();
        void sedgelengths(vector<vector<double>>);
        vector<int> gvorovector();
        void svorovector(vector<int>);
        vector<vector<double>> gvertexpositions();
        void svertexpositions(vector<vector<double>>);


        //-------------------------------------------------------
        // Angle related properties
        //-------------------------------------------------------
        double angular;
        double avg_angular;
        vector<int> chiparams;

        //-------------------------------------------------------
        // CNA parameters
        //-------------------------------------------------------
        vector<vector<int>> cna;
        vector<vector<int>> common;
        vector<vector<int>> bonds;
        int nn1[4];

        //-------------------------------------------------------
        // Other order parameters
        //-------------------------------------------------------
        vector<double> sro;
        double centrosymmetry;

        //these are the cutoffs for cna
        double lcutsmall;
        double lcutlarge;
        int lneigh;

        double sigma;
        double rho;
        double rstart;
        double rstop;
        double h;
        double kb;
        double gmr(double);
        double entropy_integrand(double);
        void trapezoid_integration();
        //results
        double entropy;
        double avg_entropy;
        double energy;
        double avg_energy;



};
