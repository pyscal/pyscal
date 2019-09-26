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

using namespace std;


const int MAXNUMBEROFNEIGHBORS = 100;
const int NILVALUE = 333333;

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
        Atom(vector<double>, int, int);
        virtual ~Atom();

        //basic atom properties
        double posx,posy,posz;
        int id;
        int loc;

        //neighbor related properties
        int neighbors[MAXNUMBEROFNEIGHBORS];
        double neighbordist[MAXNUMBEROFNEIGHBORS];
        double neighborweight[MAXNUMBEROFNEIGHBORS];
        int facevertices[MAXNUMBEROFNEIGHBORS];
        int faceverticenumbers[MAXNUMBEROFNEIGHBORS];
        double faceperimeters[MAXNUMBEROFNEIGHBORS];
        double n_diffx[MAXNUMBEROFNEIGHBORS];
        double n_diffy[MAXNUMBEROFNEIGHBORS];
        double n_diffz[MAXNUMBEROFNEIGHBORS];
        double n_r[MAXNUMBEROFNEIGHBORS];
        double n_phi[MAXNUMBEROFNEIGHBORS];
        double n_theta[MAXNUMBEROFNEIGHBORS];
        int n_neighbors;

        //vertex vectors
        vector<double> vertex_vectors;
        vector<int> vertex_numbers;
        vector<double> gvertex_vectors();
        void svertex_vectors(vector<double>);
        vector<int> gvertex_numbers();
        void svertex_numbers(vector<int>);



        double realQ6[13],imgQ6[13];

        int frenkelnumber;
        double avq6q6;
        //volume calculated by voronoi tesselation
        double volume;
        double avgvolume;
        double gvolume();
        double gavgvolume();
        void svolume(double);
        void savgvolume(double);
        double gasij();
        void sasij(double);

        int belongsto;
        int lcluster;
        int issurface;
        int issolid;
        int structure;
        int type;
        int condition;


        //indicator which is 1 if neighbors are already provided
        int isneighborset;

        //we need some functions to fetch atom properties
        vector<double> gx();
        void sx(vector<double>);
        //int* gx();
        //double gy();
        //double gz();
        //probably wont work
        //in that case will have to return a vector
        //its probably expensive
        //but that doesnt matter because we wont use it regularl

        //function to set neighbors
        void sneighbors(vector<int> nns);
        vector<int> gneighbors();
        vector<int> gfacevertices();
        void sfacevertices(vector<int>);
        vector<double> gfaceperimeters();
        void sfaceperimeters(vector<double>);

        int gnneighbors();
        void sneighborweights(vector<double> nns);
        vector<double> gneighborweights();
        //vector<double> gx();
        vector<int> gcluster();
        void scluster(vector<int>);

        //variables for storing q2-12
        //invidual variables or arrays - individual ones are easier!
        double q[11];
        double aq[11];
        double realq[11][25];
        double imgq[11][25];
        double arealq[11][25];
        double aimgq[11][25];

        vector<double> gallq();
        vector<double> gallaq();
        void sallq(vector<double>);
        void sallaq(vector<double>);

        double gq(int);
        void sq(int, double);
        vector <vector<double>> gqlm(int);

        double gaq(int);
        int gid();
        void sid(int);
        int gloc();
        void sloc(int);
        int gtype();
        void stype(int);
        void saq(int, double);
        vector <vector<double>> gaqlm(int);

        int gsolid();
        void ssolid(int);
        int gstructure();
        void sstructure(int);

        //for custom values
        vector <double> custom;
        void scustom(vector <double>);
        vector<double> gcustom();

        void scondition(int);
        int gcondition();

        int gfrenkelnumber();
        void sfrenkelnumber(int);

        //vector<vector<vector<double>>> gallqcomps();
        //void sallqcomps(vector<vector<vector<double>>>);

};
