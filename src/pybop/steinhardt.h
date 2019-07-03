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


using namespace std;


const int MAXNUMBEROFNEIGHBORS = 100;
const double PI = 3.141592653589793;
const int NILVALUE = 33333333;

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
        Atom();
        virtual ~Atom();
        
        //basic atom properties
        double posx,posy,posz;
        int id;
        int loc;

        //neighbor related properties
        int neighbors[MAXNUMBEROFNEIGHBORS];
        double neighbordist[MAXNUMBEROFNEIGHBORS];
        double neighborweight[MAXNUMBEROFNEIGHBORS];
        double n_diffx[MAXNUMBEROFNEIGHBORS];
        double n_diffy[MAXNUMBEROFNEIGHBORS];
        double n_diffz[MAXNUMBEROFNEIGHBORS];
        double n_r[MAXNUMBEROFNEIGHBORS];
        double n_phi[MAXNUMBEROFNEIGHBORS];
        double n_theta[MAXNUMBEROFNEIGHBORS];
        int n_neighbors;


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

        int belongsto;
        int lcluster;
        int issurface;
        int issolid;
        int structure;
        int type;
        
        
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

        //for custom values
        vector <double> custom;
        void scustom(vector <double>);
        vector<double> gcustom();

        //for vorocell identification
        int vorovector[4];
        vector<int> gvorovector();
        void svorovector(vector<int>);

};


class System{
  
    public:
        
        double dfactorial(int ,int );
        void convert_to_spherical_coordinates(double , double , double , double &, double &, double &);
        double PLM(int, int, double);
        void YLM(int , int , double , double , double &, double &);
        void QLM(int ,int ,double ,double ,double &, double & );
        void get_all_neighbors_normal();
        void get_all_neighbors_adaptive(int, double);
        void get_all_neighbors_voronoi();
        void reset_all_neighbors();
        void calculate_complexQLM_6();
        double get_number_from_bond(int,int);
        double get_number_from_bond(Atom ,Atom );
        void calculate_frenkel_numbers();
        double get_abs_distance(int,int,double&,double&,double&);
        double get_abs_distance(Atom , Atom );
        System();
        ~System();

        //Atom* atoms;
        vector<Atom> atoms;
    
        void read_particle_file(string);
        //void read_particle_instance(int,int);
        int calculate_nucsize();	//variant of function above
        int cluster_criteria(int,int );
        void find_solids();
        void find_clusters();
        void harvest_cluster(const int, const int);
        void find_clusters_recursive();
        int largest_cluster();
        void set_nucsize_parameters(double,int,double,double);
        //void set_inputfile(string);
        void set_neighbordistance(double);
        void assign_particles( vector<Atom>, vector<vector<double>>);
        void assign_triclinic_params(vector<vector<double>>, vector<vector<double>>);
        void get_largest_cluster_atoms();

        //average volumes
        void find_average_volume();
        //functions to set the list of reqd qs
        //again, error checking would be amazing here.
        void set_reqd_qs(vector<int>);
        void set_reqd_aqs(vector<int>);
        void calculate_q(vector <int>);
        void calculate_aq(vector <int>);
        int *reqdqs;
        int lenqs;
        int *reqdaqs;
        int lenaqs;
        //old params
        int nop;
        int baseunit;
        int minfrenkel;
        double boxx, boxy, boxz;
        //array for box
        double boxdims[3][2];
        string inputfile;
        vector<double> get_pairdistances();
        
        double neighbordistance;

        double threshold;
        double avgthreshold;
        int maxclusterid;

        double rot[3][3];
        double rotinv[3][3];
        int triclinic;

        //variables for a filter
        int filter;
        void sfilter(int);
        //int apply_filter(int, int);

        //some access functions for system
        Atom gatom(int);
        void satom(Atom);
        int glargestclusterid();
        int gnop();
        vector<Atom> gallatoms();
        vector<double> gqvals(int qq);
        vector<double> gaqvals(int qq);
        vector<int> rq_backup;
        vector<double> gbox();
        vector<vector<double>> gboxvecs();
        void sbox(vector<vector<double>>);
        vector<double> gboxdims();
        
        //system flags
        int neighborsfound;
        int qsfound;
        int fileread;

};