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
    //documentation is available on the binding module    
    public:
        Atom();
        virtual ~Atom();
        double posx,posy,posz;
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

        int belongsto;
        int lcluster;
        int issurface;
        int issolid;
        int structure;
        int id;
        int loc;
        //indicator which is 1 if neighbors are already provided
        int isneighborset;

        //we need some functions to fetch atom properties
        vector<double> gx();
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
        void sneighborweights(vector<double> nns);
        vector<double> gneighborweights();
        //vector<double> gx();
        vector<int> gcluster(); 

        //variables for storing q2-12
        //invidual variables or arrays - individual ones are easier!
        double q[11];
        double aq[11];
        double realq[11][25];
        double imgq[11][25];
        double arealq[11][25];
        double aimgq[11][25];

        
        double gq(int);
        void sq(int, double);
        vector <vector<double>> gqlm(int);

        double gaq(int);
        int gid();
        void saq(int, double);
        vector <vector<double>> gaqlm(int);

        //for custom values
        vector <double> custom;
        void scustom(vector <double>);
        vector<double> gcustom();

        //for vorocell identification
        vector<int> vorovector;
        vector<int> gvorovector();

};


class System{
  
    public:
        
        double dfactorial(int ,int );
        void convert_to_spherical_coordinates(double , double , double , double &, double &, double &);
        double PLM(int, int, double);
        void YLM(int , int , double , double , double &, double &);
        void QLM(int ,int ,double ,double ,double &, double & );
        void get_all_neighbors();
        void get_all_neighbors(string &);
        void reset_all_neighbors();
        void calculate_complexQLM_6();
        double get_number_from_bond(int,int);
        double get_number_from_bond(Atom ,Atom );
        void calculate_frenkel_numbers();
        double get_abs_distance(int,int,double&,double&,double&);
        double get_abs_distance(Atom , Atom );
        System();
        virtual ~System();

        Atom* atoms;
    
        void read_particle_file();
        void read_particle_instance(int,int);
        int calculate_nucsize();	//variant of function above
        int cluster_criteria(int,int );
        void find_solids();
        void find_clusters();
        int largest_cluster();
        void set_nucsize_parameters(int,double,double);
        void set_inputfile(string);
        void set_neighbordistance(double);
        void assign_particles( vector<Atom>, vector<double>);
        void get_largest_cluster_atoms();
        //functions to set the list of reqd qs
        //again, error checking would be amazing here.
        void set_reqd_qs(vector<int>);
        void set_reqd_aqs(vector<int>);
        void calculate_q();
        void calculate_aq();
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
        double neighbordistance;
        double threshold;
        double avgthreshold;
        int maxclusterid;

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
        vector<double> gboxdims();
        
        //system flags
        int neighborsfound;
        int qsfound;
        int fileread;

};