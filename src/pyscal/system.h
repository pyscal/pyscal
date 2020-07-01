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
#include "atom.h"

namespace py = pybind11;
using namespace std;

struct cell{
  vector<int> members;
  vector<int> neighbor_cells;
};


class System{

    public:

        //-----------------------------------------------------
        // Constructor, Destructor and Access functions
        //-----------------------------------------------------
        System();
        ~System();
        int nop;

        int baseunit;
        string inputfile;
        int fileread;


        //-----------------------------------------------------
        // Simulation box related methods
        //-----------------------------------------------------
        double rot[3][3];
        double rotinv[3][3];
        int triclinic;
        double boxx, boxy, boxz;
        double boxdims[3][2];
        void assign_triclinic_params(vector<vector<double>>, vector<vector<double>>);
        vector<vector<double>> get_triclinic_params();
        vector<vector<double>> gboxvecs();
        void sbox(vector<vector<double>>);
        vector<vector<double>> gbox();
        vector<double> gboxdims();

        //-----------------------------------------------------
        // Atom related methods
        //-----------------------------------------------------
        vector<Atom> atoms;
        void assign_particles( vector<Atom>);
        void read_particle_file(string);    // TBDep
        void set_atoms( vector<Atom>);
        vector<Atom> get_atoms();
        Atom gatom(int);
        void satom(Atom);


        //----------------------------------------------------
        // Neighbor methods
        //----------------------------------------------------
        void get_all_neighbors_normal();
        void process_neighbor(int, int);
        int get_all_neighbors_sann(double);
        int get_all_neighbors_adaptive(double, int, double);
        void get_all_neighbors_voronoi();
        void reset_all_neighbors();        
        double get_abs_distance(int,int,double&,double&,double&);
        double get_abs_distance(Atom , Atom );
        vector<double> get_distance_vector(Atom , Atom);
        void set_neighbordistance(double);
        vector<double> get_pairdistances();
        //variables for a filter
        int filter;
        void sfilter(int);
        int gfilter();
        double neighbordistance;
        int neighborsfound;
        int usecells;
        void susecells(int);
        int gusecells();
        cell *cells;
        int nx, ny, nz;
        int total_cells;
        int cell_index(int, int, int);
        void set_up_cells();
        vector<int> cell_periodic(int, int, int);
        void get_all_neighbors_cells();
        void get_temp_neighbors_cells();
        void get_temp_neighbors_brute();


        //---------------------------------------------------
        // Methods for q calculation
        //---------------------------------------------------
        void set_reqd_qs(vector<int>);
        void set_reqd_aqs(vector<int>);
        void calculate_q(vector <int>);
        void calculate_aq(vector <int>);
        int *reqdqs;
        int lenqs;
        int *reqdaqs;
        int lenaqs;
        double dfactorial(int ,int );
        void convert_to_spherical_coordinates(double , double , double , double &, double &, double &);
        double PLM(int, int, double);
        void YLM(int , int , double , double , double &, double &);
        void QLM(int ,int ,double ,double ,double &, double & );
        void calculate_complexQLM_6();
        double get_number_from_bond(int,int);
        double get_number_from_bond(Atom ,Atom );
        void calculate_frenkel_numbers();
        vector<double> gqvals(int qq);
        vector<double> gaqvals(int qq);
        vector<int> rq_backup;
        int qsfound;
        //disorder vars
        void calculate_disorder();
        void find_average_disorder();        

        //-----------------------------------------------------
        // Solids and Clustering methods
        //-----------------------------------------------------
        double minfrenkel;
        double threshold;
        double avgthreshold;
        int maxclusterid;
        void find_solid_atoms();
        void find_clusters(double);
        void harvest_cluster(const int, const int);
        void find_clusters_recursive(double);
        int largest_cluster();
        void set_nucsize_parameters(double,double,double);
        void get_largest_cluster_atoms();
        int glargestclusterid();
        void slargestclusterid(int);
        int solidq;
        int gsolidq();
        void ssolidq( int);
        int criteria;
        int gcriteria();
        void scriteria( int);

        //-----------------------------------------------------
        // Voronoi based methods
        //-----------------------------------------------------
        int alpha;
        void salpha(int);
        int galpha();
        void find_average_volume();
        int voronoiused;
        double face_cutoff;
        void set_face_cutoff(double);


};
