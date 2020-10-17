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

        //-----------------------------------------------------
        // Simulation box related methods
        //-----------------------------------------------------
        double rot[3][3];
        double rotinv[3][3];
        int triclinic;
        double boxx, boxy, boxz;
        double boxdims[3][2];
        double box[3][3];
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
        int filter;
        int usecells;
        int nx, ny, nz;
        int total_cells;
        cell *cells;
        double neighbordistance;
        void get_all_neighbors_normal();
        void process_neighbor(int, int);
        int get_all_neighbors_sann(double);
        int get_all_neighbors_bynumber(double, int, int);
        int get_neighbors_from_temp(int);
        int get_all_neighbors_adaptive(double, int, double);
        void get_all_neighbors_voronoi();
        void reset_all_neighbors();
        void reset_main_neighbors();        
        double get_abs_distance(int,int,double&,double&,double&);
        double get_abs_distance(Atom , Atom );
        vector<double> get_distance_vector(Atom , Atom);
        void set_neighbordistance(double);
        vector<double> get_pairdistances();
        //variables for a filter
        void susecells(int);
        int gusecells();
        int cell_index(int, int, int);
        void set_up_cells();
        vector<int> cell_periodic(int, int, int);
        void get_all_neighbors_cells();
        void get_temp_neighbors_cells();
        void get_temp_neighbors_brute();
        void store_neighbor_info();
        void get_diamond_neighbors();
        void set_atom_cutoff(double);

        //---------------------------------------------------
        // Methods for q calculation
        //---------------------------------------------------
        int *reqdqs;
        int lenqs;
        int *reqdaqs;
        int lenaqs;
        vector<double> gqvals(int qq);
        vector<double> gaqvals(int qq);
        vector<int> rq_backup;
        void set_reqd_qs(vector<int>);
        void set_reqd_aqs(vector<int>);
        void calculate_q(vector <int>);
        void calculate_aq(vector <int>);
        double dfactorial(int ,int );
        void convert_to_spherical_coordinates(double , double , double , double &, double &, double &);
        double PLM(int, int, double);
        void YLM(int , int , double , double , double &, double &);
        void QLM(int ,int ,double ,double ,double &, double & );
        void calculate_complexQLM_6();
        double get_number_from_bond(int,int);
        double get_number_from_bond(Atom ,Atom );
        void calculate_frenkel_numbers();
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
        int solidq;
        int criteria;
        int comparecriteria;
        void find_solid_atoms();
        void find_clusters(double);
        void harvest_cluster(const int, const int);
        void find_clusters_recursive(double);
        int largest_cluster();
        void set_nucsize_parameters(double,double,double);
        void get_largest_cluster_atoms();

        //-----------------------------------------------------
        // Voronoi based methods
        //-----------------------------------------------------
        int alpha;
        void find_average_volume();
        int voronoiused;
        double face_cutoff;
        void set_face_cutoff(double);

        //-------------------------------------------------------
        // Other order parameters
        //-------------------------------------------------------
        vector<int> calculate_acna();
        vector<int> calculate_cna();
        double switching_fn(double, double, int, int);
        void average_entropy();
        void average_entropy_switch(double, int, int);
        void entropy(double, double, double, double, double, double);

};
