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

using namespace std;

const int MAXNUMBEROFNEIGHBORS = 100;
const double PI = 3.14159265;
const int nilvalue = 33333333;
const double pi = 3.141592653589793;
const long int FACTORIALS[17] = {1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600,6227020800,87178291200,1307674368000,20922789888000};


class Atom {
  private:
    //The position of the particle
  public:
    Atom();
    virtual ~Atom();
    double posx,posy,posz;
    
    int neighbors[MAXNUMBEROFNEIGHBORS];
    double neighbordist[MAXNUMBEROFNEIGHBORS];
    //JR: also save the distance in cartesian and spherical coordinates
    double n_diffx[MAXNUMBEROFNEIGHBORS];
    double n_diffy[MAXNUMBEROFNEIGHBORS];
    double n_diffz[MAXNUMBEROFNEIGHBORS];
    double n_r[MAXNUMBEROFNEIGHBORS];
    double n_phi[MAXNUMBEROFNEIGHBORS];
    double n_theta[MAXNUMBEROFNEIGHBORS];

    int n_neighbors;

    double realQ4[9],imgQ4[9];
    double realQ6[13],imgQ6[13];
    double arealQ4[9],aimgQ4[9];
    double arealQ6[13],aimgQ6[13];
    
    double frenkelnumber;
    double avq6q6;

    int belongsto;
    int issolid;
    int structure;
    int id;
};


class System {
  public:
    void convert_SphericalCoordinates(double , double , double , double &, double &, double &);
    double PLM(int, int, double);
    void YLM(int , int , double , double , double &, double &);
    void QLM(int ,int ,double ,double ,double &, double & );
    void get_AllNeighborsAndDistances();
    void calculate_complexQLM_6();
    double get_NumberFromBond(int,int);
    void calculate_frenkelNumbers();
    double get_absDistance(int,int,double&,double&,double&);
    System();
    virtual ~System();
    //Properties of one single molecule
    Atom* molecules;
    
    //Init the system
    //void initializeMolecules(int);

    void readParticleFile();
    int calculate_largestClusterparameter_Full();	//variant of function above
    int clusterCriterium(int,int );
    void find_solids();
    void find_clusters();
    int largest_cluster();
    void set_minfrenkel(int);
    void set_inputfile(string);
    void set_neighbordistance(double);
    //old params
    int nop;
    int baseunit;
    int minfrenkel;
    double boxx, boxy, boxz;
    string inputfile;
    double neighbordistance;
    double threshold;
    double avgthreshold;

};
