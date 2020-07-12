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
Atom::Atom( vector<double> pos, int idd, int typ){

    posx = pos[0];
    posy = pos[1];
    posz = pos[2];
    id = idd;
    type = typ;

    //assign other values - the default ones
    belongsto = -1;
    issolid = 0;
    issurface = 0;
    loc = 0;
    isneighborset = 0;
    n_neighbors = 0;
    lcluster = 0;
    head = -1;


    for (int tn = 0; tn<MAXNUMBEROFNEIGHBORS; tn++){
        neighbors[tn] = -1;
        neighbordist[tn] = -1.0;
        neighborweight[tn] = -1.0;
        facevertices[tn] = -1;
        faceverticenumbers[tn] = -1;
        faceperimeters[tn] = -1.0;
        sij[tn] = -1.0;
        //edgelengths[tn] = -1.0;

    }

    for (int tn = 0; tn<11; tn++){
        q[tn] = -1;
        aq[tn] = -1;

        for (int tnn =0; tnn<25; tnn++){
            realq[tn][tnn] = -1;
            imgq[tn][tnn] = -1;
            arealq[tn][tnn] = -1;
            aimgq[tn][tnn] = -1;
        }
    }

}

Atom::~Atom(){ }

//-------------------------------------------------------
// Basic Atom properties
//-------------------------------------------------------


//aceesss funcs
vector<double> Atom::gx(){
    vector<double> pos;
    pos.emplace_back(posx);
    pos.emplace_back(posy);
    pos.emplace_back(posz);
    return pos;
}

void Atom::sx(vector<double> rls){
    posx = rls[0];
    posy = rls[1];
    posz = rls[2];
}



//-------------------------------------------------------
// Neighbor related properties
//-------------------------------------------------------
vector<int> Atom::gneighbors(){
    vector<int> nn;
    nn.reserve(n_neighbors);
    for(int i=0;i<n_neighbors;i++){
        nn.emplace_back(neighbors[i]);
    }
    return nn;
}

void Atom::sneighdist(vector<double> dd){
}

vector<double> Atom::gneighdist(){
  vector<double> neighdist;
  for(int i=0; i<n_neighbors; i++){
    neighdist.emplace_back(neighbordist[i]);
  }
  return neighdist;
}

void Atom::sneighbors(vector<int> nns){

    //first reset all neighbors
    for (int i = 0;i<MAXNUMBEROFNEIGHBORS;i++){
        neighbors[i] = NILVALUE;
        neighbordist[i] = -1.0;
    }

    //now assign the neighbors
    for(int i=0; i<nns.size(); i++){
        neighbors[i] = nns[i];
        //auto assign weight to 1
        neighborweight[i] = 1.00;
    }

    n_neighbors = nns.size();
    isneighborset = 1;

}

void Atom::sneighborweights(vector<double> nss){
    for(int i=0; i<nss.size(); i++){
        neighborweight[i] = nss[i];
    }
}

vector<double> Atom::gneighborweights(){
    vector <double> rqlms;
    for(int i=0; i<n_neighbors; i++){
        rqlms.emplace_back(neighborweight[i]);
    }
    return rqlms;
}

void Atom::sdistvecs(vector<vector<double>> nss){

}

vector<vector<double>> Atom::gdistvecs(){
    vector<vector<double>> m1;
    vector <double> m2;

    for(int i=0; i<n_neighbors; i++){
        m2.clear();
        m2.emplace_back(n_diffx[i]);
        m2.emplace_back(n_diffy[i]);
        m2.emplace_back(n_diffz[i]);

        m1.emplace_back(m2);
    }
    return m1;
}

void Atom::slocalangles(vector<vector<double>> nss){

}

vector<vector<double>> Atom::glocalangles(){
    vector<vector<double>> m1;
    vector <double> m2;

    for(int i=0; i<n_neighbors; i++){
        m2.clear();
        m2.emplace_back(n_phi[i]);
        m2.emplace_back(n_theta[i]);

        m1.emplace_back(m2);
    }
    return m1;
}

//-------------------------------------------------------
// Q parameter properties
//-------------------------------------------------------
void Atom::ssij(vector<double> dd){
}

vector<double> Atom::gsij(){
  vector<double> ss;
  for(int i=0; i<n_neighbors; i++){
    ss.emplace_back(sij[i]);
  }
  return ss;
}

vector<double> Atom::gallq(){
    vector<double> allq;
    for(int i=0; i<11; i++){
        allq.emplace_back(q[i]);
    }
    return allq;
}

vector<double> Atom::gallaq(){
    vector<double> allq;
    for(int i=0; i<11; i++){
        allq.emplace_back(aq[i]);
    }
    return allq;
}

void Atom::sallq(vector<double> allq){
    for(int i=0; i<11; i++){
        q[i] = allq[i];
    }
}

void Atom::sallaq(vector<double> allaq){
    for(int i=0; i<11; i++){
        aq[i] = allaq[i];
    }
}
double Atom::gq(int qq){ return q[qq-2]; }
void Atom::sq(int qq, double qval){ q[qq-2] = qval; }

double Atom::gaq(int qq){ return aq[qq-2]; }
void Atom::saq(int qq, double qval){ aq[qq-2] = qval; }

double Atom::gq_big(int qval, bool averaged){

    if ((qval < 2) || (qval > 12)){
        throw invalid_argument("q value should be between 2-12");
    }
    if(averaged == true) { return gaq(qval);}
    else {return gq(qval);}
}

void Atom::sq_big(int qval, double val,  bool averaged){

    if ((qval < 2) || (qval > 12)){
        throw invalid_argument("q value should be between 2-12");
    }
    if(averaged == true) { saq(qval, val);}
    else { sq(qval, val);}
}

//overloaded version which takes a vector
vector<double> Atom::gq_big(vector<int> qval, bool averaged ){
    int d;
    if(averaged == true) {
        vector<double> retvals;
        for(int i=0; i<qval.size(); i++){
            if ((qval[i] < 2) || (qval[i] > 12)){
                throw invalid_argument("q value should be between 2-12");
            }
            retvals.push_back(gaq(qval[i]));
        }
        return retvals;
    }
    else {
        vector<double> retvals;
        for(int i=0; i<qval.size(); i++){
            if ((qval[i] < 2) || (qval[i] > 12)){
                throw invalid_argument("q value should be between 2-12");
            }
            retvals.push_back(gq(qval[i]));
        }
        return retvals;
    }
}

//overloaded version which takes a vector
void Atom::sq_big(vector<int> qval, vector<double> setvals, bool averaged){

    if(averaged == true) {

        for(int i=0; i<qval.size(); i++){
            if ((qval[i] < 2) || (qval[i] > 12)){
                throw invalid_argument("q value should be between 2-12");
            }

            saq(qval[i], setvals[i]);
        }
    }
    else {

        for(int i=0; i<qval.size(); i++){
            if ((qval[i] < 2) || (qval[i] > 12)){
                throw invalid_argument("q value should be between 2-12");
            }
            sq(qval[i], setvals[i]);
        }
    }
}


vector<complex<double>> Atom::get_qcomps(int qq, bool averaged){

  vector<complex<double>> qlms;
  qlms.reserve(2*qq+1);
  if (!averaged){
    for(int i=0;i<(2*qq+1);i++){
        complex<double> lmval(realq[qq-2][i], imgq[qq-2][i]);
        qlms.emplace_back(lmval);
    }
  }
  else{
    for(int i=0;i<(2*qq+1);i++){
        complex<double> lmval(arealq[qq-2][i], aimgq[qq-2][i]);
        qlms.emplace_back(lmval);
    }
  }
  return qlms;
}


//-------------------------------------------------------
// Solid related properties
//-------------------------------------------------------




//-------------------------------------------------------
// Voronoi related properties
//-------------------------------------------------------
void Atom::sfacevertices(vector<int> nss){
    for(int i=0; i<nss.size(); i++){
        facevertices[i] = nss[i];
    }
}

vector<int> Atom::gfacevertices(){
    vector <int> rqlms;
    for(int i=0; i<n_neighbors; i++){
        rqlms.emplace_back(facevertices[i]);
    }
    return rqlms;
}

void Atom::sfaceperimeters(vector<double> nss){
    for(int i=0; i<nss.size(); i++){
        faceperimeters[i] = nss[i];
    }
}

vector<double> Atom::gfaceperimeters(){
    vector <double> rqlms;
    for(int i=0; i<n_neighbors; i++){
        rqlms.emplace_back(faceperimeters[i]);
    }
    return rqlms;
}

void Atom::sedgelengths(vector<vector<double>> nss){
    edgelengths.clear();
    edgelengths = nss;
}

vector<vector<double>> Atom::gedgelengths(){
    return edgelengths;
}

vector<int> Atom::gvorovector(){
    vector<int> voro;
    voro.emplace_back(n3);
    voro.emplace_back(n4);
    voro.emplace_back(n5);
    voro.emplace_back(n6);
    return voro;
}

void Atom::svorovector(vector<int> voro){

    n3 = voro[0];
    n4 = voro[1];
    n5 = voro[2];
    n6 = voro[3];
}

//-------------------------------------------------------
// Angle related properties
//-------------------------------------------------------


//-------------------------------------------------------
// Other order parameters
//-------------------------------------------------------
void Atom::find_common_neighbors(){
    
    common_neighbors.clear();
    common_neighbor_count.clear();
    common_neighbors.resize(n_neighbors);
    common_neighbor_count.resize(n_neighbors);

    int main_n;
    int sub_n;

    //for each neighbor
    for(int i=0; i< n_neighbors; i++){
        //loop over that neighbor
        common_neighbor_count[i] = 0;

        for(int j=0; j < next_neighbor_counts[i]; j++){
            main_n = next_neighbors[i][j];

            //now loop over self neighbors again
            for(int k=0; k < n_neighbors; k++){
                sub_n = neighbors[k];
                if (main_n == sub_n){
                    common_neighbor_count[i]++;
                    common_neighbors[i].emplace_back(k);
                }                  
            }
        }
    }

}

   
void Atom::find_bonds_of_common_neighbors(){

    common_neighbor_bonds.clear();
    common_neighbor_bond_count.clear();
    common_neighbor_bonds.resize(n_neighbors);
    common_neighbor_bond_count.resize(n_neighbors);

    int main_n, sub_n;
    int indexmain, indexsub;

    vector<int> pair;
    pair.resize(2);

    for(int i=0; i< n_neighbors; i++){
        //now take the corresponding common neighbors
        common_neighbor_bond_count[i] = 0;
        for(int j=0; j<common_neighbor_count[i]; j++){
            //hcovert it to index
            indexmain = common_neighbors[i][j];
            main_n = neighbors[indexmain];
            //this is to avoid double counting
            for(int k=j+1; k<common_neighbor_count[i]; k++){
                //if (j==k) continue;
                //convert to index
                indexsub = common_neighbors[i][k];

                //check if they are each others neighbors
                //as in, if j appears in k's neighbors
                for(int l=0; l < next_neighbor_counts[indexsub]; l++){
                    sub_n = next_neighbors[indexsub][l];
                    if (main_n == sub_n){
                        common_neighbor_bond_count[i] ++;
                        pair[0] = main_n;
                        pair[1] = neighbors[indexsub];
                        common_neighbor_bonds[i].emplace_back(pair);
                        break;
                    }
                }
            }
        }
    }    
}

void Atom::find_bond_chains(){
    
    bond_chain_count.clear();
    bond_chain_count.resize(n_neighbors);

    int max_length;
    int count;
    vector<int> indicator;
    vector<int> path;
    int pathlength;
    int chains;
    int finished;

    //we need a new algo, which loops over
    //first loop is always over neighbors
    for(int n=0; n<n_neighbors; n++){
        //now - for each neighbor, we have a bond chain
        //we need a replica array of ones
        //we need a count which we will reduce
        //max length of chain
        max_length = common_neighbor_bond_count[n];
        //chain number
        chains = 0;
        //set the count, which we will reduce
        count = max_length;
        //now make array of ones
        indicator.clear();
        for(int i=0; i<max_length; i++){
            indicator.emplace_back(1);
        }
        //we are reading - start the first while loop
        while(true){
            //clear the path
            path.clear();
            pathlength = 0;
            //if no more terms are to be added, break
            if (count == 0) break;
            //add first elements to path
            if (indicator[chains]==0) chains++;

            path.emplace_back(common_neighbor_bonds[n][chains][0]);
            path.emplace_back(common_neighbor_bonds[n][chains][1]);
            pathlength++;
            //reduce count and set this term to added
            count--;
            indicator[chains] = 0;

            //start the next loop
            while(true){
                //set finished to true
                finished = 1;
                //now loop over bonds
                for(int i=0; i<max_length; i++){
                    //is this term valid?
                    if(!indicator[i]) continue;
                    //otherwise check if it is in path
                    //loop over path
                    int rr = pathlength*2;
                    for(int j=0; j<rr; j++){
                        if((common_neighbor_bonds[n][i][0] == path[j]) or (common_neighbor_bonds[n][i][1] == path[j])){
                            //add the values
                            path.emplace_back(common_neighbor_bonds[n][i][0]);
                            path.emplace_back(common_neighbor_bonds[n][i][1]);
                            pathlength++;
                            finished = 0;
                            indicator[i] = 0;
                            count--;
                        }    
                    }
                }
                if (finished == 1) break;
            }
            if (count <= pathlength) break;

        }
        bond_chain_count[n] = pathlength-1;

    }

}