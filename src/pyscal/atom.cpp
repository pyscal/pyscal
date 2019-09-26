#include "atom.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include "string.h"
#include <chrono>

//functions for atoms
//-------------------------------------------------------------------------------------------------------------------------
Atom::Atom( vector<double> pos, int idd, int typ){

    posx = pos[0];
    posy = pos[1];
    posz = pos[2];
    id = idd;
    type = typ;
}

Atom::~Atom(){ }


vector<int> Atom::gneighbors(){
    vector<int> nn;
    nn.reserve(n_neighbors);
    for(int i=0;i<n_neighbors;i++){
        nn.emplace_back(neighbors[i]);
    }
    return nn;
}

int Atom::gnneighbors(){
    return n_neighbors;
}


int Atom::gid(){ return id; }
int Atom::gfrenkelnumber(){ return frenkelnumber; }
void Atom::sfrenkelnumber(int nn){ frenkelnumber=nn; }
void Atom::sid(int idd){ id=idd; }
int Atom::gloc(){ return loc; }
void Atom::sloc(int idd){ loc=idd; }
int Atom::gtype(){ return type; }
void Atom::stype(int idd){ type=idd; }
double Atom::gvolume(){ return volume; }
void Atom::svolume(double vv){ volume = vv; }
double Atom::gasij(){ return avq6q6; }
void Atom::sasij(double vv){ avq6q6 = vv; }
double Atom::gavgvolume(){ return avgvolume; }
void Atom::savgvolume(double vv){ avgvolume = vv; }
int Atom::gsolid(){ return issolid; }
void Atom::ssolid(int idd){ issolid=idd; }
int Atom::gstructure(){ return structure; }
void Atom::sstructure(int idd){ structure=idd; }
void Atom::scondition(int idd){ condition=idd; }
int Atom::gcondition(){ return condition; }

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

vector<int> Atom::gcluster(){
    vector<int> cl;
    cl.emplace_back(issolid);
    cl.emplace_back(issurface);
    cl.emplace_back(lcluster);
    cl.emplace_back(belongsto);
    return cl;
}

void Atom::scluster(vector<int> c1){
    issolid = c1[0];
    issurface = c1[1];
    lcluster = c1[2];
    belongsto = c1[3];
}




vector<vector <double>> Atom::gqlm(int qq) {

    vector< vector <double>> qlms;
    vector <double> rqlms;
    vector <double> iqlms;
    qlms.reserve(2);
    rqlms.reserve(2*qq+1);
    iqlms.reserve(2*qq+1);

    for(int i=0;i<(2*qq+1);i++){
        rqlms.emplace_back(realq[qq-2][i]);
        iqlms.emplace_back(imgq[qq-2][i]);
    }

    qlms.emplace_back(rqlms);
    qlms.emplace_back(iqlms);

    return qlms;

}

//for q vals
double Atom::gq(int qq){ return q[qq-2]; }
void Atom::sq(int qq, double qval){ q[qq-2] = qval; }

//vector<double> Atom::gqvec(vector<int> qq){

//  return q[qq-2];
//}
//void Atom::sq(int qq, double qval){ q[qq-2] = qval; }

//takes vector<int> and returns vector<double>
//vector<double> Atom::gaq(int qq){ return aq[qq-2]; }
//void Atom::saq(int qq, double qval){ aq[qq-2] = qval; }


//takes int and returns double - one value
double Atom::gaq(int qq){ return aq[qq-2]; }
void Atom::saq(int qq, double qval){ aq[qq-2] = qval; }

vector<vector <double>> Atom::gaqlm(int qq) {

    vector< vector <double>> qlms;
    vector <double> rqlms;
    vector <double> iqlms;
    qlms.reserve(2);
    rqlms.reserve(2*qq+1);
    iqlms.reserve(2*qq+1);

    for(int i=0;i<(2*qq+1);i++){
        rqlms.emplace_back(arealq[qq-2][i]);
        iqlms.emplace_back(aimgq[qq-2][i]);
    }

    qlms.emplace_back(rqlms);
    qlms.emplace_back(iqlms);

    return qlms;

}

void Atom::scustom(vector<double> cvals) {
    for(int i=0; i<cvals.size(); i++){
        custom.emplace_back(cvals[i]);
    }
    //custom = cvals;
}

vector<double> Atom::gcustom() {
    vector <double> rqlms;
    for(int i=0; i<custom.size(); i++){
        rqlms.emplace_back(custom[i]);
    }
    //rqlms = custom;
    return rqlms;
}

//functions to set the neighbors for each atoms
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

void Atom::svertex_numbers(vector<int> nss){
    vertex_numbers = nss;
}

vector<int> Atom::gvertex_numbers(){
    return vertex_numbers;
}

void Atom::svertex_vectors(vector<double> nss){
    vertex_vectors = nss;
}

vector<double> Atom::gvertex_vectors(){
    return vertex_vectors;
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
