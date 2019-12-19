#include "atom.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include "string.h"
#include <chrono>
#include <pybind11/stl.h>

//functions for atoms
//-------------------------------------------------------------------------------------------------------------------------
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
void Atom::snneighbors(int dd){

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
void Atom::ssolid(int idd){
  if (!((idd == 0) || (idd == 1))){
    throw invalid_argument("surface should be 1 or 0");
  }
  issolid=idd; }

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


//structure properties
int Atom::gsurface() {return issurface; }
int Atom::gcluster() {return belongsto; }
void Atom::ssurface( int val) {
  if (!((val == 0) || (val == 1))){
    throw invalid_argument("surface should be 1 or 0");
  }
  issurface = val; }

void Atom::scluster( int val) {
  belongsto = val; }

int Atom::glcluster() {return lcluster; }
void Atom::slcluster( int val) {
  if (!((val == 0) || (val == 1))){
    throw invalid_argument("largest_cluster should be 1 or 0");
  }

  lcluster = val; }


py::dict Atom::gcustom(){

    return custom;
}

void Atom::scustom(py::dict cc){
    custom = cc;
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


//takes int and returns double - one value

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

    //loop over face vertices
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

double Atom::gangular(){
    return angular;
}

void Atom::sangular(double dd){
    angular = dd;
}

double Atom::gavgangular(){
    return avg_angular;
}

void Atom::savgangular(double dd){
    avg_angular = dd;
}


vector<int> Atom::gchiparams(){
  return chiparams;
}

void Atom::schiparams(vector<int> nns){
  chiparams.clear();
  chiparams = nns;
}

double Atom::gdisorder(){
    return disorder;
}

void Atom::sdisorder(double dd){
    disorder = dd;
}

double Atom::gavgdisorder(){
    return avgdisorder;
}

void Atom::savgdisorder(double dd){
    avgdisorder = dd;
}

vector<double> Atom::gsro(){
    return sro;
}

void Atom::ssro(vector<double> dd){
    sro = dd;
}
