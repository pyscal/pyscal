#include "system.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include "voro++.hh"
#include "string.h"
#include <chrono>
#include <pybind11/stl.h>

using namespace voro;

//-----------------------------------------------------
// Constructor, Destructor and Access functions
//-----------------------------------------------------
System::System(){

    nop = -1;
    maxclusterid = -1;
    neighborsfound = 0;
    qsfound = 0;
    fileread = 0;
    filter = 0;
    triclinic = 0;
    alpha = 1;
    voronoiused = 0;
    solidq = 6;
    criteria = 0;
    usecells = 0;
    neighbordistance = 0;

}

System::~System(){

}


void System::read_particle_file(string nn){

    fileread = 1;

}



//-----------------------------------------------------
// Simulation box related methods
//-----------------------------------------------------
void System::assign_triclinic_params(vector<vector<double>> drot, vector<vector<double>> drotinv){

    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            rot[i][j] = drot[i][j];
            rotinv[i][j] = drotinv[i][j];
        }
    }

    triclinic = 1;
}

vector<vector<double>> System::get_triclinic_params(){

    vector<vector<double>> drot;
    vector<double> dummydrot;
    for(int i=0; i<3; i++){
        dummydrot.clear();
        for(int j=0; j<3; j++){
            dummydrot.emplace_back(rot[i][j]);
        }
        drot.emplace_back(dummydrot);
    }
    return drot;
}

void System::sbox(vector<vector <double>> boxd) {

    boxdims[0][0] = boxd[0][0];
    boxdims[0][1] = boxd[0][1];
    boxdims[1][0] = boxd[1][0];
    boxdims[1][1] = boxd[1][1];
    boxdims[2][0] = boxd[2][0];
    boxdims[2][1] = boxd[2][1];

    boxx = boxd[0][1] - boxd[0][0];
    boxy = boxd[1][1] - boxd[1][0];
    boxz = boxd[2][1] - boxd[2][0];
}

vector<vector<double>> System::gbox(){
    vector<vector<double>> qres;
    vector<double> qd;

    for(int i=0;i<3;i++){
        qd.clear();
        for(int j=0;j<2;j++){
            qd.emplace_back(boxdims[i][j]);
        }
        qres.emplace_back(qd);
    }
    return qres;
}

vector<vector<double>> System::gboxvecs(){
    vector<vector<double>> qres;
    vector<double> dqres;
    if (triclinic==1){
        for(int i=0; i<3; i++){
            dqres.clear();
            for(int j=0; j<3; j++){
                dqres.emplace_back(rot[j][i]);
            }
            qres.emplace_back(dqres);
        }
    }
    else{
        for(int i=0; i<3; i++){
            dqres.clear();
            for(int j=0; j<3; j++){
                if(i==j){
                    dqres.emplace_back(boxdims[i][1]-boxdims[i][0]);
                }
                else{
                    dqres.emplace_back(0.0);
                }
            }
            qres.emplace_back(dqres);
        }
    }
    return qres;
}


//-----------------------------------------------------
// Atom related methods
//-----------------------------------------------------
//this function allows for handling custom formats of atoms and so on
void System::set_atoms( vector<Atom> atomitos){

    atoms.clear();
    nop = atomitos.size();
    atoms.reserve(nop);
    atoms.assign(atomitos.begin(), atomitos.end());

}


//this function allows for handling custom formats of atoms and so on
vector<Atom> System::get_atoms( ){
    return atoms;

}

Atom System::gatom(int i) { return atoms[i]; }
void System::satom(Atom atom1) {
    int idd = atom1.loc;
    atoms[idd] = atom1;
}

//----------------------------------------------------
// Neighbor methods
//----------------------------------------------------
void System::susecells(int a){

    usecells = a;
}

int System::gusecells(){

    return usecells ;
}

double System::get_abs_distance(int ti ,int tj,double &diffx ,double &diffy,double &diffz){

    double abs, ax, ay, az;
    diffx = atoms[tj].posx - atoms[ti].posx;
    diffy = atoms[tj].posy - atoms[ti].posy;
    diffz = atoms[tj].posz - atoms[ti].posz;

    if (triclinic == 1){

        //convert to the triclinic system
        ax = rotinv[0][0]*diffx + rotinv[0][1]*diffy + rotinv[0][2]*diffz;
        ay = rotinv[1][0]*diffx + rotinv[1][1]*diffy + rotinv[1][2]*diffz;
        az = rotinv[2][0]*diffx + rotinv[2][1]*diffy + rotinv[2][2]*diffz;

        //scale to match the triclinic box size
        diffx = ax*boxx;
        diffy = ay*boxy;
        diffz = az*boxz;

        //now check pbc
        //nearest image
        if (diffx> boxx/2.0) {diffx-=boxx;};
        if (diffx<-boxx/2.0) {diffx+=boxx;};
        if (diffy> boxy/2.0) {diffy-=boxy;};
        if (diffy<-boxy/2.0) {diffy+=boxy;};
        if (diffz> boxz/2.0) {diffz-=boxz;};
        if (diffz<-boxz/2.0) {diffz+=boxz;};

        //now divide by box vals - scale down the size
        diffx = diffx/boxx;
        diffy = diffy/boxy;
        diffz = diffz/boxz;

        //now transform back to normal system
        ax = rot[0][0]*diffx + rot[0][1]*diffy + rot[0][2]*diffz;
        ay = rot[1][0]*diffx + rot[1][1]*diffy + rot[1][2]*diffz;
        az = rot[2][0]*diffx + rot[2][1]*diffy + rot[2][2]*diffz;

        //now assign to diffs and calculate distnace
        diffx = ax;
        diffy = ay;
        diffz = az;

        //finally distance
        abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

    }
    else{
        //nearest image
        if (diffx> boxx/2.0) {diffx-=boxx;};
        if (diffx<-boxx/2.0) {diffx+=boxx;};
        if (diffy> boxy/2.0) {diffy-=boxy;};
        if (diffy<-boxy/2.0) {diffy+=boxy;};
        if (diffz> boxz/2.0) {diffz-=boxz;};
        if (diffz<-boxz/2.0) {diffz+=boxz;};
        abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
    }
    return abs;
}

//function for binding
double System::get_abs_distance(Atom atom1 , Atom atom2 ){

    double abs, ax, ay, az;
    double diffx = atom1.posx - atom2.posx;
    double diffy = atom1.posy - atom2.posy;
    double diffz = atom1.posz - atom2.posz;

    if (triclinic == 1){

        //convert to the triclinic system
        ax = rotinv[0][0]*diffx + rotinv[0][1]*diffy + rotinv[0][2]*diffz;
        ay = rotinv[1][0]*diffx + rotinv[1][1]*diffy + rotinv[1][2]*diffz;
        az = rotinv[2][0]*diffx + rotinv[2][1]*diffy + rotinv[2][2]*diffz;

        //scale to match the triclinic box size
        diffx = ax*boxx;
        diffy = ay*boxy;
        diffz = az*boxz;

        //now check pbc
        //nearest image
        if (diffx> boxx/2.0) {diffx-=boxx;};
        if (diffx<-boxx/2.0) {diffx+=boxx;};
        if (diffy> boxy/2.0) {diffy-=boxy;};
        if (diffy<-boxy/2.0) {diffy+=boxy;};
        if (diffz> boxz/2.0) {diffz-=boxz;};
        if (diffz<-boxz/2.0) {diffz+=boxz;};

        //now divide by box vals - scale down the size
        diffx = diffx/boxx;
        diffy = diffy/boxy;
        diffz = diffz/boxz;

        //now transform back to normal system
        ax = rot[0][0]*diffx + rot[0][1]*diffy + rot[0][2]*diffz;
        ay = rot[1][0]*diffx + rot[1][1]*diffy + rot[1][2]*diffz;
        az = rot[2][0]*diffx + rot[2][1]*diffy + rot[2][2]*diffz;

        //now assign to diffs and calculate distnace
        diffx = ax;
        diffy = ay;
        diffz = az;

        //finally distance
        abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

    }
    else{

        //nearest image
        if (diffx> boxx/2.0) {diffx-=boxx;};
        if (diffx<-boxx/2.0) {diffx+=boxx;};
        if (diffy> boxy/2.0) {diffy-=boxy;};
        if (diffy<-boxy/2.0) {diffy+=boxy;};
        if (diffz> boxz/2.0) {diffz-=boxz;};
        if (diffz<-boxz/2.0) {diffz+=boxz;};
        abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
    }

    return abs;
}

//function for binding
vector<double> System::get_distance_vector(Atom atom1 , Atom atom2 ){

    double ax, ay, az;
    double diffx = atom1.posx - atom2.posx;
    double diffy = atom1.posy - atom2.posy;
    double diffz = atom1.posz - atom2.posz;

    if (triclinic == 1){

        //convert to the triclinic system
        ax = rotinv[0][0]*diffx + rotinv[0][1]*diffy + rotinv[0][2]*diffz;
        ay = rotinv[1][0]*diffx + rotinv[1][1]*diffy + rotinv[1][2]*diffz;
        az = rotinv[2][0]*diffx + rotinv[2][1]*diffy + rotinv[2][2]*diffz;

      double dummy;
        //scale to match the triclinic box size
        diffx = ax*boxx;
        diffy = ay*boxy;
        diffz = az*boxz;

        //now check pbc
        //nearest image
        if (diffx> boxx/2.0) {diffx-=boxx;};
        if (diffx<-boxx/2.0) {diffx+=boxx;};
        if (diffy> boxy/2.0) {diffy-=boxy;};
        if (diffy<-boxy/2.0) {diffy+=boxy;};
        if (diffz> boxz/2.0) {diffz-=boxz;};
        if (diffz<-boxz/2.0) {diffz+=boxz;};

        //now divide by box vals - scale down the size
        diffx = diffx/boxx;
        diffy = diffy/boxy;
        diffz = diffz/boxz;

        //now transform back to normal system
        ax = rot[0][0]*diffx + rot[0][1]*diffy + rot[0][2]*diffz;
        ay = rot[1][0]*diffx + rot[1][1]*diffy + rot[1][2]*diffz;
        az = rot[2][0]*diffx + rot[2][1]*diffy + rot[2][2]*diffz;

        //now assign to diffs and calculate distnace
        diffx = ax;
        diffy = ay;
        diffz = az;

    }
    else{

        //nearest image
        if (diffx> boxx/2.0) {diffx-=boxx;};
        if (diffx<-boxx/2.0) {diffx+=boxx;};
        if (diffy> boxy/2.0) {diffy-=boxy;};
        if (diffy<-boxy/2.0) {diffy+=boxy;};
        if (diffz> boxz/2.0) {diffz-=boxz;};
        if (diffz<-boxz/2.0) {diffz+=boxz;};

    }

    vector<double> abs;
    abs.emplace_back(diffx);
    abs.emplace_back(diffy);
    abs.emplace_back(diffz);

    return abs;
}


void System::reset_all_neighbors(){
    for (int ti = 0;ti<nop;ti++){

        atoms[ti].n_neighbors=0;
        atoms[ti].temp_neighbors.clear();

        for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++){

            atoms[ti].neighbors[tn] = NILVALUE;
            atoms[ti].neighbordist[tn] = -1.0;
        }
    }
}

void System::sfilter(int fno){

    filter = fno;
}
int System::gfilter(){
  return filter;
}

vector<double> System::get_pairdistances(){

    vector<double> res;
    double d;
    double diffx,diffy,diffz;

    for (int ti=0; ti<nop; ti++){
        for (int tj=ti; tj<nop; tj++){
            if(ti==tj) { continue; }
            d = get_abs_distance(ti,tj,diffx,diffy,diffz);
            res.emplace_back(d);

        }
    }
    return res;
}

//function to create cell lists
//snmall function that returns cell index when provided with cx, cy, cz
int System::cell_index(int cx, int cy, int cz){
    return cx*ny*nz + cy*nz + cz;
}


//if number of particles are small, use brute force
//if box is triclinic, use brute force
void System::set_up_cells(){

      int si,sj,sk, maincell, subcell;
      vector<int> cc;
      //find of all find the number of cells in each direction
      nx = boxx/neighbordistance;
      ny = boxy/neighbordistance;
      nz = boxz/neighbordistance;

      //now use this to find length of cell in each direction
      double lx = boxx/nx;
      double ly = boxy/ny;
      double lz = boxz/nz;

      //find the total number of cells
      total_cells = nx*ny*nz;

      //create a vector of cells
      cells = new cell[total_cells];

      //now run over and for each cell create its neighbor cells
      //all neighbor cells are also added
      for(int i=0; i<nx; i++){
         for(int j=0; j<ny; j++){
           for(int k=0; k<nz; k++){
              maincell = cell_index(i, j, k);
              for(int si=i-1; si<=i+1; si++){
                  for(int sj=j-1; sj<=j+1; sj++){
                      for(int sk=k-1; sk<=k+1; sk++){
                         cc = cell_periodic(si, sj, sk);
                         subcell = cell_index(cc[0], cc[1], cc[2]);
                         //add this to the list of neighbors
                         cells[maincell].neighbor_cells.emplace_back(subcell);

                      }
                  }
              }
           }
         }
      }

      int cx, cy, cz;
      double dx, dy, dz;
      int ind;

      //now loop over all atoms and assign cells
      for(int ti=0; ti<nop; ti++){

          //calculate c indices for the atom
          dx = atoms[ti].posx;
          dy = atoms[ti].posy;
          dz = atoms[ti].posz;

          //now apply boxdims
          if (dx < 0) dx+=boxx;
          else if (dx > boxx) dx-=boxx;
          if (dy < 0) dy+=boxy;
          else if (dy > boxy) dy-=boxy;
          if (dz < 0) dz+=boxz;
          else if (dz > boxz) dz-=boxz;

          //now find c vals
          cx = dx/lx;
          cy = dy/ly;
          cz = dz/lz;

          //now get cell index
          ind = cell_index(cx, cy, cz);

          //got cell index
          //now add the atom to the corresponding cells
          cells[ind].members.emplace_back(ti);

      }

      //end of loop - all cells, the member atoms and neighboring cells are added
}

vector<int> System::cell_periodic(int i, int j, int k){
    vector<int> ci;
    //apply periodic conditions
    if (i<0) i = i + nx;
    else if (i>nx-1) i = i -nx;
    ci.emplace_back(i);
    if (j<0) j = j + ny;
    else if (j>ny-1) j = j -ny;
    ci.emplace_back(j);
    if (k<0) k = k + nz;
    else if (k>nz-1) k = k -nz;
    ci.emplace_back(k);
    return ci;

}

//get all neighbor info but using cell lists
void System::get_all_neighbors_cells(){

    voronoiused = 0;

    double d;
    double diffx,diffy,diffz;
    double r,theta,phi;
    int ti, tj;

    //first create cells
    set_up_cells();
    int maincell, subcell;

    //now loop to find distance
    for(int i=0; i<total_cells; i++){
        //now go over the neighbor cells
        //for each member in cell i
        for(int mi=0; mi<cells[i].members.size(); mi++){
            //now go through the neighbors
            ti = cells[i].members[mi];
            for(int j=0 ; j<cells[i].neighbor_cells.size(); j++){
               //loop through members of j
               subcell = cells[i].neighbor_cells[j];
               for(int mj=0; mj<cells[subcell].members.size(); mj++){
                  //now we have mj -> members/compare with
                  tj = cells[subcell].members[mj];
                  //compare ti and tj and add
                  if (ti < tj){
                      d = get_abs_distance(ti,tj,diffx,diffy,diffz);
                      if (d < neighbordistance){

                        if ((filter == 1) && (atoms[ti].type != atoms[tj].type)){
                            continue;
                        }
                        //process_neighbor(ti, tj);

                        atoms[ti].neighbors[atoms[ti].n_neighbors] = tj;
                        atoms[ti].neighbordist[atoms[ti].n_neighbors] = d;
                        //weight is set to 1.0, unless manually reset
                        atoms[ti].neighborweight[atoms[ti].n_neighbors] = 1.00;
                        atoms[ti].n_diffx[atoms[ti].n_neighbors] = diffx;
                        atoms[ti].n_diffy[atoms[ti].n_neighbors] = diffy;
                        atoms[ti].n_diffz[atoms[ti].n_neighbors] = diffz;
                        convert_to_spherical_coordinates(diffx, diffy, diffz, r, phi, theta);
                        atoms[ti].n_r[atoms[ti].n_neighbors] = r;
                        atoms[ti].n_phi[atoms[ti].n_neighbors] = phi;
                        atoms[ti].n_theta[atoms[ti].n_neighbors] = theta;
                        atoms[ti].n_neighbors += 1;
                        atoms[ti].cutoff = neighbordistance;

                        atoms[tj].neighbors[atoms[tj].n_neighbors] = ti;
                        atoms[tj].neighbordist[atoms[tj].n_neighbors] = d;
                        //weight is set to 1.0, unless manually reset
                        atoms[tj].neighborweight[atoms[tj].n_neighbors] = 1.00;
                        atoms[tj].n_diffx[atoms[tj].n_neighbors] = -diffx;
                        atoms[tj].n_diffy[atoms[tj].n_neighbors] = -diffy;
                        atoms[tj].n_diffz[atoms[tj].n_neighbors] = -diffz;
                        convert_to_spherical_coordinates(-diffx, -diffy, -diffz, r, phi, theta);
                        atoms[tj].n_r[atoms[tj].n_neighbors] = r;
                        atoms[tj].n_phi[atoms[tj].n_neighbors] = phi;
                        atoms[tj].n_theta[atoms[tj].n_neighbors] = theta;
                        atoms[tj].n_neighbors +=1;
                        atoms[tj].cutoff = neighbordistance;
                      }
                  }
               }

            }

        }
    }

    neighborsfound = 1;

}



void System::get_all_neighbors_normal(){


    //reset voronoi flag
    voronoiused = 0;

    double d;
    double diffx,diffy,diffz;
    double r,theta,phi;

    if (!fileread) { read_particle_file(inputfile); }

    for (int ti=0; ti<nop; ti++){
        for (int tj=ti; tj<nop; tj++){
            if(ti==tj) { continue; }
            d = get_abs_distance(ti,tj,diffx,diffy,diffz);
            if (d < neighbordistance){
                if ((filter == 1) && (atoms[ti].type != atoms[tj].type)){
                    continue;
                }
                //process_neighbor(ti, tj);

                atoms[ti].neighbors[atoms[ti].n_neighbors] = tj;
                atoms[ti].neighbordist[atoms[ti].n_neighbors] = d;
                //weight is set to 1.0, unless manually reset
                atoms[ti].neighborweight[atoms[ti].n_neighbors] = 1.00;
                atoms[ti].n_diffx[atoms[ti].n_neighbors] = diffx;
                atoms[ti].n_diffy[atoms[ti].n_neighbors] = diffy;
                atoms[ti].n_diffz[atoms[ti].n_neighbors] = diffz;
                convert_to_spherical_coordinates(diffx, diffy, diffz, r, phi, theta);
                atoms[ti].n_r[atoms[ti].n_neighbors] = r;
                atoms[ti].n_phi[atoms[ti].n_neighbors] = phi;
                atoms[ti].n_theta[atoms[ti].n_neighbors] = theta;
                atoms[ti].n_neighbors += 1;
                atoms[ti].cutoff = neighbordistance;

                atoms[tj].neighbors[atoms[tj].n_neighbors] = ti;
                atoms[tj].neighbordist[atoms[tj].n_neighbors] = d;
                //weight is set to 1.0, unless manually reset
                atoms[tj].neighborweight[atoms[tj].n_neighbors] = 1.00;
                atoms[tj].n_diffx[atoms[tj].n_neighbors] = -diffx;
                atoms[tj].n_diffy[atoms[tj].n_neighbors] = -diffy;
                atoms[tj].n_diffz[atoms[tj].n_neighbors] = -diffz;
                convert_to_spherical_coordinates(-diffx, -diffy, -diffz, r, phi, theta);
                atoms[tj].n_r[atoms[tj].n_neighbors] = r;
                atoms[tj].n_phi[atoms[tj].n_neighbors] = phi;
                atoms[tj].n_theta[atoms[tj].n_neighbors] = theta;
                atoms[tj].n_neighbors +=1;
                atoms[tj].cutoff = neighbordistance;
            }
        }

    }

    //mark end of neighbor calc
    neighborsfound = 1;

}

void System::process_neighbor(int ti, int tj){
    /*
    Calculate all info and add it to list
    ti - loc of host atom
    tj - loc of the neighbor
     d - interatomic distance
     */

    double d, diffx, diffy, diffz;
    double r, phi, theta;

    d = get_abs_distance(ti, tj, diffx,diffy,diffz);

    atoms[ti].neighbors[atoms[ti].n_neighbors] = tj;
    atoms[ti].neighbordist[atoms[ti].n_neighbors] = d;
    //weight is set to 1.0, unless manually reset
    atoms[ti].neighborweight[atoms[ti].n_neighbors] = 1.00;
    atoms[ti].n_diffx[atoms[ti].n_neighbors] = diffx;
    atoms[ti].n_diffy[atoms[ti].n_neighbors] = diffy;
    atoms[ti].n_diffz[atoms[ti].n_neighbors] = diffz;
    convert_to_spherical_coordinates(diffx, diffy, diffz, r, phi, theta);
    atoms[ti].n_r[atoms[ti].n_neighbors] = r;
    atoms[ti].n_phi[atoms[ti].n_neighbors] = phi;
    atoms[ti].n_theta[atoms[ti].n_neighbors] = theta;
    atoms[ti].n_neighbors += 1;

}

/*
To increase the speed of the other methods, we need some functions using cells and
otheriwse which adds atoms to the temp_neighbors list
*/
void System::get_temp_neighbors_brute(){

    //reset voronoi flag

    double d;
    double diffx,diffy,diffz;

    for (int ti=0; ti<nop; ti++){
        for (int tj=ti; tj<nop; tj++){
            if(ti==tj) { continue; }
            d = get_abs_distance(ti,tj,diffx,diffy,diffz);
            if (d <= neighbordistance){
                datom x = {d, tj};
                atoms[ti].temp_neighbors.emplace_back(x);
                datom y = {d, ti};
                atoms[tj].temp_neighbors.emplace_back(y);
            }
        }
    }

}

/*
Cells should only be used when the system has a minimum size - in this case,
about 2000 atoms.
*/
void System::get_temp_neighbors_cells(){

    //first create cells
    set_up_cells();

    int maincell, subcell;
    int ti, tj;
    double d;
    double diffx,diffy,diffz;

    //now loop to find distance
    for(int i=0; i<total_cells; i++){
        //now go over the neighbor cells
        //for each member in cell i
        for(int mi=0; mi<cells[i].members.size(); mi++){
            //now go through the neighbors
            ti = cells[i].members[mi];
            for(int j=0 ; j<cells[i].neighbor_cells.size(); j++){
               //loop through members of j
               subcell = cells[i].neighbor_cells[j];
               for(int mj=0; mj<cells[subcell].members.size(); mj++){
                  //now we have mj -> members/compare with
                  tj = cells[subcell].members[mj];
                  //compare ti and tj and add
                  if (ti < tj){
                      d = get_abs_distance(ti,tj,diffx,diffy,diffz);
                      if (d < neighbordistance){
                        datom x = {d, tj};
                        atoms[ti].temp_neighbors.emplace_back(x);
                        datom y = {d, ti};
                        atoms[tj].temp_neighbors.emplace_back(y);
                      }
                  }
               }

            }

        }
    }

}


int System::get_all_neighbors_sann(double prefactor){
    /*
    A new adaptive algorithm. Similar to the old ones, we guess a basic distance with padding,
    and sort them up.
    After that, we use the algorithm by in J. Chem. Phys. 136, 234107 (2012) to find the list of
    neighbors.
     */

    //reset voronoi flag
    voronoiused = 0;

    double d, dcut;
    double diffx,diffy,diffz;
    double r,theta,phi;
    int m, maxneighs, finished;
    finished = 1;

    vector<int> nids;
    vector<double> dists, sorted_dists;

    //double prefactor = 1.21;
    double summ;
    double boxvol;

    //some guesswork here
    //find the box volumes
    if (triclinic==1){
        double a1, a2, a3, b1, b2, b3, c1, c2, c3;
        //rot is the cell vectors transposed
        a1 = rot[0][0];
        a2 = rot[1][0];
        a3 = rot[2][0];
        b1 = rot[0][1];
        b2 = rot[1][1];
        b3 = rot[2][1];
        c1 = rot[0][2];
        c2 = rot[1][2];
        c3 = rot[2][2];
        boxvol = c1*(a2*b3-a3*b2) - c2*(a1*b3-b1*a3) + c3*(a1*b2-a2*b1);
    }
    else{
        boxvol = boxx*boxy*boxz;
    }


    //now find the volume per particle
    double guessvol = boxvol/float(nop);

    //guess the side of a cube that is occupied by an atom - this is a guess distance
    double guessdist = cbrt(guessvol);

    //now add some safe padding - this is the prefactor which we will read in
    guessdist = prefactor*guessdist;
    neighbordistance = guessdist;


    //if file is not read - read it in at this point
    //we have to work on indicator functions
    if (!fileread) { read_particle_file(inputfile); }

    if (usecells){
        get_temp_neighbors_cells();
    }
    else{
        get_temp_neighbors_brute();
    }

    for (int ti=0; ti<nop; ti++){
        if (atoms[ti].temp_neighbors.size() < 3){
            return 0;
        }

        sort(atoms[ti].temp_neighbors.begin(), atoms[ti].temp_neighbors.end(), by_dist());

        //start with initial routine
        m = 3;
        summ = 0;
        for(int i=0 ; i<m; i++){
            summ += atoms[ti].temp_neighbors[i].dist;
            int tj = atoms[ti].temp_neighbors[i].index;
            process_neighbor(ti, tj);
        }

        //find cutoff
        dcut = summ/float(m-2);
        maxneighs = atoms[ti].temp_neighbors.size();

        while( (m < maxneighs) && (dcut >= atoms[ti].temp_neighbors[m].dist)){
            //increase m
            m = m+1;

            //here now we can add this to the list neighbors and process things
            int tj = atoms[ti].temp_neighbors[m].index;
            process_neighbor(ti, tj);

            //find new dcut
            summ = summ + atoms[ti].temp_neighbors[m].dist;
            dcut = summ/float(m-2);
            atoms[ti].cutoff = dcut;
        }

        //find if there was an error
        if (m==maxneighs){
            finished = 0;
            break;
        }
        else{
            finished = 1;
        }

    }

    //mark end of neighbor calc
    neighborsfound = 1;
    return finished;


}



int System::get_all_neighbors_adaptive(double prefactor, int nlimit, double padding){

    double d, dcut;
    double diffx,diffy,diffz;
    double r,theta,phi;
    int m, maxneighs, finished;

    double summ;
    double boxvol;
    //some guesswork here
    //find the box volumes
    if (triclinic==1){
        double a1, a2, a3, b1, b2, b3, c1, c2, c3;
        //rot is the cell vectors transposed
        a1 = rot[0][0];
        a2 = rot[1][0];
        a3 = rot[2][0];
        b1 = rot[0][1];
        b2 = rot[1][1];
        b3 = rot[2][1];
        c1 = rot[0][2];
        c2 = rot[1][2];
        c3 = rot[2][2];
        boxvol = c1*(a2*b3-a3*b2) - c2*(a1*b3-b1*a3) + c3*(a1*b2-a2*b1);
    }
    else{
        boxvol = boxx*boxy*boxz;
    }

    //now find the volume per particle
    double guessvol = boxvol/float(nop);

    //guess the side of a cube that is occupied by an atom - this is a guess distance
    double guessdist = cbrt(guessvol);

    //now add some safe padding - this is the prefactor which we will read in
    guessdist = prefactor*guessdist;
    neighbordistance = guessdist;

    //if file is not read - read it in at this point
    //we have to work on indicator functions
    if (!fileread) { read_particle_file(inputfile); }


    //introduce cell lists here - instead of looping over all neighbors
    //use cells

    if (usecells){
        get_temp_neighbors_cells();
    }
    else{
        get_temp_neighbors_brute();
    }

    //end of callstructural competition
    //subatoms would now be populated
    //now starts the main loop
    for (int ti=0; ti<nop; ti++){
        //check if its zero size
        if (atoms[ti].temp_neighbors.size() < nlimit){
            return 0;
        }

        sort(atoms[ti].temp_neighbors.begin(), atoms[ti].temp_neighbors.end(), by_dist());

        summ = 0;
        for(int i=0; i<nlimit; i++){
            summ += atoms[ti].temp_neighbors[i].dist;
        }
        dcut = padding*(1.0/float(nlimit))*summ;

        //now we are ready to loop over again, but over the lists
        for(int j=0; j<atoms[ti].temp_neighbors.size(); j++){
            int tj = atoms[ti].temp_neighbors[j].index;
            if (atoms[ti].temp_neighbors[j].dist < dcut){

                if ((filter == 1) && (atoms[ti].type != atoms[tj].type)){
                    continue;
                }

                d = get_abs_distance(ti,tj,diffx,diffy,diffz);
                atoms[ti].neighbors[atoms[ti].n_neighbors] = tj;
                atoms[ti].neighbordist[atoms[ti].n_neighbors] =d;
                //weight is set to 1.0, unless manually reset
                atoms[ti].neighborweight[atoms[ti].n_neighbors] = 1.00;
                atoms[ti].n_diffx[atoms[ti].n_neighbors] = diffx;
                atoms[ti].n_diffy[atoms[ti].n_neighbors] = diffy;
                atoms[ti].n_diffz[atoms[ti].n_neighbors] = diffz;
                convert_to_spherical_coordinates(diffx, diffy, diffz, r, phi, theta);
                atoms[ti].n_r[atoms[ti].n_neighbors] = r;
                atoms[ti].n_phi[atoms[ti].n_neighbors] = phi;
                atoms[ti].n_theta[atoms[ti].n_neighbors] = theta;
                atoms[ti].n_neighbors += 1;
                atoms[ti].cutoff = dcut;

            }
        }

    }

    //mark end of neighbor calc
    neighborsfound = 1;
    return 1;

}

void System::set_neighbordistance(double nn) { neighbordistance = nn; }


//---------------------------------------------------
// Methods for q calculation
//---------------------------------------------------
double System::dfactorial(int l,int m){

    double fac = 1.00;
    for(int i=0;i<2*m;i++){
        fac*=double(l+m-i);
    }
    return (1.00/fac);
}

void System::set_reqd_qs(vector <int> qs){

    lenqs = qs.size();
    reqdqs = new int[lenqs];
    for(int i=0;i<lenqs;i++){
        reqdqs[i] = qs[i];
    }

    rq_backup = qs;
}


void System::set_reqd_aqs(vector <int> qs){

    lenaqs = qs.size();
    reqdaqs = new int[lenaqs];
    for(int i=0;i<lenaqs;i++){
        for(int j=0;j<lenqs;j++){
            if(qs[i]==reqdqs[j]) { reqdaqs[i] = qs[i]; }
        }
    }
    //only qvlaues in the normal set will be included in the aq list
    //check here if its in the qlist
    //cout<<"corresponding q value should also be set."<<endl;

}

double System::PLM(int l, int m, double x){

    double fact,pll,pmm,pmmp1,somx2;
    int i,ll;
    pll = 0.0;
    if (m < 0 || m > l || fabs(x) > 1.0)
        cerr << "impossible combination of l and m" << "\n";
    pmm=1.0;
    if (m > 0){
        somx2=sqrt((1.0-x)*(1.0+x));
        fact=1.0;
        for (i=1;i<=m;i++){
            pmm *= -fact*somx2;
            fact += 2.0;
        }
    }

    if (l == m)
        return pmm;
    else{
        pmmp1=x*(2*m+1)*pmm;
        if (l == (m+1))
            return pmmp1;
        else{
            for (ll=m+2;ll<=l;ll++){
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
            pmm=pmmp1;
            pmmp1=pll;
            }
        return pll;
        }
    }
}

void System::convert_to_spherical_coordinates(double x, double y, double z, double &r, double &phi, double &theta){
    r = sqrt(x*x+y*y+z*z);
    theta = acos(z/r);
    phi = atan2(y,x);
}


void System::YLM(int l, int m, double theta, double phi, double &realYLM, double &imgYLM){

    double factor;
    double m_PLM;
    m_PLM = PLM(l,m,cos(theta));
    factor = ((2.0*double(l) + 1.0)/ (4.0*PI))*dfactorial(l,m);
    realYLM = sqrt(factor) * m_PLM * cos(double(m)*phi);
    imgYLM  = sqrt(factor) * m_PLM * sin(double(m)*phi);
}


void System::QLM(int l,int m,double theta,double phi,double &realYLM, double &imgYLM ){

    realYLM = 0.0;
    imgYLM = 0.0;
    if (m < 0) {
        YLM(l, abs(m), theta, phi, realYLM, imgYLM);
        realYLM = pow(-1.0,m)*realYLM;
        imgYLM = pow(-1.0,m)*imgYLM;
    }
    else{
        YLM(l, m, theta, phi, realYLM, imgYLM);
    }
}

void System::calculate_complexQLM_6(){

    //nn = number of neighbors
    int nn;
    double realti,imgti;
    double realYLM,imgYLM;

    // nop = parameter.nop;
    for (int ti= 0;ti<nop;ti++){

        nn = atoms[ti].n_neighbors;
        for (int mi = -6;mi < 7;mi++){

            realti = 0.0;
            imgti = 0.0;
            for (int ci = 0;ci<nn;ci++){

                QLM(6,mi,atoms[ti].n_theta[ci],atoms[ti].n_phi[ci],realYLM, imgYLM);
                realti += atoms[ti].neighborweight[ci]*realYLM;
                imgti += atoms[ti].neighborweight[ci]*imgYLM;
            }

            realti = realti/(double(nn));
            imgti = imgti/(double(nn));
            atoms[ti].realq[4][mi+6] = realti;
            atoms[ti].imgq[4][mi+6] = imgti;
        }
    }
}

//calculation of any complex qval
void System::calculate_q(vector <int> qs){

    //set_reqd_qs(qs);

    //nn = number of neighbors
    int nn;
    double realti,imgti;
    double realYLM,imgYLM;
    int q;
    double summ;

    //first make space in atoms for the number of qs needed - assign with null values
    for(int ti=0;ti<nop;ti++){
        for(int tj=0;tj<11;tj++){

            atoms[ti].q[tj] = -1;
            atoms[ti].aq[tj] = -1;
            for(int tk=0;tk<25;tk++){
                atoms[ti].realq[tj][tk] = 0;
                atoms[ti].imgq[tj][tk] = 0;
                atoms[ti].arealq[tj][tk] = 0;
                atoms[ti].aimgq[tj][tk] = 0;
            }
        }
    }

    //now check if neighbors are found
    if (!neighborsfound){
        get_all_neighbors_normal();
    }

    //note that the qvals will be in -2 pos
    //q2 will be in q0 pos and so on
    double weightsum;
    // nop = parameter.nop;
    for (int ti= 0;ti<nop;ti++){

        nn = atoms[ti].n_neighbors;
        //for(int tq=0;tq<lenqs;tq++){
        for(int tq=0;tq<qs.size();tq++){
            //find which q?
            q = qs[tq];
            //cout<<q<<endl;
            summ = 0;
            for (int mi = -q;mi < q+1;mi++){
                realti = 0.0;
                imgti = 0.0;
                weightsum = 0;
                for (int ci = 0;ci<nn;ci++){

                    QLM(q,mi,atoms[ti].n_theta[ci],atoms[ti].n_phi[ci],realYLM, imgYLM);
                    realti += atoms[ti].neighborweight[ci]*realYLM;
                    imgti += atoms[ti].neighborweight[ci]*imgYLM;
                    weightsum += atoms[ti].neighborweight[ci];
                }

            //the weights are not normalised,
            if(!voronoiused){
                realti = realti/float(weightsum);
                imgti = imgti/float(weightsum);
            }


            atoms[ti].realq[q-2][mi+q] = realti;
            atoms[ti].imgq[q-2][mi+q] = imgti;

            summ+= realti*realti + imgti*imgti;
            //summ+= realti;
            }
            //normalise summ
            summ = pow(((4.0*PI/(2*q+1)) * summ),0.5);
            atoms[ti].q[q-2] = summ;

        }

    }

    qsfound = 1;
}


//calculation of any complex aqvalb
void System::calculate_aq(vector <int> qs){

    //nn = number of neighbors
    int nn;
    double realti,imgti;
    //double realYLM,imgYLM;
    int q;
    double summ, weightsum;

    //if (!qsfound) { calculate_q(qs); }
    //note that the qvals will be in -2 pos
    //q2 will be in q0 pos and so on

    // nop = parameter.nop;
    for (int ti= 0;ti<nop;ti++){

        nn = atoms[ti].n_neighbors;

        for(int tq=0;tq<qs.size();tq++){
            //find which q?
            q = qs[tq];
            //cout<<q<<endl;
            summ = 0;
            for (int mi = 0;mi < 2*q+1;mi++){
                realti = atoms[ti].realq[q-2][mi];
                imgti = atoms[ti].imgq[q-2][mi];
                weightsum = 0;
                for (int ci = 0;ci<nn;ci++){

                    realti += atoms[atoms[ti].neighbors[ci]].realq[q-2][mi];
                    imgti += atoms[atoms[ti].neighbors[ci]].imgq[q-2][mi];
                    weightsum += atoms[ti].neighborweight[ci];
                }

            //realti = realti/(1.0+weightsum);
            //realti = realti/(1.0+weightsum);

            realti = realti/(double(nn+1));
            imgti = imgti/(double(nn+1));

            atoms[ti].arealq[q-2][mi] = realti;
            atoms[ti].aimgq[q-2][mi] = imgti;

            summ+= realti*realti + imgti*imgti;
            }
            //normalise summ
            summ = pow(((4.0*PI/(2*q+1)) * summ),0.5);
            atoms[ti].aq[q-2] = summ;

        }

    }
}

vector<double> System::gqvals(int qq){
    vector<double> qres;
    qres.reserve(nop);
    for(int i=0;i<nop;i++){
        qres.emplace_back(atoms[i].q[qq-2]);
    }

    return qres;
}

vector<double> System::gaqvals(int qq){
    vector<double> qres;
    qres.reserve(nop);
    for(int i=0;i<nop;i++){
        qres.emplace_back(atoms[i].aq[qq-2]);
    }

    return qres;
}

void System::calculate_disorder(){

    //for disorder we need sjj which is dot product with itself, self dot prouct of neighbors
    //and cross dot product
    double sumSquareti,sumSquaretj;
    double realdotproduct,imgdotproduct;
    double connection;
    double dis;

    for(int ti=0; ti<nop; ti++){

        sumSquareti = 0.0;
        realdotproduct = 0.0;
        imgdotproduct = 0.0;

        for (int mi = 0;mi < 2*solidq+1 ;mi++){
            sumSquareti += atoms[ti].realq[solidq-2][mi]*atoms[ti].realq[solidq-2][mi] + atoms[ti].imgq[solidq-2][mi] *atoms[ti].imgq[solidq-2][mi];
            realdotproduct += atoms[ti].realq[solidq-2][mi]*atoms[ti].realq[solidq-2][mi];
            imgdotproduct  += atoms[ti].imgq[solidq-2][mi] *atoms[ti].imgq[solidq-2][mi];
        }
        connection = (realdotproduct+imgdotproduct)/(sqrt(sumSquareti)*sqrt(sumSquareti));
        atoms[ti].sii = connection;

    }

    //first round is over
    //now find cross terms
    for(int ti=0; ti<nop; ti++){

        sumSquareti = 0.0;
        sumSquaretj = 0.0;
        realdotproduct = 0.0;
        imgdotproduct = 0.0;
        dis = 0;

        for(int tj=0; tj<atoms[ti].n_neighbors; tj++){
            for (int mi = 0;mi < 2*solidq+1 ;mi++){
                sumSquareti += atoms[ti].realq[solidq-2][mi]*atoms[ti].realq[solidq-2][mi] + atoms[ti].imgq[solidq-2][mi] *atoms[ti].imgq[solidq-2][mi];
                sumSquaretj += atoms[tj].realq[solidq-2][mi]*atoms[tj].realq[solidq-2][mi] + atoms[tj].imgq[solidq-2][mi] *atoms[tj].imgq[solidq-2][mi];
                realdotproduct += atoms[ti].realq[solidq-2][mi]*atoms[tj].realq[solidq-2][mi];
                imgdotproduct  += atoms[ti].imgq[solidq-2][mi] *atoms[tj].imgq[solidq-2][mi];
            }
            connection = (realdotproduct+imgdotproduct)/(sqrt(sumSquaretj)*sqrt(sumSquareti));
            dis += (atoms[ti].sii + atoms[tj].sii - 2*connection);
        }
        atoms[ti].disorder = dis/float(atoms[ti].n_neighbors);

    }
}

void System::find_average_disorder(){
    double vv;
    int nn;

    for (int ti= 0;ti<nop;ti++){
        nn = atoms[ti].n_neighbors;
        vv = atoms[ti].disorder;
        for (int ci = 0; ci<nn; ci++){
            vv += atoms[atoms[ti].neighbors[ci]].disorder;
        }
        vv = vv/(double(nn+1));
        atoms[ti].avgdisorder = vv;
    }
}
//-----------------------------------------------------
// Solids and Clustering methods
//-----------------------------------------------------
//also has to be overloaded - could be a useful function
double System::get_number_from_bond(int ti,int tj){

    double sumSquareti,sumSquaretj;
    double realdotproduct,imgdotproduct;
    double connection;
    sumSquareti = 0.0;
    sumSquaretj = 0.0;
    realdotproduct = 0.0;
    imgdotproduct = 0.0;

    for (int mi = 0;mi < 2*solidq+1 ;mi++){

        sumSquareti += atoms[ti].realq[solidq-2][mi]*atoms[ti].realq[solidq-2][mi] + atoms[ti].imgq[solidq-2][mi] *atoms[ti].imgq[solidq-2][mi];
        sumSquaretj += atoms[tj].realq[solidq-2][mi]*atoms[tj].realq[solidq-2][mi] + atoms[tj].imgq[solidq-2][mi] *atoms[tj].imgq[solidq-2][mi];
        realdotproduct += atoms[ti].realq[solidq-2][mi]*atoms[tj].realq[solidq-2][mi];
        imgdotproduct  += atoms[ti].imgq[solidq-2][mi] *atoms[tj].imgq[solidq-2][mi];
    }

    connection = (realdotproduct+imgdotproduct)/(sqrt(sumSquaretj)*sqrt(sumSquareti));
    //cout<<connection<<endl;
    return connection;
}

//overloaded version
double System::get_number_from_bond(Atom atom1,Atom atom2){

    double sumSquareti,sumSquaretj;
    double realdotproduct,imgdotproduct;
    double connection;
    sumSquareti = 0.0;
    sumSquaretj = 0.0;
    realdotproduct = 0.0;
    imgdotproduct = 0.0;

    for (int mi = 0;mi < 2*solidq+1 ;mi++){

        sumSquareti += atom1.realq[solidq-2][mi]*atom1.realq[solidq-2][mi] + atom1.imgq[solidq-2][mi] *atom1.imgq[solidq-2][mi];
        sumSquaretj += atom2.realq[solidq-2][mi]*atom2.realq[solidq-2][mi] + atom2.imgq[solidq-2][mi] *atom2.imgq[solidq-2][mi];
        realdotproduct += atom1.realq[solidq-2][mi]*atom2.realq[solidq-2][mi];
        imgdotproduct  += atom1.imgq[solidq-2][mi] *atom2.imgq[solidq-2][mi];
    }

    connection = (realdotproduct+imgdotproduct)/(sqrt(sumSquaretj)*sqrt(sumSquareti));
    return connection;
}

void System::calculate_frenkel_numbers(){

    int frenkelcons;
    double scalar;

    for (int ti= 0;ti<nop;ti++){

        frenkelcons = 0;
        atoms[ti].avq6q6 = 0.0;
        for (int c = 0;c<atoms[ti].n_neighbors;c++){

            scalar = get_number_from_bond(ti,atoms[ti].neighbors[c]);
            atoms[ti].sij[c] = scalar;
            if (scalar > threshold) frenkelcons += 1;
            atoms[ti].avq6q6 += scalar;
        }

        atoms[ti].frenkelnumber = frenkelcons;
        atoms[ti].avq6q6 /= atoms[ti].n_neighbors;

    }
}


void System::find_solid_atoms(){

    int tfrac;
    if (criteria == 0){
        for (int ti= 0;ti<nop;ti++){
          atoms[ti].issolid = ( (atoms[ti].frenkelnumber > minfrenkel) && (atoms[ti].avq6q6 > avgthreshold) );
        }
    }
    else if (criteria == 1){
        for (int ti= 0;ti<nop;ti++){
            tfrac = ((atoms[ti].frenkelnumber/double(atoms[ti].n_neighbors)) > minfrenkel);
            atoms[ti].issolid = (tfrac && (atoms[ti].avq6q6 > avgthreshold));
        }
    }

}


void System::find_clusters(double clustercutoff){
        if (clustercutoff != 0){
          for(int ti=0; ti<nop;ti++){
              atoms[ti].cutoff = clustercutoff;
          }
        }
        for(int ti=0; ti<nop;ti++){
            atoms[ti].belongsto = -1;
        }

        for (int ti= 0;ti<nop;ti++){

            if (!atoms[ti].condition) continue;

            if (atoms[ti].belongsto==-1) {atoms[ti].belongsto = atoms[ti].id; }
            for (int c = 0;c<atoms[ti].n_neighbors;c++){

                if(!atoms[atoms[ti].neighbors[c]].condition) continue;
                if(!(atoms[ti].neighbordist[atoms[ti].neighbors[c]] <= atoms[ti].cutoff)) continue;
                if (atoms[atoms[ti].neighbors[c]].belongsto==-1){
                    atoms[atoms[ti].neighbors[c]].belongsto = atoms[ti].belongsto;
                }
                else{
                    atoms[ti].belongsto = atoms[atoms[ti].neighbors[c]].belongsto;
                }
            }
        }
}

//we have to test with a recursive algorithm - to match the values that is presented
//in Grisells code.
void System::harvest_cluster(const int ti, const int clusterindex){

    int neigh;
    for(int i=0; i<atoms[ti].n_neighbors; i++){
        neigh = atoms[ti].neighbors[i];
        if(!atoms[neigh].condition) continue;
        if(!(atoms[ti].neighbordist[i] <= atoms[ti].cutoff)) continue;
        if (atoms[neigh].belongsto==-1){
            atoms[neigh].belongsto = clusterindex;
            harvest_cluster(neigh, clusterindex);
        }
    }
}

void System::find_clusters_recursive(double clustercutoff){

  if (clustercutoff != 0){
    for(int ti=0; ti<nop;ti++){
        atoms[ti].cutoff = clustercutoff;
    }
  }

    int clusterindex;
    clusterindex = 0;

    //reset belongsto indices
    for(int ti=0; ti<nop;ti++){
        atoms[ti].belongsto = -1;
    }

    for (int ti= 0;ti<nop;ti++){
        if (!atoms[ti].condition) continue;
        if (atoms[ti].belongsto==-1){
            clusterindex += 1;
            atoms[ti].belongsto = clusterindex;
            harvest_cluster(ti, clusterindex);
        }

    }
}



int System::largest_cluster(){

        int *freq = new int[nop];
        for(int ti=0;ti<nop;ti++){
            freq[ti] = 0;
        }

        for (int ti= 0;ti<nop;ti++)
        {
            if (atoms[ti].belongsto==-1) continue;
            freq[atoms[ti].belongsto-1]++;
        }

        int max=0;
        for (int ti= 0;ti<nop;ti++)
        {
            if (freq[ti]>max){
                max=freq[ti];
                maxclusterid = ti+1;
            }

        }

        get_largest_cluster_atoms();

        return max;
}

void System::get_largest_cluster_atoms(){
        for(int ti=0; ti<nop; ti++){
            atoms[ti].issurface = 1;
            atoms[ti].lcluster = 0;
            //if its in same cluster as max cluster assign it as one
            if(atoms[ti].belongsto == maxclusterid){
                atoms[ti].lcluster = 1;
            }
           //if its solid- identfy if it has liquid
            if(atoms[ti].issolid == 1){
                atoms[ti].issurface = 0;
                for(int tj=0; tj<atoms[ti].n_neighbors; tj++){
                    if(atoms[atoms[ti].neighbors[tj]].issolid == 0){
                        atoms[ti].issurface = 1;
                        break;
                    }
                }
            }
        }
}

void System::set_nucsize_parameters(double n1, double n2, double n3 ) { minfrenkel = n1; threshold = n2; avgthreshold = n3; }

int System::gsolidq() { return solidq; }
void System::ssolidq( int n) { solidq = n; }
//access function for criteria
int System::gcriteria() { return criteria; }
void System::scriteria( int n) { criteria = n; }

int System::glargestclusterid() { return maxclusterid; }
void System::slargestclusterid(int idd) { }

//-----------------------------------------------------
// Voronoi based methods
//-----------------------------------------------------
void System::salpha(int a){

    alpha = a;
}

int System::galpha(){

    return alpha ;
}

void System::set_face_cutoff(double fcut){
    face_cutoff = fcut;
}

//overloaded function; would be called
//if neighbor method voronoi is selected.
void System::get_all_neighbors_voronoi(){

    //reset voronoi flag
    voronoiused = 1;

    double d;
    double diffx,diffy,diffz;
    double r,theta,phi;
    int i;
    int ti,id,tnx,tny,tnz;

    double rx,ry,rz,tsum, fa, x, y, z, vol;
    vector<int> neigh,f_vert, vert_nos;
    vector<double> facearea, v, faceperimeters;
    voronoicell_neighbor c;
    vector< vector<double> > nweights;
    vector< vector<int> > nneighs;
    vector<int> idss;
    //vector<int> nvector;
    double weightsum;


    if (!fileread) { read_particle_file(inputfile); }

    //pre_container pcon(boxdims[0][0],boxdims[1][1],boxdims[1][0],boxdims[1][1],boxdims[2][0],boxdims[2][1],true,true,true);
    pre_container pcon(0.00, boxx, 0.00, boxy, 0.0, boxz, true, true, true);
    for(int i=0; i<nop; i++){
        pcon.put(i, atoms[i].posx-boxdims[0][0], atoms[i].posy-boxdims[1][0], atoms[i].posz-boxdims[2][0]);
    }
    pcon.guess_optimal(tnx,tny,tnz);
    //container con(boxdims[0][0],boxdims[1][1],boxdims[1][0],boxdims[1][1],boxdims[2][0],boxdims[2][1],tnx,tny,tnz,true,true,true, nop);
    container con(0.00, boxx, 0.00, boxy, 0.0, boxz, tnx, tny, tnz, true, true, true, nop);
    pcon.setup(con);

    c_loop_all cl(con);
    if (cl.start()) do if(con.compute_cell(c,cl)) {
            ti=cl.pid();
            c.face_areas(facearea);
            c.neighbors(neigh);
            c.face_orders(f_vert);
            c.face_vertices(vert_nos);
            c.vertices(x,y,z,v);
            c.face_perimeters(faceperimeters);

            vol = c.volume();
            tsum = 0;
            vector <double> dummyweights;
            vector <int> dummyneighs;

            //only loop over neighbors
            weightsum = 0.0;
            for (int i=0; i<facearea.size(); i++){
                weightsum += pow(facearea[i], alpha);
            }


            //assign to nvector
            atoms[ti].volume = vol;
            atoms[ti].vertex_vectors = v;
            atoms[ti].vertex_numbers = vert_nos;
            atoms[ti].cutoff = cbrt(3*vol/(4*3.141592653589793));
            //assign to the atom
            //atoms[ti].vorovector = nvector;

            //only loop over neighbors
            //weightsum = 0.0;
            //for (int i=0; i<facearea.size(); i++){
            //    weightsum += facearea[i];
            //}
            for (int tj=0; tj<neigh.size(); tj++){

                //if filter doesnt work continue
                if ((filter == 1) && (atoms[ti].type != atoms[tj].type)){
                    continue;
                }

                atoms[ti].neighbors[tj] = neigh[tj];
                atoms[ti].n_neighbors += 1;
                d = get_abs_distance(ti,neigh[tj],diffx,diffy,diffz);
                atoms[ti].neighbordist[tj] = d;
                //weight is set to 1.0, unless manually reset
                atoms[ti].neighborweight[tj] = pow(facearea[tj], alpha)/weightsum;
                atoms[ti].facevertices[tj] = f_vert[tj];
                atoms[ti].faceperimeters[tj] = faceperimeters[tj];
                atoms[ti].n_diffx[tj] = diffx;
                atoms[ti].n_diffy[tj] = diffy;
                atoms[ti].n_diffz[tj] = diffz;
                convert_to_spherical_coordinates(diffx, diffy, diffz, r, phi, theta);
                atoms[ti].n_r[tj] = r;
                atoms[ti].n_phi[tj] = phi;
                atoms[ti].n_theta[tj] = theta;

            }

    } while (cl.inc());

    //mark end of neighbor calc
    neighborsfound = 1;

    //now calculate the averged volume
    find_average_volume();


}


void System::find_average_volume(){
    double vv;
    int nn;

    for (int ti= 0;ti<nop;ti++){
        nn = atoms[ti].n_neighbors;
        vv = atoms[ti].volume;
        for (int ci = 0;ci<nn;ci++){
            vv += atoms[atoms[ti].neighbors[ci]].volume;
        }
        vv = vv/(double(nn+1));
        atoms[ti].avgvolume = vv;
    }
}






