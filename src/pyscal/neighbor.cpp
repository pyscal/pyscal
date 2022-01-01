#include "modsystem.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include "string.h"
#include <chrono>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <map>
#include <string>
#include <any>

double get_abs_distance(vector<double> pos1, vector<double> pos2, 
	const int& triclinic, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
	const vector<double>& box,
    double& diffx,
    double& diffy,
    double& diffz){
    /*
    Get absolute distance between two atoms
    */

    double abs, ax, ay, az;
    diffx = pos1[0] - pos2[0];
    diffy = pos1[1] - pos2[1];
    diffz = pos1[2] - pos2[2];


    if (triclinic == 1){

        //convert to the triclinic system
        ax = rotinv[0][0]*diffx + rotinv[0][1]*diffy + rotinv[0][2]*diffz;
        ay = rotinv[1][0]*diffx + rotinv[1][1]*diffy + rotinv[1][2]*diffz;
        az = rotinv[2][0]*diffx + rotinv[2][1]*diffy + rotinv[2][2]*diffz;

        //scale to match the triclinic box size
        diffx = ax*box[0];
        diffy = ay*box[1];
        diffz = az*box[2];

        //now check pbc
        //nearest image
        if (diffx> box[0]/2.0) {diffx-=box[0];};
        if (diffx<-box[0]/2.0) {diffx+=box[0];};
        if (diffy> box[1]/2.0) {diffy-=box[1];};
        if (diffy<-box[1]/2.0) {diffy+=box[1];};
        if (diffz> box[2]/2.0) {diffz-=box[2];};
        if (diffz<-box[2]/2.0) {diffz+=box[2];};

        //now divide by box vals - scale down the size
        diffx = diffx/box[0];
        diffy = diffy/box[1];
        diffz = diffz/box[2];

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
        if (diffx> box[0]/2.0) {diffx-=box[0];};
        if (diffx<-box[0]/2.0) {diffx+=box[0];};
        if (diffy> box[1]/2.0) {diffy-=box[1];};
        if (diffy<-box[1]/2.0) {diffy+=box[1];};
        if (diffz> box[2]/2.0) {diffz-=box[2];};
        if (diffz<-box[2]/2.0) {diffz+=box[2];};
        abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
    }
    return abs;
}

void reset_all_neighbors(py::dict& atoms){
    /*
    Reset all neighbors
    */   
    vector<vector<int>> temp2int;
    vector<double> tempdouble;
    vector<vector<double>> temp2double;
    vector<vector<vector<double>>> temp3double;

    for(unsigned int ti=0; ti<atoms.size(); ti++){
        atoms[py::str("neighbors")] = temp2int;
        atoms[py::str("neighbordist")] = temp2double;
        atoms[py::str("temp_neighbors")] = temp2int;
        atoms[py::str("temp_neighbordist")] = temp2double;
        atoms[py::str("neighborweight")] = temp2double;
        atoms[py::str("diff")] = temp3double;
        atoms[py::str("r")] = temp2double;
        atoms[py::str("phi")] = temp2double;
        atoms[py::str("theta")] = temp2double;
        atoms[py::str("cutoff")] = tempdouble;
    }
}

void convert_to_spherical_coordinates(double x, 
    double y, 
    double z, 
    double &r, 
    double &phi, 
    double &theta){

    r = sqrt(x*x+y*y+z*z);
    theta = acos(z/r);
    phi = atan2(y,x);
}

void get_all_neighbors_normal(py::dict& atoms,
    const double neighbordistance,
    const int triclinic,
    const int filter, 
    const vector<vector<double>> rot, 
    const vector<vector<double>> rotinv,
    const vector<double> box,
    const int filter)
    {
    
    double d;
    double diffx,diffy,diffz;
    double tempr,temptheta,tempphi;
    vector<double> diffi, diffj;

    //access positions and put it in an array
    vector<vector<double>> positions = atoms[py::str("positions")].cast<vector<vector<double>>>();
    vector<int> types = atoms[py::str("types")].cast<vector<int>>();
    //auto positions = atoms[py::str("positions")].cast<py::array_t<double>>();

    int nop = positions.size();
    vector<vector<int>> neighbors(nop);
    vector<vector<double>> neighbordist(nop);
    vector<vector<double>> neighborweight(nop);
    vector<vector<vector<double>>> diff(nop);
    vector<vector<double>> r(nop);
    vector<vector<double>> phi(nop);
    vector<vector<double>> theta(nop);
    vector<double> cutoff(nop); 

    //now loop and calculate things
    for (int ti=0; ti<nop; ti++){
        for (int tj=ti+1; tj<nop; tj++){
            d = get_abs_distance(positions[ti], positions[tj],
                triclinic, rot, rotinv, box, 
                diffx, diffy, diffz);
            if (d < neighbordistance){
                //if ((filter == 1) && (atoms[ti].type != atoms[tj].type)){
                //    continue;
                //}
                //else if ((filter == 2) && (atoms[ti].type == atoms[tj].type)){
                 //   continue;
                //}
                
                neighbors[ti].emplace_back(tj);
                neighbors[tj].emplace_back(ti);

                neighbordist[ti].emplace_back(d);
                neighbordist[tj].emplace_back(d);

                neighborweight[ti].emplace_back(1.00);
                neighborweight[tj].emplace_back(1.00);

                diffi.clear();
                diffi.emplace_back(diffx);
                diffi.emplace_back(diffy);
                diffi.emplace_back(diffz);

                diffj.clear();
                diffj.emplace_back(-diffx);
                diffj.emplace_back(-diffy);
                diffj.emplace_back(-diffz);

                diff[ti].emplace_back(diffi);
                diff[tj].emplace_back(diffj);
                
                convert_to_spherical_coordinates(diffx, diffy, diffz, tempr, tempphi, temptheta);
                
                r[ti].emplace_back(tempr);
                phi[ti].emplace_back(tempphi);
                theta[ti].emplace_back(temptheta);
                //n_neighbors += 1;
                cutoff[ti] = neighbordistance;

                convert_to_spherical_coordinates(-diffx, -diffy, -diffz, tempr, tempphi, temptheta);

                r[tj].emplace_back(tempr);
                phi[tj].emplace_back(tempphi);
                theta[tj].emplace_back(temptheta);
                //atoms[tj].n_neighbors += 1;
                cutoff[tj] = neighbordistance;

            }
        }
    }

    //calculation over lets assign
    atoms[py::str("neighbors")] = neighbors;
    atoms[py::str("neighbordist")] = neighbordist;
    atoms[py::str("neighborweight")] = neighborweight;
    atoms[py::str("diff")] = diff;
    atoms[py::str("r")] = r;
    atoms[py::str("theta")] = theta;
    atoms[py::str("phi")] = phi;
    atoms[py::str("cutoff")] = cutoff;
}


int cell_index(int cx, int cy, int cz, int nx, int ny, int nz){
    return cx*ny*nz + cy*nz + cz;
}

vector<int> cell_periodic(int i, int j, int k, int nx, int ny, int nz){
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

vector<cell> set_up_cells(const vector<vector<double>>& positions,
    const vector<double>& box,
    const double neighbordistance){

      int maincell, subcell, total_cells;
      int nx, ny, nz;
      double lx, ly, lz;
      vector<int> cc;
      vector<cell> cells;
      
      //find of all find the number of cells in each direction
      nx = box[0]/neighbordistance;
      ny = box[1]/neighbordistance;
      nz = box[2]/neighbordistance;
      
      //now use this to find length of cell in each direction
      lx = box[0]/nx;
      ly = box[1]/ny;
      lz = box[2]/nz;

      //find the total number of cells
      total_cells = nx*ny*nz;
      cells.resize(total_cells);
      
      //all neighbor cells are also added
      for(int i=0; i<nx; i++){
         for(int j=0; j<ny; j++){
           for(int k=0; k<nz; k++){
              maincell = cell_index(i, j, k, nx, ny, nz);
              for(int si=i-1; si<=i+1; si++){
                  for(int sj=j-1; sj<=j+1; sj++){
                      for(int sk=k-1; sk<=k+1; sk++){
                         cc = cell_periodic(si, sj, sk, nx, ny, nz);
                         subcell = cell_index(cc[0], cc[1], cc[2], nx, ny, nz);
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
      int nop = positions.size();

      //now loop over all atoms and assign cells
      for(int ti=0; ti<nop; ti++){

          //calculate c indices for the atom
          dx = positions[ti][0];
          dy = positions[ti][1];
          dz = positions[ti][2];
          
          //now apply boxdims
          if( abs(dx-0) < 1E-6)
              dx = 0;
          if( abs(dy-0) < 1E-6)
              dy = 0;
          if( abs(dz-0) < 1E-6)
              dz = 0;

          if (dx < 0) dx+=box[0];
          else if (dx >= box[0]) dx-=box[0];
          if (dy < 0) dy+=box[1];
          else if (dy >= box[1]) dy-=box[1];
          if (dz < 0) dz+=box[2];
          else if (dz >= box[2]) dz-=box[2];

          //now find c vals
          cx = dx/lx;
          cy = dy/ly;
          cz = dz/lz;

          //now get cell index
          ind = cell_index(cx, cy, cz, nx, ny, nz);
          
          //now add the atom to the corresponding cells
          cells[ind].members.emplace_back(ti);

      }
      return cells;
}


void get_all_neighbors_cells(py::dict& atoms,
    const double& neighbordistance,
    const int& triclinic,
    const int& filter, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box,
    const int filter){

    double d;
    double diffx,diffy,diffz;
    double tempr,temptheta,tempphi;
    vector<double> diffi, diffj;
    int ti, tj;

    //access positions and put it in an array
    vector<vector<double>> positions = atoms[py::str("positions")].cast<vector<vector<double>>>();

    int nop = positions.size();
    vector<vector<int>> neighbors(nop);
    vector<vector<double>> neighbordist(nop);
    vector<vector<double>> neighborweight(nop);
    vector<vector<vector<double>>> diff(nop);
    vector<vector<double>> r(nop);
    vector<vector<double>> phi(nop);
    vector<vector<double>> theta(nop);
    vector<double> cutoff(nop); 
    
    vector<cell> cells = set_up_cells(positions, box, neighbordistance);
    int total_cells = cells.size();
    int subcell;

    //now loop to find distance
    for(int i=0; i<total_cells; i++){
        //for each member in cell i
        for(size_t mi=0; mi<cells[i].members.size(); mi++){
            //now go through the neighbors
            ti = cells[i].members[mi];
            for(size_t j=0 ; j<cells[i].neighbor_cells.size(); j++){
               //loop through members of j
               subcell = cells[i].neighbor_cells[j];
               for(size_t mj=0; mj<cells[subcell].members.size(); mj++){
                    //now we have mj -> members/compare with
                    tj = cells[subcell].members[mj];
                    //compare ti and tj and add
                    if (ti < tj){
                        d = get_abs_distance(positions[ti], positions[tj],
                            triclinic, rot, rotinv, box, 
                            diffx, diffy, diffz);
                        if (d < neighbordistance){
                            //if ((filter == 1) && (atoms[ti].type != atoms[tj].type)){
                            //    continue;
                            //}
                            //else if ((filter == 2) && (atoms[ti].type == atoms[tj].type)){
                            //    continue;
                            //}
                            neighbors[ti].emplace_back(tj);
                            neighbors[tj].emplace_back(ti);

                            neighbordist[ti].emplace_back(d);
                            neighbordist[tj].emplace_back(d);

                            neighborweight[ti].emplace_back(1.00);
                            neighborweight[tj].emplace_back(1.00);

                            diffi.clear();
                            diffi.emplace_back(diffx);
                            diffi.emplace_back(diffy);
                            diffi.emplace_back(diffz);

                            diffj.clear();
                            diffj.emplace_back(-diffx);
                            diffj.emplace_back(-diffy);
                            diffj.emplace_back(-diffz);

                            diff[ti].emplace_back(diffi);
                            diff[tj].emplace_back(diffj);
                            
                            convert_to_spherical_coordinates(diffx, diffy, diffz, tempr, tempphi, temptheta);
                            
                            r[ti].emplace_back(tempr);
                            phi[ti].emplace_back(tempphi);
                            theta[ti].emplace_back(temptheta);
                            //n_neighbors += 1;
                            cutoff[ti] = neighbordistance;

                            convert_to_spherical_coordinates(-diffx, -diffy, -diffz, tempr, tempphi, temptheta);

                            r[tj].emplace_back(tempr);
                            phi[tj].emplace_back(tempphi);
                            theta[tj].emplace_back(temptheta);
                            //atoms[tj].n_neighbors += 1;
                            cutoff[tj] = neighbordistance;
                        }
                    }
                }
            }
        }
    }

    //calculation over lets assign
    atoms[py::str("neighbors")] = neighbors;
    atoms[py::str("neighbordist")] = neighbordist;
    atoms[py::str("neighborweight")] = neighborweight;
    atoms[py::str("diff")] = diff;
    atoms[py::str("r")] = r;
    atoms[py::str("theta")] = theta;
    atoms[py::str("phi")] = phi;
    atoms[py::str("cutoff")] = cutoff;
}


void get_temp_neighbors_brute(const vector<vector<double>>& positions,
    vector<vector<datom>>& temp_neighbors,
    const int& triclinic,
    const double neighbordistance, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box)
{
    double d;
    double diffx,diffy,diffz;
    int nop = positions.size();
    temp_neighbors.resize(nop);

    for (int ti=0; ti<nop; ti++){
        for (int tj=ti+1; tj<nop; tj++){
            d = get_abs_distance(positions[ti], positions[tj],
                triclinic, rot, rotinv, box, 
                diffx, diffy, diffz);
            if (d <= neighbordistance){
                datom x = {d, tj};
                temp_neighbors[ti].emplace_back(x);
                datom y = {d, ti};
                temp_neighbors[tj].emplace_back(y);
            }
        }
    }
}


void get_temp_neighbors_cells(const vector<vector<double>>& positions,
    vector<vector<datom>>& temp_neighbors,
    const int& triclinic,
    const double neighbordistance, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box){

    //first create cells
    vector<cell> cells = set_up_cells(positions, box, neighbordistance);
    int total_cells = cells.size();
    int subcell;
    int ti, tj;
    double d;
    double diffx,diffy,diffz;
    int nop = positions.size();
    temp_neighbors.resize(nop);

    //now loop to find distance
    for(int i=0; i<total_cells; i++){
        //now go over the neighbor cells
        //for each member in cell i
        for(size_t mi=0; mi<cells[i].members.size(); mi++){
            //now go through the neighbors
            ti = cells[i].members[mi];
            for(size_t j=0 ; j<cells[i].neighbor_cells.size(); j++){
               //loop through members of j
               subcell = cells[i].neighbor_cells[j];
               for(size_t mj=0; mj<cells[subcell].members.size(); mj++){
                    //now we have mj -> members/compare with
                    tj = cells[subcell].members[mj];
                    //compare ti and tj and add
                    if (ti < tj){
                        d = get_abs_distance(positions[ti], positions[tj],
                            triclinic, rot, rotinv, box, 
                            diffx, diffy, diffz);
                        if (d <= neighbordistance){
                            datom x = {d, tj};
                            temp_neighbors[ti].emplace_back(x);
                            datom y = {d, ti};
                            temp_neighbors[tj].emplace_back(y);
                      }
                  }
               }
            }
        }
    }
}


int get_all_neighbors_bynumber(py::dict& atoms,
    double& neighbordistance,
    const int& triclinic,
    const int& filter, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box,
    double prefactor,
    int nns, 
    int usecells,
    int assign){

    double d;
    double diffx,diffy,diffz;
    double tempr,temptheta,tempphi;
    vector<double> diffi, diffj;
    int finished;
    finished = 1;

    //access positions and put it in an array
    vector<vector<double>> positions = atoms[py::str("positions")].cast<vector<vector<double>>>();

    int nop = positions.size();
    vector<vector<int>> neighbors(nop);
    vector<vector<double>> neighbordist(nop);
    vector<vector<double>> neighborweight(nop);
    vector<vector<vector<double>>> diff(nop);
    vector<vector<double>> r(nop);
    vector<vector<double>> phi(nop);
    vector<vector<double>> theta(nop);
    vector<double> cutoff(nop);     
    vector<vector<int>> atom_temp_neighbors(nop);
    vector<vector<double>> atom_temp_neighbordist(nop);

    vector<int> nids;
    vector<double> dists, sorted_dists;
    double boxvol;
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
        boxvol = box[0]*box[1]*box[2];
    }
    //now find the volume per particle
    double guessvol = boxvol/float(nop);

    //guess the side of a cube that is occupied by an atom - this is a guess distance
    double guessdist = cbrt(guessvol);

    //now add some safe padding - this is the prefactor which we will read in
    guessdist = prefactor*guessdist;
    neighbordistance = guessdist;
    vector<vector<datom>> temp_neighbors;
    if (usecells){
        get_temp_neighbors_cells(positions, temp_neighbors, triclinic, neighbordistance, rot, rotinv, box);
    }
    else{
        get_temp_neighbors_brute(positions, temp_neighbors, triclinic, neighbordistance, rot, rotinv, box);
    }
    for (int ti=0; ti<nop; ti++){
        if (int(temp_neighbors[ti].size()) < nns){
            return 0;
        }

        sort(temp_neighbors[ti].begin(), temp_neighbors[ti].end(), by_dist());

        for(size_t i=0; i<temp_neighbors[ti].size(); i++){
            atom_temp_neighbors[ti].emplace_back(temp_neighbors[ti][i].index);
            atom_temp_neighbordist[ti].emplace_back(temp_neighbors[ti][i].dist);
        }

        if(assign == 1){
            //assign the neighbors
            for(int i=0; i<nns; i++){
                int tj = temp_neighbors[ti][i].index;
                d = get_abs_distance(positions[ti], positions[tj],
                    triclinic, rot, rotinv, box, 
                    diffx, diffy, diffz);

                neighbors[ti].emplace_back(tj);
                neighbordist[ti].emplace_back(d);
                neighborweight[ti].emplace_back(1.00);
                
                diffi.clear();
                diffi.emplace_back(diffx);
                diffi.emplace_back(diffy);
                diffi.emplace_back(diffz);
                diff[ti].emplace_back(diffi);
                
                convert_to_spherical_coordinates(diffx, diffy, diffz, tempr, tempphi, temptheta);
                
                r[ti].emplace_back(tempr);
                phi[ti].emplace_back(tempphi);
                theta[ti].emplace_back(temptheta);
                //n_neighbors += 1;
                cutoff[ti] = neighbordistance;
            }
        }

        finished = 1;            
    }
    if (assign==1){
        atoms[py::str("neighbors")] = neighbors;
        atoms[py::str("neighbordist")] = neighbordist;
        atoms[py::str("neighborweight")] = neighborweight;
        atoms[py::str("diff")] = diff;
        atoms[py::str("r")] = r;
        atoms[py::str("theta")] = theta;
        atoms[py::str("phi")] = phi;
        atoms[py::str("cutoff")] = cutoff;
    }
    //in any case add the temp neighbors
    atoms[py::str("temp_neighbors")] = atom_temp_neighbors;
    atoms[py::str("temp_neighbordist")] = atom_temp_neighbordist;

    return finished;
}


int get_all_neighbors_sann(py::dict& atoms,
    double& neighbordistance,
    const int& triclinic,
    const int& filter, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box,
    double prefactor,
    int usecells){

    double d, dcut;
    double diffx,diffy,diffz;
    double tempr,temptheta,tempphi;
    vector<double> diffi, diffj;
    int m, maxneighs, finished;
    finished = 1;

    //access positions and put it in an array
    vector<vector<double>> positions = atoms[py::str("positions")].cast<vector<vector<double>>>();

    int nop = positions.size();
    vector<vector<int>> neighbors(nop);
    vector<vector<double>> neighbordist(nop);
    vector<vector<double>> neighborweight(nop);
    vector<vector<vector<double>>> diff(nop);
    vector<vector<double>> r(nop);
    vector<vector<double>> phi(nop);
    vector<vector<double>> theta(nop);
    vector<double> cutoff(nop);     
    vector<vector<int>> atom_temp_neighbors(nop);
    vector<vector<double>> atom_temp_neighbordist(nop);

    vector<int> nids;
    vector<double> dists, sorted_dists;
    double summ;
    double boxvol;

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
        boxvol = box[0]*box[1]*box[2];
    }

    //now find the volume per particle
    double guessvol = boxvol/float(nop);

    //guess the side of a cube that is occupied by an atom - this is a guess distance
    double guessdist = cbrt(guessvol);

    //now add some safe padding - this is the prefactor which we will read in
    guessdist = prefactor*guessdist;
    neighbordistance = guessdist;
    vector<vector<datom>> temp_neighbors;

    if (usecells){
        get_temp_neighbors_cells(positions, temp_neighbors, triclinic, neighbordistance, rot, rotinv, box);
    }
    else{
        get_temp_neighbors_brute(positions, temp_neighbors, triclinic, neighbordistance, rot, rotinv, box);
    }
    
    for (int ti=0; ti<nop; ti++){
        if (temp_neighbors[ti].size() < 3){
            return 0;
        }

        sort(temp_neighbors[ti].begin(), temp_neighbors[ti].end(), by_dist());

        for(size_t i=0; i<temp_neighbors[ti].size(); i++){
            atom_temp_neighbors[ti].emplace_back(temp_neighbors[ti][i].index);
            atom_temp_neighbordist[ti].emplace_back(temp_neighbors[ti][i].dist);
        }

        m = 3;
        summ = 0;

        for(int i=0; i<m; i++){
            summ += temp_neighbors[ti][i].dist;
            int tj = temp_neighbors[ti][i].index;
            d = get_abs_distance(positions[ti], positions[tj],
                triclinic, rot, rotinv, box, 
                diffx, diffy, diffz);

            neighbors[ti].emplace_back(tj);
            neighbordist[ti].emplace_back(d);
            neighborweight[ti].emplace_back(1.00);
            
            diffi.clear();
            diffi.emplace_back(diffx);
            diffi.emplace_back(diffy);
            diffi.emplace_back(diffz);
            diff[ti].emplace_back(diffi);
            
            convert_to_spherical_coordinates(diffx, diffy, diffz, tempr, tempphi, temptheta);
            
            r[ti].emplace_back(tempr);
            phi[ti].emplace_back(tempphi);
            theta[ti].emplace_back(temptheta);
            //n_neighbors += 1;
        }

        dcut = summ/float(m-2);
        maxneighs = temp_neighbors[ti].size();

        while( (m < maxneighs) && (dcut >= temp_neighbors[ti][m].dist)){
            //increase m
            m = m+1;

            //here now we can add this to the list neighbors and process things
            int tj = temp_neighbors[ti][m].index;
            d = get_abs_distance(positions[ti], positions[tj],
                triclinic, rot, rotinv, box, 
                diffx, diffy, diffz);

            neighbors[ti].emplace_back(tj);
            neighbordist[ti].emplace_back(d);
            neighborweight[ti].emplace_back(1.00);
            
            diffi.clear();
            diffi.emplace_back(diffx);
            diffi.emplace_back(diffy);
            diffi.emplace_back(diffz);
            diff[ti].emplace_back(diffi);
            
            convert_to_spherical_coordinates(diffx, diffy, diffz, tempr, tempphi, temptheta);
            
            r[ti].emplace_back(tempr);
            phi[ti].emplace_back(tempphi);
            theta[ti].emplace_back(temptheta);
            //n_neighbors += 1;        

            //find new dcut
            summ = summ + temp_neighbors[ti][m].dist;
            dcut = summ/float(m-2);
            cutoff[ti] = dcut;
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
    atoms[py::str("neighbors")] = neighbors;
    atoms[py::str("neighbordist")] = neighbordist;
    atoms[py::str("neighborweight")] = neighborweight;
    atoms[py::str("diff")] = diff;
    atoms[py::str("r")] = r;
    atoms[py::str("theta")] = theta;
    atoms[py::str("phi")] = phi;
    atoms[py::str("cutoff")] = cutoff;

    //in any case add the temp neighbors
    atoms[py::str("temp_neighbors")] = atom_temp_neighbors;
    atoms[py::str("temp_neighbordist")] = atom_temp_neighbordist;

    return finished;
}


int get_all_neighbors_adaptive(py::dict& atoms,
    double& neighbordistance,
    const int& triclinic,
    const int& filter, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box,
    double prefactor,
    int nlimit,
    double padding, 
    int usecells){

    double d, dcut;
    double diffx,diffy,diffz;
    double tempr,temptheta,tempphi;
    vector<double> diffi, diffj;
    int finished = 1;

    //access positions and put it in an array
    vector<vector<double>> positions = atoms[py::str("positions")].cast<vector<vector<double>>>();

    int nop = positions.size();
    vector<vector<int>> neighbors(nop);
    vector<vector<double>> neighbordist(nop);
    vector<vector<double>> neighborweight(nop);
    vector<vector<vector<double>>> diff(nop);
    vector<vector<double>> r(nop);
    vector<vector<double>> phi(nop);
    vector<vector<double>> theta(nop);
    vector<double> cutoff(nop);     
    vector<vector<int>> atom_temp_neighbors(nop);
    vector<vector<double>> atom_temp_neighbordist(nop);

    vector<int> nids;
    vector<double> dists, sorted_dists;
    double summ;
    double boxvol;

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
        boxvol = box[0]*box[1]*box[2];
    }

    //now find the volume per particle
    double guessvol = boxvol/float(nop);

    //guess the side of a cube that is occupied by an atom - this is a guess distance
    double guessdist = cbrt(guessvol);

    //now add some safe padding - this is the prefactor which we will read in
    guessdist = prefactor*guessdist;
    neighbordistance = guessdist;
    vector<vector<datom>> temp_neighbors;

    if (usecells){
        get_temp_neighbors_cells(positions, temp_neighbors, triclinic, neighbordistance, rot, rotinv, box);
    }
    else{
        get_temp_neighbors_brute(positions, temp_neighbors, triclinic, neighbordistance, rot, rotinv, box);
    }
    
    for (int ti=0; ti<nop; ti++){
        if (int(temp_neighbors[ti].size()) < nlimit){
            return 0;
        }

        sort(temp_neighbors[ti].begin(), temp_neighbors[ti].end(), by_dist());

        for(size_t i=0; i<temp_neighbors[ti].size(); i++){
            atom_temp_neighbors[ti].emplace_back(temp_neighbors[ti][i].index);
            atom_temp_neighbordist[ti].emplace_back(temp_neighbors[ti][i].dist);
        }

        summ = 0;
        for(int i=0; i<nlimit; i++){
            summ += temp_neighbors[ti][i].dist;
        }
        dcut = padding*(1.0/float(nlimit))*summ;

        for(size_t j=0; j<temp_neighbors[ti].size(); j++){
            int tj = temp_neighbors[ti][j].index;
            if (temp_neighbors[ti][j].dist < dcut){

                //if ((filter == 1) && (atoms[ti].type != atoms[tj].type)){
                //    continue;
                //}
                //else if ((filter == 2) && (atoms[ti].type == atoms[tj].type)){
                //    continue;
                //}
                d = get_abs_distance(positions[ti], positions[tj],
                    triclinic, rot, rotinv, box, 
                    diffx, diffy, diffz);
                neighbors[ti].emplace_back(tj);
                neighbordist[ti].emplace_back(d);
                neighborweight[ti].emplace_back(1.00);
                
                diffi.clear();
                diffi.emplace_back(diffx);
                diffi.emplace_back(diffy);
                diffi.emplace_back(diffz);
                diff[ti].emplace_back(diffi);
                
                convert_to_spherical_coordinates(diffx, diffy, diffz, tempr, tempphi, temptheta);
                
                r[ti].emplace_back(tempr);
                phi[ti].emplace_back(tempphi);
                theta[ti].emplace_back(temptheta);
                //n_neighbors += 1;        

            }
        }        

    }
    atoms[py::str("neighbors")] = neighbors;
    atoms[py::str("neighbordist")] = neighbordist;
    atoms[py::str("neighborweight")] = neighborweight;
    atoms[py::str("diff")] = diff;
    atoms[py::str("r")] = r;
    atoms[py::str("theta")] = theta;
    atoms[py::str("phi")] = phi;
    atoms[py::str("cutoff")] = cutoff;

    //in any case add the temp neighbors
    atoms[py::str("temp_neighbors")] = atom_temp_neighbors;
    atoms[py::str("temp_neighbordist")] = atom_temp_neighbordist;

    return 1;
}
