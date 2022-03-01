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

void get_cna_neighbors(py::dict& atoms,
    const int& triclinic,
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box,
	double lattice_constant,
	int style){

    double d;
    double diffx,diffy,diffz;
    double tempr,temptheta,tempphi;
    vector<double> diffi, diffj;

    vector<vector<double>> positions = atoms[py::str("positions")].cast<vector<vector<double>>>();
    vector<vector<int>> atom_temp_neighbors = atoms[py::str("temp_neighbors")].cast<vector<vector<int>>>();;
    vector<vector<double>> atom_temp_neighbordist = atoms[py::str("temp_neighbordist")].cast<vector<vector<double>>>();;

    int nop = positions.size();
    vector<vector<int>> neighbors(nop);
    vector<vector<double>> neighbordist(nop);
    vector<vector<double>> neighborweight(nop);
    vector<vector<vector<double>>> diff(nop);
    vector<vector<double>> r(nop);
    vector<vector<double>> phi(nop);
    vector<vector<double>> theta(nop);
    vector<double> cutoff(nop);     

    double factor, ncount;

    if (style == 1){
        factor = 0.854;
        ncount = 12;
    }
    else{
        factor = 1.207;
        ncount = 14;
    }

    for (int ti=0; ti<nop; ti++){
        cutoff[ti] = factor*lattice_constant;
        for(int i=0 ; i<ncount; i++){
            int tj = atom_temp_neighbors[ti][i];
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
}

void get_acna_neighbors_cn12(py::dict& atoms,
    const int& triclinic,
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box){

    double d;
    double diffx,diffy,diffz;
    double tempr,temptheta,tempphi;
    vector<double> diffi, diffj;

    vector<vector<double>> positions = atoms[py::str("positions")].cast<vector<vector<double>>>();
    vector<vector<int>> atom_temp_neighbors = atoms[py::str("temp_neighbors")].cast<vector<vector<int>>>();;
    vector<vector<double>> atom_temp_neighbordist = atoms[py::str("temp_neighbordist")].cast<vector<vector<double>>>();;

    int nop = positions.size();
    vector<vector<int>> neighbors(nop);
    vector<vector<double>> neighbordist(nop);
    vector<vector<double>> neighborweight(nop);
    vector<vector<vector<double>>> diff(nop);
    vector<vector<double>> r(nop);
    vector<vector<double>> phi(nop);
    vector<vector<double>> theta(nop);
    vector<double> cutoff(nop);

    for (int ti=0; ti<nop; ti++){
        if (atom_temp_neighbors[ti].size() > 11){
            double ssum = 0;
            for(int i=0 ; i<12; i++){
                ssum += atom_temp_neighbordist[ti][i];
                int tj = atom_temp_neighbors[ti][i];
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
            }
            cutoff[ti] = 1.207*ssum/12.00;
        }
    }     

}


void get_acna_neighbors_cn14(py::dict& atoms,
    const int& triclinic,
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box){

    double d;
    double diffx,diffy,diffz;
    double tempr,temptheta,tempphi;
    vector<double> diffi, diffj;

    vector<vector<double>> positions = atoms[py::str("positions")].cast<vector<vector<double>>>();
    vector<vector<int>> atom_temp_neighbors = atoms[py::str("temp_neighbors")].cast<vector<vector<int>>>();;
    vector<vector<double>> atom_temp_neighbordist = atoms[py::str("temp_neighbordist")].cast<vector<vector<double>>>();;

    int nop = positions.size();
    vector<vector<int>> neighbors(nop);
    vector<vector<double>> neighbordist(nop);
    vector<vector<double>> neighborweight(nop);
    vector<vector<vector<double>>> diff(nop);
    vector<vector<double>> r(nop);
    vector<vector<double>> phi(nop);
    vector<vector<double>> theta(nop);
    vector<double> cutoff(nop);

    for (int ti=0; ti<nop; ti++){
        if (atom_temp_neighbors[ti].size() > 11){
            double ssum = 0;
            for(int i=0 ; i<8; i++){
                ssum += 1.1547*atom_temp_neighbordist[ti][i];
            }
            for(int i=8 ; i<14; i++){
                ssum += atom_temp_neighbordist[ti][i];
            }
            for(int i=0 ; i<14; i++){
                int tj = atom_temp_neighbors[ti][i];
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
            }
            cutoff[ti] = 1.207*ssum/14.00;
        }
    }     
}


void get_common_neighbors(const int& triclinic,
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box,
    const int ti,
    const vector<vector<double>>& positions,
    const vector<double>& cutoff,
    const vector<vector<int>>& neighbors,
    vector<vector<vector<int>>>& cna,
    vector<vector<vector<int>>>& common){
    
    int m, n;
    double d, diffx, diffy, diffz;
   
    cna[ti].clear();
    cna[ti].resize(neighbors[ti].size());
    common[ti].clear();
    common[ti].resize(neighbors[ti].size());

    for(int i=0; i<neighbors[ti].size(); i++){
        for(int j=0; j<4; j++){
            cna[ti][i].emplace_back(0);
        }
    }
    
    for(int i=0; i<neighbors[ti].size()-1; i++){
        m = neighbors[ti][i];
        for(int j=i+1; j<neighbors[ti].size(); j++){
            n = neighbors[ti][j];
            d = get_abs_distance(positions[m], positions[n],
                triclinic, rot, rotinv, box, 
                diffx, diffy, diffz);

            if (d <= cutoff[ti]){
                cna[ti][i][0]++;
                common[ti][i].emplace_back(n);
                cna[ti][j][0]++;
                common[ti][j].emplace_back(m);
            }
        }
    }
}


void get_common_bonds(const int& triclinic,
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box,
    const int ti,
    const vector<vector<double>>& positions,
    const vector<double>& cutoff,
    const vector<vector<int>>& neighbors,
    vector<vector<vector<int>>>& cna,
    vector<vector<vector<int>>>& common,
    vector<vector<vector<int>>>& bonds){
    
    int c1, c2, maxbonds, minbonds;
    double d, diffx, diffy, diffz;

    bonds[ti].clear();
    bonds[ti].resize(neighbors[ti].size());

    for(int k=0; k<neighbors[ti].size(); k++){
        for(int l=0; l<cna[ti][k][0]; l++){
            bonds[ti][k].emplace_back(0);
        }
        
        for(int l=0; l<cna[ti][k][0]-1; l++){
            for(int m=l+1; m<cna[ti][k][0]; m++){
                
                c1 = common[ti][k][l];
                c2 = common[ti][k][m];
                
                d = get_abs_distance(positions[c1], positions[c2],
                    triclinic, rot, rotinv, box, 
                    diffx, diffy, diffz);
                
                if(d <= cutoff[ti]){
                    cna[ti][k][1]++;
                    bonds[ti][k][l]++;
                    bonds[ti][k][m]++;
                }
            }
        }

        maxbonds = 0;
        minbonds = 8;
        
        for(int l=0; l<cna[ti][k][0]; l++){
            maxbonds = max(bonds[ti][k][l], maxbonds);
            minbonds = min(bonds[ti][k][l], minbonds);
        }
        cna[ti][k][2] = maxbonds;
        cna[ti][k][3] = minbonds;    
    }
}


void identify_cn12(py::dict& atoms,
    const int& triclinic,
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box){

    int c1, c2, c3, c4;
    int nfcc, nhcp, nico;

    vector<vector<double>> positions = atoms[py::str("positions")].cast<vector<vector<double>>>();
    vector<double> cutoff = atoms[py::str("cutoff")].cast<vector<double>>();
    vector<int> structure = atoms[py::str("structure")].cast<vector<int>>();
    vector<vector<int>> neighbors = atoms[py::str("neighbors")].cast<vector<vector<int>>>();

    int nop = positions.size();
    vector<vector<vector<int>>> cna(nop);
    vector<vector<vector<int>>> common(nop);
    vector<vector<vector<int>>> bonds(nop);

    for(int ti=0; ti<nop; ti++){
        if(structure[ti] == 0){
            get_common_neighbors(triclinic, rot, rotinv, box, ti, positions,
                cutoff, neighbors, cna, common);
            get_common_bonds(triclinic, rot, rotinv, box, ti, positions,
                cutoff, neighbors, cna, common, bonds);

            nfcc = 0;
            nhcp = 0;
            nico = 0;

            for(int k=0; k<neighbors[ti].size(); k++){
                
                c1 = cna[ti][k][0];
                c2 = cna[ti][k][1];
                c3 = cna[ti][k][2];
                c4 = cna[ti][k][3];

                if((c1==4) && (c2==2) && (c3==1) && (c4==1)){
                    nfcc++;
                }
                else if ((c1==4) && (c2==2) && (c3==2) && (c4==0)){
                    nhcp++;
                }
                else if ((c1==5) && (c2==5) && (c3==2) && (c4==2)){
                    nico++;
                }

            }
            if(nfcc==12){
                structure[ti] = 1;
            }
            else if((nfcc==6) && (nhcp==6)){
                structure[ti] = 2;   
            }
            else if (nico==12){
                structure[ti] = 4;   
            }
        }
    }
    atoms[py::str("structure")] = structure;
}


void identify_cn14(py::dict& atoms,
    const int& triclinic,
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box){

    int c1, c2, c3, c4;
    int nbcc1, nbcc2;

    vector<vector<double>> positions = atoms[py::str("positions")].cast<vector<vector<double>>>();
    vector<double> cutoff = atoms[py::str("cutoff")].cast<vector<double>>();
    vector<int> structure = atoms[py::str("structure")].cast<vector<int>>();
    vector<vector<int>> neighbors = atoms[py::str("neighbors")].cast<vector<vector<int>>>();

    int nop = positions.size();
    vector<vector<vector<int>>> cna(nop);
    vector<vector<vector<int>>> common(nop);
    vector<vector<vector<int>>> bonds(nop);

    for(int ti=0; ti<nop; ti++){
        if(structure[ti] == 0){
            get_common_neighbors(triclinic, rot, rotinv, box, ti, positions,
                cutoff, neighbors, cna, common);
            get_common_bonds(triclinic, rot, rotinv, box, ti, positions,
                cutoff, neighbors, cna, common, bonds);

            nbcc1 = 0;
            nbcc2 = 0;

            for(int k=0; k<neighbors[ti].size(); k++){
                
                c1 = cna[ti][k][0];
                c2 = cna[ti][k][1];
                c3 = cna[ti][k][2];
                c4 = cna[ti][k][3];

                if((c1==4) && (c2==4) && (c3==2) && (c4==2)){
                    nbcc1++;
                }
                else if ((c1==6) && (c2==6) && (c3==2) && (c4==2)){
                    nbcc2++;
                }

            }
            cout<<nbcc1<<" "<<nbcc2<<endl;
            if((nbcc1==6) && (nbcc2==8)){
                structure[ti] = 3;   
            }
        }
    }
    atoms[py::str("structure")] = structure;
}