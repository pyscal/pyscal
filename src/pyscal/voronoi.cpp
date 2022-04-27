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
#include "voro++.hh"

using namespace voro;

void get_all_neighbors_voronoi(py::dict& atoms,
    const double neighbordistance,
    const int triclinic,
    const vector<vector<double>> rot, 
    const vector<vector<double>> rotinv,
    const vector<double> box,
    const double face_area_exponent)
    {

    double d;
    double diffx,diffy,diffz;
    double tempr,temptheta,tempphi;
    vector<double> diffi, diffj, pos;
    int tnx,tny,tnz, ti, tj, nverts;
    double rx,ry,rz,tsum, fa, x, y, z, vol, weightsum;

    vector<int> neigh,f_vert, vert_nos;
    vector<double> facearea, v, faceperimeters;
    voronoicell_neighbor c;


    vector<vector<double>> positions = atoms[py::str("positions")].cast<vector<vector<double>>>();
    vector<bool> mask_1 = atoms[py::str("mask_1")].cast<vector<bool>>();
    vector<bool> mask_2 = atoms[py::str("mask_2")].cast<vector<bool>>();

    int nop = positions.size();
    vector<vector<int>> neighbors(nop);
    vector<vector<double>> neighbordist(nop);
    vector<vector<double>> neighborweight(nop);
    vector<vector<vector<double>>> diff(nop);
    vector<vector<double>> r(nop);
    vector<vector<double>> phi(nop);
    vector<vector<double>> theta(nop);
    vector<double> cutoff(nop);

    //specific properties related to  Voronoi
    vector<double> volume(nop);
    vector<vector<int>> face_vertices(nop);
    vector<vector<double>> face_perimeters(nop);
    vector<vector<double>> vertex_vectors(nop);
    vector<vector<int>> vertex_numbers(nop);
    vector<vector<vector<double>>> vertex_positions(nop);

    pre_container pcon(0.00, box[0], 0.00, box[1], 0.0, box[2], true, true, true);
    for(int i=0; i<nop; i++){
        pos = positions[i];
        pos = remap_atom_into_box(pos, triclinic, rot, rotinv, box);
        pcon.put(i, pos[0], pos[1], pos[2]);
    }
    pcon.guess_optimal(tnx, tny, tnz);
    //container con(boxdims[0][0],boxdims[1][1],boxdims[1][0],boxdims[1][1],boxdims[2][0],boxdims[2][1],tnx,tny,tnz,true,true,true, nop);
    container con(0.00, box[0], 0.00, box[1], 0.0, box[2], tnx, tny, tnz, true, true, true, nop);
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

        weightsum = 0.0;
        for (int i=0; i<facearea.size(); i++){
            weightsum += pow(facearea[i], face_area_exponent);
        }

        volume[ti] = vol;
        vertex_vectors[ti] = v;
        vertex_numbers[ti] = vert_nos;
        cutoff[ti] = cbrt(3*vol/(4*3.141592653589793));

        //clean up and add vertex positions
        nverts = int(v.size())/3;
        pos = positions[ti];
        for(int si=0; si<nverts; si++){
            vector<double> temp;
            int li=0;
            for(int vi=si*3; vi<(si*3+3); vi++){
                temp.emplace_back(v[vi]+pos[li]);
                li++;
            }
            vertex_positions[ti].emplace_back(temp);
        }


        for (int tj=0; tj<neigh.size(); tj++){
            d = get_abs_distance(positions[ti], positions[tj],
                triclinic, rot, rotinv, box, 
                diffx, diffy, diffz);
            neighbors[ti].emplace_back(tj);
            neighbordist[ti].emplace_back(d);
            neighborweight[ti].emplace_back(pow(facearea[tj], face_area_exponent)/weightsum);

            face_vertices[ti].emplace_back(f_vert[tj]);
            face_perimeters[ti].emplace_back(faceperimeters[tj]);

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

    } while (cl.inc());


    //get unique vertex positions


    //calculation over lets assign
    atoms[py::str("neighbors")] = neighbors;
    atoms[py::str("neighbordist")] = neighbordist;
    atoms[py::str("neighborweight")] = neighborweight;
    atoms[py::str("diff")] = diff;
    atoms[py::str("r")] = r;
    atoms[py::str("theta")] = theta;
    atoms[py::str("phi")] = phi;
    atoms[py::str("cutoff")] = cutoff;
    atoms[py::str("voronoi_volume")] = volume;
    atoms[py::str("face_vertices")] = face_vertices;
    atoms[py::str("face_perimeters")] = face_perimeters;
    atoms[py::str("vertex_vectors")] = vertex_vectors;
    atoms[py::str("vertex_numbers")] = vertex_numbers;
    atoms[py::str("vertex_positions")] = vertex_positions;
} 


void clean_voronoi_vertices(py::dict& atoms,
    const double neighbordistance,
    const int triclinic,
    const vector<vector<double>> rot, 
    const vector<vector<double>> rotinv,
    const vector<double> box,
    const double distance_cutoff){

    vector<vector<vector<double>>> positions = atoms[py::str("vertex_positions")].cast<vector<vector<vector<double>>>>();
    vector<bool> ghost = atoms[py::str("ghost")].cast<vector<bool>>();
    int nop = positions.size();

    vector<vector<double>> vertex_positions_all;
    vector<int> to_remove;

    double d, diffx, diffy, diffz;

    for(int ti=0; ti<nop; ti++){
        if (ghost[ti]==true) continue;
        for(int vi=0; vi<positions[ti].size(); vi++){
            auto pos = remap_atom_into_box(positions[ti][vi], triclinic, rot, rotinv, box);
            if ((pos[0] < 0) || (pos[0] > box[0])) continue;
            if ((pos[1] < 0) || (pos[1] > box[1])) continue;
            if ((pos[2] < 0) || (pos[2] > box[2])) continue;
            vertex_positions_all.emplace_back(pos);
        }        

    }

    for(int ii=0; ii<vertex_positions_all.size(); ii++){
        auto pos1 = vertex_positions_all[ii];
        for(int jj=ii+1; jj<vertex_positions_all.size(); jj++){
            auto pos2 = vertex_positions_all[jj];
                d = get_abs_distance(pos1, pos2,
                    triclinic, rot, rotinv, box, 
                    diffx, diffy, diffz);
                if (d < distance_cutoff){
                    to_remove.emplace_back(jj);            
                }
        }
    }
    
    atoms[py::str("vertex_positions_all_skipcheck")] = vertex_positions_all;
    atoms[py::str("to_remove_skipcheck")] = to_remove;
    
}	