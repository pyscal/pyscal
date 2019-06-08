#include "steinhardt.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include "voro++.hh"
#include "string.h"

using namespace voro;

/*
Constructor for the system.
*/
System::System(){
    
    nop = -1;
    maxclusterid = -1;
    neighborsfound = 0;
    qsfound = 0;
    fileread = 0;
    filter = 0;
    triclinic = 0;

}

/*
Destructor of the system class
 */
System::~System(){
    
    //delete [] atoms;
}

/*
Calculate factorial of a number
 */
double System::dfactorial(int l,int m){

    double fac = 1.00;
    for(int i=0;i<2*m;i++){
        fac*=double(l+m-i);
    }
    return (1.00/fac);

}

void System::assign_triclinic_params(vector<vector<double>> drot, vector<vector<double>> drotinv){

    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            rot[i][j] = drot[i][j];
            rotinv[i][j] = drotinv[i][j];
        }
    }
    
    triclinic = 1;
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

void System::read_particle_file(string nn){
    
    inputfile = nn;

    double posx,posy,posz;
    int id;
    double xsizeinf,ysizeinf,zsizeinf,xsizesup,ysizesup,zsizesup;
    double dummy;                        //dummy variable
    char dummy_char[256];                //dummy line
    ifstream confFile;
  
    confFile.open(inputfile.c_str(),ifstream::in);

    if (confFile.is_open()){ 
        
        confFile.getline(dummy_char,256);
        confFile.getline(dummy_char,256);
        confFile.getline(dummy_char,256);
        confFile >> nop;
        atoms = new Atom[nop];
    
        confFile.getline(dummy_char,256);
        confFile.getline(dummy_char,256);
        confFile >> xsizeinf;
        confFile >> xsizesup;
        confFile >> ysizeinf;
        confFile >> ysizesup;
        confFile >> zsizeinf;
        confFile >> zsizesup;;
        confFile.getline(dummy_char,256);
        confFile.getline(dummy_char,256);
  
        boxx = xsizesup - xsizeinf;
        boxy = ysizesup - ysizeinf;
        boxz = zsizesup - zsizeinf;
        boxdims[0][0] = xsizeinf;
        boxdims[0][1] = xsizesup;
        boxdims[1][0] = ysizeinf;
        boxdims[1][1] = ysizesup;
        boxdims[2][0] = zsizeinf;
        boxdims[2][1] = zsizesup;

        //so lets read the particles positions
        for (int ti = 0;ti<nop;ti++){
            confFile>>id;
            confFile>>dummy;
            confFile>>dummy;
            confFile>>posx;
            confFile>>posy;
            confFile>>posz;
            confFile>>dummy;
            confFile>>dummy;
            confFile>>dummy;
      
            atoms[ti].posx = posx;
            atoms[ti].posy = posy;
            atoms[ti].posz = posz;
            atoms[ti].id = id;
            atoms[ti].belongsto = -1;
            atoms[ti].issolid = 0; 
            atoms[ti].loc = ti;
            atoms[ti].isneighborset = 0;


            atoms[ti].n_neighbors=0;
            
            for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++){          
                atoms[ti].neighbors[tn] = NILVALUE;
                atoms[ti].neighbordist[tn] = -1.0;
            }
            
        }
  
    }

    fileread = 1;
  
 }



//this function allows for handling custom formats of atoms and so on
void System::assign_particles( vector<Atom> atomitos, vector<vector<double>> boxd ){
    //atomitos are just a list of Atom objects
    //boxd is a vector of 6 values - [xlow, xhigh, ylow, yhigh, zlow, zhigh]
    nop = atomitos.size();
    atoms = new Atom[nop];
    //triclinic = 0;

    boxdims[0][0] = boxd[0][0];
    boxdims[0][1] = boxd[0][1];
    boxdims[1][0] = boxd[1][0];
    boxdims[1][1] = boxd[1][1];
    boxdims[2][0] = boxd[2][0];
    boxdims[2][1] = boxd[2][1];
    
    boxx = boxd[0][1] - boxd[0][0];
    boxy = boxd[1][1] - boxd[1][0];
    boxz = boxd[2][1] - boxd[2][0];

    for(int ti=0; ti<nop; ti++){
        
        atoms[ti].posx = atomitos[ti].posx;
        atoms[ti].posy = atomitos[ti].posy;
        atoms[ti].posz = atomitos[ti].posz;
        atoms[ti].id = atomitos[ti].id;
        atoms[ti].type = atomitos[ti].type;
        atoms[ti].belongsto = -1;
        atoms[ti].issolid = 0;
        atoms[ti].loc = ti;
        atoms[ti].isneighborset = 0;
        atoms[ti].custom = atomitos[ti].custom;
        atoms[ti].n_neighbors=0;
        atoms[ti].isneighborset = 0;
        
        for (int tn = 0; tn<MAXNUMBEROFNEIGHBORS; tn++){          
            atoms[ti].neighbors[tn] = NILVALUE;
            atoms[ti].neighbordist[tn] = -1.0;
        }

        for (int tn = 0; tn<11; tn++){
            atoms[ti].q[tn] = -1;
            atoms[ti].aq[tn] = -1;
            for (int tnn =0; tnn<25; tnn++){
                atoms[ti].realq[tn][tnn] = -1;
                atoms[ti].imgq[tn][tnn] = -1;
                atoms[ti].arealq[tn][tnn] = -1;
                atoms[ti].aimgq[tn][tnn] = -1;
            } 
        }

        for (int tn = 0; tn<4; tn++){
            atoms[ti].vorovector[tn] = -1;
        }

        
    }

    //atomitos.shrink_to_fit();

    fileread = 1;
}


//needs two version of the function; one for fast inbuilt calculation.
//the other for being accessed to the python interface

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


void System::reset_all_neighbors(){
    for (int ti = 0;ti<nop;ti++){
        
        atoms[ti].n_neighbors=0;
        for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++){
                        
            atoms[ti].neighbors[tn] = NILVALUE;
            atoms[ti].neighbordist[tn] = -1.0;
        }
    }
}

void System::sfilter(int fno){

    filter = fno;
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

void System::get_all_neighbors_normal(){

    
    double d;
    double diffx,diffy,diffz;
    double r,theta,phi;

    if (!fileread) { read_particle_file(inputfile); }


    for (int ti=0; ti<nop; ti++){
        if (atoms[ti].isneighborset == 0){
            for (int tj=ti; tj<nop; tj++){
                if(ti==tj) { continue; }
                d = get_abs_distance(ti,tj,diffx,diffy,diffz); 
                if (d < neighbordistance){
                    if ((filter == 1) && (atoms[ti].type != atoms[tj].type)){
                        continue;
                    }
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
                }
            }
        }
        //if neighbors are already read in
        else {
            //only loop over neighbors
            for (int tj=0; tj<atoms[ti].n_neighbors; tj++){
        
                d = get_abs_distance(ti,atoms[ti].neighbors[tj],diffx,diffy,diffz); 
                //atoms[ti].neighbors[atoms[ti].n_neighbors] = tj; 
                atoms[ti].neighbordist[tj] = d;
                //weight is set to 1.0, unless manually reset
                //atoms[ti].neighborweight[atoms[ti].n_neighbors] = 1.00; 
                atoms[ti].n_diffx[tj] = diffx;
                atoms[ti].n_diffy[tj] = diffy;
                atoms[ti].n_diffz[tj] = diffz;
                convert_to_spherical_coordinates(diffx, diffy, diffz, r, phi, theta);
                atoms[ti].n_r[tj] = r;
                atoms[ti].n_phi[tj] = phi;
                atoms[ti].n_theta[tj] = theta;
                //atoms[ti].n_neighbors += 1;   
                
            }   
        }
    }

    //mark end of neighbor calc
    neighborsfound = 1;

}

//overloaded function; would be called
//if neighbor method voronoi is selected.
void System::get_all_neighbors_voronoi(){


    double d;
    double diffx,diffy,diffz;
    double r,theta,phi;
    int i;
    int ti,id,tnx,tny,tnz;
    int n3, n4, n5, n6;
    double rx,ry,rz,tsum, fa, x, y, z, vol;
    vector<int> neigh,f_vert;
    vector<double> facearea, v;
    voronoicell_neighbor c;
    vector< vector<double> > nweights;
    vector< vector<int> > nneighs;
    vector<int> idss;
    //vector<int> nvector;
    double weightsum;

    if (!fileread) { read_particle_file(inputfile); }

    pre_container pcon(boxdims[0][0],boxdims[1][1],boxdims[1][0],boxdims[1][1],boxdims[2][0],boxdims[2][1],true,true,true);
    for(int i=0; i<nop; i++){
        pcon.put(i, atoms[i].posx, atoms[i].posy, atoms[i].posz);
    }
    pcon.guess_optimal(tnx,tny,tnz);        
    container con(boxdims[0][0],boxdims[1][1],boxdims[1][0],boxdims[1][1],boxdims[2][0],boxdims[2][1],tnx,tny,tnz,true,true,true, nop);
    pcon.setup(con);

    c_loop_all cl(con);
    if (cl.start()) do if(con.compute_cell(c,cl)) {
            ti=cl.pid();
            c.face_areas(facearea);
            c.neighbors(neigh);
            c.face_orders(f_vert);
            c.vertices(x,y,z,v);
            vol = c.volume();
            tsum = 0;
            vector <double> dummyweights;
            vector <int> dummyneighs;

            //only loop over neighbors
            weightsum = 0.0;
            for (int i=0; i<facearea.size(); i++){
            	weightsum += facearea[i];
            }
            
            //find vertices index
            n3 = 0;
            n4 = 0;
            n5 = 0;
            n6 = 0;

            for(int i=0; i<f_vert.size(); i++){
                if ((facearea[i]/weightsum) < 0.002){
                    continue;
                }

                if (f_vert[i] == 3){
                    n3++;
                }
                else if (f_vert[i] == 4){
                    n4++;
                }
                else if (f_vert[i] == 5){
                    n5++;
                }
                else if (f_vert[i] == 6){
                    n6++;
                }
            }

            //assign to nvector
            //nvector.clear();
            atoms[ti].vorovector[0] = n3;
            atoms[ti].vorovector[1] = n4;
            atoms[ti].vorovector[2] = n5;
            atoms[ti].vorovector[3] = n6;
            atoms[ti].volume = vol;
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
                atoms[ti].neighborweight[tj] = facearea[tj]/weightsum;
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

}


void System::get_all_neighbors_adaptive(int nlimit, double prefactor){

    
    double d, dcut;
    double diffx,diffy,diffz;
    double r,theta,phi;
    vector<int> nids;
    vector<double> dists, sorted_dists;

    //double prefactor = 1.21;
    double summ;

    if (!fileread) { read_particle_file(inputfile); }

    for (int ti=0; ti<nop; ti++){
        if (atoms[ti].isneighborset == 0){
            //clear vectors
            nids.clear();
            dists.clear();
            sorted_dists.clear();

            //start looping over every other particle
            for (int tj=0; tj<nop; tj++){
                if(ti==tj) { continue; }
                d = get_abs_distance(ti,tj,diffx,diffy,diffz);
                //add the ids of atom, and distance to the list
                nids.emplace_back(tj);
                dists.emplace_back(d);
                sorted_dists.emplace_back(d);
            }

            //we have all the info now. Pick the top six
            //first sort distances
            
            sort(sorted_dists.begin(), sorted_dists.end());
            summ = 0;
            for(int i=0; i<nlimit; i++){
                summ+=sorted_dists[i];
            }
            dcut = prefactor*(1.0/float(nlimit))*summ;

            //now we are ready to loop over again, but over the lists        
            for(int j=0; j<dists.size(); j++){
                int tj = nids[j];
                if (dists[j] < dcut){

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

                }                
            }
        }
        
        //if neighbors are already read in
        else {
            //only loop over neighbors
            for (int tj=0; tj<atoms[ti].n_neighbors; tj++){
        
                d = get_abs_distance(ti,atoms[ti].neighbors[tj],diffx,diffy,diffz); 
                //atoms[ti].neighbors[atoms[ti].n_neighbors] = tj; 
                atoms[ti].neighbordist[tj] = d;
                //weight is set to 1.0, unless manually reset
                //atoms[ti].neighborweight[atoms[ti].n_neighbors] = 1.00; 
                atoms[ti].n_diffx[tj] = diffx;
                atoms[ti].n_diffy[tj] = diffy;
                atoms[ti].n_diffz[tj] = diffz;
                convert_to_spherical_coordinates(diffx, diffy, diffz, r, phi, theta);
                atoms[ti].n_r[tj] = r;
                atoms[ti].n_phi[tj] = phi;
                atoms[ti].n_theta[tj] = theta;
                //atoms[ti].n_neighbors += 1;   
                
            }   
        }
    }

    //mark end of neighbor calc
    neighborsfound = 1;

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
            atoms[ti].realQ6[mi+6] = realti;
            atoms[ti].imgQ6[mi+6] = imgti;
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
            if(weightsum>1.01){
                realti = realti/weightsum;
                imgti = imgti/weightsum;                
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

//also has to be overloaded - could be a useful function
double System::get_number_from_bond(int ti,int tj){

    double sumSquareti,sumSquaretj;
    double realdotproduct,imgdotproduct;
    double connection;
    sumSquareti = 0.0;
    sumSquaretj = 0.0;
    realdotproduct = 0.0;
    imgdotproduct = 0.0;

    for (int mi = 0;mi < 13;mi++){
                
        sumSquareti += atoms[ti].realQ6[mi]*atoms[ti].realQ6[mi] + atoms[ti].imgQ6[mi] *atoms[ti].imgQ6[mi];
        sumSquaretj += atoms[tj].realQ6[mi]*atoms[tj].realQ6[mi] + atoms[tj].imgQ6[mi] *atoms[tj].imgQ6[mi];
        realdotproduct += atoms[ti].realQ6[mi]*atoms[tj].realQ6[mi];
        imgdotproduct  += atoms[ti].imgQ6[mi] *atoms[tj].imgQ6[mi];
    }
        
    connection = (realdotproduct+imgdotproduct)/(sqrt(sumSquaretj)*sqrt(sumSquareti));
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

    for (int mi = 0;mi < 13;mi++){
                
        sumSquareti += atom1.realQ6[mi]*atom1.realQ6[mi] + atom1.imgQ6[mi] *atom1.imgQ6[mi];
        sumSquaretj += atom2.realQ6[mi]*atom2.realQ6[mi] + atom2.imgQ6[mi] *atom2.imgQ6[mi];
        realdotproduct += atom1.realQ6[mi]*atom2.realQ6[mi];
        imgdotproduct  += atom1.imgQ6[mi] *atom2.imgQ6[mi];
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
            if (scalar > threshold) frenkelcons += 1;
            atoms[ti].avq6q6 += scalar;
        }
        
        atoms[ti].frenkelnumber = frenkelcons;
        atoms[ti].avq6q6 /= atoms[ti].n_neighbors;
    
    }
}

//again to be overloaded?
//maybe not now-im lazy
int System::cluster_criteria(int ti,int criterium){
        
    int value=0;
          
    if (criterium == 0){
                
        if ( (atoms[ti].frenkelnumber > minfrenkel) && (atoms[ti].avq6q6 > avgthreshold) ){ 
            value = 1; 
        }  
        else{
            value = 0;
        }
    }
    return value;
}


void System::find_solids(){

    int criteria = 0;

    for (int ti= 0;ti<nop;ti++){
                
        atoms[ti].issolid = cluster_criteria(ti,criteria);
    }
}


void System::find_clusters(){

        for (int ti= 0;ti<nop;ti++){
                
            if (!atoms[ti].issolid) continue;
            if (atoms[ti].belongsto==-1) {atoms[ti].belongsto = atoms[ti].id; }
            for (int c = 0;c<atoms[ti].n_neighbors;c++){

                if(!atoms[atoms[ti].neighbors[c]].issolid) continue;
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
        if(!atoms[neigh].issolid) continue;
        if (atoms[neigh].belongsto==-1){
            atoms[neigh].belongsto = clusterindex;
            harvest_cluster(neigh, clusterindex);
        }
    }
}

void System::find_clusters_recursive(){

    int clusterindex;
    clusterindex = 0;

    for (int ti= 0;ti<nop;ti++){
        if (!atoms[ti].issolid) continue;
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

int System::calculate_nucsize()
{

        int greatestbelongsto;
        //Find all particles within a radius of neighbourdistancess
        //read_particle_file();
        get_all_neighbors_normal();
        //Get Q6 values
        //cout<<"step 1"<<endl;
        calculate_complexQLM_6();
        //cout<<"step 2"<<endl;
        //and the number of bonds to find the largest cluster
        calculate_frenkel_numbers();
        //cout<<"step 3"<<endl;

        find_solids();
        //cout<<"step 4"<<endl;
        //find_clusters();
        find_clusters_recursive();
        //cout<<"step 5"<<endl;
        greatestbelongsto = largest_cluster();
        //cout<<"step 6"<<endl;
        get_largest_cluster_atoms();
        //cout<<"step 7"<<endl;
        return greatestbelongsto;
}




//access functions for system
//------------------------------------------------------------------------------------------------------------------------
//void System::set_inputfile(string nn) { inputfile = nn; }
void System::set_neighbordistance(double nn) { neighbordistance = nn; }
void System::set_nucsize_parameters(double cutoff, int n1, double n2, double n3 ) { neighbordistance = cutoff; minfrenkel = n1; threshold = n2; avgthreshold = n3; }
Atom System::gatom(int i) { return atoms[i]; }
void System::satom(Atom atom1) { 
    int idd = atom1.loc;
    atoms[idd] = atom1;
}

//add function to return nop
int System::gnop() { return nop; }
//int System::gnop() { return nop; } 
//add function to pack and return the whole set of atoms
vector<Atom> System::gallatoms(){
    vector<Atom> allatoms;
    allatoms.reserve(nop);
    for(int i=0;i<nop;i++){
        allatoms.emplace_back(atoms[i]);
    }

    return allatoms;
}

int System::glargestclusterid() { return maxclusterid; }

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

vector<double> System::gbox(){
    vector<double> qres;
    qres.reserve(3);
    qres.emplace_back(boxx);
    qres.emplace_back(boxy);
    qres.emplace_back(boxz);
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

vector<double> System::gboxdims(){
    vector<double> qres;
    qres.reserve(6);
    for(int i=0;i<3;i++){
        for(int j=0;j<2;j++){
            qres.emplace_back(boxdims[i][j]);
        }
    }
    return qres; 
}

//functions for atoms
//-------------------------------------------------------------------------------------------------------------------------
Atom::Atom(){ }
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

double Atom::gq(int qq){ return q[qq-2]; }
int Atom::gid(){ return id; }
void Atom::sid(int idd){ id=idd; }
int Atom::gloc(){ return loc; }
void Atom::sloc(int idd){ loc=idd; }
int Atom::gtype(){ return type; }
void Atom::stype(int idd){ type=idd; }
double Atom::gvolume(){ return volume; }
void Atom::svolume(double vv){ volume = vv; }

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


void Atom::sq(int qq, double qval){ q[qq-2] = qval; }

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

vector<int> Atom::gvorovector(){
    vector <int> voro;
    for(int i=0; i<4; i++){
        voro.emplace_back(vorovector[i]);
    }
    return voro;
}

void Atom::svorovector(vector<int> voro){
    for(int i=0; i<4; i++){
        vorovector[i] = voro[i];
    }
}