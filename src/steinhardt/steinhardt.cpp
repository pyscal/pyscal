#include "steinhardt.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include "voro++.hh"
#include "string.h"

using namespace voro;

System::System(){
    
    nop = -1;
    maxclusterid = -1;
    neighborsfound = 0;
    qsfound = 0;
    fileread = 0;

}

System::~System(){
    
    delete [] atoms;
}


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

void System::read_particle_file(){
    
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

            for (int ti = 0;ti<nop;ti++){
                atoms[ti].n_neighbors=0;
                for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++){          
                    atoms[ti].neighbors[tn] = NILVALUE;
                    atoms[ti].neighbordist[tn] = -1.0;
                }
            }
        }
  
    }

    fileread = 1;
  
 }

void System::read_particle_instance(int startblock,int natoms){
    

    string line,str;
    //stringstream ss;
    int count = 1;
    int minc = 0;
    double xsizeinf,ysizeinf,zsizeinf,xsizesup,ysizesup,zsizesup;
    nop = natoms;
    atoms = new Atom[nop];
    int block = natoms +9;
    double posx,posy,posz;
    int id;
    double dummy;
    int idummy;
    //cout<<startblock*block<<endl;
    ifstream infile(inputfile.c_str());
    
    if (infile.is_open()){
        while(getline(infile,line)){
            //now we have to skip everything until the lines we need
            if (count>startblock*block){
                //this is in the reading range
                if (count==6+startblock*block){
                    sscanf(line.c_str(),"%lf %lf", &xsizeinf, &xsizesup);
                    boxx = xsizesup - xsizeinf;
                    boxdims[0][0] = xsizeinf;
                    boxdims[0][1] = xsizesup;

                }
                else if (count==7+startblock*block){
                    sscanf(line.c_str(),"%lf %lf", &ysizeinf, &ysizesup);
                    boxy = ysizesup - ysizeinf;
                    boxdims[1][0] = ysizeinf;
                    boxdims[1][1] = ysizesup;

                }
                else if (count==8+startblock*block){
                    sscanf(line.c_str(),"%lf %lf", &zsizeinf, &zsizesup);
                    boxz = zsizesup - zsizeinf;
                    boxdims[2][0] = zsizeinf;
                    boxdims[2][1] = zsizesup;
                }

                else if (count>9+startblock*block){
                    sscanf(line.c_str(),"%d %d %lf %lf %lf %lf %lf %lf %lf", &id, &idummy, &dummy, &posx, &posy, &posz, &dummy, &dummy, &dummy);
                    atoms[minc].posx = posx;
                    atoms[minc].posy = posy;
                    atoms[minc].posz = posz;
                    atoms[minc].id = id;
                    atoms[minc].belongsto = -1;
                    atoms[minc].issolid = 0; 
                    atoms[minc].loc = minc-9;
                    atoms[minc].isneighborset = 0;
                    minc++;

                }

                

                if (count==block+startblock*block) { break; }


            }
            //cout<<"count "<<count<<endl;
            //break loop if exceeded
            count++;
        }
    }  
}

//this is the old - actual function
//now an overloaded version would be provided, the aim of this is to provide two infos-
//1. The list of neighbors for each atom. This would be position of neighbor on the list of atoms
//2. The list of weights for each neighbors. This is used to scale the YLM contributions.
////wait - its better to just assign a variable weights to Atom and functions to set these and the weights

void System::assign_particles( vector<Atom> atomitos, vector<double> boxd ){
//dont know if this will be faster-
//what it does would be to recieve a 
//vector of atoms and then process it
//add it to the normal list we have.
    nop = atomitos.size();
    atoms = new Atom[nop];

    boxx = boxd[0];
    boxy = boxd[1];
    boxz = boxd[2];

    for(int ti=0;ti<nop;ti++){
        atoms[ti].posx = atomitos[ti].posx;
        atoms[ti].posy = atomitos[ti].posy;
        atoms[ti].posz = atomitos[ti].posz;
        atoms[ti].id = atomitos[ti].id;
        atoms[ti].belongsto = -1;
        atoms[ti].issolid = 0;
        atoms[ti].isneighborset = 0;
        atoms[ti].custom = atomitos[ti].custom;
    }

    atomitos.shrink_to_fit();

    fileread = 1;
}


//needs two version of the function; one for fast inbuilt calculation.
//the other for being accessed to the python interface

double System::get_abs_distance(int ti ,int tj,double &diffx ,double &diffy,double &diffz){
  
    double abs;
    diffx = atoms[tj].posx - atoms[ti].posx;
    diffy = atoms[tj].posy - atoms[ti].posy;
    diffz = atoms[tj].posz - atoms[ti].posz;
        
    //nearest image
    if (diffx> boxx/2.0) {diffx-=boxx;};
    if (diffx<-boxx/2.0) {diffx+=boxx;};
    if (diffy> boxy/2.0) {diffy-=boxy;};
    if (diffy<-boxy/2.0) {diffy+=boxy;};
    if (diffz> boxz/2.0) {diffz-=boxz;};
    if (diffz<-boxz/2.0) {diffz+=boxz;};
    abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
    return abs;
}

//function for binding
double System::get_abs_distance(Atom atom1 , Atom atom2 ){
  
    double abs;
    double diffx = atom1.posx - atom2.posx;
    double diffy = atom1.posy - atom2.posy;
    double diffz = atom1.posz - atom2.posz;
        
    //nearest image
    if (diffx> boxx/2.0) {diffx-=boxx;};
    if (diffx<-boxx/2.0) {diffx+=boxx;};
    if (diffy> boxy/2.0) {diffy-=boxy;};
    if (diffy<-boxy/2.0) {diffy+=boxy;};
    if (diffz> boxz/2.0) {diffz-=boxz;};
    if (diffz<-boxz/2.0) {diffz+=boxz;};
    abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
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
void System::get_all_neighbors(){

    
    double d;
    double diffx,diffy,diffz;
    double r,theta,phi;

    if (!fileread) { read_particle_file(); }

    //for (int ti = 0;ti<nop;ti++){
    //    
    //    atoms[ti].n_neighbors=0;
    //    for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++){
    //                    
    //        atoms[ti].neighbors[tn] = NILVALUE;
    //        atoms[ti].neighbordist[tn] = -1.0;
    //    }
    //}


    for (int ti=0; ti<nop; ti++){
        if (atoms[ti].isneighborset == 0){
            for (int tj=ti; tj<nop; tj++){
                if(ti==tj) { continue; }
                d = get_abs_distance(ti,tj,diffx,diffy,diffz); 
                if (d < neighbordistance){

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
void System::get_all_neighbors( string &jkl){

    if(strcmp(jkl.c_str(),"voronoi")!=0){
        return;
    }

    double d;
    double diffx,diffy,diffz;
    double r,theta,phi;
    int i;
    int ti,id,tnx,tny,tnz;
    double rx,ry,rz,tsum, fa;
    vector<int> neigh,f_vert;
    vector<double> facearea;
    voronoicell_neighbor c;
    vector< vector<double> > nweights;
    vector< vector<int> > nneighs;
    vector<int> idss;
    double weightsum;

    if (!fileread) { read_particle_file(); }

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
            tsum = 0;
            vector <double> dummyweights;
            vector <int> dummyneighs;

            //only loop over neighbors
            weightsum = 0.0;
            for (int i=0; i<facearea.size(); i++){
                weightsum += facearea[i];
            }
            for (int tj=0; tj<neigh.size(); tj++){
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
void System::calculate_q(){
        
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
        get_all_neighbors();
    }

    //note that the qvals will be in -2 pos
    //q2 will be in q0 pos and so on
    double weightsum;  
    // nop = parameter.nop;
    for (int ti= 0;ti<nop;ti++){
        
        nn = atoms[ti].n_neighbors;
        for(int tq=0;tq<lenqs;tq++){
            //find which q?
            q = reqdqs[tq];
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
void System::calculate_aq(){
        
    //nn = number of neighbors
    int nn;
    double realti,imgti;
    //double realYLM,imgYLM;
    int q;
    double summ, weightsum;

    if (!qsfound) { set_reqd_qs(rq_backup); calculate_q(); }
    //note that the qvals will be in -2 pos
    //q2 will be in q0 pos and so on
        
    // nop = parameter.nop;
    for (int ti= 0;ti<nop;ti++){
        
        nn = atoms[ti].n_neighbors;
        
        for(int tq=0;tq<lenqs;tq++){
            //find which q?
            q = reqdqs[tq];
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
        get_all_neighbors();
        //Get Q6 values
        calculate_complexQLM_6();
        //and the number of bonds to find the largest cluster
        calculate_frenkel_numbers();

        find_solids();
        find_clusters();
        greatestbelongsto = largest_cluster();
        get_largest_cluster_atoms();
        return greatestbelongsto;
}




//access functions for system
//------------------------------------------------------------------------------------------------------------------------
void System::set_inputfile(string nn) { inputfile = nn; }
void System::set_neighbordistance(double nn) { neighbordistance = nn; }
void System::set_nucsize_parameters(int n1, double n2, double n3 ) { minfrenkel = n1; threshold = n2; avgthreshold = n3; }
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

double Atom::gq(int qq){ return q[qq-2]; }
int Atom::gid(){ return id; }

//aceesss funcs 
vector<double> Atom::gx(){ 
    vector<double> pos;
    pos.emplace_back(posx);
    pos.emplace_back(posy);
    pos.emplace_back(posz);
    return pos; 
}

vector<int> Atom::gcluster(){
    vector<int> cl;
    cl.emplace_back(issolid);
    cl.emplace_back(issurface);
    cl.emplace_back(lcluster);
    cl.emplace_back(belongsto);
    return cl;
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
