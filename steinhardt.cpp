#include "steinhardt.h"
#include <iostream>
#include <iomanip>


System::System(){
    
    nop = -1;
}

System::~System(){
    
    delete [] atoms;
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
        }
  
    }
  
 }

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


void System::get_all_neighbors(){

    double nd,d;
    double diffx,diffy,diffz;
    double r,theta,phi;
    int ccount[nop];

    for (int ti = 0;ti<nop;ti++){
        
        atoms[ti].n_neighbors=0;
        for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++){
                        
            atoms[ti].neighbors[tn] = NILVALUE;
            atoms[ti].neighbordist[tn] = -1.0;
        }
    }

    for (int ti=0; ti<(nop-1); ti++){
        for (int tj=ti+1; tj<nop; tj++){
                        
            d = get_abs_distance(ti,tj,diffx,diffy,diffz); 
            if (d < neighbordistance){

                atoms[ti].neighbors[atoms[ti].n_neighbors] = tj; 
                atoms[ti].neighbordist[atoms[ti].n_neighbors] = d; 
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
    factor = ((2.0*double(l) + 1.0)*FACTORIALS[l-m]) / (4.0*PI*FACTORIALS[l+m]);
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
                realti += realYLM;
                imgti += imgYLM;
            }
            
            realti = realti/(double(nn));
            imgti = imgti/(double(nn));
            atoms[ti].realQ6[mi+6] = realti;
            atoms[ti].imgQ6[mi+6] = imgti;
        }
    }
}


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


void System::calculate_frenkel_numbers(){
        
    int frenkelcons;
    double scalar;
        
    for (int ti= 0;ti<nop;ti++){
                
        frenkelcons = 0;
        atoms[ti].avq6q6 = 0.0;
        for (int c = 0;c<atoms[ti].n_neighbors;c++){
                        
            scalar = get_number_from_bond(ti,atoms[ti].neighbors[c]);
            if (scalar > 0.5) frenkelcons += 1;
            atoms[ti].avq6q6 += scalar;
        }
        
        atoms[ti].frenkelnumber = frenkelcons;
        atoms[ti].avq6q6 /= atoms[ti].n_neighbors;
    
    }
}


int System::cluster_criteria(int ti,int criterium){
        
    int value=0;
          
    if (criterium == 0){
                
        if ( (atoms[ti].frenkelnumber > minfrenkel) && (atoms[ti].avq6q6 > 0.5) ){ 
            value = 1; 
        }  
        else{
            value = 0;
        }
    }
    return value;
}


void System::find_solids(){

    int val;
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
            if (freq[ti]>max) max=freq[ti];
        }

        return max;
} 


int System::calculate_nucsize()
{
        int criterium;
        int greatestbelongsto;
        //Find all particles within a radius of neighbourdistancess
        read_particle_file();
        get_all_neighbors();
        //Get Q6 values
        calculate_complexQLM_6();
        //and the number of bonds to find the largest cluster
        calculate_frenkel_numbers();

        find_solids();
        find_clusters();
        greatestbelongsto = largest_cluster();
        return greatestbelongsto;
}



Atom::Atom(){

}

Atom::~Atom(){
        
}

void System::set_minfrenkel(int nn) { minfrenkel = nn; }
void System::set_inputfile(string nn) { inputfile = nn; }
void System::set_neighbordistance(double nn) { neighbordistance = nn; }