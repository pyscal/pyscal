#include "steinhardt.h"
#include <iostream>
#include <iomanip>

System::System()
{
        //CParameter parameter = param;

        nop = -1;
}

System::~System()
{
        //delete parameter;
        delete [] molecules;
}


void System::readParticleFile()
{
  double posx,posy,posz;
  
  int id;

  double xsizeinf,ysizeinf,zsizeinf,xsizesup,ysizesup,zsizesup;
  double dummy;                        //dummy variable
  char dummy_char[256];                //dummy line
  ifstream confFile;
  
  confFile.open(inputfile.c_str(),ifstream::in);

  if (confFile.is_open())
  { 
    
   confFile.getline(dummy_char,256);
   confFile.getline(dummy_char,256);
   confFile.getline(dummy_char,256);
   confFile >> nop;
   molecules   = new Atom[nop];
    
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
  for (int ti = 0;ti<nop;ti++)
    {
      confFile>>id;
      confFile>>dummy;
      confFile>>dummy;
      confFile>>posx;
      confFile>>posy;
      confFile>>posz;
      confFile>>dummy;
      confFile>>dummy;
      confFile>>dummy;
     
      //cout<<posx<<"\n"<<posy<<"\n"<<posz<<"\n";
      
      molecules[ti].posx = posx;
      molecules[ti].posy = posy;
      molecules[ti].posz = posz;
      molecules[ti].id = id;
      molecules[ti].belongsto = -1;
      molecules[ti].issolid = 0; 
    }
  
  }
  
 }

double System::get_absDistance(int ti ,int tj,double &diffx ,double &diffy,double &diffz)
{
  
        double abs;
        diffx = molecules[tj].posx - molecules[ti].posx;
        diffy = molecules[tj].posy - molecules[ti].posy;
        diffz = molecules[tj].posz - molecules[ti].posz;
        //nearest image
        if (diffx >  boxx/2.0) {diffx = diffx - boxx;};
        if (diffx < -boxx/2.0) {diffx = diffx + boxx;};
        if (diffy >  boxy/2.0) {diffy = diffy - boxy;};
        if (diffy < -boxy/2.0) {diffy = diffy + boxy;};
        if (diffz >  boxz/2.0) {diffz = diffz - boxz;};
        if (diffz < -boxz/2.0) {diffz = diffz + boxz;};
        abs = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);
        return abs;
}


 void System::get_AllNeighborsAndDistances()
{

        double nd,d;
        double diffx,diffy,diffz;
        double r,theta,phi;
        nd = neighbordistance;
        //int nop;
        //nop = parameter.nop;
        int ccount[nop];

        for (int ti = 0;ti<nop;ti++)
        {
                ccount[ti]=0;
                for (int tn = 0;tn<MAXNUMBEROFNEIGHBORS;tn++)
                {
                        molecules[ti].neighbors[tn] = nilvalue;
                        molecules[ti].neighbordist[tn] = -1.0;
                }
        }


        for (int ti=0; ti<(nop-1); ti++)
        {
                for (int tj=ti+1; tj<nop; tj++)
                {
                        d = get_absDistance(ti,tj,diffx,diffy,diffz); 
                        if (d < nd) {
                                molecules[ti].neighbors[ccount[ti]] = tj; 
                                molecules[ti].neighbordist[ccount[ti]] = d; 
                                molecules[ti].n_diffx[ccount[ti]] = diffx;
                                molecules[ti].n_diffy[ccount[ti]] = diffy;
                                molecules[ti].n_diffz[ccount[ti]] = diffz;
                                convert_SphericalCoordinates(diffx, diffy, diffz, r, phi, theta);
                                molecules[ti].n_r[ccount[ti]] = r;
                                molecules[ti].n_phi[ccount[ti]] = phi;
                                molecules[ti].n_theta[ccount[ti]] = theta;
                                ccount[ti] += 1;   

                                molecules[tj].neighbors[ccount[tj]] = ti;
                                molecules[tj].neighbordist[ccount[tj]] = d;
                                molecules[tj].n_diffx[ccount[tj]] = -diffx;
                                molecules[tj].n_diffy[ccount[tj]] = -diffy;
                                molecules[tj].n_diffz[ccount[tj]] = -diffz;
                                convert_SphericalCoordinates(-diffx, -diffy, -diffz, r, phi, theta);
                                molecules[tj].n_r[ccount[tj]] = r;
                                molecules[tj].n_phi[ccount[tj]] = phi;
                                molecules[tj].n_theta[ccount[tj]] = theta;
                                ccount[tj] +=1;
                        }
                }
        }
        for (int ti=0; ti<nop; ti++){
                molecules[ti].n_neighbors = ccount[ti];
        }
}


double System::PLM(int l, int m, double x)
{
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

void System::convert_SphericalCoordinates(double x, double y, double z, double &r, double &phi, double &theta)
{
        //r = sqrt(x*x + y*y + z*z);
        //if (y >= 0.0) { phi = acos(x/sqrt(x*x + y*y)); } else {phi = 2.0*pi - acos(x/sqrt(x*x + y*y)); }
        //if (fabs(x)<0.00000001 && fabs(y)<0.000000001) phi = 0.0;
        //theta = pi/2 - atan(z/sqrt(x*x + y*y));
        r = sqrt(x*x+y*y+z*z);
        theta = acos(z/r);
        phi = atan2(y,x);
}


void System::YLM(int l, int m, double theta, double phi, double &realYLM, double &imgYLM)
{
        double factor;
        double m_PLM;
        m_PLM = PLM(l,m,cos(theta));
        factor = ((2.0*double(l) + 1.0)*FACTORIALS[l-m]) / (4.0*pi*FACTORIALS[l+m]);
        realYLM = sqrt(factor) * m_PLM * cos(double(m)*phi);
        imgYLM  = sqrt(factor) * m_PLM * sin(double(m)*phi);
}


void System::QLM(int l,int m,double theta,double phi,double &realYLM, double &imgYLM )
{
        realYLM = 0.0;
        imgYLM = 0.0;
        if (m < 0) {
                YLM(l, abs(m), theta, phi, realYLM, imgYLM);
                realYLM = pow(-1.0,m)*realYLM;
                imgYLM = pow(-1.0,m)*imgYLM;
        }
        else
        {
                YLM(l, m, theta, phi, realYLM, imgYLM);
        }
}

void System::calculate_complexQLM_6()
{
        //nn = number of neighbors
        int nn;
        double realti,imgti;
        double realYLM,imgYLM;
        
       // nop = parameter.nop;
        for (int ti= 0;ti<nop;ti++)
        {
                nn = molecules[ti].n_neighbors;
                for (int mi = -6;mi < 7;mi++)
                {
                        realti = 0.0;
                        imgti = 0.0;
                        for (int ci = 0;ci<nn;ci++)
                        {
                                QLM(6,mi,molecules[ti].n_theta[ci],molecules[ti].n_phi[ci],realYLM, imgYLM);
                                realti += realYLM;
                                imgti += imgYLM;
                        }
                        realti = realti/(double(nn));
                        imgti = imgti/(double(nn));
                        molecules[ti].realQ6[mi+6] = realti;
                        molecules[ti].imgQ6[mi+6] = imgti;
                }
        }
}


double System::get_NumberFromBond(int ti,int tj)
{
        double sumSquareti,sumSquaretj;
        double realdotproduct,imgdotproduct;
        double connection;
        sumSquareti = 0.0;
        sumSquaretj = 0.0;
        realdotproduct = 0.0;
        imgdotproduct = 0.0;

        for (int mi = 0;mi < 13;mi++)
        {
                sumSquareti += molecules[ti].realQ6[mi]*molecules[ti].realQ6[mi] + molecules[ti].imgQ6[mi] *molecules[ti].imgQ6[mi];
                sumSquaretj += molecules[tj].realQ6[mi]*molecules[tj].realQ6[mi] + molecules[tj].imgQ6[mi] *molecules[tj].imgQ6[mi];
                realdotproduct += molecules[ti].realQ6[mi]*molecules[tj].realQ6[mi];
                imgdotproduct  += molecules[ti].imgQ6[mi] *molecules[tj].imgQ6[mi];
        }
        connection = (realdotproduct+imgdotproduct)/(sqrt(sumSquaretj)*sqrt(sumSquareti));
        return connection;
}


void System::calculate_frenkelNumbers()
{
        int frenkelcons;
        double scalar;
        ofstream fout("out.dat");
        for (int ti= 0;ti<nop;ti++)
        {
                frenkelcons = 0;
                molecules[ti].avq6q6 = 0.0;
                for (int c = 0;c<molecules[ti].n_neighbors;c++)
                {
                        scalar = get_NumberFromBond(ti,molecules[ti].neighbors[c]);
                        if (scalar > 0.5) frenkelcons += 1;
      molecules[ti].avq6q6 += scalar;
                }
                molecules[ti].frenkelnumber = frenkelcons;
    molecules[ti].avq6q6 /= molecules[ti].n_neighbors;
    fout<<molecules[ti].avq6q6<<endl;
                //cout<<this->molecules[ti].frenkelnumber<<""<<"\n";
        }
  fout.close();

 

}


int System::clusterCriterium(int ti,int criterium)
{
        int value;
        value = 0; 
        //cout<< parameter->minfrenkel<<""<<"\n"; 
        if (criterium == 0)
        {
                // cout<<parameter->minfrenkel<<"\t"<<this->molecules[ti].frenkelnumber<<"\n";
                if ( (molecules[ti].frenkelnumber > minfrenkel) && (molecules[ti].avq6q6 > 0.5) ){ 
                        value = 1; 
                }  
                else {
                        value = 0;
                }
        }
        return value;
}


void System::find_solids(){

        int val;
        int criteria = 0;

        for (int ti= 0;ti<nop;ti++)
        {
                molecules[ti].issolid = clusterCriterium(ti,criteria);
        }
}


void System::find_clusters(){

        for (int ti= 0;ti<nop;ti++)
        {
                if (!molecules[ti].issolid) continue;
                if (molecules[ti].belongsto==-1) {molecules[ti].belongsto = molecules[ti].id; }
                for (int c = 0;c<molecules[ti].n_neighbors;c++)
                {   
                        if(!molecules[molecules[ti].neighbors[c]].issolid) continue;
                        if (molecules[molecules[ti].neighbors[c]].belongsto==-1){
                                molecules[molecules[ti].neighbors[c]].belongsto = molecules[ti].belongsto;
                        }
                        else{
                                molecules[ti].belongsto = molecules[molecules[ti].neighbors[c]].belongsto;  
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
                if (molecules[ti].belongsto==-1) continue;
                freq[molecules[ti].belongsto-1]++;
        }

        int max=0;
        for (int ti= 0;ti<nop;ti++)
        {
                if (freq[ti]>max) max=freq[ti];
        }

        return max;
} 


int System::calculate_largestClusterparameter_Full()
{
        int criterium;
        int greatestbelongsto;
        //Find all particles within a radius of neighbourdistancess
        readParticleFile();
        get_AllNeighborsAndDistances();
        //Get Q6 values
        calculate_complexQLM_6();
        //and the number of bonds to find the largest cluster
        calculate_frenkelNumbers();

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