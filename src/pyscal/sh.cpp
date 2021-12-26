#include <iostream>
#include <cmath>
#include <vector>
#include "modsystem.h"

void calculate_factors(const int lm, 
	vector<vector<double>>& alm,
	vector<vector<double>>& blm, 
	vector<vector<double>>& clm,
	vector<double>& dl,
	vector<double>& el)
	{
	
	double t1, t2, t3;
	double temp1, temp2, temp3, temp4;

	//resize and set everything to zero
	alm.resize(lm+1);
	blm.resize(lm+1);
	clm.resize(lm+1);
	dl.resize(lm+1);
	el.resize(lm+1);

	for(int l=0; l<=lm; l++){
		for(int m=0; m<=l; m++){
			alm[l].emplace_back(0);
			blm[l].emplace_back(0);
			clm[l].emplace_back(0);
		}
		dl[l] = 0;
		el[l] = 0;
	}


	for(int l=0; l<=lm; l++){
		
		temp1 = 4*l*l-1;
		temp2 = l*l-2*l+1;
		temp3 = 2*l+1;
		temp4 = 2*l-1;

		for(int m=0; m<l-1; m++){
			
			t1 = sqrt((temp1)/(l*l-m*m));
			t2 = -sqrt((temp2-m*m)/(4*temp2-1));
			t3 = sqrt(((l-m)*(l+m)*temp3)/(temp4));

			alm[l][m] = t1;
			blm[l][m] = t2;
			clm[l][m] = t3;

		}
	}

	for(int l=2; l<=lm; l++){
		dl[l] = sqrt(2*(l-1)+3);
		el[l] = sqrt(1+0.5/double(l));
	}
}


vector<vector<double>> calculate_plm(const int lm,
	const double costheta, 
	const double sintheta)
	{

	vector<vector<double>> alm;
	vector<vector<double>> blm; 
	vector<vector<double>> clm;
	vector<double> dl;
	vector<double> el;

	calculate_factors(lm, alm, blm, clm, dl, el);

	double a1 = sqrt(0.5/3.141592653589793);
	double a2 = sqrt(3.0);
	double a3 = sqrt(1.5);

	vector<vector<double>> plm;
	plm.resize(lm+1);
	for(int l=0; l<=lm; l++){
		for(int m=0; m<=l; m++){
			plm[l].emplace_back(0);
		}
	}


	plm[0][0] = a1;
	plm[1][0] = a1*costheta*a2;
	a1 = a1*sintheta*-a3;
	plm[1][1] = a1;


	for(int l=2; l<=lm; l++){

		int m = 0;
		for(int m=0; m<l-1; m++){
			plm[l][m] = alm[l][m]*(costheta*plm[l-1][m]+blm[l][m]*plm[l-2][m]);
		}

		plm[l][l-1] = dl[l]*costheta*a1;
		a1 *= -el[l]*sintheta;
		plm[l][l] = a1;
	}

	return plm;
}


double dfactorial(int l,
	int m){

    double fac = 1.00;
    for(int i=0;i<2*m;i++){
        fac*=double(l+m-i);
    }
    return (1.00/fac);
}

vector<vector<vector<double>>> calculate_ylm(const int lm,
	const double costheta,
	const double sintheta,
	const double cosphi,
	const double sinphi)
	{
	vector<vector<vector<double>>> ylm;
	vector<vector<double>> plm;

	ylm.resize(lm+1);
	
	for(int l=0; l<=lm; l++){
		ylm[l].resize(2*l+1);
		for(int m=0; m<(2*l+1); m++){
			ylm[l][m].emplace_back(0.0);
			ylm[l][m].emplace_back(0.0);
		}
	}
	plm = calculate_plm(lm, costheta, sintheta);
	double zv = 1.0/sqrt(2.0);
	double cosi = 1.0;
	double cosf = cosphi;
	double sini = 0.0;
	double sinf = -sinphi; 
	double ci, si;
	double fa, fb, factor;
	for(int l=0; l<=lm; l++){

		ylm[l][l][0] = plm[l][0]*zv;

	}
	for(int m=1; m<=lm; m++){

		ci = 2*cosphi*cosi-cosf;
		si = 2*cosphi*sini-sinf;
		sinf = sini;
		sini = si;
		cosf = cosi;
		cosi = ci;

		for(int l=m; l<=lm; l++){
			
			fa = plm[l][m]*ci*zv;
			fb = plm[l][m]*si*zv;
			//factor = sqrt(((2.0*double(l) + 1.0)/ (2.0*PI))*dfactorial(l,m));
			factor = 1.00;

			ylm[l][l+m][0] = factor*fa;
			ylm[l][l+m][1] = factor*fb;
			ylm[l][l-m][0] = factor*fa*pow(-1.0,-m);
			ylm[l][l-m][1] = factor*fb*pow(-1.0,-m);
		}
	}

	return ylm;

}

vector<vector<vector<vector<double>>>> calculate_q_atom(const int lm,
	const vector<double>& theta,
	const vector<double>& phi)
{
	vector<vector<vector<vector<double>>>> ylm_atom;
	//this is like a loop over neighbors
	for(int i=0; i<theta.size(); i++){
		ylm_atom.emplace_back(calculate_ylm(lm, cos(theta[i]), sin(theta[i]), cos(phi[i]), sin(phi[i])));
	}
	return ylm_atom;
}

void calculate_q(py::dict& atoms,
	const int lm)
{
	//we need theta and pi
    vector<vector<double>> theta = atoms[py::str("theta")].cast<vector<vector<double>>>();
    vector<vector<double>> phi = atoms[py::str("phi")].cast<vector<vector<double>>>();
    vector<vector<double>> weights = atoms[py::str("neighborweight")].cast<vector<vector<double>>>();
    int nop = theta.size();
    vector<vector<vector<double>>> qlm_real(nop);
    vector<vector<vector<double>>> qlm_img(nop);
    vector<vector<double>> q(nop);

    int nn;
    double summ, weightsum;
    double realti, imgti;
	for (int ti=0; ti<nop; ti++){
		//calculate ylm first
		auto ylm_atom = calculate_q_atom(lm, theta[ti], phi[ti]);	
		qlm_real[ti].resize(lm+1);
		qlm_img[ti].resize(lm+1);
		
		for(int l=0; l<=lm; l++){
			summ = 0;
			for(int m=0; m<(2*l+1); m++){
				realti = 0;
				imgti = 0;
				weightsum = 0;
				for(int ci=0; ci<theta[ti].size(); ci++){
					//TODO: add condition
					realti += weights[ti][ci]*ylm_atom[ci][l][m][0];
					imgti += weights[ti][ci]*ylm_atom[ci][l][m][1];					
					weightsum += weights[ti][ci];
				}
				//TODO: turn off for Voronoi
				realti = realti/float(weightsum);
				imgti = imgti/float(weightsum);

				qlm_real[ti][l].emplace_back(realti);
				qlm_img[ti][l].emplace_back(imgti);

				summ += realti*realti + imgti*imgti;
			}
			summ = pow(((4.0*PI/(2*l+1)) * summ),0.5);
			q[ti].emplace_back(summ);
		}
	}
    atoms[py::str("qlm_real")] = qlm_real;
    atoms[py::str("qlm_img")] = qlm_img;
    atoms[py::str("q")] = q;
}