#include <iostream>
#include <cmath>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace std;

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
		ylm[l].resize(lm+1);
		for(int m=0; m<=l; m++){
			ylm[l][m].emplace_back(0);
			ylm[l][m].emplace_back(0);
		}
	}

	plm = calculate_plm(lm, costheta, sintheta);

	double zv = 1.0/sqrt(2.0);
	double cosi = 1.0;
	double cosf = cosphi;
	double sini = 0.0;
	double sinf = -sinphi; 
	double ci, si;


	for(int l=0; l<=lm; l++){

		ylm[l][0][0] = plm[l][0]*zv;

	}

	for(int m=1; m<=lm; m++){

		ci = 2*cosphi*cosi-cosf;
		si = 2*cosphi*sini-sinf;
		sinf = sini;
		sini = si;
		cosf = cosi;
		cosi = ci;

		for(int l=m; l<=lm; l++){
			ylm[l][m][0] = plm[l][m]*ci*zv;
			ylm[l][m][1] = plm[l][m]*si*zv;
		}
	}

	return ylm;

}


PYBIND11_MODULE(csh, m) {
    py::options options;
    options.disable_function_signatures();
    m.def("calculate_factors", &calculate_factors);
    m.def("calculate_plm", &calculate_plm);
    m.def("calculate_ylm", &calculate_ylm);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}