#include "steinhardt.h"
#include <ctime>
#include <sstream>
#include <string>
#include <cstdlib>

int main(int argc, char *argv[])
{

        string filename="";
        if (argc > 1){
             filename = string(argv[1]);
        }
        else{
            cout<<"Must provide input file name to program."<<endl;
            cout<<"Exiting program..."<<endl;
            exit(1);
        }
        // The main obejct is created. It hold all the functions and data used in the analysis.
        
        //CParameter *param = new CParameter(3.63,7);
        //CParameter param(3.63,7);
        System m_MolSys;
        m_MolSys.set_minfrenkel(7);
        m_MolSys.set_inputfile(filename);
        m_MolSys.set_neighbordistance(3.63);
        m_MolSys.set_threshold(0.5);
        m_MolSys.set_avgthreshold(0.5);

        int res;

        res = m_MolSys.calculate_nucsize();
        cout<<res<<endl;
}
