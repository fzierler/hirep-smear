#include <cmath>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <tuple>
#include <map>
#include <vector>
#include "../colors.h"
// functions for the test
#include "lib.h"

void testGenerators(){
    smatrix *TSUN,*TSPN;
    group::init(group::N,group::TYPESPN,TSPN);
    group::init(group::N,group::TYPESUN,TSUN);

    using group::N;
    smatrix tmp;

    for(int isu = 0; isu<N*N-1;++isu){
        complex ccheck;
        tmp.mult(TSUN[isu],TSUN[isu]);
        trace(ccheck,tmp);
        ccheck -= 0.5;
        double check = ccheck.re*ccheck.re + ccheck.im*ccheck.im; 
        if(fabs(check) > 1e-14 ){
            cout << BOLDRED "Problem with SUN generator no. " << isu << RESET << endl;
            cout << check<< endl;
            exit(1);
        }
    }
    cout << BOLDGREEN "SUN CHECK OK" RESET << endl;
    for(int isp = 0; isp<N*(N-1)/2;++isp){
        complex ccheck;
        tmp.mult(TSPN[isp],TSPN[isp]);
        trace(ccheck,tmp);
        ccheck -= 0.5;
        double check = ccheck.re*ccheck.re + ccheck.im*ccheck.im; 
        if(fabs(check) > 1e-14 ){
            cout << BOLDRED "Problem with SPN generator no. " << isp << RESET << endl;
            cout << check<< endl;
            exit(1);
        }
    }
    cout << BOLDGREEN "SPN CHECK OK" RESET << endl;
    delete[] TSUN;
    delete[] TSPN;

}


int main(int argc, char* argv[]){

    string acfname("spn_sun_algconv.h");// algebra converter file
    string mtpname_adj("spnalgmacros_adj.cpp"); // macro test program
    string mtpname_asym("spnalgmacros_asym.cpp"); // macro test program

	int N = atoi(argv[1]);
	group::N = N;
    testGenerators();

    write_acf(string(argv[0]), acfname, N);
    write_mtp_adj(string(argv[0]), acfname, N, mtpname_adj);
    write_mtp_asym(string(argv[0]), acfname, N, mtpname_asym);
    system("g++ -o spnalgmacros_adj spnalgmacros_adj.cpp && ./spnalgmacros_adj");
    system("g++ -o spnalgmacros_asym spnalgmacros_asym.cpp && ./spnalgmacros_asym");

    return 0;

}

