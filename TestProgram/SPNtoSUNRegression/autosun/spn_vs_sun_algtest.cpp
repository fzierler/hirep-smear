#include <cmath>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <tuple>
#include <map>
#include <vector>
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
            cout << "Problem with SUN generator no. " << isu << endl;
            cout << check<< endl;
            exit(1);
        }
    }
    cout << "SUN CHECK OK" << endl;
    for(int isp = 0; isp<N*(N-1)/2;++isp){
        complex ccheck;
        tmp.mult(TSPN[isp],TSPN[isp]);
        trace(ccheck,tmp);
        ccheck -= 0.5;
        double check = ccheck.re*ccheck.re + ccheck.im*ccheck.im; 
        if(fabs(check) > 1e-14 ){
            cout << "Problem with SPN generator no. " << isp << endl;
            cout << check<< endl;
            exit(1);
        }
    }
    cout << "SPN CHECK OK" << endl;
    delete[] TSUN;
    delete[] TSPN;

}


int main(int argc, char* argv[]){

    string acfname("spn_sun_algconv.h");// algebra converter file
    string mtpname("spnalgmacros.cpp"); // macro test program

	int N = atoi(argv[1]);
	group::N = N;
    testGenerators();

    write_acf(string(argv[0]), acfname, N);
    write_mtp(string(argv[0]), acfname, N, mtpname);
    system("g++ -o spnalgmacros spnalgmacros.cpp && ./spnalgmacros");

    return 0;

}

