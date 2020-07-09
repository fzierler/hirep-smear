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

int main(int argc, char* argv[]){

   	int N = atoi(argv[1]);
	group::N = N;

    string acfname("spn_sun_algconv.h");// algebra converter file
    write_acf(string(argv[0]), acfname, N);

    return 0;

}

