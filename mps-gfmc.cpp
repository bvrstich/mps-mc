#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::complex;
using std::vector;

#include "include.h"

int main(int argc,char *argv[]){

   cout.precision(15);

   int L = atoi(argv[1]);
   int d = atoi(argv[2]);
   int D = atoi(argv[3]);
   double f = atof(argv[4]);

   //initialize the dimensions of the problem, set the trial
   global::init(D,d,L,f);

   double dtau = 0.01;
   int Nw = 10000;

   GFMC gfmc(dtau,Nw);
   gfmc.walk(1000000);

}
