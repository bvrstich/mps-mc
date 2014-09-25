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

   //initialize the dimensions of the problem, set the trial
   global::init(D,d,L);

   int Nw = 100000;

   VMC vmc(Nw);
   vmc.walk(3000);

   char filename[200];
   sprintf(filename,"output/VMC/L=%d/D=%d.walk",L,D);

   vmc.dump(filename);

}
