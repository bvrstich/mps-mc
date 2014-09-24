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
/*
   double dtau = 0.01;
   int Nw = 132893;

   vector< Walker > walker(Nw);

   ifstream in("/home/bright/bestanden/programmas/jastrow_gfmc/output/L=100/f=0.300000.walk");

   for(int i = 0;i < Nw;++i){

      for(int j = 0;j < walker[i].size();++j){

         bool tmp;

         in >> tmp;
         walker[i][j] = tmp;

      }

   }

   vector<double> overlap(Nw);

   ifstream in_over("/home/bright/bestanden/programmas/jastrow_gfmc/output/L=100/overlap.out");

   int ind;

   for(int i = 0;i < Nw;++i)
      in_over >> ind >> overlap[i];
*/   
   cout << global::mps.energy() << endl;

}
