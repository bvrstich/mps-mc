#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using std::ofstream;

#include "include.h"

namespace global{

   int DT;

   int L;

   int d;

   MPS<double> mps;

   vector< SL_MPS > U;
   vector< SL_MPS > I;

   vector< Walker > backup_walker;

   int omp_num_threads;

   Random RN;

   /**
    * @param DT_in virtual dimension of the trial
    * @param d_in physical dimension
    * @param L_in dimension of the lattice
    */
   void init(int DT_in,int d_in,int L_in){

      L = L_in;

      d = d_in;

      DT = DT_in;

#ifdef _OPENMP
      omp_num_threads = omp_get_max_threads();
#else
      omp_num_threads = 1;
#endif

      mps.resize(L);

      U.resize(omp_num_threads);
      I.resize(omp_num_threads);

      backup_walker.resize(omp_num_threads);
      
      char filename[200];
      sprintf(filename,"/home/bright/bestanden/results/mps-mc/trial/Heis_1D/L=%d/Psi0/seba_D=%d.mps",L,DT);

      mps.load(filename);

      Walker walker;

      U[0] = SL_MPS(mps,walker);

      for(int thr = 1;thr < omp_num_threads;++thr)
         U[thr] = U[0];

      for(int thr = 0;thr < omp_num_threads;++thr)
         I[thr] = U[thr];

   }

   //!function which generates random complex numbers uniformly on a square of side 2 [(-1,1):(-1,1)]
   template<>
      complex<double> rgen(){ 

         return complex<double>(2.0*RN() - 1.0,2.0*RN() - 1.0); 

      }

   //!function which generates uniform random numbers between [-1:1]
   template<>
      double rgen(){ 

         return 2.0*RN() - 1.0;

      }

   //!function which generates uniform random numbers between [0:1]
   double rgen_pos(){ 

      return RN();

   }



}
