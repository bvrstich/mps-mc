#ifndef GLOBAL_H
#define GLOBAL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

using std::ostream;
using std::vector;
using std::complex;

#include "MPS.h"
#include "SL_MPS.h"

using namespace btas;

namespace global {

   extern Random RN;

   //!x dimension of the lattice
   extern int L;

   //!physical dimension of sites
   extern int d;

   //!virtual dimension of the trial
   extern int DT;

   //!single layer mps for energy calculation
   extern vector< SL_MPS > U;

   //!single layer mps for energy calculation
   extern vector< SL_MPS > I;

   //!number of omp threads
   extern int omp_num_threads;

   //!trial wavefunction
   extern MPS<double> mps;

   void init(int,int,int);

   template<typename T>
      T rgen();

   //!test value which sets the minimal overlap
   extern double cutoff;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
