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
using std::ifstream;

#include "include.h"

using namespace global;

/** 
 * empty constructor:
 */
SL_MPS::SL_MPS() : vector< TArray<double,2> >(L) { }

/**
 * construct and fill the SL_MPS object by contracting input MPS with input Walker
 * @param mps input MPS
 * @param walker input Walker
 */
SL_MPS::SL_MPS(const MPS<double> &mps,const Walker &walker) : vector< TArray<double,2> >(L) { 

   for(int i = 0;i < L;++i){

      (*this)[i].resize(mps[i].shape(1),mps[i].shape(2));

      int dim = mps[i].shape(1) * mps[i].shape(2);

      blas::copy(dim, mps[i].data() + dim * walker[i], 1, (*this)[i].data(), 1);

   }

}

/**
 * copy constructor
 */
SL_MPS::SL_MPS(const SL_MPS &slmps_copy) : vector< TArray<double,2> >(slmps_copy) { }

/**
 * empty destructor
 */
SL_MPS::~SL_MPS(){ }

/**
 * fill the SL_MPS object by contracting input MPS with input Walker
 * @param inverse if true take the inverse
 * @param mps input MPS
 * @param walker input Walker
 */
void SL_MPS::fill(bool inverse,const MPS<double> &mps,const Walker &walker) { 

   for(int i = 0;i < L;++i){

      int dim = mps[i].shape(1) * mps[i].shape(2);

      if(inverse)
         blas::copy(dim, mps[i].data() + dim * (!walker[i]), 1, (*this)[i].data(), 1);
      else
         blas::copy(dim, mps[i].data() + dim * walker[i], 1, (*this)[i].data(), 1);

   }

}
