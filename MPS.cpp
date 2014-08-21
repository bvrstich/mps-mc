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
template<typename T>
MPS<T>::MPS() : vector< TArray<T,3> >(L) { }

/** 
 * standard constructor: just takes in
 * @param D_in virtual max bond dimension
 * allocates the tensors and fills them randomly
 */
template<typename T>
MPS<T>::MPS(int D_in) : vector< TArray<T,3> >(L) {

   D = D_in;

   vector<int> vdim(L + 1);

   vdim[0] = 1;

   for(int i = 1;i < L;++i){

      int tmp = vdim[i - 1] * d;

      if(tmp < D)
         vdim[i] = tmp;
      else 
         vdim[i] = D;

   }

   vdim[L] = 1;

   for(int i = L - 1;i > 0;--i){

      int tmp = vdim[i + 1] * d;

      if(tmp < vdim[i])
         vdim[i] = tmp;

   }

   for(int i = 0;i < this->size();++i){

      (*this)[i].resize(d,vdim[i],vdim[i+1]);
      (*this)[i].generate(rgen<T>);

   }

}

/** 
 * @param filename inputfile
 */
template<typename T>
void MPS<T>::load(const char *filename){

   ifstream in(filename);

   int a,b,c;

   in >> a >> b >> c;

   vector<int> vdim(L + 1);

   for(int i = 0;i <= L;++i)
      in >> i >> vdim[i];

   int dim;

   for(int i = 0;i < L;++i){

      (*this)[i].resize(d,vdim[i],vdim[i+1]);

      int teller;

      in >> a >> b;

      for(int s = 0;s < d;++s)
         for(int k = 0;k < vdim[i + 1];++k)
            for(int j = 0;j < vdim[i];++j){

               in >> i >> teller >> (*this)[i](s,j,k);

            }

   }

}


/**
 * copy constructor
 */
template<typename T>
MPS<T>::MPS(const MPS<T> &mps_copy) : vector< TArray<T,3> >(mps_copy) {

   D = mps_copy.gD();

}

/**
 * empty destructor
 */
template<typename T>
MPS<T>::~MPS(){ }

/**
 * @return virtual dimension of the MPS
 */
template<typename T>
int MPS<T>::gD() const {

   return D;

}

template MPS<double>::MPS();
template MPS< complex<double> >::MPS();

template MPS<double>::MPS(int);
template MPS< complex<double> >::MPS(int);

template MPS<double>::MPS(const MPS<double> &);
template MPS< complex<double> >::MPS(const MPS< complex<double> > &);

template MPS<double>::~MPS();
template MPS< complex<double> >::~MPS();

template int MPS<double>::gD() const;
template int MPS< complex<double> >::gD() const;

template void MPS<double>::load(const char *);
template void MPS< complex<double> >::load(const char *);
