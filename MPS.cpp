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
 * copy constructor
 */
template<typename T>
MPS<T>::MPS(const Walker &walker) : vector< TArray<T,3> >(L) {

   D = 1;

   vector<int> vdim(L + 1);


   for(int i = 0;i <= L;++i)
      vdim[i] = 1;

   for(int i = 0;i < this->size();++i){

      (*this)[i].resize(d,vdim[i],vdim[i+1]);
      (*this)[i] = 0.0;

      if(walker[i])
         (*this)[i](1,0,0) = 1.0;
      else
         (*this)[i](0,0,0) = 1.0;

   }


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

/**
 * calculate the 1D Heisenberg model energy expectation value for this MPS
 */

template<>
double MPS<double>::energy() const{

#ifdef _OPENMP
   int myID = omp_get_thread_num();
#else
   int myID = 0;
#endif

   DArray<2> Sz(d,d);
   DArray<2> Sp(d,d);
   DArray<2> Sm(d,d);

   //Sz
   Sz(0,0) = -0.5;
   Sz(0,1) = 0.0;
   Sz(1,0) = 0.0;
   Sz(1,1) = 0.5;
   
   //S+
   Sp(0,0) = 0.0;
   Sp(0,1) = 1.0;
   Sp(1,0) = 0.0;
   Sp(1,1) = 0.0;
   
   //S-
   Sm(0,0) = 0.0;
   Sm(0,1) = 0.0;
   Sm(1,0) = 1.0;
   Sm(1,1) = 0.0;

   double ener = 0.0;

   enum {j,k,l,m,n,o};

   vector< DArray<4> > R(L-1);

   DArray<5> I;

   Contract(1.0,mps[L-1],shape(j,k,l),mps[L-1],shape(j,m,n),0.0,R[L-2],shape(k,l,m,n));

   for(int i = L-2;i > 0;--i){

      I.clear();
      Contract(1.0,mps[i],shape(j,k,l),R[i],shape(l,m,n,o),0.0,I,shape(j,k,m,n,o));

      R[i-1].clear();
      Contract(1.0,mps[i],shape(j,l,n),I,shape(j,k,m,n,o),0.0,R[i-1],shape(k,m,l,o));

   }

   //first site
   DArray<4> LU;

   Contract(1.0,mps[0],shape(j,k,l),mps[0],shape(j,m,n),0.0,LU,shape(k,l,m,n));

   //operators
   DArray<4> Lz;
   DArray<4> Lp;
   DArray<4> Lm;

   DArray<3> tmp;

   //Sz
   Contract(1.0,mps[0],shape(j,k,l),Sz,shape(j,m),0.0,tmp,shape(m,k,l));
   Contract(1.0,tmp,shape(j,k,l),mps[0],shape(j,m,n),0.0,Lz,shape(k,l,m,n));

   //Sp
   Contract(1.0,mps[0],shape(j,k,l),Sp,shape(j,m),0.0,tmp,shape(m,k,l));
   Contract(1.0,tmp,shape(j,k,l),mps[0],shape(j,m,n),0.0,Lp,shape(k,l,m,n));

   //Sm
   Contract(1.0,mps[0],shape(j,k,l),Sm,shape(j,m),0.0,tmp,shape(m,k,l));
   Contract(1.0,tmp,shape(j,k,l),mps[0],shape(j,m,n),0.0,Lm,shape(k,l,m,n));

   for(int i = 1;i < L - 1;++i){

      //Sz
      tmp.clear();
      Contract(1.0,mps[i],shape(j,k,l),Sz,shape(j,m),0.0,tmp,shape(m,k,l));

      I.clear();
      Contract(1.0,Lz,shape(k,l,m,n),tmp,shape(j,l,o),0.0,I,shape(j,k,o,m,n));

      Lz.clear();
      Contract(1.0,I,shape(j,k,o,m,n),mps[i],shape(j,n,l),0.0,Lz,shape(k,o,m,l));

      ener += blas::dot(Lz.size(),Lz.data(),1,R[i].data(),1);

      //create new Lz
      I.clear();
      Contract(1.0,LU,shape(k,l,m,n),tmp,shape(j,l,o),0.0,I,shape(j,k,o,m,n));

      Lz.clear();
      Contract(1.0,I,shape(j,k,o,m,n),mps[i],shape(j,n,l),0.0,Lz,shape(k,o,m,l));
 
      //Sp
      tmp.clear();
      Contract(1.0,mps[i],shape(j,k,l),Sm,shape(j,m),0.0,tmp,shape(m,k,l));

      I.clear();
      Contract(1.0,Lp,shape(k,l,m,n),tmp,shape(j,l,o),0.0,I,shape(j,k,o,m,n));

      Lp.clear();
      Contract(1.0,I,shape(j,k,o,m,n),mps[i],shape(j,n,l),0.0,Lp,shape(k,o,m,l));

      ener += 0.5 * blas::dot(Lp.size(),Lp.data(),1,R[i].data(),1);

      //Sm
      tmp.clear();
      Contract(1.0,mps[i],shape(j,k,l),Sp,shape(j,m),0.0,tmp,shape(m,k,l));

      I.clear();
      Contract(1.0,Lm,shape(k,l,m,n),tmp,shape(j,l,o),0.0,I,shape(j,k,o,m,n));

      Lm.clear();
      Contract(1.0,I,shape(j,k,o,m,n),mps[i],shape(j,n,l),0.0,Lm,shape(k,o,m,l));

      ener += 0.5 * blas::dot(Lm.size(),Lm.data(),1,R[i].data(),1);
      
      //create new Lm
      tmp.clear();
      Contract(1.0,mps[i],shape(j,k,l),Sm,shape(j,m),0.0,tmp,shape(m,k,l));

      I.clear();
      Contract(1.0,LU,shape(k,l,m,n),tmp,shape(j,l,o),0.0,I,shape(j,k,o,m,n));

      Lm.clear();
      Contract(1.0,I,shape(j,k,o,m,n),mps[i],shape(j,n,l),0.0,Lm,shape(k,o,m,l));

      //create new Lp
      tmp.clear();
      Contract(1.0,mps[i],shape(j,k,l),Sp,shape(j,m),0.0,tmp,shape(m,k,l));

      I.clear();
      Contract(1.0,LU,shape(k,l,m,n),tmp,shape(j,l,o),0.0,I,shape(j,k,o,m,n));

      Lp.clear();
      Contract(1.0,I,shape(j,k,o,m,n),mps[i],shape(j,n,l),0.0,Lp,shape(k,o,m,l));

      //finally create new LU
      I.clear();
      Contract(1.0,LU,shape(k,l,m,n),mps[i],shape(j,l,o),0.0,I,shape(j,k,o,m,n));

      LU.clear();
      Contract(1.0,I,shape(j,k,o,m,n),mps[i],shape(j,n,l),0.0,LU,shape(k,o,m,l));

   }
   
   //Sz
   tmp.clear();
   Contract(1.0,mps[L-1],shape(j,k,l),Sz,shape(j,m),0.0,tmp,shape(m,k,l));

   LU.clear();
   Contract(1.0,tmp,shape(j,k,l),mps[L-1],shape(j,m,n),0.0,LU,shape(k,l,m,n));

   ener += blas::dot(Lz.size(),Lz.data(),1,LU.data(),1);

   //Sp
   tmp.clear();
   Contract(1.0,mps[L-1],shape(j,k,l),Sp,shape(j,m),0.0,tmp,shape(m,k,l));

   LU.clear();
   Contract(1.0,tmp,shape(j,k,l),mps[L-1],shape(j,m,n),0.0,LU,shape(k,l,m,n));

   ener += 0.5 * blas::dot(Lm.size(),Lm.data(),1,LU.data(),1);

   //Sm
   tmp.clear();
   Contract(1.0,mps[L-1],shape(j,k,l),Sm,shape(j,m),0.0,tmp,shape(m,k,l));

   LU.clear();
   Contract(1.0,tmp,shape(j,k,l),mps[L-1],shape(j,m,n),0.0,LU,shape(k,l,m,n));

   ener += 0.5 * blas::dot(Lp.size(),Lp.data(),1,LU.data(),1);

   return ener;

}

/**
 * calculate the overlap of the trial with the walker.
 * @param mps_in input MPS
 */
template<>
double MPS<double>::dot(const MPS<double> &mps_in) const{

#ifdef _OPENMP
   int myID = omp_get_thread_num();
#else
   int myID = 0;
#endif

   enum {j,k,l,m,n,o};

   DArray<4> E;
   DArray<5> I;

   Contract(1.0,mps[L-1],shape(j,k,l),mps_in[L-1],shape(j,m,n),0.0,E,shape(k,l,m,n));

   for(int i = L-1;i > 0;--i){

      I.clear();
      Contract(1.0,mps[i-1],shape(j,k,l),E,shape(l,m,n,o),0.0,I,shape(j,k,m,n,o));

      E.clear();
      Contract(1.0,mps_in[i-1],shape(j,l,n),I,shape(j,k,m,n,o),0.0,E,shape(k,m,l,o));

   }

   return E(0,0,0,0);

}

template MPS<double>::MPS();
template MPS< complex<double> >::MPS();

template MPS<double>::MPS(int);
template MPS< complex<double> >::MPS(int);

template MPS<double>::MPS(const MPS<double> &);
template MPS< complex<double> >::MPS(const MPS< complex<double> > &);

template MPS<double>::MPS(const Walker &);
template MPS< complex<double> >::MPS(const Walker &);

template MPS<double>::~MPS();
template MPS< complex<double> >::~MPS();

template int MPS<double>::gD() const;
template int MPS< complex<double> >::gD() const;

template void MPS<double>::load(const char *);
template void MPS< complex<double> >::load(const char *);
