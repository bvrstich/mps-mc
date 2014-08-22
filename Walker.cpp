#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

using std::cout;
using std::ifstream;
using std::ofstream;
using std::endl;

#include "include.h"

using namespace global;

/**
 * construct a Walker object: initialize on AF state
 * @param option decides wich AF state, start with up or down: standard == 0
 * @param n_trot_in number of trotter terms
 */
Walker::Walker(int option) : std::vector< bool >( L ){

   weight = 1.0;

   sign = 1;

   for(int site = 0;site < L;++site){

      if( (site + option) % 2 == 0)
         (*this)[site] = true;
      else
         (*this)[site] = false;

   }

}

/**
 * copy constructor
 * @param walker input Walker object to be copied
 */
Walker::Walker(const Walker &walker) : std::vector< bool >(walker) {

   this->weight = walker.gWeight();
   this->nn_over = walker.gnn_over();
   this->EL = walker.gEL();

   this->sign = walker.gsign();

}

/**
 * destructor
 */
Walker::~Walker(){ }

/**
 * @return the sign of the walker
 */
int Walker::gsign() const {

   return sign;

}

/** 
 * @return the weight corresponding to the walker
 */
double Walker::gWeight() const{

   return weight; 

}

/**
 * muliply the weight by a factor
 */
void Walker::multWeight(double factor){

   weight *= factor; 

}

/**
 * set new weight
 */
void Walker::sWeight(double new_weight){

   weight = new_weight;

}

/** 
 * @return the overlap of the walker with the Trial
 */
double Walker::gOverlap() const{

   return nn_over[0]; 

}

/** 
 * @return the vector containing the overlaps of all the neighbouring walker states with the trial
 */
const vector<double> &Walker::gnn_over() const {

   return nn_over; 

}

/** 
 * @param index of the neighbour walker
 * @return the overlap of the walker with index 'index' with the trial
 */
double Walker::gnn_over(int index) const {

   return nn_over[index]; 

}

/** 
 * @return the local energy
 */
double Walker::gEL() const{

   return EL; 

}

ostream &operator<<(ostream &output,const Walker &walker_p){

   for(int site = 0;site < L;++site)
      output << walker_p[site];

   return output;

}

/**
 * flip the sign of the walker
 */
void Walker::sign_flip(){

   sign *= -1;

}

/**
 * get the 'potential' energy of the walker: < \sum_i Sz_i Sz_{i+1} >
 */
double Walker::pot_en() const {

   double tmp = 0.0;

   for(int site = 0;site < L - 1;++site){

      //Sz Sz
      if( (*this)[site] == (*this)[site + 1] )//up up or down down
         tmp += 0.25;
      else //up down or down up
         tmp -= 0.25;

   }

   return tmp;

}

/**
 * calculate the local energy expectation value and overlap with the accesible states
 * @param mps trial wave function represented as a matrix product state
 */
void Walker::calc_EL(const MPS< double > &mps){

#ifdef _OPENMP
   int myID = omp_get_thread_num();
#else
   int myID = 0;
#endif

   nn_over.clear();

   U[myID].fill(false,mps,*this);
   I[myID].fill(true,mps,*this);

   //construct right renormalized operator
   vector< TArray<double,1> > R(L - 2);

   int m,n;
   double ward;

   //rightmost site
   int dim = U[myID][L-1].shape(0);

   R[L-3].resize(dim);

   blas::copy(dim, U[myID][L-1].data(), 1, R[L-3].data(), 1);

   for(int i = L - 2;i > 1;--i){

      m = U[myID][i].shape(0);
      n = U[myID][i].shape(1);

      R[i-2].resize(m);

      blas::gemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, U[myID][i].data(), n, R[i - 1].data(), 1, 0.0, R[i - 2].data(), 1);

   }

   //now the left going operators: regular and inverse
   TArray<double,1> LU;
   TArray<double,1> LI;

   TArray<double,1> tmp;

   LU.resize(dim);
   LI.resize(dim);

   blas::copy(dim, U[myID][0].data(), 1, LU.data(), 1);
   blas::copy(dim, I[myID][0].data(), 1, LI.data(), 1);

   //put the overlap here later!
   nn_over.push_back(0);

   //now calculate the energy!
   for(int i = 1;i < L - 1;++i){

      //contribution to the energy
      if( (*this)[i - 1] != (*this)[i]){

         m = I[myID][i].shape(0);
         n = I[myID][i].shape(1);

         tmp.resize(n);

         blas::gemv(CblasRowMajor, CblasTrans, m, n, 1.0, I[myID][i].data(), n, LI.data(), 1, 0.0, tmp.data(), 1);

         nn_over.push_back(blas::dot(n, tmp.data(), 1, R[i-1].data(), 1));

      }

      //new left sides

      //if it contributes make LI:
      if( (*this)[i] != (*this)[i + 1]){

         m = I[myID][i].shape(0);
         n = I[myID][i].shape(1);

         LI.resize(n);

         blas::gemv(CblasRowMajor, CblasTrans, m, n, 1.0, I[myID][i].data(), n, LU.data(), 1, 0.0, LI.data(), 1);

      }

      //new left unity
      m = U[myID][i].shape(0);
      n = U[myID][i].shape(1);

      tmp.resize(n);
      blas::gemv(CblasRowMajor, CblasTrans, m, n, 1.0, U[myID][i].data(), n, LU.data(), 1, 0.0, tmp.data(), 1);

      LU.resize(n);
      blas::copy(n, tmp.data(), 1, LU.data(), 1);

   }

   //final contribution to the energy
   if( (*this)[L - 2] != (*this)[L - 1])
      nn_over.push_back(blas::dot(LI.size(), LI.data(), 1, I[myID][L-1].data(), 1));

   //overlap
   nn_over[0] = blas::dot(LU.size(), LU.data(), 1, U[myID][L-1].data(), 1);

   //calculate the local energy
   EL = this->pot_en();

   for(int i = 1;i < nn_over.size();++i)
      EL += 0.5 * nn_over[i]/nn_over[0];

}

/**
 * save the walker to file
 * @param filename name of the file
 */
void Walker::save(const char *filename){

   ofstream fout(filename);

   for(int site = 0;site < L;++site)
      fout << site << "\t" << (*this)[site] << endl;

}

/**
 * count the number of differences between two walkers
 * @param walker_i input Walker
 */
int Walker::num_diff(const Walker &walker_i) const {

   int num_diff = 0;

   for(int site = 0;site < L;++site)
      if( (*this)[site] != walker_i[site])
         num_diff++;

   return num_diff;

}

/**
 * load the walker
 */
void Walker::load(const char *filename){

   ifstream fin(filename);

   bool tmp;

   for(int site = 0;site < L;++site){

      fin >> site >> tmp; 
      (*this)[site] = tmp;

   }

}
