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

}

/**
 * destructor
 */
Walker::~Walker(){ }

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

   U[myID].fill(false,mps,*this);
   I[myID].fill(true,mps,*this);

   nn_over.clear();

   //construct right renormalized operator for real and inverted lattice
   vector< TArray<double,1> > RU(L - 2);
   vector< TArray<double,1> > RI(L - 2);

   int m,n;
   double value;

   //rightmost site
   int dim = U[myID][L-1].shape(0);

   RU[L-3].resize(dim);
   RI[L-3].resize(dim);

   blas::copy(dim, U[myID][L-1].data(), 1, RU[L-3].data(), 1);
   blas::copy(dim, I[myID][L-1].data(), 1, RI[L-3].data(), 1);

   for(int i = L - 2;i > 1;--i){

      m = U[myID][i].shape(0);
      n = U[myID][i].shape(1);

      RU[i-2].resize(m);
      RI[i-2].resize(m);

      blas::gemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, U[myID][i].data(), n, RU[i - 1].data(), 1, 0.0, RU[i - 2].data(), 1);
      blas::gemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, I[myID][i].data(), n, RI[i - 1].data(), 1, 0.0, RI[i - 2].data(), 1);

   }

   //now the left going operators: regular and inverse
   TArray<double,1> LUU;
   TArray<double,1> LUI;
   TArray<double,1> LIU;
   TArray<double,1> LII;

   TArray<double,1> tmp;

   LUU.resize(dim);
   LUI.resize(dim);

   LIU.resize(dim);
   LII.resize(dim);

   blas::copy(dim, U[myID][0].data(), 1, LUU.data(), 1);
   blas::copy(dim, I[myID][0].data(), 1, LUI.data(), 1);

   blas::copy(dim, I[myID][0].data(), 1, LIU.data(), 1);
   blas::copy(dim, U[myID][0].data(), 1, LII.data(), 1);

   //put the overlap here later!
   nn_over.push_back(0);

   //now calculate the energy!
   for(int i = 1;i < L - 1;++i){

      //contribution to the energy
      if( (*this)[i - 1] != (*this)[i]){

         m = I[myID][i].shape(0);
         n = I[myID][i].shape(1);

         tmp.resize(n);

         blas::gemv(CblasRowMajor, CblasTrans, m, n, 1.0, I[myID][i].data(), n, LUI.data(), 1, 0.0, tmp.data(), 1);
         value = blas::dot(n, tmp.data(), 1, RU[i-1].data(), 1);

         blas::gemv(CblasRowMajor, CblasTrans, m, n, 1.0, U[myID][i].data(), n, LII.data(), 1, 0.0, tmp.data(), 1);
         value += blas::dot(n, tmp.data(), 1, RI[i-1].data(), 1);

         if(fabs(value) > cutoff )
            nn_over.push_back(value);
         else
            nn_over.push_back(cutoff);

      }

      //new left sides

      //if it contributes make LUI and LII:
      if( (*this)[i] != (*this)[i + 1]){

         m = I[myID][i].shape(0);
         n = I[myID][i].shape(1);

         LUI.resize(n);
         LII.resize(n);

         blas::gemv(CblasRowMajor, CblasTrans, m, n, 1.0, I[myID][i].data(), n, LUU.data(), 1, 0.0, LUI.data(), 1);
         blas::gemv(CblasRowMajor, CblasTrans, m, n, 1.0, U[myID][i].data(), n, LIU.data(), 1, 0.0, LII.data(), 1);

      }

      //new left unity: LUU and LIU
      m = U[myID][i].shape(0);
      n = U[myID][i].shape(1);

      tmp.resize(n);
      blas::gemv(CblasRowMajor, CblasTrans, m, n, 1.0, U[myID][i].data(), n, LUU.data(), 1, 0.0, tmp.data(), 1);

      LUU.resize(n);
      blas::copy(n, tmp.data(), 1, LUU.data(), 1);

      tmp.resize(n);
      blas::gemv(CblasRowMajor, CblasTrans, m, n, 1.0, I[myID][i].data(), n, LIU.data(), 1, 0.0, tmp.data(), 1);

      LIU.resize(n);
      blas::copy(n, tmp.data(), 1, LIU.data(), 1);

   }

   //final contribution to the energy
   if( (*this)[L - 2] != (*this)[L - 1]){

      value = blas::dot(LUI.size(), LUI.data(), 1, I[myID][L-1].data(), 1);
      value += blas::dot(LII.size(), LII.data(), 1, U[myID][L-1].data(), 1);

      if(fabs(value) > cutoff )
         nn_over.push_back(value);
      else
         nn_over.push_back(cutoff);

   }

   //overlap
   value = blas::dot(LUU.size(), LUU.data(), 1, U[myID][L-1].data(), 1);
   value += blas::dot(LII.size(), LIU.data(), 1, I[myID][L-1].data(), 1);

   if( fabs(value) > cutoff )
      nn_over[0] = value;
   else
      nn_over[0] = cutoff;

   //calculate the local energy
   EL = this->pot_en();

   for(int i = 1;i < nn_over.size();++i)
      EL -= 0.5 * fabs(nn_over[i] / nn_over[0]);

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
