#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <vector>

#include "include.h"

using std::cout;
using std::endl;
using std::complex;
using std::vector;
using std::ofstream;
using std::ifstream;
using std::ios;

using namespace global;

/**
 * constructor of the GFMC object, takes input parameters that define the QMC walk.
 * @param dtau_in timestep
 * @param Nw_in number of Walker states
 */
GFMC::GFMC(double dtau_in,int Nw_in){

   this->dtau = dtau_in;
   this->Nw = Nw_in;

   dist.resize(omp_num_threads);

   SetupWalkers();

}

/**
 * unnecessary destructor...
 */
GFMC::~GFMC(){ }

/**
 * initialize the walkers: read in distribution from walkers dir
 */
void GFMC::SetupWalkers(){

   walker.resize(Nw);

   char filename[200];
   sprintf(filename,"output/VMC/L=%d/D=%d.walk",L,DT);

   ifstream input(filename);

   for(int i = 0;i < Nw;++i){

      for(unsigned int j = 0;j < walker[i].size();++j){

         bool tmp;

         input >> tmp;
         walker[i][j] = tmp;

      }

      walker[i].calc_EL(global::mps);
      walker[i].calc_overlap_jastrow();

   }

}

void GFMC::walk(const int n_steps){

   //set projected energy
   sEP();

   char filename[200];
   sprintf(filename,"output/GFMC/L=%d/D=%d.txt",L,DT);

   ofstream output(filename,ios::trunc);

   output << "#Step\t\tE_P\t\tE_T\t" << endl;
   output.close();

#ifdef _DEBUG
   cout << "Energy at start = " << EP << "\t" << endl;
   cout << "---------------------------------------------------------" << endl;
#endif

   for(int step = 0;step < n_steps;step++){

      //Propagate the walkers of each rank separately
      double wsum = propagate();

      double scaling = Nw / wsum;

      ET = 1.0/dtau * ( 1 - 1.0/scaling);

      //calculate the energy
      sEP();

#ifdef _DEBUG
      cout << "        Step = " << step << endl;
      cout << "   # walkers = " << walker.size() << endl;
      cout << "         E_P = " << EP << endl;
      cout << "         E_T = " << ET << endl;
      cout << "---------------------------------------------------------" << endl;
#endif

      write(step);

      //Based on scaling, first control the population on each rank separately, and then balance the population over the ranks (uses MPI)
      PopulationControl(scaling);

      double max_ov = 0.0;
      double min_ov = 1.0;

      double max_en = -100.0;
      double min_en = 100.0;

      for(unsigned int i = 0;i < walker.size();++i){

         if(max_ov < std::abs(walker[i].gOverlap()))
            max_ov = std::abs(walker[i].gOverlap());

         if(min_ov > std::abs(walker[i].gOverlap()))
            min_ov = std::abs(walker[i].gOverlap());

         if(max_en < walker[i].gEL())
            max_en = walker[i].gEL();

         if(min_en > walker[i].gEL())
            min_en = walker[i].gEL();

      }

#ifdef _DEBUG
      cout << endl;
      cout << "Minimal Overlap:\t" << min_ov << endl;
      cout << "Maximal Overlap:\t" << max_ov << endl;
      cout << endl;
      cout << "Minimal Energy:\t" << min_en << endl;
      cout << "Maximal Energy:\t" << max_en << endl;
      cout << endl;
#endif

   }

}

/**
 * Here the trotter terms, propagator terms are applied to every walker individually.
 */
double GFMC::propagate(){

   double total_weight = 0.0;

#pragma omp parallel for reduction(+:total_weight)
   for(unsigned int i = 0;i < walker.size();i++){

#ifdef _OPENMP
      int myID = omp_get_thread_num();
#else
      int myID = 0;
#endif

      //construct distribution
      dist[myID].construct(walker[i],dtau,0.0);
      double nrm = dist[myID].normalize();

      //draw new walker
      int pick = dist[myID].draw();

      //copy the correct walker
      walker[i] = dist[myID].gwalker(pick);

      //multiply weight
      walker[i].multWeight(nrm);

      walker[i].calc_EL(global::mps);

      total_weight += walker[i].gWeight();

   }

   return total_weight;

}

/**
 * set total projected energy of the walkers at a certain timestep
 */
void GFMC::sEP(){

   double projE_num = 0.0;
   double projE_den = 0.0;

   double max_ratio = 0.0;
   double min_ratio = 100;

   for(unsigned int wi = 0;wi < walker.size();wi++){

      double over_ratio = fabs(walker[wi].gnn_over(0)/walker[wi].gOverlap_jastrow());

#ifdef _DEBUG
      if(over_ratio > max_ratio)
         max_ratio = over_ratio;

      if(over_ratio < min_ratio)
         min_ratio = over_ratio;
#endif

      double w_loc_en = walker[wi].gEL(); // <Psi_T | H | walk > / <Psi_T | walk >

      //For the projected energy
      projE_num += over_ratio * walker[wi].gWeight() * w_loc_en;
      projE_den += over_ratio * walker[wi].gWeight();

   }

#ifdef _DEBUG
   cout << endl;
   cout << "Max overlap ratio =\t" << max_ratio << endl;
   cout << "Min overlap ratio =\t" << min_ratio << endl;
   cout << endl;
#endif

   EP = projE_num / projE_den;

}

/**
 * write output to file
 */
void GFMC::write(const int step){

   char filename[200];
   sprintf(filename,"output/GFMC/L=%d/D=%d.txt",L,DT);

   ofstream output(filename,ios::app);
   output.precision(16);
   output << step << "\t\t" << walker.size() << "\t" << EP << "\t" << ET << endl;
   output.close();

}

/**
 * redistribute the weights to stabilize the walk, keep the population in check
 */
void GFMC::PopulationControl(double scaling){

   double minw = 1.0;
   double maxw = 1.0;

   for(unsigned int i = 0;i < walker.size();i++)
      walker[i].multWeight(scaling);

   for(unsigned int i = 0;i < walker.size();i++){

      double weight = walker[i].gWeight();

      if(weight < minw)
         minw = weight;

      if(weight > maxw)
         maxw = weight;

      if (weight < 0.1){ //Energy doesn't change statistically

         int nCopies = (int) ( weight + RN());

         if(nCopies == 0){

            walker.erase(walker.begin() + i);

         }
         else
            walker[i].sWeight(1.0);

      }

      if(weight > 2.0){ //statically energy doesn't change

         int nCopies =(int) ( weight + RN());
         double new_weight = weight / (double) nCopies;

         walker[i].sWeight(new_weight);

         for(int n = 1;n < nCopies;++n){

            Walker nw = walker[i];
            walker.push_back(nw);

         }

      }

   }

   double sum = 0.0;

   for(unsigned int i = 0;i < walker.size();++i)
      sum += walker[i].gWeight();
   
   //rescale the weights to unity for correct ET estimate in next iteration
   for(unsigned int i = 0;i < walker.size();++i)
      walker[i].multWeight((double)Nw/sum);


#ifdef _DEBUG
   cout << endl;
   cout << "total weight:\t" << sum << endl;
   cout << endl;

   cout << "The min. encountered weight is " << minw << " ." << endl;
   cout << "The max. encountered weight is " << maxw << " ." << endl;
   cout << endl;
#endif

}
