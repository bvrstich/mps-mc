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
 * constructor of the VMC object, takes input parameters that define the QMC walk.
 * @param Nw_in number of Walker states
 */
VMC::VMC(int Nw_in){

   this->Nw = Nw_in;

   dist.resize(omp_num_threads);

   SetupWalkers();

}

/**
 * unnecessary destructor...
 */
VMC::~VMC(){ }

/**
 * initialize the walkers: read in distribution from walkers dir
 */
void VMC::SetupWalkers(){

   walker.resize(Nw);

   walker[0].calc_EL(mps);

   walker[1] = Walker(1);
   walker[1].calc_EL(mps);

   for(int i = 2;i < Nw;++i){

      if(i % 2 == 0)
         walker[i] = walker[0];
      else
         walker[i] = walker[1];

   }

}

void VMC::walk(const int n_steps){

   //set projected energy
   sEP();

   char filename[200];
   sprintf(filename,"output/VMC/L=%d/D=%d.txt",L,DT);

   ofstream output(filename,ios::trunc);

   output << "#Step\t\tE_P\t\tE_T\t" << endl;
   output.close();

#ifdef _DEBUG
   cout << "Energy at start = " << EP << "\t" << endl;
   cout << "---------------------------------------------------------" << endl;
#endif

   for(int step = 0;step < n_steps;step++){

      //Propagate the walkers of each rank separately
      propagate();

      //calculate the energy
      sEP();

#ifdef _DEBUG
      cout << "        Step = " << step << endl;
      cout << "   # walkers = " << walker.size() << endl;
      cout << "         E_P = " << EP << endl;
      cout << "---------------------------------------------------------" << endl;
#endif

      write(step);

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
void VMC::propagate(){

#pragma omp parallel for
   for(unsigned int i = 0;i < walker.size();i++){

#ifdef _OPENMP
      int myID = omp_get_thread_num();
#else
      int myID = 0;
#endif

      //construct distribution
      dist[myID].construct(walker[i]);

      //draw new walker
      int pick = dist[myID].metropolis();

      walker[i] = dist[myID].gwalker(pick);

      //calculate new properties
      walker[i].calc_EL(mps);

   }

}

/**
 * set total projected energy of the walkers at a certain timestep
 */
void VMC::sEP(){

   double projE_num = 0.0;
   double projE_den = 0.0;

   for(unsigned int wi = 0;wi < walker.size();wi++){

      double w_loc_en = walker[wi].gEL(); // <Psi_T | H | walk > / <Psi_T | walk >

      //For the projected energy
      projE_num += walker[wi].gWeight() * w_loc_en;
      projE_den += walker[wi].gWeight();

   }

   EP = projE_num / projE_den;

}

/**
 * write output to file
 */
void VMC::write(const int step){

   char filename[200];
   sprintf(filename,"output/VMC/L=%d/D=%d.txt",L,DT);

   ofstream output(filename,ios::app);
   output.precision(16);
   output << step << "\t\t" << walker.size() << "\t" << EP << endl;
   output.close();

}


/**
 * dump the walkers to a single file
 */
void VMC::dump(const char *filename){

   ofstream out(filename);
   out.precision(16);

   for(unsigned int i = 0;i < walker.size();++i)
      for(int j = 0;j < L;++j)
         out << walker[i][j] << " ";
   out << endl;

}
