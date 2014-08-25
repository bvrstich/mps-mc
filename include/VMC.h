#ifndef VMC_H
#define VMC_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using std::complex;
using std::vector;

class VMC {

   public:
   
      //constructor with input trialwavefunction
      VMC(int);
      
      //Destructor
      virtual ~VMC();
      
      //Let the walkers propagate for n_steps steps
      void walk(int);

      //Propagate my population of walkers for 1 timestep. Return the sum of the coeff of my walkers.
      void propagate();
      
      //Calculate the single walker projected energies
      void sEP();

      //Write the projected energy, target energy
      void write(int);

      //Setup the walkers
      void SetupWalkers();

   private:
      
      //!The total desired number of walkers
      int Nw;

      //!projected energy at current timestep
      double EP;

      double cum_num;
      double cum_den;

      //!The walkers
      std::vector<Walker> walker;

      //!Distribution of possible final states given a walker state
      std::vector< Distribution > dist;
      
};

#endif
