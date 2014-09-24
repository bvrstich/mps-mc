#ifndef WALKER_H
#define WALKER_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

using std::complex;
using std::vector;

using namespace btas;

#include "MPS.h"

/**
 * class definition of Walker, made to describe the product state walkers. An array of L*L booleans vector. Each site represents a spin up (0) or spin down (1)
 */
class Walker : public vector< bool > {

   friend ostream &operator<<(ostream &output,const Walker &walker_p);

   public:

      //empty contstructor
      Walker(int option = 0);
   
      //Constructor copying an entire Walker
      Walker(const Walker &walker);
      
      //Destructor
      virtual ~Walker();

      double gWeight() const;

      void sWeight(double);

      void multWeight(double);

      double gOverlap() const;

      const vector<double> &gnn_over() const;

      double gnn_over(int) const;

      double gEL() const;

      double pot_en() const;

      void calc_EL(const MPS<double> &);

      void save(const char *filename);

      void load(const char *filename);

      int num_diff(const Walker &) const;

  private:

      //!The walker weight
      double weight;

      //!local energy
      double EL;

      //!overlap of nn configurations with trial
      vector<double> nn_over;

      //!overlap with trial
      double overlap;

};

#endif
