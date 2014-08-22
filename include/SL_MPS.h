#ifndef SL_MPS_H
#define SL_MPS_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

using namespace btas;

#include "Walker.h"

/**
 * @author Brecht Verstichel
 * @date 26-03-2014\n\n
 * This class SL_MPS is a class written for the construction of single-layer matrix product states
 * i.e. MPS contracted with a walker
 */
class SL_MPS : public vector< DArray<2> > {

   public:

      SL_MPS();

      SL_MPS(const MPS<double> &,const Walker &);

      //copy constructor
      SL_MPS(const SL_MPS &);

      //destructor
      virtual ~SL_MPS();

      void fill(bool,const MPS<double> &,const Walker &);

   private:

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
