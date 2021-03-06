#ifndef MPS_H
#define MPS_H

#include <iostream>
#include <fstream>
#include <vector>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>
#include <btas/DENSE/TArray.h>

using std::ostream;
using std::vector;

using namespace btas;

class Walker;

/**
 * @author Brecht Verstichel
 * @date 26-03-2014\n\n
 * This class MPS is a class written for the construction of matrix product states without symmetry.
 */
template<typename T>
class MPS : public vector< TArray<T,3> > {

   public:

      MPS();

      MPS(int D);

      MPS(double f);

      MPS(const Walker &);

      //copy constructor
      MPS(const MPS &);

      //destructor
      virtual ~MPS();

      int gD() const;

      void load(const char *);

      double dot(const MPS<double> &) const;

      double energy() const;

      void scal(double);

   private:

      //!dimension of the bonds
      int D;

};

/**
 * output stream operator overloaded for MPS<T> 
 */
template<typename T>
ostream &operator<<(ostream &output,const MPS<T> &mps_p){

   for(int s = 0;s < mps_p.size();++s){

         output << std::endl;
         output << "Tensor on site (" << s << ")\t" << std::endl;
         output << std::endl;
         output << mps_p[s] << std::endl;

      }

   return output;

}

#endif

/* vim: set ts=3 sw=3 expandtab :*/
