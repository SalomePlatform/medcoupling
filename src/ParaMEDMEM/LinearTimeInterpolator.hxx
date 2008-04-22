#ifndef LINEARTIMEINTERPOLATOR_HXX_
#define LINEARTIMEINTERPOLATOR_HXX_

#include <map>
#include <iostream>

#include "TimeInterpolator.hxx"

namespace ParaMEDMEM {

  class DEC;

  class LinearTimeInterpolator:public TimeInterpolator {

    public:  
      LinearTimeInterpolator( double InterpPrecision=0, int nStepBefore=1,
                              int nStepAfter=1 ) ;
      virtual ~LinearTimeInterpolator();

      void DoInterp( double time0, double time1, double time, int recvcount,
                     int nbuff0, int nbuff1,
                     int **recvbuff0, int **recvbuff1, int *result ) ;
      void DoInterp( double time0, double time1, double time, int recvcount,
                     int nbuff0, int nbuff1,
                     double **recvbuff0, double **recvbuff1, double *result ) ;

    private :
  };
}

#endif
