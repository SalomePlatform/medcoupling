#ifndef TIMEINTERPOLATOR_HXX_
#define TIMEINTERPOLATOR_HXX_

#include <map>
#include <iostream>

#include "ProcessorGroup.hxx"

namespace ParaMEDMEM {

  class TimeInterpolator {

    public:  
      TimeInterpolator( double InterpPrecision, int nStepBefore=1, int nStepAfter=1 ) ;
      virtual ~TimeInterpolator();

      void SetInterpParams( double InterpPrecision, int nStepBefore=1, int nStepAfter=1 ) {
           _InterpPrecision = InterpPrecision ;
           _nStepBefore = nStepBefore ;
           _nStepAfter = nStepAfter ; } ;
      void Steps( int &nStepBefore, int &nStepAfter ) {
           nStepBefore = _nStepBefore ;
           nStepAfter = _nStepAfter ; } ;
      virtual void DoInterp( double time0, double time1, double time, int recvcount ,
                             int nbuff0, int nbuff1,
                             int **recvbuff0, int **recvbuff1, int *result )= 0 ;
      virtual void DoInterp( double time0, double time1, double time, int recvcount ,
                             int nbuff0, int nbuff1,
                             double **recvbuff0, double **recvbuff1, double *result )= 0 ;

    protected :
      double     _InterpPrecision ;
      int        _nStepBefore ;
      int        _nStepAfter ;
  };
}

#endif
