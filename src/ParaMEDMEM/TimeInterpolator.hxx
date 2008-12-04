//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
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
