// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#ifndef __TIMEINTERPOLATOR_HXX__
#define __TIMEINTERPOLATOR_HXX__

#include "ProcessorGroup.hxx"

#include <map>
#include <iostream>

namespace MEDCoupling
{

  /*!
   * Internal class, not part of the public API.
   *
   * Abstract class for all time-related interpolation in a parallel context.
   */
  class TimeInterpolator
  {
  public:  
    TimeInterpolator( double InterpPrecision, int nStepBefore=1, int nStepAfter=1 );
    virtual ~TimeInterpolator();

    void setInterpParams( double InterpPrecision, int nStepBefore=1, int nStepAfter=1 ) { _interp_precision=InterpPrecision; _n_step_before=nStepBefore; _n_step_after=nStepAfter; }
    void steps( int &nStepBefore, int &nStepAfter ) { nStepBefore=_n_step_before; nStepAfter=_n_step_after ; }
    virtual void doInterp( double time0, double time1, double time, int recvcount ,
                           int nbuff0, int nbuff1,
                           int **recvbuff0, int **recvbuff1, int *result ) = 0;
    virtual void doInterp( double time0, double time1, double time, int recvcount ,
                           int nbuff0, int nbuff1,
                           double **recvbuff0, double **recvbuff1, double *result ) = 0;
  protected :
    double _interp_precision;
    int _n_step_before;
    int _n_step_after;
  };
}

#endif
