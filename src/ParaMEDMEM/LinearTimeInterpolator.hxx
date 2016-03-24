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

#ifndef __LINEARTIMEINTERPOLATOR_HXX__
#define __LINEARTIMEINTERPOLATOR_HXX__

#include "TimeInterpolator.hxx"

#include <map>
#include <iostream>

namespace MEDCoupling
{
  class DEC;
  
  /*!
   * Internal class, not part of the public API.
   *
   * Linear interpolation of a block of data between two given times.
   */
  class LinearTimeInterpolator : public TimeInterpolator
  {
    public:  
      LinearTimeInterpolator( double InterpPrecision=0, int nStepBefore=1,
                              int nStepAfter=1 ) ;
      virtual ~LinearTimeInterpolator();
      void doInterp( double time0, double time1, double time, int recvcount,
                     int nbuff0, int nbuff1,
                     int **recvbuff0, int **recvbuff1, int *result );
      void doInterp( double time0, double time1, double time, int recvcount,
                     int nbuff0, int nbuff1,
                     double **recvbuff0, double **recvbuff1, double *result );
  };
}

#endif
