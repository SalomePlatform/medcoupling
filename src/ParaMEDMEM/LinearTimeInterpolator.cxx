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

#include "LinearTimeInterpolator.hxx"

using namespace std;

namespace MEDCoupling
{    

  LinearTimeInterpolator::LinearTimeInterpolator( double InterpPrecision, int nStepBefore,
                                                  int nStepAfter ):
    TimeInterpolator( InterpPrecision, nStepBefore, nStepAfter )
  {
  }
  
  LinearTimeInterpolator::~LinearTimeInterpolator()
  {
  } 
  
  void LinearTimeInterpolator::doInterp( double time0, double time1, double time,
                                         int recvcount , int nbuff0, int nbuff1,
                                         int **recvbuff0, int **recvbuff1, int *result )
  {
    for(int i = 0 ; i < recvcount ; i++ )
        result[i] = (int) ((recvbuff0[0][i]*(time1 - time) + recvbuff1[0][i]*(time - time0))/(time1 - time0) + _interp_precision);
  }

  void LinearTimeInterpolator::doInterp( double time0, double time1, double time,
                                         int recvcount , int nbuff0, int nbuff1,
                                         double **recvbuff0, double **recvbuff1,
                                         double *result )
  {
    for(int i = 0 ; i < recvcount ; i++ )
      result[i] = (recvbuff0[0][i]*(time1 - time) + recvbuff1[0][i]*(time - time0))/(time1 - time0);
  }
  
}
