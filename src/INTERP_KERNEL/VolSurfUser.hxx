// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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
// Author : Anthony Geay (CEA/DEN)

#ifndef __VOLSURFUSER_HXX__
#define __VOLSURFUSER_HXX__

#include "NormalizedUnstructuredMesh.hxx"

namespace INTERP_KERNEL
{
  template<class ConnType, NumberingPolicy numPolConn, int SPACEDIM>
  double computeVolSurfOfCell(NormalizedCellType type, const ConnType *connec, int lgth, const double *coords);

  template<class ConnType, NumberingPolicy numPolConn>
  double computeVolSurfOfCell2(NormalizedCellType type, const ConnType *connec, int lgth, const double *coords, int spaceDim);

  template<class ConnType, NumberingPolicy numPolConn, int SPACEDIM>
  void computeBarycenter(NormalizedCellType type, const ConnType *connec, int lgth, const double *coords, double *res);

  template<class ConnType, NumberingPolicy numPolConn>
  void computeBarycenter2(NormalizedCellType type, const ConnType *connec, int lgth, const double *coords, int spaceDim, double *res);
}

#endif
