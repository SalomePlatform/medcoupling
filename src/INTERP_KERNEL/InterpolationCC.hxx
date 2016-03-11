// Copyright (C) 2009-2016  OPEN CASCADE
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

// File      : InterpolationCC.hxx
// Created   : Fri Aug 14 11:33:17 2009
// Author    : Edward AGAPOV (eap)
//
#ifndef __InterpolationCC_HXX__
#define __InterpolationCC_HXX__

#include "Interpolation.hxx"

namespace INTERP_KERNEL
{
  /*!
   * \brief Interpolator of cartesian/cartesian meshes
   */
  class InterpolationCC : public Interpolation<InterpolationCC>
  {
//     static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
//     static const int MESHDIM=MyMeshType::MY_MESHDIM;
//     typedef typename MyMeshType::MyConnType ConnType;
//     static const NumberingPolicy numPol=MyMeshType::My_numPol;
  public:
    InterpolationCC();
    InterpolationCC(const InterpolationOptions& io);
    template<class MyMeshType, class MatrixType>
    int interpolateMeshes(const MyMeshType& srcMesh, const MyMeshType& targetMesh, MatrixType& result, const char *method);

  private:
  };
}


#endif
