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
#ifndef __INTERPOLATION_HXX__
#define __INTERPOLATION_HXX__

/**
 * \mainpage
 * Status : documentation of 3D - part of intersection matrix calculation more or less complete
 *
 *
 */
#include "INTERPKERNELDefines.hxx"
#include "InterpolationOptions.hxx"

namespace INTERP_KERNEL
{
  template<class TrueMainInterpolator>
  class Interpolation : public InterpolationOptions
  {
  public:
    Interpolation() { }
    Interpolation(const InterpolationOptions& io) :InterpolationOptions(io){}
    //interpolation of two triangular meshes.
    template<class MatrixType, class MyMeshType>
    int interpolateMeshes(const MyMeshType& mesh1, const MyMeshType& mesh2, MatrixType& result)
    { return asLeaf().interpolateMeshes(mesh1,mesh2,result); }
  protected:
    TrueMainInterpolator& asLeaf() { return static_cast<TrueMainInterpolator&>(*this); }
  };
}

#endif
