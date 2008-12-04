//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D, OPEN CASCADE
//
//  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
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
#ifndef _INTERPOLATION_HXX_
#define _INTERPOLATION_HXX_

/**
 * \mainpage
 * Status : documentation of 3D - part of intersection matrix calculation more or less complete
 *
 *
 */
#include "INTERPKERNEL_defines.hxx"
#include "InterpolationOptions.hxx"

namespace INTERP_KERNEL
{
  template<class TrueMainInterpolator>
  class INTERPKERNEL_EXPORT Interpolation : public InterpolationOptions
  {
  public:
    Interpolation() { 
	// 	double InterpolationOptions::_precision=1e-12;
// 		int InterpolationOptions::_printLevel=0;
// 		InterpolationOptions::setIntersectionType(Triangulation);
// 		InterpolationOptions::setMedianPlane(0.5);
// 		InterpolationOptions::setDoRotate(true);
// 		InterpolationOptions::setBoundingBoxAdjustment(0.1);
// 		InterpolationOptions::setSplittingPolicy(GENERAL_48);
		}
		Interpolation(const InterpolationOptions& io) :InterpolationOptions(io){}
    //interpolation of two triangular meshes.
    template<class MatrixType, class MyMeshType>
    void interpolateMeshes(const MyMeshType& mesh1, const MyMeshType& mesh2, MatrixType& result)
    { return asLeaf().interpolateMeshes(mesh1,mesh2,result); }
  protected:
    TrueMainInterpolator& asLeaf() { return static_cast<TrueMainInterpolator&>(*this); }
  };
}

#endif
