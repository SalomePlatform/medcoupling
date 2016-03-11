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
// Author : Anthony Geay (CEA/DEN)

#ifndef __INTERPOLATION3DSURF_HXX__
#define __INTERPOLATION3DSURF_HXX__

#include "InterpolationPlanar.txx"
#include "INTERPKERNELDefines.hxx"
#include "InterpolationOptions.hxx"

namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT Interpolation3DSurf : public InterpolationPlanar<Interpolation3DSurf>
  {
  public:
    Interpolation3DSurf();
    Interpolation3DSurf(const InterpolationOptions& io);
    void setOptions(double precision, int printLevel, double medianPlane, 
                    IntersectionType intersectionType, bool doRotate, int orientation=0);
  public:
    template<class MyMeshType, class MyMatrixRow>
    void performAdjustmentOfBB(PlanarIntersector<MyMeshType,MyMatrixRow>* intersector, std::vector<double>& bbox) const
    { intersector->adjustBoundingBoxes(bbox,InterpolationPlanar<Interpolation3DSurf>::getBoundingBoxAdjustment(),InterpolationPlanar<Interpolation3DSurf>::getBoundingBoxAdjustmentAbs()); }
  };
}

#endif
