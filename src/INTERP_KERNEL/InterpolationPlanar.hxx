// Copyright (C) 2007-2024  CEA, EDF
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

#ifndef __INTERPOLATIONPLANAR_HXX__
#define __INTERPOLATIONPLANAR_HXX__

#include "Interpolation.hxx"
#include "PlanarIntersector.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpolationOptions.hxx"

namespace INTERP_KERNEL
{
  /**
   * \defgroup interpolationPlanar InterpolationPlanar
   *
   * \class InterpolationPlanar
   * \brief Class used to compute the coefficients of the interpolation matrix between
   * two local meshes in two dimensions. Meshes can contain mixed triangular and quadrangular elements.
   */

  template<class RealPlanar>
  class InterpolationPlanar : public Interpolation< InterpolationPlanar<RealPlanar> >
  {
  private:
    double _dim_caracteristic;
  public:
    InterpolationPlanar();
    InterpolationPlanar(const InterpolationOptions & io);

    // geometric precision, debug print level, choice of the median plane, intersection etc ...
    void setOptions(double precision, int printLevel,
                    IntersectionType intersectionType, int orientation=0);

    // Main function to interpolate triangular and quadratic meshes
    template<class MyMeshType, class MatrixType>
    typename MyMeshType::MyConnType interpolateMeshes(const MyMeshType& meshS, const MyMeshType& meshT, MatrixType& result, const std::string& method);
  public:
    bool doRotate() const { return asLeafInterpPlanar().doRotate(); }
    double medianPlane() const { return asLeafInterpPlanar().medianPlane(); }
    template<class MyMeshType, class MyMatrixRow>
      void performAdjustmentOfBB(PlanarIntersector<MyMeshType,MyMatrixRow>* intersector, std::vector<double>& bbox) const
    { return asLeafInterpPlanar().performAdjustmentOfBB(intersector,bbox); }
  protected:
    RealPlanar& asLeafInterpPlanar() { return static_cast<RealPlanar&>(*this); }
    const RealPlanar& asLeafInterpPlanar() const { return static_cast< const RealPlanar& >(*this); }
  };
}

#endif
