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

#include "Interpolation3DSurf.hxx"
#include "InterpolationPlanar.txx"

namespace INTERP_KERNEL
{
  Interpolation3DSurf::Interpolation3DSurf()
  {
  }

  Interpolation3DSurf::Interpolation3DSurf(const InterpolationOptions& io):InterpolationPlanar<Interpolation3DSurf>(io)
  {
  }

  
  /**
     \brief  Function used to set the options for the intersection calculation
     \details The following options can be modified:
     -# intersectionType: the type of algorithm to be used in the computation of the cell-cell intersections.
     - Values: Triangle, Convex.
     - Default: Triangle.
     -# medianPlan: Position of the median plane where both cells will be projected
     - Values: between 0 and 1.
     - Default: 0.5.
     -# doRotat: rotate the coordinate system such that the target cell is in the Oxy plane.
     - Values: true (necessarilly if Intersection_type=Triangle), false.
     - Default: true (as default Intersection_type=Triangle)
     -# precision: Level of precision of the computations is precision times the characteristic size of the mesh.
     - Values: positive real number.
     - Default: 1.0E-12.
     -# printLevel: Level of verboseness during the computations.
     - Values: interger between 0 and 3.
     - Default: 0.
  */
  void Interpolation3DSurf::setOptions(double precision, int printLevel, double medianPlan, 
                                       IntersectionType intersectionType, bool doRotat, int orientation)
  {
    InterpolationPlanar<Interpolation3DSurf>::setOptions(precision,printLevel,intersectionType, orientation);
    InterpolationPlanar<Interpolation3DSurf>::setDoRotate(doRotat);
    InterpolationPlanar<Interpolation3DSurf>::setMedianPlane(medianPlan);
  }
}
