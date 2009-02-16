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
#ifndef __INTERPOLATION3DSURF_TXX__
#define __INTERPOLATION3DSURF_TXX__

#include "Interpolation3DSurf.hxx"
#include "InterpolationPlanar.txx"

namespace INTERP_KERNEL
{
  const double Interpolation3DSurf::DFT_MEDIAN_PLANE=0.5;
  const double Interpolation3DSurf::DFT_SURF3D_ADJ_EPS=1e-4;
  
  Interpolation3DSurf::Interpolation3DSurf():_do_rotate(true)
                                            ,_median_plane(DFT_MEDIAN_PLANE)
                                            ,_surf_3D_adjustment_eps(DFT_SURF3D_ADJ_EPS)
  {
  }

  Interpolation3DSurf::Interpolation3DSurf(const InterpolationOptions& io):InterpolationPlanar<Interpolation3DSurf>(io)
  {
  }

  
  /**
     \brief  Function used to set the options for the intersection calculation
     \details The following options can be modified:
     -# Intersection_type: the type of algorithm to be used in the computation of the cell-cell intersections.
     - Values: Triangle, Convex.
     - Default: Triangle.
     -# MedianPlane: Position of the median plane where both cells will be projected
     - Values: between 0 and 1.
     - Default: 0.5.
     -# DoRotate: rotate the coordinate system such that the target cell is in the Oxy plane.
     - Values: true (necessarilly if Intersection_type=Triangle), false.
     - Default: true (as default Intersection_type=Triangle)
     -# Precision: Level of precision of the computations is precision times the characteristic size of the mesh.
     - Values: positive real number.
     - Default: 1.0E-12.
     -# PrintLevel: Level of verboseness during the computations.
     - Values: interger between 0 and 3.
     - Default: 0.
  */
  void Interpolation3DSurf::setOptions(double precision, int printLevel, double medianPlane, 
                                       IntersectionType intersectionType, bool doRotate, int orientation)
  {
    InterpolationPlanar<Interpolation3DSurf>::setOptions(precision,printLevel,intersectionType, orientation);
    _do_rotate=doRotate;
    _median_plane=medianPlane;
  }
}

#endif
