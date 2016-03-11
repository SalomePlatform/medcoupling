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

#ifndef __INTERPOLATION2D1D_HXX__
#define __INTERPOLATION2D1D_HXX__

#include "Interpolation.hxx"
#include "Planar2D1DIntersectorP0P0.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpolationOptions.hxx"

namespace INTERP_KERNEL
{
  /*!
   * Contrary to its name this class deals with 1D mesh in source and 2D mesh in target.
   * The meshdim of 'MyMeshType' in input is ignored that's why 'meshS' and 'meshT'
   * have the same type.
   * '_duplicate_faces' attribute stores duplicated faces in the following format.
   * The key of '_duplicate_faces' represents the 1D cellId that is shared by
   * more than one 2D target cell, and the value of '_duplicate_faces'
   * the 2D target cells. The size of the value of '_duplicate_faces' is more than or equal to 2.
   */
  class Interpolation2D1D : public Interpolation<Interpolation2D1D>
  {
  public:
    typedef std::map<int,std::set<int> > DuplicateFacesType;

    Interpolation2D1D() { setOrientation(2); }
    Interpolation2D1D(const InterpolationOptions& io):Interpolation<Interpolation2D1D>(io) { }
  public:

    // Main function to interpolate triangular and quadratic meshes
    template<class MyMeshType, class MatrixType>
    int interpolateMeshes(const MyMeshType& meshS, const MyMeshType& meshT, MatrixType& result, const std::string& method);
    DuplicateFacesType retrieveDuplicateFaces() const
    {
      return _duplicate_faces;
    }
  private:
    DuplicateFacesType _duplicate_faces;
  private:
    double _dim_caracteristic;
  };
}

#endif
