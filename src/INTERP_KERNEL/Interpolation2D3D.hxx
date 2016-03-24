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

#ifndef __INTERPOLATION2D3D_HXX__
#define __INTERPOLATION2D3D_HXX__

#include <set>
#include <map>

#include "INTERPKERNELDefines.hxx"
#include "Interpolation.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpolationOptions.hxx"

namespace INTERP_KERNEL
{
  /*!
   * Contrary to its name this class deals with 2D mesh in source and 3D mesh in target.
   * The meshdim of 'MyMeshType' in input is ignored that's why 'meshS' and 'meshT'
   * have the same type.
   * '_duplicate_faces' attribute stores duplicated faces in the following format.
   * The key of '_duplicate_faces' represents the 2D cellId that is shared by
   * more than one 3D target cell, and the value of '_duplicate_faces'
   * the 3D target cells. The size of the value of '_duplicate_faces' is more than or equal to 2.
   */
  class Interpolation2D3D : public Interpolation<Interpolation2D3D>
  {
  public:
    typedef std::map<int,std::set<int> > DuplicateFacesType;

    INTERPKERNEL_EXPORT Interpolation2D3D();
    INTERPKERNEL_EXPORT Interpolation2D3D(const InterpolationOptions& io);
    template<class MyMeshType, class MyMatrixType>
    int interpolateMeshes(const MyMeshType& srcMesh,
                          const MyMeshType& targetMesh,
                          MyMatrixType& matrix,
                          const std::string& method);
    INTERPKERNEL_EXPORT DuplicateFacesType retrieveDuplicateFaces() const { return _duplicate_faces; }
  private:
    SplittingPolicy _splitting_policy;
    DuplicateFacesType _duplicate_faces;
  };
}

#endif
