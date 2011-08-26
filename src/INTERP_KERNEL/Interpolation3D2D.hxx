// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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

#ifndef __INTERPOLATION3D2D_HXX__
#define __INTERPOLATION3D2D_HXX__

#include <set>
#include <map>

#include "Interpolation.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpolationOptions.hxx"

namespace INTERP_KERNEL
{
  class Interpolation3D2D : public Interpolation<Interpolation3D2D>
  {
  public:
    typedef std::map<int,std::set<int> > DuplicateFacesType;

    Interpolation3D2D();
    Interpolation3D2D(const InterpolationOptions& io);
    template<class MyMeshType, class MyMatrixType>
    int interpolateMeshes(const MyMeshType& srcMesh,
                          const MyMeshType& targetMesh,
                          MyMatrixType& matrix,
                          const char *method);
    DuplicateFacesType retrieveDuplicateFaces() const
    {
      return _duplicate_faces;
    }
  private:
    SplittingPolicy _splitting_policy;
    DuplicateFacesType _duplicate_faces;
  };
}

#endif
