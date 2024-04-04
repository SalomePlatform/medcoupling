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

#ifndef __INTERPOLATIONCURVE_HXX__
#define __INTERPOLATIONCURVE_HXX__

#include "Interpolation.hxx"
#include "InterpolationOptions.hxx"

#include "BBTree.txx"

#include <functional>
#include <string>
#include <vector>

namespace INTERP_KERNEL
{
  /**
   * \defgroup interpolationCurve InterpolationCurve
   *
   * \class InterpolationCurve
   * \brief Class used to compute the coefficients of the interpolation matrix between
   * two local meshes in two dimensions.
   */

  template<class RealCurve>
  class InterpolationCurve : public Interpolation< InterpolationCurve<RealCurve> >
  {
  public:
    InterpolationCurve();
    InterpolationCurve(const InterpolationOptions & io);

    // Main function to interpolate
    template<class MyMeshType, class MatrixType>
    typename MyMeshType::MyConnType interpolateMeshesInternal(const MyMeshType& meshS, const MyMeshType& meshT,
                                                              MatrixType& result, const std::string& method,
                                                              std::function< void(const BBTree< MyMeshType::MY_SPACEDIM , typename MyMeshType::MyConnType>&, const double*, std::vector<typename MyMeshType::MyConnType>&) > bbtreeMethod);
    template<class MyMeshType, class MatrixType>
    typename MyMeshType::MyConnType interpolateMeshes(const MyMeshType& meshS, const MyMeshType& meshT,
                                                      MatrixType& result, const std::string& method)
    {
      std::function< void(const BBTree< MyMeshType::MY_SPACEDIM , typename MyMeshType::MyConnType>&, const double*, std::vector<typename MyMeshType::MyConnType>&) > bbtreeMethod =
        [](const BBTree< MyMeshType::MY_SPACEDIM , typename MyMeshType::MyConnType>& bbtree, const double *bb, std::vector<typename MyMeshType::MyConnType>& intersecting_elems)
        { bbtree.getIntersectingElems(bb, intersecting_elems); };
      return this->interpolateMeshesInternal(meshS,meshT,result,method,bbtreeMethod);
    }
    
    template<class MyMeshType, class MatrixType>
    typename MyMeshType::MyConnType interpolateMeshes0D(const MyMeshType& meshS, const MyMeshType& meshT,
                                                        MatrixType& result, const std::string& method)
    {
      std::function< void(const BBTree< MyMeshType::MY_SPACEDIM , typename MyMeshType::MyConnType>&, const double*, std::vector<typename MyMeshType::MyConnType>&) > bbtreeMethod =
        [](const BBTree< MyMeshType::MY_SPACEDIM , typename MyMeshType::MyConnType>& bbtree, const double *bb, std::vector<typename MyMeshType::MyConnType>& intersecting_elems)
        {
          double TMP[MyMeshType::MY_SPACEDIM];
          for(int i=0;i<MyMeshType::MY_SPACEDIM;++i)
            TMP[i] = bb[2*i];
          bbtree.getElementsAroundPoint(TMP, intersecting_elems);
        };
      return this->interpolateMeshesInternal(meshS,meshT,result,method,bbtreeMethod);
    }
  };
}

#endif
