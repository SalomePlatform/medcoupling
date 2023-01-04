// Copyright (C) 2022  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (EDF R&D)

#pragma once

#include "BBTreeStandAlone.txx"
#include "MeshElement.txx"
#include "Log.hxx"

#include <memory>
#include <functional>

namespace INTERP_KERNEL
{
  template<class MyMeshType>
  BBTreeStandAlone<3,typename MyMeshType::MyConnType> BuildBBTreeWithAdjustment(const MyMeshType& srcMesh, std::function<void(double *,typename MyMeshType::MyConnType)> bboxAdjuster)
  {
    using ConnType = typename MyMeshType::MyConnType;
    const ConnType numSrcElems = srcMesh.getNumberOfElements();
    LOG(2, "Source mesh has " << numSrcElems << " elements");
    // create BBTree structure
    // - get bounding boxes
    const ConnType nbElts = 6 * numSrcElems;
    std::unique_ptr<double[]> bboxes( new double[nbElts] );
    for(ConnType i = 0; i < numSrcElems ; ++i)
      {
        MeshElement<ConnType> srcElem(i,srcMesh);
        // get source bboxes in right order
        const BoundingBox *box( srcElem.getBoundingBox() );
        box->fillInXMinXmaxYminYmaxZminZmaxFormat(bboxes.get()+6*i);
      }
    bboxAdjuster(bboxes.get(),nbElts);
    return BBTreeStandAlone<3,ConnType>(std::move(bboxes),numSrcElems);
  }

  template<class MyMeshType>
  BBTreeStandAlone<3,typename MyMeshType::MyConnType> BuildBBTree(const MyMeshType& srcMesh)
  {
    return BuildBBTreeWithAdjustment(srcMesh,[](double *,typename MyMeshType::MyConnType){});
  }
}
