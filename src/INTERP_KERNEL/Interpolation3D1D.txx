// Copyright (C) 2007-2020  CEA/DEN, EDF R&D
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

#include "Interpolation3D1D.hxx"
#include "Interpolation.txx"
#include "MeshElement.txx"
#include "PointLocator3DIntersectorP0P0.txx"
#include "PointLocator3DIntersectorP0P1.txx"
#include "PointLocator3DIntersectorP1P0.txx"
#include "PointLocator3DIntersectorP1P1.txx"
#include "Log.hxx"

#include "BBTree.txx"

#include <memory>

namespace INTERP_KERNEL
{
  /**
   *  Very similar to Interpolation3D::interpolateMeshes, except for the bounding boxes that can be
   *  adjusted in a similar fashion as in InterpolationPlanar::performAdjustmentOfBB()
   **/
  template<class MyMeshType, class MatrixType>
  typename MyMeshType::MyConnType Interpolation3D1D::interpolateMeshes(const MyMeshType& srcMesh, const MyMeshType& targetMesh, MatrixType& result, const std::string& method)
  {
    if(InterpolationOptions::getIntersectionType() != PointLocator)
      INTERP_KERNEL::Exception("Invalid 3D/1D-0D intersection type specified : must be PointLocator.");

    typedef typename MyMeshType::MyConnType ConnType;
    // create MeshElement objects corresponding to each element of the two meshes
    const ConnType numSrcElems = srcMesh.getNumberOfElements();
    const ConnType numTargetElems = targetMesh.getNumberOfElements();

    LOG(2, "Source mesh has " << numSrcElems << " elements and target mesh has " << numTargetElems << " elements ");

    std::vector< std::unique_ptr< MeshElement<ConnType> > > srcElems(numSrcElems);
    std::vector< std::unique_ptr< MeshElement<ConnType> > > targetElems(numTargetElems);

    std::map<MeshElement<ConnType>*, ConnType> indices;

    for(ConnType i = 0 ; i < numSrcElems ; ++i)
      srcElems[i].reset( new MeshElement<ConnType>(i, srcMesh) );

    for(ConnType i = 0 ; i < numTargetElems ; ++i)
      targetElems[i].reset( new MeshElement<ConnType>(i, targetMesh) );

    std::unique_ptr< Intersector3D<MyMeshType,MatrixType> > intersector;
    std::string methC = InterpolationOptions::filterInterpolationMethod(method);
    if(methC=="P0P0")
      { intersector.reset( new PointLocator3DIntersectorP0P0<MyMeshType,MatrixType>(targetMesh, srcMesh, getPrecision()) );
      }
    else if(methC=="P0P1")
      {  intersector.reset( new PointLocator3DIntersectorP0P1<MyMeshType,MatrixType>(targetMesh, srcMesh, getPrecision()) );
      }
    else if(methC=="P1P0")
      {  intersector.reset( new PointLocator3DIntersectorP1P0<MyMeshType,MatrixType>(targetMesh, srcMesh, getPrecision()) );
      }
    else if(methC=="P1P1")
      {  intersector.reset( new PointLocator3DIntersectorP1P1<MyMeshType,MatrixType>(targetMesh, srcMesh, getPrecision()) );
      }
    else
      throw Exception("Invalid method chosen must be in \"P0P0\", \"P0P1\", \"P1P0\" or \"P1P1\".");
    // create empty maps for all source elements
    result.resize(intersector->getNumberOfRowsOfResMatrix());

    // create BBTree structure
    // - get bounding boxes
    std::vector<double> bboxes(6*numSrcElems);
    std::unique_ptr<ConnType[]> srcElemIdx{ new ConnType[numSrcElems] };
    for(ConnType i = 0; i < numSrcElems ; ++i)
      {
        // get source bboxes in right order
        const BoundingBox* box = srcElems[i]->getBoundingBox();
        bboxes[6*i+0] = box->getCoordinate(BoundingBox::XMIN);
        bboxes[6*i+1] = box->getCoordinate(BoundingBox::XMAX);
        bboxes[6*i+2] = box->getCoordinate(BoundingBox::YMIN);
        bboxes[6*i+3] = box->getCoordinate(BoundingBox::YMAX);
        bboxes[6*i+4] = box->getCoordinate(BoundingBox::ZMIN);
        bboxes[6*i+5] = box->getCoordinate(BoundingBox::ZMAX);

        srcElemIdx[i] = srcElems[i]->getIndex();
      }

    adjustBoundingBoxes(bboxes);
    const double *bboxPtr = nullptr;
    if(numSrcElems>0)
      bboxPtr=bboxes.data();
    BBTree<3,ConnType> tree(bboxPtr, srcElemIdx.get(), 0, numSrcElems);

    // for each target element, get source elements with which to calculate intersection
    // - calculate intersection by calling intersectCells
    for(ConnType i = 0; i < numTargetElems; ++i)
      {
        const BoundingBox* box = targetElems[i]->getBoundingBox();
        const ConnType targetIdx = targetElems[i]->getIndex();

        // get target bbox in right order
        double targetBox[6];
        targetBox[0] = box->getCoordinate(BoundingBox::XMIN);
        targetBox[1] = box->getCoordinate(BoundingBox::XMAX);
        targetBox[2] = box->getCoordinate(BoundingBox::YMIN);
        targetBox[3] = box->getCoordinate(BoundingBox::YMAX);
        targetBox[4] = box->getCoordinate(BoundingBox::ZMIN);
        targetBox[5] = box->getCoordinate(BoundingBox::ZMAX);

        std::vector<ConnType> intersectElems;

        tree.getIntersectingElems(targetBox, intersectElems);

        if ( !intersectElems.empty() )
          intersector->intersectCells(targetIdx,intersectElems,result);
      }
    return intersector->getNumberOfColsOfResMatrix();
  }
}
