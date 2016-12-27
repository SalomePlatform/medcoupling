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

#ifndef __INTERPOLATION3D1D_TXX__
#define __INTERPOLATION3D1D_TXX__

#include "Interpolation3D1D.hxx"
#include "Interpolation.txx"
#include "MeshElement.txx"
#include "PointLocator3DIntersectorP0P0.txx"
#include "PointLocator3DIntersectorP0P1.txx"
#include "PointLocator3DIntersectorP1P0.txx"
#include "PointLocator3DIntersectorP1P1.txx"
#include "Log.hxx"

#include "BBTree.txx"

namespace INTERP_KERNEL
{
  /**
   *  Very similar to Interpolation3D::interpolateMeshes, except for the bounding boxes that can be
   *  adjusted in a similar fashion as in InterpolationPlanar::performAdjustmentOfBB()
   **/
  template<class MyMeshType, class MatrixType>
  int Interpolation3D1D::interpolateMeshes(const MyMeshType& srcMesh, const MyMeshType& targetMesh, MatrixType& result, const std::string& method)
  {
    if(InterpolationOptions::getIntersectionType() != PointLocator)
      INTERP_KERNEL::Exception("Invalid 3D/1D intersection type specified : must be PointLocator.");

    typedef typename MyMeshType::MyConnType ConnType;
    // create MeshElement objects corresponding to each element of the two meshes
    const unsigned long numSrcElems = srcMesh.getNumberOfElements();
    const unsigned long numTargetElems = targetMesh.getNumberOfElements();

    LOG(2, "Source mesh has " << numSrcElems << " elements and target mesh has " << numTargetElems << " elements ");

    std::vector<MeshElement<ConnType>*> srcElems(numSrcElems);
    std::vector<MeshElement<ConnType>*> targetElems(numTargetElems);

    std::map<MeshElement<ConnType>*, int> indices;

    for(unsigned long i = 0 ; i < numSrcElems ; ++i)
      srcElems[i] = new MeshElement<ConnType>(i, srcMesh);       

    for(unsigned long i = 0 ; i < numTargetElems ; ++i)
      targetElems[i] = new MeshElement<ConnType>(i, targetMesh);

    Intersector3D<MyMeshType,MatrixType>* intersector=0;
    std::string methC = InterpolationOptions::filterInterpolationMethod(method);
    if(methC=="P0P0")
      { intersector=new PointLocator3DIntersectorP0P0<MyMeshType,MatrixType>(targetMesh, srcMesh, getPrecision());
      }
    else if(methC=="P0P1")
      {  intersector=new PointLocator3DIntersectorP0P1<MyMeshType,MatrixType>(targetMesh, srcMesh, getPrecision());
      }
    else if(methC=="P1P0")
      {  intersector=new PointLocator3DIntersectorP1P0<MyMeshType,MatrixType>(targetMesh, srcMesh, getPrecision());
      }
    else if(methC=="P1P1")
      {  intersector=new PointLocator3DIntersectorP1P1<MyMeshType,MatrixType>(targetMesh, srcMesh, getPrecision());
      }
    else
      throw Exception("Invalid method chosen must be in \"P0P0\", \"P0P1\", \"P1P0\" or \"P1P1\".");
    // create empty maps for all source elements
    result.resize(intersector->getNumberOfRowsOfResMatrix());

    // create BBTree structure
    // - get bounding boxes
    std::vector<double> bboxes(6*numSrcElems);
    int* srcElemIdx = new int[numSrcElems];
    for(unsigned long i = 0; i < numSrcElems ; ++i)
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
    const double *bboxPtr=0;
    if(numSrcElems>0)
      bboxPtr=&bboxes[0];
    BBTree<3,ConnType> tree(bboxPtr, srcElemIdx, 0, numSrcElems);

    // for each target element, get source elements with which to calculate intersection
    // - calculate intersection by calling intersectCells
    for(unsigned long i = 0; i < numTargetElems; ++i)
      {
        const BoundingBox* box = targetElems[i]->getBoundingBox();
        const int targetIdx = targetElems[i]->getIndex();

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

    // free allocated memory
    delete [] srcElemIdx;

    int ret=intersector->getNumberOfColsOfResMatrix();

    delete intersector;

    for(unsigned long i = 0 ; i < numSrcElems ; ++i)
      {
        delete srcElems[i];
      }
    for(unsigned long i = 0 ; i < numTargetElems ; ++i)
      {
        delete targetElems[i];
      }
    return ret;

  }
}

#endif
