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
#ifndef __INTERPOLATION2D3D_TXX__
#define __INTERPOLATION2D3D_TXX__

#include "Interpolation2D3D.hxx"
#include "Interpolation.txx"
#include "MeshElement.txx"
#include "TransformedTriangle.hxx"
#include "Polyhedron3D2DIntersectorP0P0.txx"
#include "PointLocator3DIntersectorP0P0.txx"
#include "PolyhedronIntersectorP0P1.txx"
#include "PointLocator3DIntersectorP0P1.txx"
#include "PolyhedronIntersectorP1P0.txx"
#include "PolyhedronIntersectorP1P0Bary.txx"
#include "PointLocator3DIntersectorP1P0.txx"
#include "PolyhedronIntersectorP1P1.txx"
#include "PointLocator3DIntersectorP1P1.txx"
#include "Log.hxx"

#include "BBTree.txx"

namespace INTERP_KERNEL
{
  /**
   * Calculates the matrix of volumes of intersection between the elements of srcMesh and the elements of targetMesh.
   * The calculation is done in two steps. First a filtering process reduces the number of pairs of elements for which the
   * calculation must be carried out by eliminating pairs that do not intersect based on their bounding boxes. Then, the 
   * volume of intersection is calculated by an object of type Intersector3D for the remaining pairs, and entered into the
   * intersection matrix. 
   * 
   * The matrix is partially sparse : it is a vector of maps of integer - double pairs. 
   * It can also be an INTERP_KERNEL::Matrix object.
   * The length of the vector is equal to the number of target elements - for each target element there is a map, regardless
   * of whether the element intersects any source elements or not. But in the maps there are only entries for those source elements
   * which have a non-zero intersection volume with the target element. The vector has indices running from 
   * 0 to (nb target elements - 1), meaning that the map for target element i is stored at index i - 1. In the maps, however,
   * the indexing is more natural : the intersection volume of the target element i with source element j is found at matrix[i-1][j].
   * 

   * @param srcMesh     3DSurf source mesh (meshDim=2,spaceDim=3)
   * @param targetMesh  3D target mesh, containing only tetraedra
   * @param matrix      matrix in which the result is stored
   *
   */
  template<class MyMeshType, class MyMatrixType>
  int Interpolation2D3D::interpolateMeshes(const MyMeshType& srcMesh,
                                           const MyMeshType& targetMesh,
                                           MyMatrixType& matrix,
                                           const std::string& method)
  {
    typedef typename MyMeshType::MyConnType ConnType;
    // create MeshElement objects corresponding to each element of the two meshes
    const unsigned long numSrcElems = srcMesh.getNumberOfElements();
    const unsigned long numTargetElems = targetMesh.getNumberOfElements();

    LOG(2, "Source mesh has " << numSrcElems << " elements and target mesh has " << numTargetElems << " elements ");

    std::vector<MeshElement<ConnType>*> srcElems(numSrcElems);
    std::vector<MeshElement<ConnType>*> targetElems(numTargetElems);

    std::map<MeshElement<ConnType>*, int> indices;
    DuplicateFacesType intersectFaces;

    for(unsigned long i = 0 ; i < numSrcElems ; ++i)
      srcElems[i] = new MeshElement<ConnType>(i, srcMesh);       

    for(unsigned long i = 0 ; i < numTargetElems ; ++i)
      targetElems[i] = new MeshElement<ConnType>(i, targetMesh);

    Intersector3D<MyMeshType,MyMatrixType>* intersector=0;
    std::string methC = InterpolationOptions::filterInterpolationMethod(method);
    const double dimCaracteristic = CalculateCharacteristicSizeOfMeshes(srcMesh, targetMesh, InterpolationOptions::getPrintLevel());
    if(methC=="P0P0")
      {
        switch(InterpolationOptions::getIntersectionType())
          {
          case Triangulation:
            intersector=new Polyhedron3D2DIntersectorP0P0<MyMeshType,MyMatrixType>(targetMesh,
                                                                                   srcMesh,
                                                                                   dimCaracteristic,
                                                                                   getPrecision(),
                                                                                   intersectFaces,
                                                                                   getSplittingPolicy());
            break;
          case PointLocator:// switch target and source
            intersector=new PointLocator3DIntersectorP0P0<MyMeshType,MyMatrixType>(srcMesh,targetMesh,getPrecision());
            break;
          default:
            throw INTERP_KERNEL::Exception("Invalid 3D to 2D intersection type for P0P0 interp specified : must be Triangulation or PointLocator.");
          }
      }
    else
      throw Exception("Invalid method chosen must be in \"P0P0\".");
    // create empty maps for all source elements
    matrix.resize(intersector->getNumberOfRowsOfResMatrix());

    // create BBTree structure
    // - get bounding boxes
    double* bboxes = new double[6 * numSrcElems];
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

        // source indices have to begin with zero for BBox, I think
        srcElemIdx[i] = srcElems[i]->getIndex();
      }

    BBTree<3,ConnType> tree(bboxes, srcElemIdx, 0, numSrcElems, 0.);

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
            intersector->intersectCells(targetIdx, intersectElems, matrix);

      }

    delete [] bboxes;
    delete [] srcElemIdx;

    DuplicateFacesType::iterator iter;
    for (iter = intersectFaces.begin(); iter != intersectFaces.end(); ++iter)
      {
        if (iter->second.size() > 1)
          {
            _duplicate_faces.insert(std::make_pair(iter->first, iter->second));
          }
      }

    // free allocated memory
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
