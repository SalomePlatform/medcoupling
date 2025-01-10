// Copyright (C) 2007-2025  CEA, EDF
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

#ifndef __INTERPOLATION3D_TXX__
#define __INTERPOLATION3D_TXX__

#include "Interpolation3D.hxx"
#include "Interpolation.txx"
#include "MeshElement.txx"
#include "TransformedTriangle.hxx"
#include "PolyhedronIntersectorP0P0.txx"
#include "PointLocator3DIntersectorP0P0.txx"
#include "PolyhedronIntersectorP0P1.txx"
#include "PointLocator3DIntersectorP0P1.txx"
#include "PolyhedronIntersectorP1P0.txx"
#include "PolyhedronIntersectorP1P0Bary.txx"
#include "PointLocator3DIntersectorP1P0.txx"
#include "PolyhedronIntersectorP1P1.txx"
#include "PointLocator3DIntersectorP1P1.txx"
#include "Barycentric3DIntersectorP1P1.txx"
#include "MappedBarycentric3DIntersectorP1P1.txx"
#include "Log.hxx"
// If defined, use recursion to traverse the binary search tree, else use the BBTree class
//#define USE_RECURSIVE_BBOX_FILTER

#ifdef USE_RECURSIVE_BBOX_FILTER
#include "MeshRegion.txx"
#include "RegionNode.hxx"
#include <stack>

#else // use BBTree class

#include "InterpolationHelper.txx"

#endif

#include <memory>

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

   * @param srcMesh     3-dimensional source mesh
   * @param targetMesh  3-dimesional target mesh, containing only tetraedra
   * @param result      matrix in which the result is stored 
   *
   */
  template<class MyMeshType, class MatrixType>
  typename MyMeshType::MyConnType Interpolation3D::interpolateMeshes(const MyMeshType& srcMesh, const MyMeshType& targetMesh, MatrixType& result, const std::string& method)
  {
    using ConnType = typename MyMeshType::MyConnType;
    // create MeshElement objects corresponding to each element of the two meshes
    const ConnType numTargetElems = targetMesh.getNumberOfElements();

    LOG(2, "Target mesh has " << numTargetElems << " elements ");

    std::unique_ptr<Intersector3D<MyMeshType,MatrixType>> intersector;
    std::string methC = InterpolationOptions::filterInterpolationMethod(method);
    if(methC=="P0P0")
      {
        switch(InterpolationOptions::getIntersectionType())
          {
          case Triangulation:
            intersector.reset( new PolyhedronIntersectorP0P0<MyMeshType,MatrixType>(targetMesh, srcMesh, getSplittingPolicy()) );
            break;
          case PointLocator:
            intersector.reset( new PointLocator3DIntersectorP0P0<MyMeshType,MatrixType>(targetMesh, srcMesh, getPrecision()) );
            break;
          default:
            throw INTERP_KERNEL::Exception("Invalid 3D intersection type for P0P0 interp specified : must be Triangle or PointLocator.");
          }
      }
    else if(methC=="P0P1")
      {
        switch(InterpolationOptions::getIntersectionType())
          {
          case Triangulation:
            intersector.reset( new PolyhedronIntersectorP0P1<MyMeshType,MatrixType>(targetMesh, srcMesh, getSplittingPolicy()) );
            break;
          case PointLocator:
            intersector.reset( new PointLocator3DIntersectorP0P1<MyMeshType,MatrixType>(targetMesh, srcMesh, getPrecision()) );
            break;
          default:
            throw INTERP_KERNEL::Exception("Invalid 3D intersection type for P0P1 interp specified : must be Triangle or PointLocator.");
          }
      }
    else if(methC=="P1P0")
      {
        switch(InterpolationOptions::getIntersectionType())
          {
          case Triangulation:
            intersector.reset( new PolyhedronIntersectorP1P0<MyMeshType,MatrixType>(targetMesh, srcMesh, getSplittingPolicy()) );
            break;
          case PointLocator:
            intersector.reset( new PointLocator3DIntersectorP1P0<MyMeshType,MatrixType>(targetMesh, srcMesh, getPrecision()) );
            break;
          case Barycentric:
            intersector.reset( new PolyhedronIntersectorP1P0Bary<MyMeshType,MatrixType>(targetMesh, srcMesh, getSplittingPolicy()) );
            break;
          default:
            throw INTERP_KERNEL::Exception("Invalid 3D intersection type for P1P0 interp specified : must be Triangle, PointLocator or Barycentric.");
          }
      }
    else if(methC=="P1P1")
      {
        switch(InterpolationOptions::getIntersectionType())
          {
          case Triangulation:
            intersector.reset( new PolyhedronIntersectorP1P1<MyMeshType,MatrixType>(targetMesh, srcMesh, getSplittingPolicy()) );
            break;
          case PointLocator:
            intersector.reset( new PointLocator3DIntersectorP1P1<MyMeshType,MatrixType>(targetMesh, srcMesh, getPrecision()) );
            break;
          case Barycentric:
            intersector.reset( new Barycentric3DIntersectorP1P1<MyMeshType,MatrixType>(targetMesh, srcMesh, getPrecision()) );
            break;
          case MappedBarycentric:
            intersector.reset( new MappedBarycentric3DIntersectorP1P1<MyMeshType,MatrixType>(targetMesh, srcMesh, getPrecision()) );
            break;
          default:
            throw INTERP_KERNEL::Exception("Invalid 3D intersection type for P1P1 interp specified : must be Triangle, PointLocator, Barycentric or MappedBarycentric.");
          }
      }
    else
      throw Exception("Invalid method chosen must be in \"P0P0\", \"P0P1\", \"P1P0\" or \"P1P1\".");
    // create empty maps for all source elements
    result.resize(intersector->getNumberOfRowsOfResMatrix());

#ifdef USE_RECURSIVE_BBOX_FILTER

    /*
     * Performs a depth-first search over srcMesh, using bounding boxes to recursively eliminate the elements of targetMesh
     * which cannot intersect smaller and smaller regions of srcMesh. At each level, each region is divided in two, forming
     * a binary search tree with leaves consisting of only one element of the source mesh together with the elements of the
     * target mesh that can intersect it. The recursion is implemented with a stack of RegionNodes, each one containing a 
     * source region and a target region. Each region has an associated bounding box and a vector of pointers to the elements 
     * that belong to it. Each MeshElement contains a bounding box and the global number of the corresponding element in the mesh.
     */

    // create initial RegionNode and fill up its source region with all the source mesh elements and
    // its target region with all the target mesh elements whose bounding box
    // intersects that of the source region

    RegionNode<ConnType>* firstNode = new RegionNode<ConnType>();

    MeshRegion<ConnType>& srcRegion = firstNode->getSrcRegion();

    for(ConnType i = 0 ; i < numSrcElems ; ++i)
      {
        srcRegion.addElement(srcElems[i], srcMesh);
      }

    MeshRegion<ConnType>& targetRegion = firstNode->getTargetRegion();

    for(ConnType i = 0 ; i < numTargetElems ; ++i)
      {
        if(!srcRegion.isDisjointWithElementBoundingBox( *(targetElems[i]) ))
          {
            targetRegion.addElement(targetElems[i], targetMesh);
          }
      }

    // Using a stack, descend recursively, creating at each step two new RegionNodes having as source region the left and
    // right part of the source region of the current node (created using MeshRegion::split()) and as target region all the 
    // elements of the target mesh whose bounding box intersects the corresponding part
    // Continue until the source region contains only one element, at which point the intersection volumes are
    // calculated with all the remaining target mesh elements and stored in the matrix if they are non-zero.

    std::stack< RegionNode<ConnType>* > nodes;
    nodes.push(firstNode);

    while(!nodes.empty())
      {
        RegionNode<ConnType>* currNode = nodes.top();
        nodes.pop();
        LOG(4, "Popping node ");

        if(currNode->getTargetRegion().getNumberOfElements() == 1)
          {
            // calculate volumes
            LOG(4, " - One element");

            MeshElement<ConnType>* targetElement = *(currNode->getTargetRegion().getBeginElements());
            std::vector<ConnType> intersectElems;
            for(typename std::vector< MeshElement<ConnType>* >::const_iterator iter = currNode->getSrcRegion().getBeginElements();iter != currNode->getSrcRegion().getEndElements();++iter)
              intersectElems.push_back((*iter)->getIndex());
            intersector->intersectCells(targetElement->getIndex(),intersectElems,result);
          }
        else // recursion 
          {

            LOG(4, " - Recursion");

            RegionNode<ConnType>* leftNode = new RegionNode<ConnType>();
            RegionNode<ConnType>* rightNode = new RegionNode<ConnType>();

            // split current source region
            //} decide on axis
            static BoundingBox::BoxCoord axis = BoundingBox::XMAX;

            currNode->getTargetRegion().split(leftNode->getTargetRegion(), rightNode->getTargetRegion(), axis, targetMesh);

            LOG(5, "After split, left target region has " << leftNode->getTargetRegion().getNumberOfElements()
                << " elements and right target region has " << rightNode->getTargetRegion().getNumberOfElements() 
                << " elements");

            // ugly hack to avoid problem with enum which does not start at 0
            // I guess I ought to implement ++ for it instead ...
            // Anyway, it basically chooses the next axis, cyclically
            axis = (axis != BoundingBox::ZMAX) ? static_cast<BoundingBox::BoxCoord>(axis + 1) : BoundingBox::XMAX;

            // add source elements of current node that overlap the target regions of the new nodes
            LOG(5, " -- Adding source elements");
            ConnType numLeftElements = 0;
            ConnType numRightElements = 0;
            for(typename std::vector<MeshElement<ConnType>*>::const_iterator iter = currNode->getSrcRegion().getBeginElements() ; 
                iter != currNode->getSrcRegion().getEndElements() ; ++iter)
              {
                LOG(6, " --- New target node");

                if(!leftNode->getTargetRegion().isDisjointWithElementBoundingBox(**iter))
                  {
                    leftNode->getSrcRegion().addElement(*iter, srcMesh);
                    ++numLeftElements;
                  }

                if(!rightNode->getTargetRegion().isDisjointWithElementBoundingBox(**iter))
                  {
                    rightNode->getSrcRegion().addElement(*iter, srcMesh);
                    ++numRightElements;
                  }

              }

            LOG(5, "Left src region has " << numLeftElements << " elements and right src region has " 
                << numRightElements << " elements");

            // push new nodes on stack
            if(numLeftElements != 0)
              {
                nodes.push(leftNode);
              }
            else
              {
                delete leftNode;
              }

            if(numRightElements != 0)
              {
                nodes.push(rightNode);
              }
            else
              {
                delete rightNode;
              }
          }

        // all nodes are deleted here
        delete currNode;

        LOG(4, "Next iteration. Nodes left : " << nodes.size());
      }

#else // Use BBTree
    // create BBTree structure
    BBTreeStandAlone<3,ConnType> tree( BuildBBTree(srcMesh) );

    // for each target element, get source elements with which to calculate intersection
    // - calculate intersection by calling intersectCells
    for(ConnType i = 0; i < numTargetElems; ++i)
      {
        MeshElement<ConnType> trgMeshElem(i, targetMesh);

        const BoundingBox *box = trgMeshElem.getBoundingBox();

        // get target bbox in right order
        double targetBox[6];
        box->fillInXMinXmaxYminYmaxZminZmaxFormat(targetBox);

        std::vector<ConnType> intersectElems;

        tree.getIntersectingElems(targetBox, intersectElems);

        if ( !intersectElems.empty() )
          intersector->intersectCells(i,intersectElems,result);
      }

#endif
    return intersector->getNumberOfColsOfResMatrix();
  }
}

#endif
