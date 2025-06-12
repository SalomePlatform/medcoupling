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
#include "InterpolationHelper.txx"

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
 * The length of the vector is equal to the number of target elements - for each target element there is a map,
 regardless
 * of whether the element intersects any source elements or not. But in the maps there are only entries for those source
 elements
 * which have a non-zero intersection volume with the target element. The vector has indices running from
 * 0 to (nb target elements - 1), meaning that the map for target element i is stored at index i - 1. In the maps,
 however,
 * the indexing is more natural : the intersection volume of the target element i with source element j is found at
 matrix[i-1][j].
 *

 * @param srcMesh     3DSurf source mesh (meshDim=2,spaceDim=3)
 * @param targetMesh  3D target mesh, containing only tetraedra
 * @param matrix      matrix in which the result is stored
 *
 */
template <class MyMeshType, class MyMatrixType>
typename MyMeshType::MyConnType
Interpolation2D3D::interpolateMeshes(
    const MyMeshType &srcMesh, const MyMeshType &targetMesh, MyMatrixType &matrix, const std::string &method
)
{
    typedef typename MyMeshType::MyConnType ConnType;
    // create MeshElement objects corresponding to each element of the two meshes
    const ConnType numTargetElems = targetMesh.getNumberOfElements();

    LOG(2, "Target mesh has " << numTargetElems << " elements ");

    DuplicateFacesType intersectFaces;

    std::unique_ptr<Intersector3D<MyMeshType, MyMatrixType> > intersector;
    std::string methC = InterpolationOptions::filterInterpolationMethod(method);
    const double dimCaracteristic =
        CalculateCharacteristicSizeOfMeshes(srcMesh, targetMesh, InterpolationOptions::getPrintLevel());
    if (methC == "P0P0")
    {
        switch (InterpolationOptions::getIntersectionType())
        {
            case Triangulation:
                intersector.reset(new Polyhedron3D2DIntersectorP0P0<MyMeshType, MyMatrixType>(
                    targetMesh, srcMesh, dimCaracteristic, getPrecision(), intersectFaces, getSplittingPolicy()
                ));
                break;
            default:
                throw INTERP_KERNEL::Exception(
                    "Invalid 2D to 3D intersection type for P0P0 interp specified : must be Triangulation."
                );
        }
    }
    else
        throw Exception("Invalid method chosen must be in \"P0P0\".");
    // create empty maps for all source elements
    matrix.resize(intersector->getNumberOfRowsOfResMatrix());

    // create BBTree structure
    // [ABN] Adjust 2D bounding box (those might be flat in the cases where the 2D surf are perfectly aligned with the
    // axis)
    BBTreeStandAlone<3, ConnType> tree(BuildBBTreeWithAdjustment(
        srcMesh,
        [this, &intersector](double *bbox, typename MyMeshType::MyConnType sz)
        { this->performAdjustmentOfBB(intersector.get(), bbox, sz); }
    ));

    // for each target element, get source elements with which to calculate intersection
    // - calculate intersection by calling intersectCells
    for (ConnType i = 0; i < numTargetElems; ++i)
    {
        MeshElement<ConnType> trgMeshElem(i, targetMesh);

        const BoundingBox *box = trgMeshElem.getBoundingBox();

        // get target bbox in right order
        double targetBox[6];
        box->fillInXMinXmaxYminYmaxZminZmaxFormat(targetBox);

        std::vector<ConnType> intersectElems;

        tree.getIntersectingElems(targetBox, intersectElems);

        if (!intersectElems.empty())
            intersector->intersectCells(i, intersectElems, matrix);
    }

    DuplicateFacesType::iterator iter;
    for (iter = intersectFaces.begin(); iter != intersectFaces.end(); ++iter)
    {
        if (iter->second.size() > 1)
        {
            _duplicate_faces.insert(std::make_pair(iter->first, iter->second));
        }
    }

    // free allocated memory
    ConnType ret = intersector->getNumberOfColsOfResMatrix();

    return ret;
}
}  // namespace INTERP_KERNEL

#endif
