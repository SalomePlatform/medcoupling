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
#ifndef __POLYHEDRON3D2DINTERSECTORP0P0_TXX__
#define __POLYHEDRON3D2DINTERSECTORP0P0_TXX__

#include "Polyhedron3D2DIntersectorP0P0.hxx"
#include "Intersector3DP0P0.txx"
#include "MeshUtils.hxx"

#include "SplitterTetra.txx"

namespace INTERP_KERNEL
{

  /**
   * Constructor creating object from target cell global number 
   * The constructor first calculates the necessary nodes, 
   * (depending on the splitting policy) and then splits the hexahedron into 
   * tetrahedra, placing these in the internal vector _tetra.
   * 
   * @param targetMesh  mesh containing the target elements
   * @param srcMesh     mesh containing the source elements
   * @param policy      splitting policy to be used
   */
  template<class MyMeshType, class MyMatrixType>
  Polyhedron3D2DIntersectorP0P0<MyMeshType,MyMatrixType>::Polyhedron3D2DIntersectorP0P0(const MyMeshType& targetMesh,
                                                                                        const MyMeshType& srcMesh,
                                                                                        const double dimCaracteristic,
                                                                                        const double precision,
                                                                                        DuplicateFacesType& intersectFaces,
                                                                                        SplittingPolicy policy)
    : Intersector3DP0P0<MyMeshType,MyMatrixType>(targetMesh,srcMesh),
      _split(targetMesh,srcMesh,policy),
      _dim_caracteristic(dimCaracteristic),
      _precision(precision),
      _intersect_faces(intersectFaces)
  {
  }

  /**
   * Destructor.
   * Liberates the SplitterTetra objects and potential sub-node points that have been allocated.
   *
   */
  template<class MyMeshType, class MyMatrixType>
  Polyhedron3D2DIntersectorP0P0<MyMeshType,MyMatrixType>::~Polyhedron3D2DIntersectorP0P0()
  {
    releaseArrays();
  }
    
  template<class MyMeshType, class MyMatrixType>
  void Polyhedron3D2DIntersectorP0P0<MyMeshType,MyMatrixType>::releaseArrays()
  {
    for(typename std::vector< SplitterTetra<MyMeshType>* >::iterator iter = _tetra.begin(); iter != _tetra.end(); ++iter)
      delete *iter;
    _split.releaseArrays();
    _tetra.clear();
  }

  /**
   * Calculates the volume of intersection of an element in the source mesh and the target element
   * represented by the object.
   * The calculation is performed by calling the corresponding method for
   * each SplitterTetra object created by the splitting.
   * 
   * @param targetCell in C mode.
   * @param srcCells in C mode.
   */
  template<class MyMeshType, class MyMatrixType>
  void Polyhedron3D2DIntersectorP0P0<MyMeshType,MyMatrixType>::intersectCells(ConnType targetCell,
                                                                              const std::vector<ConnType>& srcCells,
                                                                              MyMatrixType& matrix)
  {
    int nbOfNodesT=Intersector3D<MyMeshType,MyMatrixType>::_target_mesh.getNumberOfNodesOfElement(OTT<ConnType,numPol>::indFC(targetCell));
    releaseArrays();
    _split.splitTargetCell(targetCell,nbOfNodesT,_tetra);

    for(typename std::vector<ConnType>::const_iterator iterCellS=srcCells.begin();iterCellS!=srcCells.end();iterCellS++)
      {
        double surface = 0.;
        std::multiset<TriangleFaceKey> listOfTetraFacesTreated;
        std::set<TriangleFaceKey> listOfTetraFacesColinear;

        // calculate the coordinates of the nodes
        typename MyMeshType::MyConnType cellSrc = *iterCellS;
        int cellSrcIdx = OTT<ConnType,numPol>::indFC(cellSrc);
        NormalizedCellType normCellType=Intersector3D<MyMeshType,MyMatrixType>::_src_mesh.getTypeOfElement(cellSrcIdx);
        const CellModel& cellModelCell=CellModel::GetCellModel(normCellType);
        const MyMeshType& src_mesh = Intersector3D<MyMeshType,MyMatrixType>::_src_mesh;
        unsigned nbOfNodes4Type=cellModelCell.isDynamic() ? src_mesh.getNumberOfNodesOfElement(cellSrcIdx) : cellModelCell.getNumberOfNodes();
        int *polyNodes=new int[nbOfNodes4Type];
        double **polyCoords = new double*[nbOfNodes4Type];
        for(int i = 0;i<(int)nbOfNodes4Type;++i)
          {
            // we could store mapping local -> global numbers too, but not sure it is worth it
            const int globalNodeNum = getGlobalNumberOfNode(i, OTT<ConnType,numPol>::indFC(*iterCellS), src_mesh);
            polyNodes[i] = globalNodeNum;
            polyCoords[i] = const_cast<double*>(src_mesh.getCoordinatesPtr()+MyMeshType::MY_SPACEDIM*globalNodeNum);
          }

        for(typename std::vector<SplitterTetra<MyMeshType>*>::iterator iter = _tetra.begin(); iter != _tetra.end(); ++iter)
            surface += (*iter)->intersectSourceFace(normCellType,
                                                    nbOfNodes4Type,
                                                    polyNodes,
                                                    polyCoords,
                                                    _dim_caracteristic,
                                                    _precision,
                                                    listOfTetraFacesTreated,
                                                    listOfTetraFacesColinear);

        if(surface!=0.)
          {
            
            matrix[targetCell].insert(std::make_pair(cellSrcIdx, surface));
            
            bool isSrcFaceColinearWithFaceOfTetraTargetCell = false;
            std::set<TriangleFaceKey>::iterator iter;
            for (iter = listOfTetraFacesColinear.begin(); iter != listOfTetraFacesColinear.end(); ++iter)
              {
                if (listOfTetraFacesTreated.count(*iter) != 1)
                  {
                    isSrcFaceColinearWithFaceOfTetraTargetCell = false;
                    break;
                  }
                else
                  {
                    isSrcFaceColinearWithFaceOfTetraTargetCell = true;
                  }
              }
            
            if (isSrcFaceColinearWithFaceOfTetraTargetCell)
              {
                DuplicateFacesType::iterator intersectFacesIter = _intersect_faces.find(cellSrcIdx);
                if (intersectFacesIter != _intersect_faces.end())
                  {
                    intersectFacesIter->second.insert(targetCell);
                  }
                else
                  {
                    std::set<int> targetCellSet;
                    targetCellSet.insert(targetCell);
                    _intersect_faces.insert(std::make_pair(cellSrcIdx, targetCellSet));
                  }
              }
          }
        delete [] polyNodes;
        delete [] polyCoords;
      }
    _split.releaseArrays();
  }
}

#endif
