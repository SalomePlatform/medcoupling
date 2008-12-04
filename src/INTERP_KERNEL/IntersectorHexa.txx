//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D, OPEN CASCADE
//
//  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#ifndef __INTERSECTORHEXA_TXX__
#define __INTERSECTORHEXA_TXX__

#include "IntersectorHexa.hxx"

#include "MeshUtils.hxx"

#include "IntersectorTetra.hxx"

namespace INTERP_KERNEL
{

  /**
   * Constructor creating object from target cell global number 
   * The constructor first calculates the necessary nodes, 
   * (depending on the splitting policy) and then splits the hexahedron into 
   * tetrahedra, placing these in the internal vector _tetra.
   * 
   * @param srcMesh     mesh containing the source elements
   * @param targetMesh  mesh containing the target elements
   * @param targetCell  global number of the target cell
   * @param policy      splitting policy to be used
   */
  template<class MyMeshType>
  IntersectorHexa<MyMeshType>::IntersectorHexa(const MyMeshType& srcMesh, const MyMeshType& targetMesh,
                                               typename MyMeshType::MyConnType targetCell, SplittingPolicy policy)
  {
    const int numTetra = static_cast<int>(policy);

    // pre-calculate nodes
    calculateSubNodes(targetMesh, targetCell, policy);

    _tetra.reserve(numTetra);
    _nodes.reserve(30); // we never have more than this

    switch(policy)
      {
      case PLANAR_FACE_5:
        {
          const int subZone[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
          fiveSplit(srcMesh, subZone);
        }
        break;

      case PLANAR_FACE_6:
        {
          const int subZone[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
          sixSplit(srcMesh, subZone);
        }
        break;

      case GENERAL_24:
        {
          calculateGeneral24Tetra(srcMesh);
        }
        break;

      case GENERAL_48:
        {
          calculateGeneral48Tetra(srcMesh);
        }
        break;
      default:
        assert(false);
      }
  }
    
  /**
   * Destructor.
   * Liberates the IntersectorTetra objects and potential sub-node points that have been allocated.
   *
   */
  template<class MyMeshType>
  IntersectorHexa<MyMeshType>::~IntersectorHexa()
  {
    for(typename std::vector< IntersectorTetra<MyMeshType>* >::iterator iter = _tetra.begin(); iter != _tetra.end(); ++iter)
      delete *iter;
    
    // free potential sub-mesh nodes that have been allocated
    std::vector<const double*>::iterator iter = _nodes.begin() + 8;
    while(iter != _nodes.end())
      {
        delete[] *iter;
        ++iter;
      }
  }

  /**
   * Splits the hexahedron into five tetrahedra.
   * This method adds five IntersectorTetra objects to the vector _tetra. 
   *
   * @param srcMesh  the source mesh
   * @param subZone  the local node numbers corresponding to the hexahedron corners - these are mapped onto {0,..,7}. Providing this allows the 
   *                 splitting to be reused on the subzones of the GENERAL_* types of splitting
   */
  template<class MyMeshType>
  void IntersectorHexa<MyMeshType>::fiveSplit(const MyMeshType& srcMesh, const int* const subZone)
  {
    // Schema according to which the splitting is performed.
    // Each line represents one tetrahedron. The numbering is as follows :
    //
    //          7 ------ 6
    //         /|       /|
    //        / |      / |
    //       3 ------ 2  |
    //       |  |     |  |
    //       |  |     |  |
    //       |  4-----|- 5
    //       | /      | /
    //       0 ------ 1


    static const int SPLIT_NODES_5[20] = 
      {
        0, 1, 5, 2,
        0, 4, 5, 7,
        0, 3, 7, 2,
        5, 6, 7, 2,
        0, 2, 5, 7
      };
    
    // create tetrahedra
    for(int i = 0; i < 5; ++i)
      {
        const double* nodes[4];
        for(int j = 0; j < 4; ++j)
          {
            nodes[j] = getCoordsOfSubNode(subZone[ SPLIT_NODES_5[4*i+j] ]);
          }
        IntersectorTetra<MyMeshType>* t = new IntersectorTetra<MyMeshType>(srcMesh, nodes);
        _tetra.push_back(t);
      }
  }

  /**
   * Splits the hexahedron into six tetrahedra.
   * This method adds six IntersectorTetra objects to the vector _tetra. 
   *
   * @param srcMesh  the source mesh
   * @param subZone  the local node numbers corresponding to the hexahedron corners - these are mapped onto {0,..,7}. Providing this allows the 
   *                 splitting to be reused on the subzones of the GENERAL_* types of splitting
   */
  template<class MyMeshType>
  void IntersectorHexa<MyMeshType>::sixSplit(const MyMeshType& srcMesh, const int* const subZone)
  {
    // Schema according to which the splitting is performed.
    // Each line represents one tetrahedron. The numbering is as follows :
    //
    //          7 ------ 6
    //         /|       /|
    //        / |      / |
    //       3 ------ 2  |
    //       |  |     |  |
    //       |  |     |  |
    //       |  4-----|- 5
    //       | /      | /
    //       0 ------ 1

    static const int SPLIT_NODES_6[24] = 
      {
        0, 1, 5, 6,
        0, 2, 1, 6,
        0, 5, 4, 6,
        0, 4, 7, 6,
        0, 3, 2, 6,
        0, 7, 3, 6
      };

    for(int i = 0; i < 6; ++i)
      {
        const double* nodes[4];
        for(int j = 0; j < 4; ++j)
          {
            nodes[j] = getCoordsOfSubNode(subZone[ SPLIT_NODES_6[4*i+j] ]);
          }
        IntersectorTetra<MyMeshType>* t = new IntersectorTetra<MyMeshType>(srcMesh, nodes);
        _tetra.push_back(t);
      }
  }

  /**
   * Splits the hexahedron into 24 tetrahedra.
   * The splitting is done by combining the barycenter of the tetrahedron, the barycenter of each face 
   * and the nodes of each edge of the face. This creates 6 faces * 4 edges / face = 24 tetrahedra.
   * The submesh nodes introduced are the barycenters of the faces and the barycenter of the cell.
   *
   * @param srcMesh  the source mesh
   * 
   */
  template<class MyMeshType>
  void IntersectorHexa<MyMeshType>::calculateGeneral24Tetra(const MyMeshType& srcMesh)
  {
    // The two nodes of the original mesh cell used in each tetrahedron.
    // The tetrahedra all have nodes (cellCenter, faceCenter, edgeNode1, edgeNode2)
    // For the correspondance of the nodes, see the GENERAL_48_SUB_NODES table in calculateSubNodes
    static const int TETRA_EDGES[48] = 
      {
        // face with center 9
        1, 2,
        2, 6,
        6, 5,
        5, 1,
        // face with center 10
        1, 2,
        2, 3,
        3, 4, 
        4, 1,
        // face with center 11
        1, 5,
        5, 8,
        8, 4,
        4, 1,
        // face with center 12
        2, 6,
        6, 7,
        7, 3,
        3, 2,
        // face with center 13
        6, 7,
        7, 8,
        8, 5,
        5, 6,
        // face with center 14
        3, 7,
        7, 8, 
        8, 4,
        4, 3
      };
    
    // nodes to use for tetrahedron
    const double* nodes[4];
    
    // get the cell center
    nodes[0] = getCoordsOfSubNode(15);
    
    for(int faceCenterNode = 9; faceCenterNode <= 14; ++faceCenterNode)
      {
        // get the face center
        nodes[1] = getCoordsOfSubNode(faceCenterNode);
        for(int j = 0; j < 4; ++j)
          {
            const int row = 4*(faceCenterNode - 9) + j;
            nodes[2] = getCoordsOfSubNode(TETRA_EDGES[2*row]);
            nodes[3] = getCoordsOfSubNode(TETRA_EDGES[2*row + 1]);
           
            IntersectorTetra<MyMeshType>* t = new IntersectorTetra<MyMeshType>(srcMesh, nodes);
            _tetra.push_back(t);
          }
      }
  }


  /**
   * Splits the hexahedron into 48 tetrahedra.
   * The splitting is done by introducing the midpoints of all the edges 
   * and the barycenter of the element as submesh nodes. The 8 hexahedral subzones thus defined
   * are then split into 6 tetrahedra each, as in Grandy, p. 449. The division of the subzones 
   * is done by calling sixSplit().
   *
   * @param srcMesh  the source mesh
   * 
   */
  template<class MyMeshType>
  void IntersectorHexa<MyMeshType>::calculateGeneral48Tetra(const MyMeshType& srcMesh)
  {
    // Define 8 hexahedral subzones as in Grandy, p449
    // the values correspond to the nodes that correspond to nodes 1,2,3,4,5,6,7,8 in the subcell
    // For the correspondance of the nodes, see the GENERAL_48_SUB_NODES table in calculateSubNodes
    static const int subZones[64] = 
      {
        1, 9, 22, 13, 10, 21, 27, 23,
        9, 2, 14, 22, 21, 11, 24, 27,
        13, 22, 17, 4, 23, 27, 26, 18,
        22, 14, 3, 17, 27, 24, 19, 26,
        10, 21, 27, 23, 5, 12, 25, 15,
        21, 11, 24, 27, 12, 6, 16, 25,
        23, 27, 26, 18, 15, 25, 20, 8,
        27, 24, 19, 26, 25, 16, 7, 20
      };
    
    for(int i = 0; i < 8; ++i)
      {
        sixSplit(srcMesh, &subZones[8*i]);
      }
  }
  
  /**
   * Precalculates all the nodes.
   * Retrieves the mesh nodes and allocates the necessary sub-mesh 
   * nodes according to the splitting policy used.
   * This method is meant to be called once by the constructor.
   *
   * @param targetMesh  the target mesh
   * @param targetCell  the global number of the cell that the object represents
   * @param policy      the splitting policy of the object
   *
   */
  template<class MyMeshType>
  void IntersectorHexa<MyMeshType>::calculateSubNodes(const MyMeshType& targetMesh, typename MyMeshType::MyConnType targetCell, SplittingPolicy policy)
  {
    // retrieve real mesh nodes
    for(int node = 0; node < 8 ; ++node)
      {
        // calculate only normal nodes
        _nodes.push_back(getCoordsOfNode(node, targetCell, targetMesh));
      }

    // create sub-mesh nodes if needed
    switch(policy)
      {
      case GENERAL_24:
        {
          // Each sub-node is the barycenter of 4 other nodes.
          // For the faces, these are on the orignal mesh.
          // For the barycenter, the four face sub-nodes are used.
          static const int GENERAL_24_SUB_NODES[28] = 
            {
              1, 2, 5, 6, // sub-node 9 (face)
              1, 2, 3, 4, // sub-node 10 (face)
              1, 4, 5, 8, // sub-node 11 (face)
              2, 3, 6, 7, // sub-node 12 (face)
              5, 6, 7, 8, // sub-node 13 (face)
              3, 4, 7, 8, // sub-node 14 (face)
              10, 11, 12, 13 // sub-node 15 (cell)
            };
         
          for(int i = 0; i < 7; ++i)
            {
              double* barycenter = new double[3];
              calcBarycenter<4>(barycenter, &GENERAL_24_SUB_NODES[4*i]);
              _nodes.push_back(barycenter);
            }
        }
        break;
       
      case GENERAL_48:
        {
          // Each sub-node is the barycenter of two other nodes.
          // For the edges, these lie on the original mesh.
          // For the faces, these are the edge sub-nodes.
          // For the cell these are two face sub-nodes.
          static const int GENERAL_48_SUB_NODES[38] = 
            {
              1, 2,  // sub-node 9 (edge)
              1, 5,  // sub-node 10 (edge)
              2, 6,  // sub-node 11 (edge)
              5, 6,  // sub-node 12 (edge)
              1, 4,  // sub-node 13 (edge)
              2, 3,  // sub-node 14 (edge)
              5, 8,  // sub-node 15 (edge)
              6, 7,  // sub-node 16 (edge)
              3, 4,  // sub-node 17 (edge)
              4, 8,  // sub-node 18 (edge)
              3, 7,  // sub-node 19 (edge)
              7, 8,  // sub-node 20 (edge)
              9, 12, // sub-node 21 (face)
              13, 14,// sub-node 22 (face)
              10, 18,// sub-node 23 (face)
              11, 19,// sub-node 24 (face)
              15, 16,// sub-node 25 (face)
              17, 20,// sub-node 26 (face)
              21, 26 // sub-node 27 (cell)
            };

          for(int i = 0; i < 19; ++i)
            {
              double* barycenter = new double[3];
              calcBarycenter<2>(barycenter, &GENERAL_48_SUB_NODES[2*i]);
              _nodes.push_back(barycenter);
            }
        }
        break;
       
      default:
        break;
      }
  }       

  /**
   * Calculates the volume of intersection of an element in the source mesh and the target element
   * represented by the object.
   * The calculation is performed by calling the corresponding method for
   * each IntersectorTetra object created by the splitting.
   * 
   * @param srcCell   global number of the source element (1 <= srcCell < # source cells)
   *
   */
  template<class MyMeshType>
  double IntersectorHexa<MyMeshType>::intersectSourceCell(typename MyMeshType::MyConnType srcCell)
  {
    double volume = 0.0;
    for(typename std::vector<IntersectorTetra<MyMeshType>*>::iterator iter = _tetra.begin(); iter != _tetra.end(); ++iter)
      volume += (*iter)->intersectSourceCell(srcCell);
    return volume;
  }
}

#endif
