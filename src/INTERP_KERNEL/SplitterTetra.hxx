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

#ifndef __SPLITTERTETRA_HXX__
#define __SPLITTERTETRA_HXX__

#include "INTERPKERNELDefines.hxx"
#include "TransformedTriangle.hxx"
#include "TetraAffineTransform.hxx"
#include "InterpolationOptions.hxx"
#include "InterpKernelException.hxx"
#include "InterpKernelHashMap.hxx"
#include "VectorUtils.hxx"

#include <functional>
#include <vector>
#include <cassert>
#include <map>
#include <set>

namespace INTERP_KERNEL
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

  static const int SPLIT_NODES_5[20] = /* WHY not all well oriented ???? */
    {
      0, 1, 5, 2,
      0, 4, 5, 7,
      0, 3, 7, 2,
      5, 6, 7, 2,
      0, 2, 5, 7
    };

  static const int SPLIT_NODES_5_WO[20] = /* WO for well oriented !!! normals of 3 first points are OUTSIDE the TETRA4 */
    {
      0, 5, 1, 2,
      0, 4, 5, 7,
      0, 3, 7, 2,
      5, 7, 6, 2,
      0, 5, 2, 7
    };

  static const int SPLIT_NODES_6[24] = /* WHY all badly oriented ???? */
    {
      0, 1, 5, 6,
      0, 2, 1, 6,
      0, 5, 4, 6,
      0, 4, 7, 6,
      0, 3, 2, 6,
      0, 7, 3, 6
    };
  
  static const int SPLIT_NODES_6_WO[24] = /* WO for well oriented !!! normals of 3 first points are OUTSIDE the TETRA4 */
    {
      0, 5, 1, 6,
      0, 1, 2, 6,
      0, 4, 5, 6,
      0, 7, 4, 6,
      0, 2, 3, 6,
      0, 3, 7, 6
    };
  
  // Each sub-node is the barycenter of 4 other nodes.
  // For the faces, these are on the orignal mesh.
  // For the barycenter, the four face sub-nodes are used.
  static const int GENERAL_24_SUB_NODES[28] = 
    {
      0,1,4,5,// sub-node 8  (face)
      0,1,2,3,// sub-node 9  (face)
      0,3,4,7,// sub-node 10 (face)
      1,2,5,6,// sub-node 11 (face)
      4,5,6,7,// sub-node 12 (face)
      2,3,6,7,// sub-node 13 (face)
      8,9,10,11// sub-node 14 (cell)
    };

  static const int GENERAL_24_SUB_NODES_WO[28] = 
    {
      0,4,5,1,// sub-node 8  (face)
      0,1,2,3,// sub-node 9  (face)
      0,3,7,4,// sub-node 10 (face)
      1,5,6,2,// sub-node 11 (face)
      4,7,6,5,// sub-node 12 (face)
      2,6,7,3,// sub-node 13 (face)
      8,9,10,11// sub-node 14 (cell)
    };
  
  static const int TETRA_EDGES_GENERAL_24[48] = 
    {
      // face with center 8
      0,1,
      1,5,
      5,4,
      4,0,
      // face with center 9
      0,1,
      1,2,
      2,3,
      3,0,
      // face with center 10
      0,4,
      4,7,
      7,3,
      3,0,
      // face with center 11
      1,5,
      5,6,
      6,2,
      2,1,
      // face with center 12
      5,6,
      6,7,
      7,4,
      4,5,
      // face with center 13
      2,6,
      6,7,
      7,3,
      3,2
    };
  
  // Each sub-node is the barycenter of two other nodes.
  // For the edges, these lie on the original mesh.
  // For the faces, these are the edge sub-nodes.
  // For the cell these are two face sub-nodes.
  static const int GENERAL_48_SUB_NODES[38] = 
    {
      0,1,   // sub-node 8 (edge)
      0,4,   // sub-node 9 (edge)
      1,5,   // sub-node 10 (edge)
      4,5,   // sub-node 11 (edge)
      0,3,   // sub-node 12 (edge)
      1,2,   // sub-node 13 (edge)
      4,7,   // sub-node 14 (edge)
      5,6,   // sub-node 15 (edge)
      2,3,   // sub-node 16 (edge)
      3,7,   // sub-node 17 (edge)
      2,6,   // sub-node 18 (edge)
      6,7,   // sub-node 19 (edge)
      8,11,  // sub-node 20 (face)
      12,13, // sub-node 21 (face)
      9,17,  // sub-node 22 (face)
      10,18, // sub-node 23 (face)
      14,15, // sub-node 24 (face)
      16,19, // sub-node 25 (face)
      20,25  // sub-node 26 (cell)
    };

  // Define 8 hexahedral subzones as in Grandy, p449
  // the values correspond to the nodes that correspond to nodes 1,2,3,4,5,6,7,8 in the subcell
  // For the correspondance of the nodes, see the GENERAL_48_SUB_NODES table in calculateSubNodes
  static const int GENERAL_48_SUBZONES[64] = 
    {
      0,8,21,12,9,20,26,22,
      8,1,13,21,20,10,23,26,
      12,21,16,3,22,26,25,17,
      21,13,2,16,26,23,18,25,
      9,20,26,22,4,11,24,14,
      20,10,23,26,11,5,15,24,
      22,26,25,17,14,24,19,7,
      26,23,18,25,24,15,6,19
    };

  static const int GENERAL_48_SUBZONES_2[64] = 
    {
      0,-1,-14,-5,-2,-13,-19,-15,
      -1,1,-6,-14,-13,-3,-16,-19,
      -5,-14,-9,3,-15,-19,-18,-10,
      -14,-6,2,-9,-19,-16,-11,-18,
      -2,-13,-19,-15,4,-4,-17,-7,
      -13,-3,-16,-19,-4,5,-8,-17,
      -15,-19,-18,-10,-7,-17,-12,7,
      -19,-16,-11,-18,-17,-8,6,-12};

  void SplitHexa8IntoTetras(SplittingPolicy policy, const int *nodalConnBg, const int *nodalConnEnd, const double *coords,
                            std::vector<int>& tetrasNodalConn, std::vector<double>& addCoords);
  
  INTERPKERNEL_EXPORT void SplitIntoTetras(SplittingPolicy policy, NormalizedCellType gt, const int *nodalConnBg, const int *nodalConnEnd, const double *coords,
                                           std::vector<int>& tetrasNodalConn, std::vector<double>& addCoords);
  
  /**
   * \brief Class representing a triangular face, used as key in caching hash map in SplitterTetra.
   *
   */
  class TriangleFaceKey
  {
  public:
    
    /**
     * Constructor
     * Sorts the given nodes (so that the order in which they are passed does not matter) and 
     * calculates a hash value for the key.
     *
     * @param node1  global number of the first node of the face
     * @param node2  global number of the second node of the face
     * @param node3  global number of the third node of the face
     */
    TriangleFaceKey(int node1, int node2, int node3)
    {
      Sort3Ints(_nodes, node1, node2, node3);
      _hashVal = ( _nodes[0] + _nodes[1] + _nodes[2] ) % 29;
    }

    /**
     * Equality comparison operator.
     * Compares this TriangleFaceKey object to another and determines if they represent the same face.
     * 
     * @param   key  TriangleFaceKey with which to compare
     * @return  true if key has the same three nodes as this object, false if not
     */ 
    bool operator==(const TriangleFaceKey& key) const
    {
      return _nodes[0] == key._nodes[0] && _nodes[1] == key._nodes[1] && _nodes[2] == key._nodes[2];
    }

    /**
     * Less than operator.
     *
     * @param   key  TriangleFaceKey with which to compare
     * @return  true if this object has the three nodes less than the nodes of the key object, false if not
     */
    bool operator<(const TriangleFaceKey& key) const
    {
      for (int i = 0; i < 3; ++i)
        {
          if (_nodes[i] < key._nodes[i])
            {
              return true;
            }
          else if (_nodes[i] > key._nodes[i])
            {
              return false;
            }
        }
      return false;
    }

    /**
     * Returns a hash value for the object, based on its three nodes.
     * This value is not unique for each face.
     *
     * @return  a hash value for the object
     */
    int hashVal() const
    {
      return _hashVal;
    }
     
    inline static void Sort3Ints(int* sorted, int node1, int node2, int node3);

  private:
    /// global numbers of the three nodes, sorted in ascending order
    int _nodes[3];
    
    /// hash value for the object, calculated in the constructor
    int _hashVal;
  };
  
  /**
   * Method to sort three integers in ascending order
   *
   * @param sorted  int[3] array in which to store the result
   * @param x1   first integer
   * @param x2   second integer
   * @param x3   third integer
   */
  inline void TriangleFaceKey::Sort3Ints(int* sorted, int x1, int x2, int x3)
  {
    if(x1 < x2)
      {
        if(x1 < x3)
          {
            // x1 is min
            sorted[0] = x1;
            sorted[1] = x2 < x3 ? x2 : x3;
            sorted[2] = x2 < x3 ? x3 : x2;
          }
        else
          {
            // x3, x1, x2
            sorted[0] = x3;
            sorted[1] = x1;
            sorted[2] = x2;
          }
      }
    else // x2 < x1
      {
        if(x2 < x3)
          {
            // x2 is min
            sorted[0] = x2;
            sorted[1] = x1 < x3 ? x1 : x3;
            sorted[2] = x1 < x3 ? x3 : x1;
          } 
        else
          {
            // x3, x2, x1
            sorted[0] = x3;
            sorted[1] = x2;
            sorted[2] = x1;
          }
      }
  }

  /**
   * \brief Template specialization of INTERP_KERNEL::hash<T> function object for use with a 
   * with TriangleFaceKey as key class.
   * 
   */
  template<> class hash<INTERP_KERNEL::TriangleFaceKey>
  {
  public:
    /**
     * Operator() that returns the precalculated hashvalue of a TriangleFaceKey object.
     *
     * @param key  a TriangleFaceKey object
     * @return an integer hash value for key
     */
    int operator()(const INTERP_KERNEL::TriangleFaceKey& key) const
    {
      return key.hashVal();
    }
  };
}

namespace INTERP_KERNEL
{
  /** 
   * \brief Class calculating the volume of intersection between a tetrahedral target element and
   * source elements with triangular or quadratilateral faces.
   *
   */
  template<class MyMeshType>
  class SplitterTetra
  {
  public: 
    
    SplitterTetra(const MyMeshType& srcMesh, const double** tetraCorners, const typename MyMeshType::MyConnType *nodesId);

    SplitterTetra(const MyMeshType& srcMesh, const double tetraCorners[12], const int *conn = 0);

    ~SplitterTetra();

    double intersectSourceCell(typename MyMeshType::MyConnType srcCell, double* baryCentre=0);
    double intersectSourceFace(const NormalizedCellType polyType,
                               const int polyNodesNbr,
                               const int *const polyNodes,
                               const double *const *const polyCoords,
                               const double dimCaracteristic,
                               const double precision,
                               std::multiset<TriangleFaceKey>& listOfTetraFacesTreated,
                               std::set<TriangleFaceKey>& listOfTetraFacesColinear);

    double intersectTetra(const double** tetraCorners);

    typename MyMeshType::MyConnType getId(int id) { return _conn[id]; }
    
    void splitIntoDualCells(SplitterTetra<MyMeshType> **output);

    void splitMySelfForDual(double* output, int i, typename MyMeshType::MyConnType& nodeId);

    void clearVolumesCache();

  private:
    inline static void CheckIsOutside(const double* pt, bool* isOutside, const double errTol = DEFAULT_ABS_TOL);
    inline static void CheckIsStrictlyOutside(const double* pt, bool* isStrictlyOutside, const double errTol = DEFAULT_ABS_TOL);
    inline void calculateNode(typename MyMeshType::MyConnType globalNodeNum);
    inline void calculateNode2(typename MyMeshType::MyConnType globalNodeNum, const double* node);
    inline void calculateVolume(TransformedTriangle& tri, const TriangleFaceKey& key);
    inline void calculateSurface(TransformedTriangle& tri, const TriangleFaceKey& key);

    static inline bool IsFacesCoplanar(const double *const planeNormal, const double planeConstant,
                                const double *const *const coordsFace, const double precision);
    static inline double CalculateIntersectionSurfaceOfCoplanarTriangles(const double *const planeNormal,
                                                                         const double planeConstant,
                                                                         const double *const p1, const double *const p2, const double *const p3,
                                                                         const double *const p4, const double *const p5, const double *const p6,
                                                                         const double dimCaracteristic, const double precision);

    /// disallow copying
    SplitterTetra(const SplitterTetra& t);
    
    /// disallow assignment
    SplitterTetra& operator=(const SplitterTetra& t);

    // member variables
    /// affine transform associated with this target element
    TetraAffineTransform* _t;
    
    /// HashMap relating node numbers to transformed nodes, used for caching
    HashMap< int , double* > _nodes;
    
    /// HashMap relating triangular faces to calculated volume contributions, used for caching
    HashMap< TriangleFaceKey, double > _volumes;

    /// reference to the source mesh
    const MyMeshType& _src_mesh;
                
    // node id of the first node in target mesh in C mode.
    typename MyMeshType::MyConnType _conn[4];

    double _coords[12];
    
    /// Smallest volume of the intersecting elements in the transformed space that will be returned as non-zero. 
    /// Since the scale is always the same in the transformed space (the target tetrahedron is unitary), this number is independent of the scale of the meshes.
    static const double SPARSE_TRUNCATION_LIMIT;
  };

  /**
   * Function used to filter out elements by checking if they belong to one of the halfspaces
   * x <= 0, x >= 1, y <= 0, y >= 1, z <= 0, z >= 1, (indexed 0 - 7). The function updates an array of boolean variables
   * which indicates whether the points that have already been checked are all in a halfspace. For each halfspace, 
   * the corresponding array element will be true if and only if it was true when the method was called and pt lies in the halfspace.
   * 
   * @param pt        double[3] containing the coordiantes of a transformed point
   * @param isOutside bool[8] which indicate the results of earlier checks. 
   */
  template<class MyMeshType>
  inline void SplitterTetra<MyMeshType>::CheckIsOutside(const double* pt, bool* isOutside, const double errTol)
  {
    isOutside[0] = isOutside[0] && (pt[0] < errTol);
    isOutside[1] = isOutside[1] && (pt[0] > (1.0-errTol) );
    isOutside[2] = isOutside[2] && (pt[1] < errTol);
    isOutside[3] = isOutside[3] && (pt[1] > (1.0-errTol));
    isOutside[4] = isOutside[4] && (pt[2] < errTol);
    isOutside[5] = isOutside[5] && (pt[2] > (1.0-errTol));
    isOutside[6] = isOutside[6] && (1.0 - pt[0] - pt[1] - pt[2] < errTol);
    isOutside[7] = isOutside[7] && (1.0 - pt[0] - pt[1] - pt[2] > (1.0-errTol) );
  }
  
  template<class MyMeshType>
  inline void SplitterTetra<MyMeshType>::CheckIsStrictlyOutside(const double* pt, bool* isStrictlyOutside, const double errTol)
  {
    isStrictlyOutside[0] = isStrictlyOutside[0] && (pt[0] < -errTol);
    isStrictlyOutside[1] = isStrictlyOutside[1] && (pt[0] > (1.0 + errTol));
    isStrictlyOutside[2] = isStrictlyOutside[2] && (pt[1] < -errTol);
    isStrictlyOutside[3] = isStrictlyOutside[3] && (pt[1] > (1.0 + errTol));
    isStrictlyOutside[4] = isStrictlyOutside[4] && (pt[2] < -errTol);
    isStrictlyOutside[5] = isStrictlyOutside[5] && (pt[2] > (1.0 + errTol));
    isStrictlyOutside[6] = isStrictlyOutside[6] && (1.0 - pt[0] - pt[1] - pt[2] < -errTol);
    isStrictlyOutside[7] = isStrictlyOutside[7] && (1.0 - pt[0] - pt[1] - pt[2] > (1.0 + errTol));
  }

  /**
   * Calculates the transformed node with a given global node number.
   * Gets the coordinates for the node in _src_mesh with the given global number and applies TetraAffineTransform
   * _t to it. Stores the result in the cache _nodes. The non-existance of the node in _nodes should be verified before
   * calling.
   *
   * @param globalNodeNum  global node number of the node in the mesh _src_mesh
   *
   */
  template<class MyMeshType>
  inline void SplitterTetra<MyMeshType>::calculateNode(typename MyMeshType::MyConnType globalNodeNum)
  {  
    const double* node = _src_mesh.getCoordinatesPtr()+MyMeshType::MY_SPACEDIM*globalNodeNum;
    double* transformedNode = new double[MyMeshType::MY_SPACEDIM];
    assert(transformedNode != 0);
    _t->apply(transformedNode, node);
    _nodes[globalNodeNum] = transformedNode;
  }


  /**
   * Calculates the transformed node with a given global node number.
   * Applies TetraAffineTransform * _t to it.
   * Stores the result in the cache _nodes. The non-existence of the node in _nodes should be verified before * calling.
   * The only difference with the previous method calculateNode is that the coordinates of the node are passed in arguments
   * and are not recalculated in order to optimize the method.
   *
   * @param globalNodeNum  global node number of the node in the mesh _src_mesh
   *
   */
  template<class MyMeshType>
  inline void SplitterTetra<MyMeshType>::calculateNode2(typename MyMeshType::MyConnType globalNodeNum, const double* node)
  {
    double* transformedNode = new double[MyMeshType::MY_SPACEDIM];
    assert(transformedNode != 0);
    _t->apply(transformedNode, node);
    _nodes[globalNodeNum] = transformedNode;
  }

  /**
   * Calculates the volume contribution from the given TransformedTriangle and stores it with the given key in .
   * Calls TransformedTriangle::calculateIntersectionVolume to perform the calculation.
   *
   * @param tri    triangle for which to calculate the volume contribution
   * @param key    key associated with the face
   */
  template<class MyMeshType>
  inline void SplitterTetra<MyMeshType>::calculateVolume(TransformedTriangle& tri, const TriangleFaceKey& key)
  {
    const double vol = tri.calculateIntersectionVolume();
    _volumes.insert(std::make_pair(key, vol));
  }

  /**
   * Calculates the surface contribution from the given TransformedTriangle and stores it with the given key in.
   * Calls TransformedTriangle::calculateIntersectionSurface to perform the calculation.
   *
   * @param tri    triangle for which to calculate the surface contribution
   * @param key    key associated with the face
   */
  template<class MyMeshType>
  inline void SplitterTetra<MyMeshType>::calculateSurface(TransformedTriangle& tri, const TriangleFaceKey& key)
  {
    const double surf = tri.calculateIntersectionSurface(_t);
    _volumes.insert(std::make_pair(key, surf));
  }

  template<class MyMeshTypeT, class MyMeshTypeS=MyMeshTypeT>
  class SplitterTetra2
  {
  public:
    SplitterTetra2(const MyMeshTypeT& targetMesh, const MyMeshTypeS& srcMesh, SplittingPolicy policy);
    ~SplitterTetra2();
    void releaseArrays();
    void splitTargetCell2(typename MyMeshTypeT::MyConnType targetCell, typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra);
    void splitTargetCell(typename MyMeshTypeT::MyConnType targetCell, typename MyMeshTypeT::MyConnType nbOfNodesT,
                         typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra);//to suppress
    void fiveSplit(const int* const subZone, typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra);//to suppress
    void sixSplit(const int* const subZone, typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra);//to suppress
    void calculateGeneral24Tetra(typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra);//to suppress
    void calculateGeneral48Tetra(typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra);//to suppress
    void splitPyram5(typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra);//to suppress
    void splitConvex(typename MyMeshTypeT::MyConnType                     targetCell,//to suppress
                     typename std::vector< SplitterTetra<MyMeshTypeS>* >& tetra);//to suppress
    void calculateSubNodes(const MyMeshTypeT& targetMesh, typename MyMeshTypeT::MyConnType targetCell);//to suppress
    inline const double* getCoordsOfSubNode(typename MyMeshTypeT::MyConnType node);//to suppress
    inline const double* getCoordsOfSubNode2(typename MyMeshTypeT::MyConnType node, typename MyMeshTypeT::MyConnType& nodeId);//to suppress
    //template<int n>
    inline void calcBarycenter(int n, double* barycenter, const typename MyMeshTypeT::MyConnType* pts);//to suppress
  private:
    const MyMeshTypeT& _target_mesh;
    const MyMeshTypeS& _src_mesh;
    SplittingPolicy _splitting_pol;
    /// vector of pointers to double[3] containing the coordinates of the
    /// (sub) - nodes of split target cell
    std::vector<const double*> _nodes;
    std::vector<typename MyMeshTypeT::MyConnType> _node_ids;
  };

  /**
   * Calculates the barycenter of n (sub) - nodes
   *
   * @param  n  number of nodes for which to calculate barycenter
   * @param  barycenter  pointer to double[3] array in which to store the result
   * @param  pts pointer to int[n] array containing the (sub)-nodes for which to calculate the barycenter
   */
  template<class MyMeshTypeT, class MyMeshTypeS>
  //template<int n>
  inline void SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::calcBarycenter(int n, double* barycenter, const typename MyMeshTypeT::MyConnType* pts)
  {
    barycenter[0] = barycenter[1] = barycenter[2] = 0.0;
    for(int i = 0; i < n ; ++i)
      {
       const double* pt = getCoordsOfSubNode(pts[i]);
       barycenter[0] += pt[0];
       barycenter[1] += pt[1];
       barycenter[2] += pt[2];
      }
    
    barycenter[0] /= n;
    barycenter[1] /= n;
    barycenter[2] /= n;
  }

  /**
   * Accessor to the coordinates of a given (sub)-node
   *
   * @param  node  local number of the (sub)-node 0,..,7 are the elements nodes, sub-nodes are numbered from 8,..
   * @return pointer to double[3] containing the coordinates of the nodes
   */
  template<class MyMeshTypeT, class MyMeshTypeS>
  inline const double* SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::getCoordsOfSubNode(typename MyMeshTypeT::MyConnType node)
  {
    // replace "at()" with [] for unsafe but faster access
    return _nodes.at(node);
  }

  /**
   * Accessor to the coordinates of a given (sub)-node
   *
   * @param  node  local number of the (sub)-node 0,..,7 are the elements nodes, sub-nodes are numbered from 8,..
   * @param nodeId is an output that is node id in target whole mesh in C mode.
   * @return pointer to double[3] containing the coordinates of the nodes
   */
  template<class MyMeshTypeT, class MyMeshTypeS>
  const double* SplitterTetra2<MyMeshTypeT, MyMeshTypeS>::getCoordsOfSubNode2(typename MyMeshTypeT::MyConnType node, typename MyMeshTypeT::MyConnType& nodeId)
  {
    const double *ret=_nodes.at(node);
    if(node<8)
      nodeId=_node_ids[node];
    else
      nodeId=-1;
    return ret;    
  }
}

#endif
