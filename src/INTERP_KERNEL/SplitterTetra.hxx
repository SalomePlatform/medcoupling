//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#ifndef __SPLITTERTETRA_HXX__
#define __SPLITTERTETRA_HXX__

#include "TransformedTriangle.hxx"
#include "TetraAffineTransform.hxx"
#include "InterpolationOptions.hxx"

#include <assert.h>
#include <vector>
#include <functional>
#include <map>
#ifdef WIN32
# include <hash_map>
#else
# include <ext/hash_map>
#endif

#ifndef WIN32
using __gnu_cxx::hash_map;
#else
using stdext::hash_map;
using stdext::hash_compare;
#endif

namespace INTERP_KERNEL
{
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
      sort3Ints(_nodes, node1, node2, node3);
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
     * Returns a hash value for the object, based on its three nodes.
     * This value is not unique for each face.
     *
     * @return  a hash value for the object
     */
    int hashVal() const
    {
      return _hashVal;
    }

#ifdef WIN32
    operator size_t () const
    {
      return _hashVal;
    }
#endif
     
    inline void sort3Ints(int* sorted, int node1, int node2, int node3);

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
  inline void TriangleFaceKey::sort3Ints(int* sorted, int x1, int x2, int x3)
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
}
#ifndef WIN32
namespace __gnu_cxx
{


  /**
   * \brief Template specialization of __gnu_cxx::hash<T> function object for use with a __gnu_cxx::hash_map 
   * with TriangleFaceKey as key class.
   * 
   */
  template<>
  class hash<INTERP_KERNEL::TriangleFaceKey>

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
#else
  struct TriangleFaceKeyComparator
  {
    bool operator()(const INTERP_KERNEL::TriangleFaceKey& key1,
                    const INTERP_KERNEL::TriangleFaceKey& key2 ) const
    {
      return key1.hashVal() < key2.hashVal();
    }
  };
#endif

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

    ~SplitterTetra();

    double intersectSourceCell(typename MyMeshType::MyConnType srcCell, double* baryCentre=0);

    double intersectTetra(const double** tetraCorners);

    typename MyMeshType::MyConnType getId(int id) { return _conn[id]; }
    
    void splitIntoDualCells(SplitterTetra<MyMeshType> **output);

    void splitMySelfForDual(double* output, int i, typename MyMeshType::MyConnType& nodeId);

    void clearVolumesCache();

  private:
    // member functions
    inline void createAffineTransform(const double** corners);
    inline void checkIsOutside(const double* pt, bool* isOutside) const;
    inline void calculateNode(typename MyMeshType::MyConnType globalNodeNum);
    inline void calculateVolume(TransformedTriangle& tri, const TriangleFaceKey& key);
        

    /// disallow copying
    SplitterTetra(const SplitterTetra& t);
    
    /// disallow assignment
    SplitterTetra& operator=(const SplitterTetra& t);

    // member variables
    /// affine transform associated with this target element
    TetraAffineTransform* _t;
    
    /// hash_map relating node numbers to transformed nodes, used for caching
    hash_map< int , double* > _nodes;
    
    /// hash_map relating triangular faces to calculated volume contributions, used for caching
    hash_map< TriangleFaceKey, double
#ifdef WIN32
        , hash_compare<TriangleFaceKey,TriangleFaceKeyComparator> 
#endif
    > _volumes;

    /// reference to the source mesh
    const MyMeshType& _src_mesh;
                
    // node id of the first node in target mesh in C mode.
    typename MyMeshType::MyConnType _conn[4];

    double _coords[12];
  };

  /**
   * Creates the affine transform _t from the corners of the tetrahedron. Used by the constructors
   *
   * @param corners  double*[4] array containing pointers to four double[3] arrays with the 
   *                 coordinates of the corners of the tetrahedron
   */
  template<class MyMeshType>
  inline void SplitterTetra<MyMeshType>::createAffineTransform(const double** corners)
  {
    // create AffineTransform from tetrahedron
    _t = new TetraAffineTransform( corners );
  }

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
  inline void SplitterTetra<MyMeshType>::checkIsOutside(const double* pt, bool* isOutside) const
  {
    isOutside[0] = isOutside[0] && (pt[0] <= 0.0);
    isOutside[1] = isOutside[1] && (pt[0] >= 1.0);
    isOutside[2] = isOutside[2] && (pt[1] <= 0.0);
    isOutside[3] = isOutside[3] && (pt[1] >= 1.0);
    isOutside[4] = isOutside[4] && (pt[2] <= 0.0);
    isOutside[5] = isOutside[5] && (pt[2] >= 1.0);
    isOutside[6] = isOutside[6] && (1.0 - pt[0] - pt[1] - pt[2] <= 0.0);
    isOutside[7] = isOutside[7] && (1.0 - pt[0] - pt[1] - pt[2] >= 1.0);
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

  template<class MyMeshType>
  class SplitterTetra2
  {
  public:
    SplitterTetra2(const MyMeshType& targetMesh, const MyMeshType& srcMesh, SplittingPolicy policy);
    ~SplitterTetra2();
    void releaseArrays();
    void splitTargetCell(typename MyMeshType::MyConnType targetCell, typename MyMeshType::MyConnType nbOfNodesT,
                         typename std::vector< SplitterTetra<MyMeshType>* >& tetra);
    void fiveSplit(const int* const subZone, typename std::vector< SplitterTetra<MyMeshType>* >& tetra);
    void sixSplit(const int* const subZone, typename std::vector< SplitterTetra<MyMeshType>* >& tetra);
    void calculateGeneral24Tetra(typename std::vector< SplitterTetra<MyMeshType>* >& tetra);
    void calculateGeneral48Tetra(typename std::vector< SplitterTetra<MyMeshType>* >& tetra);
    void calculateSubNodes(const MyMeshType& targetMesh, typename MyMeshType::MyConnType targetCell);
    inline const double* getCoordsOfSubNode(typename MyMeshType::MyConnType node);
    inline const double* getCoordsOfSubNode2(typename MyMeshType::MyConnType node, typename MyMeshType::MyConnType& nodeId);
    template<int n>
    inline void calcBarycenter(double* barycenter, const typename MyMeshType::MyConnType* pts);
  private:
    const MyMeshType& _target_mesh;
    const MyMeshType& _src_mesh;
    SplittingPolicy _splitting_pol;
    /// vector of pointers to double[3] containing the coordinates of the
    /// (sub) - nodes of split target cell
    std::vector<const double*> _nodes;
    std::vector<typename MyMeshType::MyConnType> _node_ids;
  };

  /**
   * Calculates the barycenter of n (sub) - nodes
   *
   * @param  n  number of nodes for which to calculate barycenter
   * @param  barycenter  pointer to double[3] array in which to store the result
   * @param  pts pointer to int[n] array containing the (sub)-nodes for which to calculate the barycenter
   */
  template<class MyMeshType>
  template<int n>
  inline void SplitterTetra2<MyMeshType>::calcBarycenter(double* barycenter, const typename MyMeshType::MyConnType* pts)
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
  template<class MyMeshType>
  inline const double* SplitterTetra2<MyMeshType>::getCoordsOfSubNode(typename MyMeshType::MyConnType node)
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
  template<class MyMeshType>
  const double* SplitterTetra2<MyMeshType>::getCoordsOfSubNode2(typename MyMeshType::MyConnType node, typename MyMeshType::MyConnType& nodeId)
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
