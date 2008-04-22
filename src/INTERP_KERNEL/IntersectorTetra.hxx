#ifndef __INTERSECTORTETRA_HXX__
#define __INTERSECTORTETRA_HXX__

#include "TargetIntersector.hxx"
#include "TransformedTriangle.hxx"
#include "TetraAffineTransform.hxx"
#include "NormalizedUnstructuredMesh.hxx"

#include <assert.h>
#include <vector>
#include <functional>
#include <map>
#include <ext/hash_map>

#include "MEDMEM_define.hxx"
#include "MEDMEM_CellModel.hxx"

using __gnu_cxx::hash_map;

namespace INTERP_KERNEL
{
  /**
   * \brief Class representing a triangular face, used as key in caching hash map in IntersectorTetra.
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

namespace INTERP_KERNEL
{

  /** 
   * \brief Class calculating the volume of intersection between a tetrahedral target element and
   * source elements with triangular or quadratilateral faces.
   *
   */
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  class IntersectorTetra : public TargetIntersector<ConnType>
  {

  public: 

    IntersectorTetra(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& srcMesh,
                     const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& targetMesh, ConnType targetCell);
    
    IntersectorTetra(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& srcMesh,
                     const double** tetraCorners);
    
    ~IntersectorTetra();

    double intersectSourceCell(ConnType srcCell);

  private:
    
    // member functions
    inline void createAffineTransform(const double** corners);
    inline void checkIsOutside(const double* pt, bool* isOutside) const;
    inline void calculateNode(ConnType globalNodeNum);
    inline void calculateVolume(TransformedTriangle& tri, const TriangleFaceKey& key);
	

    /// disallow copying
    IntersectorTetra(const IntersectorTetra& t);
    
    /// disallow assignment
    IntersectorTetra& operator=(const IntersectorTetra& t);

    // member variables
    /// affine transform associated with this target element
    TetraAffineTransform* _t;
    
    /// hash_map relating node numbers to transformed nodes, used for caching
    hash_map< int, double* > _nodes;
    
    /// hash_map relating triangular faces to calculated volume contributions, used for caching
    hash_map< TriangleFaceKey, double > _volumes;

    /// reference to the source mesh
    const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& _srcMesh;
		
		
  };

  /**
   * Creates the affine transform _t from the corners of the tetrahedron. Used by the constructors
   *
   * @param corners  double*[4] array containing pointers to four double[3] arrays with the 
   *                 coordinates of the corners of the tetrahedron
   */
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  inline void IntersectorTetra<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::createAffineTransform(const double** corners)
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
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  inline void IntersectorTetra<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::checkIsOutside(const double* pt, bool* isOutside) const
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
   * Gets the coordinates for the node in _srcMesh with the given global number and applies TetraAffineTransform
   * _t to it. Stores the result in the cache _nodes. The non-existance of the node in _nodes should be verified before
   * calling.
   *
   * @param globalNodeNum  global node number of the node in the mesh _srcMesh
   *
   */
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  inline void IntersectorTetra<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::calculateNode(ConnType globalNodeNum)
  {  
    const double* node = _srcMesh.getCoordinatesPtr()+SPACEDIM*globalNodeNum;
    double* transformedNode = new double[SPACEDIM];
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
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  inline void IntersectorTetra<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::calculateVolume(TransformedTriangle& tri, const TriangleFaceKey& key)
  {
    const double vol = tri.calculateIntersectionVolume();
    _volumes.insert(std::make_pair(key, vol));
  }
}

#endif
