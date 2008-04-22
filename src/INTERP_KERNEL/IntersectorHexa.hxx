#ifndef __INTERSECTORHEXA_HXX__
#define __INTERSECTORHEXA_HXX__

#include "TargetIntersector.hxx"
#include "IntersectorTetra.hxx"
#include "NormalizedUnstructuredMesh.hxx"

namespace INTERP_KERNEL
{
  /** 
   * \brief Class responsible for calculating intersection between a hexahedron target element and  
   * the source elements.
   *
   */
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  class IntersectorHexa : public TargetIntersector<ConnType>
  {

  public:

    /// Type describing the different ways in which the hexahedron can be split into tetrahedra.
    /// The PLANAR_* policies persume that each face is to be considered planar, while the general
    /// policies make no such hypothesis. The integer at the end gives the number of tetrahedra
    /// that result from the split.
    enum SplittingPolicy { PLANAR_FACE_5 = 5, PLANAR_FACE_6 = 6, GENERAL_24 = 24, GENERAL_48 = 48 };

    IntersectorHexa(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& srcMesh,
                    const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& targetMesh,
                    ConnType targetCell, SplittingPolicy policy = GENERAL_24);
    
    ~IntersectorHexa();

    virtual double intersectSourceCell(ConnType srcCell);

  private:

    void fiveSplit(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& srcMesh,
                   const int* const subZone);
    
    void sixSplit(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& srcMesh,
                  const int* const subZone);

    void calculateGeneral24Tetra(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& srcMesh);

    void calculateGeneral48Tetra(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& srcMesh);
    
    void calculateSubNodes(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& targetMesh, ConnType targetCell, SplittingPolicy policy);
    
    inline const double* getCoordsOfSubNode(ConnType node);
    
    template<int n>
    inline void calcBarycenter(double* barycenter, const int* const pts);

    /// pointers to the IntersectorTetra objects representing the tetrahedra 
    /// that result from the splitting of the hexahedron
    vector< IntersectorTetra<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>* > _tetra;
    
    /// vector of pointers to double[3] containing the coordinates of the
    /// (sub) - nodes
    vector<const double*> _nodes;
    
  };

  /**
   * Accessor to the coordinates of a given (sub)-node
   *
   * @param  node  local number of the (sub)-node 1,..,8 are the elements nodes, sub-nodes are numbered from 9,..
   * @return pointer to double[3] containing the coordinates of the nodes
   */
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  inline const double* IntersectorHexa<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::getCoordsOfSubNode(ConnType node)
  {
    // replace "at()" with [] for unsafe but faster access
    return _nodes.at(node-1);//OTT<ConnType,numPol>::coo2C(node));
  }
  
  /**
   * Calculates the barycenter of n (sub) - nodes
   *
   * @param  n  number of nodes for which to calculate barycenter
   * @param  barycenter  pointer to double[3] array in which to store the result
   * @param  pts pointer to int[n] array containing the (sub)-nodes for which to calculate the barycenter
   */
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  template<int n>
  inline void IntersectorHexa<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::calcBarycenter(double* barycenter, const int* const pts)
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

}

#endif
