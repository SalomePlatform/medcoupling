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
#ifndef __INTERSECTORHEXA_HXX__
#define __INTERSECTORHEXA_HXX__

#include "TargetIntersector.hxx"
#include "IntersectorTetra.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpolationOptions.hxx"

namespace INTERP_KERNEL
{


  /** 
   * \brief Class responsible for calculating intersection between a hexahedron target element and  
   * the source elements.
   *
   */
  template<class MyMeshType>
  class INTERPKERNEL_EXPORT IntersectorHexa : public TargetIntersector<typename MyMeshType::MyConnType>,
                                              public InterpolationOptions
  {
  public:

    IntersectorHexa(const MyMeshType& srcMesh, const MyMeshType& targetMesh,
                    typename MyMeshType::MyConnType targetCell, SplittingPolicy policy = GENERAL_24);

    ~IntersectorHexa();

    virtual double intersectSourceCell(typename MyMeshType::MyConnType srcCell);

  private:

    void fiveSplit(const MyMeshType& srcMesh, const int* const subZone);

    void sixSplit(const MyMeshType& srcMesh, const int* const subZone);

    void calculateGeneral24Tetra(const MyMeshType& srcMesh);

    void calculateGeneral48Tetra(const MyMeshType& srcMesh);

    void calculateSubNodes(const MyMeshType& targetMesh, typename MyMeshType::MyConnType targetCell, SplittingPolicy policy);
    
    inline const double* getCoordsOfSubNode(typename MyMeshType::MyConnType node);
    
    template<int n>
    inline void calcBarycenter(double* barycenter, const typename MyMeshType::MyConnType* pts);

    /// pointers to the IntersectorTetra objects representing the tetrahedra 
    /// that result from the splitting of the hexahedron
    std::vector< IntersectorTetra<MyMeshType>* > _tetra;
    
    /// vector of pointers to double[3] containing the coordinates of the
    /// (sub) - nodes
    std::vector<const double*> _nodes;
    
  };

  /**
   * Accessor to the coordinates of a given (sub)-node
   *
   * @param  node  local number of the (sub)-node 1,..,8 are the elements nodes, sub-nodes are numbered from 9,..
   * @return pointer to double[3] containing the coordinates of the nodes
   */
  template<class MyMeshType>
  inline const double* IntersectorHexa<MyMeshType>::getCoordsOfSubNode(typename MyMeshType::MyConnType node)
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
  template<class MyMeshType>
  template<int n>
  inline void IntersectorHexa<MyMeshType>::calcBarycenter(double* barycenter, const typename MyMeshType::MyConnType* pts)
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
