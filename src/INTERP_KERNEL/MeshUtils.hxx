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

#ifndef __MESHUTILS_HXX__
#define __MESHUTILS_HXX__

#include "InterpolationUtils.hxx"

namespace INTERP_KERNEL
{
  /**
   * Returns the global number of the node of an element.
   *
   * @param      node       the node for which the global number is sought (ALWAYS in C mode)
   * @param      element    an element of the mesh (in numPol policy)
   * @param      mesh       a mesh
   * @return    the node's global number so that (its coordinates in the coordinates array are at [SPACEDIM*globalNumber, SPACEDIM*globalNumber + SPACEDIM]
   */
  template<class MyMeshType>
  inline typename MyMeshType::MyConnType getGlobalNumberOfNode(typename MyMeshType::MyConnType node, typename MyMeshType::MyConnType element, const MyMeshType& mesh)
  {
    typedef typename MyMeshType::MyConnType ConnType;
    const NumberingPolicy numPol=MyMeshType::My_numPol;
    const ConnType elemIdx=OTT<ConnType,numPol>::conn2C(mesh.getConnectivityIndexPtr()[OTT<ConnType,numPol>::ind2C(element)]);
    if(mesh.getTypeOfElement(element)!=INTERP_KERNEL::NORM_POLYHED)
      return OTT<ConnType,numPol>::coo2C(mesh.getConnectivityPtr()[elemIdx + node]);
    else
      {
        const ConnType *startNodalConnOfElem=mesh.getConnectivityPtr()+elemIdx;
        ConnType ptr=0,ret=0;
        while(startNodalConnOfElem[ret]==-1 || ptr!=node)
          {
            ret++;
            if(startNodalConnOfElem[ret]!=-1)
              ptr++;
          }
        return OTT<ConnType,numPol>::coo2C(startNodalConnOfElem[ret]);
      }
  }

  /**
   * Returns the coordinates of a node of an element
   *
   * @param      node       the node for which the coordinates are sought. In C mode.
   * @param      element    an element of the mesh. In mesh policy.
   * @param      mesh       a mesh
   * @return    pointer to an array of 3 doubles containing the coordinates of the node
   */
  template<class MyMeshType>
  inline const double* getCoordsOfNode(typename MyMeshType::MyConnType node, typename MyMeshType::MyConnType element, const MyMeshType& mesh)
  {
    typedef typename MyMeshType::MyConnType ConnType;
    const ConnType connIdx = getGlobalNumberOfNode(node, element, mesh);
    const double *ret=mesh.getCoordinatesPtr()+MyMeshType::MY_SPACEDIM*connIdx;
    return ret;
  }

  /**
   * Returns the coordinates of a node of an element
   *
   * @param      node       the node for which the coordinates are sought. In C mode.
   * @param      element    an element of the mesh. In mesh policy.
   * @param      mesh       a mesh
   * @param      nodeId     globale nodeId in whole mesh point of view in C mode.
   * @return    pointer to an array of 3 doubles containing the coordinates of the node
   */
  template<class MyMeshType>
  inline const double* getCoordsOfNode2(typename MyMeshType::MyConnType node, typename MyMeshType::MyConnType element, const MyMeshType& mesh, typename MyMeshType::MyConnType& nodeId)
  {
    nodeId= getGlobalNumberOfNode(node, element, mesh);
    return mesh.getCoordinatesPtr()+MyMeshType::MY_SPACEDIM*nodeId;
  }

  /**
   * Returns the barycentric coordinates of a point within a triangle or tetrahedron
   *
   * @param      point       the point for which the barycentric coordinates are sought
   * @param      element    an element of the mesh
   * @param      mesh       a mesh
   * @param     barycentricCoords an array of 3 doubles containing the coordinates of the node
   */
  template<class MyMeshType, int NB_NODES>
  inline void getBarycentricCoordinates(const double*                   point,
                                        typename MyMeshType::MyConnType element,
                                        const MyMeshType&               mesh,
                                        double*                         barycentricCoords)
  {
    std::vector<const double*> nodes( NB_NODES );
    for ( int node = 0; node < NB_NODES; ++node )
      {
        nodes[ node ] = getCoordsOfNode( node, element, mesh );
      }
    barycentric_coords( nodes, point, barycentricCoords );
  }
    
}

#endif
