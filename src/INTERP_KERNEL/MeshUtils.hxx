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
#ifndef __MESHUTILS_HXX__
#define __MESHUTILS_HXX__

#include "InterpolationUtils.hxx"

namespace INTERP_KERNEL
{
  /**
   * Returns the global number of the node of an element.
   * (1 <= node <= no. nodes of element)
   *
   * @param      node       the node for which the global number is sought (ALWAYS in C mode)
   * @param      element    an element of the mesh (in numPol policy)
   * @param      mesh       a mesh
   * @return    the node's global number so that (its coordinates in the coordinates array are at [SPACEDIM*globalNumber, SPACEDIM*globalNumber + 2]
   */
  template<class MyMeshType>
  inline typename MyMeshType::MyConnType getGlobalNumberOfNode(typename MyMeshType::MyConnType node, typename MyMeshType::MyConnType element, const MyMeshType& mesh)
  {
    typedef typename MyMeshType::MyConnType ConnType;
    const NumberingPolicy numPol=MyMeshType::My_numPol;
    const ConnType elemIdx=OTT<ConnType,numPol>::conn2C(mesh.getConnectivityIndexPtr()[OTT<ConnType,numPol>::ind2C(element)]);
    return OTT<ConnType,numPol>::coo2C(mesh.getConnectivityPtr()[elemIdx + node]);
  }

  /**
   * Returns the coordinates of a node of an element
   * (1 <= node <= no. nodes)
   *
   * @param      node       the node for which the coordinates are sought
   * @param      element    an element of the mesh
   * @param      mesh       a mesh
   * @return    pointer to an array of 3 doubles containing the coordinates of the node
   */
  template<class MyMeshType>
  inline const double* getCoordsOfNode(int node, int element, const MyMeshType& mesh)
  {
    typedef typename MyMeshType::MyConnType ConnType;
    const ConnType connIdx = getGlobalNumberOfNode(node, element, mesh);
    return mesh.getCoordinatesPtr()+MyMeshType::MY_SPACEDIM*connIdx;
  }
    
}

#endif
