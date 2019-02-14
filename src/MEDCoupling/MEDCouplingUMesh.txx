// Copyright (C) 2018-2019  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (EDF R&D)

#pragma once
#include "MEDCouplingUMesh.hxx"

#include <sstream>

template<class MAPCLS>
void MEDCoupling::MEDCouplingUMesh::renumberNodesInConnT(const MAPCLS& newNodeNumbersO2N)
{
  checkConnectivityFullyDefined();
  int *conn(getNodalConnectivity()->getPointer());
  const int *connIndex(getNodalConnectivityIndex()->getConstPointer());
  int nbOfCells(getNumberOfCells());
  for(int i=0;i<nbOfCells;i++)
    for(int iconn=connIndex[i]+1;iconn!=connIndex[i+1];iconn++)
      {
        int& node=conn[iconn];
        if(node>=0)//avoid polyhedron separator
          {
            auto it(newNodeNumbersO2N.find(node));
            if(it!=newNodeNumbersO2N.end())
              {
                node=(*it).second;
              }
            else
              {
                std::ostringstream oss; oss << "MEDCouplingUMesh::renumberNodesInConn(map) : presence in connectivity for cell #" << i << " of node #" << node << " : Not in map !";
                throw INTERP_KERNEL::Exception(oss.str());
              }
          }
      }
  _nodal_connec->declareAsNew();
  updateTime();
}
