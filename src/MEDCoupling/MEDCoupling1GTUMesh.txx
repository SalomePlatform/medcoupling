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
#include "MEDCoupling1GTUMesh.hxx"

#include <sstream>

template<class MAPCLS>
void MEDCoupling::MEDCoupling1SGTUMesh::renumberNodesInConnT(const MAPCLS& newNodeNumbersO2N)
{
  getNumberOfCells();//only to check that all is well defined.
  int *begPtr(_conn->getPointer());
  int nbElt(_conn->getNumberOfTuples());
  int *endPtr(begPtr+nbElt);
  for(int *it=begPtr;it!=endPtr;it++)
    {
      auto it2(newNodeNumbersO2N.find(*it));
      if(it2!=newNodeNumbersO2N.end())
        {
          *it=(*it2).second;
        }
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::renumberNodesInConn : At pos #" << std::distance(begPtr,it) << " of nodal connectivity value is " << *it << ". Not in map !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  updateTime();
}

template<class MAPCLS>
void MEDCoupling::MEDCoupling1DGTUMesh::renumberNodesInConnT(const MAPCLS& newNodeNumbersO2N)
{
  getNumberOfCells();//only to check that all is well defined.
  //
  int nbOfTuples(_conn->getNumberOfTuples());
  int *pt(_conn->getPointer());
  for(int i=0;i<nbOfTuples;i++,pt++)
    {
      if(*pt==-1) continue;
      if(*pt>=0)
        {
          auto it(newNodeNumbersO2N.find(*pt));
          if(it!=newNodeNumbersO2N.end())
            *pt=(*it).second;
          else
            {
              std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::renumberNodesInConn : At pos #" << i << " of connectivity, node id is " << *pt << ". Not in keys of input map !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::renumberNodesInConn : error on tuple #" << i << " value is " << *pt << " ! Should be >=0 !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  //
  updateTime();
}
