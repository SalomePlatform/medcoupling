// Copyright (C) 2020  CEA/DEN, EDF R&D
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
#include "ProcessorGroup.hxx"
#include "MEDCouplingMemArray.hxx"

#include <string>
#include <vector>

namespace MEDCoupling
{
  /*!
   * Parallel representation of an unstructured mesh.
   *
   * This class is very specific to the requirement of parallel code computations.
   */
  class ParaUMesh
  {
  public:
    ParaUMesh(MEDCouplingUMesh *mesh, DataArrayIdType *globalCellIds, DataArrayIdType *globalNodeIds);
    MCAuto<DataArrayIdType> getCellIdsLyingOnNodes(const DataArrayIdType *globalNodeIds, bool fullyIn) const;
    MCAuto<ParaUMesh> redistributeCells(const DataArrayIdType *globalCellIds) const;
    virtual ~ParaUMesh() { }
  private:
    MCAuto<MEDCouplingUMesh> _mesh;
    MCAuto<DataArrayIdType> _cell_global;
    MCAuto<DataArrayIdType> _node_global;
  };
}
