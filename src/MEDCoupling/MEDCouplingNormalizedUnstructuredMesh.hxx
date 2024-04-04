// Copyright (C) 2007-2024  CEA, EDF
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
// Author : Anthony Geay (CEA/DEN)

#ifndef __PARAMEDMEM_MEDCOUPLINGNORMALIZEDUNSTRUCTUREDMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGNORMALIZEDUNSTRUCTUREDMESH_HXX__

#include "MCIdType.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "NormalizedGeometricTypes"

namespace MEDCoupling
{
  class MEDCouplingPointSet;
}

template<int SPACEDIM,int MESHDIM>
class MEDCouplingNormalizedUnstructuredMesh
{
public:
  static const int MY_SPACEDIM=SPACEDIM;
  static const int MY_MESHDIM=MESHDIM;
  using MyConnType = mcIdType;
  static const INTERP_KERNEL::NumberingPolicy My_numPol=INTERP_KERNEL::ALL_C_MODE;
public:
  MEDCouplingNormalizedUnstructuredMesh(const MEDCoupling::MEDCouplingPointSet *mesh);
  void getBoundingBox(double *boundingBox) const;
  INTERP_KERNEL::NormalizedCellType getTypeOfElement(mcIdType eltId) const;
  mcIdType getNumberOfNodesOfElement(mcIdType eltId) const;
  mcIdType getNumberOfElements() const;
  mcIdType getNumberOfNodes() const;
  const mcIdType *getConnectivityPtr() const;
  const double *getCoordinatesPtr() const;
  const mcIdType *getConnectivityIndexPtr() const;
  void releaseTempArrays();
  ~MEDCouplingNormalizedUnstructuredMesh();
private:
  void prepare();
private:
  const MEDCoupling::MEDCouplingPointSet *_mesh;
  mcIdType *_conn_for_interp;
  mcIdType *_conn_index_for_interp;
};

#endif
