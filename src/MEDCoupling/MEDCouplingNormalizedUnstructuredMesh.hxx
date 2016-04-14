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
// Author : Anthony Geay (CEA/DEN)

#ifndef __PARAMEDMEM_MEDCOUPLINGNORMALIZEDUNSTRUCTUREDMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGNORMALIZEDUNSTRUCTUREDMESH_HXX__

#include "NormalizedUnstructuredMesh.hxx"

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
  typedef int MyConnType;
  static const INTERP_KERNEL::NumberingPolicy My_numPol=INTERP_KERNEL::ALL_C_MODE;
public:
  MEDCouplingNormalizedUnstructuredMesh(const MEDCoupling::MEDCouplingPointSet *mesh);
  void getBoundingBox(double *boundingBox) const;
  INTERP_KERNEL::NormalizedCellType getTypeOfElement(int eltId) const;
  int getNumberOfNodesOfElement(int eltId) const;
  int getNumberOfElements() const;
  int getNumberOfNodes() const;
  const int *getConnectivityPtr() const;
  const double *getCoordinatesPtr() const;
  const int *getConnectivityIndexPtr() const;
  void releaseTempArrays();
  ~MEDCouplingNormalizedUnstructuredMesh();
private:
  void prepare();
private:
  const MEDCoupling::MEDCouplingPointSet *_mesh;
  int *_conn_for_interp;
  int *_conn_index_for_interp;
};

#endif
