// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#ifndef __PARAMEDMEM_MEDCOUPLINGNORMALIZEDUNSTRUCTUREDMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGNORMALIZEDUNSTRUCTUREDMESH_HXX__

#include "NormalizedUnstructuredMesh.hxx"

namespace ParaMEDMEM
{
  class MEDCouplingUMesh;
}

template<int SPACEDIM,int MESHDIM>
class MEDCouplingNormalizedUnstructuredMesh : public INTERP_KERNEL::GenericMesh
{
public:
  static const int MY_SPACEDIM=SPACEDIM;
  static const int MY_MESHDIM=MESHDIM;
  typedef int MyConnType;
  static const INTERP_KERNEL::NumberingPolicy My_numPol=INTERP_KERNEL::ALL_C_MODE;
public:
  MEDCouplingNormalizedUnstructuredMesh(const ParaMEDMEM::MEDCouplingUMesh *mesh);
  void getBoundingBox(double *boundingBox) const;
  INTERP_KERNEL::NormalizedCellType getTypeOfElement(int eltId) const;
  unsigned char getNumberOfNodesOfElement(int eltId) const;
  unsigned long getNumberOfElements() const;
  unsigned long getNumberOfNodes() const;
  const int *getConnectivityPtr() const;
  const double *getCoordinatesPtr() const;
  const int *getConnectivityIndexPtr() const;
  void releaseTempArrays();
  ~MEDCouplingNormalizedUnstructuredMesh();
private:
  void prepare();
private:
  const ParaMEDMEM::MEDCouplingUMesh *_mesh;
  int *_conn_for_interp;
  int *_conn_index_for_interp;
};

#endif
