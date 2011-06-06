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

#ifndef __PARAMEDMEM_MEDCOUPLINGNORMALIZEDCARTESIANMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGNORMALIZEDCARTESIANMESH_HXX__

#include "NormalizedUnstructuredMesh.hxx"

namespace ParaMEDMEM
{
  class MEDCouplingCMesh;
}

template<int SPACEDIM>
class MEDCouplingNormalizedCartesianMesh : public INTERP_KERNEL::GenericMesh
{
public:
  static const int MY_SPACEDIM=SPACEDIM;
  static const int MY_MESHDIM=SPACEDIM;
  typedef int MyConnType;
  static const INTERP_KERNEL::NumberingPolicy My_numPol=INTERP_KERNEL::ALL_C_MODE;
public:
  MEDCouplingNormalizedCartesianMesh(ParaMEDMEM::MEDCouplingCMesh *mesh);
  //void getBoundingBox(double *boundingBox) const;
  //INTERP_KERNEL::NormalizedCellType getTypeOfElement(int eltId) const;
  //unsigned char getNumberOfNodesOfElement(int eltId) const;
  //unsigned long getNumberOfNodes() const;
  unsigned long getNumberOfElements() const;
  unsigned long nbCellsAlongAxis(int axis) const;
  const double * getCoordsAlongAxis(int axis) const;
  ~MEDCouplingNormalizedCartesianMesh();
private:
  ParaMEDMEM::MEDCouplingCMesh *_mesh;
};

#endif
