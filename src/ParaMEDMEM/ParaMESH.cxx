//
// Copyright (C) 2007-2019  CEA/DEN, EDF R&D
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

#include "ParaMESH.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "Topology.hxx"
#include "BlockTopology.hxx"
#include "MEDCouplingMemArray.hxx"

#include <fstream>
#include <vector>

//inclusion for the namespaces
using namespace std;

namespace MEDCoupling
{
  ParaMESH::ParaMESH( MEDCouplingPointSet *subdomain_mesh, MEDCouplingPointSet *subdomain_face,
            DataArrayIdType *CorrespElt_local2global, DataArrayIdType *CorrespFace_local2global,
            DataArrayIdType *CorrespNod_local2global, const ProcessorGroup& proc_group ):
    _my_domain_id(proc_group.myRank()),
    _block_topology(new BlockTopology(proc_group, subdomain_mesh->getNumberOfCells())),
    _explicit_topology(nullptr)
  {
    _cell_mesh.takeRef(subdomain_mesh);
    _face_mesh.takeRef(subdomain_face);
    _node_global.takeRef(CorrespNod_local2global);
    _face_global.takeRef(CorrespFace_local2global);
    _cell_global.takeRef(CorrespElt_local2global);
  }

  ParaMESH::ParaMESH( MEDCouplingPointSet *mesh, const ProcessorGroup& proc_group, const std::string& name):
    _my_domain_id(proc_group.myRank()),
    _block_topology(new BlockTopology(proc_group, mesh->getNumberOfCells()))
  {
    _cell_mesh.takeRef(mesh);
    mcIdType nb_elem=mesh->getNumberOfCells();
    _explicit_topology=new BlockTopology(proc_group,nb_elem);
    mcIdType nbOfCells=mesh->getNumberOfCells();
    _cell_global = DataArrayIdType::New();
    _cell_global->alloc(nbOfCells,1);
    mcIdType *cellglobal=_cell_global->getPointer();
    mcIdType offset = _block_topology->localToGlobal(make_pair(_my_domain_id,0));
    for (mcIdType i=0; i<nbOfCells; i++)
      {
        cellglobal[i]=offset+i;
      }
  }

  void ParaMESH::setNodeGlobal(DataArrayIdType *nodeGlobal)
  {
    _node_global.takeRef(nodeGlobal);
  }

  void ParaMESH::setCellGlobal(DataArrayIdType *cellGlobal)
  {
    _cell_global.takeRef(cellGlobal);
  }

  ParaMESH::~ParaMESH()
  {
    delete _block_topology;
    delete _explicit_topology;
  }

}
