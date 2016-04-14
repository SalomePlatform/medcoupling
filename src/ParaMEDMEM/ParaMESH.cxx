//
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
            DataArrayInt *CorrespElt_local2global, DataArrayInt *CorrespFace_local2global,
            DataArrayInt *CorrespNod_local2global, const ProcessorGroup& proc_group ):
    _cell_mesh(subdomain_mesh),
    _face_mesh(subdomain_face),
    _my_domain_id(proc_group.myRank()),
    _block_topology (new BlockTopology(proc_group, subdomain_mesh->getNumberOfCells())),
    _explicit_topology(0),
    _node_global(CorrespNod_local2global),
    _face_global(CorrespFace_local2global),
    _cell_global(CorrespElt_local2global)
  {
    if(_cell_mesh)
      _cell_mesh->incrRef();
    if(_face_mesh)
      _face_mesh->incrRef();
    if(CorrespElt_local2global)
      CorrespElt_local2global->incrRef();
    if(CorrespFace_local2global)
      CorrespFace_local2global->incrRef();
    if(CorrespNod_local2global)
      CorrespNod_local2global->incrRef();
  }

  ParaMESH::ParaMESH( MEDCouplingPointSet *mesh, const ProcessorGroup& proc_group, const std::string& name):
    _cell_mesh(mesh),
    _face_mesh(0),
    _my_domain_id(proc_group.myRank()),
    _block_topology (new BlockTopology(proc_group, mesh->getNumberOfCells())),
    _node_global(0),
    _face_global(0)
  {
    if(_cell_mesh)
      _cell_mesh->incrRef();
    int nb_elem=mesh->getNumberOfCells();
    _explicit_topology=new BlockTopology(proc_group,nb_elem);
    int nbOfCells=mesh->getNumberOfCells();
    _cell_global = DataArrayInt::New();
    _cell_global->alloc(nbOfCells,1);
    int *cellglobal=_cell_global->getPointer();
    int offset = _block_topology->localToGlobal(make_pair(_my_domain_id,0));
    for (int i=0; i<nbOfCells; i++)
      {
        cellglobal[i]=offset+i;
      }
  }

  void ParaMESH::setNodeGlobal(DataArrayInt *nodeGlobal)
  {
    if(nodeGlobal!=_node_global)
      {
        if(_node_global)
          _node_global->decrRef();
        _node_global=nodeGlobal;
        if(_node_global)
          _node_global->incrRef();
      }
  }

  void ParaMESH::setCellGlobal(DataArrayInt *cellGlobal)
  {
    if(cellGlobal!=_cell_global)
      {
        if(_cell_global)
          _cell_global->decrRef();
        _cell_global=cellGlobal;
        if(_cell_global)
          _cell_global->incrRef();
      }
  }

  ParaMESH::~ParaMESH()
  {
    if(_cell_mesh)
      _cell_mesh->decrRef();
    if(_face_mesh)
      _face_mesh->decrRef();
    delete _block_topology;
    if(_node_global)
      _node_global->decrRef();
    if(_cell_global)
      _cell_global->decrRef();
    if(_face_global)
      _face_global->decrRef();
    delete _explicit_topology;
  }

}
