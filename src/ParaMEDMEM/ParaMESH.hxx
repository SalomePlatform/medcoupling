// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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

#ifndef __PARAMESH_HXX__
#define __PARAMESH_HXX__

#include "MEDCouplingPointSet.hxx"
#include "ProcessorGroup.hxx"
#include "MEDCouplingMemArray.hxx"

#include <string>
#include <vector>

namespace ParaMEDMEM
{
  class Topology;
  class BlockTopology;
  class DataArrayInt;

  class ParaMESH
  {
  public:
    ParaMESH( MEDCouplingPointSet *subdomain_mesh,
              MEDCouplingPointSet *subdomain_face,
              DataArrayInt *CorrespElt_local2global,
              DataArrayInt *CorrespFace_local2global,
              DataArrayInt *CorrespNod_local2global,
              const ProcessorGroup& proc_group ) ;
    ParaMESH( MEDCouplingPointSet *mesh,
              const ProcessorGroup& proc_group, const std::string& name);

    virtual ~ParaMESH();
    void setNodeGlobal(DataArrayInt *nodeGlobal);
    void setCellGlobal(DataArrayInt *cellGlobal);
    Topology* getTopology() const { return _explicit_topology; }
    bool isStructured() const { return _cell_mesh->isStructured(); }
    MEDCouplingPointSet *getCellMesh() const { return _cell_mesh; }
    MEDCouplingPointSet *getFaceMesh() const { return _face_mesh; }
    BlockTopology* getBlockTopology() const { return _block_topology; }

    DataArrayInt* getGlobalNumberingNodeDA() const { if(_node_global) _node_global->incrRef(); return _node_global; }
    DataArrayInt* getGlobalNumberingFaceDA() const { if(_face_global) _face_global->incrRef(); return _face_global; }
    DataArrayInt* getGlobalNumberingCellDA() const { if(_cell_global) _cell_global->incrRef(); return _cell_global; }
    const int* getGlobalNumberingNode() const { if(_node_global) return _node_global->getConstPointer(); return 0; }
    const int* getGlobalNumberingFace() const { if(_face_global) return _face_global->getConstPointer(); return 0; }
    const int* getGlobalNumberingCell() const { if(_cell_global) return _cell_global->getConstPointer(); return 0; }

  private:
    //mesh object underlying the ParaMESH object
    MEDCouplingPointSet *_cell_mesh ;
    MEDCouplingPointSet *_face_mesh ;

    //id of the local grid
    int _my_domain_id;

    //global topology of the cells
    ParaMEDMEM::BlockTopology* _block_topology;
    Topology*  _explicit_topology;
    // pointers to global numberings
    DataArrayInt* _node_global;
    DataArrayInt* _face_global;
    DataArrayInt* _cell_global;
  };
}

#endif
