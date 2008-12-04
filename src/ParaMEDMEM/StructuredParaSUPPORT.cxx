//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#include "Topology.hxx"
#include "BlockTopology.hxx"
#include "ParaGRID.hxx"
#include "ParaMESH.hxx"
#include "StructuredParaSUPPORT.hxx"
#include "MEDMEM_Support.hxx"

namespace ParaMEDMEM 
{
	
/*! Constructor on all elements from a GRID */
StructuredParaSUPPORT::StructuredParaSUPPORT(const ParaGRID* const grid, const MED_EN::medEntityMesh entity):
_block_topology(grid->getBlockTopology()),
_grid(grid), 
_entity(entity)
{
  _support=new SUPPORT(grid->getGrid(), "support on all entities", entity);
}
/*! Constructor on all elements from a GRID */
StructuredParaSUPPORT::StructuredParaSUPPORT(const ParaMESH* const mesh, const MED_EN::medEntityMesh entity):
_block_topology(mesh->getBlockTopology()),
_grid(0),
_entity(entity)
{
  _mesh=mesh;
  _support=new SUPPORT(mesh->getMesh(), "support on all entities", entity);
}

StructuredParaSUPPORT::~StructuredParaSUPPORT()
{
	delete _support;
}

}//end of namespace ParaMEDMEM
