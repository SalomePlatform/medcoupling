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
