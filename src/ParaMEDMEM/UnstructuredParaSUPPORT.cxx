#include "Topology.hxx"
#include "ExplicitTopology.hxx"
#include "BlockTopology.hxx"
#include "ParaMESH.hxx"
#include "UnstructuredParaSUPPORT.hxx"
#include "MEDMEM_Support.hxx"

namespace ParaMEDMEM 
{
	
/*! Constructor on from a ParaMESH and a local support*/
UnstructuredParaSUPPORT::UnstructuredParaSUPPORT(const ParaMESH* const mesh, const MEDMEM::SUPPORT* support):
ParaSUPPORT(mesh,support),
_entity(support->getEntity())
{
  //_mesh=mesh;
  //_support=support;
  _explicit_topology=new ExplicitTopology(*(ParaSUPPORT*)this);
}

/*! Constructor on from a ProcessorGroup and a local support*/
UnstructuredParaSUPPORT::UnstructuredParaSUPPORT(const MEDMEM::SUPPORT* support, const ProcessorGroup& proc_group):
	ParaSUPPORT(*support,proc_group),
  _entity(support->getEntity())
{
  //ostringstream name ("mesh associated with support ");
 // name << support->getName();
 // _mesh=new ParaMESH(*support->getMesh(), proc_group, name.str() );
  //_support=support;
  int nb_elem=_support->getNumberOfElements(MED_EN::MED_ALL_ELEMENTS);
  //  _explicit_topology=new ExplicitTopology(*(ParaSUPPORT*)this);
  _explicit_topology=new BlockTopology(proc_group,nb_elem);
}

UnstructuredParaSUPPORT::UnstructuredParaSUPPORT(const ParaMESH* const mesh, const MED_EN::medEntityMesh entity):
  ParaSUPPORT(mesh, new MEDMEM::SUPPORT(mesh->getMesh(), "support on all entities", entity)),
  _entity(entity),
  _explicit_topology(new ExplicitTopology(*this))
{
}

UnstructuredParaSUPPORT::~UnstructuredParaSUPPORT()
{
	//	if (_has_support_ownership)
	//		delete _support;
	delete _explicit_topology;
}

 

}//end of namespace ParaMEDMEM
