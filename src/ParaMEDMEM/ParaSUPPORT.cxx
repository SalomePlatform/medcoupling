#include "ParaSUPPORT.hxx"
#include "ParaMESH.hxx"
#include "MEDMEM_Support.hxx"

namespace ParaMEDMEM
{

  ParaSUPPORT::ParaSUPPORT()
  {
  }
  
  ParaSUPPORT::ParaSUPPORT(const MEDMEM::SUPPORT& support, const ProcessorGroup& proc_group):
  _support(&support), 
  _has_support_ownership(false),
  _has_mesh_ownership(true)
   {
    _mesh = new ParaMESH(*(support.getMesh()),  proc_group, "mesh from support");
  } 

  ParaSUPPORT::~ParaSUPPORT()
  {
		if (_has_support_ownership)
			delete _support;
    if (_has_mesh_ownership)
      delete _mesh;
  }

	const int* ParaSUPPORT::getGlobalNumbering() const
	{
		if (! _support->isOnAllElements())
			throw MEDMEM::MEDEXCEPTION("GlobalNumbering can only be retrieved on supports on all elements");
		return _mesh->getGlobalNumbering(_support->getEntity());
	}

}

