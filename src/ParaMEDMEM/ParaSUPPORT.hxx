#ifndef PARASUPPORT_HXX_
#define PARASUPPORT_HXX_

#include "MEDMEM_Exception.hxx"
namespace MEDMEM
{
  class SUPPORT;
}
namespace ParaMEDMEM
{
  class Topology;
  class ParaMESH;
  class ProcessorGroup;
  class ParaSUPPORT
  {
  public:
    ParaSUPPORT();
    ParaSUPPORT(const ParaMESH* mesh, const MEDMEM::SUPPORT* support):
      _support(support), _mesh(mesh), _has_support_ownership(false), _has_mesh_ownership(false){};
		ParaSUPPORT(const ParaMESH* mesh):_has_mesh_ownership(false){};
		ParaSUPPORT(const MEDMEM::SUPPORT&, const ProcessorGroup& proc_group);
		virtual ~ParaSUPPORT();
    virtual const Topology* getTopology() const =0;
    virtual const MEDMEM::SUPPORT* getSupport() const {return _support;}
    virtual const ParaMESH* getMesh() const {return _mesh;}
		virtual const int* getGlobalNumbering() const;

  protected :
		
    const MEDMEM::SUPPORT* _support;
    const ParaMESH* _mesh;
		bool _has_support_ownership;
    bool _has_mesh_ownership;
  };
  
}


#endif /*PARASUPPORT_HXX_*/
