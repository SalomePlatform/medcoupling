
#ifndef _ICOCOMEDFIELD_HXX_
#define _ICOCOMEDFIELD_HXX_

#include <vector>
#include "ICoCoField.hxx"
#include "MEDMEM_Field.hxx"

namespace ParaMEDMEM
{
	class ParaMESH;
	class ParaFIELD;
	class ParaSUPPORT;
	class ComponentTopology;
	class ProcessorGroup;
}
namespace MEDMEM
{
	class CONNECTZONE;
	class MESHING;
	class SUPPORT;
}
namespace ICoCo
{
	class TrioField;
	
  class MEDField:public ICoCo::Field
    {
    public:
      
      MEDField(){};
      MEDField(ParaMEDMEM::ParaMESH* mesh, ParaMEDMEM::ParaFIELD* field);
			MEDField(TrioField& , const ParaMEDMEM::ProcessorGroup& group);
			virtual ~MEDField();
      ParaMEDMEM::ParaFIELD* getField() const  {return _field;}
      ParaMEDMEM::ParaMESH* getMesh()const {return _mesh;}

    private:
      ParaMEDMEM::ParaMESH* _mesh;
      ParaMEDMEM::ParaFIELD* _field;
      MEDMEM::FIELD<double>* _local_field;
			bool _has_field_ownership;
			MEDMEM::MESHING* _local_mesh;
			MEDMEM::SUPPORT* _local_support;
			ParaMEDMEM::ParaSUPPORT* _support;
			ParaMEDMEM::ComponentTopology* _comp_topology;
			
      std::vector<MEDMEM::CONNECTZONE*> _connectzones;
    };
};

#endif
