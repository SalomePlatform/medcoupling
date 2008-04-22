#ifndef PARAFIELD_HXX_
#define PARAFIELD_HXX_

#include "MEDMEM_define.hxx"
#include "MEDMEM_GenDriver.hxx"
#include "MEDMEM_Field.hxx"
#include "ComponentTopology.hxx"

namespace MEDMEM{
	class MEDEXCEPTION;
  template <class T> class FIELD<T>;
}


namespace ParaMEDMEM
{

class ParaSUPPORT;
class ProcessorGroup;

class ParaFIELD
{
public:

	ParaFIELD(const ParaSUPPORT* support, const ComponentTopology& component_topology); 

	ParaFIELD(MEDMEM::driverTypes driver_type, const string& file_name, 
		const string& driver_name, const ComponentTopology& component_topology) 
		throw (MEDMEM::MEDEXCEPTION);
  ParaFIELD(MEDMEM::FIELD<double>* field, const ProcessorGroup& group);
  
	virtual ~ParaFIELD();
	void write(MEDMEM::driverTypes driverType, const string& fileName="", const string& meshName="");
	void synchronizeTarget( ParaMEDMEM::ParaFIELD* source_field);
	void synchronizeSource( ParaMEDMEM::ParaFIELD* target_field);
	MEDMEM::FIELD<double>* getField() const {return _field;}
  const ParaSUPPORT* getSupport() const {return _support;}
	Topology* getTopology() const {return _topology;}
	int nbComponents() const {return _component_topology.nbComponents();}
	double getVolumeIntegral(int icomp) const;
	double getL2norm()const{return -1;}
private:
	MEDMEM::FIELD<double>* getSupportVolumes(const MEDMEM::SUPPORT& support) const;
	MEDMEM::FIELD<double>* _field;
	const  ParaMEDMEM::ComponentTopology& _component_topology;
	 ParaMEDMEM::Topology* _topology; 

	const  ParaMEDMEM::ParaSUPPORT* _support;
	bool _has_field_ownership;
	bool _has_support_ownership;
};

}

#endif /*PARAFIELD_HXX_*/
