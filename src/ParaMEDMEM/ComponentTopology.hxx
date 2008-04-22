#ifndef COMPONENTTOPOLOGY_HXX_
#define COMPONENTTOPOLOGY_HXX_

#include <vector>
#include "ProcessorGroup.hxx"
#include "Topology.hxx"


namespace ParaMEDMEM
{
class ComponentTopology
{
public:
	ComponentTopology(int nb_comp, ProcessorGroup* group);
	ComponentTopology(int nb_comp, int nb_blocks);
	ComponentTopology(int nb_comp);
	ComponentTopology();
	virtual ~ComponentTopology();
	//!returns the number of MED components in the topology
	int nbComponents() const {return component_array[component_array.size()-1];}
	//!returns the number of MED components on local processor
	int nbLocalComponents() const ;
	//!returns the number of the first MED component on local processor
	int firstLocalComponent() const ;
	//!returns the number of blocks in the topology
	int nbBlocks()const {return component_array.size()-1;}
	//!returns the block structure
	const vector <int> * getBlockIndices() const {return &component_array;}
	const ProcessorGroup* getProcGroup()const {return _proc_group;} 
private:
	std::vector<int> component_array;
	ProcessorGroup* _proc_group;
};

}

#endif /*COMPONENTTOPOLOGY_HXX_*/
