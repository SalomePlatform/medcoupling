#include "ComponentTopology.hxx"
#include "MEDMEM_Exception.hxx"

namespace ParaMEDMEM
{
	
/* Generic constructor for \a nb_comp components equally parted
 * in \a nb_blocks blocks
 */

ComponentTopology::ComponentTopology(int nb_comp, ProcessorGroup* group):_proc_group(group)
{
	int nb_blocks=group->size();
	
	if (nb_blocks>nb_comp) throw MEDMEM::MEDEXCEPTION(
	LOCALIZED("ComponentTopology Number of components must be larger than number of blocks"));
		
	component_array.resize(nb_blocks+1);
	component_array[0]=0;
	for (int i=1; i<=nb_blocks; i++)
	{
		component_array[i]=component_array[i-1]+nb_comp/nb_blocks;
		if (i<=nb_comp%nb_blocks)
			component_array[i]++;
	}
	
}
	
/* Generic constructor for \a nb_comp components equally parted
 * in \a nb_blocks blocks
 */

ComponentTopology::ComponentTopology(int nb_comp, int nb_blocks):_proc_group(0)
{
	if (nb_blocks>nb_comp) throw MEDMEM::MEDEXCEPTION(
	LOCALIZED("ComponentTopology Number of components must be larger than number of blocks"));
		
	component_array.resize(nb_blocks+1);
	component_array[0]=0;
	for (int i=1; i<=nb_blocks; i++)
	{
		component_array[i]=component_array[i-1]+nb_comp/nb_blocks;
		if (i<=nb_comp%nb_blocks)
			component_array[i]++;
	}
	
}
//!Constructor for one block of \a nb_comp components
ComponentTopology::ComponentTopology(int nb_comp):_proc_group(0)
{
		
	component_array.resize(2);
	component_array[0]=0;
	component_array[1]=nb_comp;
	
}

//! Constructor for one component
ComponentTopology::ComponentTopology():_proc_group(0)
{
	component_array.resize(2);
	component_array[0]=0;
	component_array[1]=1;
	
}
ComponentTopology::~ComponentTopology()
{
}

int ComponentTopology::nbLocalComponents() const{
	if (_proc_group==0) return nbComponents();
	
	int nbcomp;
	int myrank = _proc_group->myRank();
	//cout << "nbLocalComp "<<myrank<<" "<< component_array[myrank+1]<< " " <<component_array[myrank]<<endl;
	if (myrank!=-1)
		nbcomp = component_array[myrank+1]-component_array[myrank];
	else 
		nbcomp=0;
	return nbcomp;
}

int ComponentTopology::firstLocalComponent() const{
		if (_proc_group==0) return 0;
	
	int icomp;
	int myrank = _proc_group->myRank();
	if (myrank!=-1)
		icomp = component_array[myrank];
	else 
		icomp=-1;
	return icomp;
}
}
