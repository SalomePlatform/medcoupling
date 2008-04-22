#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Support.hxx"
#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "ParaSUPPORT.hxx"
#include "ParaMESH.hxx"
#include "Topology.hxx"
#include "ExplicitTopology.hxx"
#include "BlockTopology.hxx"
#include "ComponentTopology.hxx"

#include <vector>
#include <algorithm>

using namespace std;
using namespace MEDMEM;
namespace ParaMEDMEM
{

ExplicitTopology::ExplicitTopology(const ParaSUPPORT& parasupport ):
_proc_group(parasupport.getMesh()->getBlockTopology()->getProcGroup()),
_nb_components(1)
{
	_nb_elems=parasupport.getSupport()->getNumberOfElements(MED_EN::MED_ALL_ELEMENTS);
	MED_EN::medEntityMesh entity= parasupport.getSupport()->getEntity();
	const int* global=parasupport.getMesh()->getGlobalNumbering(entity);
	_loc2glob=new int[_nb_elems]; 
	
	if (parasupport.getSupport()->isOnAllElements())
	{
		for (int i=0; i<_nb_elems; i++)
		{
			_loc2glob[i]=global[i];
			_glob2loc[global[i]]=i;
		}
	}
	else
	{
		const int* number= parasupport.getSupport()->getNumber(MED_EN::MED_ALL_ELEMENTS);	
		for (int i=0; i<_nb_elems; i++)
		{
			int local=number[i];
			_loc2glob[i]=global[local];
			_glob2loc[global[local]]=i;
		}
	}
}

ExplicitTopology::ExplicitTopology(const ExplicitTopology& topo, int nb_components)
{
  _proc_group = topo._proc_group;
  _nb_elems = topo._nb_elems;
  _nb_components = nb_components;
  _loc2glob=new int[_nb_elems];
  for (int i=0; i<_nb_elems; i++)
    {
      _loc2glob[i]=topo._loc2glob[i];
    }
  _glob2loc=topo._glob2loc;
}


ExplicitTopology::~ExplicitTopology()
{
  if (_loc2glob != 0) delete[] _loc2glob;
}


/*! Serializes the data contained in the Explicit Topology
 * for communication purposes*/
void ExplicitTopology::serialize(int* & serializer, int& size) const 
{
	vector <int> buffer;
	
	buffer.push_back(_nb_elems);
	for (int i=0; i<_nb_elems; i++)
	{
	  buffer.push_back(_loc2glob[i]);
	}
		
	serializer=new int[buffer.size()];
	size=	buffer.size();
	copy(buffer.begin(), buffer.end(), serializer);
	
}
/*! Unserializes the data contained in the Explicit Topology
 * after communication. Uses the same structure as the one used for serialize()
 * 
 * */
void ExplicitTopology::unserialize(const int* serializer,const CommInterface& comm_interface)
{
	const int* ptr_serializer=serializer;
	cout << "unserialize..."<<endl;
	_nb_elems=*ptr_serializer++;
	cout << "nbelems "<<_nb_elems<<endl;
	_loc2glob=new int[_nb_elems];
	for (int i=0; i<_nb_elems; i++)
	{
	  _loc2glob[i]=*ptr_serializer;
	  _glob2loc[*ptr_serializer]=i;
	  ptr_serializer++;
	  
	}

}

}
