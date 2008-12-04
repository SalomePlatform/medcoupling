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
#ifndef BLOCKTOPOLOGY_HXX_
#define BLOCKTOPOLOGY_HXX_
#include <vector>
#include <utility>
#include <iostream>
#include "ProcessorGroup.hxx"

using namespace std;
namespace MEDMEM
{
	class GRID;
}

namespace ParaMEDMEM
{
class Topology;
class ComponentTopology;
typedef enum{Block,Cycle} CYCLE_TYPE; 

class BlockTopology: public Topology
{
public:
	BlockTopology(){};
	BlockTopology(const ProcessorGroup& group, const MEDMEM::GRID& grid); 
	BlockTopology(const BlockTopology& geom_topo, const ComponentTopology& comp_topo);
	BlockTopology(const ProcessorGroup& group, int nb_elem);
	virtual ~BlockTopology();
	inline int getNbElements()const;
	inline int getNbLocalElements() const;
	const ProcessorGroup* getProcGroup()const {return _proc_group;};
	inline std::pair<int,int> globalToLocal (const int) const ;
	inline int localToGlobal (const std::pair<int,int>) const;
	std::vector<std::pair<int,int> > getLocalArrayMinMax() const ;
	int getDimension() const {return _dimension;};
	void serialize(int* & serializer, int& size) const ;
	void unserialize(const int* serializer, const CommInterface& comm_interface);
private:
	//dimension : 2 or 3
	int _dimension;
	//proc array
	std::vector<int> _nb_procs_per_dim;
	//stores the offsets vector  
	std::vector<std::vector<int> > _local_array_indices;
	//stores the cycle type (block or cyclic)
	std::vector<CYCLE_TYPE> _cycle_type;
	//Processor group
	const ProcessorGroup* _proc_group;
	//nb of elements
	int _nb_elems;
	bool _owns_processor_group;
};

//!converts a pair <subdomainid,local> to a global number 
inline std::pair<int,int> BlockTopology::globalToLocal(const int global) const {
	int subdomain_id=0;
	//int local=global;
	int position=global;
	int size=_nb_elems;
	int size_procs=_proc_group->size();
	int increment=size;
	vector<int>axis_position(_dimension);
	vector<int>axis_offset(_dimension);
	for (int idim=0; idim<_dimension; idim++)
	{
		int axis_size=_local_array_indices[idim].size()-1;
		int axis_nb_elem=_local_array_indices[idim][axis_size];
		increment=increment/axis_nb_elem;
		//cout << "increment "<<increment<<endl;
		int proc_increment = size_procs/(axis_size);
		int axis_pos=position/increment;
		position=position%increment;
//		if (_cycle_type[idim]==Block)
	//	{	
			int iaxis=1;
		//	cout << "local array "<<_local_array_indices[idim][iaxis]<<" "<<axis_pos<<endl;
			while (_local_array_indices[idim][iaxis]<=axis_pos)
			{
				subdomain_id+=proc_increment;
					iaxis++;
			}
			axis_position[idim]=axis_pos-_local_array_indices[idim][iaxis-1];
			axis_offset[idim]=iaxis;
			
//		}
		
	//	else
		//{
//			int size = axis_nb_elem/axis_size;
//			if ((position%axis_size)<(axis_nb_elem%axis_size))
//				size+=1;	
//			subdomain_id+=proc_increment*(position%axis_size);
//			local -= (axis_nb_elem-size)*increment;
		//}
	}
	int local=0;
	int local_increment=1;
	for (int idim=_dimension-1; idim>=0; idim--)
	{
		local+=axis_position[idim]*local_increment;
		local_increment*=_local_array_indices[idim][axis_offset[idim]]-_local_array_indices[idim][axis_offset[idim]-1];
	}	
	
	return make_pair(subdomain_id,local);
}

//!converts local number to a global number
inline int BlockTopology::localToGlobal(const pair<int,int> local) const {
	
	int subdomain_id=local.first;
	int global=0;
	int loc=local.second;
	int increment=_nb_elems;
	int proc_increment=_proc_group->size();
	int local_increment=getNbLocalElements();
	for (int idim=0; idim < _dimension; idim++)
	{
		int axis_size=_local_array_indices[idim].size()-1;
		int axis_nb_elem=_local_array_indices[idim][axis_size];
		increment=increment/axis_nb_elem;
		proc_increment = proc_increment/(axis_size);
		int proc_axis=subdomain_id/proc_increment;
		subdomain_id=subdomain_id%proc_increment;
		int local_axis_nb_elem=_local_array_indices[idim][proc_axis+1]-_local_array_indices[idim][proc_axis];
		local_increment = local_increment/local_axis_nb_elem;
	//	if (_cycle_type[idim]==Block)
		//{	
			int iaxis=loc/local_increment+_local_array_indices[idim][proc_axis];
			global+=increment*iaxis;
			loc = loc%local_increment;
		//}
		//else
		//{
			//cout << "cyclic Not implemented yet"<<endl;
			//exit (2); 
		//}
	}
	return global;
	
}

//!Retrieves the number of elements for a given topology
inline int BlockTopology::getNbElements()const {return _nb_elems;}

//Retrieves the local number of elements 
inline int BlockTopology::getNbLocalElements()const 
{
	int position=_proc_group->myRank();
	int nb_elem = 1;
	int increment=1;
	for (int i=_dimension-1; i>=0; i--)
		{	
			increment *=_nb_procs_per_dim[i];
			int idim=position%increment;
			position=position/increment;
				//cout << "i idim dimension"<<i<<" "<<idim<<" "<<_dimension<<endl;
			int imin=_local_array_indices[i][idim];
			int imax=_local_array_indices[i][idim+1];
//			cout << "position imax imin "<<position<<" "<< imax <<" "<< imin<< " "<<idim<<endl;
			nb_elem*=(imax-imin);
		}
	return nb_elem;
}
}


#endif /*BLOCKTOPOLOGY_HXX_*/
