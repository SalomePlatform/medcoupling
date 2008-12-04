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
#ifndef ExplicitTOPOLOGY_HXX_
#define ExplicitTOPOLOGY_HXX_
#include <vector>
#include <utility>
#include <iostream>
#include "ProcessorGroup.hxx"
#include <ext/hash_map>

using namespace std;
using namespace __gnu_cxx;

namespace MEDMEM
{
	class GRID;
}

namespace ParaMEDMEM
{
class Topology;
class ComponentTopology;
  class ParaSUPPORT;

class ExplicitTopology: public Topology
{
public:
	ExplicitTopology(){};
  ExplicitTopology( const ExplicitTopology& topo, int nbcomponents);
	ExplicitTopology(const ParaSUPPORT&);
	virtual ~ExplicitTopology();
	
	inline int getNbElements()const;
	inline int getNbLocalElements() const;
	const ProcessorGroup* getProcGroup()const {return _proc_group;};
// 	inline std::pair<int,int> globalToLocal (const int global) const {
// 	pair <int,int>local;
// 	local.first=_proc_group->myRank();
// 	local.second=globalToLocal(global);}
    int localToGlobal (const std::pair<int,int> local) const {return localToGlobal(local.second);}
	inline int localToGlobal(int) const;
	inline int globalToLocal(int) const;
	void serialize(int* & serializer, int& size) const ;
	void unserialize(const int* serializer, const CommInterface& comm_interface);
  int getNbComponents () const {return _nb_components;}
private:
  //Processor group
  const ProcessorGroup* _proc_group;
  //nb of elements
  int _nb_elems;
  //nb of components
  int _nb_components;
  //mapping local to global
  int* _loc2glob;
  //mapping global to local
  hash_map<int,int> _glob2loc;
};

//!converts a pair <subdomainid,local> to a global number 
inline int ExplicitTopology::globalToLocal(const int global) const {
	return (_glob2loc.find(global))->second;;
	}

//!converts local number to a global number
int ExplicitTopology::localToGlobal(int local) const {
	return _loc2glob[local];
	}

//!Retrieves the number of elements for a given topology
inline int ExplicitTopology::getNbElements()const {return _nb_elems;}

//Retrieves the local number of elements 
inline int ExplicitTopology::getNbLocalElements()const 
{
	return _glob2loc.size();
}
}


#endif /*ExplicitTOPOLOGY_HXX_*/
