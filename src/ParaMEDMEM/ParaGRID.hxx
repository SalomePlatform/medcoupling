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
#ifndef PARAGRID_HXX_
#define PARAGRID_HXX_

#include "CommInterface.hxx"

#include <vector>

#include "MEDMEM_Exception.hxx"
#include "MEDMEM_define.hxx"
#include "MEDMEM_GenDriver.hxx"
#include "MEDMEM_Grid.hxx"
#include "MEDMEM_ConnectZone.hxx"
#include "ProcessorGroup.hxx"
#include "Topology.hxx"
#include "BlockTopology.hxx"

using namespace MEDMEM;
namespace ParaMEDMEM
{


const int MYRANK_ID=-100;
class ParaGRID
{
public:
	ParaGRID(MEDMEM::GRID* global_grid, Topology* topology)throw (MEDMEM::MEDEXCEPTION);
	ParaGRID(MEDMEM::driverTypes driver_type, const string& file_name, 
		const string& driver_name, int domain_id=MYRANK_ID)
	throw (MEDMEM::MEDEXCEPTION);
	
	void write(MEDMEM::driverTypes driverType, const string& fileName="")
	throw (MEDMEM::MEDEXCEPTION);
	
	ParaMEDMEM::BlockTopology * getBlockTopology() const {return _block_topology;}
	virtual ~ParaGRID();
	MEDMEM::GRID* getGrid() const {return _grid;} 
	
private:
	MEDMEM::GRID* _grid;
	//grid name 
	const string _name;
	// structured grid topology
	ParaMEDMEM::BlockTopology* _block_topology;
	// stores the x,y,z axes on the global grid
	std::vector<std::vector<double> > _global_axis;
	//id of the local grid
	int _my_domain_id;

};

}

#endif /*PARAGRID_H_*/
