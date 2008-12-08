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
#ifndef PARAMESH_HXX_
#define PARAMESH_HXX_

#include "CommInterface.hxx"

#include <string>
#include <vector>

#include "MEDMEM_Exception.hxx"
#include "MEDMEM_define.hxx"
#include "MEDMEM_GenDriver.hxx"
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_ConnectZone.hxx"
#include "ProcessorGroup.hxx"

namespace ParaMEDMEM
{
class BlockTopology;

class ParaMESH
{
public:
	ParaMESH(MEDMEM::driverTypes driver_type, const std::string& file_name, 
		const ProcessorGroup& group)
	throw (MEDMEM::MEDEXCEPTION);
  ParaMESH(MEDMEM::MESH& subdomain_mesh, const ProcessorGroup& proc_group, const string& name);
	void write(MEDMEM::driverTypes driverType, const std::string& fileName="")
	throw (MEDMEM::MEDEXCEPTION);
	virtual ~ParaMESH();
	MEDMEM::MESH* getMesh() const {return _mesh;}
	ParaMEDMEM::BlockTopology* getBlockTopology()const {return _block_topology;}
//	const string& getFilename() const {return _medfilename;}
	const int* getGlobalNumbering(MED_EN::medEntityMesh)const; 
private:
	//mesh object underlying the ParaMESH object
	MEDMEM::MESH* _mesh;
  
  bool _has_mesh_ownership;
	//name of the mesh
	string _name;
	//connect zone
	std::vector<MEDMEM::CONNECTZONE*> _connect_zone;
	//id of the local grid
	int _my_domain_id;
	//global topology of the cells
	ParaMEDMEM::BlockTopology* _block_topology;
	// pointers to global numberings
	int* _nodeglobal;
	int* _edgeglobal;
	int* _faceglobal;
	int* _cellglobal;
 
};

}

#endif /*PARAMESH_H_*/
