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
#include "ParaGRID.hxx"

#include <fstream>
#include <vector>

using namespace std;

namespace ParaMEDMEM
{
	
ParaGRID::ParaGRID(MEDMEM::GRID* global_grid, Topology* topology)throw (MEDMEM::MEDEXCEPTION):
_name(global_grid->getName())
{
	
	_block_topology = dynamic_cast<BlockTopology*>(topology);
	if (_block_topology==0) throw MEDEXCEPTION(LOCALIZED(
	"ParaGRID::ParaGRID topology must be block topology"));
	
	if (!_block_topology->getProcGroup()->containsMyRank()) return;
	
	int dimension=_block_topology->getDimension() ;
	if (dimension != global_grid->getSpaceDimension())
		throw MEDEXCEPTION(LOCALIZED("ParaGrid::ParaGrid incompatible topology"));
	
	vector<vector<double> > xyz_array(dimension);
	vector<pair<int,int> > local_indices = _block_topology->getLocalArrayMinMax();
	int myrank=_block_topology->getProcGroup()->myRank();
	vector <string> coordinates_names;
	vector <string> coordinates_units;
	for (int idim=0; idim<dimension ; idim++)
	{
		cout << " Indices "<< local_indices[idim].first <<" "<<local_indices[idim].second<<endl;
		for (int i=(local_indices)[idim].first; i<(local_indices)[idim].second; i++)
			xyz_array[idim].push_back(global_grid->getArrayValue(idim+1,i));
		coordinates_names.push_back(global_grid->getCoordinatesNames()[idim]);
		coordinates_units.push_back(global_grid->getCoordinatesUnits()[idim]);
	}
	_grid=new MEDMEM::GRID(xyz_array,
						   coordinates_names,
  					       coordinates_units);
	_grid->setName(global_grid->getName());
	_grid->setDescription(global_grid->getDescription());
	
}

ParaGRID::ParaGRID(MEDMEM::driverTypes driver_type, const string& file_name, 
	const string& driver_name, int domain_id)
throw (MEDMEM::MEDEXCEPTION){};
ParaGRID::~ParaGRID(){if (_grid !=0) delete _grid;};


/*! method for writing a distributed grid
 * 
 * \param driverType type of driver used (MED_DRIVER,VTK_DRIVER)
 * \param master_filename name of the master file
 */
 
void ParaGRID::write(MEDMEM::driverTypes driverType, const string& master_filename)
throw (MEDMEM::MEDEXCEPTION){
	
  const char* LOC = "ParaMEDMEM::ParaGRID::write()";
  BEGIN_OF_MED(LOC);
	 
	if (!_block_topology->getProcGroup()->containsMyRank()) return;
	 
	int myrank=_block_topology->getProcGroup()->myRank();
	
	ofstream file(master_filename.c_str());
	int nbdomains= _block_topology->getProcGroup()->size();	
	vector<string> filename(nbdomains);
		
		
	//loop on the domains
	for (int i=0; i<nbdomains;i++)
	{
		char distfilename[256];
	
		ostringstream suffix;
		
		suffix << master_filename<< i+1 <<".med";
		
		strcpy(distfilename,suffix.str().c_str());
		filename[i]=string(distfilename);
		MESSAGE_MED("File name "<<string(distfilename));
	}	
	
	//creation of the master file by proc 0 on ProgGroup	
	if (myrank==0)
	{
		file <<"#MED Fichier V 2.3"<<" "<<endl;
		file <<"#"<<" "<<endl;
		file<<nbdomains<<" "<<endl;
		for (int i=0; i<nbdomains;i++)
		{
			//updating the ascii description file
			file << _name <<" "<< i+1 << " "<< _name << " localhost " << filename[i] << " "<<endl;
		}
	}
		
	int id=_grid->addDriver(MEDMEM::MED_DRIVER,filename[myrank],_name);
		
	MESSAGE_MED("Start writing");
	_grid->write(id);
			
	END_OF_MED("ParaMEDMEM::ParaGRID::write()");
};

	
	
}
