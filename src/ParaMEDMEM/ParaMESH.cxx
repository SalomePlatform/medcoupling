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
#include "ParaMESH.hxx"

#include <fstream>
#include <vector>

#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "Topology.hxx"
#include "BlockTopology.hxx"
#include "MEDMEM_ConnectZone.hxx"

//inclusion for the namespaces
#include "MEDMEM_Field.hxx"
using namespace std;

namespace ParaMEDMEM
{
	

ParaMESH::ParaMESH(MEDMEM::driverTypes driver_type, const string& filename, 
	const ProcessorGroup& group)
throw (MEDMEM::MEDEXCEPTION) :_has_mesh_ownership(true) {

  const char* LOC = "MEDSPLITTER::MESHCollectionDriver::read()";
  BEGIN_OF_MED(LOC);
  
 
  string meshstring;
  char file[256];
  char meshname[MED_TAILLE_NOM];
  int domain_id=group.myRank();
  
  // reading ascii master file
  try{
  	MESSAGE_MED("Start reading");
    ifstream asciiinput(filename.c_str());
    char charbuffer[512];
    asciiinput.getline(charbuffer,512);
    
    while (charbuffer[0]=='#')
      {
	asciiinput.getline(charbuffer,512);
     }

    //reading number of domains
    int nbdomain=atoi(charbuffer);
    cout << "nb domain"<<nbdomain<<endl;
    //    asciiinput>>nbdomain;
   
    string mesh;
	int idomain;
	string host;
	string medfilename;
  
    for (int i=0; i<=domain_id;i++)
      {
	//reading information about the domain
		
      
		asciiinput >> mesh >> idomain >> meshstring >> host >> medfilename;
		
		if (idomain!=i+1)
		  {
		    cerr<<"Error : domain must be written from 1 to N in asciifile descriptor"<<endl;
		 	throw (MEDEXCEPTION("Error : domain must be written from 1 to N in asciifile descriptor"));
		  }
		strcpy(meshname,meshstring.c_str());
		strcpy(file,medfilename.c_str());
      }
    _name=meshstring;
	///////////////////////////////////////////
	// treatment of the domain that corresponds
	// to the local id
	///////////////////////////////////////////
	_mesh=new MEDMEM::MESH(driver_type,file, meshname);
		
	//reading MEDSPLITTER::CONNECTZONEs NODE/NODE and CELL/CELL
	med_2_3::med_idt fid = med_2_3::MEDouvrir(file,med_2_3::MED_LECTURE);
 	med_2_3::med_int njoint = med_2_3::MEDnJoint(fid, meshname);
 	for (int ijoint=1; ijoint<=njoint; ijoint++)
 	{
		int distant;
		char joint_description[MED_TAILLE_DESC];
	    char name[MED_TAILLE_NOM];
	    char name_distant[MED_TAILLE_NOM];
	    // cout << "arguments"<< fid<<" "<<file<<" "<<ijoint<<" "<<name<<" "<<joint_description<<" "<<distant<<" "<<name_distant<<endl;
	    int ncorr = med_2_3::MEDjointInfo(fid,meshname, ijoint, name, 
	    	joint_description,
		       &distant, name_distant);
		cout << "Found " << ncorr <<"correspondances in joint "<<ijoint<<endl;

    	for (int ic=1; ic<=ncorr; ic++)
	    {
	      med_2_3::med_entite_maillage cor_typent_local;
	      med_2_3::med_geometrie_element cor_typgeo_local;
	      med_2_3::med_entite_maillage cor_typent_dist;
	      med_2_3::med_geometrie_element cor_typgeo_dist;
	  
	   
	      int ncouples;
	      ncouples = med_2_3::MEDjointTypeCorres(fid, meshname, name, ic,
					    &cor_typent_local,  &cor_typgeo_local,
					    &cor_typent_dist, &cor_typgeo_dist
					    );
	      int* node_corresp=new int[ncouples];
	      if (cor_typent_local == med_2_3::MED_NOEUD && cor_typent_dist == med_2_3::MED_NOEUD)
			{
	  
			  med_2_3::MEDjointLire(fid, meshname, name,
		       node_corresp,ncouples,
		       cor_typent_local,  cor_typgeo_local,
		       cor_typent_dist, cor_typgeo_dist
		       );
			}
			//constructing the connect zone and adding it to the connect zone list
			MEDMEM::CONNECTZONE* cz = new MEDMEM::CONNECTZONE();
			cz->setName(string(name));
			cz->setDescription(joint_description);
			cz->setLocalDomainNumber(domain_id);
			cz->setDistantDomainNumber(distant);
			//cz->setLocalMesh((m_collection->getMesh())[i]);
			//cz->setDistantMesh((m_collection->getMesh())[distant]);
			cz->setNodeCorresp(node_corresp,ncouples);
			_connect_zone.push_back(cz);
			
	    }//loop on correspondances
	}//loop on joints	
	
	//
	// Reading global numbering
	// 
	int ncell=_mesh->getNumberOfElements(MED_EN::MED_CELL,MED_EN::MED_ALL_ELEMENTS);
	int * array=new int[ncell];
	int offset=0;
	MESSAGE_MED("Reading cell global numbering for mesh "<< domain_id);
	MED_EN::MESH_ENTITIES::const_iterator currentEntity;
	list<MED_EN::medGeometryElement>::const_iterator iter;
	currentEntity  = MED_EN::meshEntities.find(MED_EN::MED_CELL);
	char meshchar[MED_TAILLE_NOM];
	strcpy(meshchar,_mesh->getName().c_str());
	for (iter = (*currentEntity).second.begin();iter != (*currentEntity).second.end(); iter++)
	{
		MED_EN::medGeometryElement type=*iter;
		if (type/100 != _mesh->getMeshDimension()) continue;
		int ntype = _mesh->getNumberOfElements(MED_EN::MED_CELL,type);
		if (ntype==0) continue;
		med_2_3::MEDglobalNumLire(fid,meshname, array+offset, ntype,
	        med_2_3::MED_MAILLE, (med_2_3::med_geometrie_element)type);
        offset+=ntype;
	}
	_cellglobal=array;
	
	MESSAGE_MED("Reading node global numbering");
	int nnode= _mesh->getNumberOfNodes();
	array=new int[nnode];
	med_2_3::MEDglobalNumLire(fid,meshname, array, nnode,
		        med_2_3::MED_NOEUD, med_2_3::MED_POINT1); 
	_nodeglobal=array;
	
	MESSAGE_MED("Reading face global numbering for mesh "<<domain_id);
	int nbface=_mesh->getNumberOfElements(MED_EN::MED_FACE,MED_EN::MED_ALL_ELEMENTS);
	array=new int[nbface];
	currentEntity  = MED_EN::meshEntities.find(MED_EN::MED_FACE);
	offset=0;
	for (iter = (*currentEntity).second.begin();iter != (*currentEntity).second.end(); iter++)
	{
		MED_EN::medGeometryElement type=*iter;
		if (type/100 != _mesh->getMeshDimension()-1) continue;
		int ntype = _mesh->getNumberOfElements(MED_EN::MED_FACE,type);
		if (ntype==0) continue;
		med_2_3::MEDglobalNumLire(fid,meshname, array+offset, ntype,
	        med_2_3::MED_FACE, (med_2_3::med_geometrie_element)type);
        offset+=ntype;
	}
	_faceglobal=array;
//		faceglobal[i]=0;
    med_2_3::MEDfermer(fid);
	
    _block_topology=new BlockTopology(group,ncell); 
  
     MESSAGE_MED("end of read");
      
  }//of try
  catch(...)
    {
      cerr << "I/O error reading parallel MED file"<<endl;
      throw;
    }
   
  
  //creation of topology from mesh and connect zones
  //m_collection->setTopology(
  //	new ParallelTopology((m_collection->getMesh()),(m_collection->getCZ()),cellglobal,nodeglobal,faceglobal)
   // );
    
    END_OF_MED("MEDSPLITTER::MESHCollectionDriver::read()")
};

/*! Constructor for creating a ParaMESH from a local mesh and
 * a processor group. Constructor must be called by all the processors
 * in the group. */

ParaMESH::ParaMESH(MEDMEM::MESH& subdomain_mesh, const ProcessorGroup& proc_group, const string& name):
_mesh(&subdomain_mesh),
_my_domain_id(proc_group.myRank()),
_block_topology (new BlockTopology(proc_group, subdomain_mesh.getNumberOfElements(MED_EN::MED_CELL,MED_EN::MED_ALL_ELEMENTS))),
_has_mesh_ownership(false)
{
	ostringstream stream;
	stream<<name<<"_"<<_my_domain_id+1;
	_name=stream.str();

  _cellglobal = new int[subdomain_mesh.getNumberOfElements(MED_EN::MED_CELL, MED_EN::MED_ALL_ELEMENTS)];
  int offset = _block_topology->localToGlobal(make_pair(_my_domain_id,0));
  for (int i=0; i<subdomain_mesh.getNumberOfElements(MED_EN::MED_CELL, MED_EN::MED_ALL_ELEMENTS); i++)
    {
      _cellglobal[i]=offset+i;
    }
}



ParaMESH::~ParaMESH()
{
  if (_has_mesh_ownership) delete _mesh;
  delete _block_topology;
  delete[] _cellglobal;
}

/*! method for writing a distributed MESH
 * 
 * \param driverType type of driver used (MED_DRIVER,VTK_DRIVER)
 * \param master_filename name of the master file
 */
 
void ParaMESH::write(MEDMEM::driverTypes driverType, const string& master_filename)
throw (MEDMEM::MEDEXCEPTION){
	
  const char* LOC = "ParaMEDMEM::ParaMESH::write()";
  BEGIN_OF_MED(LOC);
	 
	if (!_block_topology->getProcGroup()->containsMyRank()) return;
	 
	int myrank=_block_topology->getProcGroup()->myRank();
	cout << "Myrank in write " << myrank<<endl;
	
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
		MESSAGE_MED("Master File Name "<<master_filename);
		ofstream file(master_filename.c_str());
		if (!file) throw (MEDEXCEPTION("Unable to create master file"));
		file <<"#MED Fichier V 2.3"<<" "<<endl;
		file <<"#"<<" "<<endl;
		file<<nbdomains<<" "<<endl;
		for (int i=0; i<nbdomains;i++)
		{
			//updating the ascii description file
		  file << _name <<" "<< i+1 << " "<< _name <<"_"<<i+1<< " localhost " << filename[i] << " "<<endl;
		}
		
	}
		
	int id=_mesh->addDriver(MEDMEM::MED_DRIVER,filename[myrank],_name);
		
	MESSAGE_MED("Start writing");
	_mesh->write(id);
	_mesh->rmDriver(id);		
	END_OF_MED("ParaMEDMEM::ParaMESH::write()");
};

const int* ParaMESH::getGlobalNumbering(const MED_EN::medEntityMesh entity)const
{
	switch (entity)
	{
		case MED_CELL:
		 	return _cellglobal;
		case MED_FACE :
			return _faceglobal;
		case MED_EDGE :
			return _edgeglobal;
		case MED_NODE:
			return _nodeglobal;
    default :
      return 0;
	}
}	
	
}
