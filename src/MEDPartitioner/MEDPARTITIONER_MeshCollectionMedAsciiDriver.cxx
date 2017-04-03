// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#include "MEDPARTITIONER_ParallelTopology.hxx"
#include "MEDPARTITIONER_MeshCollectionDriver.hxx"
#include "MEDPARTITIONER_MeshCollection.hxx"
#include "MEDPARTITIONER_MeshCollectionMedAsciiDriver.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_Utils.hxx"

#include "MEDCouplingUMesh.hxx"
#include "MEDLoader.hxx"

#include <map>
#include <set>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

using namespace MEDPARTITIONER;

MeshCollectionMedAsciiDriver::MeshCollectionMedAsciiDriver(MeshCollection* collection):MeshCollectionDriver(collection)
{
}

/*!reads a MED File v>=2.3
 * and mounts the corresponding meshes in memory
 * the connect zones are created from the joints
 * 
 *\param filename ascii file containing the list of MED v2.3 files
 * */

int MeshCollectionMedAsciiDriver::read(MEDCoupling::MEDFileData* filedata)
{
  readMEDFileData(filedata);

  std::vector<MEDPARTITIONER::ConnectZone*> cz; // to fill from filedata
  std::vector<int*> cellglobal;
  std::vector<int*> nodeglobal;
  std::vector<int*> faceglobal;
  int size = (_collection->getMesh()).size();
  cellglobal.resize(size);
  nodeglobal.resize(size);
  faceglobal.resize(size);
  for ( int idomain = 0; idomain < size; ++idomain )
    {
      cellglobal[idomain]=0;
      faceglobal[idomain]=0;
      nodeglobal[idomain]=0;
      if ( (_collection->getMesh())[idomain] && (_collection->getMesh())[idomain]->getNumberOfNodes() > 0 )
        _collection->setNonEmptyMesh(idomain);
    }
  //creation of topology from mesh and connect zones
  ParallelTopology* aPT = new ParallelTopology((_collection->getMesh()), cz, cellglobal, nodeglobal, faceglobal);
  _collection->setTopology(aPT,true);

  return 0;
}

/*!reads a MED File v>=2.3
 * and mounts the corresponding meshes in memory
 * the connect zones are created from the joints
 *
 *\param filename ascii file containing the list of MED v2.3 files
 * */

int MeshCollectionMedAsciiDriver::read(const char* filename, ParaDomainSelector* domainSelector)
{
  //distributed meshes
  std::vector<int*> cellglobal;
  std::vector<int*> nodeglobal;
  std::vector<int*> faceglobal;
  int nbdomain;

  //reading ascii master file
  try
    {
      std::ifstream asciiinput(filename);
      if (!asciiinput)
        throw INTERP_KERNEL::Exception("Master ASCII File does not exist");
      char charbuffer[512];
      asciiinput.getline(charbuffer,512);

      while (charbuffer[0]=='#')
        {
          asciiinput.getline(charbuffer,512);
        }

      //reading number of domains
      nbdomain=atoi(charbuffer);
      MyGlobals::_File_Names.resize(nbdomain);
      MyGlobals::_Mesh_Names.resize(nbdomain);
      (_collection->getMesh()).resize(nbdomain);
      cellglobal.resize(nbdomain);
      nodeglobal.resize(nbdomain);
      faceglobal.resize(nbdomain);

      if (nbdomain == 0)
        throw INTERP_KERNEL::Exception("Empty ASCII master file");
      for (int i=0; i<nbdomain;i++)
        {
          //reading information about the domain
          std::string mesh,host;
          int idomain;
          cellglobal[i]=0;
          faceglobal[i]=0;
          nodeglobal[i]=0;

          asciiinput >> mesh >> idomain >> MyGlobals::_Mesh_Names[i] >> host >> MyGlobals::_File_Names[i];

          //Setting the name of the global mesh (which should be is the same for all the subdomains)
          if (i==0)
            _collection->setName(mesh);

          if (idomain!=i+1)
            {
              throw INTERP_KERNEL::Exception("domain must be written from 1 to N in ASCII file descriptor");
            }
          if ( !domainSelector || domainSelector->isMyDomain(i))
            readSubdomain(i);

        } //loop on domains
    } //of try
  catch(...)
    {
      throw INTERP_KERNEL::Exception("I/O error reading parallel MED file");
    }

  //creation of topology from mesh and connect zones
  ParallelTopology* aPT = new ParallelTopology((_collection->getMesh()), (_collection->getCZ()), cellglobal, nodeglobal, faceglobal);
  _collection->setTopology(aPT, true);

  for (int i=0; i<nbdomain; i++)
    {
      delete [] cellglobal[i];
      delete [] nodeglobal[i];
      delete [] faceglobal[i];
    }
  return 0;
}

/*! writes the collection of meshes in a MED v2.3 file
 * with the connect zones being written as joints
 * \param filename name of the ascii file containing the meshes description
 */
void MeshCollectionMedAsciiDriver::write(const char* filename, ParaDomainSelector* domainSelector) const
{
  int nbdomains=_collection->getMesh().size();
  std::vector<std::string> filenames;
  filenames.resize(nbdomains);

  //loop on the domains
  for (int idomain=0; idomain<nbdomains; idomain++)
    {
      std::string distfilename;
      std::ostringstream suffix;
      suffix << filename << idomain+1 << ".med";
      distfilename=suffix.str();
      filenames[idomain]=distfilename;

      if ( !domainSelector || domainSelector->isMyDomain( idomain ) )
        {
          // [ABN] spurious test in 8.2 - fixed as I think it should be:
          if ( _collection->getMesh()[idomain]->getNumberOfCells() == 0 ) continue;
          WriteUMesh(distfilename.c_str(),(_collection->getMesh())[idomain],true);
          //writeSubdomain(idomain, nbdomains, distfilename.c_str(), domainSelector);
        }
    }

  //write master file
  if ( !domainSelector || domainSelector->rank() == 0 )
    {
      std::ofstream file(filename);
      file << "#MED Fichier V 2.3"<<" " << std::endl;
      file << "#" << " " << std::endl;
      file << _collection->getMesh().size() << " " << std::endl;

      for (int idomain=0; idomain<nbdomains; idomain++)
        file << _collection->getName() <<" "<< idomain+1 << " "
             << (_collection->getMesh())[idomain]->getName() << " localhost "
             << filenames[idomain] << " "<< std::endl;
    }

}
