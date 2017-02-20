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
// Author : Anthony Geay (CEA/DEN)

#include "ParaMEDLoader.hxx"
#include "MEDLoader.hxx"
#include "ParaMESH.hxx"
#include "BlockTopology.hxx"
#include "MEDCouplingUMesh.hxx"

#include <fstream>
#include <sstream>

using namespace MEDCoupling;

ParaMEDLoader::ParaMEDLoader()
{
}

void ParaMEDLoader::WriteParaMesh(const char *fileName, MEDCoupling::ParaMESH *mesh)
{
  if(!mesh->getBlockTopology()->getProcGroup()->containsMyRank())
    return ;
  int myRank=mesh->getBlockTopology()->getProcGroup()->myRank();
  int nbDomains=mesh->getBlockTopology()->getProcGroup()->size();
  std::vector<std::string> fileNames(nbDomains);
  for(int i=0;i<nbDomains;i++)
    {
      std::ostringstream sstr;
      sstr << fileName << i+1 << ".med";
      fileNames[i]=sstr.str();
    }
  if(myRank==0)
    WriteMasterFile(fileName,fileNames,mesh->getCellMesh()->getName().c_str());
  WriteUMesh(fileNames[myRank].c_str(),dynamic_cast<MEDCouplingUMesh *>(mesh->getCellMesh()),true);
}

/*!
 * This method builds the master file 'fileName' of a parallel MED file defined in 'fileNames'.
 */
void ParaMEDLoader::WriteMasterFile(const char *fileName, const std::vector<std::string>& fileNames, const char *meshName)
{
  int nbOfDom=fileNames.size();
  std::ofstream fs(fileName);
  fs << "#MED Fichier V 2.3" << " " << std::endl;
  fs << "#"<<" " << std::endl;
  fs << nbOfDom <<" " << std::endl;
  for(int i=0;i<nbOfDom;i++)
    fs << meshName << " " << i+1 << " " << meshName << "_" << i+1 << " localhost " << fileNames[i] << " " << std::endl;
}
