// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_Utils.hxx"

#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#include "MEDFileMesh.hxx"

#include <map>
#include <set>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

#include "med.h"

using namespace MEDPARTITIONER;

MeshCollectionDriver::MeshCollectionDriver(MeshCollection* collection):_collection(collection)
{
}

/*!reads a unique MED File v>=2.1
 * and mounts the corresponding mesh in memory
 *\param filename binary file
 *\param meshname mesh name in the MED file
 * */
int MeshCollectionDriver::readSeq(const char* filename, const char* meshname)
{
  std::cout << "readSeq" << std::endl;
  MyGlobals::_File_Names.resize(1);
  MyGlobals::_File_Names[0]=std::string(filename);

  ParaMEDMEM::MEDFileUMesh* mfm=ParaMEDMEM::MEDFileUMesh::New(filename,meshname);
  //puts the only mesh in the mesh vector
  (_collection->getMesh()).push_back(mfm->getLevel0Mesh(false));
  (_collection->getFaceMesh()).push_back(mfm->getLevelM1Mesh(false));

  //reading family ids
  ParaMEDMEM::DataArrayInt* cellIds(mfm->getFamilyFieldAtLevel(0)->deepCpy());
  ParaMEDMEM::DataArrayInt* faceIds(mfm->getFamilyFieldAtLevel(-1)->deepCpy());
  (_collection->getCellFamilyIds()).push_back(cellIds);
  (_collection->getFaceFamilyIds()).push_back(faceIds); 

  //reading groups
  (_collection->getFamilyInfo())=mfm->getFamilyInfo();
  (_collection->getGroupInfo())=mfm->getGroupInfo();

  (_collection->getCZ()).clear();
  
  ParallelTopology* aPT = new ParallelTopology((_collection->getMesh()));
  _collection->setTopology(aPT);
  _collection->setName(meshname);
  _collection->setDomainNames(meshname);
  return 0;
}


//================================================================================
/*!
 * \brief Return mesh dimension from distributed med file had being read
 */
//================================================================================

void MeshCollectionDriver::readSubdomain(std::vector<int*>& cellglobal,
                                         std::vector<int*>& faceglobal,
                                         std::vector<int*>& nodeglobal, int idomain)
{
  std::string meshname=MyGlobals::_Mesh_Names[idomain];
  std::string file=MyGlobals::_File_Names[idomain];

  ParaMEDMEM::MEDFileUMesh* mfm=ParaMEDMEM::MEDFileUMesh::New(file,meshname);
  std::vector<int> nonEmpty=mfm->getNonEmptyLevels();
  
  try 
    { 
      (_collection->getMesh())[idomain]=mfm->getLevel0Mesh(false); 
      //reading families groups
      ParaMEDMEM::DataArrayInt* cellIds(mfm->getFamilyFieldAtLevel(0)->deepCpy());
      (_collection->getCellFamilyIds())[idomain]=cellIds;
    }
  catch(...)
    { 
      (_collection->getMesh())[idomain]=CreateEmptyMEDCouplingUMesh(); // or 0 if you want tests;
      ParaMEDMEM::DataArrayInt* empty=ParaMEDMEM::DataArrayInt::New();
      empty->alloc(0,1);
      (_collection->getCellFamilyIds())[idomain]=empty;
      std::cout << "\nNO Level0Mesh (Cells)\n";
    }
  try 
    { 
      if (nonEmpty.size()>1 && nonEmpty[1]==-1)
        {
          (_collection->getFaceMesh())[idomain]=mfm->getLevelM1Mesh(false);
          //reading families groups
          ParaMEDMEM::DataArrayInt* faceIds(mfm->getFamilyFieldAtLevel(-1)->deepCpy());
          (_collection->getFaceFamilyIds())[idomain]=faceIds;
          if (MyGlobals::_Verbose>10)
            std::cout << "proc " << MyGlobals::_Rank << " : WITH Faces\n";

        }
      else
        {
          throw INTERP_KERNEL::Exception("no faces");
        }
    }
  catch(...)
    {
      (_collection->getFaceMesh())[idomain]=CreateEmptyMEDCouplingUMesh(); // or 0 if you want test;
      ParaMEDMEM::DataArrayInt* empty=ParaMEDMEM::DataArrayInt::New();
      (_collection->getFaceFamilyIds())[idomain]=empty;
      if (MyGlobals::_Verbose>10)
        std::cout << "proc " << MyGlobals::_Rank << " : WITHOUT Faces\n";
    }
  
  //reading groups
  _collection->getFamilyInfo()=mfm->getFamilyInfo();
  _collection->getGroupInfo()=mfm->getGroupInfo();

  mfm->decrRef();
  
  std::vector<std::string> localInformation;
  std::string str;
  localInformation.push_back(str+"ioldDomain="+IntToStr(idomain));
  localInformation.push_back(str+"meshName="+meshname);
  MyGlobals::_General_Informations.push_back(SerializeFromVectorOfString(localInformation));
  std::vector<std::string> localFields=BrowseAllFieldsOnMesh(file, meshname, idomain);
  if (localFields.size()>0) 
    MyGlobals::_Field_Descriptions.push_back(SerializeFromVectorOfString(localFields));
}


void MeshCollectionDriver::readSubdomain(int idomain)
{
  std::string meshname=MyGlobals::_Mesh_Names[idomain];
  std::string file=MyGlobals::_File_Names[idomain];

  ParaMEDMEM::MEDFileUMesh* mfm=ParaMEDMEM::MEDFileUMesh::New(file,meshname);
  std::vector<int> nonEmpty=mfm->getNonEmptyLevels();
  
  try 
    { 
      (_collection->getMesh())[idomain]=mfm->getLevel0Mesh(false); 
      //reading families groups
      ParaMEDMEM::DataArrayInt* cellIds(mfm->getFamilyFieldAtLevel(0)->deepCpy());
      (_collection->getCellFamilyIds())[idomain]=cellIds;
    }
  catch(...)
    { 
      (_collection->getMesh())[idomain]=CreateEmptyMEDCouplingUMesh(); // or 0 if you want tests;
      ParaMEDMEM::DataArrayInt* empty=ParaMEDMEM::DataArrayInt::New();
      empty->alloc(0,1);
      (_collection->getCellFamilyIds())[idomain]=empty;
      std::cout<<"\nNO Level0Mesh (Cells)\n";
    }
  try 
    { 
      if (nonEmpty.size()>1 && nonEmpty[1]==-1)
        {
          (_collection->getFaceMesh())[idomain]=mfm->getLevelM1Mesh(false);
          //reading families groups
          ParaMEDMEM::DataArrayInt* faceIds(mfm->getFamilyFieldAtLevel(-1)->deepCpy());
          (_collection->getFaceFamilyIds())[idomain]=faceIds;
          if (MyGlobals::_Verbose>10)
            std::cout << "proc " << MyGlobals::_Rank << " : WITH Faces\n";
        }
      else
        {
          throw INTERP_KERNEL::Exception("no faces");
        }
    }
  catch(...)
    {
      (_collection->getFaceMesh())[idomain]=CreateEmptyMEDCouplingUMesh(); // or 0 if you want test;
      ParaMEDMEM::DataArrayInt* empty=ParaMEDMEM::DataArrayInt::New();
      (_collection->getFaceFamilyIds())[idomain]=empty;
      if (MyGlobals::_Verbose>10)
        std::cout << "proc " << MyGlobals::_Rank << " : WITHOUT Faces\n";
    }
  
  //reading groups
  _collection->getFamilyInfo()=mfm->getFamilyInfo();
  _collection->getGroupInfo()=mfm->getGroupInfo();

  mfm->decrRef();
  
  std::vector<std::string> localInformation;
  std::string str;
  localInformation.push_back(str+"ioldDomain="+IntToStr(idomain));
  localInformation.push_back(str+"meshName="+meshname);
  MyGlobals::_General_Informations.push_back(SerializeFromVectorOfString(localInformation));
  std::vector<std::string> localFields=BrowseAllFieldsOnMesh(file, meshname, idomain);
  if (localFields.size()>0) 
    MyGlobals::_Field_Descriptions.push_back(SerializeFromVectorOfString(localFields));
}


void MeshCollectionDriver::writeMedFile(int idomain, const std::string& distfilename) const
{
  std::vector<const ParaMEDMEM::MEDCouplingUMesh*> meshes;
  ParaMEDMEM::MEDCouplingUMesh* cellMesh=_collection->getMesh(idomain);
  ParaMEDMEM::MEDCouplingUMesh* faceMesh=_collection->getFaceMesh(idomain);
  //ParaMEDMEM::MEDCouplingUMesh* faceMeshFilter=0;
  
  std::string finalMeshName=ExtractFromDescription(MyGlobals::_General_Informations[0], "finalMeshName=");
  // std::string cleFilter=Cle1ToStr("filterFaceOnCell",idomain);
  // ParaMEDMEM::DataArrayInt* filter=0;
  // if (_collection->getMapDataArrayInt().find(cleFilter)!=_collection->getMapDataArrayInt().end())
  //   {
  //     filter=_collection->getMapDataArrayInt().find(cleFilter)->second;
  //     int* index=filter->getPointer();
  //     faceMeshFilter=(ParaMEDMEM::MEDCouplingUMesh *) faceMesh->buildPartOfMySelf(index,index+filter->getNbOfElems(),true);
  //     faceMesh=faceMeshFilter;
  //   }
  cellMesh->setName(finalMeshName);
  meshes.push_back(cellMesh);
  
  faceMesh->checkCoherency();
  if (faceMesh->getNumberOfCells()>0)
    {
      faceMesh->tryToShareSameCoordsPermute(*cellMesh, 1e-10);
      meshes.push_back(faceMesh);
    }
  
  //ParaMEDMEM::MEDCouplingUMesh* boundaryMesh=0;
  // if (MyGlobals::_Creates_Boundary_Faces>0)
  //   {
  //     //try to write Boundary meshes
  //     bool keepCoords=false; //TODO or true
  //     boundaryMesh=(ParaMEDMEM::MEDCouplingUMesh *) cellMesh->buildBoundaryMesh(keepCoords);
  //     boundaryMesh->setName("boundaryMesh");
  //   }

  MEDLoader::WriteUMeshes(distfilename, meshes, true);
  // if (faceMeshFilter!=0)
  //   faceMeshFilter->decrRef();

  // if (boundaryMesh!=0)
  //   {
  //     //doing that testMesh becomes second mesh sorted by alphabetical order of name
  //     MEDLoader::WriteUMesh(distfilename, boundaryMesh, false);
  //     boundaryMesh->decrRef();
  //   }
  ParaMEDMEM::MEDFileUMesh* mfm=ParaMEDMEM::MEDFileUMesh::New(distfilename, _collection->getMesh(idomain)->getName());

  mfm->setFamilyInfo(_collection->getFamilyInfo());
  mfm->setGroupInfo(_collection->getGroupInfo());

  std::string key=Cle1ToStr("faceFamily_toArray",idomain);
  if ( meshes.size() == 2 &&
      _collection->getMapDataArrayInt().find(key)!=_collection->getMapDataArrayInt().end())
    {
      ParaMEDMEM::DataArrayInt *fam=_collection->getMapDataArrayInt().find(key)->second;
      mfm->setFamilyFieldArr(-1,fam);
    }

  key=Cle1ToStr("cellFamily_toArray",idomain);
  if (_collection->getMapDataArrayInt().find(key)!=_collection->getMapDataArrayInt().end())
    mfm->setFamilyFieldArr(0,_collection->getMapDataArrayInt().find(key)->second);

  mfm->write(distfilename,0);
  key="/inewFieldDouble="+IntToStr(idomain)+"/";

  std::map<std::string,ParaMEDMEM::DataArrayDouble*>::iterator it;
  int nbfFieldFound=0;
  for (it=_collection->getMapDataArrayDouble().begin() ; it!=_collection->getMapDataArrayDouble().end(); it++)
    {
      std::string desc=(*it).first;
      size_t found=desc.find(key);
      if (found==std::string::npos)
        continue;
      if (MyGlobals::_Verbose>20)
        std::cout << "proc " << MyGlobals::_Rank << " : write field " << desc << std::endl;
      std::string meshName, fieldName;
      int typeField, DT, IT, entity;
      FieldShortDescriptionToData(desc, fieldName, typeField, entity, DT, IT);
      double time=StrToDouble(ExtractFromDescription(desc, "time="));
      int typeData=StrToInt(ExtractFromDescription(desc, "typeData="));
      std::string entityName=ExtractFromDescription(desc, "entityName=");
      ParaMEDMEM::MEDCouplingFieldDouble* field=0;
      if (typeData!=6)
        {
          std::cout << "WARNING : writeMedFile : typeData " << typeData << " not implemented for fields\n";
          continue;
        }
      if (entityName=="MED_CELL")
        {
          //there is a field of idomain to write
          field=ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_CELLS,ParaMEDMEM::ONE_TIME);
        }
      if (entityName=="MED_NODE_ELEMENT")
        {
          //there is a field of idomain to write
          field=ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_GAUSS_NE,ParaMEDMEM::ONE_TIME);
        }
      if (!field)
        {
          std::cout << "WARNING : writeMedFile : entityName " << entityName << " not implemented for fields\n";
          continue;
        }
      nbfFieldFound++;
      field->setName(fieldName);
      field->setMesh(mfm->getLevel0Mesh(false));
      ParaMEDMEM::DataArrayDouble *da=(*it).second;
    
      //get information for components etc..
      std::vector<std::string> r1;
      r1=SelectTagsInVectorOfString(MyGlobals::_General_Informations,"fieldName="+fieldName);
      r1=SelectTagsInVectorOfString(r1,"typeField="+IntToStr(typeField));
      r1=SelectTagsInVectorOfString(r1,"DT="+IntToStr(DT));
      r1=SelectTagsInVectorOfString(r1,"IT="+IntToStr(IT));
      //not saved in file? field->setDescription(ExtractFromDescription(r1[0], "fieldDescription="));
      int nbc=StrToInt(ExtractFromDescription(r1[0], "nbComponents="));
      if (nbc==da->getNumberOfComponents())
        {
          for (int i=0; i<nbc; i++) 
            da->setInfoOnComponent(i,ExtractFromDescription(r1[0], "componentInfo"+IntToStr(i)+"="));
        }
      else
        {
          std::cerr << "Problem On field " << fieldName << " : number of components unexpected " << da->getNumberOfComponents() << std::endl;
        }
    
      field->setArray(da);
      field->setTime(time,DT,IT);
      field->checkCoherency();
      try
        {
          MEDLoader::WriteField(distfilename,field,false);
        }
      catch(INTERP_KERNEL::Exception& e)
        {
          //cout trying rewrite all data, only one field defined
          std::string tmp,newName=distfilename;
          tmp+="_"+fieldName+"_"+IntToStr(nbfFieldFound)+".med";
          newName.replace(newName.find(".med"),4,tmp);
          std::cout << "WARNING : writeMedFile : create a new file name with only one field because MEDLoader::WriteField throw:" << newName << std::endl;
          MEDLoader::WriteField(newName,field,true);
        }
    }
  mfm->decrRef();
}
