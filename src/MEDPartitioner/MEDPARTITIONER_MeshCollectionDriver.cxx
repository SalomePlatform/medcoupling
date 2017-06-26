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

#include "MEDPARTITIONER_MeshCollectionDriver.hxx"

#include "MEDPARTITIONER_ConnectZone.hxx"
#include "MEDPARTITIONER_MeshCollection.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_ParallelTopology.hxx"
#include "MEDPARTITIONER_Utils.hxx"

#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "MEDCouplingSkyLineArray.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDFileData.hxx"
#include "MEDFileField.hxx"
#include "MEDFileJoint.hxx"
#include "MEDFileMesh.hxx"
#include "MEDLoader.hxx"

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

  MEDCoupling::MEDFileUMesh* mfm=MEDCoupling::MEDFileUMesh::New(filename,meshname);
  //puts the only mesh in the mesh vector
  (_collection->getMesh()).push_back(mfm->getLevel0Mesh(false));
  (_collection->getFaceMesh()).push_back(mfm->getLevelM1Mesh(false));

  //reading family ids
  MEDCoupling::DataArrayInt* cellIds(mfm->getFamilyFieldAtLevel(0)->deepCopy());
  MEDCoupling::DataArrayInt* faceIds(mfm->getFamilyFieldAtLevel(-1)->deepCopy());
  (_collection->getCellFamilyIds()).push_back(cellIds);
  (_collection->getFaceFamilyIds()).push_back(faceIds); 

  //reading groups
  (_collection->getFamilyInfo())=mfm->getFamilyInfo();
  (_collection->getGroupInfo())=mfm->getGroupInfo();

  (_collection->getCZ()).clear();
  
  ParallelTopology* aPT = new ParallelTopology((_collection->getMesh()));
  _collection->setTopology(aPT, true);
  _collection->setName(meshname);
  _collection->setDomainNames(meshname);
  return 0;
}


void MeshCollectionDriver::readMEDFileData(const MEDCoupling::MEDFileData* filedata)
{
  const int nbDomains = filedata->getMeshes()->getNumberOfMeshes();
  _collection->getMesh()         .resize( nbDomains, 0 );
  _collection->getFaceMesh()     .resize( nbDomains, 0 );
  _collection->getCellFamilyIds().resize( nbDomains, 0 );
  _collection->getFaceFamilyIds().resize( nbDomains, 0 );

  for (int i=0; i<nbDomains; i++)
    {
      MEDCoupling::MEDFileUMesh *mfm = dynamic_cast<MEDCoupling::MEDFileUMesh *>(filedata->getMeshes()->getMeshAtPos(i));
      readData(mfm,i);
      if ( mfm && mfm->getMeshDimension() > 0 )
        _collection->setNonEmptyMesh( i );
    }

  ParallelTopology* aPT = new ParallelTopology(_collection->getMesh());
  _collection->setTopology(aPT, true);
  if ( nbDomains > 0 )
    {
      _collection->setName( filedata->getMeshes()->getMeshAtPos(0)->getName() );
      _collection->setDomainNames( _collection->getName() );
    }
  if ( ParaDomainSelector* domainSelector = _collection->getParaDomainSelector() )
    if ( _collection->isParallelMode() )
      {
        //to know nb of cells on each proc to compute global cell ids from locally global
        domainSelector->gatherNbOf(_collection->getMesh());
      }
}

void MeshCollectionDriver::readFileData(std::string file,std::string meshname,int idomain) const
{
  MEDCoupling::MEDFileUMesh* mfm=MEDCoupling::MEDFileUMesh::New(file,meshname);
  readData(mfm,idomain);
  mfm->decrRef();
}

void MeshCollectionDriver::readData(MEDCoupling::MEDFileUMesh* mfm, int idomain) const
{
  std::vector<int> nonEmpty=mfm->getNonEmptyLevels();
  try
    {
      (_collection->getMesh())[idomain]=mfm->getLevel0Mesh(false);
      //reading families groups
      MEDCoupling::DataArrayInt* cellIds(mfm->getFamilyFieldAtLevel(0)->deepCopy());
      (_collection->getCellFamilyIds())[idomain]=cellIds;
    }
  catch(...)
    {
      (_collection->getMesh())[idomain]=CreateEmptyMEDCouplingUMesh(); // or 0 if you want tests;
      MEDCoupling::DataArrayInt* empty=MEDCoupling::DataArrayInt::New();
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
          MEDCoupling::DataArrayInt* faceIds(mfm->getFamilyFieldAtLevel(-1)->deepCopy());
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
      MEDCoupling::DataArrayInt* empty=MEDCoupling::DataArrayInt::New();
      (_collection->getFaceFamilyIds())[idomain]=empty;
      if (MyGlobals::_Verbose>10)
        std::cout << "proc " << MyGlobals::_Rank << " : WITHOUT Faces\n";
    }
  //reading groups
  _collection->getFamilyInfo()=mfm->getFamilyInfo();
  _collection->getGroupInfo()=mfm->getGroupInfo();
}

void MeshCollectionDriver::readSubdomain(int idomain)
{
  std::string meshname=MyGlobals::_Mesh_Names[idomain];
  std::string file=MyGlobals::_File_Names[idomain];
  readFileData(file,meshname,idomain);
  
  std::vector<std::string> localInformation;
  std::string str;
  localInformation.push_back(str+"ioldDomain="+IntToStr(idomain));
  localInformation.push_back(str+"meshName="+meshname);
  MyGlobals::_General_Informations.push_back(SerializeFromVectorOfString(localInformation));
  std::vector<std::string> localFields=BrowseAllFieldsOnMesh(file, meshname, idomain);
  if (localFields.size()>0)
    MyGlobals::_Field_Descriptions.push_back(SerializeFromVectorOfString(localFields));
}

MEDCoupling::MEDFileMesh* MeshCollectionDriver::getMesh(int idomain) const
{
  MEDCoupling::MEDFileUMesh* mfm = MEDCoupling::MEDFileUMesh::New();

  MEDCoupling::MEDCouplingUMesh* cellMesh=_collection->getMesh(idomain);
  MEDCoupling::MEDCouplingUMesh* faceMesh=_collection->getFaceMesh(idomain);
  // std::string cleFilter=Cle1ToStr("filterFaceOnCell",idomain);
  // MEDCoupling::DataArrayInt* filter=0;
  // if (_collection->getMapDataArrayInt().find(cleFilter)!=_collection->getMapDataArrayInt().end())
  //   {
  //     filter=_collection->getMapDataArrayInt().find(cleFilter)->second;
  //     int* index=filter->getPointer();
  //     faceMeshFilter=(MEDCoupling::MEDCouplingUMesh *) faceMesh->buildPartOfMySelf(index,index+filter->getNbOfElems(),true);
  //     faceMesh=faceMeshFilter;
  //   }
  // if (faceMeshFilter!=0)
  //   faceMeshFilter->decrRef();
  std::string finalMeshName="";
  if (MyGlobals::_General_Informations.size()!=0)
    {
      std::size_t found=MyGlobals::_General_Informations[0].find("finalMeshName=");
      if ((found!=std::string::npos) && (found>0))
        {
          finalMeshName=ExtractFromDescription(MyGlobals::_General_Informations[0], "finalMeshName=");
        }
    }
  if (finalMeshName.empty())
    {
      finalMeshName=_collection->getName();
    }
  cellMesh->setName(finalMeshName);
  mfm->setMeshAtLevel( 0, cellMesh );

  faceMesh->checkConsistencyLight();
  if (faceMesh->getNumberOfCells()>0)
    {
      faceMesh->tryToShareSameCoordsPermute(*cellMesh, 1e-10);
      faceMesh->setName(finalMeshName);
      mfm->setMeshAtLevel( -1, faceMesh );
    }

  // MEDCoupling::MEDCouplingUMesh* boundaryMesh=0;
  // if (MyGlobals::_Create_Boundary_Faces>0)
  //   {
  //     //try to write Boundary meshes
  //     bool keepCoords=false; //TODO or true
  //     boundaryMesh=(MEDCoupling::MEDCouplingUMesh *) cellMesh->buildBoundaryMesh(keepCoords);
  //     boundaryMesh->setName("boundaryMesh");
  // if (boundaryMesh!=0)
  //   {
  //     //doing that testMesh becomes second mesh sorted by alphabetical order of name
  //     WriteUMesh(distfilename, boundaryMesh, false);
  //     boundaryMesh->decrRef();
  //   }

  mfm->setFamilyInfo(_collection->getFamilyInfo());
  mfm->setGroupInfo(_collection->getGroupInfo());
  std::string key=Cle1ToStr("faceFamily_toArray",idomain);
  if ( faceMesh->getNumberOfCells()>0 && _collection->getMapDataArrayInt().find(key)!=_collection->getMapDataArrayInt().end())
    mfm->setFamilyFieldArr(-1,_collection->getMapDataArrayInt().find(key)->second);
  key=Cle1ToStr("cellFamily_toArray",idomain);
  if (_collection->getMapDataArrayInt().find(key)!=_collection->getMapDataArrayInt().end())
    mfm->setFamilyFieldArr(0,_collection->getMapDataArrayInt().find(key)->second);

  // add joints

  using MEDCoupling::MCAuto;
  using MEDCoupling::MEDCouplingSkyLineArray;
  using MEDCoupling::MEDFileJoint;
  using MEDCoupling::MEDFileJointCorrespondence;
  using MEDCoupling::MEDFileJointOneStep;
  using MEDCoupling::MEDFileJoints;
  using MEDCoupling::MEDFileJoints;

  if ( _collection->getCZ().size() > 0 )
    {
      MCAuto< MEDFileJoints > joints = MEDFileJoints::New();

      for ( size_t i = 0; i < _collection->getCZ().size(); ++i )
        {
          ConnectZone* cz = _collection->getCZ()[i];
          if ( !cz ||
               cz->getLocalDomainNumber() != idomain )
            continue;
          {
            std::ostringstream oss;
            oss << "joint_" << cz->getDistantDomainNumber();
            cz->setName( oss.str() );
          }
          {
            std::ostringstream oss;
            oss << "connect_zone_" << i;
            cz->setDescription( oss.str() );
          }

          MCAuto< MEDFileJoint>
            joint = MEDFileJoint::New( cz->getName(), finalMeshName,
                                       finalMeshName, cz->getDistantDomainNumber() );
          joint->setDescription( cz->getDescription() );
          joints->pushJoint( joint );

          MCAuto< MEDFileJointOneStep> j1st = MEDFileJointOneStep::New();
          joint->pushStep( j1st );

          const MEDCouplingSkyLineArray * nodeCorr = cz->getNodeCorresp();
          if ( nodeCorr )
            {
              MCAuto< MEDFileJointCorrespondence >
                corr = MEDFileJointCorrespondence::New( nodeCorr->getValuesArray() );
              j1st->pushCorrespondence( corr );
            }

          std::vector< std::pair< int,int > > types = cz->getEntities();
          INTERP_KERNEL::NormalizedCellType t1, t2;
          for ( size_t it = 0; it < types.size(); ++it )
            {
              const MEDCouplingSkyLineArray * cellCorr =
                cz->getEntityCorresp( types[it].first, types[it].second );
              if ( cellCorr && cellCorr->getNumberOf() > 0 )
                {
                  t1 = INTERP_KERNEL::NormalizedCellType( types[it].first );
                  t2 = INTERP_KERNEL::NormalizedCellType( types[it].second );
                  MCAuto< MEDFileJointCorrespondence>
                    corr = MEDFileJointCorrespondence::New( cellCorr->getValuesArray(), t1, t2 );
                  j1st->pushCorrespondence( corr );
                }
            }
        }
      mfm->setJoints( joints );
    }

  return mfm;
}

MEDCoupling::MEDCouplingFieldDouble* MeshCollectionDriver::getField(std::string key, std::string description, MEDCoupling::DataArrayDouble* data, MEDCoupling::MEDFileMesh* mfm, int idomain) const
{
  std::string desc=description;
  if (MyGlobals::_Verbose>20)
    std::cout << "proc " << MyGlobals::_Rank << " : write field " << desc << std::endl;
  std::string meshName, fieldName;
  int typeField, DT, IT, entity;
  FieldShortDescriptionToData(desc, fieldName, typeField, entity, DT, IT);
  double time=StrToDouble(ExtractFromDescription(desc, "time="));
  int typeData=StrToInt(ExtractFromDescription(desc, "typeData="));
  std::string entityName=ExtractFromDescription(desc, "entityName=");
  MEDCoupling::MEDCouplingFieldDouble* field=0;
  if (typeData!=6)
    {
      std::cout << "WARNING : writeMedFile : typeData " << typeData << " not implemented for fields\n";
    }
  if (entityName=="MED_CELL")
    {
      //there is a field of idomain to write
      field=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::ONE_TIME);
    }
  if (entityName=="MED_NODE_ELEMENT")
    {
      //there is a field of idomain to write
      field=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_GAUSS_NE,MEDCoupling::ONE_TIME);
    }
  if (!field)
    {
      std::cout << "WARNING : writeMedFile : entityName " << entityName << " not implemented for fields\n";
    }
  if (field && typeData==6)
    {
      field->setName(fieldName);
      field->setMesh(mfm->getMeshAtLevel(0));
      MEDCoupling::DataArrayDouble *da=data;
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
      field->checkConsistencyLight();
    }
  return field;
}

void MeshCollectionDriver::writeMedFile(int idomain, const std::string& distfilename) const
{
  MEDCoupling::MEDFileMesh* mfm = getMesh( idomain );
  mfm->write(distfilename,2);

  std::string key="/inewFieldDouble="+IntToStr(idomain)+"/";
  std::map<std::string,MEDCoupling::DataArrayDouble*>::iterator it;
  int nbfFieldFound=0;
  for (it=_collection->getMapDataArrayDouble().begin() ; it!=_collection->getMapDataArrayDouble().end(); it++)
    {
      size_t found=(*it).first.find(key);
      if (found==std::string::npos)
        continue;
      MEDCoupling::MEDCouplingFieldDouble* field=0;
      field = getField(key, (*it).first, (*it).second, mfm, idomain);
      nbfFieldFound++;
      try
        {
          WriteField(distfilename,field,false);
        }
      catch(INTERP_KERNEL::Exception& e)
        {
          //cout trying rewrite all data, only one field defined
          std::string tmp,newName=distfilename;
          std::string fieldName;
          fieldName=field->getName();
          tmp+="_"+fieldName+"_"+IntToStr(nbfFieldFound)+".med";
          newName.replace(newName.find(".med"),4,tmp);
          std::cout << "WARNING : writeMedFile : create a new file name with only one field because WriteField throw:" << newName << std::endl;
          WriteField(newName,field,true);
        }
    }
  mfm->decrRef();
}

MEDCoupling::MEDFileData* MeshCollectionDriver::getMEDFileData()
{
  MEDCoupling::MEDFileData* newdata = MEDCoupling::MEDFileData::New();

  MEDCoupling::MCAuto<MEDCoupling::MEDFileMeshes> meshes;
  MEDCoupling::MCAuto<MEDCoupling::MEDFileFields> fields;
  meshes = MEDCoupling::MEDFileMeshes::New();
  fields = MEDCoupling::MEDFileFields::New();

  for (size_t i=0; i<_collection->getMesh().size(); i++)
    {
      MEDCoupling::MEDFileMesh* mfm = getMesh( i );
      meshes->pushMesh(mfm);

      std::string key="/inewFieldDouble="+IntToStr(i)+"/";
      std::map<std::string,MEDCoupling::DataArrayDouble*>::iterator it;
      MEDCoupling::MEDFileFieldMultiTS* fieldsMTS = MEDCoupling::MEDFileFieldMultiTS::New();
      for (it=_collection->getMapDataArrayDouble().begin() ; it!=_collection->getMapDataArrayDouble().end(); it++)
        {
          size_t found=(*it).first.find(key);
          if (found==std::string::npos)
            continue;
          MEDCoupling::MEDCouplingFieldDouble* field=0;
          field=getField(key, (*it).first, (*it).second, mfm, i);
          MEDCoupling::MEDFileField1TS* f1ts = MEDCoupling::MEDFileField1TS::New();
          f1ts->setFieldNoProfileSBT(field);
          fieldsMTS->pushBackTimeStep(f1ts);

          field->decrRef();
          f1ts->decrRef();
        }
      fields->pushField(fieldsMTS);

      fieldsMTS->decrRef();
      mfm->decrRef();
    }
  newdata->setMeshes(meshes);
  newdata->setFields(fields);
  return newdata;
}
