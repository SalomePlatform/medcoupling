// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#include <vector>
#include <string>
#include <map>
#include <set>
#include <iostream>
#include <fstream>

#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

#include <sys/time.h>

extern "C" {
#include "med.h"
}

using namespace std;
using namespace ParaMEDMEM;
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
  cout<<"readSeq"<<endl;
  MyGlobals::_File_Names.resize(1);
  MyGlobals::_File_Names[0]=string(filename);

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

void MeshCollectionDriver::readSubdomain(vector<int*>& cellglobal, //cvwat03
                                         vector<int*>& faceglobal,
                                         vector<int*>& nodeglobal, int idomain)
{
  string meshname=MyGlobals::_Mesh_Names[idomain];
  string file=MyGlobals::_File_Names[idomain];
  //cout << "Reading "<<meshname<<" in "<<file<<endl; //cvw

  ParaMEDMEM::MEDFileUMesh* mfm=ParaMEDMEM::MEDFileUMesh::New(file.c_str(),meshname.c_str());
  vector<int> nonEmpty=mfm->getNonEmptyLevels();
  
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
      cout<<"\nNO Level0Mesh (Cells)\n";
    }
  try 
    { 
      if (nonEmpty.size()>1 && nonEmpty[1]==-1)
        {
          (_collection->getFaceMesh())[idomain]=mfm->getLevelM1Mesh(false);
          //reading families groups
          ParaMEDMEM::DataArrayInt* faceIds(mfm->getFamilyFieldAtLevel(-1)->deepCpy());
          (_collection->getFaceFamilyIds())[idomain]=faceIds;
        }
      else
        {
          throw "no faces";
        }
    }
  catch(...)
    {
      (_collection->getFaceMesh())[idomain]=CreateEmptyMEDCouplingUMesh(); // or 0 if you want test;
      ParaMEDMEM::DataArrayInt* empty=ParaMEDMEM::DataArrayInt::New();
      (_collection->getFaceFamilyIds())[idomain]=empty;
      if (MyGlobals::_Verbose>10) cout<<"proc "<<MyGlobals::_Rank<<" : NO LevelM1Mesh (Faces)\n";
    }
  
  //reading groups
  _collection->getFamilyInfo()=mfm->getFamilyInfo();
  _collection->getGroupInfo()=mfm->getGroupInfo();

  mfm->decrRef();
  
  vector<string> localInformation;
  string str;
  localInformation.push_back(str+"ioldDomain="+IntToStr(idomain));
  localInformation.push_back(str+"meshName="+meshname);
  MyGlobals::_General_Informations.push_back(SerializeFromVectorOfString(localInformation));
  vector<string> localFields=BrowseAllFieldsOnMesh(file, meshname, idomain); //cvwat07
  if (localFields.size()>0) 
    MyGlobals::_Field_Descriptions.push_back(SerializeFromVectorOfString(localFields));
}


void MeshCollectionDriver::readSubdomain(int idomain)
{
  string meshname=MyGlobals::_Mesh_Names[idomain];
  string file=MyGlobals::_File_Names[idomain];
  //cout << "Reading "<<meshname<<" in "<<file<<endl; //cvw

  ParaMEDMEM::MEDFileUMesh* mfm=ParaMEDMEM::MEDFileUMesh::New(file.c_str(),meshname.c_str());
  vector<int> nonEmpty=mfm->getNonEmptyLevels();
  
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
      cout<<"\nNO Level0Mesh (Cells)\n";
    }
  try 
    { 
      if (nonEmpty.size()>1 && nonEmpty[1]==-1)
        {
          (_collection->getFaceMesh())[idomain]=mfm->getLevelM1Mesh(false);
          //reading families groups
          ParaMEDMEM::DataArrayInt* faceIds(mfm->getFamilyFieldAtLevel(-1)->deepCpy());
          (_collection->getFaceFamilyIds())[idomain]=faceIds;
        }
      else
        {
          throw "no faces";
        }
    }
  catch(...)
    {
      (_collection->getFaceMesh())[idomain]=CreateEmptyMEDCouplingUMesh(); // or 0 if you want test;
      ParaMEDMEM::DataArrayInt* empty=ParaMEDMEM::DataArrayInt::New();
      (_collection->getFaceFamilyIds())[idomain]=empty;
      if (MyGlobals::_Verbose>10) cout<<"proc "<<MyGlobals::_Rank<<" : NO LevelM1Mesh (Faces)\n";
    }
  
  //reading groups
  _collection->getFamilyInfo()=mfm->getFamilyInfo();
  _collection->getGroupInfo()=mfm->getGroupInfo();

  mfm->decrRef();
  
  vector<string> localInformation;
  string str;
  localInformation.push_back(str+"ioldDomain="+IntToStr(idomain));
  localInformation.push_back(str+"meshName="+meshname);
  MyGlobals::_General_Informations.push_back(SerializeFromVectorOfString(localInformation));
  vector<string> localFields=BrowseAllFieldsOnMesh(file, meshname, idomain); //cvwat07
  if (localFields.size()>0) 
    MyGlobals::_Field_Descriptions.push_back(SerializeFromVectorOfString(localFields));
}


void MeshCollectionDriver::writeMedFile(int idomain, const string& distfilename)
{
  vector<const ParaMEDMEM::MEDCouplingUMesh*> meshes;
  ParaMEDMEM::MEDCouplingUMesh* cellMesh=_collection->getMesh(idomain);
  ParaMEDMEM::MEDCouplingUMesh* faceMesh=_collection->getFaceMesh(idomain);
  ParaMEDMEM::MEDCouplingUMesh* faceMeshFilter=0;
  
  string finalMeshName=ExtractFromDescription(MyGlobals::_General_Informations[0], "finalMeshName=");
  string cleFilter=Cle1ToStr("filterFaceOnCell",idomain);
  DataArrayInt* filter=0;
  if (_collection->getMapDataArrayInt().find(cleFilter)!=_collection->getMapDataArrayInt().end())
    {
      filter=_collection->getMapDataArrayInt().find(cleFilter)->second;
      int* index=filter->getPointer();
      faceMeshFilter=(MEDCouplingUMesh *) faceMesh->buildPartOfMySelf(index,index+filter->getNbOfElems(),true);
      faceMesh=faceMeshFilter;
    }
  cellMesh->setName(finalMeshName.c_str());
  meshes.push_back(cellMesh);
  //cellMesh->zipCoords();
  //faceMesh->zipCoords();
  
  faceMesh->checkCoherency();
  if (faceMesh->getNumberOfCells()>0)
    {
      faceMesh->tryToShareSameCoordsPermute(*cellMesh, 1e-10);
      meshes.push_back(faceMesh);
    }
  
  /*do not work
    ParaMEDMEM::MEDFileUMesh* mfm2=ParaMEDMEM::MEDFileUMesh::New();
    MEDFileUMesh* mfm2 = static_cast<MEDFileUMesh*>(cellMesh->getMeshes()->getMeshAtPos(0));
    MEDFileUMesh* mfm2 = ParaMEDMEM::MEDFileUMesh::New(cellMesh);
    string fname="FUM_"+distfilename;
    mfm2->setMeshAtLevel(0, cellMesh );
    mfm2->setMeshAtLevel(-1, faceMesh );
    mfm2->write(fname.c_str(),0);
    mfm2->decrRef();
  */
  
  ParaMEDMEM::MEDCouplingUMesh* boundaryMesh=0;
  //ParaMEDMEM::MEDCouplingUMesh* boundaryMesh1=0;
  //ParaMEDMEM::MEDCouplingUMesh* finalboundaryMesh=0;
  if (MyGlobals::_Creates_Boundary_Faces>0)
    {
      //try to write Boundary meshes
      bool keepCoords=false; //TODO or true
      boundaryMesh=(MEDCouplingUMesh *) cellMesh->buildBoundaryMesh(keepCoords);
      boundaryMesh->setName("boundaryMesh");
      //cout<<"boundaryMesh "<<boundaryMesh->getNumberOfCells()<<endl;
      //do not work if faceMesh present yet //The mesh dimension of meshes must be different each other!
      //boundaryMesh->checkCoherency();
      //boundaryMesh->tryToShareSameCoordsPermute(*cellMesh, 1e-10);
      //meshes.push_back(boundaryMesh);
      //string boundary="boundary_"+distfilename;
    
      /*try to find joint do no work
        int rang=MyGlobals::_Rank;
        if (rang==1) (_collection->getParaDomainSelector())->sendMesh(*(boundaryMesh),0);
        if (rang==0) 
        {
        (_collection->getParaDomainSelector())->recvMesh(boundaryMesh1,1);
        //vector<const ParaMEDMEM::MEDCouplingUMesh*> meshes;
        //vector<DataArrayInt* > corr;
        //meshes.push_back(boundaryMesh);
        //meshes.push_back(boundaryMesh1);
        //need share the same coords
        //boundaryMesh1->tryToShareSameCoordsPermute(*boundaryMesh, 1e-10);
        //finalboundaryMesh=MEDCouplingUMesh::FuseUMeshesOnSameCoords(meshes,2, corr);
        //boundaryMesh=finalboundaryMesh;
      
        boundaryMesh->zipCoords();
        boundaryMesh1->zipCoords();
        finalboundaryMesh=MEDCouplingUMesh::MergeUMeshes(boundaryMesh,boundaryMesh1);
        DataArrayInt* commonNodes=0;
        commonNodes=finalboundaryMesh->zipCoordsTraducer();
        boundaryMesh=finalboundaryMesh;
        cout<<"zipcoords"<<commonNodes->repr()<<endl;
        }
      */
    }
  
  MEDLoader::WriteUMeshes(distfilename.c_str(), meshes, true);
  if (faceMeshFilter!=0) faceMeshFilter->decrRef();
  
  
  if (boundaryMesh!=0)
    {
      //doing that testMesh becomes second mesh sorted by alphabetical order of name
      MEDLoader::WriteUMesh(distfilename.c_str(), boundaryMesh, false);
      boundaryMesh->decrRef();
    }

  //cout<<"familyInfo :\n"<<ReprMapOfStringInt(_collection->getFamilyInfo())<<endl;
  //cout<<"groupInfo :\n"<<ReprMapOfStringVectorOfString(_collection->getGroupInfo())<<endl;
  ParaMEDMEM::MEDFileUMesh* mfm=ParaMEDMEM::MEDFileUMesh::New(distfilename.c_str(), _collection->getMesh(idomain)->getName());
        
  /*example of adding new family
    (_collection->getFamilyInfo())["FaceNotOnCell"]=-500;
    vector<string> FaceNotOnCell;
    FaceNotOnCell.push_back("FaceNotOnCell");
    (_collection->getGroupInfo())["FaceNotOnCell"]=FaceNotOnCell;
  */
  
  mfm->setFamilyInfo(_collection->getFamilyInfo());
  mfm->setGroupInfo(_collection->getGroupInfo());
  
  //without filter mfm->setFamilyFieldArr(-1,(_collection->getFaceFamilyIds())[idomain]);
  
  string cle=Cle1ToStr("faceFamily_toArray",idomain);
  if (_collection->getMapDataArrayInt().find(cle)!=_collection->getMapDataArrayInt().end())
    {
      DataArrayInt* fam=_collection->getMapDataArrayInt().find(cle)->second;
      DataArrayInt* famFilter=0;
      if (filter!=0)
        {
          int* index=filter->getPointer();
          int nbTuples=filter->getNbOfElems();
          //not the good one...buildPartOfMySelf do not exist for DataArray 
          //Filter=fam->renumberAndReduce(index, filter->getNbOfElems());
          famFilter=DataArrayInt::New();
          famFilter->alloc(nbTuples,1);
          int* pfamFilter=famFilter->getPointer();
          int* pfam=fam->getPointer();
          for (int i=0; i<nbTuples; i++) pfamFilter[i]=pfam[index[i]];
          fam=famFilter;
          mfm->setFamilyFieldArr(-1,fam);
          famFilter->decrRef();
        }
      //cout<<"proc "<<MyGlobals::_Rank<<"cvw111 "<<nbTuples<<endl;
      //mfm->setFamilyFieldArr(-1,fam);
      //if (famFilter!=0) famFilter->decrRef();
    }
  
  /*example visualisation of filter
    if (_collection->getMapDataArrayInt().find(cle)!=_collection->getMapDataArrayInt().end())
    {
    DataArrayInt* fam=_collection->getMapDataArrayInt().find(cle)->second;
    string cle2=Cle1ToStr("filterNotFaceOnCell",idomain);
    if (_collection->getMapDataArrayInt().find(cle2)!=_collection->getMapDataArrayInt().end())
    {
    DataArrayInt* filter=_collection->getMapDataArrayInt().find(cle2)->second;
    int* index=filter->getPointer();
    int* pfam=fam->getPointer();
    for (int i=0; i<filter->getNbOfElems(); i++) pfam[index[i]]=-500;
    }
    mfm->setFamilyFieldArr(-1,fam);
    //mfm->setFamilyFieldArr(-1,_collection->getMapDataArrayInt().find(cle)->second);
    }
  */
  
  cle=Cle1ToStr("cellFamily_toArray",idomain);
  if (_collection->getMapDataArrayInt().find(cle)!=_collection->getMapDataArrayInt().end())
    mfm->setFamilyFieldArr(0,_collection->getMapDataArrayInt().find(cle)->second);
  
  mfm->write(distfilename.c_str(),0);
  cle="/inewFieldDouble="+IntToStr(idomain)+"/";
    
  map<string,ParaMEDMEM::DataArrayDouble*>::iterator it;
  int nbfFieldFound=0;
  for (it=_collection->getMapDataArrayDouble().begin() ; it!=_collection->getMapDataArrayDouble().end(); it++)
    {
      string desc=(*it).first;
      size_t found=desc.find(cle);
      if (found==string::npos) continue;
      if (MyGlobals::_Verbose>20) cout<<"proc "<<MyGlobals::_Rank<<" : write field "<<desc<<endl;
      string meshName, fieldName;
      int typeField, DT, IT, entity;
      FieldShortDescriptionToData(desc, fieldName, typeField, entity, DT, IT);
      double time=StrToDouble(ExtractFromDescription(desc, "time="));
      int typeData=StrToInt(ExtractFromDescription(desc, "typeData="));
      //int nbPtGauss=StrToInt(ExtractFromDescription(desc, "nbPtGauss="));
      string entityName=ExtractFromDescription(desc, "entityName=");
      MEDCouplingFieldDouble* field=0;
      if (typeData!=6)
        {
          cout<<"WARNING : writeMedFile : typeData "<<typeData<<" not implemented for fields\n";
          continue;
        }
      if (entityName=="MED_CELL")
        {
          //there is a field of idomain to write
          field=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
        }
      if (entityName=="MED_NODE_ELEMENT")
        {
          //there is a field of idomain to write
          field=MEDCouplingFieldDouble::New(ON_GAUSS_NE,ONE_TIME);
        }
      if (!field)
        {
          cout<<"WARNING : writeMedFile : entityName "<<entityName<<" not implemented for fields\n";
          continue;
        }
      nbfFieldFound++;
      field->setName(fieldName.c_str());
      field->setMesh(mfm->getLevel0Mesh(false));
      DataArrayDouble* da=(*it).second;
    
      //get information for components etc..
      vector<string> r1;
      r1=SelectTagsInVectorOfString(MyGlobals::_General_Informations,"fieldName="+fieldName);
      r1=SelectTagsInVectorOfString(r1,"typeField="+IntToStr(typeField));
      r1=SelectTagsInVectorOfString(r1,"DT="+IntToStr(DT));
      r1=SelectTagsInVectorOfString(r1,"IT="+IntToStr(IT));
      //not saved in file? field->setDescription(ExtractFromDescription(r1[0], "fieldDescription=").c_str());
      int nbc=StrToInt(ExtractFromDescription(r1[0], "nbComponents="));
      //double time=StrToDouble(ExtractFromDescription(r1[0], "time="));
      if (nbc==da->getNumberOfComponents())
        {
          for (int i=0; i<nbc; i++) 
            da->setInfoOnComponent(i,ExtractFromDescription(r1[0], "componentInfo"+IntToStr(i)+"=").c_str());
        }
      else
        {
          cerr<<"Problem On field "<<fieldName<<" : number of components unexpected "<<da->getNumberOfComponents()<<endl;
        }
    
      field->setArray(da);
      field->setTime(time,DT,IT);
      field->checkCoherency();
      try
        {
          MEDLoader::WriteField(distfilename.c_str(),field,false);
          //if entityName=="MED_NODE_ELEMENT"
          //AN INTERP_KERNEL::EXCEPTION HAS BEEN THROWN : Not implemented other profile fitting from already written mesh for fields than on NODES and on CELLS.**********
          //modification MEDLoader.cxx done
        }
      catch(INTERP_KERNEL::Exception& e)
        {
          //cout trying rewrite all data, only one field defined
          string tmp,newName=distfilename;
          tmp+="_"+fieldName+"_"+IntToStr(nbfFieldFound)+".med";
          newName.replace(newName.find(".med"),4,tmp);
          cout<<"WARNING : writeMedFile : new file name with only one field :"<<newName<<endl;
          MEDLoader::WriteField(newName.c_str(),field,true);
        }
      //cout<<"proc "<<MyGlobals::_Rank<<" : write field "<<cle<<" done"<<endl;
    }
  mfm->decrRef();
}
