//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#include "MEDFileMesh.hxx"

extern "C" {
#include "med.h"
}
//MEDPARTITIONER includes
#include "MEDPARTITIONER_ParallelTopology.hxx"
#include "MEDPARTITIONER_MESHCollectionDriver.hxx"
#include "MEDPARTITIONER_MESHCollection.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_utils.hxx"

using namespace MEDPARTITIONER;
using namespace std;

//template inclusion
//#include "MEDPARTITIONER_MESHCollectionDriver.H"

// med_geometrie_element typmai[MED_NBR_GEOMETRIE_MAILLE+2] = { MED_POINT1,
//                                                              MED_SEG2,
//                                                              MED_SEG3,
//                                                              MED_TRIA3,
//                                                              MED_TRIA6,
//                                                              MED_QUAD4,
//                                                              MED_QUAD8,
//                                                              MED_TETRA4,
//                                                              MED_TETRA10,
//                                                              MED_HEXA8,
//                                                              MED_HEXA20,
//                                                              MED_PENTA6,
//                                                              MED_PENTA15,
//                                                              MED_PYRA5,
//                                                              MED_PYRA13,
//                                                              MED_POLYGONE,
//                                                              MED_POLYEDRE };

MESHCollectionDriver::MESHCollectionDriver(MESHCollection* collection):_collection(collection)
{
}


/*!reads a unique MED File v>=2.1
 * and mounts the corresponding mesh in memory
 *\param filename binary file
 *\param meshname mesh name in the MED file
 * */
int MESHCollectionDriver::readSeq(const char* filename, const char* meshname)
{
  cout<<"readSeq"<<endl;
  MyGlobals::_fileNames.resize(1);
  MyGlobals::_fileNames[0]=string(filename);

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
  
// (_collection->getGroupMeshes()).resize(groupNames.size());
  
//   for (int i=0; i< groupNames.size();i++)
//     {
//       vector<string> myGroup;
//       myGroup.push_back(groupNames[i]);
//       (_collection->getGroupMeshes())[i].push_back(MEDLoader::ReadUMeshFromGroups(filename,meshname,-1,myGroup));
//     }


  (_collection->getCZ()).clear();
  /*cvw 
  vector<int*> cellglobal,nodeglobal,faceglobal;
  cellglobal.resize(1);
  nodeglobal.resize(1);
  faceglobal.resize(1);
  cellglobal[0]=0;
  nodeglobal[0]=0;
  faceglobal[0]=0;
  //creation of topology from mesh 
  //connectzone argument is 0
  ParallelTopology* aPT = new ParallelTopology
    ((_collection->getMesh()), (_collection->getCZ()), cellglobal, nodeglobal, faceglobal);
  */
  
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

void MESHCollectionDriver::readSubdomain(vector<int*>& cellglobal, //cvwat03
                                         vector<int*>& faceglobal,
                                         vector<int*>& nodeglobal, int idomain)
{
  string meshname=MyGlobals::_meshNames[idomain];
  string file=MyGlobals::_fileNames[idomain];

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
    (_collection->getMesh())[idomain]=createEmptyMEDCouplingUMesh(); // or 0 if you want tests;
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
    //ParaMEDMEM::MEDCouplingUMesh *umesh=ParaMEDMEM::MEDCouplingUMesh::New(); //empty one
    //umesh->setMeshDimension(3);
    //umesh->allocateCells(0);
    //int nb=umesh->getNumberOfCells(); //no use if no allocateCells(0)! because thrown exception
    //cout<<"\nempty mesh"<<nb<<endl;
    
    (_collection->getFaceMesh())[idomain]=createEmptyMEDCouplingUMesh(); // or 0 if you want test;
    ParaMEDMEM::DataArrayInt* empty=ParaMEDMEM::DataArrayInt::New();
    (_collection->getFaceFamilyIds())[idomain]=empty;
    if (MyGlobals::_verbose>10) cout<<"proc "<<MyGlobals::_rank<<" : NO LevelM1Mesh (Faces)\n";
  }
  
  //reading groups
  _collection->getFamilyInfo()=mfm->getFamilyInfo();
  _collection->getGroupInfo()=mfm->getGroupInfo();

  mfm->decrRef();
  
  vector<string> localInformation;
  string str;
  localInformation.push_back(str+"ioldDomain="+intToStr(idomain));
  localInformation.push_back(str+"meshName="+meshname);
  MyGlobals::_generalInformations.push_back(serializeFromVectorOfString(localInformation));
  vector<string> localFields=browseAllFieldsOnMesh(file, meshname, idomain); //cvwat07
  if (localFields.size()>0) 
    MyGlobals::_fieldDescriptions.push_back(serializeFromVectorOfString(localFields));
  //cout<< "End Reading "<<meshname<<" in "<<file<<endl;
}


void MESHCollectionDriver::readSubdomain(int idomain)
{
  string meshname=MyGlobals::_meshNames[idomain];
  string file=MyGlobals::_fileNames[idomain];

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
    (_collection->getMesh())[idomain]=createEmptyMEDCouplingUMesh(); // or 0 if you want tests;
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
    (_collection->getFaceMesh())[idomain]=createEmptyMEDCouplingUMesh(); // or 0 if you want test;
    ParaMEDMEM::DataArrayInt* empty=ParaMEDMEM::DataArrayInt::New();
    (_collection->getFaceFamilyIds())[idomain]=empty;
    if (MyGlobals::_verbose>10) cout<<"proc "<<MyGlobals::_rank<<" : NO LevelM1Mesh (Faces)\n";
  }
  
  //reading groups
  _collection->getFamilyInfo()=mfm->getFamilyInfo();
  _collection->getGroupInfo()=mfm->getGroupInfo();

  mfm->decrRef();
  
  vector<string> localInformation;
  string str;
  localInformation.push_back(str+"ioldDomain="+intToStr(idomain));
  localInformation.push_back(str+"meshName="+meshname);
  MyGlobals::_generalInformations.push_back(serializeFromVectorOfString(localInformation));
  vector<string> localFields=browseAllFieldsOnMesh(file, meshname, idomain); //cvwat07
  if (localFields.size()>0) 
    MyGlobals::_fieldDescriptions.push_back(serializeFromVectorOfString(localFields));
  //cout<< "End Reading "<<meshname<<" in "<<file<<endl;
}


void MESHCollectionDriver::writeMedFile(int idomain, const string& distfilename)
{
  vector<const ParaMEDMEM::MEDCouplingUMesh*> meshes;
  ParaMEDMEM::MEDCouplingUMesh* cellMesh=_collection->getMesh(idomain);
  ParaMEDMEM::MEDCouplingUMesh* faceMesh=_collection->getFaceMesh(idomain);
  ParaMEDMEM::MEDCouplingUMesh* faceMeshFilter=0;
  
  string finalMeshName=extractFromDescription(MyGlobals::_generalInformations[0], "finalMeshName=");
  string cleFilter=cle1ToStr("filterFaceOnCell",idomain);
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
  if (MyGlobals::_creates_boundary_faces>0)
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
    int rang=MyGlobals::_rank;
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

  //cout<<"familyInfo :\n"<<reprMapOfStringInt(_collection->getFamilyInfo())<<endl;
  //cout<<"groupInfo :\n"<<reprMapOfStringVectorOfString(_collection->getGroupInfo())<<endl;
  
  ParaMEDMEM::MEDFileUMesh* mfm=ParaMEDMEM::MEDFileUMesh::New(distfilename.c_str(), _collection->getMesh(idomain)->getName());
	
  /*example of adding new family
  (_collection->getFamilyInfo())["FaceNotOnCell"]=-500;
  vector<string> FaceNotOnCell;
  FaceNotOnCell.push_back("FaceNotOnCell");
  (_collection->getGroupInfo())["FaceNotOnCell"]=FaceNotOnCell;
  */
  
  mfm->setFamilyInfo(_collection->getFamilyInfo());
  mfm->setGroupInfo(_collection->getGroupInfo());
  
  //cvwat08 
  //without filter mfm->setFamilyFieldArr(-1,(_collection->getFaceFamilyIds())[idomain]);
  
  string cle=cle1ToStr("faceFamily_toArray",idomain);
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
    //cout<<"proc "<<MyGlobals::_rank<<"cvw111 "<<nbTuples<<endl;
    //mfm->setFamilyFieldArr(-1,fam);
   // if (famFilter!=0) famFilter->decrRef();
  }
  
  /*example visualisation of filter
  if (_collection->getMapDataArrayInt().find(cle)!=_collection->getMapDataArrayInt().end())
  {
    DataArrayInt* fam=_collection->getMapDataArrayInt().find(cle)->second;
    string cle2=cle1ToStr("filterNotFaceOnCell",idomain);
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
  
  cle=cle1ToStr("cellFamily_toArray",idomain);
  if (_collection->getMapDataArrayInt().find(cle)!=_collection->getMapDataArrayInt().end())
    mfm->setFamilyFieldArr(0,_collection->getMapDataArrayInt().find(cle)->second);
  
  mfm->write(distfilename.c_str(),0);
  cle="/inewFieldDouble="+intToStr(idomain)+"/";
    
  map<string,ParaMEDMEM::DataArrayDouble*>::iterator it;
  int nbfFieldFound=0;
  for (it=_collection->getMapDataArrayDouble().begin() ; it!=_collection->getMapDataArrayDouble().end(); it++)
  {
    string desc=(*it).first;
    size_t found=desc.find(cle);
    if (found==string::npos) continue;
    if (MyGlobals::_verbose>20) cout<<"proc "<<MyGlobals::_rank<<" : write field "<<desc<<endl;
    string meshName, fieldName;
    int typeField, DT, IT, entity;
    fieldShortDescriptionToData(desc, fieldName, typeField, entity, DT, IT);
    double time=strToDouble(extractFromDescription(desc, "time="));
    int typeData=strToInt(extractFromDescription(desc, "typeData="));
    //int nbPtGauss=strToInt(extractFromDescription(desc, "nbPtGauss="));
    string entityName=extractFromDescription(desc, "entityName=");
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
    r1=selectTagsInVectorOfString(MyGlobals::_generalInformations,"fieldName="+fieldName);
    r1=selectTagsInVectorOfString(r1,"typeField="+intToStr(typeField));
    r1=selectTagsInVectorOfString(r1,"DT="+intToStr(DT));
    r1=selectTagsInVectorOfString(r1,"IT="+intToStr(IT));
    //not saved in file? field->setDescription(extractFromDescription(r1[0], "fieldDescription=").c_str());
    int nbc=strToInt(extractFromDescription(r1[0], "nbComponents="));
    //double time=strToDouble(extractFromDescription(r1[0], "time="));
    if (nbc==da->getNumberOfComponents())
    {
      for (int i=0; i<nbc; i++) 
        da->setInfoOnComponent(i,extractFromDescription(r1[0], "componentInfo"+intToStr(i)+"=").c_str());
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
      tmp+="_"+fieldName+"_"+intToStr(nbfFieldFound)+".med";
      newName.replace(newName.find(".med"),4,tmp);
      cout<<"WARNING : writeMedFile : new file name with only one field :"<<newName<<endl;
      MEDLoader::WriteField(newName.c_str(),field,true);
    }
    //cout<<"proc "<<MyGlobals::_rank<<" : write field "<<cle<<" done"<<endl;
  }

  mfm->decrRef();

}

/*
void writeFieldNodeCellTryingToFitExistingMesh(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f)
{
  med_int numdt,numo;
  med_float dt;
  int nbComp=f->getNumberOfComponents();
  med_idt fid=appendFieldSimpleAtt(fileName,f,numdt,numo,dt);
  std::list<MEDLoader::MEDFieldDoublePerCellType> split;
  prepareCellFieldDoubleForWriting(f,thisMeshCellIdsPerType,split);
  const double *pt=f->getArray()->getConstPointer();
  int number=0;
  for(std::list<MEDLoader::MEDFieldDoublePerCellType>::const_iterator iter=split.begin();iter!=split.end();iter++)
    {
      INTERP_KERNEL::AutoPtr<char> nommaa=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      MEDLoaderBase::safeStrCpy(f->getMesh()->getName(),MED_NAME_SIZE,nommaa,MEDLoader::_TOO_LONG_STR);
      INTERP_KERNEL::AutoPtr<char> profileName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      std::ostringstream oss; oss << "Pfl" << f->getName() << "_" << number++;
      MEDLoaderBase::safeStrCpy(oss.str().c_str(),MED_NAME_SIZE,profileName,MEDLoader::_TOO_LONG_STR);
      const std::vector<int>& ids=(*iter).getCellIdPerType();
      int *profile=new int [ids.size()];
      std::transform(ids.begin(),ids.end(),profile,std::bind2nd(std::plus<int>(),1));
      MEDprofileWr(fid,profileName,ids.size(),profile);
      delete [] profile;
      MEDfieldValueWithProfileWr(fid,f->getName(),numdt,numo,dt,MED_NODE_CELL,typmai3[(int)(*iter).getType()],MED_COMPACT_PFLMODE,profileName,
                                 MED_NO_LOCALIZATION,MED_FULL_INTERLACE,MED_ALL_CONSTITUENT,(*iter).getNbOfTuple(),(const unsigned char*)pt);
      pt+=(*iter).getNbOfTuple()*nbComp;
    }
  MEDfileClose(fid);
}*/

