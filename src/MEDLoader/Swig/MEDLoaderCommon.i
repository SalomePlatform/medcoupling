// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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

%module MEDLoader

#define MEDCOUPLING_EXPORT
#define MEDLOADER_EXPORT

#ifdef WITH_DOCSTRINGS
%include "MEDLoader_doc.i"
#endif

%include "MEDCouplingCommon.i"

%{
#include "MEDLoader.hxx"
#include "MEDFileJoint.hxx"
#include "MEDFileMesh.hxx"
#include "MEDFileField.hxx"
#include "MEDFileParameter.hxx"
#include "MEDFileData.hxx"
#include "MEDFileMeshReadSelector.hxx"
#include "MEDFileFieldOverView.hxx"
#include "MEDLoaderTypemaps.i"
#include "SauvReader.hxx"
#include "SauvWriter.hxx"

using namespace ParaMEDMEM;
%}

#if SWIG_VERSION >= 0x010329
%template()  std::vector<std::string>;
#endif

%typemap(out) ParaMEDMEM::MEDFileMesh*
{
  $result=convertMEDFileMesh($1,$owner);
}

%typemap(out) ParaMEDMEM::MEDFileParameter1TS*
{
  $result=convertMEDFileParameter1TS($1,$owner);
}

%typemap(out) ParaMEDMEM::MEDFileAnyTypeFieldMultiTS*
{
  $result=convertMEDFileFieldMultiTS($1,$owner);
}

%typemap(out) ParaMEDMEM::MEDFileAnyTypeField1TS*
{
  $result=convertMEDFileField1TS($1,$owner);
}

%typemap(out) ParaMEDMEM::MEDMeshMultiLev*
{
  $result=convertMEDMeshMultiLev($1,$owner);
}

%newobject MEDLoader::ReadUMeshFromFamilies;
%newobject MEDLoader::ReadUMeshFromGroups;
%newobject MEDLoader::ReadUMeshFromFile;
%newobject MEDLoader::ReadMeshFromFile;
%newobject MEDLoader::ReadField;
%newobject MEDLoader::ReadFieldCell;
%newobject MEDLoader::ReadFieldNode;
%newobject MEDLoader::ReadFieldGauss;
%newobject MEDLoader::ReadFieldGaussNE;
%newobject ParaMEDMEM::MEDFileMesh::New;
%newobject ParaMEDMEM::MEDFileMesh::createNewEmpty;
%newobject ParaMEDMEM::MEDFileMesh::deepCpy;
%newobject ParaMEDMEM::MEDFileMesh::shallowCpy;
%newobject ParaMEDMEM::MEDFileMesh::getGenMeshAtLevel;
%newobject ParaMEDMEM::MEDFileMesh::__getitem__;
%newobject ParaMEDMEM::MEDFileMesh::getGroupArr;
%newobject ParaMEDMEM::MEDFileMesh::getGroupsArr;
%newobject ParaMEDMEM::MEDFileMesh::getFamilyArr;
%newobject ParaMEDMEM::MEDFileMesh::getFamiliesArr;
%newobject ParaMEDMEM::MEDFileMesh::getNodeGroupArr;
%newobject ParaMEDMEM::MEDFileMesh::getNodeGroupsArr;
%newobject ParaMEDMEM::MEDFileMesh::getNodeFamilyArr;
%newobject ParaMEDMEM::MEDFileMesh::getNodeFamiliesArr;
%newobject ParaMEDMEM::MEDFileMesh::getAllFamiliesIdsReferenced;
%newobject ParaMEDMEM::MEDFileMesh::computeAllFamilyIdsInUse;
%newobject ParaMEDMEM::MEDFileData::getJoints;
%newobject ParaMEDMEM::MEDFileStructuredMesh::getImplicitFaceMesh;
%newobject ParaMEDMEM::MEDFileUMesh::New;
%newobject ParaMEDMEM::MEDFileUMesh::LoadPartOf;
%newobject ParaMEDMEM::MEDFileUMesh::getCoords;
%newobject ParaMEDMEM::MEDFileUMesh::getPartDefAtLevel;
%newobject ParaMEDMEM::MEDFileUMesh::getGroup;
%newobject ParaMEDMEM::MEDFileUMesh::getGroups;
%newobject ParaMEDMEM::MEDFileUMesh::getFamily;
%newobject ParaMEDMEM::MEDFileUMesh::getFamilies;
%newobject ParaMEDMEM::MEDFileUMesh::getMeshAtLevel;
%newobject ParaMEDMEM::MEDFileUMesh::getLevel0Mesh;
%newobject ParaMEDMEM::MEDFileUMesh::getLevelM1Mesh;
%newobject ParaMEDMEM::MEDFileUMesh::getLevelM2Mesh;
%newobject ParaMEDMEM::MEDFileUMesh::getLevelM3Mesh;
%newobject ParaMEDMEM::MEDFileUMesh::getDirectUndergroundSingleGeoTypeMesh;
%newobject ParaMEDMEM::MEDFileUMesh::extractFamilyFieldOnGeoType;
%newobject ParaMEDMEM::MEDFileUMesh::extractNumberFieldOnGeoType;
%newobject ParaMEDMEM::MEDFileUMesh::zipCoords;
%newobject ParaMEDMEM::MEDFileUMesh::buildExtrudedMesh;
%newobject ParaMEDMEM::MEDFileUMesh::linearToQuadratic;
%newobject ParaMEDMEM::MEDFileUMesh::quadraticToLinear;
%newobject ParaMEDMEM::MEDFileCMesh::New;
%newobject ParaMEDMEM::MEDFileCurveLinearMesh::New;
%newobject ParaMEDMEM::MEDFileMeshMultiTS::New;
%newobject ParaMEDMEM::MEDFileMeshMultiTS::deepCpy;
%newobject ParaMEDMEM::MEDFileMeshMultiTS::getOneTimeStep;
%newobject ParaMEDMEM::MEDFileMeshes::New;
%newobject ParaMEDMEM::MEDFileMeshes::deepCpy;
%newobject ParaMEDMEM::MEDFileMeshes::getMeshAtPos;
%newobject ParaMEDMEM::MEDFileMeshes::getMeshWithName;
%newobject ParaMEDMEM::MEDFileMeshes::__getitem__;
%newobject ParaMEDMEM::MEDFileMeshes::__iter__;

%newobject ParaMEDMEM::MEDFileFields::New;
%newobject ParaMEDMEM::MEDFileFields::LoadPartOf;
%newobject ParaMEDMEM::MEDFileFields::LoadSpecificEntities;
%newobject ParaMEDMEM::MEDFileFields::deepCpy;
%newobject ParaMEDMEM::MEDFileFields::shallowCpy;
%newobject ParaMEDMEM::MEDFileFields::getFieldWithName;
%newobject ParaMEDMEM::MEDFileFields::getFieldAtPos;
%newobject ParaMEDMEM::MEDFileFields::partOfThisLyingOnSpecifiedMeshName;
%newobject ParaMEDMEM::MEDFileFields::partOfThisLyingOnSpecifiedTimeSteps;
%newobject ParaMEDMEM::MEDFileFields::partOfThisNotLyingOnSpecifiedTimeSteps;
%newobject ParaMEDMEM::MEDFileFields::__iter__;

%newobject ParaMEDMEM::MEDFileAnyTypeFieldMultiTS::New;
%newobject ParaMEDMEM::MEDFileAnyTypeFieldMultiTS::deepCpy;
%newobject ParaMEDMEM::MEDFileAnyTypeFieldMultiTS::shallowCpy;
%newobject ParaMEDMEM::MEDFileAnyTypeFieldMultiTS::getTimeStepAtPos;
%newobject ParaMEDMEM::MEDFileAnyTypeFieldMultiTS::getTimeStep;
%newobject ParaMEDMEM::MEDFileAnyTypeFieldMultiTS::getTimeStepGivenTime;
%newobject ParaMEDMEM::MEDFileAnyTypeFieldMultiTS::__iter__;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::New;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::LoadSpecificEntities;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::getFieldAtLevel;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::getFieldAtTopLevel;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::getFieldOnMeshAtLevel;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::getFieldAtLevelOld;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::getUndergroundDataArray;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::convertToInt;
%newobject ParaMEDMEM::MEDFileIntFieldMultiTS::New;
%newobject ParaMEDMEM::MEDFileIntFieldMultiTS::LoadSpecificEntities;
%newobject ParaMEDMEM::MEDFileIntFieldMultiTS::getUndergroundDataArray;
%newobject ParaMEDMEM::MEDFileIntFieldMultiTS::convertToDouble;

%newobject ParaMEDMEM::MEDFileAnyTypeField1TS::New;
%newobject ParaMEDMEM::MEDFileAnyTypeField1TS::shallowCpy;
%newobject ParaMEDMEM::MEDFileAnyTypeField1TS::deepCpy;
%newobject ParaMEDMEM::MEDFileField1TS::New;
%newobject ParaMEDMEM::MEDFileField1TS::getFieldAtLevel;
%newobject ParaMEDMEM::MEDFileField1TS::getFieldAtTopLevel;
%newobject ParaMEDMEM::MEDFileField1TS::getFieldOnMeshAtLevel;
%newobject ParaMEDMEM::MEDFileField1TS::getFieldAtLevelOld;
%newobject ParaMEDMEM::MEDFileField1TS::getUndergroundDataArray;
%newobject ParaMEDMEM::MEDFileField1TS::convertToInt;
%newobject ParaMEDMEM::MEDFileIntField1TS::New;
%newobject ParaMEDMEM::MEDFileIntField1TS::getUndergroundDataArray;
%newobject ParaMEDMEM::MEDFileIntField1TS::convertToDouble;

%newobject ParaMEDMEM::MEDFileData::New;
%newobject ParaMEDMEM::MEDFileData::deepCpy;
%newobject ParaMEDMEM::MEDFileData::getMeshes;
%newobject ParaMEDMEM::MEDFileData::getFields;
%newobject ParaMEDMEM::MEDFileData::getParams;

%newobject ParaMEDMEM::MEDFileParameterDouble1TS::New;
%newobject ParaMEDMEM::MEDFileParameterDouble1TS::deepCpy;
%newobject ParaMEDMEM::MEDFileParameterMultiTS::New;
%newobject ParaMEDMEM::MEDFileParameterMultiTS::deepCpy;
%newobject ParaMEDMEM::MEDFileParameterMultiTS::getTimeStepAtPos;
%newobject ParaMEDMEM::MEDFileParameterMultiTS::__getitem__;
%newobject ParaMEDMEM::MEDFileParameters::New;
%newobject ParaMEDMEM::MEDFileParameters::deepCpy;
%newobject ParaMEDMEM::MEDFileParameters::getParamAtPos;
%newobject ParaMEDMEM::MEDFileParameters::getParamWithName;
%newobject ParaMEDMEM::MEDFileParameters::__getitem__;

%newobject ParaMEDMEM::MEDFileJointCorrespondence::New;
%newobject ParaMEDMEM::MEDFileJointCorrespondence::deepCpy;
%newobject ParaMEDMEM::MEDFileJointCorrespondence::shallowCpy;
%newobject ParaMEDMEM::MEDFileJointOneStep::New;
%newobject ParaMEDMEM::MEDFileJointOneStep::deepCpy;
%newobject ParaMEDMEM::MEDFileJointOneStep::shallowCpy;
%newobject ParaMEDMEM::MEDFileJoint::New;
%newobject ParaMEDMEM::MEDFileJoint::deepCpy;
%newobject ParaMEDMEM::MEDFileJoint::shallowCpy;
%newobject ParaMEDMEM::MEDFileJoints::New;
%newobject ParaMEDMEM::MEDFileJoints::deepCpy;
%newobject ParaMEDMEM::MEDFileJoints::getJointAtPos;
%newobject ParaMEDMEM::MEDFileJoints::getJointWithName;
%newobject ParaMEDMEM::MEDFileJoints::__getitem__;

%newobject ParaMEDMEM::SauvWriter::New;
%newobject ParaMEDMEM::SauvReader::New;
%newobject ParaMEDMEM::SauvReader::loadInMEDFileDS;

%newobject ParaMEDMEM::MEDFileMeshStruct::New;
%newobject ParaMEDMEM::MEDMeshMultiLev::prepare;
%newobject ParaMEDMEM::MEDMeshMultiLev::buildDataArray;
%newobject ParaMEDMEM::MEDFileFastCellSupportComparator::New;
%newobject ParaMEDMEM::MEDFileFastCellSupportComparator::buildFromScratchDataSetSupport;

%feature("unref") MEDFileMesh "$this->decrRef();"
%feature("unref") MEDFileUMesh "$this->decrRef();"
%feature("unref") MEDFileCMesh "$this->decrRef();"
%feature("unref") MEDFileMeshMultiTS "$this->decrRef();"
%feature("unref") MEDFileMeshes "$this->decrRef();"
%feature("unref") MEDFileFieldLoc "$this->decrRef();"
%feature("unref") MEDFileAnyTypeField1TS "$this->decrRef();"
%feature("unref") MEDFileField1TS "$this->decrRef();"
%feature("unref") MEDFileIntField1TS "$this->decrRef();"
%feature("unref") MEDFileAnyTypeFieldMultiTS "$this->decrRef();"
%feature("unref") MEDFileFieldMultiTS "$this->decrRef();"
%feature("unref") MEDFileIntFieldMultiTS "$this->decrRef();"
%feature("unref") MEDFileFields "$this->decrRef();"
%feature("unref") MEDFileParameter1TS "$this->decrRef();"
%feature("unref") MEDFileParameterDouble1TSWTI "$this->decrRef();"
%feature("unref") MEDFileParameterDouble1TS "$this->decrRef();"
%feature("unref") MEDFileParameterMultiTS "$this->decrRef();"
%feature("unref") MEDFileParameters "$this->decrRef();"
%feature("unref") MEDFileJointCorrespondence "$this->decrRef();"
%feature("unref") MEDFileJointOneStep "$this->decrRef();"
%feature("unref") MEDFileJoint "$this->decrRef();"
%feature("unref") MEDFileJoints "$this->decrRef();"
%feature("unref") MEDFileData "$this->decrRef();"
%feature("unref") SauvReader "$this->decrRef();"
%feature("unref") SauvWriter "$this->decrRef();"
%feature("unref") MEDFileFastCellSupportComparator "$this->decrRef();"
%feature("unref") MEDMeshMultiLev "$this->decrRef();"
%feature("unref") MEDUMeshMultiLev "$this->decrRef();"
%feature("unref") MEDCMeshMultiLev "$this->decrRef();"
%feature("unref") MEDCurveLinearMeshMultiLev "$this->decrRef();"
%feature("unref") MEDFileMeshStruct "$this->decrRef();"

class MEDLoader
{
public:
  static bool HasXDR();
  static std::string MEDFileVersionStr();
  static void SetEpsilonForNodeComp(double val) throw(INTERP_KERNEL::Exception);
  static void SetCompPolicyForCell(int val) throw(INTERP_KERNEL::Exception);
  static void SetTooLongStrPolicy(int val) throw(INTERP_KERNEL::Exception);
  static void CheckFileForRead(const std::string& fileName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshNames(const std::string& fileName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshNamesOnField(const std::string& fileName, const std::string& fieldName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshGroupsNames(const std::string& fileName, const std::string& meshName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshFamiliesNames(const std::string& fileName, const std::string& meshName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshFamiliesNamesOnGroup(const std::string& fileName, const std::string& meshName, const std::string& grpName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshGroupsNamesOnFamily(const std::string& fileName, const std::string& meshName, const std::string& famName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetAllFieldNamesOnMesh(const std::string& fileName, const std::string& meshName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetAllFieldNames(const std::string& fileName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetFieldNamesOnMesh(ParaMEDMEM::TypeOfField type, const std::string& fileName, const std::string& meshName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetCellFieldNamesOnMesh(const std::string& fileName, const std::string& meshName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetNodeFieldNamesOnMesh(const std::string& fileName, const std::string& meshName) throw(INTERP_KERNEL::Exception);
  static double GetTimeAttachedOnFieldIteration(const std::string& fileName, const std::string& fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static void AssignStaticWritePropertiesTo(ParaMEDMEM::MEDFileWritable& obj) throw(INTERP_KERNEL::Exception);
  %extend
     {
       static PyObject *MEDFileVersion()
       {
         int major,minor,release;
         MEDLoader::MEDFileVersion(major,minor,release);
         PyObject *ret(PyTuple_New(3));
         PyTuple_SetItem(ret,0,SWIG_From_int(major));
         PyTuple_SetItem(ret,1,SWIG_From_int(minor));
         PyTuple_SetItem(ret,2,SWIG_From_int(release));
         return ret;
       }

       static PyObject *GetFieldIterations(ParaMEDMEM::TypeOfField type, const std::string& fileName, const std::string& meshName, const std::string& fieldName) throw(INTERP_KERNEL::Exception)
       {
         std::vector< std::pair<int,int> > res=MEDLoader::GetFieldIterations(type,fileName,meshName,fieldName);
         PyObject *ret=PyList_New(res.size());
         int rk=0;
         for(std::vector< std::pair<int,int> >::const_iterator iter=res.begin();iter!=res.end();iter++,rk++)
           {
             PyObject *elt=PyTuple_New(2);
             PyTuple_SetItem(elt,0,SWIG_From_int((*iter).first));
             PyTuple_SetItem(elt,1,SWIG_From_int((*iter).second));
             PyList_SetItem(ret,rk,elt);
           }
         return ret;
       }

       static PyObject *GetAllFieldIterations(const std::string& fileName, const std::string& fieldName) throw(INTERP_KERNEL::Exception)
       {
         std::vector< std::pair< std::pair<int,int>, double> > res=MEDLoader::GetAllFieldIterations(fileName,fieldName);
         PyObject *ret=PyList_New(res.size());
         int rk=0;
         for(std::vector< std::pair< std::pair<int,int>, double> >::const_iterator iter=res.begin();iter!=res.end();iter++,rk++)
           {
             PyObject *elt=PyTuple_New(3);
             PyTuple_SetItem(elt,0,SWIG_From_int((*iter).first.first));
             PyTuple_SetItem(elt,1,SWIG_From_int((*iter).first.second));
             PyTuple_SetItem(elt,2,SWIG_From_double((*iter).second));
             PyList_SetItem(ret,rk,elt);
           }
         return ret;
       }

       static PyObject *GetCellFieldIterations(const std::string& fileName, const std::string& meshName, const std::string& fieldName) throw(INTERP_KERNEL::Exception)
       {
         std::vector< std::pair<int,int> > res=MEDLoader::GetCellFieldIterations(fileName,meshName,fieldName);
         PyObject *ret=PyList_New(res.size());
         int rk=0;
         for(std::vector< std::pair<int,int> >::const_iterator iter=res.begin();iter!=res.end();iter++,rk++)
           {
             PyObject *elt=PyTuple_New(2);
             PyTuple_SetItem(elt,0,SWIG_From_int((*iter).first));
             PyTuple_SetItem(elt,1,SWIG_From_int((*iter).second));
             PyList_SetItem(ret,rk,elt);
           }
         return ret;
       }
       static PyObject *GetNodeFieldIterations(const std::string& fileName, const std::string& meshName, const std::string& fieldName) throw(INTERP_KERNEL::Exception)
       {
         std::vector< std::pair<int,int> > res=MEDLoader::GetNodeFieldIterations(fileName,meshName,fieldName);
         PyObject *ret=PyList_New(res.size());
         int rk=0;
         for(std::vector< std::pair<int,int> >::const_iterator iter=res.begin();iter!=res.end();iter++,rk++)
           {
             PyObject *elt=PyTuple_New(2);
             PyTuple_SetItem(elt,0,SWIG_From_int((*iter).first));
             PyTuple_SetItem(elt,1,SWIG_From_int((*iter).second));
             PyList_SetItem(ret,rk,elt);
           }
         return ret;
       }
       static PyObject *GetComponentsNamesOfField(const std::string& fileName, const std::string& fieldName) throw(INTERP_KERNEL::Exception)
       {
         std::vector< std::pair<std::string,std::string> > res=MEDLoader::GetComponentsNamesOfField(fileName,fieldName);
         PyObject *ret=PyList_New(res.size());
         int rk=0;
         for(std::vector< std::pair<std::string,std::string> >::const_iterator iter=res.begin();iter!=res.end();iter++,rk++)
           {
             PyObject *elt=PyTuple_New(2);
             PyTuple_SetItem(elt,0,PyString_FromString((*iter).first.c_str()));
             PyTuple_SetItem(elt,1,PyString_FromString((*iter).second.c_str()));
             PyList_SetItem(ret,rk,elt);
           }
         return ret;
       }
       static PyObject *GetUMeshGlobalInfo(const std::string& fileName, const std::string& meshName) throw(INTERP_KERNEL::Exception)
       {
         int meshDim,spaceDim,numberOfNodes;
         std::vector< std::vector< std::pair<INTERP_KERNEL::NormalizedCellType,int> > > res=MEDLoader::GetUMeshGlobalInfo(fileName,meshName,meshDim,spaceDim,numberOfNodes);
         PyObject *ret=PyTuple_New(4);
         PyObject *elt0=PyList_New(res.size());
         int i=0;
         for(std::vector< std::vector< std::pair<INTERP_KERNEL::NormalizedCellType,int> > >::const_iterator it=res.begin();it!=res.end();it++,i++)
           {
             const std::vector< std::pair<INTERP_KERNEL::NormalizedCellType,int> >&obj2=(*it);
             int j=0;
             PyObject *elt1=PyList_New(obj2.size());
             for(std::vector< std::pair<INTERP_KERNEL::NormalizedCellType,int> >::const_iterator it2=obj2.begin();it2!=obj2.end();it2++,j++)
               {
                 PyObject *elt2=PyTuple_New(2);
                 PyTuple_SetItem(elt2,0,SWIG_From_int((int)(*it2).first));
                 PyTuple_SetItem(elt2,1,SWIG_From_int((*it2).second));
                 PyList_SetItem(elt1,j,elt2);
               }
             PyList_SetItem(elt0,i,elt1);
           }
         PyTuple_SetItem(ret,0,elt0);
         PyTuple_SetItem(ret,1,SWIG_From_int(meshDim));
         PyTuple_SetItem(ret,2,SWIG_From_int(spaceDim));
         PyTuple_SetItem(ret,3,SWIG_From_int(numberOfNodes));
         return ret;
       }
       static PyObject *ReadFieldsOnSameMesh(ParaMEDMEM::TypeOfField type, const std::string& fileName, const std::string& meshName, int meshDimRelToMax,
                                             const std::string& fieldName, PyObject *liIts) throw(INTERP_KERNEL::Exception)
       {
         std::vector<std::pair<int,int> > its=convertTimePairIdsFromPy(liIts);
         std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> res=MEDLoader::ReadFieldsOnSameMesh(type,fileName,meshName,meshDimRelToMax,fieldName,its);
         return convertFieldDoubleVecToPy(res);
       }
       static void WriteUMeshesPartition(const std::string& fileName, const std::string& meshName, PyObject *li, bool writeFromScratch) throw(INTERP_KERNEL::Exception)
       {
         std::vector<const ParaMEDMEM::MEDCouplingUMesh *> v;
         convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingUMesh *>(li,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,"MEDCouplingUMesh",v);
         MEDLoader::WriteUMeshesPartition(fileName,meshName,v,writeFromScratch);
       }
       static void WriteUMeshesPartitionDep(const std::string& fileName, const std::string& meshName, PyObject *li, bool writeFromScratch) throw(INTERP_KERNEL::Exception)
       {
         std::vector<const ParaMEDMEM::MEDCouplingUMesh *> v;
         convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingUMesh *>(li,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,"MEDCouplingUMesh",v);
         MEDLoader::WriteUMeshesPartitionDep(fileName,meshName,v,writeFromScratch);
       }
       static void WriteUMeshes(const std::string& fileName, PyObject *li, bool writeFromScratch) throw(INTERP_KERNEL::Exception)
       {
         std::vector<const ParaMEDMEM::MEDCouplingUMesh *> v;
         convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingUMesh *>(li,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,"MEDCouplingUMesh",v);
         MEDLoader::WriteUMeshes(fileName,v,writeFromScratch);
       }
       static PyObject *GetTypesOfField(const std::string& fileName, const std::string& meshName, const std::string& fieldName) throw(INTERP_KERNEL::Exception)
       {
         std::vector< ParaMEDMEM::TypeOfField > v=MEDLoader::GetTypesOfField(fileName,meshName,fieldName);
         int size=v.size();
         PyObject *ret=PyList_New(size);
         for(int i=0;i<size;i++)
           PyList_SetItem(ret,i,PyInt_FromLong((int)v[i]));
         return ret;
       }
       static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromGroups(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, PyObject *li) throw(INTERP_KERNEL::Exception)
       {
         std::vector<std::string> grps;
         converPyListToVecString(li,grps);
         return MEDLoader::ReadUMeshFromGroups(fileName,meshName,meshDimRelToMax,grps);
       }
       static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFamilies(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, PyObject *li) throw(INTERP_KERNEL::Exception)
       {
         std::vector<std::string> fams;
         converPyListToVecString(li,fams);
         return MEDLoader::ReadUMeshFromFamilies(fileName,meshName,meshDimRelToMax,fams);
       }
     }
  static ParaMEDMEM::MEDCouplingMesh *ReadMeshFromFile(const std::string& fileName, const std::string& meshName, int meshDimRelToMax=0) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingMesh *ReadMeshFromFile(const std::string& fileName, int meshDimRelToMax=0) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const std::string& fileName, const std::string& meshName, int meshDimRelToMax=0) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const std::string& fileName, int meshDimRelToMax=0) throw(INTERP_KERNEL::Exception);
  static int ReadUMeshDimFromFile(const std::string& fileName, const std::string& meshName) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadField(ParaMEDMEM::TypeOfField type, const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldCell(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldNode(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldGauss(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldGaussNE(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static void WriteMesh(const std::string& fileName, const ParaMEDMEM::MEDCouplingMesh *mesh, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteUMesh(const std::string& fileName, const ParaMEDMEM::MEDCouplingUMesh *mesh, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteUMeshDep(const std::string& fileName, const ParaMEDMEM::MEDCouplingUMesh *mesh, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteField(const std::string& fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteFieldDep(const std::string& fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteFieldUsingAlreadyWrittenMesh(const std::string& fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f) throw(INTERP_KERNEL::Exception);
};

namespace ParaMEDMEM
{
  class MEDFileWritable
  {
  public:
    void copyOptionsFrom(const MEDFileWritable& other) const;
    int getTooLongStrPolicy() const throw(INTERP_KERNEL::Exception);
    void setTooLongStrPolicy(int newVal) throw(INTERP_KERNEL::Exception);
    int getZipConnPolicy() throw(INTERP_KERNEL::Exception);
    void setZipConnPolicy(int newVal) throw(INTERP_KERNEL::Exception);
  };

  class MEDFileMeshReadSelector
  {
  public:
    MEDFileMeshReadSelector();
    MEDFileMeshReadSelector(unsigned int code);
    unsigned int getCode() const;
    void setCode(unsigned int newCode);
    bool isCellFamilyFieldReading() const;
    bool isNodeFamilyFieldReading() const;
    bool isCellNameFieldReading() const;
    bool isNodeNameFieldReading() const;
    bool isCellNumFieldReading() const;
    bool isNodeNumFieldReading() const;
    void setCellFamilyFieldReading(bool b);
    void setNodeFamilyFieldReading(bool b);
    void setCellNameFieldReading(bool b);
    void setNodeNameFieldReading(bool b);
    void setCellNumFieldReading(bool b);
    void setNodeNumFieldReading(bool b);
    %extend
    {
      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        std::ostringstream oss;
        self->reprAll(oss);
        return oss.str();
      }
      
      std::string __repr__() const throw(INTERP_KERNEL::Exception)
      {
        std::ostringstream oss; oss << "MEDFileMeshReadSelector C++ instance at " << self << " (with code=" << self->getCode() << ").";
        return oss.str();
      }
    }
  };
  class MEDFileJointCorrespondence : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileJointCorrespondence *New() throw(INTERP_KERNEL::Exception);
    static MEDFileJointCorrespondence *New(DataArrayInt* correspondence) // nodes
      throw(INTERP_KERNEL::Exception);
    static MEDFileJointCorrespondence *New(DataArrayInt* correspondence,  // cells
                                           INTERP_KERNEL::NormalizedCellType loc_geo_type,
                                           INTERP_KERNEL::NormalizedCellType rem_geo_type)
      throw(INTERP_KERNEL::Exception);
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDFileJointCorrespondence *deepCpy() const;
    MEDFileJointCorrespondence *shallowCpy() const;
    void setIsNodal(bool isNodal);
    bool getIsNodal() const;
    bool isEqual(const MEDFileJointCorrespondence *other) const;
    void setLocalGeometryType(INTERP_KERNEL::NormalizedCellType type);
    INTERP_KERNEL::NormalizedCellType getLocalGeometryType() const;
    void setRemoteGeometryType(INTERP_KERNEL::NormalizedCellType type);
    INTERP_KERNEL::NormalizedCellType getRemoteGeometryType() const;
    void setCorrespondence(DataArrayInt *corr) throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getCorrespondence() const throw(INTERP_KERNEL::Exception);
    void write(const std::string& fileName, int mode, const std::string& localMeshName, const std::string& jointName, int order, int iteration) const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const throw(INTERP_KERNEL::Exception);
    %extend
    {
      MEDFileJointCorrespondence()
      {
        return MEDFileJointCorrespondence::New();
      }
      MEDFileJointCorrespondence(DataArrayInt* correspondence) throw(INTERP_KERNEL::Exception)
      {
        return MEDFileJointCorrespondence::New(correspondence);
      }
      MEDFileJointCorrespondence(DataArrayInt* correspondence,  // cells
                                 INTERP_KERNEL::NormalizedCellType loc_geo_type,
                                 INTERP_KERNEL::NormalizedCellType rem_geo_type) throw(INTERP_KERNEL::Exception)
      {
        return MEDFileJointCorrespondence::New(correspondence, loc_geo_type, rem_geo_type);
      }

      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->simpleRepr();
      }
    }
  };

  class MEDFileJointOneStep : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileJointOneStep *New(int dt=-1, int it=-1) throw(INTERP_KERNEL::Exception);
    static MEDFileJointOneStep *New(const std::string& fileName, const std::string& mName, const std::string& jointName, int number=1) throw(INTERP_KERNEL::Exception);
    MEDFileJointOneStep *deepCpy() const;
    MEDFileJointOneStep *shallowCpy() const;
    bool isEqual(const MEDFileJointOneStep *other) const;
    void setOrder(int order);
    int getOrder() const;
    void setIteration(int it);
    int getIteration() const;
    void pushCorrespondence(MEDFileJointCorrespondence* correspondence);
    int getNumberOfCorrespondences() const;
    MEDFileJointCorrespondence *getCorrespondenceAtPos(int i) const;
    void write(const std::string& fileName, int mode, const std::string& localMeshName, const std::string& jointName) const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const throw(INTERP_KERNEL::Exception);
    %extend
    {
      MEDFileJointOneStep()
      {
        return MEDFileJointOneStep::New();
      }

      MEDFileJointOneStep(const std::string& fileName, const std::string& mName, const std::string& jointName, int number) throw(INTERP_KERNEL::Exception)
      {
        return MEDFileJointOneStep::New(fileName,mName,jointName,number);
      }

      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->simpleRepr();
      }
    }
  };
  class MEDFileJoint : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileJoint *New() throw(INTERP_KERNEL::Exception);
    static MEDFileJoint *New(const std::string& fileName, const std::string& mName, int num) throw(INTERP_KERNEL::Exception);
    static MEDFileJoint *New(const std::string& jointName, const std::string& locMeshName, const std::string& remoteMeshName, int remoteMeshNum ) throw(INTERP_KERNEL::Exception);
    MEDFileJoint *deepCpy() const;
    MEDFileJoint *shallowCpy() const;
    bool isEqual(const MEDFileJoint *other) const;
    void setLocalMeshName(const std::string& name);
    std::string getLocalMeshName() const;
    void setRemoteMeshName(const std::string& name);
    std::string getRemoteMeshName() const;
    void setDescription(const std::string& name);
    std::string getDescription() const;
    void setJointName(const std::string& name);
    std::string getJointName() const;
    bool changeJointNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    void setDomainNumber(const int& number);
    int getDomainNumber() const;
    void pushStep(MEDFileJointOneStep* step);
    int getNumberOfSteps() const;
    MEDFileJointOneStep *getStepAtPos(int i) const;
    void write(const std::string& fileName, int mode) const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    %extend
    {
      MEDFileJoint()
      {
        return MEDFileJoint::New();
      }
      
      MEDFileJoint(const std::string& fileName, const std::string& mName, int num) throw(INTERP_KERNEL::Exception)
      {
        return MEDFileJoint::New(fileName,mName,num);
      }

      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->simpleRepr();
      }
    }
  };

  class MEDFileJoints : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileJoints *New() throw(INTERP_KERNEL::Exception);
    static MEDFileJoints *New(const std::string& fileName, const std::string& meshName) throw(INTERP_KERNEL::Exception);
    MEDFileJoints *deepCpy() const;
    std::string simpleRepr() const;
    void write(const std::string& fileName, int mode) const throw(INTERP_KERNEL::Exception);
    std::string getMeshName() const;
    int getNumberOfJoints() const;
    std::vector<std::string> getJointsNames() const;
    bool changeJointNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    void resize(int newSize) throw(INTERP_KERNEL::Exception);
    void pushJoint(MEDFileJoint *joint);
    void setJointAtPos(int i, MEDFileJoint *joint) throw(INTERP_KERNEL::Exception);
    void destroyJointAtPos(int i) throw(INTERP_KERNEL::Exception);
    %extend
    {
      MEDFileJoints()
      {
        return MEDFileJoints::New();
      }
      
      MEDFileJoints(const std::string& fileName, const std::string& meshName) throw(INTERP_KERNEL::Exception)
      {
        return MEDFileJoints::New(fileName,meshName);
      }

      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->simpleRepr();
      }

      MEDFileJoint *__getitem__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        if(PyInt_Check(obj))
          {
            MEDFileJoint *ret=self->getJointAtPos(InterpreteNegativeInt((int)PyInt_AS_LONG(obj),self->getNumberOfJoints()));
            if(ret)
              ret->incrRef();
            return ret;
          }
        else if(PyString_Check(obj))
          {
            MEDFileJoint *ret=self->getJointWithName(PyString_AsString(obj));
            if(ret)
              ret->incrRef();
            return ret;
          }
        else
          throw INTERP_KERNEL::Exception("MEDFileJoints::__getitem__ : only integer or string with meshname supported !");
      }

      int __len__() const throw(INTERP_KERNEL::Exception)
      {
        return self->getNumberOfJoints();
      }

      MEDFileJoint *getJointAtPos(int i) const throw(INTERP_KERNEL::Exception)
      {
        MEDFileJoint *ret=self->getJointAtPos(i);
        if(ret)
          ret->incrRef();
        return ret;
      }

      MEDFileJoint *getJointWithName(const std::string& paramName) const throw(INTERP_KERNEL::Exception)
      {
        MEDFileJoint *ret=self->getJointWithName(paramName);
        if(ret)
          ret->incrRef();
        return ret;
      }
    }
  };

  class MEDFileMesh : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileMesh *New(const std::string& fileName, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    static MEDFileMesh *New(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    virtual MEDFileMesh *createNewEmpty() const throw(INTERP_KERNEL::Exception);
    virtual MEDFileMesh *deepCpy() const throw(INTERP_KERNEL::Exception);
    virtual MEDFileMesh *shallowCpy() const throw(INTERP_KERNEL::Exception);
    virtual void clearNonDiscrAttributes() const throw(INTERP_KERNEL::Exception);
    void setName(const std::string& name);
    std::string getName();
    std::string getUnivName() const;
    bool getUnivNameWrStatus() const;
    void setUnivNameWrStatus(bool newStatus);
    void setDescription(const std::string& name);
    std::string getDescription() const;
    void setOrder(int order);
    int getOrder() const;
    void setIteration(int it);
    int getIteration();
    void setTimeValue(double time);
    void setTime(int dt, int it, double time);
    double getTimeValue() const;
    void setTimeUnit(const std::string& unit);
    std::string getTimeUnit() const;
    virtual int getNumberOfNodes() const throw(INTERP_KERNEL::Exception);
    virtual int getNumberOfCellsAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    virtual bool hasImplicitPart() const throw(INTERP_KERNEL::Exception);
    virtual int buildImplicitPartIfAny(INTERP_KERNEL::NormalizedCellType gt) const throw(INTERP_KERNEL::Exception);
    virtual void releaseImplicitPartIfAny() const throw(INTERP_KERNEL::Exception);
    virtual std::vector<int> getFamArrNonEmptyLevelsExt() const throw(INTERP_KERNEL::Exception);
    virtual std::vector<int> getNumArrNonEmptyLevelsExt() const throw(INTERP_KERNEL::Exception);
    virtual std::vector<int> getNameArrNonEmptyLevelsExt() const throw(INTERP_KERNEL::Exception);
    virtual std::vector<int> getDistributionOfTypes(int meshDimRelToMax) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getNonEmptyLevels() const throw(INTERP_KERNEL::Exception);
    std::vector<int> getNonEmptyLevelsExt() const throw(INTERP_KERNEL::Exception);
    void write(const std::string& fileName, int mode) const throw(INTERP_KERNEL::Exception);
    int getSizeAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    //
    bool existsGroup(const std::string& groupName) const throw(INTERP_KERNEL::Exception);
    bool existsFamily(int famId) const throw(INTERP_KERNEL::Exception);
    bool existsFamily(const std::string& familyName) const throw(INTERP_KERNEL::Exception);
    void setFamilyId(const std::string& familyName, int id) throw(INTERP_KERNEL::Exception);
    void setFamilyIdUnique(const std::string& familyName, int id) throw(INTERP_KERNEL::Exception);
    void addFamily(const std::string& familyName, int id) throw(INTERP_KERNEL::Exception);
    void addFamilyOnGrp(const std::string& grpName, const std::string& famName) throw(INTERP_KERNEL::Exception);
    virtual void createGroupOnAll(int meshDimRelToMaxExt, const std::string& groupName) throw(INTERP_KERNEL::Exception);
    virtual bool keepFamIdsOnlyOnLevs(const std::vector<int>& famIds, const std::vector<int>& levs) throw(INTERP_KERNEL::Exception);
    void copyFamGrpMapsFrom(const MEDFileMesh& other) throw(INTERP_KERNEL::Exception);
    void clearGrpMap() throw(INTERP_KERNEL::Exception);
    void clearFamMap() throw(INTERP_KERNEL::Exception);
    void clearFamGrpMaps() throw(INTERP_KERNEL::Exception);
    const std::map<std::string,int>& getFamilyInfo() const throw(INTERP_KERNEL::Exception);
    const std::map<std::string, std::vector<std::string> >& getGroupInfo() const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getFamiliesOnGroup(const std::string& name) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getFamiliesOnGroups(const std::vector<std::string>& grps) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamiliesIdsOnGroup(const std::string& name) const throw(INTERP_KERNEL::Exception);
    void setFamiliesOnGroup(const std::string& name, const std::vector<std::string>& fams) throw(INTERP_KERNEL::Exception);
    void setFamiliesIdsOnGroup(const std::string& name, const std::vector<int>& famIds) throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getGroupsOnFamily(const std::string& name) const throw(INTERP_KERNEL::Exception);
    void setGroupsOnFamily(const std::string& famName, const std::vector<std::string>& grps) throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getGroupsNames() const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getFamiliesNames() const throw(INTERP_KERNEL::Exception);
    void assignFamilyNameWithGroupName() throw(INTERP_KERNEL::Exception);
    std::vector<std::string> removeEmptyGroups() throw(INTERP_KERNEL::Exception);
    void removeGroup(const std::string& name) throw(INTERP_KERNEL::Exception);
    void removeFamily(const std::string& name) throw(INTERP_KERNEL::Exception);
    std::vector<std::string> removeOrphanGroups() throw(INTERP_KERNEL::Exception);
    std::vector<std::string> removeOrphanFamilies() throw(INTERP_KERNEL::Exception);
    void removeFamiliesReferedByNoGroups() throw(INTERP_KERNEL::Exception);
    void rearrangeFamilies() throw(INTERP_KERNEL::Exception);
    void checkOrphanFamilyZero() const throw(INTERP_KERNEL::Exception);
    void changeGroupName(const std::string& oldName, const std::string& newName) throw(INTERP_KERNEL::Exception);
    void changeFamilyName(const std::string& oldName, const std::string& newName) throw(INTERP_KERNEL::Exception);
    void changeFamilyId(int oldId, int newId) throw(INTERP_KERNEL::Exception);
    void changeAllGroupsContainingFamily(const std::string& familyNameToChange, const std::vector<std::string>& newFamiliesNames) throw(INTERP_KERNEL::Exception);
    void setFamilyInfo(const std::map<std::string,int>& info);
    void setGroupInfo(const std::map<std::string, std::vector<std::string> >&info);
    int getFamilyId(const std::string& name) const throw(INTERP_KERNEL::Exception);
    int getMaxAbsFamilyId() const throw(INTERP_KERNEL::Exception);
    int getMaxFamilyId() const throw(INTERP_KERNEL::Exception);
    int getMinFamilyId() const throw(INTERP_KERNEL::Exception);
    int getTheMaxAbsFamilyId() const throw(INTERP_KERNEL::Exception);
    int getTheMaxFamilyId() const throw(INTERP_KERNEL::Exception);
    int getTheMinFamilyId() const throw(INTERP_KERNEL::Exception);
    virtual int getMaxAbsFamilyIdInArrays() const throw(INTERP_KERNEL::Exception);
    virtual int getMaxFamilyIdInArrays() const throw(INTERP_KERNEL::Exception);
    virtual int getMinFamilyIdInArrays() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getAllFamiliesIdsReferenced() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *computeAllFamilyIdsInUse() const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamiliesIds(const std::vector<std::string>& famNames) const throw(INTERP_KERNEL::Exception);
    std::string getFamilyNameGivenId(int id) const throw(INTERP_KERNEL::Exception);
    bool ensureDifferentFamIdsPerLevel() throw(INTERP_KERNEL::Exception);
    void normalizeFamIdsTrio() throw(INTERP_KERNEL::Exception);
    void normalizeFamIdsMEDFile() throw(INTERP_KERNEL::Exception);
    virtual int getMeshDimension() const throw(INTERP_KERNEL::Exception);
    virtual std::string simpleRepr() const throw(INTERP_KERNEL::Exception);
    virtual std::string advancedRepr() const throw(INTERP_KERNEL::Exception);
    //
    virtual MEDCouplingMesh *getGenMeshAtLevel(int meshDimRelToMax, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual void setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr) throw(INTERP_KERNEL::Exception);
    virtual void setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr) throw(INTERP_KERNEL::Exception);
    virtual void setNameFieldAtLevel(int meshDimRelToMaxExt, DataArrayAsciiChar *nameArr) throw(INTERP_KERNEL::Exception);
    virtual void addNodeGroup(const DataArrayInt *ids) throw(INTERP_KERNEL::Exception);
    virtual void addGroup(int meshDimRelToMaxExt, const DataArrayInt *ids) throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getGroupsArr(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getGroupArr(int meshDimRelToMaxExt, const std::string& grp, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getFamilyArr(int meshDimRelToMaxExt, const std::string& fam, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeGroupArr(const std::string& grp, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeGroupsArr(const std::vector<std::string>& grps, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeFamilyArr(const std::string& fam, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeFamiliesArr(const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception);
    int getNumberOfJoints();
    MEDFileJoints *getJoints();
    void           setJoints( MEDFileJoints* joints );
    %extend
       {
         std::string __str__() const throw(INTERP_KERNEL::Exception)
         {
           return self->simpleRepr();
         }

         MEDCouplingMesh *__getitem__(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
         {
           return self->getGenMeshAtLevel(meshDimRelToMaxExt,false);
         }

         PyObject *getTime() throw(INTERP_KERNEL::Exception)
         {
           int tmp1,tmp2;
           double tmp0=self->getTime(tmp1,tmp2);
           PyObject *res = PyList_New(3);
           PyList_SetItem(res,0,SWIG_From_int(tmp1));
           PyList_SetItem(res,1,SWIG_From_int(tmp2));
           PyList_SetItem(res,2,SWIG_From_double(tmp0));
           return res;
         }

         virtual PyObject *isEqual(const MEDFileMesh *other, double eps) const throw(INTERP_KERNEL::Exception)
         {
           std::string what;
           bool ret0=self->isEqual(other,eps,what);
           PyObject *res=PyList_New(2);
           PyObject *ret0Py=ret0?Py_True:Py_False;
           Py_XINCREF(ret0Py);
           PyList_SetItem(res,0,ret0Py);
           PyList_SetItem(res,1,PyString_FromString(what.c_str()));
           return res;
         }

         void setGroupsAtLevel(int meshDimRelToMaxExt, PyObject *li, bool renum=false) throw(INTERP_KERNEL::Exception)
         {
           std::vector<const DataArrayInt *> grps;
           convertFromPyObjVectorOfObj<const ParaMEDMEM::DataArrayInt *>(li,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,"DataArrayInt",grps);
           self->setGroupsAtLevel(meshDimRelToMaxExt,grps,renum);
         }
         
         PyObject *areFamsEqual(const MEDFileMesh *other) const throw(INTERP_KERNEL::Exception)
         {
           std::string what;
           bool ret0=self->areFamsEqual(other,what);
           PyObject *res=PyList_New(2);
           PyObject *ret0Py=ret0?Py_True:Py_False;
           Py_XINCREF(ret0Py);
           PyList_SetItem(res,0,ret0Py);
           PyList_SetItem(res,1,PyString_FromString(what.c_str()));
           return res;
         }

         PyObject *areGrpsEqual(const MEDFileMesh *other) const throw(INTERP_KERNEL::Exception)
         {
           std::string what;
           bool ret0=self->areGrpsEqual(other,what);
           PyObject *res=PyList_New(2);
           PyObject *ret0Py=ret0?Py_True:Py_False;
           Py_XINCREF(ret0Py);
           PyList_SetItem(res,0,ret0Py);
           PyList_SetItem(res,1,PyString_FromString(what.c_str()));
           return res;
         }

         PyObject *getAllGeoTypes() const throw(INTERP_KERNEL::Exception)
         {
           std::vector<INTERP_KERNEL::NormalizedCellType> result(self->getAllGeoTypes());
           std::vector<INTERP_KERNEL::NormalizedCellType>::const_iterator iL=result.begin();
           PyObject *res=PyList_New(result.size());
           for(int i=0;iL!=result.end(); i++, iL++)
             PyList_SetItem(res,i,PyInt_FromLong(*iL));
           return res;
         }

         PyObject *getGeoTypesAtLevel(int meshDimRelToMax) const throw(INTERP_KERNEL::Exception)
         {
           std::vector<INTERP_KERNEL::NormalizedCellType> result(self->getGeoTypesAtLevel(meshDimRelToMax));
           std::vector<INTERP_KERNEL::NormalizedCellType>::const_iterator iL=result.begin();
           PyObject *res=PyList_New(result.size());
           for(int i=0;iL!=result.end(); i++, iL++)
             PyList_SetItem(res,i,PyInt_FromLong(*iL));
           return res;
         }

         PyObject *getFamilyFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
         {
           const DataArrayInt *tmp=self->getFamilyFieldAtLevel(meshDimRelToMaxExt);
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
         }

         PyObject *getOrCreateAndGetFamilyFieldAtLevel(int meshDimRelToMaxExt) throw(INTERP_KERNEL::Exception)
         {
           const DataArrayInt *tmp=self->getOrCreateAndGetFamilyFieldAtLevel(meshDimRelToMaxExt);
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
         }

         PyObject *getNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
         {
           const DataArrayInt *tmp=self->getNumberFieldAtLevel(meshDimRelToMaxExt);
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
         }

         PyObject *getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
         {
           const DataArrayInt *tmp=self->getRevNumberFieldAtLevel(meshDimRelToMaxExt);
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
         }
         
         PyObject *getNameFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
         {
           const DataArrayAsciiChar *tmp=self->getNameFieldAtLevel(meshDimRelToMaxExt);
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_ParaMEDMEM__DataArrayAsciiChar, SWIG_POINTER_OWN | 0 );
         }

         PyObject *findOrCreateAndGiveFamilyWithId(int id, bool& created) throw(INTERP_KERNEL::Exception)
         {
           bool ret1;
           std::string ret0=self->findOrCreateAndGiveFamilyWithId(id,ret1);
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,PyString_FromString(ret0.c_str()));
           PyTuple_SetItem(ret,1,SWIG_From_bool(ret1));
           return ret;
         }
         
         PyObject *unPolyze() throw(INTERP_KERNEL::Exception)
         {
           DataArrayInt *ret3=0;
           std::vector<int> ret1,ret2;
           bool ret0=self->unPolyze(ret1,ret2,ret3);
           PyObject *ret=PyTuple_New(4);
           PyTuple_SetItem(ret,0,SWIG_From_bool(ret0));
           //
           PyObject *retLev1_0=PyList_New((int)ret1.size()/3);
           for(int j=0;j<(int)ret1.size()/3;j++)
             {
               PyObject *retLev2=PyList_New(3);
               PyList_SetItem(retLev2,0,SWIG_From_int(ret1[3*j]));
               PyList_SetItem(retLev2,1,SWIG_From_int(ret1[3*j+1]));
               PyList_SetItem(retLev2,2,SWIG_From_int(ret1[3*j+2]));
               PyList_SetItem(retLev1_0,j,retLev2);
             }
           PyTuple_SetItem(ret,1,retLev1_0);
           //
           PyObject *retLev1_1=PyList_New((int)ret2.size()/3);
           for(int j=0;j<(int)ret2.size()/3;j++)
             {
               PyObject *retLev2=PyList_New(3);
               PyList_SetItem(retLev2,0,SWIG_From_int(ret2[3*j]));
               PyList_SetItem(retLev2,1,SWIG_From_int(ret2[3*j+1]));
               PyList_SetItem(retLev2,2,SWIG_From_int(ret2[3*j+2]));
               PyList_SetItem(retLev1_1,j,retLev2);
             }
           PyTuple_SetItem(ret,2,retLev1_1);
           //
           PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(ret3),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           return ret;
         }
       }
  };

  class MEDFileUMesh : public MEDFileMesh
  {
  public:
    static MEDFileUMesh *New(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    static MEDFileUMesh *New(const std::string& fileName, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    static MEDFileUMesh *New();
    ~MEDFileUMesh();
    int getSpaceDimension() const throw(INTERP_KERNEL::Exception);
    int getRelativeLevOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const throw(INTERP_KERNEL::Exception);
    //
    std::vector<int> getGrpNonEmptyLevels(const std::string& grp) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getGrpNonEmptyLevelsExt(const std::string& grp) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamNonEmptyLevels(const std::string& fam) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamNonEmptyLevelsExt(const std::string& fam) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getGrpsNonEmptyLevels(const std::vector<std::string>& grps) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getGrpsNonEmptyLevelsExt(const std::vector<std::string>& grps) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamsNonEmptyLevels(const std::vector<std::string>& fams) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamsNonEmptyLevelsExt(const std::vector<std::string>& fams) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getGroupsOnSpecifiedLev(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getGroup(int meshDimRelToMaxExt, const std::string& grp, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getGroups(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getFamily(int meshDimRelToMaxExt, const std::string& fam, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getFamilies(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getNodeGroupsArr(const std::vector<std::string>& grps, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getMeshAtLevel(int meshDimRelToMaxExt, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getLevel0Mesh(bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getLevelM1Mesh(bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getLevelM2Mesh(bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getLevelM3Mesh(bool renum=false) const throw(INTERP_KERNEL::Exception);
    void forceComputationOfParts() const throw(INTERP_KERNEL::Exception);
    //
    void setFamilyNameAttachedOnId(int id, const std::string& newFamName) throw(INTERP_KERNEL::Exception);
    void setCoords(DataArrayDouble *coords) throw(INTERP_KERNEL::Exception);
    void eraseGroupsAtLevel(int meshDimRelToMaxExt) throw(INTERP_KERNEL::Exception);
    void removeMeshAtLevel(int meshDimRelToMax) throw(INTERP_KERNEL::Exception);
    void setMeshAtLevel(int meshDimRelToMax, MEDCoupling1GTUMesh *m) throw(INTERP_KERNEL::Exception);
    void setMeshAtLevel(int meshDimRelToMax, MEDCouplingUMesh *m, bool newOrOld=false) throw(INTERP_KERNEL::Exception);
    void optimizeFamilies() throw(INTERP_KERNEL::Exception);
    DataArrayInt *zipCoords() throw(INTERP_KERNEL::Exception);
    DataArrayInt *extractFamilyFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *extractNumberFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const throw(INTERP_KERNEL::Exception);
    MEDFileUMesh *buildExtrudedMesh(const MEDCouplingUMesh *m1D, int policy) const throw(INTERP_KERNEL::Exception);
    MEDFileUMesh *linearToQuadratic(int conversionType=0, double eps=1e-12) const throw(INTERP_KERNEL::Exception);
    MEDFileUMesh *quadraticToLinear(double eps=1e-12) const throw(INTERP_KERNEL::Exception);
    %extend
       { 
         MEDFileUMesh(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileUMesh::New(fileName,mName,dt,it,mrs);
         }

         MEDFileUMesh(const std::string& fileName, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileUMesh::New(fileName,mrs);
         }

         MEDFileUMesh()
         {
           return MEDFileUMesh::New();
         }

         // serialization
         static PyObject *___new___(PyObject *cls, PyObject *args) throw(INTERP_KERNEL::Exception)
         {
           return NewMethWrapCallInitOnlyIfEmptyDictInInput(cls,args,"MEDFileUMesh");
         }

         static MEDFileUMesh *LoadPartOf(const std::string& fileName, const std::string& mName, PyObject *types, const std::vector<int>& slicPerTyp, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception)
         {
           std::vector<int> typesCpp1;
           convertPyToNewIntArr3(types,typesCpp1);
           std::size_t sz(typesCpp1.size());
           std::vector<INTERP_KERNEL::NormalizedCellType> typesCpp2(sz);
           for(std::size_t ii=0;ii<sz;ii++)
             typesCpp2[ii]=(INTERP_KERNEL::NormalizedCellType)typesCpp1[ii];
           return MEDFileUMesh::LoadPartOf(fileName,mName,typesCpp2,slicPerTyp,dt,it,mrs);
         }

         PyObject *__getnewargs__() throw(INTERP_KERNEL::Exception)
         {// put an empty dict in input to say to __new__ to call __init__...
           PyObject *ret(PyTuple_New(1));
           PyObject *ret0(PyDict_New());
           PyTuple_SetItem(ret,0,ret0);
           return ret;
         }

         PyObject *__getstate__() throw(INTERP_KERNEL::Exception)
         {
           std::vector<double> a0;
           std::vector<int> a1;
           std::vector<std::string> a2;
           std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > a3;
           MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> a4;
           self->serialize(a0,a1,a2,a3,a4);
           PyObject *ret(PyTuple_New(5));
           PyTuple_SetItem(ret,0,convertDblArrToPyList2(a0));
           PyTuple_SetItem(ret,1,convertIntArrToPyList2(a1));
           int sz(a2.size());
           PyObject *ret2(PyList_New(sz));
           for(int i=0;i<sz;i++)
             PyList_SetItem(ret2,i,PyString_FromString(a2[i].c_str()));
           PyTuple_SetItem(ret,2,ret2);
           sz=a3.size();
           PyObject *ret3(PyList_New(sz));
           for(int i=0;i<sz;i++)
             {
               DataArrayInt *elt(a3[i]);
               if(elt)
                 elt->incrRef();
               PyList_SetItem(ret3,i,SWIG_NewPointerObj(SWIG_as_voidptr(elt),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             }
           PyTuple_SetItem(ret,3,ret3);
           DataArrayDouble *ret4(a4);
           if(ret4)
             ret4->incrRef();
           PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(ret4),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
           return ret;
         }

         void __setstate__(PyObject *inp) throw(INTERP_KERNEL::Exception)
         {
           static const char MSG[]="MEDFileUMesh.__setstate__ : expected input is a tuple of size 4 !";
           if(!PyTuple_Check(inp))
             throw INTERP_KERNEL::Exception(MSG);
           int sz(PyTuple_Size(inp));
           if(sz!=5)
             throw INTERP_KERNEL::Exception(MSG);
           std::vector<double> a0;
           std::vector<int> a1;
           std::vector<std::string> a2;
           std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > a3;
           MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> a4;
           //
           PyObject *a0py(PyTuple_GetItem(inp,0)),*a1py(PyTuple_GetItem(inp,1)),*a2py(PyTuple_GetItem(inp,2));
           int tmp(-1);
           fillArrayWithPyListDbl3(a0py,tmp,a0);
           convertPyToNewIntArr3(a1py,a1);
           fillStringVector(a2py,a2);
           //
           PyObject *b0py(PyTuple_GetItem(inp,3)),*b1py(PyTuple_GetItem(inp,4));
           void *argp(0);
           int status(SWIG_ConvertPtr(b1py,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble,0|0));
           if(!SWIG_IsOK(status))
             throw INTERP_KERNEL::Exception(MSG);
           a4=reinterpret_cast<DataArrayDouble *>(argp);
           if((DataArrayDouble *)a4)
             a4->incrRef();
           {
             std::vector< DataArrayInt * > a3Tmp;
             convertFromPyObjVectorOfObj<ParaMEDMEM::DataArrayInt *>(b0py,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,"DataArrayInt",a3Tmp);
             std::size_t sz(a3Tmp.size());
             a3.resize(sz);
             for(std::size_t i=0;i<sz;i++)
               {
                 a3[i]=a3Tmp[i];
                 if(a3Tmp[i])
                   a3Tmp[i]->incrRef();
               }
             self->unserialize(a0,a1,a2,a3,a4);
           }
         }

         void __setitem__(int meshDimRelToMax, MEDCouplingPointSet *mesh) throw(INTERP_KERNEL::Exception)
         {
           if(!mesh)
             throw INTERP_KERNEL::Exception("MEDFileUMesh::__setitem__ : Input mesh is NULL !");
           MEDCouplingUMesh *m0(dynamic_cast<MEDCouplingUMesh *>(mesh));
           if(m0)
             {
               self->setMeshAtLevel(meshDimRelToMax,m0,false);
               return ;
             }
           MEDCoupling1GTUMesh *m1(dynamic_cast<MEDCoupling1GTUMesh *>(mesh));
           if(m1)
             {
               self->setMeshAtLevel(meshDimRelToMax,m1);
               return ;
             }
           throw INTERP_KERNEL::Exception("MEDFileUMesh::__setitem__ : Not recognized input mesh !");
         }

         void __delitem__(int meshDimRelToMax) throw(INTERP_KERNEL::Exception)
         {
           self->removeMeshAtLevel(meshDimRelToMax);
         }

         void setMeshes(PyObject *li, bool renum=false) throw(INTERP_KERNEL::Exception)
         {
           std::vector<const MEDCouplingUMesh *> ms;
           convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingUMesh *>(li,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,"MEDCouplingUMesh",ms);
           self->setMeshes(ms,renum);
         }

         void setGroupsFromScratch(int meshDimRelToMax, PyObject *li, bool renum=false) throw(INTERP_KERNEL::Exception)
         {
           std::vector<const MEDCouplingUMesh *> ms;
           convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingUMesh *>(li,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,"MEDCouplingUMesh",ms);
           self->setGroupsFromScratch(meshDimRelToMax,ms,renum);
         }
         
         void setGroupsOnSetMesh(int meshDimRelToMax, PyObject *li, bool renum=false) throw(INTERP_KERNEL::Exception)
         {
           std::vector<const MEDCouplingUMesh *> ms;
           convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingUMesh *>(li,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,"MEDCouplingUMesh",ms);
           self->setGroupsOnSetMesh(meshDimRelToMax,ms,renum);
         }

         DataArrayDouble *getCoords() const throw(INTERP_KERNEL::Exception)
         {
           DataArrayDouble *ret=self->getCoords();
           if(ret)
             ret->incrRef();
           return ret;
         }

         PartDefinition *getPartDefAtLevel(int meshDimRelToMaxExt, INTERP_KERNEL::NormalizedCellType gt=INTERP_KERNEL::NORM_ERROR) const throw(INTERP_KERNEL::Exception)
         {
           const PartDefinition *ret(self->getPartDefAtLevel(meshDimRelToMaxExt,gt));
           if(ret)
             ret->incrRef();
           return const_cast<PartDefinition *>(ret);
         }

         PyObject *buildInnerBoundaryAlongM1Group(const std::string& grpNameM1) throw(INTERP_KERNEL::Exception)
         {
           DataArrayInt *ret0=0,*ret1=0,*ret2=0;
           self->buildInnerBoundaryAlongM1Group(grpNameM1,ret0,ret1,ret2);
           PyObject *ret=PyTuple_New(3);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(ret2),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           return ret;
         }
         
         MEDCoupling1GTUMesh *getDirectUndergroundSingleGeoTypeMesh(INTERP_KERNEL::NormalizedCellType gt) const throw(INTERP_KERNEL::Exception)
         {
           MEDCoupling1GTUMesh *ret(self->getDirectUndergroundSingleGeoTypeMesh(gt));
           if(ret)
             ret->incrRef();
           return ret;
         }

         PyObject *getDirectUndergroundSingleGeoTypeMeshes(int meshDimRelToMax) const throw(INTERP_KERNEL::Exception)
         {
           std::vector<MEDCoupling1GTUMesh *> tmp(self->getDirectUndergroundSingleGeoTypeMeshes(meshDimRelToMax));
           std::size_t sz(tmp.size());
           PyObject *ret=PyList_New(sz);
           for(std::size_t i=0;i<sz;i++)
             {
               if(tmp[i])
                 tmp[i]->incrRef();
               PyList_SetItem(ret,i,convertMesh(tmp[i], SWIG_POINTER_OWN | 0 ));
             }
           return ret;
         }
       }
  };

  class MEDFileStructuredMesh : public MEDFileMesh
  {
  public:
    %extend
    {
      MEDCoupling1SGTUMesh *getImplicitFaceMesh() const throw(INTERP_KERNEL::Exception)
      {
        MEDCoupling1SGTUMesh *ret(self->getImplicitFaceMesh());
        if(ret)
          ret->incrRef();
        return ret;
      }
    }
  };

  class MEDFileCMesh : public MEDFileStructuredMesh
  {
  public:
    static MEDFileCMesh *New();
    static MEDFileCMesh *New(const std::string& fileName, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    static MEDFileCMesh *New(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    void setMesh(MEDCouplingCMesh *m) throw(INTERP_KERNEL::Exception);
    int getSpaceDimension() const throw(INTERP_KERNEL::Exception);
    %extend
       {
         MEDFileCMesh()
         {
           return MEDFileCMesh::New();
         }

         MEDFileCMesh(const std::string& fileName, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileCMesh::New(fileName,mrs);
         }

         MEDFileCMesh(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileCMesh::New(fileName,mName,dt,it,mrs);
         }
         
         PyObject *getMesh() const throw(INTERP_KERNEL::Exception)
         {
           const MEDCouplingCMesh *tmp=self->getMesh();
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_ParaMEDMEM__MEDCouplingCMesh, SWIG_POINTER_OWN | 0 );
         }
       }
  };

  class MEDFileCurveLinearMesh : public MEDFileStructuredMesh
  {
  public:
    static MEDFileCurveLinearMesh *New();
    static MEDFileCurveLinearMesh *New(const std::string& fileName, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    static MEDFileCurveLinearMesh *New(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    void setMesh(MEDCouplingCurveLinearMesh *m) throw(INTERP_KERNEL::Exception);
    %extend
       {
         MEDFileCurveLinearMesh()
         {
           return MEDFileCurveLinearMesh::New();
         }

         MEDFileCurveLinearMesh(const std::string& fileName, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileCurveLinearMesh::New(fileName,mrs);
         }

         MEDFileCurveLinearMesh(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileCurveLinearMesh::New(fileName,mName,dt,it,mrs);
         }
         
         PyObject *getMesh() const throw(INTERP_KERNEL::Exception)
         {
           const MEDCouplingCurveLinearMesh *tmp=self->getMesh();
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_ParaMEDMEM__MEDCouplingCurveLinearMesh, SWIG_POINTER_OWN | 0 );
         }
       }
  };

  class MEDFileMeshMultiTS : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileMeshMultiTS *New();
    static MEDFileMeshMultiTS *New(const std::string& fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileMeshMultiTS *New(const std::string& fileName, const std::string& mName) throw(INTERP_KERNEL::Exception);
    MEDFileMeshMultiTS *deepCpy() const throw(INTERP_KERNEL::Exception);
    std::string getName() const throw(INTERP_KERNEL::Exception);
    void write(const std::string& fileName, int mode) const throw(INTERP_KERNEL::Exception);
    void setOneTimeStep(MEDFileMesh *mesh1TimeStep) throw(INTERP_KERNEL::Exception);
    %extend
       { 
         MEDFileMeshMultiTS()
         {
           return MEDFileMeshMultiTS::New();
         }

         MEDFileMeshMultiTS(const std::string& fileName) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileMeshMultiTS::New(fileName);
         }

         MEDFileMeshMultiTS(const std::string& fileName, const std::string& mName) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileMeshMultiTS::New(fileName,mName);
         }

         MEDFileMesh *getOneTimeStep() const throw(INTERP_KERNEL::Exception)
           {
             MEDFileMesh *ret=self->getOneTimeStep();
             if(ret)
               ret->incrRef();
             return ret;
           }
       }
  };

  class MEDFileMeshesIterator
  {
  public:
    %extend
    {
      PyObject *next() throw(INTERP_KERNEL::Exception)
      {
        MEDFileMesh *ret=self->nextt();
        if(ret)
          {
            ret->incrRef();
            return convertMEDFileMesh(ret,SWIG_POINTER_OWN | 0 );
          }
        else
          {
            PyErr_SetString(PyExc_StopIteration,"No more data.");
            return 0;
          }
      }
    }
  };

  class MEDFileMeshes : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileMeshes *New();
    static MEDFileMeshes *New(const std::string& fileName) throw(INTERP_KERNEL::Exception);
    MEDFileMeshes *deepCpy() const throw(INTERP_KERNEL::Exception);
    void write(const std::string& fileName, int mode) const throw(INTERP_KERNEL::Exception);
    int getNumberOfMeshes() const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getMeshesNames() const throw(INTERP_KERNEL::Exception);
    //
    void resize(int newSize) throw(INTERP_KERNEL::Exception);
    void pushMesh(MEDFileMesh *mesh) throw(INTERP_KERNEL::Exception);
    void setMeshAtPos(int i, MEDFileMesh *mesh) throw(INTERP_KERNEL::Exception);
    void destroyMeshAtPos(int i) throw(INTERP_KERNEL::Exception);
    %extend
       {
         MEDFileMeshes()
         {
           return MEDFileMeshes::New();
         }

         MEDFileMeshes(const std::string& fileName) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileMeshes::New(fileName);
         }

         std::string __str__() const throw(INTERP_KERNEL::Exception)
           {
             return self->simpleRepr();
           }

         MEDFileMesh *__getitem__(PyObject *obj) throw(INTERP_KERNEL::Exception)
         {
           if(PyInt_Check(obj))
             {
               MEDFileMesh *ret=self->getMeshAtPos(InterpreteNegativeInt((int)PyInt_AS_LONG(obj),self->getNumberOfMeshes()));
               if(ret)
                 ret->incrRef();
               return ret;
             }
           else if(PyString_Check(obj))
             {
               MEDFileMesh *ret=self->getMeshWithName(PyString_AsString(obj));
               if(ret)
                 ret->incrRef();
               return ret;
             }
           else
             throw INTERP_KERNEL::Exception("MEDFileMeshes::__getitem__ : only integer or string with meshname supported !");
         }

         MEDFileMeshes *__setitem__(int obj, MEDFileMesh *mesh) throw(INTERP_KERNEL::Exception)
         {
           self->setMeshAtPos(obj,mesh);
           return self;
         }

         MEDFileMeshesIterator *__iter__() throw(INTERP_KERNEL::Exception)
         {
           return self->iterator();
         }

         int __len__() const throw(INTERP_KERNEL::Exception)
         {
           return self->getNumberOfMeshes();
         }
         
         MEDFileMesh *getMeshAtPos(int i) const throw(INTERP_KERNEL::Exception)
           {
             MEDFileMesh *ret=self->getMeshAtPos(i);
             if(ret)
               ret->incrRef();
             return ret;
           }
         MEDFileMesh *getMeshWithName(const std::string& mname) const throw(INTERP_KERNEL::Exception)
           {
             MEDFileMesh *ret=self->getMeshWithName(mname);
             if(ret)
               ret->incrRef();
             return ret;
           }
       }
  };

  class MEDFileFieldLoc : public RefCountObject
  {
  public:
    std::string getName() const;
    int getDimension() const;
    int getNumberOfGaussPoints() const;
    int getNumberOfPointsInCells() const;
    const std::vector<double>& getRefCoords() const;
    const std::vector<double>& getGaussCoords() const;
    const std::vector<double>& getGaussWeights() const;
    bool isEqual(const MEDFileFieldLoc& other, double eps) const throw(INTERP_KERNEL::Exception);
  %extend
    {
      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->repr();
      }
    }
  };

  class MEDFileFieldGlobsReal
  {
  public:
    void resetContent();
    void shallowCpyGlobs(const MEDFileFieldGlobsReal& other) throw(INTERP_KERNEL::Exception);
    void deepCpyGlobs(const MEDFileFieldGlobsReal& other) throw(INTERP_KERNEL::Exception);
    void shallowCpyOnlyUsedGlobs(const MEDFileFieldGlobsReal& other) throw(INTERP_KERNEL::Exception);
    void deepCpyOnlyUsedGlobs(const MEDFileFieldGlobsReal& other) throw(INTERP_KERNEL::Exception);
    void appendGlobs(const MEDFileFieldGlobsReal& other, double eps) throw(INTERP_KERNEL::Exception);
    void checkGlobsCoherency() const throw(INTERP_KERNEL::Exception);
    void checkGlobsPflsPartCoherency() const throw(INTERP_KERNEL::Exception);
    void checkGlobsLocsPartCoherency() const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getPfls() const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getLocs() const throw(INTERP_KERNEL::Exception);
    bool existsPfl(const std::string& pflName) const throw(INTERP_KERNEL::Exception);
    bool existsLoc(const std::string& locName) const throw(INTERP_KERNEL::Exception);
    std::string createNewNameOfPfl() const throw(INTERP_KERNEL::Exception);
    std::string createNewNameOfLoc() const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<int> > whichAreEqualProfiles() const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<int> > whichAreEqualLocs(double eps) const throw(INTERP_KERNEL::Exception);
    virtual std::vector<std::string> getPflsReallyUsed() const throw(INTERP_KERNEL::Exception);
    virtual std::vector<std::string> getLocsReallyUsed() const throw(INTERP_KERNEL::Exception);
    virtual std::vector<std::string> getPflsReallyUsedMulti() const throw(INTERP_KERNEL::Exception);
    virtual std::vector<std::string> getLocsReallyUsedMulti() const throw(INTERP_KERNEL::Exception);
    void killProfileIds(const std::vector<int>& pflIds) throw(INTERP_KERNEL::Exception);
    void killLocalizationIds(const std::vector<int>& locIds) throw(INTERP_KERNEL::Exception);
    void changePflName(const std::string& oldName, const std::string& newName) throw(INTERP_KERNEL::Exception);
    void changeLocName(const std::string& oldName, const std::string& newName) throw(INTERP_KERNEL::Exception);
    int getNbOfGaussPtPerCell(int locId) const throw(INTERP_KERNEL::Exception);
    int getLocalizationId(const std::string& loc) const throw(INTERP_KERNEL::Exception);
  %extend
     {
       PyObject *getProfile(const std::string& pflName) const throw(INTERP_KERNEL::Exception)
       {
         const DataArrayInt *ret=self->getProfile(pflName);
         if(ret)
           ret->incrRef();
         return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
       }

       PyObject *getProfileFromId(int pflId) const throw(INTERP_KERNEL::Exception)
       {
         const DataArrayInt *ret=self->getProfileFromId(pflId);
         if(ret)
           ret->incrRef();
         return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
       }

       PyObject *getLocalizationFromId(int locId) const throw(INTERP_KERNEL::Exception)
       {
         const MEDFileFieldLoc *loc=&self->getLocalizationFromId(locId);
         if(loc)
           loc->incrRef();
         return SWIG_NewPointerObj(SWIG_as_voidptr(loc),SWIGTYPE_p_ParaMEDMEM__MEDFileFieldLoc, SWIG_POINTER_OWN | 0 );
       }
       
       PyObject *getLocalization(const std::string& locName) const throw(INTERP_KERNEL::Exception)
       {
         const MEDFileFieldLoc *loc=&self->getLocalization(locName);
         if(loc)
           loc->incrRef();
         return SWIG_NewPointerObj(SWIG_as_voidptr(loc),SWIGTYPE_p_ParaMEDMEM__MEDFileFieldLoc, SWIG_POINTER_OWN | 0 );
       }
       
       PyObject *zipPflsNames() throw(INTERP_KERNEL::Exception)
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > ret=self->zipPflsNames();
         return convertVecPairVecStToPy(ret);
       }

       PyObject *zipLocsNames(double eps) throw(INTERP_KERNEL::Exception)
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > ret=self->zipLocsNames(eps);
         return convertVecPairVecStToPy(ret);
       }

       void changePflsNames(PyObject *li) throw(INTERP_KERNEL::Exception)
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > v=convertVecPairVecStFromPy(li);
         self->changePflsNames(v);
       }

       void changePflsRefsNamesGen(PyObject *li) throw(INTERP_KERNEL::Exception)
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > v=convertVecPairVecStFromPy(li);
         self->changePflsRefsNamesGen(v);
       }

       void changePflsNamesInStruct(PyObject *li) throw(INTERP_KERNEL::Exception)
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > v=convertVecPairVecStFromPy(li);
         self->changePflsNamesInStruct(v);
       }

       void changeLocsNames(PyObject *li) throw(INTERP_KERNEL::Exception)
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > v=convertVecPairVecStFromPy(li);
         self->changeLocsNames(v);
       }

       void changeLocsRefsNamesGen(PyObject *li) throw(INTERP_KERNEL::Exception)
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > v=convertVecPairVecStFromPy(li);
         self->changeLocsRefsNamesGen(v);
       }
       
       void changeLocsNamesInStruct(PyObject *li) throw(INTERP_KERNEL::Exception)
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > v=convertVecPairVecStFromPy(li);
         self->changeLocsNamesInStruct(v);
       }

       std::string simpleReprGlobs() const throw(INTERP_KERNEL::Exception)
       {
         std::ostringstream oss;
         self->simpleReprGlobs(oss);
         return oss.str();
       }
     }
  };

  class MEDFileAnyTypeField1TS : public RefCountObject, public MEDFileFieldGlobsReal, public MEDFileWritable
  {
  public:
    static MEDFileAnyTypeField1TS *New(const std::string& fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileAnyTypeField1TS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileAnyTypeField1TS *New(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    void write(const std::string& fileName, int mode) const throw(INTERP_KERNEL::Exception);
    void loadArrays() throw(INTERP_KERNEL::Exception);
    void loadArraysIfNecessary() throw(INTERP_KERNEL::Exception);
    void unloadArrays() throw(INTERP_KERNEL::Exception);
    void unloadArraysWithoutDataLoss() throw(INTERP_KERNEL::Exception);
    int getDimension() const throw(INTERP_KERNEL::Exception);
    int getIteration() const throw(INTERP_KERNEL::Exception);
    int getOrder() const throw(INTERP_KERNEL::Exception);
    std::string getName() throw(INTERP_KERNEL::Exception);
    void setName(const std::string& name) throw(INTERP_KERNEL::Exception);
    std::string getMeshName() throw(INTERP_KERNEL::Exception);
    void setMeshName(const std::string& newMeshName) throw(INTERP_KERNEL::Exception);
    int getMeshIteration() const throw(INTERP_KERNEL::Exception);
    int getMeshOrder() const throw(INTERP_KERNEL::Exception);
    int getNumberOfComponents() const throw(INTERP_KERNEL::Exception);
    bool isDealingTS(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    void setInfo(const std::vector<std::string>& infos) throw(INTERP_KERNEL::Exception);
    const std::vector<std::string>& getInfo() const throw(INTERP_KERNEL::Exception);
    bool presenceOfMultiDiscPerGeoType() const throw(INTERP_KERNEL::Exception);
    void setTime(int iteration, int order, double val) throw(INTERP_KERNEL::Exception);
    virtual MEDFileAnyTypeField1TS *shallowCpy() const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TS *deepCpy() const throw(INTERP_KERNEL::Exception);
    std::string getDtUnit() const throw(INTERP_KERNEL::Exception);
    void setDtUnit(const std::string& dtUnit) throw(INTERP_KERNEL::Exception);
    %extend
    {
      PyObject *getTime() throw(INTERP_KERNEL::Exception)
      {
        int tmp1,tmp2;
        double tmp0=self->getTime(tmp1,tmp2);
        PyObject *res = PyList_New(3);
        PyList_SetItem(res,0,SWIG_From_int(tmp1));
        PyList_SetItem(res,1,SWIG_From_int(tmp2));
        PyList_SetItem(res,2,SWIG_From_double(tmp0));
        return res;
      }

      PyObject *getDtIt() const throw(INTERP_KERNEL::Exception)
      {
        std::pair<int,int> res=self->getDtIt();
        PyObject *elt=PyTuple_New(2);
        PyTuple_SetItem(elt,0,SWIG_From_int(res.first));
        PyTuple_SetItem(elt,1,SWIG_From_int(res.second));
        return elt;
      }

      void setProfileNameOnLeaf(INTERP_KERNEL::NormalizedCellType typ, int locId, const std::string& newPflName, bool forceRenameOnGlob=false) throw(INTERP_KERNEL::Exception)
      {
        self->setProfileNameOnLeaf(0,typ,locId,newPflName,forceRenameOnGlob);
      }
      
      void setLocNameOnLeaf(INTERP_KERNEL::NormalizedCellType typ, int locId, const std::string& newLocName, bool forceRenameOnGlob=false) throw(INTERP_KERNEL::Exception)
      {
        self->setLocNameOnLeaf(0,typ,locId,newLocName,forceRenameOnGlob);
      }

      bool changeMeshNames(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<std::string,std::string> > modifTab=convertVecPairStStFromPy(li);
        return self->changeMeshNames(modifTab);
      }
      
      PyObject *getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception)
      {
        std::vector<TypeOfField> ret=self->getTypesOfFieldAvailable();
        PyObject *ret2=PyList_New(ret.size());
        for(int i=0;i<(int)ret.size();i++)
          PyList_SetItem(ret2,i,SWIG_From_int(ret[i]));
        return ret2;
      }

      PyObject *getNonEmptyLevels(const std::string& mname=std::string()) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> ret1;
        int ret0=self->getNonEmptyLevels(mname,ret1);
        PyObject *elt=PyTuple_New(2);
        PyTuple_SetItem(elt,0,SWIG_From_int(ret0));
        PyTuple_SetItem(elt,1,convertIntArrToPyList2(ret1));
        return elt;
      }

      PyObject *getFieldSplitedByType(const std::string& mname=std::string()) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<INTERP_KERNEL::NormalizedCellType> types;
        std::vector< std::vector<TypeOfField> > typesF;
        std::vector< std::vector<std::string> > pfls;
        std::vector< std::vector<std::string> > locs;
        std::vector< std::vector< std::pair<int,int> > > ret=self->getFieldSplitedByType(mname,types,typesF,pfls,locs);
        int sz=ret.size();
        PyObject *ret2=PyList_New(sz);
           for(int i=0;i<sz;i++)
             {
               const std::vector< std::pair<int,int> >& dadsI=ret[i];
               const std::vector<TypeOfField>& typesFI=typesF[i];
               const std::vector<std::string>& pflsI=pfls[i];
               const std::vector<std::string>& locsI=locs[i];
               PyObject *elt=PyTuple_New(2);
               PyTuple_SetItem(elt,0,SWIG_From_int(types[i]));
               int sz2=ret[i].size();
               PyObject *elt2=PyList_New(sz2);
               for(int j=0;j<sz2;j++)
                 {
                   PyObject *elt3=PyTuple_New(4);
                   PyTuple_SetItem(elt3,0,SWIG_From_int(typesFI[j]));
                   PyObject *elt4=PyTuple_New(2); PyTuple_SetItem(elt4,0,SWIG_From_int(dadsI[j].first)); PyTuple_SetItem(elt4,1,SWIG_From_int(dadsI[j].second));
                   PyTuple_SetItem(elt3,1,elt4);
                   PyTuple_SetItem(elt3,2,PyString_FromString(pflsI[j].c_str()));
                   PyTuple_SetItem(elt3,3,PyString_FromString(locsI[j].c_str()));
                   PyList_SetItem(elt2,j,elt3);
                 }
               PyTuple_SetItem(elt,1,elt2);
               PyList_SetItem(ret2,i,elt);
             }
           return ret2;
      }

      PyObject *splitComponents() const throw(INTERP_KERNEL::Exception)
      {
        std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TS > > ret=self->splitComponents();
        std::size_t sz=ret.size();
        PyObject *retPy=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(retPy,i,convertMEDFileField1TS(ret[i].retn(), SWIG_POINTER_OWN | 0 ));
        return retPy;
      }

      PyObject *splitDiscretizations() const throw(INTERP_KERNEL::Exception)
      {
        std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TS > > ret=self->splitDiscretizations();
        std::size_t sz=ret.size();
        PyObject *retPy=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(retPy,i,convertMEDFileField1TS(ret[i].retn(), SWIG_POINTER_OWN | 0 ));
        return retPy;
      }

      PyObject *splitMultiDiscrPerGeoTypes() const throw(INTERP_KERNEL::Exception)
      {
        std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TS > > ret=self->splitMultiDiscrPerGeoTypes();
        std::size_t sz=ret.size();
        PyObject *retPy=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(retPy,i,convertMEDFileField1TS(ret[i].retn(), SWIG_POINTER_OWN | 0 ));
        return retPy;
      }
    }
  };

  class MEDFileField1TS : public MEDFileAnyTypeField1TS
  {
  public:
    static MEDFileField1TS *New(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileField1TS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileField1TS *New(const std::string& fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileField1TS *New();
    ParaMEDMEM::MEDFileIntField1TS *convertToInt(bool isDeepCpyGlobs=true) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtLevelOld(TypeOfField type, const std::string& mname, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    //
    void setFieldNoProfileSBT(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception);
    void setFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception);
    void setProfileNameOnLeaf(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const std::string& newPflName, bool forceRenameOnGlob=false) throw(INTERP_KERNEL::Exception);
    void setLocNameOnLeaf(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const std::string& newLocName, bool forceRenameOnGlob=false) throw(INTERP_KERNEL::Exception);
    %extend
       {
         MEDFileField1TS(const std::string& fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileField1TS::New(fileName,loadAll);
         }
         
         MEDFileField1TS(const std::string& fileName, const std::string& fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileField1TS::New(fileName,fieldName,loadAll);
         }

         MEDFileField1TS(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileField1TS::New(fileName,fieldName,iteration,order,loadAll);
         }

         MEDFileField1TS()
         {
           return MEDFileField1TS::New();
         }

         void copyTinyInfoFrom(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception)
         {
           const DataArrayDouble *arr=0;
           if(field)
             arr=field->getArray();
           self->copyTinyInfoFrom(field,arr);
         }
         
         std::string __str__() const throw(INTERP_KERNEL::Exception)
         {
           return self->simpleRepr();
         }
         
         PyObject *getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh) const throw(INTERP_KERNEL::Exception)
         {
           DataArrayInt *ret1=0;
           DataArrayDouble *ret0=self->getFieldWithProfile(type,meshDimRelToMax,mesh,ret1);
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           return ret;
         }

         PyObject *getFieldSplitedByType2(const std::string& mname=std::string()) const throw(INTERP_KERNEL::Exception)
         {
           std::vector<INTERP_KERNEL::NormalizedCellType> types;
           std::vector< std::vector<TypeOfField> > typesF;
           std::vector< std::vector<std::string> > pfls;
           std::vector< std::vector<std::string> > locs;
           std::vector< std::vector<DataArrayDouble *> > ret=self->getFieldSplitedByType2(mname,types,typesF,pfls,locs);
           int sz=ret.size();
           PyObject *ret2=PyList_New(sz);
           for(int i=0;i<sz;i++)
             {
               const std::vector<DataArrayDouble *>& dadsI=ret[i];
               const std::vector<TypeOfField>& typesFI=typesF[i];
               const std::vector<std::string>& pflsI=pfls[i];
               const std::vector<std::string>& locsI=locs[i];
               PyObject *elt=PyTuple_New(2);
               PyTuple_SetItem(elt,0,SWIG_From_int(types[i]));
               int sz2=ret[i].size();
               PyObject *elt2=PyList_New(sz2);
               for(int j=0;j<sz2;j++)
                 {
                   PyObject *elt3=PyTuple_New(4);
                   PyTuple_SetItem(elt3,0,SWIG_From_int(typesFI[j]));
                   PyTuple_SetItem(elt3,1,SWIG_NewPointerObj(SWIG_as_voidptr(dadsI[j]),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
                   PyTuple_SetItem(elt3,2,PyString_FromString(pflsI[j].c_str()));
                   PyTuple_SetItem(elt3,3,PyString_FromString(locsI[j].c_str()));
                   PyList_SetItem(elt2,j,elt3);
                 }
               PyTuple_SetItem(elt,1,elt2);
               PyList_SetItem(ret2,i,elt);
             }
           return ret2;
         }

         DataArrayDouble *getUndergroundDataArray() const throw(INTERP_KERNEL::Exception)
         {
           DataArrayDouble *ret=self->getUndergroundDataArray();
           if(ret)
             ret->incrRef();
           return ret;
         }

         PyObject *getUndergroundDataArrayExt() const throw(INTERP_KERNEL::Exception)
         {
           std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > > elt1Cpp;
           DataArrayDouble *elt0=self->getUndergroundDataArrayExt(elt1Cpp);
           if(elt0)
             elt0->incrRef();
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(elt0),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
           std::size_t sz=elt1Cpp.size();
           PyObject *elt=PyList_New(sz);
           for(std::size_t i=0;i<sz;i++)
             {
               PyObject *elt1=PyTuple_New(2);
               PyObject *elt2=PyTuple_New(2);
               PyTuple_SetItem(elt2,0,SWIG_From_int((int)elt1Cpp[i].first.first));
               PyTuple_SetItem(elt2,1,SWIG_From_int(elt1Cpp[i].first.second));
               PyObject *elt3=PyTuple_New(2);
               PyTuple_SetItem(elt3,0,SWIG_From_int(elt1Cpp[i].second.first));
               PyTuple_SetItem(elt3,1,SWIG_From_int(elt1Cpp[i].second.second));
               PyTuple_SetItem(elt1,0,elt2);
               PyTuple_SetItem(elt1,1,elt3);
               PyList_SetItem(elt,i,elt1);
             }
           PyTuple_SetItem(ret,1,elt);
           return ret;
         }
       }
  };

  class MEDFileIntField1TS : public MEDFileAnyTypeField1TS
  {
  public:
    static MEDFileIntField1TS *New();
    static MEDFileIntField1TS *New(const std::string& fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileIntField1TS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileIntField1TS *New(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    ParaMEDMEM::MEDFileField1TS *convertToDouble(bool isDeepCpyGlobs=true) const throw(INTERP_KERNEL::Exception);
    //
    void setFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals) throw(INTERP_KERNEL::Exception);
    void setFieldProfile(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception);
    %extend
    {
      MEDFileIntField1TS() throw(INTERP_KERNEL::Exception)
      {
        return MEDFileIntField1TS::New();
      }

      MEDFileIntField1TS(const std::string& fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception)
      {
        return MEDFileIntField1TS::New(fileName,loadAll);
      }

      MEDFileIntField1TS(const std::string& fileName, const std::string& fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception)
      {
        return MEDFileIntField1TS::New(fileName,fieldName,loadAll);
      }

      MEDFileIntField1TS(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true) throw(INTERP_KERNEL::Exception)
      {
        return MEDFileIntField1TS::New(fileName,fieldName,iteration,order,loadAll);
      }

      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->simpleRepr();
      }

      PyObject *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0;
        MEDCouplingFieldDouble *ret0=self->getFieldAtLevel(type,meshDimRelToMax,ret1,renumPol);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *getFieldAtTopLevel(TypeOfField type, int renumPol=0) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0;
        MEDCouplingFieldDouble *ret0=self->getFieldAtTopLevel(type,ret1,renumPol);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0;
        MEDCouplingFieldDouble *ret0=self->getFieldOnMeshAtLevel(type,meshDimRelToMax,mesh,ret1,renumPol);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }
      
      PyObject *getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0;
        MEDCouplingFieldDouble *ret0=self->getFieldOnMeshAtLevel(type,mesh,ret1,renumPol);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }
      
      PyObject *getFieldAtLevelOld(TypeOfField type, const std::string& mname, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0;
        MEDCouplingFieldDouble *ret0=self->getFieldAtLevelOld(type,mname,meshDimRelToMax,ret1,renumPol);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh) const throw(INTERP_KERNEL::Exception)
      {
         DataArrayInt *ret1=0;
         DataArrayInt *ret0=self->getFieldWithProfile(type,meshDimRelToMax,mesh,ret1);
         PyObject *ret=PyTuple_New(2);
         PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
         PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
         return ret;
      }
      
      DataArrayInt *getUndergroundDataArray() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret=self->getUndergroundDataArray();
        if(ret)
          ret->incrRef();
        return ret;
      }
    }
  };

  class MEDFileAnyTypeFieldMultiTSIterator
  {
  public:
    %extend
    {
      PyObject *next() throw(INTERP_KERNEL::Exception)
      {
        MEDFileAnyTypeField1TS *ret=self->nextt();
        if(ret)
          return convertMEDFileField1TS(ret, SWIG_POINTER_OWN | 0 );
        else
          {
            PyErr_SetString(PyExc_StopIteration,"No more data.");
            return 0;
          }
      }
    }
  };

  class MEDFileAnyTypeFieldMultiTS : public RefCountObject, public MEDFileFieldGlobsReal, public MEDFileWritable
  {
  public:
    static MEDFileAnyTypeFieldMultiTS *New(const std::string& fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileAnyTypeFieldMultiTS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTS *deepCpy() const throw(INTERP_KERNEL::Exception);
    virtual MEDFileAnyTypeFieldMultiTS *shallowCpy() const throw(INTERP_KERNEL::Exception);
    std::string getName() const throw(INTERP_KERNEL::Exception);
    void setName(const std::string& name) throw(INTERP_KERNEL::Exception);
    std::string getDtUnit() const throw(INTERP_KERNEL::Exception);
    void setDtUnit(const std::string& dtUnit) throw(INTERP_KERNEL::Exception);
    std::string getMeshName() const throw(INTERP_KERNEL::Exception);
    void setMeshName(const std::string& newMeshName) throw(INTERP_KERNEL::Exception);
    const std::vector<std::string>& getInfo() const throw(INTERP_KERNEL::Exception);
    bool presenceOfMultiDiscPerGeoType() const throw(INTERP_KERNEL::Exception);
    int getNumberOfComponents() const throw(INTERP_KERNEL::Exception);
    int getNumberOfTS() const throw(INTERP_KERNEL::Exception);
    void eraseEmptyTS() throw(INTERP_KERNEL::Exception);
    int getPosOfTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    int getPosGivenTime(double time, double eps=1e-8) const throw(INTERP_KERNEL::Exception);
    void write(const std::string& fileName, int mode) const throw(INTERP_KERNEL::Exception);
    void loadArrays() throw(INTERP_KERNEL::Exception);
    void loadArraysIfNecessary() throw(INTERP_KERNEL::Exception);
    void unloadArrays() throw(INTERP_KERNEL::Exception);
    void unloadArraysWithoutDataLoss() throw(INTERP_KERNEL::Exception);
    //
    virtual MEDFileAnyTypeField1TS *getTimeStepAtPos(int pos) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TS *getTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TS *getTimeStepGivenTime(double time, double eps=1e-8) const throw(INTERP_KERNEL::Exception);
    void pushBackTimeStep(MEDFileAnyTypeField1TS *f1ts) throw(INTERP_KERNEL::Exception);
    void synchronizeNameScope() throw(INTERP_KERNEL::Exception);
    %extend
    {
      int __len__() const throw(INTERP_KERNEL::Exception)
      {
        return self->getNumberOfTS();
      }

      int getTimeId(PyObject *elt0) const throw(INTERP_KERNEL::Exception)
      {
        if(elt0 && PyInt_Check(elt0))
          {//fmts[3]
            int pos=PyInt_AS_LONG(elt0);
            return pos;
          }
        else if(elt0 && PyTuple_Check(elt0))
          {
            if(PyTuple_Size(elt0)==2)
              {
                PyObject *o0=PyTuple_GetItem(elt0,0);
                PyObject *o1=PyTuple_GetItem(elt0,1);
                if(PyInt_Check(o0) && PyInt_Check(o1))
                  {//fmts(1,-1)
                    int iter=PyInt_AS_LONG(o0);
                    int order=PyInt_AS_LONG(o1);
                    return self->getPosOfTimeStep(iter,order);
                  }
                else
                  throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::__getitem__ : invalid input param ! input is a tuple of size 2 but two integers are expected in this tuple to request a time steps !");
              }
            else
              throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::__getitem__ : invalid input param ! input is a tuple of size != 2 ! two integers are expected in this tuple to request a time steps !");
          }
        else if(elt0 && PyFloat_Check(elt0))
          {
            double val=PyFloat_AS_DOUBLE(elt0);
            return self->getPosGivenTime(val);
          }
        else
          throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::__getitem__ : invalid input params ! expected fmts[int], fmts[int,int] or fmts[double] to request time step !");
      }
      
      PyObject *getIterations() const throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > res=self->getIterations();
        PyObject *ret=PyList_New(res.size());
        int rk=0;
        for(std::vector< std::pair<int,int> >::const_iterator iter=res.begin();iter!=res.end();iter++,rk++)
          {
            PyObject *elt=PyTuple_New(2);
            PyTuple_SetItem(elt,0,SWIG_From_int((*iter).first));
            PyTuple_SetItem(elt,1,SWIG_From_int((*iter).second));
            PyList_SetItem(ret,rk,elt);
          }
        return ret;
      }
      
      PyObject *getTimeSteps() const throw(INTERP_KERNEL::Exception)
      {
        std::vector<double> ret1;
        std::vector< std::pair<int,int> > ret=self->getTimeSteps(ret1);
        std::size_t sz=ret.size();
        PyObject *ret2=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          {
            PyObject *elt=PyTuple_New(3);
            PyTuple_SetItem(elt,0,SWIG_From_int(ret[i].first));
            PyTuple_SetItem(elt,1,SWIG_From_int(ret[i].second));
            PyTuple_SetItem(elt,2,SWIG_From_double(ret1[i]));
            PyList_SetItem(ret2,i,elt);
          }
        return ret2;
      }
      
      PyObject *getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::vector<TypeOfField> > ret=self->getTypesOfFieldAvailable();
        PyObject *ret2=PyList_New(ret.size());
        for(int i=0;i<(int)ret.size();i++)
          {
            const std::vector<TypeOfField>& rett=ret[i];
            PyObject *ret3=PyList_New(rett.size());
            for(int j=0;j<(int)rett.size();j++)
              PyList_SetItem(ret3,j,SWIG_From_int(rett[j]));
            PyList_SetItem(ret2,i,ret3);
          }
        return ret2;
      }
      
      PyObject *getNonEmptyLevels(int iteration, int order, const std::string& mname=std::string()) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> ret1;
        int ret0=self->getNonEmptyLevels(iteration,order,mname,ret1);
        PyObject *elt=PyTuple_New(2);
        PyTuple_SetItem(elt,0,SWIG_From_int(ret0));
        PyTuple_SetItem(elt,1,convertIntArrToPyList2(ret1));
        return elt;
      }
      
      PyObject *getFieldSplitedByType(int iteration, int order, const std::string& mname=std::string()) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<INTERP_KERNEL::NormalizedCellType> types;
        std::vector< std::vector<TypeOfField> > typesF;
        std::vector< std::vector<std::string> > pfls;
        std::vector< std::vector<std::string> > locs;
        std::vector< std::vector< std::pair<int,int> > > ret=self->getFieldSplitedByType(iteration,order,mname,types,typesF,pfls,locs);
        int sz=ret.size();
        PyObject *ret2=PyList_New(sz);
        for(int i=0;i<sz;i++)
          {
            const std::vector< std::pair<int,int> >& dadsI=ret[i];
            const std::vector<TypeOfField>& typesFI=typesF[i];
            const std::vector<std::string>& pflsI=pfls[i];
            const std::vector<std::string>& locsI=locs[i];
            PyObject *elt=PyTuple_New(2);
            PyTuple_SetItem(elt,0,SWIG_From_int(types[i]));
            int sz2=ret[i].size();
            PyObject *elt2=PyList_New(sz2);
            for(int j=0;j<sz2;j++)
              {
                PyObject *elt3=PyTuple_New(4);
                PyTuple_SetItem(elt3,0,SWIG_From_int(typesFI[j]));
                PyObject *elt4=PyTuple_New(2); PyTuple_SetItem(elt4,0,SWIG_From_int(dadsI[j].first)); PyTuple_SetItem(elt4,1,SWIG_From_int(dadsI[j].second));
                PyTuple_SetItem(elt3,1,elt4);
                PyTuple_SetItem(elt3,2,PyString_FromString(pflsI[j].c_str()));
                PyTuple_SetItem(elt3,3,PyString_FromString(locsI[j].c_str()));
                PyList_SetItem(elt2,j,elt3);
              }
            PyTuple_SetItem(elt,1,elt2);
            PyList_SetItem(ret2,i,elt);
          }
        return ret2;
      }

      std::vector<int> getTimeIds(PyObject *elts) const throw(INTERP_KERNEL::Exception)
      {
        if(PyList_Check(elts))
          {
            int sz=PyList_Size(elts);
            std::vector<int> ret(sz);
            for(int i=0;i<sz;i++)
              {
                PyObject *elt=PyList_GetItem(elts,i);
                ret[i]=ParaMEDMEM_MEDFileAnyTypeFieldMultiTS_getTimeId(self,elt);
              }
            return ret;
          }
        else
          {
            std::vector<int> ret(1);
            ret[0]=ParaMEDMEM_MEDFileAnyTypeFieldMultiTS_getTimeId(self,elts);
            return ret;
          }
      }
      
      void __delitem__(PyObject *elts) throw(INTERP_KERNEL::Exception)
      {
        if(PySlice_Check(elts))
          {
            Py_ssize_t strt=2,stp=2,step=2;
            PySliceObject *oC=reinterpret_cast<PySliceObject *>(elts);
            GetIndicesOfSlice(oC,self->getNumberOfTS(),&strt,&stp,&step,"MEDFileAnyTypeFieldMultiTS.__delitem__ : error in input slice !");
            self->eraseTimeStepIds2(strt,stp,step);
          }
        else
          {
            std::vector<int> idsToRemove=ParaMEDMEM_MEDFileAnyTypeFieldMultiTS_getTimeIds(self,elts);
            if(!idsToRemove.empty())
              self->eraseTimeStepIds(&idsToRemove[0],&idsToRemove[0]+idsToRemove.size());
          }
      }
      
      void eraseTimeStepIds(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int pos1;
        std::vector<int> pos2;
        DataArrayInt *pos3=0;
        DataArrayIntTuple *pos4=0;
        convertObjToPossibleCpp1(li,sw,pos1,pos2,pos3,pos4);
        switch(sw)
          {
          case 1:
            {
              self->eraseTimeStepIds(&pos1,&pos1+1);
              return;
            }
          case 2:
            {
              if(pos2.empty())
                return;
              self->eraseTimeStepIds(&pos2[0],&pos2[0]+pos2.size());
              return ;
            }
          case 3:
            {
              self->eraseTimeStepIds(pos3->begin(),pos3->end());
              return ;
            }
          default:
            throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::eraseTimeStepIds : unexpected input array type recognized !");
          }
      }

      MEDFileAnyTypeFieldMultiTSIterator *__iter__() throw(INTERP_KERNEL::Exception)
      {
        return self->iterator();
      }

      PyObject *__getitem__(PyObject *elt0) const throw(INTERP_KERNEL::Exception)
      {
        if(elt0 && PyList_Check(elt0))
          {
            int sz=PyList_Size(elt0);
            MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da=DataArrayInt::New(); da->alloc(sz,1);
            int *pt=da->getPointer();
            for(int i=0;i<sz;i++,pt++)
              {
                PyObject *elt1=PyList_GetItem(elt0,i);
                *pt=MEDFileAnyTypeFieldMultiTSgetitemSingleTS__(self,elt1);
              }
            return convertMEDFileFieldMultiTS(self->buildSubPart(da->begin(),da->end()),SWIG_POINTER_OWN | 0);
          }
        else if(elt0 && PySlice_Check(elt0))
          {
            Py_ssize_t strt=2,stp=2,step=2;
            PySliceObject *oC=reinterpret_cast<PySliceObject *>(elt0);
            GetIndicesOfSlice(oC,self->getNumberOfTS(),&strt,&stp,&step,"MEDFileAnyTypeFieldMultiTS.__getitem__ : error in input slice !");
            return convertMEDFileFieldMultiTS(self->buildSubPartSlice(strt,stp,step),SWIG_POINTER_OWN | 0);
          }
        else
          return convertMEDFileField1TS(self->getTimeStepAtPos(MEDFileAnyTypeFieldMultiTSgetitemSingleTS__(self,elt0)),SWIG_POINTER_OWN | 0);
      }

      bool changeMeshNames(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<std::string,std::string> > modifTab=convertVecPairStStFromPy(li);
        return self->changeMeshNames(modifTab);
      }

      PyObject *splitComponents() const throw(INTERP_KERNEL::Exception)
      {
        std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTS > > ret=self->splitComponents();
        std::size_t sz=ret.size();
        PyObject *retPy=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(retPy,i,convertMEDFileFieldMultiTS(ret[i].retn(), SWIG_POINTER_OWN | 0 ));
        return retPy;
      }

      PyObject *splitDiscretizations() const throw(INTERP_KERNEL::Exception)
      {
        std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTS > > ret=self->splitDiscretizations();
        std::size_t sz=ret.size();
        PyObject *retPy=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(retPy,i,convertMEDFileFieldMultiTS(ret[i].retn(), SWIG_POINTER_OWN | 0 ));
        return retPy;
      }

      PyObject *splitMultiDiscrPerGeoTypes() const throw(INTERP_KERNEL::Exception)
      {
        std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTS > > ret=self->splitMultiDiscrPerGeoTypes();
        std::size_t sz=ret.size();
        PyObject *retPy=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(retPy,i,convertMEDFileFieldMultiTS(ret[i].retn(), SWIG_POINTER_OWN | 0 ));
        return retPy;
      }

      void pushBackTimeSteps(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        void *argp(0);
        int status(SWIG_ConvertPtr(li,&argp,SWIGTYPE_p_ParaMEDMEM__MEDFileAnyTypeFieldMultiTS,0|0));
        if(SWIG_IsOK(status))
          {
            self->pushBackTimeSteps(reinterpret_cast<MEDFileAnyTypeFieldMultiTS *>(argp));
          }
        else
          {
            std::vector<MEDFileAnyTypeField1TS *> tmp;
            convertFromPyObjVectorOfObj<ParaMEDMEM::MEDFileAnyTypeField1TS *>(li,SWIGTYPE_p_ParaMEDMEM__MEDFileAnyTypeField1TS,"MEDFileAnyTypeField1TS",tmp);
            self->pushBackTimeSteps(tmp);
          }
      }

      static PyObject *MEDFileAnyTypeFieldMultiTS::SplitIntoCommonTimeSeries(PyObject *li) throw(INTERP_KERNEL::Exception)
      {
        std::vector<MEDFileAnyTypeFieldMultiTS *> vectFMTS;
        convertFromPyObjVectorOfObj<ParaMEDMEM::MEDFileAnyTypeFieldMultiTS *>(li,SWIGTYPE_p_ParaMEDMEM__MEDFileAnyTypeFieldMultiTS,"MEDFileAnyTypeFieldMultiTS",vectFMTS);
        std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > ret=MEDFileAnyTypeFieldMultiTS::SplitIntoCommonTimeSeries(vectFMTS);
        std::size_t sz=ret.size();
        PyObject *retPy=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          {
            std::size_t sz2=ret[i].size();
            PyObject *ret1Py=PyList_New(sz2);
            for(std::size_t j=0;j<sz2;j++)
              {
                MEDFileAnyTypeFieldMultiTS *elt(ret[i][j]);
                if(elt)
                  elt->incrRef();
                PyList_SetItem(ret1Py,j,convertMEDFileFieldMultiTS(elt,SWIG_POINTER_OWN | 0 ));
              }
            PyList_SetItem(retPy,i,ret1Py);
          }
        return retPy;
      }
      
      static PyObject *MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupport(PyObject *li, const MEDFileMesh *mesh) throw(INTERP_KERNEL::Exception)
      {
        std::vector<MEDFileAnyTypeFieldMultiTS *> vectFMTS;
        convertFromPyObjVectorOfObj<ParaMEDMEM::MEDFileAnyTypeFieldMultiTS *>(li,SWIGTYPE_p_ParaMEDMEM__MEDFileAnyTypeFieldMultiTS,"MEDFileAnyTypeFieldMultiTS",vectFMTS);
        std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFastCellSupportComparator> > ret2;
        std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > ret=MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupport(vectFMTS,mesh,ret2);
        if(ret2.size()!=ret.size())
          {
            std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupport (PyWrap) : internal error ! Size of 2 vectors must match ! (" << ret.size() << "!=" << ret2.size() << ") !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
        std::size_t sz=ret.size();
        PyObject *retPy=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          {
            std::size_t sz2=ret[i].size();
            PyObject *ret0Py=PyTuple_New(2);
            PyObject *ret1Py=PyList_New(sz2);
            for(std::size_t j=0;j<sz2;j++)
              {
                MEDFileAnyTypeFieldMultiTS *elt(ret[i][j]);
                if(elt)
                  elt->incrRef();
                PyList_SetItem(ret1Py,j,convertMEDFileFieldMultiTS(elt,SWIG_POINTER_OWN | 0 ));
              }
            PyTuple_SetItem(ret0Py,0,ret1Py);
            PyTuple_SetItem(ret0Py,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret2[i].retn()),SWIGTYPE_p_ParaMEDMEM__MEDFileFastCellSupportComparator, SWIG_POINTER_OWN | 0 ));
            PyList_SetItem(retPy,i,ret0Py);
          }
        return retPy;
      }
    }
  };

  class MEDFileFieldMultiTS : public MEDFileAnyTypeFieldMultiTS
  {
  public:
    static MEDFileFieldMultiTS *New() throw(INTERP_KERNEL::Exception);
    static MEDFileFieldMultiTS *New(const std::string& fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileFieldMultiTS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    //
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, int iteration, int order, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtLevelOld(TypeOfField type, const std::string& mname, int iteration, int order, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    //
    void appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception);
    void appendFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception);
    ParaMEDMEM::MEDFileIntFieldMultiTS *convertToInt(bool isDeepCpyGlobs=true) const throw(INTERP_KERNEL::Exception);
    %extend
       {
         MEDFileFieldMultiTS()
         {
           return MEDFileFieldMultiTS::New();
         }

         MEDFileFieldMultiTS(const std::string& fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileFieldMultiTS::New(fileName,loadAll);
         }

         MEDFileFieldMultiTS(const std::string& fileName, const std::string& fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileFieldMultiTS::New(fileName,fieldName,loadAll);
         }

         static MEDFileFieldMultiTS *LoadSpecificEntities(const std::string& fileName, const std::string& fieldName, PyObject *entities, bool loadAll=true)
         {
           std::vector<std::pair<int,int> > tmp(convertTimePairIdsFromPy(entities));
           std::size_t sz(tmp.size());
           std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > entitiesCpp(sz);
           for(std::size_t i=0;i<sz;i++)
             {
               entitiesCpp[i].first=(TypeOfField)tmp[i].first;
               entitiesCpp[i].second=(INTERP_KERNEL::NormalizedCellType)tmp[i].second;
             }
           return MEDFileFieldMultiTS::LoadSpecificEntities(fileName,fieldName,entitiesCpp,loadAll);
         }
         
         std::string __str__() const throw(INTERP_KERNEL::Exception)
         {
           return self->simpleRepr();
         }

         PyObject *getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh) const throw(INTERP_KERNEL::Exception)
         {
           DataArrayInt *ret1=0;
           DataArrayDouble *ret0=self->getFieldWithProfile(type,iteration,order,meshDimRelToMax,mesh,ret1);
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           return ret;
         }

         PyObject *getFieldSplitedByType2(int iteration, int order, const std::string& mname=std::string()) const throw(INTERP_KERNEL::Exception)
         {
           std::vector<INTERP_KERNEL::NormalizedCellType> types;
           std::vector< std::vector<TypeOfField> > typesF;
           std::vector< std::vector<std::string> > pfls;
           std::vector< std::vector<std::string> > locs;
           std::vector< std::vector<DataArrayDouble *> > ret=self->getFieldSplitedByType2(iteration,order,mname,types,typesF,pfls,locs);
           int sz=ret.size();
           PyObject *ret2=PyList_New(sz);
           for(int i=0;i<sz;i++)
             {
               const std::vector<DataArrayDouble *>& dadsI=ret[i];
               const std::vector<TypeOfField>& typesFI=typesF[i];
               const std::vector<std::string>& pflsI=pfls[i];
               const std::vector<std::string>& locsI=locs[i];
               PyObject *elt=PyTuple_New(2);
               PyTuple_SetItem(elt,0,SWIG_From_int(types[i]));
               int sz2=ret[i].size();
               PyObject *elt2=PyList_New(sz2);
               for(int j=0;j<sz2;j++)
                 {
                   PyObject *elt3=PyTuple_New(4);
                   PyTuple_SetItem(elt3,0,SWIG_From_int(typesFI[j]));
                   PyTuple_SetItem(elt3,1,SWIG_NewPointerObj(SWIG_as_voidptr(dadsI[j]),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
                   PyTuple_SetItem(elt3,2,PyString_FromString(pflsI[j].c_str()));
                   PyTuple_SetItem(elt3,3,PyString_FromString(locsI[j].c_str()));
                   PyList_SetItem(elt2,j,elt3);
                 }
               PyTuple_SetItem(elt,1,elt2);
               PyList_SetItem(ret2,i,elt);
             }
           return ret2;
         }
         DataArrayDouble *getUndergroundDataArray(int iteration, int order) const throw(INTERP_KERNEL::Exception)
         {
           DataArrayDouble *ret=self->getUndergroundDataArray(iteration,order);
           if(ret)
             ret->incrRef();
           return ret;
         }
         
         PyObject *getUndergroundDataArrayExt(int iteration, int order) const throw(INTERP_KERNEL::Exception)
         {
           std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > > elt1Cpp;
           DataArrayDouble *elt0=self->getUndergroundDataArrayExt(iteration,order,elt1Cpp);
           if(elt0)
             elt0->incrRef();
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(elt0),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
           std::size_t sz=elt1Cpp.size();
           PyObject *elt=PyList_New(sz);
           for(std::size_t i=0;i<sz;i++)
             {
               PyObject *elt1=PyTuple_New(2);
               PyObject *elt2=PyTuple_New(2);
               PyTuple_SetItem(elt2,0,SWIG_From_int(elt1Cpp[i].first.first));
               PyTuple_SetItem(elt2,1,SWIG_From_int(elt1Cpp[i].first.second));
               PyObject *elt3=PyTuple_New(2);
               PyTuple_SetItem(elt3,0,SWIG_From_int(elt1Cpp[i].second.first));
               PyTuple_SetItem(elt3,1,SWIG_From_int(elt1Cpp[i].second.second));
               PyTuple_SetItem(elt1,0,elt2);
               PyTuple_SetItem(elt1,1,elt3);
               PyList_SetItem(elt,i,elt1);
             }
           PyTuple_SetItem(ret,1,elt);
           return ret;
         }
       }
  };

  class MEDFileFieldsIterator
  {
  public:
    %extend
    {
      PyObject *next() throw(INTERP_KERNEL::Exception)
      {
        MEDFileAnyTypeFieldMultiTS *ret=self->nextt();
        if(ret)
          return convertMEDFileFieldMultiTS(ret, SWIG_POINTER_OWN | 0 );
        else
          {
            PyErr_SetString(PyExc_StopIteration,"No more data.");
            return 0;
          }
      }
    }
  };

  class MEDFileIntFieldMultiTS : public MEDFileAnyTypeFieldMultiTS
  {
  public:
    static MEDFileIntFieldMultiTS *New();
    static MEDFileIntFieldMultiTS *New(const std::string& fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileIntFieldMultiTS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    //
    void appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals) throw(INTERP_KERNEL::Exception);
    void appendFieldProfile(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception);
    ParaMEDMEM::MEDFileFieldMultiTS *convertToDouble(bool isDeepCpyGlobs=true) const throw(INTERP_KERNEL::Exception);
    %extend
    {
      MEDFileIntFieldMultiTS()
      {
        return MEDFileIntFieldMultiTS::New();
      }
      
      MEDFileIntFieldMultiTS(const std::string& fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception)
      {
        return MEDFileIntFieldMultiTS::New(fileName,loadAll);
      }
      
      MEDFileIntFieldMultiTS(const std::string& fileName, const std::string& fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception)
      {
        return MEDFileIntFieldMultiTS::New(fileName,fieldName,loadAll);
      }

      static MEDFileIntFieldMultiTS *LoadSpecificEntities(const std::string& fileName, const std::string& fieldName, PyObject *entities, bool loadAll=true)
      {
        std::vector<std::pair<int,int> > tmp(convertTimePairIdsFromPy(entities));
        std::size_t sz(tmp.size());
        std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > entitiesCpp(sz);
        for(std::size_t i=0;i<sz;i++)
          {
            entitiesCpp[i].first=(TypeOfField)tmp[i].first;
            entitiesCpp[i].second=(INTERP_KERNEL::NormalizedCellType)tmp[i].second;
          }
        return MEDFileIntFieldMultiTS::LoadSpecificEntities(fileName,fieldName,entitiesCpp,loadAll);
      }

      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->simpleRepr();
      }

      PyObject *getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0;
        MEDCouplingFieldDouble *ret0=self->getFieldAtLevel(type,iteration,order,meshDimRelToMax,ret1,renumPol);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *getFieldAtTopLevel(TypeOfField type, int iteration, int order, int renumPol=0) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0;
        MEDCouplingFieldDouble *ret0=self->getFieldAtTopLevel(type,iteration,order,ret1,renumPol);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0;
        MEDCouplingFieldDouble *ret0=self->getFieldOnMeshAtLevel(type,iteration,order,meshDimRelToMax,mesh,ret1,renumPol);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }
      
      PyObject *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0;
        MEDCouplingFieldDouble *ret0=self->getFieldOnMeshAtLevel(type,iteration,order,mesh,ret1,renumPol);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }
      
      PyObject *getFieldAtLevelOld(TypeOfField type, int iteration, int order, const std::string& mname, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0;
        MEDCouplingFieldDouble *ret0=self->getFieldAtLevelOld(type,iteration,order,mname,meshDimRelToMax,ret1,renumPol);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      PyObject *getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh) const throw(INTERP_KERNEL::Exception)
      {
         DataArrayInt *ret1=0;
         DataArrayInt *ret0=self->getFieldWithProfile(type,iteration,order,meshDimRelToMax,mesh,ret1);
         PyObject *ret=PyTuple_New(2);
         PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
         PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
         return ret;
      }

      DataArrayInt *getUndergroundDataArray(int iteration, int order) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret=self->getUndergroundDataArray(iteration,order);
        if(ret)
          ret->incrRef();
        return ret;
      }
    }
  };

  class MEDFileFields : public RefCountObject, public MEDFileFieldGlobsReal, public MEDFileWritable
  {
  public:
    static MEDFileFields *New() throw(INTERP_KERNEL::Exception);
    static MEDFileFields *New(const std::string& fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileFields *LoadPartOf(const std::string& fileName, bool loadAll=true, const MEDFileMeshes *ms=0) throw(INTERP_KERNEL::Exception);
    MEDFileFields *deepCpy() const throw(INTERP_KERNEL::Exception);
    MEDFileFields *shallowCpy() const throw(INTERP_KERNEL::Exception);
    void loadArrays() throw(INTERP_KERNEL::Exception);
    void loadArraysIfNecessary() throw(INTERP_KERNEL::Exception);
    void unloadArrays() throw(INTERP_KERNEL::Exception);
    void unloadArraysWithoutDataLoss() throw(INTERP_KERNEL::Exception);
    void write(const std::string& fileName, int mode) const throw(INTERP_KERNEL::Exception);
    int getNumberOfFields() const;
    std::vector<std::string> getFieldsNames() const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getMeshesNames() const throw(INTERP_KERNEL::Exception);
    //
    void resize(int newSize) throw(INTERP_KERNEL::Exception);
    void pushField(MEDFileAnyTypeFieldMultiTS *field) throw(INTERP_KERNEL::Exception);
    void setFieldAtPos(int i, MEDFileAnyTypeFieldMultiTS *field) throw(INTERP_KERNEL::Exception);
    int getPosFromFieldName(const std::string& fieldName) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTS *getFieldAtPos(int i) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTS *getFieldWithName(const std::string& fieldName) const throw(INTERP_KERNEL::Exception);
    MEDFileFields *partOfThisLyingOnSpecifiedMeshName(const std::string& meshName) const throw(INTERP_KERNEL::Exception);
    void destroyFieldAtPos(int i) throw(INTERP_KERNEL::Exception);
    bool removeFieldsWithoutAnyTimeStep() throw(INTERP_KERNEL::Exception);
    %extend
       {
         MEDFileFields()
         {
           return MEDFileFields::New();
         }

         MEDFileFields(const std::string& fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileFields::New(fileName,loadAll);
         }
         
         std::string __str__() const throw(INTERP_KERNEL::Exception)
         {
           return self->simpleRepr();
         }

         static MEDFileFields *LoadSpecificEntities(const std::string& fileName, PyObject *entities, bool loadAll=true) throw(INTERP_KERNEL::Exception)
         {
           std::vector<std::pair<int,int> > tmp(convertTimePairIdsFromPy(entities));
           std::size_t sz(tmp.size());
           std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > entitiesCpp(sz);
           for(std::size_t i=0;i<sz;i++)
             {
               entitiesCpp[i].first=(TypeOfField)tmp[i].first;
               entitiesCpp[i].second=(INTERP_KERNEL::NormalizedCellType)tmp[i].second;
             }
           return MEDFileFields::LoadSpecificEntities(fileName,entitiesCpp,loadAll);
         }

         PyObject *getCommonIterations() const throw(INTERP_KERNEL::Exception)
         {
           bool ret1;
           std::vector< std::pair<int,int> > ret0=self->getCommonIterations(ret1);
           PyObject *ret=PyTuple_New(2);
           PyObject *ret_0=PyList_New(ret0.size());
           int rk=0;
           for(std::vector< std::pair<int,int> >::const_iterator iter=ret0.begin();iter!=ret0.end();iter++,rk++)
             {
               PyObject *elt=PyTuple_New(2);
               PyTuple_SetItem(elt,0,SWIG_From_int((*iter).first));
               PyTuple_SetItem(elt,1,SWIG_From_int((*iter).second));
               PyList_SetItem(ret_0,rk,elt);
             }
           PyTuple_SetItem(ret,0,ret_0);
           PyObject *ret_1=ret1?Py_True:Py_False; Py_XINCREF(ret_1);
           PyTuple_SetItem(ret,1,ret_1);
           return ret;
         }

         MEDFileFields *partOfThisLyingOnSpecifiedTimeSteps(PyObject *timeSteps) const throw(INTERP_KERNEL::Exception)
         {
           std::vector< std::pair<int,int> > ts=convertTimePairIdsFromPy(timeSteps);
           return self->partOfThisLyingOnSpecifiedTimeSteps(ts);
         }

         MEDFileFields *partOfThisNotLyingOnSpecifiedTimeSteps(PyObject *timeSteps) const throw(INTERP_KERNEL::Exception)
         {
           std::vector< std::pair<int,int> > ts=convertTimePairIdsFromPy(timeSteps);
           return self->partOfThisNotLyingOnSpecifiedTimeSteps(ts);
         }
         
         PyObject *__getitem__(PyObject *obj) throw(INTERP_KERNEL::Exception)
         {
           if(obj && PyList_Check(obj))
             {
               int sz=PyList_Size(obj);
               MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da=DataArrayInt::New(); da->alloc(sz,1);
               int *pt=da->getPointer();
               for(int i=0;i<sz;i++,pt++)
                 {
                   PyObject *elt1=PyList_GetItem(obj,i);
                   *pt=MEDFileFieldsgetitemSingleTS__(self,elt1);
                 }
               return SWIG_NewPointerObj(SWIG_as_voidptr(self->buildSubPart(da->begin(),da->end())),SWIGTYPE_p_ParaMEDMEM__MEDFileFields, SWIG_POINTER_OWN | 0 );
             }
           else
             return convertMEDFileFieldMultiTS(self->getFieldAtPos(MEDFileFieldsgetitemSingleTS__(self,obj)), SWIG_POINTER_OWN | 0 );
         }

         MEDFileFields *__setitem__(int obj, MEDFileFieldMultiTS *field) throw(INTERP_KERNEL::Exception)
         {
           self->setFieldAtPos(obj,field);
           return self;
         }

         int __len__() const throw(INTERP_KERNEL::Exception)
         {
           return self->getNumberOfFields();
         }

         MEDFileFieldsIterator *__iter__() throw(INTERP_KERNEL::Exception)
         {
           return self->iterator();
         }
         
         bool changeMeshNames(PyObject *li) throw(INTERP_KERNEL::Exception)
         {
           std::vector< std::pair<std::string,std::string> > modifTab=convertVecPairStStFromPy(li);
           return self->changeMeshNames(modifTab);
         }

         int getPosOfField(PyObject *elt0) const throw(INTERP_KERNEL::Exception)
         {
           if(elt0 && PyInt_Check(elt0))
             {//fmts[3]
               return PyInt_AS_LONG(elt0);
             }
           else if(elt0 && PyString_Check(elt0))
             return self->getPosFromFieldName(PyString_AsString(elt0));
           else
             throw INTERP_KERNEL::Exception("MEDFileFields::getPosOfField : invalid input params ! expected fields[int], fields[string_of_field_name] !");
         }
         
         std::vector<int> getPosOfFields(PyObject *elts) const throw(INTERP_KERNEL::Exception)
         {
           if(PyList_Check(elts))
             {
               int sz=PyList_Size(elts);
               std::vector<int> ret(sz);
               for(int i=0;i<sz;i++)
                 {
                   PyObject *elt=PyList_GetItem(elts,i);
                   ret[i]=ParaMEDMEM_MEDFileFields_getPosOfField(self,elt);
                 }
               return ret;
             }
           else
             {
               std::vector<int> ret(1);
               ret[0]=ParaMEDMEM_MEDFileFields_getPosOfField(self,elts);
               return ret;
             }
         }

         void pushFields(PyObject *fields) throw(INTERP_KERNEL::Exception)
         {
           std::vector<MEDFileAnyTypeFieldMultiTS *> tmp;
           convertFromPyObjVectorOfObj<ParaMEDMEM::MEDFileAnyTypeFieldMultiTS *>(fields,SWIGTYPE_p_ParaMEDMEM__MEDFileAnyTypeFieldMultiTS,"MEDFileAnyTypeFieldMultiTS",tmp);
           self->pushFields(tmp);
         }
         
         void __delitem__(PyObject *elts) throw(INTERP_KERNEL::Exception)
         {
           if(elts && PySlice_Check(elts))
             {
               Py_ssize_t strt=2,stp=2,step=2;
               PySliceObject *oC=reinterpret_cast<PySliceObject *>(elts);
               GetIndicesOfSlice(oC,self->getNumberOfFields(),&strt,&stp,&step,"MEDFileFields.__delitem__ : error in input slice !");
               self->destroyFieldsAtPos2(strt,stp,step);
             }
           else
             {
               std::vector<int> idsToRemove=ParaMEDMEM_MEDFileFields_getPosOfFields(self,elts);
               if(!idsToRemove.empty())
                 self->destroyFieldsAtPos(&idsToRemove[0],&idsToRemove[0]+idsToRemove.size());
             }
         }
       }
  };

  class MEDFileParameter1TS : public RefCountObject
  {
  public:
    void setIteration(int it);
    int getIteration() const;
    void setOrder(int order);
    int getOrder() const;
    void setTimeValue(double time);
    void setTime(int dt, int it, double time);
    double getTime(int& dt, int& it);
    double getTimeValue() const;
  };

  class MEDFileParameterDouble1TSWTI : public MEDFileParameter1TS
  {
  public:
    void setValue(double val) throw(INTERP_KERNEL::Exception);
    double getValue() const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const throw(INTERP_KERNEL::Exception);
    %extend
    {
      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->simpleRepr();
      }
    }
  };

  class MEDFileParameterTinyInfo : public MEDFileWritable
  {
  public:
    void setDescription(const std::string& name);
    std::string getDescription() const;
    void setTimeUnit(const std::string& unit);
    std::string getTimeUnit() const;
  };

  class MEDFileParameterDouble1TS : public MEDFileParameterDouble1TSWTI, public MEDFileParameterTinyInfo
  {
  public:
    static MEDFileParameterDouble1TS *New();
    static MEDFileParameterDouble1TS *New(const std::string& fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileParameterDouble1TS *New(const std::string& fileName, const std::string& paramName) throw(INTERP_KERNEL::Exception);
    static MEDFileParameterDouble1TS *New(const std::string& fileName, const std::string& paramName, int dt, int it) throw(INTERP_KERNEL::Exception);
    virtual MEDFileParameter1TS *deepCpy() const throw(INTERP_KERNEL::Exception);
    virtual std::string simpleRepr() const throw(INTERP_KERNEL::Exception);
    void setName(const std::string& name) throw(INTERP_KERNEL::Exception);
    std::string getName() const throw(INTERP_KERNEL::Exception);
    void write(const std::string& fileName, int mode) const throw(INTERP_KERNEL::Exception);
    %extend
    {
      MEDFileParameterDouble1TS()
      {
        return MEDFileParameterDouble1TS::New();
      }
      
      MEDFileParameterDouble1TS(const std::string& fileName) throw(INTERP_KERNEL::Exception)
      {
        return MEDFileParameterDouble1TS::New(fileName);
      }

      MEDFileParameterDouble1TS(const std::string& fileName, const std::string& paramName) throw(INTERP_KERNEL::Exception)
      {
        return MEDFileParameterDouble1TS::New(fileName,paramName);
      }

      MEDFileParameterDouble1TS(const std::string& fileName, const std::string& paramName, int dt, int it) throw(INTERP_KERNEL::Exception)
      {
        return MEDFileParameterDouble1TS::New(fileName,paramName,dt,it);
      }

      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->simpleRepr();
      }

      PyObject *isEqual(const MEDFileParameter1TS *other, double eps) const throw(INTERP_KERNEL::Exception)
      {
        std::string what;
        bool ret0=self->isEqual(other,eps,what);
        PyObject *res=PyList_New(2);
        PyObject *ret0Py=ret0?Py_True:Py_False;
        Py_XINCREF(ret0Py);
        PyList_SetItem(res,0,ret0Py);
        PyList_SetItem(res,1,PyString_FromString(what.c_str()));
        return res;
      }
    }
  };

  class MEDFileParameterMultiTS : public RefCountObject, public MEDFileParameterTinyInfo
  {
  public:
    static MEDFileParameterMultiTS *New();
    static MEDFileParameterMultiTS *New(const std::string& fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileParameterMultiTS *New(const std::string& fileName, const std::string& paramName) throw(INTERP_KERNEL::Exception);
    std::string getName() const;
    void setName(const std::string& name);
    MEDFileParameterMultiTS *deepCpy() const throw(INTERP_KERNEL::Exception);
    void write(const std::string& fileName, int mode) const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const throw(INTERP_KERNEL::Exception);
    void appendValue(int dt, int it, double time, double val) throw(INTERP_KERNEL::Exception);
    double getDoubleValue(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    int getPosOfTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    int getPosGivenTime(double time, double eps=1e-8) const throw(INTERP_KERNEL::Exception);
    int getNumberOfTS() const throw(INTERP_KERNEL::Exception);
    %extend
    {
      MEDFileParameterMultiTS()
      {
        return MEDFileParameterMultiTS::New();
      }
      
      MEDFileParameterMultiTS(const std::string& fileName)
      {
        return MEDFileParameterMultiTS::New(fileName);
      }

      MEDFileParameterMultiTS(const std::string& fileName, const std::string& paramName)
      {
        return MEDFileParameterMultiTS::New(fileName,paramName);
      }

      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->simpleRepr();
      }
      
      PyObject *isEqual(const MEDFileParameterMultiTS *other, double eps) const throw(INTERP_KERNEL::Exception)
      {
        std::string what;
        bool ret0=self->isEqual(other,eps,what);
        PyObject *res=PyList_New(2);
        PyObject *ret0Py=ret0?Py_True:Py_False;
        Py_XINCREF(ret0Py);
        PyList_SetItem(res,0,ret0Py);
        PyList_SetItem(res,1,PyString_FromString(what.c_str()));
        return res;
      }
      
      void eraseTimeStepIds(PyObject *ids) throw(INTERP_KERNEL::Exception)
      {
        int sw;
        int pos1;
        std::vector<int> pos2;
        DataArrayInt *pos3=0;
        DataArrayIntTuple *pos4=0;
        convertObjToPossibleCpp1(ids,sw,pos1,pos2,pos3,pos4);
        switch(sw)
          {
          case 1:
            {
              self->eraseTimeStepIds(&pos1,&pos1+1);
              return;
            }
          case 2:
            {
              if(pos2.empty())
                return;
              self->eraseTimeStepIds(&pos2[0],&pos2[0]+pos2.size());
              return ;
            }
          case 3:
            {
              self->eraseTimeStepIds(pos3->begin(),pos3->end());
              return ;
            }
          default:
            throw INTERP_KERNEL::Exception("MEDFileParameterMultiTS::eraseTimeStepIds : unexpected input array type recognized !");
          }
      }

      int getTimeStepId(PyObject *elt0) const throw(INTERP_KERNEL::Exception)
      {
        if(elt0 && PyInt_Check(elt0))
          {//fmts[3]
            int pos=InterpreteNegativeInt(PyInt_AS_LONG(elt0),self->getNumberOfTS());
            return pos;
          }
        else if(elt0 && PyTuple_Check(elt0))
          {
            if(PyTuple_Size(elt0)==2)
              {
                PyObject *o0=PyTuple_GetItem(elt0,0);
                PyObject *o1=PyTuple_GetItem(elt0,1);
                if(PyInt_Check(o0) && PyInt_Check(o1))
                  {//fmts(1,-1)
                    int iter=PyInt_AS_LONG(o0);
                    int order=PyInt_AS_LONG(o1);
                    return self->getPosOfTimeStep(iter,order);
                  }
                else
                  throw INTERP_KERNEL::Exception("MEDFileParameterMultiTS::getTimeStepId : invalid input param ! input is a tuple of size 2 but two integers are expected in this tuple to request a time steps !");
              }
            else
              throw INTERP_KERNEL::Exception("MEDFileParameterMultiTS::getTimeStepId : invalid input param ! input is a tuple of size != 2 ! two integers are expected in this tuple to request a time steps !");
          }
        else if(elt0 && PyFloat_Check(elt0))
          {
            double val=PyFloat_AS_DOUBLE(elt0);
            return self->getPosGivenTime(val);
          }
        else
          throw INTERP_KERNEL::Exception("MEDFileParameterMultiTS::getTimeStepId : invalid input params ! expected fmts[int], fmts[int,int] or fmts[double] to request time step !");
      }

      MEDFileParameter1TS *__getitem__(PyObject *elt0) const throw(INTERP_KERNEL::Exception)
      {
        MEDFileParameter1TS *ret=self->getTimeStepAtPos(ParaMEDMEM_MEDFileParameterMultiTS_getTimeStepId(self,elt0));
        if(ret)
          ret->incrRef();
        return ret;
      }

      std::vector<int> getTimeStepIds(PyObject *elts) const throw(INTERP_KERNEL::Exception)
      {
        if(PyList_Check(elts))
          {
            int sz=PyList_Size(elts);
            std::vector<int> ret(sz);
            for(int i=0;i<sz;i++)
              {
                PyObject *elt=PyList_GetItem(elts,i);
                ret[i]=ParaMEDMEM_MEDFileParameterMultiTS_getTimeStepId(self,elt);
              }
            return ret;
          }
        else
          {
            std::vector<int> ret(1);
            ret[0]=ParaMEDMEM_MEDFileParameterMultiTS_getTimeStepId(self,elts);
            return ret;
          }
      }

      void __delitem__(PyObject *elts) throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> idsToRemove=ParaMEDMEM_MEDFileParameterMultiTS_getTimeStepIds(self,elts);
        if(!idsToRemove.empty())
          self->eraseTimeStepIds(&idsToRemove[0],&idsToRemove[0]+idsToRemove.size());
      }
      
      MEDFileParameter1TS *getTimeStepAtPos(int posId) const throw(INTERP_KERNEL::Exception)
      {
        MEDFileParameter1TS *ret=self->getTimeStepAtPos(posId);
        if(ret)
          ret->incrRef();
        return ret;
      }

      PyObject *getIterations() const throw(INTERP_KERNEL::Exception)
      {
        std::vector< std::pair<int,int> > res=self->getIterations();
        PyObject *ret=PyList_New(res.size());
        int rk=0;
        for(std::vector< std::pair<int,int> >::const_iterator iter=res.begin();iter!=res.end();iter++,rk++)
          {
            PyObject *elt=PyTuple_New(2);
            PyTuple_SetItem(elt,0,SWIG_From_int((*iter).first));
            PyTuple_SetItem(elt,1,SWIG_From_int((*iter).second));
            PyList_SetItem(ret,rk,elt);
          }
        return ret;
      }

      PyObject *getTimeSteps() const throw(INTERP_KERNEL::Exception)
      {
        std::vector<double> res2;
        std::vector< std::pair<int,int> > res=self->getTimeSteps(res2);
        PyObject *ret=PyList_New(res.size());
        int rk=0;
        for(std::vector< std::pair<int,int> >::const_iterator iter=res.begin();iter!=res.end();iter++,rk++)
          {
            PyObject *elt=PyTuple_New(3);
            PyTuple_SetItem(elt,0,SWIG_From_int((*iter).first));
            PyTuple_SetItem(elt,1,SWIG_From_int((*iter).second));
            PyTuple_SetItem(elt,2,SWIG_From_double(res2[rk]));
            PyList_SetItem(ret,rk,elt);
          }
        return ret;
      }
    }
  };

  class MEDFileParameters : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileParameters *New();
    static MEDFileParameters *New(const std::string& fileName) throw(INTERP_KERNEL::Exception);
    MEDFileParameters *deepCpy() const throw(INTERP_KERNEL::Exception);
    void write(const std::string& fileName, int mode) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getParamsNames() const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const throw(INTERP_KERNEL::Exception);
    void resize(int newSize) throw(INTERP_KERNEL::Exception);
    void pushParam(MEDFileParameterMultiTS *param) throw(INTERP_KERNEL::Exception);
    void setParamAtPos(int i, MEDFileParameterMultiTS *param) throw(INTERP_KERNEL::Exception);
    void destroyParamAtPos(int i) throw(INTERP_KERNEL::Exception);
    int getPosFromParamName(const std::string& paramName) const throw(INTERP_KERNEL::Exception);
    int getNumberOfParams() const throw(INTERP_KERNEL::Exception);
    %extend
    {
      MEDFileParameters()
      {
        return MEDFileParameters::New();
      }
      
      MEDFileParameters(const std::string& fileName)
      {
        return MEDFileParameters::New(fileName);
      }

      std::string __str__() const throw(INTERP_KERNEL::Exception)
      {
        return self->simpleRepr();
      }

      MEDFileParameterMultiTS *__getitem__(PyObject *obj) throw(INTERP_KERNEL::Exception)
      {
        if(PyInt_Check(obj))
          {
            MEDFileParameterMultiTS *ret=self->getParamAtPos(InterpreteNegativeInt((int)PyInt_AS_LONG(obj),self->getNumberOfParams()));
            if(ret)
              ret->incrRef();
            return ret;
          }
        else if(PyString_Check(obj))
          {
            MEDFileParameterMultiTS *ret=self->getParamWithName(PyString_AsString(obj));
            if(ret)
              ret->incrRef();
            return ret;
          }
        else
          throw INTERP_KERNEL::Exception("MEDFileParameters::__getitem__ : only integer or string with meshname supported !");
      }

      int __len__() const throw(INTERP_KERNEL::Exception)
      {
        return self->getNumberOfParams();
      }
      
      MEDFileParameterMultiTS *getParamAtPos(int i) const throw(INTERP_KERNEL::Exception)
      {
        MEDFileParameterMultiTS *ret=self->getParamAtPos(i);
        if(ret)
          ret->incrRef();
        return ret;
      }

      MEDFileParameterMultiTS *getParamWithName(const std::string& paramName) const throw(INTERP_KERNEL::Exception)
      {
        MEDFileParameterMultiTS *ret=self->getParamWithName(paramName);
        if(ret)
          ret->incrRef();
        return ret;
      }
      
      PyObject *isEqual(const MEDFileParameters *other, double eps) const throw(INTERP_KERNEL::Exception)
      {
        std::string what;
        bool ret0=self->isEqual(other,eps,what);
        PyObject *res=PyList_New(2);
        PyObject *ret0Py=ret0?Py_True:Py_False;
        Py_XINCREF(ret0Py);
        PyList_SetItem(res,0,ret0Py);
        PyList_SetItem(res,1,PyString_FromString(what.c_str()));
        return res;
      }
    }
  };

  class MEDFileData : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileData *New(const std::string& fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileData *New();
    MEDFileData *deepCpy() const throw(INTERP_KERNEL::Exception);
    void setFields(MEDFileFields *fields) throw(INTERP_KERNEL::Exception);
    void setMeshes(MEDFileMeshes *meshes) throw(INTERP_KERNEL::Exception);
    void setParams(MEDFileParameters *params) throw(INTERP_KERNEL::Exception);
    int getNumberOfFields() const throw(INTERP_KERNEL::Exception);
    int getNumberOfMeshes() const throw(INTERP_KERNEL::Exception);
    int getNumberOfParams() const throw(INTERP_KERNEL::Exception);
    //
    bool changeMeshName(const std::string& oldMeshName, const std::string& newMeshName) throw(INTERP_KERNEL::Exception);
    bool unPolyzeMeshes() throw(INTERP_KERNEL::Exception);
    //
    void write(const std::string& fileName, int mode) const throw(INTERP_KERNEL::Exception);
    %extend
       {
         MEDFileData(const std::string& fileName) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileData::New(fileName);
         }

         MEDFileData()
         {
           return MEDFileData::New();
         }

         std::string __str__() const throw(INTERP_KERNEL::Exception)
         {
           return self->simpleRepr();
         }

         MEDFileMeshes *getMeshes() const throw(INTERP_KERNEL::Exception)
         {
           MEDFileMeshes *ret=self->getMeshes();
           if(ret)
             ret->incrRef();
           return ret;
         }

         MEDFileParameters *getParams() const throw(INTERP_KERNEL::Exception)
         {
           MEDFileParameters *ret=self->getParams();
           if(ret)
             ret->incrRef();
           return ret;
         }

         MEDFileFields *getFields() const throw(INTERP_KERNEL::Exception)
         {
           MEDFileFields *ret=self->getFields();
           if(ret)
             ret->incrRef();
           return ret;
         }

         bool changeMeshNames(PyObject *li) throw(INTERP_KERNEL::Exception)
         {
           std::vector< std::pair<std::string,std::string> > modifTab=convertVecPairStStFromPy(li);
           return self->changeMeshNames(modifTab);
         }
       }
  };

  class SauvReader : public RefCountObject
  {
  public:
    static SauvReader* New(const std::string& fileName) throw(INTERP_KERNEL::Exception);
    MEDFileData * loadInMEDFileDS() throw(INTERP_KERNEL::Exception);
    %extend
    {
      SauvReader(const std::string& fileName) throw(INTERP_KERNEL::Exception)
      {
        return SauvReader::New(fileName);
      }
    }
  };

  class SauvWriter : public RefCountObject
  {
  public:
    static SauvWriter * New();
    void setMEDFileDS(const MEDFileData* medData, unsigned meshIndex = 0) throw(INTERP_KERNEL::Exception);
    void write(const std::string& fileName) throw(INTERP_KERNEL::Exception);
    void setCpyGrpIfOnASingleFamilyStatus(bool status) throw(INTERP_KERNEL::Exception);
    bool getCpyGrpIfOnASingleFamilyStatus() const throw(INTERP_KERNEL::Exception);
    %extend
    {
      SauvWriter() throw(INTERP_KERNEL::Exception)
      {
        return SauvWriter::New();
      }
    }
  };
  
  ///////////////

  class MEDFileMeshStruct;

  class MEDFileField1TSStructItem
  {
  public:
    static MEDFileField1TSStructItem BuildItemFrom(const MEDFileAnyTypeField1TS *ref, const MEDFileMeshStruct *meshSt) throw(INTERP_KERNEL::Exception);
  };

  class MEDFileMeshStruct : public RefCountObject
  {
  public:
    static MEDFileMeshStruct *New(const MEDFileMesh *mesh) throw(INTERP_KERNEL::Exception);
  protected:
    ~MEDFileMeshStruct();
  };
  
  class MEDMeshMultiLev : public RefCountObject
  {
  public:
    virtual MEDMeshMultiLev *prepare() const throw(INTERP_KERNEL::Exception);
    DataArray *buildDataArray(const MEDFileField1TSStructItem& fst, const MEDFileFieldGlobsReal *globs, const DataArray *vals) const throw(INTERP_KERNEL::Exception);
  protected:
    ~MEDMeshMultiLev();
  public:
    %extend
    {
      PyObject *retrieveFamilyIdsOnCells() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *famIds(0);
        bool isWithoutCopy(false);
        self->retrieveFamilyIdsOnCells(famIds,isWithoutCopy);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret1Py=isWithoutCopy?Py_True:Py_False;
        Py_XINCREF(ret1Py);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(famIds),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,ret1Py);
        return ret;
      }

      PyObject *retrieveNumberIdsOnCells() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *numIds(0);
        bool isWithoutCopy(false);
        self->retrieveNumberIdsOnCells(numIds,isWithoutCopy);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret1Py=isWithoutCopy?Py_True:Py_False;
        Py_XINCREF(ret1Py);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(numIds),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,ret1Py);
        return ret;
      }
      
      PyObject *retrieveFamilyIdsOnNodes() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *famIds(0);
        bool isWithoutCopy(false);
        self->retrieveFamilyIdsOnNodes(famIds,isWithoutCopy);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret1Py=isWithoutCopy?Py_True:Py_False;
        Py_XINCREF(ret1Py);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(famIds),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,ret1Py);
        return ret;
      }

      PyObject *retrieveNumberIdsOnNodes() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *numIds(0);
        bool isWithoutCopy(false);
        self->retrieveNumberIdsOnNodes(numIds,isWithoutCopy);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret1Py=isWithoutCopy?Py_True:Py_False;
        Py_XINCREF(ret1Py);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(numIds),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,ret1Py);
        return ret;
      }

      PyObject *getGeoTypes() const throw(INTERP_KERNEL::Exception)
      {
        std::vector< INTERP_KERNEL::NormalizedCellType > result(self->getGeoTypes());
        std::vector< INTERP_KERNEL::NormalizedCellType >::const_iterator iL(result.begin());
        PyObject *res(PyList_New(result.size()));
        for(int i=0;iL!=result.end(); i++, iL++)
          PyList_SetItem(res,i,PyInt_FromLong(*iL));
        return res;
      }
    }
  };

  class MEDUMeshMultiLev : public MEDMeshMultiLev
  {
  protected:
    ~MEDUMeshMultiLev();
  public:
    %extend
     {
       PyObject *buildVTUArrays() const throw(INTERP_KERNEL::Exception)
       {
         DataArrayDouble *coords(0); DataArrayByte *types(0); DataArrayInt *cellLocations(0),*cells(0),*faceLocations(0),*faces(0);
         bool ncc(self->buildVTUArrays(coords,types,cellLocations,cells,faceLocations,faces));
         PyObject *ret0Py=ncc?Py_True:Py_False;
         Py_XINCREF(ret0Py);
         PyObject *ret=PyTuple_New(7);
         PyTuple_SetItem(ret,0,ret0Py);
         PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(coords),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
         PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(types),SWIGTYPE_p_ParaMEDMEM__DataArrayByte, SWIG_POINTER_OWN | 0 ));
         PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(cellLocations),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
         PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(cells),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
         PyTuple_SetItem(ret,5,SWIG_NewPointerObj(SWIG_as_voidptr(faceLocations),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
         PyTuple_SetItem(ret,6,SWIG_NewPointerObj(SWIG_as_voidptr(faces),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
         return ret;
       }
     }
  };

  class MEDStructuredMeshMultiLev : public MEDMeshMultiLev
  {
  protected:
    ~MEDStructuredMeshMultiLev();
  };

  class MEDCMeshMultiLev : public MEDStructuredMeshMultiLev
  {
  protected:
    ~MEDCMeshMultiLev();
  public:
    %extend
    {
      PyObject *buildVTUArrays() const throw(INTERP_KERNEL::Exception)
      {
        bool isInternal;
        std::vector< DataArrayDouble * > objs(self->buildVTUArrays(isInternal));
        std::size_t sz(objs.size());
        PyObject *ret(PyTuple_New(2));
        PyObject *ret0=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(ret0,i,SWIG_NewPointerObj(SWIG_as_voidptr(objs[i]),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,0,ret0);
        PyObject *ret1Py(isInternal?Py_True:Py_False);
        Py_XINCREF(ret1Py);
        PyTuple_SetItem(ret,1,ret1Py);
        return ret;
      }
    }
  };

  class MEDCurveLinearMeshMultiLev : public MEDStructuredMeshMultiLev
  {
  protected:
    ~MEDCurveLinearMeshMultiLev();
  public:
    %extend
    {
      PyObject *buildVTUArrays() const throw(INTERP_KERNEL::Exception)
      {
        DataArrayDouble *ret0(0);
        std::vector<int> ret1;
        bool ret2;
        self->buildVTUArrays(ret0,ret1,ret2);
        std::size_t sz(ret1.size());
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
        PyObject *ret1Py=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(ret1Py,i,SWIG_From_int(ret1[i]));
        PyTuple_SetItem(ret,1,ret1Py);
        PyObject *ret2Py(ret2?Py_True:Py_False);
        Py_XINCREF(ret2Py);
        PyTuple_SetItem(ret,2,ret2Py);
        return ret;
      }
    }
  };

  class MEDFileFastCellSupportComparator : public RefCountObject
  {
  public:
    static MEDFileFastCellSupportComparator *New(const MEDFileMeshStruct *m, const MEDFileAnyTypeFieldMultiTS *ref) throw(INTERP_KERNEL::Exception);
    MEDMeshMultiLev *buildFromScratchDataSetSupport(int timeStepId, const MEDFileFieldGlobsReal *globs) const throw(INTERP_KERNEL::Exception);
    bool isDataSetSupportEqualToThePreviousOne(int timeStepId, const MEDFileFieldGlobsReal *globs) const throw(INTERP_KERNEL::Exception);
    int getNumberOfTS() const throw(INTERP_KERNEL::Exception);
  protected:
    ~MEDFileFastCellSupportComparator();
  public:
    %extend
    {
      PyObject *getGeoTypesAt(int timeStepId, const MEDFileMesh *m) const throw(INTERP_KERNEL::Exception)
      {
        std::vector< INTERP_KERNEL::NormalizedCellType > result(self->getGeoTypesAt(timeStepId,m));
        std::vector< INTERP_KERNEL::NormalizedCellType >::const_iterator iL(result.begin());
        PyObject *res(PyList_New(result.size()));
        for(int i=0;iL!=result.end(); i++, iL++)
          PyList_SetItem(res,i,PyInt_FromLong(*iL));
        return res;
      }
    }
  };
}
