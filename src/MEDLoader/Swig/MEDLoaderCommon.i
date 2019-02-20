// Copyright (C) 2017-2019  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (EDF R&D)

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
#include "MEDFileEquivalence.hxx"
#include "MEDFileEntities.hxx"
#include "MEDFileMeshReadSelector.hxx"
#include "MEDFileFieldOverView.hxx"
#include "MEDLoaderTypemaps.i"
#include "SauvReader.hxx"
#include "SauvWriter.hxx"

using namespace MEDCoupling;
%}

#if SWIG_VERSION >= 0x010329
%template()  std::vector<std::string>;
#endif

%typemap(out) MEDCoupling::MEDFileMesh*
{
  $result=convertMEDFileMesh($1,$owner);
}

%typemap(out) MEDCoupling::MEDFileParameter1TS*
{
  $result=convertMEDFileParameter1TS($1,$owner);
}

%typemap(out) MEDCoupling::MEDFileAnyTypeFieldMultiTS*
{
  $result=convertMEDFileFieldMultiTS($1,$owner);
}

%typemap(out) MEDCoupling::MEDFileAnyTypeField1TS*
{
  $result=convertMEDFileField1TS($1,$owner);
}

%typemap(out) MEDCoupling::MEDMeshMultiLev*
{
  $result=convertMEDMeshMultiLev($1,$owner);
}

%newobject ReadUMeshFromFamiliesSwig;
%newobject ReadUMeshFromGroupsSwig;
%newobject ReadFieldSwig;
%newobject MEDCoupling::ReadUMeshFromFile;
%newobject MEDCoupling::ReadMeshFromFile;
%newobject MEDCoupling::ReadFieldCell;
%newobject MEDCoupling::ReadFieldNode;
%newobject MEDCoupling::ReadFieldGauss;
%newobject MEDCoupling::ReadFieldGaussNE;
%newobject MEDCoupling::MEDFileMesh::New;
%newobject MEDCoupling::MEDFileMesh::createNewEmpty;
%newobject MEDCoupling::MEDFileMesh::deepCopy;
%newobject MEDCoupling::MEDFileMesh::shallowCpy;
%newobject MEDCoupling::MEDFileMesh::getMeshAtLevel;
%newobject MEDCoupling::MEDFileMesh::__getitem__;
%newobject MEDCoupling::MEDFileMesh::getGroupArr;
%newobject MEDCoupling::MEDFileMesh::getGroupsArr;
%newobject MEDCoupling::MEDFileMesh::getFamilyArr;
%newobject MEDCoupling::MEDFileMesh::getFamiliesArr;
%newobject MEDCoupling::MEDFileMesh::getNodeGroupArr;
%newobject MEDCoupling::MEDFileMesh::getNodeGroupsArr;
%newobject MEDCoupling::MEDFileMesh::getNodeFamilyArr;
%newobject MEDCoupling::MEDFileMesh::getNodeFamiliesArr;
%newobject MEDCoupling::MEDFileMesh::getGlobalNumFieldAtLevel;
%newobject MEDCoupling::MEDFileMesh::getAllFamiliesIdsReferenced;
%newobject MEDCoupling::MEDFileMesh::computeAllFamilyIdsInUse;
%newobject MEDCoupling::MEDFileMesh::getEquivalences;
%newobject MEDCoupling::MEDFileMesh::cartesianize;
%newobject MEDCoupling::MEDFileData::getJoints;
%newobject MEDCoupling::MEDFileStructuredMesh::getImplicitFaceMesh;
%newobject MEDCoupling::MEDFileUMesh::New;
%newobject MEDCoupling::MEDFileUMesh::LoadPartOf;
%newobject MEDCoupling::MEDFileUMesh::getCoords;
%newobject MEDCoupling::MEDFileUMesh::getPartDefAtLevel;
%newobject MEDCoupling::MEDFileUMesh::getGroup;
%newobject MEDCoupling::MEDFileUMesh::getGroups;
%newobject MEDCoupling::MEDFileUMesh::getFamily;
%newobject MEDCoupling::MEDFileUMesh::getFamilies;
%newobject MEDCoupling::MEDFileUMesh::getLevel0Mesh;
%newobject MEDCoupling::MEDFileUMesh::getLevelM1Mesh;
%newobject MEDCoupling::MEDFileUMesh::getLevelM2Mesh;
%newobject MEDCoupling::MEDFileUMesh::getLevelM3Mesh;
%newobject MEDCoupling::MEDFileUMesh::getDirectUndergroundSingleGeoTypeMesh;
%newobject MEDCoupling::MEDFileUMesh::extractFamilyFieldOnGeoType;
%newobject MEDCoupling::MEDFileUMesh::extractNumberFieldOnGeoType;
%newobject MEDCoupling::MEDFileUMesh::zipCoords;
%newobject MEDCoupling::MEDFileUMesh::deduceNodeSubPartFromCellSubPart;
%newobject MEDCoupling::MEDFileUMesh::extractPart;
%newobject MEDCoupling::MEDFileUMesh::buildExtrudedMesh;
%newobject MEDCoupling::MEDFileUMesh::linearToQuadratic;
%newobject MEDCoupling::MEDFileUMesh::quadraticToLinear;
%newobject MEDCoupling::MEDFileUMesh::symmetry3DPlane;
%newobject MEDCoupling::MEDFileUMesh::Aggregate;
%newobject MEDCoupling::MEDFileUMesh::convertToExtrudedMesh;
%newobject MEDCoupling::MEDFileCMesh::New;
%newobject MEDCoupling::MEDFileCurveLinearMesh::New;
%newobject MEDCoupling::MEDFileMeshMultiTS::New;
%newobject MEDCoupling::MEDFileMeshMultiTS::deepCopy;
%newobject MEDCoupling::MEDFileMeshMultiTS::getOneTimeStep;
%newobject MEDCoupling::MEDFileMeshes::New;
%newobject MEDCoupling::MEDFileMeshes::deepCopy;
%newobject MEDCoupling::MEDFileMeshes::getMeshAtPos;
%newobject MEDCoupling::MEDFileMeshes::getMeshWithName;
%newobject MEDCoupling::MEDFileMeshes::__getitem__;
%newobject MEDCoupling::MEDFileMeshes::__iter__;

%newobject MEDCoupling::MEDFileMeshSupports::New;
%newobject MEDCoupling::MEDFileMeshSupports::getSupMeshWithName;

%newobject MEDCoupling::MEDFileStructureElements::New;

%newobject MEDCoupling::MEDFileFields::New;
%newobject MEDCoupling::MEDFileFields::NewAdv;
%newobject MEDCoupling::MEDFileFields::NewWithDynGT;
%newobject MEDCoupling::MEDFileFields::LoadPartOf;
%newobject MEDCoupling::MEDFileFields::LoadSpecificEntities;
%newobject MEDCoupling::MEDFileFields::deepCopy;
%newobject MEDCoupling::MEDFileFields::shallowCpy;
%newobject MEDCoupling::MEDFileFields::getFieldWithName;
%newobject MEDCoupling::MEDFileFields::getFieldAtPos;
%newobject MEDCoupling::MEDFileFields::partOfThisLyingOnSpecifiedMeshName;
%newobject MEDCoupling::MEDFileFields::partOfThisLyingOnSpecifiedTimeSteps;
%newobject MEDCoupling::MEDFileFields::partOfThisNotLyingOnSpecifiedTimeSteps;
%newobject MEDCoupling::MEDFileFields::partOfThisOnStructureElements;
%newobject MEDCoupling::MEDFileFields::__iter__;
%newobject MEDCoupling::MEDFileFields::extractPart;
%newobject MEDCoupling::MEDFileFields::linearToQuadratic;

%newobject MEDCoupling::MEDFileWritableStandAlone::serialize;
%newobject MEDCoupling::MEDFileAnyTypeFieldMultiTS::New;
%newobject MEDCoupling::MEDFileAnyTypeFieldMultiTS::deepCopy;
%newobject MEDCoupling::MEDFileAnyTypeFieldMultiTS::shallowCpy;
%newobject MEDCoupling::MEDFileAnyTypeFieldMultiTS::getTimeStepAtPos;
%newobject MEDCoupling::MEDFileAnyTypeFieldMultiTS::getTimeStep;
%newobject MEDCoupling::MEDFileAnyTypeFieldMultiTS::getTimeStepGivenTime;
%newobject MEDCoupling::MEDFileAnyTypeFieldMultiTS::__iter__;
%newobject MEDCoupling::MEDFileAnyTypeFieldMultiTS::extractPart;
%newobject MEDCoupling::MEDFileAnyTypeFieldMultiTS::buildNewEmpty;
%newobject MEDCoupling::MEDFileFieldMultiTS::New;
%newobject MEDCoupling::MEDFileFieldMultiTS::LoadSpecificEntities;
%newobject MEDCoupling::MEDFileFieldMultiTS::field;
%newobject MEDCoupling::MEDFileFieldMultiTS::getFieldAtLevel;
%newobject MEDCoupling::MEDFileFieldMultiTS::getFieldAtTopLevel;
%newobject MEDCoupling::MEDFileFieldMultiTS::getFieldOnMeshAtLevel;
%newobject MEDCoupling::MEDFileFieldMultiTS::getFieldAtLevelOld;
%newobject MEDCoupling::MEDFileFieldMultiTS::getUndergroundDataArray;
%newobject MEDCoupling::MEDFileFieldMultiTS::convertToInt;

%newobject MEDCoupling::MEDFileIntFieldMultiTS::New;
%newobject MEDCoupling::MEDFileIntFieldMultiTS::field;
%newobject MEDCoupling::MEDFileIntFieldMultiTS::LoadSpecificEntities;
%newobject MEDCoupling::MEDFileIntFieldMultiTS::getUndergroundDataArray;
%newobject MEDCoupling::MEDFileIntFieldMultiTS::convertToDouble;
%newobject MEDCoupling::MEDFileIntFieldMultiTS::getFieldAtLevel;
%newobject MEDCoupling::MEDFileIntFieldMultiTS::getFieldAtTopLevel;
%newobject MEDCoupling::MEDFileIntFieldMultiTS::getFieldOnMeshAtLevel;
%newobject MEDCoupling::MEDFileIntFieldMultiTS::getFieldAtLevelOld;

%newobject MEDCoupling::MEDFileFloatFieldMultiTS::New;
%newobject MEDCoupling::MEDFileFloatFieldMultiTS::field;
%newobject MEDCoupling::MEDFileFloatFieldMultiTS::LoadSpecificEntities;
%newobject MEDCoupling::MEDFileFloatFieldMultiTS::getUndergroundDataArray;
%newobject MEDCoupling::MEDFileFloatFieldMultiTS::convertToDouble;
%newobject MEDCoupling::MEDFileFloatFieldMultiTS::getFieldAtLevel;
%newobject MEDCoupling::MEDFileFloatFieldMultiTS::getFieldAtTopLevel;
%newobject MEDCoupling::MEDFileFloatFieldMultiTS::getFieldOnMeshAtLevel;
%newobject MEDCoupling::MEDFileFloatFieldMultiTS::getFieldAtLevelOld;

%newobject MEDCoupling::MEDFileAnyTypeField1TS::New;
%newobject MEDCoupling::MEDFileAnyTypeField1TS::NewAdv;
%newobject MEDCoupling::MEDFileAnyTypeField1TS::shallowCpy;
%newobject MEDCoupling::MEDFileAnyTypeField1TS::deepCopy;
%newobject MEDCoupling::MEDFileAnyTypeField1TS::extractPart;
%newobject MEDCoupling::MEDFileField1TS::New;
%newobject MEDCoupling::MEDFileField1TS::field;
%newobject MEDCoupling::MEDFileField1TS::getFieldAtLevel;
%newobject MEDCoupling::MEDFileField1TS::getFieldAtTopLevel;
%newobject MEDCoupling::MEDFileField1TS::getFieldOnMeshAtLevel;
%newobject MEDCoupling::MEDFileField1TS::getFieldAtLevelOld;
%newobject MEDCoupling::MEDFileField1TS::getUndergroundDataArray;
%newobject MEDCoupling::MEDFileField1TS::convertToInt;

%newobject MEDCoupling::MEDFileIntField1TS::New;
%newobject MEDCoupling::MEDFileIntField1TS::field;
%newobject MEDCoupling::MEDFileIntField1TS::getFieldAtLevel;
%newobject MEDCoupling::MEDFileIntField1TS::getFieldAtTopLevel;
%newobject MEDCoupling::MEDFileIntField1TS::getFieldOnMeshAtLevel;
%newobject MEDCoupling::MEDFileIntField1TS::getFieldAtLevelOld;
%newobject MEDCoupling::MEDFileIntField1TS::getUndergroundDataArray;
%newobject MEDCoupling::MEDFileIntField1TS::convertToDouble;

%newobject MEDCoupling::MEDFileFloatField1TS::New;
%newobject MEDCoupling::MEDFileFloatField1TS::field;
%newobject MEDCoupling::MEDFileFloatField1TS::getFieldAtLevel;
%newobject MEDCoupling::MEDFileFloatField1TS::getFieldAtTopLevel;
%newobject MEDCoupling::MEDFileFloatField1TS::getFieldOnMeshAtLevel;
%newobject MEDCoupling::MEDFileFloatField1TS::getFieldAtLevelOld;
%newobject MEDCoupling::MEDFileFloatField1TS::getUndergroundDataArray;
%newobject MEDCoupling::MEDFileFloatField1TS::convertToDouble;

%newobject MEDCoupling::MEDFileData::New;
%newobject MEDCoupling::MEDFileData::deepCopy;
%newobject MEDCoupling::MEDFileData::getMeshes;
%newobject MEDCoupling::MEDFileData::getFields;
%newobject MEDCoupling::MEDFileData::getParams;
%newobject MEDCoupling::MEDFileData::Aggregate;

%newobject MEDCoupling::MEDFileEntities::BuildFrom;

%newobject MEDCoupling::MEDFileParameterDouble1TS::New;
%newobject MEDCoupling::MEDFileParameterDouble1TS::deepCopy;
%newobject MEDCoupling::MEDFileParameterMultiTS::New;
%newobject MEDCoupling::MEDFileParameterMultiTS::deepCopy;
%newobject MEDCoupling::MEDFileParameterMultiTS::getTimeStepAtPos;
%newobject MEDCoupling::MEDFileParameterMultiTS::__getitem__;
%newobject MEDCoupling::MEDFileParameters::New;
%newobject MEDCoupling::MEDFileParameters::deepCopy;
%newobject MEDCoupling::MEDFileParameters::getParamAtPos;
%newobject MEDCoupling::MEDFileParameters::getParamWithName;
%newobject MEDCoupling::MEDFileParameters::__getitem__;

%newobject MEDCoupling::MEDFileJointCorrespondence::New;
%newobject MEDCoupling::MEDFileJointCorrespondence::deepCopy;
%newobject MEDCoupling::MEDFileJointCorrespondence::shallowCpy;
%newobject MEDCoupling::MEDFileJointCorrespondence::getCorrespondence;
%newobject MEDCoupling::MEDFileJointOneStep::New;
%newobject MEDCoupling::MEDFileJointOneStep::deepCopy;
%newobject MEDCoupling::MEDFileJointOneStep::shallowCpy;
%newobject MEDCoupling::MEDFileJointOneStep::getCorrespondenceAtPos;
%newobject MEDCoupling::MEDFileJointOneStep::__getitem__;
%newobject MEDCoupling::MEDFileJoint::New;
%newobject MEDCoupling::MEDFileJoint::deepCopy;
%newobject MEDCoupling::MEDFileJoint::shallowCpy;
%newobject MEDCoupling::MEDFileJoint::getStepAtPos;
%newobject MEDCoupling::MEDFileJoint::__getitem__;
%newobject MEDCoupling::MEDFileJoints::New;
%newobject MEDCoupling::MEDFileJoints::deepCopy;
%newobject MEDCoupling::MEDFileJoints::getJointAtPos;
%newobject MEDCoupling::MEDFileJoints::getJointWithName;
%newobject MEDCoupling::MEDFileJoints::__getitem__;
%newobject MEDCoupling::MEDFileEquivalences::getEquivalence;
%newobject MEDCoupling::MEDFileEquivalences::getEquivalenceWithName;
%newobject MEDCoupling::MEDFileEquivalences::appendEmptyEquivalenceWithName;
%newobject MEDCoupling::MEDFileEquivalencePair::initCell;
%newobject MEDCoupling::MEDFileEquivalencePair::initNode;
%newobject MEDCoupling::MEDFileEquivalencePair::getCell;
%newobject MEDCoupling::MEDFileEquivalencePair::getNode;
%newobject MEDCoupling::MEDFileEquivalenceData::getArray;
%newobject MEDCoupling::MEDFileEquivalenceCell::getArray;

%newobject MEDCoupling::SauvWriter::New;
%newobject MEDCoupling::SauvReader::New;
%newobject MEDCoupling::SauvReader::loadInMEDFileDS;

%newobject MEDCoupling::MEDFileMeshStruct::New;
%newobject MEDCoupling::MEDMeshMultiLev::prepare;
%newobject MEDCoupling::MEDMeshMultiLev::buildDataArray;
%newobject MEDCoupling::MEDMeshMultiLev::retrieveGlobalNodeIdsIfAny;
%newobject MEDCoupling::MEDFileFastCellSupportComparator::New;
%newobject MEDCoupling::MEDFileFastCellSupportComparator::buildFromScratchDataSetSupport;

%feature("unref") MEDFileMesh "$this->decrRef();"
%feature("unref") MEDFileUMesh "$this->decrRef();"
%feature("unref") MEDFileCMesh "$this->decrRef();"
%feature("unref") MEDFileMeshMultiTS "$this->decrRef();"
%feature("unref") MEDFileMeshes "$this->decrRef();"
%feature("unref") MEDFileFieldLoc "$this->decrRef();"
%feature("unref") MEDFileAnyTypeField1TS "$this->decrRef();"
%feature("unref") MEDFileField1TS "$this->decrRef();"
%feature("unref") MEDFileIntField1TS "$this->decrRef();"
%feature("unref") MEDFileFloatField1TS "$this->decrRef();"
%feature("unref") MEDFileAnyTypeFieldMultiTS "$this->decrRef();"
%feature("unref") MEDFileFieldMultiTS "$this->decrRef();"
%feature("unref") MEDFileIntFieldMultiTS "$this->decrRef();"
%feature("unref") MEDFileFloatFieldMultiTS "$this->decrRef();"
%feature("unref") MEDFileMeshSupports "$this->decrRef();"
%feature("unref") MEDFileStructureElements "$this->decrRef();"
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
%feature("unref") MEDFileEquivalences "$this->decrRef();"
%feature("unref") MEDFileEquivalencePair "$this->decrRef();"
%feature("unref") MEDFileEquivalenceBase "$this->decrRef();"
%feature("unref") MEDFileEquivalenceData "$this->decrRef();"
%feature("unref") MEDFileEquivalenceCell "$this->decrRef();"
%feature("unref") MEDFileEquivalenceNode "$this->decrRef();"
%feature("unref") MEDFileData "$this->decrRef();"
%feature("unref") SauvReader "$this->decrRef();"
%feature("unref") SauvWriter "$this->decrRef();"
%feature("unref") MEDFileFastCellSupportComparator "$this->decrRef();"
%feature("unref") MEDMeshMultiLev "$this->decrRef();"
%feature("unref") MEDUMeshMultiLev "$this->decrRef();"
%feature("unref") MEDCMeshMultiLev "$this->decrRef();"
%feature("unref") MEDCurveLinearMeshMultiLev "$this->decrRef();"
%feature("unref") MEDFileMeshStruct "$this->decrRef();"

namespace MEDCoupling
{
  bool HasXDR();
  std::string MEDFileVersionStr();
  std::string MEDFileVersionOfFileStr(const std::string& fileName);
  void SetEpsilonForNodeComp(double val);
  void SetCompPolicyForCell(int val);
  void SetTooLongStrPolicy(int val);
  void CheckFileForRead(const std::string& fileName);
  std::vector<std::string> GetMeshNames(const std::string& fileName);
  std::vector<std::string> GetMeshNamesOnField(const std::string& fileName, const std::string& fieldName);
  std::vector<std::string> GetMeshGroupsNames(const std::string& fileName, const std::string& meshName);
  std::vector<std::string> GetMeshFamiliesNames(const std::string& fileName, const std::string& meshName);
  std::vector<std::string> GetMeshFamiliesNamesOnGroup(const std::string& fileName, const std::string& meshName, const std::string& grpName);
  std::vector<std::string> GetMeshGroupsNamesOnFamily(const std::string& fileName, const std::string& meshName, const std::string& famName);
  std::vector<std::string> GetAllFieldNamesOnMesh(const std::string& fileName, const std::string& meshName);
  std::vector<std::string> GetAllFieldNames(const std::string& fileName);
  std::vector<std::string> GetFieldNamesOnMesh(MEDCoupling::TypeOfField type, const std::string& fileName, const std::string& meshName);
  std::vector<std::string> GetCellFieldNamesOnMesh(const std::string& fileName, const std::string& meshName);
  std::vector<std::string> GetNodeFieldNamesOnMesh(const std::string& fileName, const std::string& meshName);
  double GetTimeAttachedOnFieldIteration(const std::string& fileName, const std::string& fieldName, int iteration, int order);
  void AssignStaticWritePropertiesTo(MEDCoupling::MEDFileWritable& obj);
  MEDCoupling::MEDCouplingMesh *ReadMeshFromFile(const std::string& fileName, const std::string& meshName, int meshDimRelToMax=0);
  MEDCoupling::MEDCouplingMesh *ReadMeshFromFile(const std::string& fileName, int meshDimRelToMax=0);
  MEDCoupling::MEDCouplingUMesh *ReadUMeshFromFile(const std::string& fileName, const std::string& meshName, int meshDimRelToMax=0);
  MEDCoupling::MEDCouplingUMesh *ReadUMeshFromFile(const std::string& fileName, int meshDimRelToMax=0);
  int ReadUMeshDimFromFile(const std::string& fileName, const std::string& meshName);
  MEDCoupling::MEDCouplingField *ReadFieldCell(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order);
  MEDCoupling::MEDCouplingField *ReadFieldNode(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order);
  MEDCoupling::MEDCouplingField *ReadFieldGauss(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order);
  MEDCoupling::MEDCouplingField *ReadFieldGaussNE(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order);
  void WriteMesh(const std::string& fileName, const MEDCoupling::MEDCouplingMesh *mesh, bool writeFromScratch);
  void WriteUMesh(const std::string& fileName, const MEDCoupling::MEDCouplingUMesh *mesh, bool writeFromScratch);
  void WriteUMeshDep(const std::string& fileName, const MEDCoupling::MEDCouplingUMesh *mesh, bool writeFromScratch);
  void WriteField(const std::string& fileName, const MEDCoupling::MEDCouplingField *f, bool writeFromScratch);
  void WriteFieldDep(const std::string& fileName, const MEDCoupling::MEDCouplingField *f, bool writeFromScratch);
  void WriteFieldUsingAlreadyWrittenMesh(const std::string& fileName, const MEDCoupling::MEDCouplingField *f);
}

%rename (MEDFileVersion) MEDFileVersionSwig;
%rename (GetFieldIterations) GetFieldIterationsSwig;
%rename (GetAllFieldIterations) GetAllFieldIterationsSwig;
%rename (GetCellFieldIterations) GetCellFieldIterationsSwig;
%rename (GetNodeFieldIterations) GetNodeFieldIterationsSwig;
%rename (GetComponentsNamesOfField) GetComponentsNamesOfFieldSwig;
%rename (GetUMeshGlobalInfo) GetUMeshGlobalInfoSwig;
%rename (ReadFieldsOnSameMesh) ReadFieldsOnSameMeshSwig;
%rename (WriteUMeshesPartition) WriteUMeshesPartitionSwig;
%rename (WriteUMeshesPartitionDep) WriteUMeshesPartitionDepSwig;
%rename (WriteUMeshes) WriteUMeshesSwig;
%rename (GetTypesOfField) GetTypesOfFieldSwig;
%rename (ReadUMeshFromGroups) ReadUMeshFromGroupsSwig;
%rename (ReadUMeshFromFamilies) ReadUMeshFromFamiliesSwig;
%rename (ReadField) ReadFieldSwig;

%inline
{
  PyObject *MEDFileVersionSwig()
  {
    int major,minor,release;
    MEDCoupling::MEDFileVersion(major,minor,release);
    PyObject *ret(PyTuple_New(3));
    PyTuple_SetItem(ret,0,SWIG_From_int(major));
    PyTuple_SetItem(ret,1,SWIG_From_int(minor));
    PyTuple_SetItem(ret,2,SWIG_From_int(release));
    return ret;
  }

  MEDCoupling::MEDCouplingField *ReadFieldSwig(const std::string& fileName)
  {
    MCAuto<MEDCoupling::MEDCouplingField> ret(MEDCoupling::ReadField(fileName));
    return ret.retn();
  }

  MEDCoupling::MEDCouplingField *ReadFieldSwig(const std::string& fileName, const std::string& fieldName)
  {
    MCAuto<MEDCoupling::MEDCouplingField> ret(MEDCoupling::ReadField(fileName,fieldName));
    return ret.retn();
  }
  
  MEDCoupling::MEDCouplingField *ReadFieldSwig(const std::string& fileName, const std::string& fieldName, int iteration, int order)
  {
    MCAuto<MEDCoupling::MEDCouplingField> ret(MEDCoupling::ReadField(fileName,fieldName,iteration,order));
    return ret.retn();
  }
  
  MEDCoupling::MEDCouplingField *ReadFieldSwig(MEDCoupling::TypeOfField type, const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order)
  {
    MCAuto<MEDCoupling::MEDCouplingField> ret(MEDCoupling::ReadField(type,fileName,meshName,meshDimRelToMax,fieldName,iteration,order));
    return ret.retn();
  }

  PyObject *GetFieldIterationsSwig(MEDCoupling::TypeOfField type, const std::string& fileName, const std::string& meshName, const std::string& fieldName)
  {
    std::vector< std::pair<int,int> > res=MEDCoupling::GetFieldIterations(type,fileName,meshName,fieldName);
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
  
  PyObject *GetAllFieldIterationsSwig(const std::string& fileName, const std::string& fieldName)
    {
      std::vector< std::pair< std::pair<int,int>, double> > res=MEDCoupling::GetAllFieldIterations(fileName,fieldName);
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
  
  PyObject *GetCellFieldIterationsSwig(const std::string& fileName, const std::string& meshName, const std::string& fieldName)
    {
      std::vector< std::pair<int,int> > res=MEDCoupling::GetCellFieldIterations(fileName,meshName,fieldName);
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

  PyObject *GetNodeFieldIterationsSwig(const std::string& fileName, const std::string& meshName, const std::string& fieldName)
    {
      std::vector< std::pair<int,int> > res=MEDCoupling::GetNodeFieldIterations(fileName,meshName,fieldName);
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

  PyObject *GetComponentsNamesOfFieldSwig(const std::string& fileName, const std::string& fieldName)
    {
      std::vector< std::pair<std::string,std::string> > res=MEDCoupling::GetComponentsNamesOfField(fileName,fieldName);
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

  PyObject *GetUMeshGlobalInfoSwig(const std::string& fileName, const std::string& meshName)
    {
      int meshDim,spaceDim,numberOfNodes;
      std::vector< std::vector< std::pair<INTERP_KERNEL::NormalizedCellType,int> > > res=MEDCoupling::GetUMeshGlobalInfo(fileName,meshName,meshDim,spaceDim,numberOfNodes);
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
  
  PyObject *ReadFieldsOnSameMeshSwig(MEDCoupling::TypeOfField type, const std::string& fileName, const std::string& meshName, int meshDimRelToMax,
                                     const std::string& fieldName, PyObject *liIts)
    {
      std::vector<std::pair<int,int> > its=convertTimePairIdsFromPy(liIts);
      std::vector<MEDCoupling::MEDCouplingFieldDouble *> res=MEDCoupling::ReadFieldsOnSameMesh(type,fileName,meshName,meshDimRelToMax,fieldName,its);
      return convertFieldDoubleVecToPy(res);
    }
  
  void WriteUMeshesPartitionSwig(const std::string& fileName, const std::string& meshName, PyObject *li, bool writeFromScratch)
  {
    std::vector<const MEDCoupling::MEDCouplingUMesh *> v;
    convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",v);
    MEDCoupling::WriteUMeshesPartition(fileName,meshName,v,writeFromScratch);
  }
  
  void WriteUMeshesPartitionDepSwig(const std::string& fileName, const std::string& meshName, PyObject *li, bool writeFromScratch)
  {
    std::vector<const MEDCoupling::MEDCouplingUMesh *> v;
    convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",v);
    MEDCoupling::WriteUMeshesPartitionDep(fileName,meshName,v,writeFromScratch);
  }
  
  void WriteUMeshesSwig(const std::string& fileName, PyObject *li, bool writeFromScratch)
  {
    std::vector<const MEDCoupling::MEDCouplingUMesh *> v;
    convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",v);
    MEDCoupling::WriteUMeshes(fileName,v,writeFromScratch);
  }
  
  PyObject *GetTypesOfFieldSwig(const std::string& fileName, const std::string& meshName, const std::string& fieldName)
    {
      std::vector< MEDCoupling::TypeOfField > v=MEDCoupling::GetTypesOfField(fileName,meshName,fieldName);
      int size=v.size();
      PyObject *ret=PyList_New(size);
      for(int i=0;i<size;i++)
        PyList_SetItem(ret,i,PyInt_FromLong((int)v[i]));
      return ret;
    }
  
  MEDCoupling::MEDCouplingUMesh *ReadUMeshFromGroupsSwig(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, PyObject *li)
    {
      std::vector<std::string> grps;
      converPyListToVecString(li,grps);
      return MEDCoupling::ReadUMeshFromGroups(fileName,meshName,meshDimRelToMax,grps);
    }

  MEDCoupling::MEDCouplingUMesh *ReadUMeshFromFamiliesSwig(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, PyObject *li)
    {
      std::vector<std::string> fams;
      converPyListToVecString(li,fams);
      return MEDCoupling::ReadUMeshFromFamilies(fileName,meshName,meshDimRelToMax,fams);
    }
}

namespace MEDCoupling
{
  class MEDFileWritable
  {
  public:
    void copyOptionsFrom(const MEDFileWritable& other) const;
    int getTooLongStrPolicy() const;
    void setTooLongStrPolicy(int newVal);
    int getZipConnPolicy();
    void setZipConnPolicy(int newVal);
  };
  
  class MEDFileWritableStandAlone : public MEDFileWritable
  {
  public:
    void write(const std::string& fileName, int mode) const;
    void write33(const std::string& fileName, int mode) const;
    %extend
       {
         DataArrayByte *serialize() const
         {
           MCAuto<DataArrayByte> ret(self->serialize());
           return ret.retn();
         }

         PyObject *__getstate__()
         {
           PyObject *ret(PyList_New(0));
           return ret;
         }

         void __setstate__(PyObject *inp)
         {
         }
       }
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
    bool isGlobalNodeNumFieldReading() const;
    void setCellFamilyFieldReading(bool b);
    void setNodeFamilyFieldReading(bool b);
    void setCellNameFieldReading(bool b);
    void setNodeNameFieldReading(bool b);
    void setCellNumFieldReading(bool b);
    void setNodeNumFieldReading(bool b);
    void setGlobalNodeNumFieldReading(bool b);
    %extend
    {
      std::string __str__() const
      {
        std::ostringstream oss;
        self->reprAll(oss);
        return oss.str();
      }
      
      std::string __repr__() const
      {
        std::ostringstream oss; oss << "MEDFileMeshReadSelector C++ instance at " << self << " (with code=" << self->getCode() << ").";
        return oss.str();
      }
    }
  };

  class MEDFileJointCorrespondence : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileJointCorrespondence *New();
    static MEDFileJointCorrespondence *New(DataArrayInt* correspondence) // nodes
     ;
    static MEDFileJointCorrespondence *New(DataArrayInt* correspondence,  // cells
                                           INTERP_KERNEL::NormalizedCellType loc_geo_type,
                                           INTERP_KERNEL::NormalizedCellType rem_geo_type)
     ;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDFileJointCorrespondence *deepCopy() const;
    MEDFileJointCorrespondence *shallowCpy() const;
    void setIsNodal(bool isNodal);
    bool getIsNodal() const;
    bool isEqual(const MEDFileJointCorrespondence *other) const;
    void setLocalGeometryType(INTERP_KERNEL::NormalizedCellType type);
    INTERP_KERNEL::NormalizedCellType getLocalGeometryType() const;
    void setRemoteGeometryType(INTERP_KERNEL::NormalizedCellType type);
    INTERP_KERNEL::NormalizedCellType getRemoteGeometryType() const;
    void setCorrespondence(DataArrayInt *corr);
    void write(const std::string& fileName, int mode, const std::string& localMeshName, const std::string& jointName, int order, int iteration) const;
    std::string simpleRepr() const;
    %extend
    {
      MEDFileJointCorrespondence()
      {
        return MEDFileJointCorrespondence::New();
      }
      MEDFileJointCorrespondence(DataArrayInt* correspondence)
      {
        return MEDFileJointCorrespondence::New(correspondence);
      }
      MEDFileJointCorrespondence(DataArrayInt* correspondence,  // cells
                                 INTERP_KERNEL::NormalizedCellType loc_geo_type,
                                 INTERP_KERNEL::NormalizedCellType rem_geo_type)
      {
        return MEDFileJointCorrespondence::New(correspondence, loc_geo_type, rem_geo_type);
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }
      
      DataArrayInt *getCorrespondence() const
      {
        const DataArrayInt *ret(self->getCorrespondence());
        if(ret)
          ret->incrRef();
        return const_cast<DataArrayInt *>(ret);
      }
    }
  };

  class MEDFileJointOneStep : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileJointOneStep *New(int dt=-1, int it=-1);
    static MEDFileJointOneStep *New(const std::string& fileName, const std::string& mName, const std::string& jointName, int number=1);
    MEDFileJointOneStep *deepCopy() const;
    MEDFileJointOneStep *shallowCpy() const;
    bool isEqual(const MEDFileJointOneStep *other) const;
    void setOrder(int order);
    int getOrder() const;
    void setIteration(int it);
    int getIteration() const;
    void pushCorrespondence(MEDFileJointCorrespondence* correspondence);
    int getNumberOfCorrespondences() const;
    void write(const std::string& fileName, int mode, const std::string& localMeshName, const std::string& jointName) const;
    std::string simpleRepr() const;
    %extend
    {
      MEDFileJointOneStep()
      {
        return MEDFileJointOneStep::New();
      }

      MEDFileJointOneStep(const std::string& fileName, const std::string& mName, const std::string& jointName, int number)
      {
        return MEDFileJointOneStep::New(fileName,mName,jointName,number);
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }
      
      MEDFileJointCorrespondence *getCorrespondenceAtPos(int i) const
      {
        MEDFileJointCorrespondence *ret(self->getCorrespondenceAtPos(i));
        if(ret)
          ret->incrRef();
        return ret;
      }

      MEDFileJointCorrespondence *__getitem__(int i) const
      {
        return MEDCoupling_MEDFileJointOneStep_getCorrespondenceAtPos(self,i);
      }
    }
  };

  class MEDFileJoint : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    static MEDFileJoint *New();
    static MEDFileJoint *New(const std::string& fileName, const std::string& mName, int num);
    static MEDFileJoint *New(const std::string& jointName, const std::string& locMeshName, const std::string& remoteMeshName, int remoteMeshNum);
    MEDFileJoint *deepCopy() const;
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
    bool changeJointNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    void setDomainNumber(const int& number);
    int getDomainNumber() const;
    void pushStep(MEDFileJointOneStep* step);
    int getNumberOfSteps() const;
    std::string simpleRepr() const;
    %extend
    {
      MEDFileJoint()
      {
        return MEDFileJoint::New();
      }
      
      MEDFileJoint(const std::string& fileName, const std::string& mName, int num)
      {
        return MEDFileJoint::New(fileName,mName,num);
      }

      MEDFileJoint(const std::string& jointName, const std::string& locMeshName, const std::string& remoteMeshName, int remoteMeshNum)
      {
        return MEDFileJoint::New(jointName,locMeshName,remoteMeshName,remoteMeshNum);
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }
      
      MEDFileJointOneStep *getStepAtPos(int i) const
      {
        MEDFileJointOneStep *ret(self->getStepAtPos(i));
        if(ret)
          ret->incrRef();
        return ret;
      }

      MEDFileJointOneStep *__getitem__(int i)
      {
        return MEDCoupling_MEDFileJoint_getStepAtPos(self,i);
      }
    }
  };

  class MEDFileJoints : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    static MEDFileJoints *New();
    static MEDFileJoints *New(const std::string& fileName, const std::string& meshName);
    MEDFileJoints *deepCopy() const;
    std::string simpleRepr() const;
    std::string getMeshName() const;
    int getNumberOfJoints() const;
    std::vector<std::string> getJointsNames() const;
    bool changeJointNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    void resize(int newSize);
    void pushJoint(MEDFileJoint *joint);
    void setJointAtPos(int i, MEDFileJoint *joint);
    void destroyJointAtPos(int i);
    %extend
    {
      MEDFileJoints()
      {
        return MEDFileJoints::New();
      }
      
      MEDFileJoints(const std::string& fileName, const std::string& meshName)
      {
        return MEDFileJoints::New(fileName,meshName);
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }

      MEDFileJoint *__getitem__(PyObject *obj)
      {
        static const char msg[]="MEDFileJoints::__getitem__ : only integer or string with meshname supported !";
        if(PyInt_Check(obj))
          {
            MEDFileJoint *ret=self->getJointAtPos(InterpreteNegativeInt((int)PyInt_AS_LONG(obj),self->getNumberOfJoints()));
            if(ret)
              ret->incrRef();
            return ret;
          }
        MEDFileJoint *ret(self->getJointWithName(convertPyObjectToStr(obj,msg)));
        if(ret)
          ret->incrRef();
        return ret;
      }

      int __len__() const
      {
        return self->getNumberOfJoints();
      }

      MEDFileJoint *getJointAtPos(int i) const
      {
        MEDFileJoint *ret=self->getJointAtPos(i);
        if(ret)
          ret->incrRef();
        return ret;
      }

      MEDFileJoint *getJointWithName(const std::string& paramName) const
      {
        MEDFileJoint *ret=self->getJointWithName(paramName);
        if(ret)
          ret->incrRef();
        return ret;
      }
    }
  };
  
  class MEDFileEquivalenceBase : public RefCountObject, public MEDFileWritableStandAlone
  {
  private:
    MEDFileEquivalenceBase();
  };

  class MEDFileEquivalenceData : public MEDFileEquivalenceBase
  {
  private:
    MEDFileEquivalenceData();
  public:
    void setArray(DataArrayInt *data);
    %extend
    {
      DataArrayInt *getArray()
      {
        DataArrayInt *ret(self->getArray());
        if(ret) ret->incrRef();
        return ret;
      }
    }
  };

  class MEDFileEquivalenceNode : public MEDFileEquivalenceData
  {
  private:
    MEDFileEquivalenceNode();
  };

  class MEDFileEquivalenceCell : public MEDFileEquivalenceBase
  {
  private:
    MEDFileEquivalenceCell();
  public:
    void clear();
    std::size_t size() const;
    void setArray(int meshDimRelToMax, DataArrayInt *da);
    void setArrayForType(INTERP_KERNEL::NormalizedCellType type, DataArrayInt *da);
    %extend
    {
      DataArrayInt *getArray(INTERP_KERNEL::NormalizedCellType type)
      {
        DataArrayInt *ret(self->getArray(type));
        if(ret) ret->incrRef();
        return ret;
      }
      
      PyObject *getTypes() const
      {
        std::vector<INTERP_KERNEL::NormalizedCellType> result(self->getTypes());
        std::vector<INTERP_KERNEL::NormalizedCellType>::const_iterator iL=result.begin();
        PyObject *res=PyList_New(result.size());
        for(int i=0;iL!=result.end(); i++, iL++)
          PyList_SetItem(res,i,PyInt_FromLong(*iL));
        return res;
      }
    }
  };

  class MEDFileEquivalencePair : public RefCountObject, public MEDFileWritableStandAlone
  {
  private:
    MEDFileEquivalencePair();
  public:
    std::string getName() const;
    void setName(const std::string& name);
    std::string getDescription() const;
    void setDescription(const std::string& descr);
    void setArray(int meshDimRelToMaxExt, DataArrayInt *da);;
    %extend
    {
      MEDFileEquivalenceCell *initCell()
      {
        MEDFileEquivalenceCell *ret(self->initCell());
        if(ret) ret->incrRef();
        return ret;
      }

      MEDFileEquivalenceNode *initNode()
      {
        MEDFileEquivalenceNode *ret(self->initNode());
        if(ret) ret->incrRef();
        return ret;
      }
      
      MEDFileEquivalenceCell *getCell()
      {
        MEDFileEquivalenceCell *ret(self->getCell());
        if(ret) ret->incrRef();
        return ret;
      }
      
      MEDFileEquivalenceNode *getNode()
      {
        MEDFileEquivalenceNode *ret(self->getNode());
        if(ret) ret->incrRef();
        return ret;
      }
    }
  };
  
  class MEDFileEquivalences : public RefCountObject, public MEDFileWritableStandAlone
  {
  private:
    MEDFileEquivalences();
  public:
    int size() const;
    std::vector<std::string> getEquivalenceNames() const;
    void killEquivalenceWithName(const std::string& name);
    void killEquivalenceAt(int i);
    void clear();
    %extend
    {
      MEDFileEquivalencePair *getEquivalence(int i)
      {
        MEDFileEquivalencePair *ret(self->getEquivalence(i));
        if(ret) ret->incrRef();
        return ret;
      }
      MEDFileEquivalencePair *getEquivalenceWithName(const std::string& name)
      {
        MEDFileEquivalencePair *ret(self->getEquivalenceWithName(name));
        if(ret) ret->incrRef();
        return ret;
      }

      MEDFileEquivalencePair *appendEmptyEquivalenceWithName(const std::string& name)
      {
        MEDFileEquivalencePair *ret(self->appendEmptyEquivalenceWithName(name));
        if(ret) ret->incrRef();
        return ret;
      }
    }
  };

  class MEDFileMesh : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    static MEDFileMesh *New(const std::string& fileName, MEDFileMeshReadSelector *mrs=0);
    static MEDFileMesh *New(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    static MEDFileMesh *New(DataArrayByte *db);
    virtual MEDFileMesh *createNewEmpty() const;
    virtual MEDFileMesh *deepCopy() const;
    virtual MEDFileMesh *shallowCpy() const;
    virtual void clearNonDiscrAttributes() const;
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
    void setAxisType(MEDCouplingAxisType at);
    MEDCouplingAxisType getAxisType() const;
    virtual int getNumberOfNodes() const;
    virtual int getNumberOfCellsAtLevel(int meshDimRelToMaxExt) const;
    virtual bool hasImplicitPart() const;
    virtual int buildImplicitPartIfAny(INTERP_KERNEL::NormalizedCellType gt) const;
    virtual void releaseImplicitPartIfAny() const;
    virtual int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType ct) const;
    virtual std::vector<int> getFamArrNonEmptyLevelsExt() const;
    virtual std::vector<int> getNumArrNonEmptyLevelsExt() const;
    virtual std::vector<int> getNameArrNonEmptyLevelsExt() const;
    virtual std::vector<int> getDistributionOfTypes(int meshDimRelToMax) const;
    virtual MEDFileMesh *cartesianize() const;
    std::vector<int> getNonEmptyLevels() const;
    std::vector<int> getNonEmptyLevelsExt() const;
    int getSizeAtLevel(int meshDimRelToMaxExt) const;
    //
    bool existsGroup(const std::string& groupName) const;
    bool existsFamily(int famId) const;
    bool existsFamily(const std::string& familyName) const;
    void setFamilyId(const std::string& familyName, int id);
    void setFamilyIdUnique(const std::string& familyName, int id);
    void addFamily(const std::string& familyName, int id);
    void addFamilyOnGrp(const std::string& grpName, const std::string& famName);
    virtual void createGroupOnAll(int meshDimRelToMaxExt, const std::string& groupName);
    virtual bool keepFamIdsOnlyOnLevs(const std::vector<int>& famIds, const std::vector<int>& levs);
    void copyFamGrpMapsFrom(const MEDFileMesh& other);
    void clearGrpMap();
    void clearFamMap();
    void clearFamGrpMaps();
    const std::map<std::string,int>& getFamilyInfo() const;
    const std::map<std::string, std::vector<std::string> >& getGroupInfo() const;
    std::vector<std::string> getFamiliesOnGroup(const std::string& name) const;
    std::vector<std::string> getFamiliesOnGroups(const std::vector<std::string>& grps) const;
    std::vector<int> getFamiliesIdsOnGroup(const std::string& name) const;
    void setFamiliesOnGroup(const std::string& name, const std::vector<std::string>& fams);
    void setFamiliesIdsOnGroup(const std::string& name, const std::vector<int>& famIds);
    std::vector<std::string> getGroupsOnFamily(const std::string& name) const;
    void setGroupsOnFamily(const std::string& famName, const std::vector<std::string>& grps);
    std::vector<std::string> getGroupsNames() const;
    std::vector<std::string> getFamiliesNames() const;
    std::vector<std::string> getGroupsOnSpecifiedLev(int meshDimRelToMaxExt) const;
    std::vector<int> getGrpNonEmptyLevelsExt(const std::string& grp) const;
    std::vector<int> getGrpNonEmptyLevels(const std::string& grp) const;
    std::vector<int> getFamsNonEmptyLevels(const std::vector<std::string>& fams) const;
    std::vector<int> getFamsNonEmptyLevelsExt(const std::vector<std::string>& fams) const;
    std::vector<int> getGrpsNonEmptyLevels(const std::vector<std::string>& grps) const;
    std::vector<int> getGrpsNonEmptyLevelsExt(const std::vector<std::string>& grps) const;
    std::vector<int> getFamNonEmptyLevels(const std::string& fam) const;
    std::vector<int> getFamNonEmptyLevelsExt(const std::string& fam) const;
    std::vector<std::string> getFamiliesNamesWithFilePointOfView() const;
    static std::string GetMagicFamilyStr();
    void assignFamilyNameWithGroupName();
    std::vector<std::string> removeEmptyGroups();
    void removeGroupAtLevel(int meshDimRelToMaxExt, const std::string& name);
    void removeGroup(const std::string& name);
    void removeFamily(const std::string& name);
    std::vector<std::string> removeOrphanGroups();
    std::vector<std::string> removeOrphanFamilies();
    void removeFamiliesReferedByNoGroups();
    void rearrangeFamilies();
    void checkOrphanFamilyZero() const;
    void changeGroupName(const std::string& oldName, const std::string& newName);
    void changeFamilyName(const std::string& oldName, const std::string& newName);
    void changeFamilyId(int oldId, int newId);
    void changeAllGroupsContainingFamily(const std::string& familyNameToChange, const std::vector<std::string>& newFamiliesNames);
    void setFamilyInfo(const std::map<std::string,int>& info);
    void setGroupInfo(const std::map<std::string, std::vector<std::string> >&info);
    int getFamilyId(const std::string& name) const;
    int getMaxAbsFamilyId() const;
    int getMaxFamilyId() const;
    int getMinFamilyId() const;
    int getTheMaxAbsFamilyId() const;
    int getTheMaxFamilyId() const;
    int getTheMinFamilyId() const;
    virtual int getMaxAbsFamilyIdInArrays() const;
    virtual int getMaxFamilyIdInArrays() const;
    virtual int getMinFamilyIdInArrays() const;
    DataArrayInt *getAllFamiliesIdsReferenced() const;
    DataArrayInt *computeAllFamilyIdsInUse() const;
    std::vector<int> getFamiliesIds(const std::vector<std::string>& famNames) const;
    std::string getFamilyNameGivenId(int id) const;
    bool ensureDifferentFamIdsPerLevel();
    void normalizeFamIdsTrio();
    void normalizeFamIdsMEDFile();
    virtual int getMeshDimension() const;
    virtual std::string simpleRepr() const;
    virtual std::string advancedRepr() const;
    //
    virtual MEDCouplingMesh *getMeshAtLevel(int meshDimRelToMax, bool renum=false) const;
    virtual void setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr);
    virtual void setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr);
    virtual void setNameFieldAtLevel(int meshDimRelToMaxExt, DataArrayAsciiChar *nameArr);
    virtual void setGlobalNumFieldAtLevel(int meshDimRelToMaxExt, DataArrayInt *globalNumArr);
    virtual void addNodeGroup(const DataArrayInt *ids);
    virtual void addGroup(int meshDimRelToMaxExt, const DataArrayInt *ids);
    virtual DataArrayInt *getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const;
    virtual DataArrayInt *getGroupsArr(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum=false) const;
    virtual DataArrayInt *getGroupArr(int meshDimRelToMaxExt, const std::string& grp, bool renum=false) const;
    virtual DataArrayInt *getFamilyArr(int meshDimRelToMaxExt, const std::string& fam, bool renum=false) const;
    virtual DataArrayInt *getNodeGroupArr(const std::string& grp, bool renum=false) const;
    virtual DataArrayInt *getNodeGroupsArr(const std::vector<std::string>& grps, bool renum=false) const;
    virtual DataArrayInt *getNodeFamilyArr(const std::string& fam, bool renum=false) const;
    virtual DataArrayInt *getNodeFamiliesArr(const std::vector<std::string>& fams, bool renum=false) const;
    int getNumberOfJoints();
    MEDFileJoints *getJoints();
    void setJoints( MEDFileJoints* joints );
    void initializeEquivalences();
    void killEquivalences();
    bool presenceOfStructureElements() const;
    void killStructureElements();
    %extend
       {
         std::string __str__() const
         {
           return self->simpleRepr();
         }

         MEDCouplingMesh *__getitem__(int meshDimRelToMaxExt) const
         {
           return self->getMeshAtLevel(meshDimRelToMaxExt,false);
         }

         PyObject *getTime()
         {
           int tmp1,tmp2;
           double tmp0=self->getTime(tmp1,tmp2);
           PyObject *res = PyList_New(3);
           PyList_SetItem(res,0,SWIG_From_int(tmp1));
           PyList_SetItem(res,1,SWIG_From_int(tmp2));
           PyList_SetItem(res,2,SWIG_From_double(tmp0));
           return res;
         }

         virtual PyObject *isEqual(const MEDFileMesh *other, double eps) const
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

         void setGroupsAtLevel(int meshDimRelToMaxExt, PyObject *li, bool renum=false)
         {
           std::vector<const DataArrayInt *> grps;
           convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayInt *>(li,SWIGTYPE_p_MEDCoupling__DataArrayInt,"DataArrayInt",grps);
           self->setGroupsAtLevel(meshDimRelToMaxExt,grps,renum);
         }
         
         PyObject *areFamsEqual(const MEDFileMesh *other) const
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

         PyObject *areGrpsEqual(const MEDFileMesh *other) const
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

         PyObject *getAllGeoTypes() const
         {
           std::vector<INTERP_KERNEL::NormalizedCellType> result(self->getAllGeoTypes());
           std::vector<INTERP_KERNEL::NormalizedCellType>::const_iterator iL=result.begin();
           PyObject *res=PyList_New(result.size());
           for(int i=0;iL!=result.end(); i++, iL++)
             PyList_SetItem(res,i,PyInt_FromLong(*iL));
           return res;
         }

         PyObject *getGeoTypesAtLevel(int meshDimRelToMax) const
         {
           std::vector<INTERP_KERNEL::NormalizedCellType> result(self->getGeoTypesAtLevel(meshDimRelToMax));
           std::vector<INTERP_KERNEL::NormalizedCellType>::const_iterator iL=result.begin();
           PyObject *res=PyList_New(result.size());
           for(int i=0;iL!=result.end(); i++, iL++)
             PyList_SetItem(res,i,PyInt_FromLong(*iL));
           return res;
         }

         PyObject *getFamilyFieldAtLevel(int meshDimRelToMaxExt) const
         {
           const DataArrayInt *tmp=self->getFamilyFieldAtLevel(meshDimRelToMaxExt);
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
         }

         PyObject *getOrCreateAndGetFamilyFieldAtLevel(int meshDimRelToMaxExt)
         {
           const DataArrayInt *tmp=self->getOrCreateAndGetFamilyFieldAtLevel(meshDimRelToMaxExt);
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
         }

         PyObject *getNumberFieldAtLevel(int meshDimRelToMaxExt) const
         {
           const DataArrayInt *tmp=self->getNumberFieldAtLevel(meshDimRelToMaxExt);
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
         }

         PyObject *getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const
         {
           const DataArrayInt *tmp=self->getRevNumberFieldAtLevel(meshDimRelToMaxExt);
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
         }
         
         PyObject *getNameFieldAtLevel(int meshDimRelToMaxExt) const
         {
           const DataArrayAsciiChar *tmp=self->getNameFieldAtLevel(meshDimRelToMaxExt);
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_MEDCoupling__DataArrayAsciiChar, SWIG_POINTER_OWN | 0 );
         }

         PyObject *findOrCreateAndGiveFamilyWithId(int id, bool& created)
         {
           bool ret1;
           std::string ret0=self->findOrCreateAndGiveFamilyWithId(id,ret1);
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,PyString_FromString(ret0.c_str()));
           PyTuple_SetItem(ret,1,SWIG_From_bool(ret1));
           return ret;
         }
         
         PyObject *unPolyze()
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
           PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(ret3),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           return ret;
         }

         MEDFileEquivalences *getEquivalences()
         {
           MEDFileEquivalences *ret(self->getEquivalences());
           if(ret) ret->incrRef();
           return ret;
         }

         virtual DataArrayInt *getGlobalNumFieldAtLevel(int meshDimRelToMaxExt) const
         {
           MCAuto<DataArrayInt> ret(self->getGlobalNumFieldAtLevel(meshDimRelToMaxExt));
           return ret.retn();
         }
       }
  };

  class MEDFileUMesh : public MEDFileMesh
  {
  public:
    static MEDFileUMesh *New(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    static MEDFileUMesh *New(const std::string& fileName, MEDFileMeshReadSelector *mrs=0);
    static MEDFileUMesh *New(const MEDCouplingMappedExtrudedMesh *mem);
    static MEDFileUMesh *New(DataArrayByte *db);
    static MEDFileUMesh *New();
    static const char *GetSpeStr4ExtMesh();
    ~MEDFileUMesh();
    int getSpaceDimension() const;
    int getRelativeLevOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const;
    void checkConsistency() const;
    void checkSMESHConsistency() const;
    void clearNodeAndCellNumbers();
    //
    MEDCouplingUMesh *getGroup(int meshDimRelToMaxExt, const std::string& grp, bool renum=false) const;
    MEDCouplingUMesh *getGroups(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum=false) const;
    MEDCouplingUMesh *getFamily(int meshDimRelToMaxExt, const std::string& fam, bool renum=false) const;
    MEDCouplingUMesh *getFamilies(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const;
    DataArrayInt *getNodeGroupsArr(const std::vector<std::string>& grps, bool renum=false) const;
    MEDCouplingUMesh *getLevel0Mesh(bool renum=false) const;
    MEDCouplingUMesh *getLevelM1Mesh(bool renum=false) const;
    MEDCouplingUMesh *getLevelM2Mesh(bool renum=false) const;
    MEDCouplingUMesh *getLevelM3Mesh(bool renum=false) const;
    void forceComputationOfParts() const;
    //
    void setFamilyNameAttachedOnId(int id, const std::string& newFamName);
    void setCoords(DataArrayDouble *coords);
    void setCoordsForced(DataArrayDouble *coords);
    void eraseGroupsAtLevel(int meshDimRelToMaxExt);
    void removeMeshAtLevel(int meshDimRelToMax);
    void setMeshAtLevel(int meshDimRelToMax, MEDCoupling1GTUMesh *m);
    void setMeshAtLevel(int meshDimRelToMax, MEDCouplingUMesh *m, bool newOrOld=false);
    void optimizeFamilies();
    DataArrayInt *zipCoords();
    DataArrayInt *extractFamilyFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const;
    DataArrayInt *extractNumberFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const;
    MEDFileUMesh *buildExtrudedMesh(const MEDCouplingUMesh *m1D, int policy) const;
    MEDFileUMesh *linearToQuadratic(int conversionType=0, double eps=1e-12) const;
    MEDFileUMesh *quadraticToLinear(double eps=1e-12) const;
    MEDCouplingMappedExtrudedMesh *convertToExtrudedMesh() const;
    %extend
       { 
         MEDFileUMesh(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0)
         {
           return MEDFileUMesh::New(fileName,mName,dt,it,mrs);
         }

         MEDFileUMesh(const std::string& fileName, MEDFileMeshReadSelector *mrs=0)
         {
           return MEDFileUMesh::New(fileName,mrs);
         }

         MEDFileUMesh(const MEDCouplingMappedExtrudedMesh *mem)
         {
           return MEDFileUMesh::New(mem);
         }

         MEDFileUMesh(DataArrayByte *db)
         {
           return MEDFileUMesh::New(db);
         }

         MEDFileUMesh()
         {
           return MEDFileUMesh::New();
         }

         static MEDFileUMesh *LoadPartOf(const std::string& fileName, const std::string& mName, PyObject *types, const std::vector<int>& slicPerTyp, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0)
         {
           std::vector<int> typesCpp1;
           convertPyToNewIntArr3(types,typesCpp1);
           std::size_t sz(typesCpp1.size());
           std::vector<INTERP_KERNEL::NormalizedCellType> typesCpp2(sz);
           for(std::size_t ii=0;ii<sz;ii++)
             typesCpp2[ii]=(INTERP_KERNEL::NormalizedCellType)typesCpp1[ii];
           return MEDFileUMesh::LoadPartOf(fileName,mName,typesCpp2,slicPerTyp,dt,it,mrs);
         }

         PyObject *__getstate__()
         {
           std::vector<double> a0;
           std::vector<int> a1;
           std::vector<std::string> a2;
           std::vector< MCAuto<DataArrayInt> > a3;
           MCAuto<DataArrayDouble> a4;
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
               PyList_SetItem(ret3,i,SWIG_NewPointerObj(SWIG_as_voidptr(elt),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
             }
           PyTuple_SetItem(ret,3,ret3);
           DataArrayDouble *ret4(a4);
           if(ret4)
             ret4->incrRef();
           PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(ret4),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
           return ret;
         }

         void __setstate__(PyObject *inp)
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
           std::vector< MCAuto<DataArrayInt> > a3;
           MCAuto<DataArrayDouble> a4;
           //
           PyObject *a0py(PyTuple_GetItem(inp,0)),*a1py(PyTuple_GetItem(inp,1)),*a2py(PyTuple_GetItem(inp,2));
           int tmp(-1);
           fillArrayWithPyListDbl3(a0py,tmp,a0);
           convertPyToNewIntArr3(a1py,a1);
           fillStringVector(a2py,a2);
           //
           PyObject *b0py(PyTuple_GetItem(inp,3)),*b1py(PyTuple_GetItem(inp,4));
           void *argp(0);
           int status(SWIG_ConvertPtr(b1py,&argp,SWIGTYPE_p_MEDCoupling__DataArrayDouble,0|0));
           if(!SWIG_IsOK(status))
             throw INTERP_KERNEL::Exception(MSG);
           a4=reinterpret_cast<DataArrayDouble *>(argp);
           if((DataArrayDouble *)a4)
             a4->incrRef();
           {
             std::vector< DataArrayInt * > a3Tmp;
             convertFromPyObjVectorOfObj<MEDCoupling::DataArrayInt *>(b0py,SWIGTYPE_p_MEDCoupling__DataArrayInt,"DataArrayInt",a3Tmp);
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

         void __setitem__(int meshDimRelToMax, MEDCouplingPointSet *mesh)
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

         void __delitem__(int meshDimRelToMax)
         {
           self->removeMeshAtLevel(meshDimRelToMax);
         }

         MEDFileUMesh *symmetry3DPlane(PyObject *point, PyObject *normalVector) const
         {
           const char msg[]="Python wrap of MEDFileUMesh::symmetry3DPlane : ";
           double val,val2;
           DataArrayDouble *a,*a2;
           DataArrayDoubleTuple *aa,*aa2;
           std::vector<double> bb,bb2;
           int sw;
           const double *centerPtr(convertObjToPossibleCpp5_Safe(point,sw,val,a,aa,bb,msg,1,3,true));
           const double *vectorPtr(convertObjToPossibleCpp5_Safe(normalVector,sw,val2,a2,aa2,bb2,msg,1,3,true));
           MCAuto<MEDFileUMesh> ret(self->symmetry3DPlane(centerPtr,vectorPtr));
           return ret.retn();
         }

         static MEDFileUMesh *Aggregate(PyObject *meshes)
         {
           std::vector<const MEDFileUMesh *> meshesCpp;
           convertFromPyObjVectorOfObj<const MEDCoupling::MEDFileUMesh *>(meshes,SWIGTYPE_p_MEDCoupling__MEDFileUMesh,"MEDFileUMesh",meshesCpp);
           MCAuto<MEDFileUMesh> ret(MEDFileUMesh::Aggregate(meshesCpp));
           return ret.retn();
         }

         PyObject *getAllDistributionOfTypes() const
         {
           std::vector< std::pair<int,int> > ret(self->getAllDistributionOfTypes());
           return convertVecPairIntToPy(ret);
         }
         
         DataArrayInt *deduceNodeSubPartFromCellSubPart(PyObject *extractDef) const
         {
           std::map<int, MCAuto<DataArrayInt> > extractDefCpp;
           convertToMapIntDataArrayInt(extractDef,extractDefCpp);
           return self->deduceNodeSubPartFromCellSubPart(extractDefCpp);
         }

         MEDFileUMesh *extractPart(PyObject *extractDef) const
         {
           std::map<int, MCAuto<DataArrayInt> > extractDefCpp;
           convertToMapIntDataArrayInt(extractDef,extractDefCpp);
           return self->extractPart(extractDefCpp);
         }

         void setMeshes(PyObject *li, bool renum=false)
         {
           std::vector<const MEDCouplingUMesh *> ms;
           convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",ms);
           self->setMeshes(ms,renum);
         }

         void setGroupsFromScratch(int meshDimRelToMax, PyObject *li, bool renum=false)
         {
           std::vector<const MEDCouplingUMesh *> ms;
           convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",ms);
           self->setGroupsFromScratch(meshDimRelToMax,ms,renum);
         }
         
         void setGroupsOnSetMesh(int meshDimRelToMax, PyObject *li, bool renum=false)
         {
           std::vector<const MEDCouplingUMesh *> ms;
           convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingUMesh *>(li,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh",ms);
           self->setGroupsOnSetMesh(meshDimRelToMax,ms,renum);
         }

         DataArrayDouble *getCoords() const
         {
           DataArrayDouble *ret=self->getCoords();
           if(ret)
             ret->incrRef();
           return ret;
         }

         PartDefinition *getPartDefAtLevel(int meshDimRelToMaxExt, INTERP_KERNEL::NormalizedCellType gt=INTERP_KERNEL::NORM_ERROR) const
         {
           const PartDefinition *ret(self->getPartDefAtLevel(meshDimRelToMaxExt,gt));
           if(ret)
             ret->incrRef();
           return const_cast<PartDefinition *>(ret);
         }

         PyObject *buildInnerBoundaryAlongM1Group(const std::string& grpNameM1)
         {
           DataArrayInt *ret0=0,*ret1=0,*ret2=0;
           self->buildInnerBoundaryAlongM1Group(grpNameM1,ret0,ret1,ret2);
           PyObject *ret=PyTuple_New(3);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(ret2),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           return ret;
         }
         
         MEDCoupling1GTUMesh *getDirectUndergroundSingleGeoTypeMesh(INTERP_KERNEL::NormalizedCellType gt) const
         {
           MEDCoupling1GTUMesh *ret(self->getDirectUndergroundSingleGeoTypeMesh(gt));
           if(ret)
             ret->incrRef();
           return ret;
         }

         PyObject *getDirectUndergroundSingleGeoTypeMeshes(int meshDimRelToMax) const
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
      MEDCoupling1SGTUMesh *getImplicitFaceMesh() const
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
    static MEDFileCMesh *New(DataArrayByte *db);
    static MEDFileCMesh *New(const std::string& fileName, MEDFileMeshReadSelector *mrs=0);
    static MEDFileCMesh *New(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    void setMesh(MEDCouplingCMesh *m);
    int getSpaceDimension() const;
    %extend
       {
         MEDFileCMesh()
         {
           return MEDFileCMesh::New();
         }

         MEDFileCMesh(const std::string& fileName, MEDFileMeshReadSelector *mrs=0)
         {
           return MEDFileCMesh::New(fileName,mrs);
         }

         MEDFileCMesh(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0)
         {
           return MEDFileCMesh::New(fileName,mName,dt,it,mrs);
         }

         MEDFileCMesh(DataArrayByte *db)
         {
           return MEDFileCMesh::New(db);
         }
         
         PyObject *getMesh() const
         {
           const MEDCouplingCMesh *tmp=self->getMesh();
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_MEDCoupling__MEDCouplingCMesh, SWIG_POINTER_OWN | 0 );
         }
       }
  };

  class MEDFileCurveLinearMesh : public MEDFileStructuredMesh
  {
  public:
    static MEDFileCurveLinearMesh *New();
    static MEDFileCurveLinearMesh *New(DataArrayByte *db);
    static MEDFileCurveLinearMesh *New(const std::string& fileName, MEDFileMeshReadSelector *mrs=0);
    static MEDFileCurveLinearMesh *New(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    void setMesh(MEDCouplingCurveLinearMesh *m);
    %extend
       {
         MEDFileCurveLinearMesh()
         {
           return MEDFileCurveLinearMesh::New();
         }

         MEDFileCurveLinearMesh(const std::string& fileName, MEDFileMeshReadSelector *mrs=0)
         {
           return MEDFileCurveLinearMesh::New(fileName,mrs);
         }

         MEDFileCurveLinearMesh(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0)
         {
           return MEDFileCurveLinearMesh::New(fileName,mName,dt,it,mrs);
         }

         MEDFileCurveLinearMesh(DataArrayByte *db)
         {
           return MEDFileCurveLinearMesh::New(db);
         }
         
         PyObject *getMesh() const
         {
           const MEDCouplingCurveLinearMesh *tmp=self->getMesh();
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_MEDCoupling__MEDCouplingCurveLinearMesh, SWIG_POINTER_OWN | 0 );
         }
       }
  };

  class MEDFileMeshMultiTS : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    static MEDFileMeshMultiTS *New();
    static MEDFileMeshMultiTS *New(const std::string& fileName);
    static MEDFileMeshMultiTS *New(const std::string& fileName, const std::string& mName);
    MEDFileMeshMultiTS *deepCopy() const;
    std::string getName() const;
    void setOneTimeStep(MEDFileMesh *mesh1TimeStep);
    void cartesianizeMe();
    %extend
       { 
         MEDFileMeshMultiTS()
         {
           return MEDFileMeshMultiTS::New();
         }

         MEDFileMeshMultiTS(const std::string& fileName)
         {
           return MEDFileMeshMultiTS::New(fileName);
         }

         MEDFileMeshMultiTS(const std::string& fileName, const std::string& mName)
         {
           return MEDFileMeshMultiTS::New(fileName,mName);
         }

         MEDFileMesh *getOneTimeStep() const
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
      PyObject *next()
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

  class MEDFileMeshes : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    static MEDFileMeshes *New();
    static MEDFileMeshes *New(DataArrayByte *db);
    MEDFileMeshes *deepCopy() const;
    int getNumberOfMeshes() const;
    std::vector<std::string> getMeshesNames() const;
    //
    void resize(int newSize);
    void pushMesh(MEDFileMesh *mesh);
    void setMeshAtPos(int i, MEDFileMesh *mesh);
    void destroyMeshAtPos(int i);
    void cartesianizeMe();
    bool presenceOfStructureElements() const;
    void killStructureElements();
    %extend
       {
         MEDFileMeshes()
         {
           return MEDFileMeshes::New();
         }

         MEDFileMeshes(const std::string& fileName)
         {
           return MEDFileMeshes::New(fileName);
         }

         MEDFileMeshes(DataArrayByte *db)
         {
           return MEDFileMeshes::New(db);
         }

         std::string __str__() const
           {
             return self->simpleRepr();
           }

         MEDFileMesh *__getitem__(PyObject *obj)
         {
           static const char msg[]="MEDFileMeshes::__getitem__ : only integer or string with meshname supported !";
             if(PyInt_Check(obj))
             {
               MEDFileMesh *ret=self->getMeshAtPos(InterpreteNegativeInt((int)PyInt_AS_LONG(obj),self->getNumberOfMeshes()));
               if(ret)
                 ret->incrRef();
               return ret;
             }
           MEDFileMesh *ret(self->getMeshWithName(convertPyObjectToStr(obj,msg)));
           if(ret)
             ret->incrRef();
           return ret;
         }

         MEDFileMeshes *__setitem__(int obj, MEDFileMesh *mesh)
         {
           self->setMeshAtPos(obj,mesh);
           return self;
         }

         MEDFileMeshesIterator *__iter__()
         {
           return self->iterator();
         }

         int __len__() const
         {
           return self->getNumberOfMeshes();
         }
         
         MEDFileMesh *getMeshAtPos(int i) const
           {
             MEDFileMesh *ret=self->getMeshAtPos(i);
             if(ret)
               ret->incrRef();
             return ret;
           }
         MEDFileMesh *getMeshWithName(const std::string& mname) const
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
    bool isEqual(const MEDFileFieldLoc& other, double eps) const;
  %extend
    {
      std::string __str__() const
      {
        return self->repr();
      }
    }
  };

  class MEDFileFieldGlobsReal
  {
  public:
    void resetContent();
    void shallowCpyGlobs(const MEDFileFieldGlobsReal& other);
    void deepCpyGlobs(const MEDFileFieldGlobsReal& other);
    void shallowCpyOnlyUsedGlobs(const MEDFileFieldGlobsReal& other);
    void deepCpyOnlyUsedGlobs(const MEDFileFieldGlobsReal& other);
    void appendGlobs(const MEDFileFieldGlobsReal& other, double eps);
    void checkGlobsCoherency() const;
    void checkGlobsPflsPartCoherency() const;
    void checkGlobsLocsPartCoherency() const;
    std::vector<std::string> getPfls() const;
    std::vector<std::string> getLocs() const;
    bool existsPfl(const std::string& pflName) const;
    bool existsLoc(const std::string& locName) const;
    std::string createNewNameOfPfl() const;
    std::string createNewNameOfLoc() const;
    std::vector< std::vector<int> > whichAreEqualProfiles() const;
    std::vector< std::vector<int> > whichAreEqualLocs(double eps) const;
    virtual std::vector<std::string> getPflsReallyUsed() const;
    virtual std::vector<std::string> getLocsReallyUsed() const;
    virtual std::vector<std::string> getPflsReallyUsedMulti() const;
    virtual std::vector<std::string> getLocsReallyUsedMulti() const;
    void killProfileIds(const std::vector<int>& pflIds);
    void killLocalizationIds(const std::vector<int>& locIds);
    void changePflName(const std::string& oldName, const std::string& newName);
    void changeLocName(const std::string& oldName, const std::string& newName);
    int getNbOfGaussPtPerCell(int locId) const;
    int getLocalizationId(const std::string& loc) const;
    int getProfileId(const std::string& pfl) const;
    void killStructureElementsInGlobs();
  %extend
     {
       PyObject *getProfile(const std::string& pflName) const
       {
         const DataArrayInt *ret=self->getProfile(pflName);
         if(ret)
           ret->incrRef();
         return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
       }

       PyObject *getProfileFromId(int pflId) const
       {
         const DataArrayInt *ret=self->getProfileFromId(pflId);
         if(ret)
           ret->incrRef();
         return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
       }

       PyObject *getLocalizationFromId(int locId) const
       {
         const MEDFileFieldLoc *loc=&self->getLocalizationFromId(locId);
         if(loc)
           loc->incrRef();
         return SWIG_NewPointerObj(SWIG_as_voidptr(loc),SWIGTYPE_p_MEDCoupling__MEDFileFieldLoc, SWIG_POINTER_OWN | 0 );
       }
       
       PyObject *getLocalization(const std::string& locName) const
       {
         const MEDFileFieldLoc *loc=&self->getLocalization(locName);
         if(loc)
           loc->incrRef();
         return SWIG_NewPointerObj(SWIG_as_voidptr(loc),SWIGTYPE_p_MEDCoupling__MEDFileFieldLoc, SWIG_POINTER_OWN | 0 );
       }
       
       PyObject *zipPflsNames()
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > ret=self->zipPflsNames();
         return convertVecPairVecStToPy(ret);
       }

       PyObject *zipLocsNames(double eps)
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > ret=self->zipLocsNames(eps);
         return convertVecPairVecStToPy(ret);
       }

       void changePflsNames(PyObject *li)
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > v=convertVecPairVecStFromPy(li);
         self->changePflsNames(v);
       }

       void changePflsRefsNamesGen(PyObject *li)
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > v=convertVecPairVecStFromPy(li);
         self->changePflsRefsNamesGen(v);
       }

       void changePflsNamesInStruct(PyObject *li)
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > v=convertVecPairVecStFromPy(li);
         self->changePflsNamesInStruct(v);
       }

       void changeLocsNames(PyObject *li)
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > v=convertVecPairVecStFromPy(li);
         self->changeLocsNames(v);
       }

       void changeLocsRefsNamesGen(PyObject *li)
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > v=convertVecPairVecStFromPy(li);
         self->changeLocsRefsNamesGen(v);
       }
       
       void changeLocsNamesInStruct(PyObject *li)
       {
         std::vector< std::pair<std::vector<std::string>, std::string > > v=convertVecPairVecStFromPy(li);
         self->changeLocsNamesInStruct(v);
       }

       std::string simpleReprGlobs() const
       {
         std::ostringstream oss;
         self->simpleReprGlobs(oss);
         return oss.str();
       }
     }
  };

  class MEDFileEntities
  {
  public:
    %extend
      {
        static MEDFileEntities *BuildFrom(PyObject *entities)
        {
          std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > inp;
          std::vector< std::pair<int,int> > inp0(convertTimePairIdsFromPy(entities));
          {
            std::size_t sz(inp0.size());
            inp.resize(sz);
            for(std::size_t i=0;i<sz;i++)
              inp[i]=std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType>((TypeOfField)inp0[i].first,(INTERP_KERNEL::NormalizedCellType)inp0[i].second);
          }
          return MEDFileEntities::BuildFrom(&inp);
        }
      }
  private:
    MEDFileEntities();
  };

  class MEDFileAnyTypeField1TS : public RefCountObject, public MEDFileFieldGlobsReal, public MEDFileWritableStandAlone
  {
  public:
    static MEDFileAnyTypeField1TS *New(const std::string& fileName, bool loadAll=true);
    static MEDFileAnyTypeField1TS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    static MEDFileAnyTypeField1TS *New(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true);
    static MEDFileAnyTypeField1TS *NewAdv(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileEntities *entities);
    void loadArrays();
    void loadArraysIfNecessary();
    void unloadArrays();
    void unloadArraysWithoutDataLoss();
    int getDimension() const;
    int getIteration() const;
    int getOrder() const;
    std::string getName();
    void setName(const std::string& name);
    std::string getMeshName();
    void setMeshName(const std::string& newMeshName);
    int getMeshIteration() const;
    int getMeshOrder() const;
    int getNumberOfComponents() const;
    bool isDealingTS(int iteration, int order) const;
    void setInfo(const std::vector<std::string>& infos);
    const std::vector<std::string>& getInfo() const;
    bool presenceOfMultiDiscPerGeoType() const;
    void setTime(int iteration, int order, double val);
    virtual MEDFileAnyTypeField1TS *shallowCpy() const;
    MEDFileAnyTypeField1TS *deepCopy() const;
    std::string getDtUnit() const;
    void setDtUnit(const std::string& dtUnit);
    %extend
    {
      PyObject *getTime()
      {
        int tmp1,tmp2;
        double tmp0=self->getTime(tmp1,tmp2);
        PyObject *res = PyList_New(3);
        PyList_SetItem(res,0,SWIG_From_int(tmp1));
        PyList_SetItem(res,1,SWIG_From_int(tmp2));
        PyList_SetItem(res,2,SWIG_From_double(tmp0));
        return res;
      }

      PyObject *getDtIt() const
      {
        std::pair<int,int> res=self->getDtIt();
        PyObject *elt=PyTuple_New(2);
        PyTuple_SetItem(elt,0,SWIG_From_int(res.first));
        PyTuple_SetItem(elt,1,SWIG_From_int(res.second));
        return elt;
      }

      void setProfileNameOnLeaf(INTERP_KERNEL::NormalizedCellType typ, int locId, const std::string& newPflName, bool forceRenameOnGlob=false)
      {
        self->setProfileNameOnLeaf(0,typ,locId,newPflName,forceRenameOnGlob);
      }
      
      void setLocNameOnLeaf(INTERP_KERNEL::NormalizedCellType typ, int locId, const std::string& newLocName, bool forceRenameOnGlob=false)
      {
        self->setLocNameOnLeaf(0,typ,locId,newLocName,forceRenameOnGlob);
      }

      bool changeMeshNames(PyObject *li)
      {
        std::vector< std::pair<std::string,std::string> > modifTab=convertVecPairStStFromPy(li);
        return self->changeMeshNames(modifTab);
      }
      
      PyObject *getTypesOfFieldAvailable() const
      {
        std::vector<TypeOfField> ret=self->getTypesOfFieldAvailable();
        PyObject *ret2=PyList_New(ret.size());
        for(int i=0;i<(int)ret.size();i++)
          PyList_SetItem(ret2,i,SWIG_From_int(ret[i]));
        return ret2;
      }

      PyObject *getNonEmptyLevels(const std::string& mname=std::string()) const
      {
        std::vector<int> ret1;
        int ret0=self->getNonEmptyLevels(mname,ret1);
        PyObject *elt=PyTuple_New(2);
        PyTuple_SetItem(elt,0,SWIG_From_int(ret0));
        PyTuple_SetItem(elt,1,convertIntArrToPyList2(ret1));
        return elt;
      }

      PyObject *getFieldSplitedByType(const std::string& mname=std::string()) const
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

      PyObject *splitComponents() const
      {
        std::vector< MCAuto< MEDFileAnyTypeField1TS > > ret=self->splitComponents();
        std::size_t sz=ret.size();
        PyObject *retPy=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(retPy,i,convertMEDFileField1TS(ret[i].retn(), SWIG_POINTER_OWN | 0 ));
        return retPy;
      }

      PyObject *splitDiscretizations() const
      {
        std::vector< MCAuto< MEDFileAnyTypeField1TS > > ret=self->splitDiscretizations();
        std::size_t sz=ret.size();
        PyObject *retPy=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(retPy,i,convertMEDFileField1TS(ret[i].retn(), SWIG_POINTER_OWN | 0 ));
        return retPy;
      }

      PyObject *splitMultiDiscrPerGeoTypes() const
      {
        std::vector< MCAuto< MEDFileAnyTypeField1TS > > ret=self->splitMultiDiscrPerGeoTypes();
        std::size_t sz=ret.size();
        PyObject *retPy=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(retPy,i,convertMEDFileField1TS(ret[i].retn(), SWIG_POINTER_OWN | 0 ));
        return retPy;
      }

      MEDFileAnyTypeField1TS *extractPart(PyObject *extractDef, MEDFileMesh *mm) const
      {
        std::map<int, MCAuto<DataArrayInt> > extractDefCpp;
        convertToMapIntDataArrayInt(extractDef,extractDefCpp);
        return self->extractPart(extractDefCpp,mm);
      }
    }
  };

  class MEDFileField1TS : public MEDFileAnyTypeField1TS
  {
  public:
    static MEDFileField1TS *New(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true);
    static MEDFileField1TS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    static MEDFileField1TS *New(const std::string& fileName, bool loadAll=true);
    static MEDFileField1TS *New(DataArrayByte *db);
    static MEDFileField1TS *New();
    MEDCoupling::MEDFileIntField1TS *convertToInt(bool isDeepCpyGlobs=true) const;
    void copyTimeInfoFrom(MEDCouplingFieldDouble *mcf);
    MEDCouplingFieldDouble *field(const MEDFileMesh *mesh) const;
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol=0) const;
    MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, int renumPol=0) const;
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, int renumPol=0) const;
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const;
    MEDCouplingFieldDouble *getFieldAtLevelOld(TypeOfField type, const std::string& mname, int meshDimRelToMax, int renumPol=0) const;
    //
    void setFieldNoProfileSBT(const MEDCouplingFieldDouble *field);
    void setFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    void setFieldProfileFlatly(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    void setProfileNameOnLeaf(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const std::string& newPflName, bool forceRenameOnGlob=false);
    void setLocNameOnLeaf(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const std::string& newLocName, bool forceRenameOnGlob=false);
    %extend
       {
         MEDFileField1TS(const std::string& fileName, bool loadAll=true)
         {
           return MEDFileField1TS::New(fileName,loadAll);
         }
         
         MEDFileField1TS(const std::string& fileName, const std::string& fieldName, bool loadAll=true)
         {
           return MEDFileField1TS::New(fileName,fieldName,loadAll);
         }

         MEDFileField1TS(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true)
         {
           return MEDFileField1TS::New(fileName,fieldName,iteration,order,loadAll);
         }

         MEDFileField1TS(DataArrayByte *db)
         {
           return MEDFileField1TS::New(db);
         }

         MEDFileField1TS()
         {
           return MEDFileField1TS::New();
         }
         
         void copyTinyInfoFrom(const MEDCouplingFieldDouble *field)
         {
           const DataArrayDouble *arr=0;
           if(field)
             arr=field->getArray();
           self->copyTinyInfoFrom(field,arr);
         }
         
         std::string __str__() const
         {
           return self->simpleRepr();
         }
         
         PyObject *getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh) const
         {
           return MEDFileField1TS_getFieldWithProfile<double>(self,type,meshDimRelToMax,mesh);
         }

         PyObject *getFieldSplitedByType2(const std::string& mname=std::string()) const
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
                   PyTuple_SetItem(elt3,1,SWIG_NewPointerObj(SWIG_as_voidptr(dadsI[j]),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
                   PyTuple_SetItem(elt3,2,PyString_FromString(pflsI[j].c_str()));
                   PyTuple_SetItem(elt3,3,PyString_FromString(locsI[j].c_str()));
                   PyList_SetItem(elt2,j,elt3);
                 }
               PyTuple_SetItem(elt,1,elt2);
               PyList_SetItem(ret2,i,elt);
             }
           return ret2;
         }

         DataArrayDouble *getUndergroundDataArray() const
         {
           DataArrayDouble *ret=self->getUndergroundDataArray();
           if(ret)
             ret->incrRef();
           return ret;
         }

         PyObject *getUndergroundDataArrayExt() const
         {
           return MEDFileField1TS_getUndergroundDataArrayExt<double>(self);
         }
       }
  };

  class MEDFileIntField1TS : public MEDFileAnyTypeField1TS
  {
  public:
    static MEDFileIntField1TS *New();
    static MEDFileIntField1TS *New(const std::string& fileName, bool loadAll=true);
    static MEDFileIntField1TS *New(DataArrayByte *db);
    static MEDFileIntField1TS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    static MEDFileIntField1TS *New(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true);
    MEDCoupling::MEDFileField1TS *convertToDouble(bool isDeepCpyGlobs=true) const;
    //
    void setFieldNoProfileSBT(const MEDCouplingFieldInt *field);
    void setFieldProfile(const MEDCouplingFieldInt *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    void setFieldProfileFlatly(const MEDCouplingFieldInt *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    void copyTimeInfoFrom(MEDCouplingFieldInt *mcf);
    MEDCouplingFieldInt *field(const MEDFileMesh *mesh) const;
    MEDCouplingFieldInt *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol=0) const;
    MEDCouplingFieldInt *getFieldAtTopLevel(TypeOfField type, int renumPol=0) const;
    MEDCouplingFieldInt *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const;
    MEDCouplingFieldInt *getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, int renumPol=0) const;
    MEDCouplingFieldInt *getFieldAtLevelOld(TypeOfField type, const std::string& mname, int meshDimRelToMax, int renumPol=0) const;
    %extend
    {
      MEDFileIntField1TS()
      {
        return MEDFileIntField1TS::New();
      }

      MEDFileIntField1TS(const std::string& fileName, bool loadAll=true)
      {
        return MEDFileIntField1TS::New(fileName,loadAll);
      }

      MEDFileIntField1TS(const std::string& fileName, const std::string& fieldName, bool loadAll=true)
      {
        return MEDFileIntField1TS::New(fileName,fieldName,loadAll);
      }

      MEDFileIntField1TS(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true)
      {
        return MEDFileIntField1TS::New(fileName,fieldName,iteration,order,loadAll);
      }

      MEDFileIntField1TS(DataArrayByte *db)
      {
        return MEDFileIntField1TS::New(db);
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }

      PyObject *getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh) const
      {
         return MEDFileField1TS_getFieldWithProfile<int>(self,type,meshDimRelToMax,mesh);
      }
      
      DataArrayInt *getUndergroundDataArray() const
      {
        DataArrayInt *ret=self->getUndergroundDataArray();
        if(ret)
          ret->incrRef();
        return ret;
      }

      PyObject *getUndergroundDataArrayExt() const
      {
        return MEDFileField1TS_getUndergroundDataArrayExt<int>(self);
      }
    }
  };

  class MEDFileFloatField1TS : public MEDFileAnyTypeField1TS
  {
  public:
    static MEDFileFloatField1TS *New();
    static MEDFileFloatField1TS *New(const std::string& fileName, bool loadAll=true);
    static MEDFileFloatField1TS *New(DataArrayByte *db);
    static MEDFileFloatField1TS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    static MEDFileFloatField1TS *New(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true);
    MEDCoupling::MEDFileField1TS *convertToDouble(bool isDeepCpyGlobs=true) const;
    //
    void setFieldNoProfileSBT(const MEDCouplingFieldFloat *field);
    void setFieldProfile(const MEDCouplingFieldFloat *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    void setFieldProfileFlatly(const MEDCouplingFieldFloat *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    void copyTimeInfoFrom(MEDCouplingFieldFloat *mcf);
    MEDCouplingFieldFloat *field(const MEDFileMesh *mesh) const;
    MEDCouplingFieldFloat *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol=0) const;
    MEDCouplingFieldFloat *getFieldAtTopLevel(TypeOfField type, int renumPol=0) const;
    MEDCouplingFieldFloat *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const;
    MEDCouplingFieldFloat *getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, int renumPol=0) const;
    MEDCouplingFieldFloat *getFieldAtLevelOld(TypeOfField type, const std::string& mname, int meshDimRelToMax, int renumPol=0) const;
    %extend
    {
      MEDFileFloatField1TS()
      {
        return MEDFileFloatField1TS::New();
      }

      MEDFileFloatField1TS(const std::string& fileName, bool loadAll=true)
      {
        return MEDFileFloatField1TS::New(fileName,loadAll);
      }

      MEDFileFloatField1TS(const std::string& fileName, const std::string& fieldName, bool loadAll=true)
      {
        return MEDFileFloatField1TS::New(fileName,fieldName,loadAll);
      }

      MEDFileFloatField1TS(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true)
      {
        return MEDFileFloatField1TS::New(fileName,fieldName,iteration,order,loadAll);
      }

      MEDFileFloatField1TS(DataArrayByte *db)
      {
        return MEDFileFloatField1TS::New(db);
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }

      PyObject *getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh) const
      {
         return MEDFileField1TS_getFieldWithProfile<float>(self,type,meshDimRelToMax,mesh);
      }
      
      DataArrayFloat *getUndergroundDataArray() const
      {
        DataArrayFloat *ret=self->getUndergroundDataArray();
        if(ret)
          ret->incrRef();
        return ret;
      }
      
      PyObject *getUndergroundDataArrayExt() const
      {
        return MEDFileField1TS_getUndergroundDataArrayExt<float>(self);
      }
    }
  };

  class MEDFileAnyTypeFieldMultiTSIterator
  {
  public:
    %extend
    {
      PyObject *next()
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

  class MEDFileAnyTypeFieldMultiTS : public RefCountObject, public MEDFileFieldGlobsReal, public MEDFileWritableStandAlone
  {
  public:
    static MEDFileAnyTypeFieldMultiTS *New(const std::string& fileName, bool loadAll=true);
    static MEDFileAnyTypeFieldMultiTS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    MEDFileAnyTypeFieldMultiTS *deepCopy() const;
    virtual MEDFileAnyTypeFieldMultiTS *shallowCpy() const;
    std::string getName() const;
    void setName(const std::string& name);
    std::string getDtUnit() const;
    void setDtUnit(const std::string& dtUnit);
    std::string getMeshName() const;
    void setMeshName(const std::string& newMeshName);
    const std::vector<std::string>& getInfo() const;
    bool presenceOfMultiDiscPerGeoType() const;
    int getNumberOfComponents() const;
    int getNumberOfTS() const;
    void eraseEmptyTS();
    int getPosOfTimeStep(int iteration, int order) const;
    int getPosGivenTime(double time, double eps=1e-8) const;
    void loadArrays();
    void loadArraysIfNecessary();
    void unloadArrays();
    void unloadArraysWithoutDataLoss();
    //
    virtual MEDFileAnyTypeField1TS *getTimeStepAtPos(int pos) const;
    MEDFileAnyTypeField1TS *getTimeStep(int iteration, int order) const;
    MEDFileAnyTypeField1TS *getTimeStepGivenTime(double time, double eps=1e-8) const;
    void pushBackTimeStep(MEDFileAnyTypeField1TS *f1ts);
    void synchronizeNameScope();
    MEDFileAnyTypeFieldMultiTS *buildNewEmpty() const;
    %extend
    {
      int __len__() const
      {
        return self->getNumberOfTS();
      }

      int getTimeId(PyObject *elt0) const
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
      
      PyObject *getIterations() const
      {
        std::vector< std::pair<int,int> > res(self->getIterations());
        return convertVecPairIntToPy(res);
      }
      
      PyObject *getTimeSteps() const
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
      
      PyObject *getTypesOfFieldAvailable() const
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
      
      PyObject *getNonEmptyLevels(int iteration, int order, const std::string& mname=std::string()) const
      {
        std::vector<int> ret1;
        int ret0=self->getNonEmptyLevels(iteration,order,mname,ret1);
        PyObject *elt=PyTuple_New(2);
        PyTuple_SetItem(elt,0,SWIG_From_int(ret0));
        PyTuple_SetItem(elt,1,convertIntArrToPyList2(ret1));
        return elt;
      }
      
      PyObject *getFieldSplitedByType(int iteration, int order, const std::string& mname=std::string()) const
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

      std::vector<int> getTimeIds(PyObject *elts) const
      {
        if(PyList_Check(elts))
          {
            int sz=PyList_Size(elts);
            std::vector<int> ret(sz);
            for(int i=0;i<sz;i++)
              {
                PyObject *elt=PyList_GetItem(elts,i);
                ret[i]=MEDCoupling_MEDFileAnyTypeFieldMultiTS_getTimeId(self,elt);
              }
            return ret;
          }
        else
          {
            std::vector<int> ret(1);
            ret[0]=MEDCoupling_MEDFileAnyTypeFieldMultiTS_getTimeId(self,elts);
            return ret;
          }
      }
      
      void __delitem__(PyObject *elts)
      {
        if(PySlice_Check(elts))
          {
            Py_ssize_t strt=2,stp=2,step=2;
            GetIndicesOfSlice(elts,self->getNumberOfTS(),&strt,&stp,&step,"MEDFileAnyTypeFieldMultiTS.__delitem__ : error in input slice !");
            self->eraseTimeStepIds2(strt,stp,step);
          }
        else
          {
            std::vector<int> idsToRemove=MEDCoupling_MEDFileAnyTypeFieldMultiTS_getTimeIds(self,elts);
            if(!idsToRemove.empty())
              self->eraseTimeStepIds(&idsToRemove[0],&idsToRemove[0]+idsToRemove.size());
          }
      }
      
      void eraseTimeStepIds(PyObject *li)
      {
        int sw;
        int pos1;
        std::vector<int> pos2;
        DataArrayInt *pos3=0;
        DataArrayIntTuple *pos4=0;
        convertIntStarLikePyObjToCpp(li,sw,pos1,pos2,pos3,pos4);
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

      MEDFileAnyTypeFieldMultiTSIterator *__iter__()
      {
        return self->iterator();
      }

      PyObject *__getitem__(PyObject *elt0) const
      {
        if(elt0 && PyList_Check(elt0))
          {
            int sz=PyList_Size(elt0);
            MCAuto<DataArrayInt> da=DataArrayInt::New(); da->alloc(sz,1);
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
            GetIndicesOfSlice(elt0,self->getNumberOfTS(),&strt,&stp,&step,"MEDFileAnyTypeFieldMultiTS.__getitem__ : error in input slice !");
            return convertMEDFileFieldMultiTS(self->buildSubPartSlice(strt,stp,step),SWIG_POINTER_OWN | 0);
          }
        else
          return convertMEDFileField1TS(self->getTimeStepAtPos(MEDFileAnyTypeFieldMultiTSgetitemSingleTS__(self,elt0)),SWIG_POINTER_OWN | 0);
      }

      bool changeMeshNames(PyObject *li)
      {
        std::vector< std::pair<std::string,std::string> > modifTab=convertVecPairStStFromPy(li);
        return self->changeMeshNames(modifTab);
      }

      PyObject *splitComponents() const
      {
        std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > ret=self->splitComponents();
        std::size_t sz=ret.size();
        PyObject *retPy=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(retPy,i,convertMEDFileFieldMultiTS(ret[i].retn(), SWIG_POINTER_OWN | 0 ));
        return retPy;
      }

      PyObject *splitDiscretizations() const
      {
        std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > ret=self->splitDiscretizations();
        std::size_t sz=ret.size();
        PyObject *retPy=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(retPy,i,convertMEDFileFieldMultiTS(ret[i].retn(), SWIG_POINTER_OWN | 0 ));
        return retPy;
      }

      PyObject *splitMultiDiscrPerGeoTypes() const
      {
        std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > ret=self->splitMultiDiscrPerGeoTypes();
        std::size_t sz=ret.size();
        PyObject *retPy=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(retPy,i,convertMEDFileFieldMultiTS(ret[i].retn(), SWIG_POINTER_OWN | 0 ));
        return retPy;
      }

      void pushBackTimeSteps(PyObject *li)
      {
        void *argp(0);
        int status(SWIG_ConvertPtr(li,&argp,SWIGTYPE_p_MEDCoupling__MEDFileAnyTypeFieldMultiTS,0|0));
        if(SWIG_IsOK(status))
          {
            self->pushBackTimeSteps(reinterpret_cast<MEDFileAnyTypeFieldMultiTS *>(argp));
          }
        else
          {
            std::vector<MEDFileAnyTypeField1TS *> tmp;
            convertFromPyObjVectorOfObj<MEDCoupling::MEDFileAnyTypeField1TS *>(li,SWIGTYPE_p_MEDCoupling__MEDFileAnyTypeField1TS,"MEDFileAnyTypeField1TS",tmp);
            self->pushBackTimeSteps(tmp);
          }
      }

      MEDFileAnyTypeFieldMultiTS *extractPart(PyObject *extractDef, MEDFileMesh *mm) const
      {
        std::map<int, MCAuto<DataArrayInt> > extractDefCpp;
        convertToMapIntDataArrayInt(extractDef,extractDefCpp);
        return self->extractPart(extractDefCpp,mm);
      }

      static PyObject *MEDFileAnyTypeFieldMultiTS::SplitIntoCommonTimeSeries(PyObject *li)
      {
        std::vector<MEDFileAnyTypeFieldMultiTS *> vectFMTS;
        convertFromPyObjVectorOfObj<MEDCoupling::MEDFileAnyTypeFieldMultiTS *>(li,SWIGTYPE_p_MEDCoupling__MEDFileAnyTypeFieldMultiTS,"MEDFileAnyTypeFieldMultiTS",vectFMTS);
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
      
      static PyObject *MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupport(PyObject *li, const MEDFileMesh *mesh)
      {
        std::vector<MEDFileAnyTypeFieldMultiTS *> vectFMTS;
        convertFromPyObjVectorOfObj<MEDCoupling::MEDFileAnyTypeFieldMultiTS *>(li,SWIGTYPE_p_MEDCoupling__MEDFileAnyTypeFieldMultiTS,"MEDFileAnyTypeFieldMultiTS",vectFMTS);
        std::vector< MCAuto<MEDFileFastCellSupportComparator> > ret2;
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
            PyTuple_SetItem(ret0Py,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret2[i].retn()),SWIGTYPE_p_MEDCoupling__MEDFileFastCellSupportComparator, SWIG_POINTER_OWN | 0 ));
            PyList_SetItem(retPy,i,ret0Py);
          }
        return retPy;
      }
    }
  };

  class MEDFileIntFieldMultiTS;
  
  class MEDFileFieldMultiTS : public MEDFileAnyTypeFieldMultiTS
  {
  public:
    static MEDFileFieldMultiTS *New();
    static MEDFileFieldMultiTS *New(const std::string& fileName, bool loadAll=true);
    static MEDFileFieldMultiTS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    static MEDFileFieldMultiTS *New(DataArrayByte *db);
    //
    MEDCouplingFieldDouble *field(int iteration, int order, const MEDFileMesh *mesh) const;
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol=0) const;
    MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, int iteration, int order, int renumPol=0) const;
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const;
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, int renumPol=0) const;
    MEDCouplingFieldDouble *getFieldAtLevelOld(TypeOfField type, int iteration, int order, const std::string& mname, int meshDimRelToMax, int renumPol=0) const;
    //
    void appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field);
    void appendFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    MEDFileIntFieldMultiTS *convertToInt(bool isDeepCpyGlobs=true) const;
    %extend
       {
         MEDFileFieldMultiTS()
         {
           return MEDFileFieldMultiTS::New();
         }

         MEDFileFieldMultiTS(const std::string& fileName, bool loadAll=true)
         {
           return MEDFileFieldMultiTS::New(fileName,loadAll);
         }

         MEDFileFieldMultiTS(const std::string& fileName, const std::string& fieldName, bool loadAll=true)
         {
           return MEDFileFieldMultiTS::New(fileName,fieldName,loadAll);
         }
         
         MEDFileFieldMultiTS(DataArrayByte *db)
         {
           return MEDFileFieldMultiTS::New(db);
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
         
         std::string __str__() const
         {
           return self->simpleRepr();
         }

         PyObject *getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh) const
         {
           DataArrayInt *ret1=0;
           DataArrayDouble *ret0=self->getFieldWithProfile(type,iteration,order,meshDimRelToMax,mesh,ret1);
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           return ret;
         }

         PyObject *getFieldSplitedByType2(int iteration, int order, const std::string& mname=std::string()) const
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
                   PyTuple_SetItem(elt3,1,SWIG_NewPointerObj(SWIG_as_voidptr(dadsI[j]),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
                   PyTuple_SetItem(elt3,2,PyString_FromString(pflsI[j].c_str()));
                   PyTuple_SetItem(elt3,3,PyString_FromString(locsI[j].c_str()));
                   PyList_SetItem(elt2,j,elt3);
                 }
               PyTuple_SetItem(elt,1,elt2);
               PyList_SetItem(ret2,i,elt);
             }
           return ret2;
         }
         DataArrayDouble *getUndergroundDataArray(int iteration, int order) const
         {
           DataArrayDouble *ret=self->getUndergroundDataArray(iteration,order);
           if(ret)
             ret->incrRef();
           return ret;
         }
         
         PyObject *getUndergroundDataArrayExt(int iteration, int order) const
         {
           std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > > elt1Cpp;
           DataArrayDouble *elt0=self->getUndergroundDataArrayExt(iteration,order,elt1Cpp);
           if(elt0)
             elt0->incrRef();
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(elt0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
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
      PyObject *next()
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
    static MEDFileIntFieldMultiTS *New(const std::string& fileName, bool loadAll=true);
    static MEDFileIntFieldMultiTS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    static MEDFileIntFieldMultiTS *New(DataArrayByte *db);
    //
    void appendFieldNoProfileSBT(const MEDCouplingFieldInt *field);
    void appendFieldProfile(const MEDCouplingFieldInt *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    MEDCoupling::MEDFileFieldMultiTS *convertToDouble(bool isDeepCpyGlobs=true) const;
    MEDCouplingFieldInt *field(int iteration, int order, const MEDFileMesh *mesh) const;
    MEDCouplingFieldInt *getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol=0) const;
    MEDCouplingFieldInt *getFieldAtTopLevel(TypeOfField type, int iteration, int order, int renumPol=0) const;
    MEDCouplingFieldInt *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const;
    MEDCouplingFieldInt *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, int renumPol=0) const;
    MEDCouplingFieldInt *getFieldAtLevelOld(TypeOfField type, int iteration, int order, const std::string& mname, int meshDimRelToMax, int renumPol=0) const;
    %extend
    {
      MEDFileIntFieldMultiTS()
      {
        return MEDFileIntFieldMultiTS::New();
      }
      
      MEDFileIntFieldMultiTS(const std::string& fileName, bool loadAll=true)
      {
        return MEDFileIntFieldMultiTS::New(fileName,loadAll);
      }
      
      MEDFileIntFieldMultiTS(const std::string& fileName, const std::string& fieldName, bool loadAll=true)
      {
        return MEDFileIntFieldMultiTS::New(fileName,fieldName,loadAll);
      }

      MEDFileIntFieldMultiTS(DataArrayByte *db)
      {
        return MEDFileIntFieldMultiTS::New(db);
      }
      
      static MEDFileIntFieldMultiTS *LoadSpecificEntities(const std::string& fileName, const std::string& fieldName, PyObject *entities, bool loadAll=true)
      {
        std::vector<std::pair<int,int> > tmp(convertTimePairIdsFromPy(entities));
        std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > entitiesCpp(convertVecPairIntToVecPairTOFCT(tmp));
        return MEDFileIntFieldMultiTS::LoadSpecificEntities(fileName,fieldName,entitiesCpp,loadAll);
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }

      PyObject *getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh) const
      {
         DataArrayInt *ret1=0;
         DataArrayInt *ret0=self->getFieldWithProfile(type,iteration,order,meshDimRelToMax,mesh,ret1);
         PyObject *ret=PyTuple_New(2);
         PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
         PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
         return ret;
      }

      DataArrayInt *getUndergroundDataArray(int iteration, int order) const
      {
        DataArrayInt *ret=self->getUndergroundDataArray(iteration,order);
        if(ret)
          ret->incrRef();
        return ret;
      }
    }
  };

  class MEDFileFloatFieldMultiTS : public MEDFileAnyTypeFieldMultiTS
  {
  public:
    static MEDFileFloatFieldMultiTS *New();
    static MEDFileFloatFieldMultiTS *New(const std::string& fileName, bool loadAll=true);
    static MEDFileFloatFieldMultiTS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    static MEDFileFloatFieldMultiTS *New(DataArrayByte *db);
    //
    void appendFieldNoProfileSBT(const MEDCouplingFieldFloat *field);
    void appendFieldProfile(const MEDCouplingFieldFloat *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    MEDCoupling::MEDFileFieldMultiTS *convertToDouble(bool isDeepCpyGlobs=true) const;
    MEDCouplingFieldFloat *field(int iteration, int order, const MEDFileMesh *mesh) const;
    MEDCouplingFieldFloat *getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol=0) const;
    MEDCouplingFieldFloat *getFieldAtTopLevel(TypeOfField type, int iteration, int order, int renumPol=0) const;
    MEDCouplingFieldFloat *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const;
    MEDCouplingFieldFloat *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, int renumPol=0) const;
    MEDCouplingFieldFloat *getFieldAtLevelOld(TypeOfField type, int iteration, int order, const std::string& mname, int meshDimRelToMax, int renumPol=0) const;
    %extend
    {
      MEDFileFloatFieldMultiTS()
      {
        return MEDFileFloatFieldMultiTS::New();
      }
      
      MEDFileFloatFieldMultiTS(const std::string& fileName, bool loadAll=true)
      {
        return MEDFileFloatFieldMultiTS::New(fileName,loadAll);
      }
      
      MEDFileFloatFieldMultiTS(const std::string& fileName, const std::string& fieldName, bool loadAll=true)
      {
        return MEDFileFloatFieldMultiTS::New(fileName,fieldName,loadAll);
      }

      MEDFileFloatFieldMultiTS(DataArrayByte *db)
      {
        return MEDFileFloatFieldMultiTS::New(db);
      }
      
      static MEDFileFloatFieldMultiTS *LoadSpecificEntities(const std::string& fileName, const std::string& fieldName, PyObject *entities, bool loadAll=true)
      {
        std::vector<std::pair<int,int> > tmp(convertTimePairIdsFromPy(entities));
        std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > entitiesCpp(convertVecPairIntToVecPairTOFCT(tmp));
        return MEDFileFloatFieldMultiTS::LoadSpecificEntities(fileName,fieldName,entitiesCpp,loadAll);
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }

      PyObject *getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh) const
      {
         DataArrayInt *ret1=0;
         DataArrayFloat *ret0=self->getFieldWithProfile(type,iteration,order,meshDimRelToMax,mesh,ret1);
         PyObject *ret=PyTuple_New(2);
         PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayFloat, SWIG_POINTER_OWN | 0 ));
         PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
         return ret;
      }

      DataArrayFloat *getUndergroundDataArray(int iteration, int order) const
      {
        DataArrayFloat *ret=self->getUndergroundDataArray(iteration,order);
        if(ret)
          ret->incrRef();
        return ret;
      }
    }
  };
  
  class MEDFileMeshSupports : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    static MEDFileMeshSupports *New(const std::string& fileName);
    std::vector<std::string> getSupMeshNames() const;
    %extend
       {
         MEDFileUMesh *getSupMeshWithName(const std::string& name) const
         {
           const MEDFileUMesh *ret(self->getSupMeshWithName(name));
           MEDFileUMesh *ret2(const_cast<MEDFileUMesh *>(ret));
           if(ret2)
             ret2->incrRef();
           return ret2;
         }
       }
  };
 
  class MEDFileStructureElements : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    static MEDFileStructureElements *New(const std::string& fileName, const MEDFileMeshSupports *ms);
  private:
    MEDFileStructureElements();
  };

  class MEDFileFields : public RefCountObject, public MEDFileFieldGlobsReal, public MEDFileWritableStandAlone
  {
  public:
    static MEDFileFields *New();
    static MEDFileFields *New(const std::string& fileName, bool loadAll=true);
    static MEDFileFields *New(DataArrayByte *db);
    static MEDFileFields *NewAdv(const std::string& fileName, bool loadAll, const MEDFileEntities *entities);
    static MEDFileFields *LoadPartOf(const std::string& fileName, bool loadAll=true, const MEDFileMeshes *ms=0);
    static MEDFileFields *NewWithDynGT(const std::string& fileName, const MEDFileStructureElements *se, bool loadAll=true);
    MEDFileFields *deepCopy() const;
    MEDFileFields *shallowCpy() const;
    void loadArrays();
    void loadArraysIfNecessary();
    void unloadArrays();
    void unloadArraysWithoutDataLoss();
    int getNumberOfFields() const;
    std::vector<std::string> getFieldsNames() const;
    std::vector<std::string> getMeshesNames() const;
    //
    void resize(int newSize);
    void pushField(MEDFileAnyTypeFieldMultiTS *field);
    void setFieldAtPos(int i, MEDFileAnyTypeFieldMultiTS *field);
    int getPosFromFieldName(const std::string& fieldName) const;
    MEDFileAnyTypeFieldMultiTS *getFieldAtPos(int i) const;
    MEDFileAnyTypeFieldMultiTS *getFieldWithName(const std::string& fieldName) const;
    MEDFileFields *partOfThisLyingOnSpecifiedMeshName(const std::string& meshName) const;
    bool presenceOfStructureElements() const;
    void aggregate(const MEDFileFields& other);
    void killStructureElements();
    void keepOnlyStructureElements();
    void keepOnlyOnMeshSE(const std::string& meshName, const std::string& seName);
    void blowUpSE(MEDFileMeshes *ms, const MEDFileStructureElements *ses);
    void destroyFieldAtPos(int i);
    bool removeFieldsWithoutAnyTimeStep();
    %extend
       {
         MEDFileFields()
         {
           return MEDFileFields::New();
         }

         MEDFileFields(const std::string& fileName, bool loadAll=true)
         {
           return MEDFileFields::New(fileName,loadAll);
         }

         MEDFileFields(DataArrayByte *db)
         {
           return MEDFileFields::New(db);
         }

         MEDFileFields(const std::string& fileName, bool loadAll, const MEDFileEntities *entities)
         {
           return MEDFileFields::NewAdv(fileName,loadAll,entities);
         }
         
         std::string __str__() const
         {
           return self->simpleRepr();
         }
         
         MEDFileFields *partOfThisOnStructureElements() const
         {
           MCAuto<MEDFileFields> ret(self->partOfThisOnStructureElements());
           return ret.retn();
         }

         MEDFileFields *partOfThisLyingOnSpecifiedMeshSEName(const std::string& meshName, const std::string& seName) const
         {
           MCAuto<MEDFileFields> ret(self->partOfThisLyingOnSpecifiedMeshSEName(meshName,seName));
           return ret.retn();
         }
         
         static MEDFileFields *LoadSpecificEntities(const std::string& fileName, PyObject *entities, bool loadAll=true)
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

         PyObject *getMeshSENames() const
         {
           std::vector< std::pair<std::string,std::string> > ps;
           self->getMeshSENames(ps);
           return convertVectPairStToPy(ps);
         }

         PyObject *getCommonIterations() const
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

         MEDFileFields *partOfThisLyingOnSpecifiedTimeSteps(PyObject *timeSteps) const
         {
           std::vector< std::pair<int,int> > ts=convertTimePairIdsFromPy(timeSteps);
           return self->partOfThisLyingOnSpecifiedTimeSteps(ts);
         }

         MEDFileFields *partOfThisNotLyingOnSpecifiedTimeSteps(PyObject *timeSteps) const
         {
           std::vector< std::pair<int,int> > ts=convertTimePairIdsFromPy(timeSteps);
           return self->partOfThisNotLyingOnSpecifiedTimeSteps(ts);
         }
         
         PyObject *__getitem__(PyObject *obj)
         {
           if(obj && PyList_Check(obj))
             {
               int sz=PyList_Size(obj);
               MCAuto<DataArrayInt> da=DataArrayInt::New(); da->alloc(sz,1);
               int *pt=da->getPointer();
               for(int i=0;i<sz;i++,pt++)
                 {
                   PyObject *elt1=PyList_GetItem(obj,i);
                   *pt=MEDFileFieldsgetitemSingleTS__(self,elt1);
                 }
               return SWIG_NewPointerObj(SWIG_as_voidptr(self->buildSubPart(da->begin(),da->end())),SWIGTYPE_p_MEDCoupling__MEDFileFields, SWIG_POINTER_OWN | 0 );
             }
           else
             return convertMEDFileFieldMultiTS(self->getFieldAtPos(MEDFileFieldsgetitemSingleTS__(self,obj)), SWIG_POINTER_OWN | 0 );
         }

         MEDFileFields *__setitem__(int obj, MEDFileFieldMultiTS *field)
         {
           self->setFieldAtPos(obj,field);
           return self;
         }

         int __len__() const
         {
           return self->getNumberOfFields();
         }

         MEDFileFieldsIterator *__iter__()
         {
           return self->iterator();
         }
         
         bool changeMeshNames(PyObject *li)
         {
           std::vector< std::pair<std::string,std::string> > modifTab=convertVecPairStStFromPy(li);
           return self->changeMeshNames(modifTab);
         }

         int getPosOfField(PyObject *elt0) const
         {
           static const char msg[]="MEDFileFields::getPosOfField : invalid input params ! expected fields[int], fields[string_of_field_name] !";
           if(!elt0)
             throw INTERP_KERNEL::Exception(msg);
           if(PyInt_Check(elt0))
             {//fmts[3]
               return PyInt_AS_LONG(elt0);
             }
           return self->getPosFromFieldName(convertPyObjectToStr(elt0,msg));
         }
         
         std::vector<int> getPosOfFields(PyObject *elts) const
         {
           if(PyList_Check(elts))
             {
               int sz=PyList_Size(elts);
               std::vector<int> ret(sz);
               for(int i=0;i<sz;i++)
                 {
                   PyObject *elt=PyList_GetItem(elts,i);
                   ret[i]=MEDCoupling_MEDFileFields_getPosOfField(self,elt);
                 }
               return ret;
             }
           else
             {
               std::vector<int> ret(1);
               ret[0]=MEDCoupling_MEDFileFields_getPosOfField(self,elts);
               return ret;
             }
         }

         void pushFields(PyObject *fields)
         {
           std::vector<MEDFileAnyTypeFieldMultiTS *> tmp;
           convertFromPyObjVectorOfObj<MEDCoupling::MEDFileAnyTypeFieldMultiTS *>(fields,SWIGTYPE_p_MEDCoupling__MEDFileAnyTypeFieldMultiTS,"MEDFileAnyTypeFieldMultiTS",tmp);
           self->pushFields(tmp);
         }
         
         void __delitem__(PyObject *elts)
         {
           if(elts && PySlice_Check(elts))
             {
               Py_ssize_t strt=2,stp=2,step=2;
               GetIndicesOfSlice(elts,self->getNumberOfFields(),&strt,&stp,&step,"MEDFileFields.__delitem__ : error in input slice !");
               self->destroyFieldsAtPos2(strt,stp,step);
             }
           else
             {
               std::vector<int> idsToRemove=MEDCoupling_MEDFileFields_getPosOfFields(self,elts);
               if(!idsToRemove.empty())
                 self->destroyFieldsAtPos(&idsToRemove[0],&idsToRemove[0]+idsToRemove.size());
             }
         }

         MEDFileFields *extractPart(PyObject *extractDef, MEDFileMesh *mm) const
         {
           std::map<int, MCAuto<DataArrayInt> > extractDefCpp;
           convertToMapIntDataArrayInt(extractDef,extractDefCpp);
           return self->extractPart(extractDefCpp,mm);
         }

         MEDFileFields *linearToQuadratic(const MEDFileMeshes *oldLin, const MEDFileMeshes *newQuad) const
         {
           MCAuto<MEDFileFields> ret(self->linearToQuadratic(oldLin,newQuad));
           return ret.retn();
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
    void setValue(double val);
    double getValue() const;
    std::string simpleRepr() const;
    %extend
    {
      std::string __str__() const
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
    static MEDFileParameterDouble1TS *New(const std::string& fileName);
    static MEDFileParameterDouble1TS *New(const std::string& fileName, const std::string& paramName);
    static MEDFileParameterDouble1TS *New(const std::string& fileName, const std::string& paramName, int dt, int it);
    virtual MEDFileParameter1TS *deepCopy() const;
    virtual std::string simpleRepr() const;
    void setName(const std::string& name);
    std::string getName() const;
    void write(const std::string& fileName, int mode) const;
    %extend
    {
      MEDFileParameterDouble1TS()
      {
        return MEDFileParameterDouble1TS::New();
      }
      
      MEDFileParameterDouble1TS(const std::string& fileName)
      {
        return MEDFileParameterDouble1TS::New(fileName);
      }

      MEDFileParameterDouble1TS(const std::string& fileName, const std::string& paramName)
      {
        return MEDFileParameterDouble1TS::New(fileName,paramName);
      }

      MEDFileParameterDouble1TS(const std::string& fileName, const std::string& paramName, int dt, int it)
      {
        return MEDFileParameterDouble1TS::New(fileName,paramName,dt,it);
      }

      std::string __str__() const
      {
        return self->simpleRepr();
      }

      PyObject *isEqual(const MEDFileParameter1TS *other, double eps) const
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
    static MEDFileParameterMultiTS *New(const std::string& fileName);
    static MEDFileParameterMultiTS *New(const std::string& fileName, const std::string& paramName);
    std::string getName() const;
    void setName(const std::string& name);
    MEDFileParameterMultiTS *deepCopy() const;
    void write(const std::string& fileName, int mode) const;
    std::string simpleRepr() const;
    void appendValue(int dt, int it, double time, double val);
    double getDoubleValue(int iteration, int order) const;
    int getPosOfTimeStep(int iteration, int order) const;
    int getPosGivenTime(double time, double eps=1e-8) const;
    int getNumberOfTS() const;
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

      std::string __str__() const
      {
        return self->simpleRepr();
      }
      
      PyObject *isEqual(const MEDFileParameterMultiTS *other, double eps) const
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
      
      void eraseTimeStepIds(PyObject *ids)
      {
        int sw;
        int pos1;
        std::vector<int> pos2;
        DataArrayInt *pos3=0;
        DataArrayIntTuple *pos4=0;
        convertIntStarLikePyObjToCpp(ids,sw,pos1,pos2,pos3,pos4);
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

      int getTimeStepId(PyObject *elt0) const
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

      MEDFileParameter1TS *__getitem__(PyObject *elt0) const
      {
        MEDFileParameter1TS *ret=self->getTimeStepAtPos(MEDCoupling_MEDFileParameterMultiTS_getTimeStepId(self,elt0));
        if(ret)
          ret->incrRef();
        return ret;
      }

      std::vector<int> getTimeStepIds(PyObject *elts) const
      {
        if(PyList_Check(elts))
          {
            int sz=PyList_Size(elts);
            std::vector<int> ret(sz);
            for(int i=0;i<sz;i++)
              {
                PyObject *elt=PyList_GetItem(elts,i);
                ret[i]=MEDCoupling_MEDFileParameterMultiTS_getTimeStepId(self,elt);
              }
            return ret;
          }
        else
          {
            std::vector<int> ret(1);
            ret[0]=MEDCoupling_MEDFileParameterMultiTS_getTimeStepId(self,elts);
            return ret;
          }
      }

      void __delitem__(PyObject *elts)
      {
        std::vector<int> idsToRemove=MEDCoupling_MEDFileParameterMultiTS_getTimeStepIds(self,elts);
        if(!idsToRemove.empty())
          self->eraseTimeStepIds(&idsToRemove[0],&idsToRemove[0]+idsToRemove.size());
      }
      
      MEDFileParameter1TS *getTimeStepAtPos(int posId) const
      {
        MEDFileParameter1TS *ret=self->getTimeStepAtPos(posId);
        if(ret)
          ret->incrRef();
        return ret;
      }

      PyObject *getIterations() const
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

      PyObject *getTimeSteps() const
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

  class MEDFileParameters : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    static MEDFileParameters *New();
    static MEDFileParameters *New(const std::string& fileName);
    static MEDFileParameters *New(DataArrayByte *db);
    MEDFileParameters *deepCopy() const;
    std::vector<std::string> getParamsNames() const;
    std::string simpleRepr() const;
    void resize(int newSize);
    void pushParam(MEDFileParameterMultiTS *param);
    void setParamAtPos(int i, MEDFileParameterMultiTS *param);
    void destroyParamAtPos(int i);
    int getPosFromParamName(const std::string& paramName) const;
    int getNumberOfParams() const;
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

      MEDFileParameters(DataArrayByte *db)
      {
        return MEDFileParameters::New(db);
      }
      
      std::string __str__() const
      {
        return self->simpleRepr();
      }

      MEDFileParameterMultiTS *__getitem__(PyObject *obj)
      {
        static const char msg[]="MEDFileParameters::__getitem__ : only integer or string with meshname supported !";
        if(PyInt_Check(obj))
          {
            MEDFileParameterMultiTS *ret=self->getParamAtPos(InterpreteNegativeInt((int)PyInt_AS_LONG(obj),self->getNumberOfParams()));
            if(ret)
              ret->incrRef();
            return ret;
          }
        MEDFileParameterMultiTS *ret(self->getParamWithName(convertPyObjectToStr(obj,msg)));
        if(ret)
          ret->incrRef();
        return ret;
      }

      int __len__() const
      {
        return self->getNumberOfParams();
      }
      
      MEDFileParameterMultiTS *getParamAtPos(int i) const
      {
        MEDFileParameterMultiTS *ret=self->getParamAtPos(i);
        if(ret)
          ret->incrRef();
        return ret;
      }

      MEDFileParameterMultiTS *getParamWithName(const std::string& paramName) const
      {
        MEDFileParameterMultiTS *ret=self->getParamWithName(paramName);
        if(ret)
          ret->incrRef();
        return ret;
      }
      
      PyObject *isEqual(const MEDFileParameters *other, double eps) const
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

  class MEDFileData : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    static MEDFileData *New(DataArrayByte *db);
    static MEDFileData *New(const std::string& fileName);
    static MEDFileData *New();
    MEDFileData *deepCopy() const;
    void setFields(MEDFileFields *fields);
    void setMeshes(MEDFileMeshes *meshes);
    void setParams(MEDFileParameters *params);
    int getNumberOfFields() const;
    int getNumberOfMeshes() const;
    int getNumberOfParams() const;
    //
    bool changeMeshName(const std::string& oldMeshName, const std::string& newMeshName);
    bool unPolyzeMeshes();
    void dealWithStructureElements();
    std::string getHeader() const;
    void setHeader(const std::string& header);
    //
    %extend
       {
         MEDFileData(const std::string& fileName)
         {
           return MEDFileData::New(fileName);
         }

         MEDFileData(DataArrayByte *db)
         {
           return MEDFileData::New(db);
         }

         MEDFileData()
         {
           return MEDFileData::New();
         }

         std::string __str__() const
         {
           return self->simpleRepr();
         }

         MEDFileMeshes *getMeshes() const
         {
           MEDFileMeshes *ret=self->getMeshes();
           if(ret)
             ret->incrRef();
           return ret;
         }

         MEDFileParameters *getParams() const
         {
           MEDFileParameters *ret=self->getParams();
           if(ret)
             ret->incrRef();
           return ret;
         }

         MEDFileFields *getFields() const
         {
           MEDFileFields *ret=self->getFields();
           if(ret)
             ret->incrRef();
           return ret;
         }

         bool changeMeshNames(PyObject *li)
         {
           std::vector< std::pair<std::string,std::string> > modifTab=convertVecPairStStFromPy(li);
           return self->changeMeshNames(modifTab);
         }

         static MEDFileData *Aggregate(PyObject *mfds)
         {
           std::vector<const MEDFileData *> mfdsCpp;
           convertFromPyObjVectorOfObj<const MEDCoupling::MEDFileData *>(mfds,SWIGTYPE_p_MEDCoupling__MEDFileData,"MEDFileData",mfdsCpp);
           MCAuto<MEDFileData> ret(MEDFileData::Aggregate(mfdsCpp));
           return ret.retn();
         }
       }
  };

  class SauvReader : public RefCountObject
  {
  public:
    static SauvReader* New(const std::string& fileName);
    MEDFileData * loadInMEDFileDS();
    %extend
    {
      SauvReader(const std::string& fileName)
      {
        return SauvReader::New(fileName);
      }
    }
  };

  class SauvWriter : public RefCountObject
  {
  public:
    static SauvWriter * New();
    void setMEDFileDS(const MEDFileData* medData, unsigned meshIndex = 0);
    void write(const std::string& fileName);
    void setCpyGrpIfOnASingleFamilyStatus(bool status);
    bool getCpyGrpIfOnASingleFamilyStatus() const;
    %extend
    {
      SauvWriter()
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
    static MEDFileField1TSStructItem BuildItemFrom(const MEDFileAnyTypeField1TS *ref, const MEDFileMeshStruct *meshSt);
  };

  class MEDFileMeshStruct : public RefCountObject
  {
  public:
    static MEDFileMeshStruct *New(const MEDFileMesh *mesh);
  protected:
    ~MEDFileMeshStruct();
  };
  
  class MEDMeshMultiLev : public RefCountObject
  {
  public:
    virtual MEDMeshMultiLev *prepare() const;
    DataArray *buildDataArray(const MEDFileField1TSStructItem& fst, const MEDFileFieldGlobsReal *globs, const DataArray *vals) const;
    DataArrayInt *retrieveGlobalNodeIdsIfAny() const;
  protected:
    ~MEDMeshMultiLev();
  public:
    %extend
    {
      PyObject *retrieveFamilyIdsOnCells() const
      {
        DataArrayInt *famIds(0);
        bool isWithoutCopy(false);
        self->retrieveFamilyIdsOnCells(famIds,isWithoutCopy);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret1Py=isWithoutCopy?Py_True:Py_False;
        Py_XINCREF(ret1Py);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(famIds),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,ret1Py);
        return ret;
      }

      PyObject *retrieveNumberIdsOnCells() const
      {
        DataArrayInt *numIds(0);
        bool isWithoutCopy(false);
        self->retrieveNumberIdsOnCells(numIds,isWithoutCopy);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret1Py=isWithoutCopy?Py_True:Py_False;
        Py_XINCREF(ret1Py);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(numIds),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,ret1Py);
        return ret;
      }
      
      PyObject *retrieveFamilyIdsOnNodes() const
      {
        DataArrayInt *famIds(0);
        bool isWithoutCopy(false);
        self->retrieveFamilyIdsOnNodes(famIds,isWithoutCopy);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret1Py=isWithoutCopy?Py_True:Py_False;
        Py_XINCREF(ret1Py);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(famIds),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,ret1Py);
        return ret;
      }

      PyObject *retrieveNumberIdsOnNodes() const
      {
        DataArrayInt *numIds(0);
        bool isWithoutCopy(false);
        self->retrieveNumberIdsOnNodes(numIds,isWithoutCopy);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret1Py=isWithoutCopy?Py_True:Py_False;
        Py_XINCREF(ret1Py);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(numIds),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,ret1Py);
        return ret;
      }

      PyObject *getGeoTypes() const
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
       PyObject *buildVTUArrays() const
       {
         DataArrayDouble *coords(0); DataArrayByte *types(0); DataArrayInt *cellLocations(0),*cells(0),*faceLocations(0),*faces(0);
         bool ncc(self->buildVTUArrays(coords,types,cellLocations,cells,faceLocations,faces));
         PyObject *ret0Py=ncc?Py_True:Py_False;
         Py_XINCREF(ret0Py);
         PyObject *ret=PyTuple_New(7);
         PyTuple_SetItem(ret,0,ret0Py);
         PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(coords),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
         PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(types),SWIGTYPE_p_MEDCoupling__DataArrayByte, SWIG_POINTER_OWN | 0 ));
         PyTuple_SetItem(ret,3,SWIG_NewPointerObj(SWIG_as_voidptr(cellLocations),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
         PyTuple_SetItem(ret,4,SWIG_NewPointerObj(SWIG_as_voidptr(cells),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
         PyTuple_SetItem(ret,5,SWIG_NewPointerObj(SWIG_as_voidptr(faceLocations),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
         PyTuple_SetItem(ret,6,SWIG_NewPointerObj(SWIG_as_voidptr(faces),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
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
      PyObject *buildVTUArrays() const
      {
        bool isInternal;
        std::vector< DataArrayDouble * > objs(self->buildVTUArrays(isInternal));
        std::size_t sz(objs.size());
        PyObject *ret(PyTuple_New(2));
        PyObject *ret0=PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(ret0,i,SWIG_NewPointerObj(SWIG_as_voidptr(objs[i]),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
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
      PyObject *buildVTUArrays() const
      {
        DataArrayDouble *ret0(0);
        std::vector<int> ret1;
        bool ret2;
        self->buildVTUArrays(ret0,ret1,ret2);
        std::size_t sz(ret1.size());
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
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
    static MEDFileFastCellSupportComparator *New(const MEDFileMeshStruct *m, const MEDFileAnyTypeFieldMultiTS *ref);
    MEDMeshMultiLev *buildFromScratchDataSetSupport(int timeStepId, const MEDFileFieldGlobsReal *globs) const;
    bool isDataSetSupportEqualToThePreviousOne(int timeStepId, const MEDFileFieldGlobsReal *globs) const;
    int getNumberOfTS() const;
  protected:
    ~MEDFileFastCellSupportComparator();
  public:
    %extend
    {
      PyObject *getGeoTypesAt(int timeStepId, const MEDFileMesh *m) const
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

%pythoncode %{
def enter1TS(self):
    self.loadArrays()
    pass
def exit1TS(self, exctype, exc, tb):
    self.unloadArrays()
    pass
MEDFileAnyTypeField1TS.__enter__=enter1TS
MEDFileAnyTypeField1TS.__exit__=exit1TS
%}
