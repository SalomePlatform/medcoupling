// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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

%module MEDLoader

#define MEDCOUPLING_EXPORT
#define MEDLOADER_EXPORT

%include "MEDCouplingCommon.i"

%{
#include "MEDLoader.hxx"
#include "MEDFileMesh.hxx"
#include "MEDFileField.hxx"
#include "MEDFileData.hxx"
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

%newobject MEDLoader::ReadUMeshFromFamilies;
%newobject MEDLoader::ReadUMeshFromGroups;
%newobject MEDLoader::ReadUMeshFromFile;
%newobject MEDLoader::ReadField;
%newobject MEDLoader::ReadFieldCell;
%newobject MEDLoader::ReadFieldNode;
%newobject MEDLoader::ReadFieldGauss;
%newobject MEDLoader::ReadFieldGaussNE;
%newobject ParaMEDMEM::MEDFileMesh::New;
%newobject ParaMEDMEM::MEDFileMesh::getGenMeshAtLevel;
%newobject ParaMEDMEM::MEDFileMesh::getGroupArr;
%newobject ParaMEDMEM::MEDFileMesh::getGroupsArr;
%newobject ParaMEDMEM::MEDFileMesh::getFamilyArr;
%newobject ParaMEDMEM::MEDFileMesh::getFamiliesArr;
%newobject ParaMEDMEM::MEDFileMesh::getNodeGroupArr;
%newobject ParaMEDMEM::MEDFileMesh::getNodeGroupsArr;
%newobject ParaMEDMEM::MEDFileMesh::getNodeFamilyArr;
%newobject ParaMEDMEM::MEDFileMesh::getNodeFamiliesArr;
%newobject ParaMEDMEM::MEDFileMesh::getAllFamiliesIdsReferenced;
%newobject ParaMEDMEM::MEDFileUMesh::New;
%newobject ParaMEDMEM::MEDFileUMesh::getCoords;
%newobject ParaMEDMEM::MEDFileUMesh::getGroup;
%newobject ParaMEDMEM::MEDFileUMesh::getGroups;
%newobject ParaMEDMEM::MEDFileUMesh::getFamily;
%newobject ParaMEDMEM::MEDFileUMesh::getFamilies;
%newobject ParaMEDMEM::MEDFileUMesh::getMeshAtLevel;
%newobject ParaMEDMEM::MEDFileUMesh::getLevel0Mesh;
%newobject ParaMEDMEM::MEDFileUMesh::getLevelM1Mesh;
%newobject ParaMEDMEM::MEDFileUMesh::getLevelM2Mesh;
%newobject ParaMEDMEM::MEDFileUMesh::getLevelM3Mesh;
%newobject ParaMEDMEM::MEDFileCMesh::New;
%newobject ParaMEDMEM::MEDFileMeshMultiTS::New;
%newobject ParaMEDMEM::MEDFileMeshMultiTS::getOneTimeStep;
%newobject ParaMEDMEM::MEDFileMeshes::New;
%newobject ParaMEDMEM::MEDFileMeshes::getMeshAtPos;
%newobject ParaMEDMEM::MEDFileMeshes::getMeshWithName;
%newobject ParaMEDMEM::MEDFileMeshes::__getitem__;
%newobject ParaMEDMEM::MEDFileMeshes::__iter__;

%newobject ParaMEDMEM::MEDFileFields::New;
%newobject ParaMEDMEM::MEDFileFields::getFieldWithName;
%newobject ParaMEDMEM::MEDFileFields::getFieldAtPos;
%newobject ParaMEDMEM::MEDFileFields::__getitem__;
%newobject ParaMEDMEM::MEDFileFields::__iter__;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::New;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::getTimeStepAtPos;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::getTimeStep;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::getTimeStepGivenTime;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::__getitem__;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::__iter__;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::getFieldAtLevel;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::getFieldAtTopLevel;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::getFieldOnMeshAtLevel;
%newobject ParaMEDMEM::MEDFileFieldMultiTS::getFieldAtLevelOld;
%newobject ParaMEDMEM::MEDFileField1TSWithoutDAS::getUndergroundDataArray;
%newobject ParaMEDMEM::MEDFileField1TS::New;
%newobject ParaMEDMEM::MEDFileField1TS::getFieldAtLevel;
%newobject ParaMEDMEM::MEDFileField1TS::getFieldAtTopLevel;
%newobject ParaMEDMEM::MEDFileField1TS::getFieldOnMeshAtLevel;
%newobject ParaMEDMEM::MEDFileField1TS::getFieldAtLevelOld;
%newobject ParaMEDMEM::MEDFileFieldMultiTSWithoutDAS::getUndergroundDataArray;

%newobject ParaMEDMEM::MEDFileData::New;
%newobject ParaMEDMEM::MEDFileData::getMeshes;
%newobject ParaMEDMEM::MEDFileData::getFields;

%newobject ParaMEDMEM::SauvWriter::New;
%newobject ParaMEDMEM::SauvReader::New;
%newobject ParaMEDMEM::SauvReader::loadInMEDFileDS;

%feature("unref") MEDFileMesh "$this->decrRef();"
%feature("unref") MEDFileUMesh "$this->decrRef();"
%feature("unref") MEDFileCMesh "$this->decrRef();"
%feature("unref") MEDFileMeshMultiTS "$this->decrRef();"
%feature("unref") MEDFileMeshes "$this->decrRef();"
%feature("unref") MEDFileFieldLoc "$this->decrRef();"
%feature("unref") MEDFileField1TS "$this->decrRef();"
%feature("unref") MEDFileFieldMultiTS "$this->decrRef();"
%feature("unref") MEDFileFields "$this->decrRef();"
%feature("unref") MEDFileData "$this->decrRef();"
%feature("unref") SauvReader "$this->decrRef();"
%feature("unref") SauvWriter "$this->decrRef();"

class MEDLoader
{
public:
  static void setEpsilonForNodeComp(double val) throw(INTERP_KERNEL::Exception);
  static void setCompPolicyForCell(int val) throw(INTERP_KERNEL::Exception);
  static void setTooLongStrPolicy(int val) throw(INTERP_KERNEL::Exception);
  static void CheckFileForRead(const char *fileName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshNames(const char *fileName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshNamesOnField(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshGroupsNames(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshFamiliesNames(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshFamiliesNamesOnGroup(const char *fileName, const char *meshName, const char *grpName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshGroupsNamesOnFamily(const char *fileName, const char *meshName, const char *famName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetAllFieldNamesOnMesh(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetAllFieldNames(const char *fileName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetFieldNamesOnMesh(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetCellFieldNamesOnMesh(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetNodeFieldNamesOnMesh(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception);
  static double GetTimeAttachedOnFieldIteration(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  %extend
     {
       static PyObject *GetFieldIterations(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, const char *fieldName) throw(INTERP_KERNEL::Exception)
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

       static PyObject *GetAllFieldIterations(const char *fileName, const char *meshName, const char *fieldName) throw(INTERP_KERNEL::Exception)
       {
         std::vector< std::pair< std::pair<int,int>, double> > res=MEDLoader::GetAllFieldIterations(fileName,meshName,fieldName);
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

       static PyObject *GetCellFieldIterations(const char *fileName, const char *meshName, const char *fieldName) throw(INTERP_KERNEL::Exception)
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
       static PyObject *GetNodeFieldIterations(const char *fileName, const char *meshName, const char *fieldName) throw(INTERP_KERNEL::Exception)
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
       static PyObject *GetComponentsNamesOfField(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
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
       static PyObject *GetUMeshGlobalInfo(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception)
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
       static PyObject *ReadFieldsOnSameMesh(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, int meshDimRelToMax,
                                             const char *fieldName, PyObject *liIts) throw(INTERP_KERNEL::Exception)
       {
         std::vector<std::pair<int,int> > its=convertTimePairIdsFromPy(liIts);
         std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> res=MEDLoader::ReadFieldsOnSameMesh(type,fileName,meshName,meshDimRelToMax,fieldName,its);
         return convertFieldDoubleVecToPy(res);
       }
       static void WriteUMeshesPartition(const char *fileName, const char *meshName, PyObject *li, bool writeFromScratch) throw(INTERP_KERNEL::Exception)
       {
         std::vector<const ParaMEDMEM::MEDCouplingUMesh *> v=convertUMeshVecFromPy(li);
         MEDLoader::WriteUMeshesPartition(fileName,meshName,v,writeFromScratch);
       }
       static void WriteUMeshesPartitionDep(const char *fileName, const char *meshName, PyObject *li, bool writeFromScratch) throw(INTERP_KERNEL::Exception)
       {
         std::vector<const ParaMEDMEM::MEDCouplingUMesh *> v=convertUMeshVecFromPy(li);
         MEDLoader::WriteUMeshesPartitionDep(fileName,meshName,v,writeFromScratch);
       }
       static void WriteUMeshes(const char *fileName, PyObject *li, bool writeFromScratch) throw(INTERP_KERNEL::Exception)
       {
         std::vector<const ParaMEDMEM::MEDCouplingUMesh *> v=convertUMeshVecFromPy(li);
         MEDLoader::WriteUMeshes(fileName,v,writeFromScratch);
       }
       static PyObject *GetTypesOfField(const char *fileName, const char *meshName, const char *fieldName) throw(INTERP_KERNEL::Exception)
       {
         std::vector< ParaMEDMEM::TypeOfField > v=MEDLoader::GetTypesOfField(fileName,meshName,fieldName);
         int size=v.size();
         PyObject *ret=PyList_New(size);
         for(int i=0;i<size;i++)
           PyList_SetItem(ret,i,PyInt_FromLong((int)v[i]));
         return ret;
       }
       static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromGroups(const char *fileName, const char *meshName, int meshDimRelToMax, PyObject *li) throw(INTERP_KERNEL::Exception)
       {
         std::vector<std::string> grps;
         converPyListToVecString(li,grps);
         return MEDLoader::ReadUMeshFromGroups(fileName,meshName,meshDimRelToMax,grps);
       }
       static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFamilies(const char *fileName, const char *meshName, int meshDimRelToMax, PyObject *li) throw(INTERP_KERNEL::Exception)
       {
         std::vector<std::string> fams;
         converPyListToVecString(li,fams);
         return MEDLoader::ReadUMeshFromFamilies(fileName,meshName,meshDimRelToMax,fams);
       }
     }
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const char *fileName, const char *meshName, int meshDimRelToMax=0) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const char *fileName, int meshDimRelToMax=0) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadField(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldCell(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldNode(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldGauss(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldGaussNE(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static void WriteUMesh(const char *fileName, const ParaMEDMEM::MEDCouplingUMesh *mesh, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteUMeshDep(const char *fileName, const ParaMEDMEM::MEDCouplingUMesh *mesh, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteField(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteFieldDep(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteFieldUsingAlreadyWrittenMesh(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f) throw(INTERP_KERNEL::Exception);
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

  class MEDFileMesh : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileMesh *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileMesh *New(const char *fileName, const char *mName, int dt=-1, int it=-1);
    virtual void clearNonDiscrAttributes() const;
    void setName(const char *name);
    const char *getName();
    void setUnivName(const char *name);
    const char *getUnivName() const;
    void setDescription(const char *name);
    const char *getDescription() const;
    void setOrder(int order);
    int getOrder() const;
    void setIteration(int it);
    int getIteration();
    void setTimeValue(double time);
    void setTime(int dt, int it, double time);
    double getTimeValue() const;
    void setTimeUnit(const char *unit);
    const char *getTimeUnit() const;
    virtual int getNumberOfNodes() const throw(INTERP_KERNEL::Exception);
    std::vector<int> getNonEmptyLevels() const;
    std::vector<int> getNonEmptyLevelsExt() const;
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    int getSizeAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    //
    bool existsGroup(const char *groupName) const;
    bool existsFamily(int famId) const;
    bool existsFamily(const char *familyName) const;
    void setFamilyId(const char *familyName, int id);
    void setFamilyIdUnique(const char *familyName, int id) throw(INTERP_KERNEL::Exception);
    void addFamily(const char *familyName, int id) throw(INTERP_KERNEL::Exception);
    void addFamilyOnGrp(const char *grpName, const char *famName) throw(INTERP_KERNEL::Exception);
    virtual void createGroupOnAll(int meshDimRelToMaxExt, const char *groupName) throw(INTERP_KERNEL::Exception);
    virtual bool keepFamIdsOnlyOnLevs(const std::vector<int>& famIds, const std::vector<int>& levs) throw(INTERP_KERNEL::Exception);
    void copyFamGrpMapsFrom(const MEDFileMesh& other);
    const std::map<std::string,int>& getFamilyInfo() const;
    const std::map<std::string, std::vector<std::string> >& getGroupInfo() const;
    std::vector<std::string> getFamiliesOnGroup(const char *name) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getFamiliesOnGroups(const std::vector<std::string>& grps) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamiliesIdsOnGroup(const char *name) const throw(INTERP_KERNEL::Exception);
    void setFamiliesOnGroup(const char *name, const std::vector<std::string>& fams) throw(INTERP_KERNEL::Exception);
    void setFamiliesIdsOnGroup(const char *name, const std::vector<int>& famIds) throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getGroupsOnFamily(const char *name) const throw(INTERP_KERNEL::Exception);
    void setGroupsOnFamily(const char *famName, const std::vector<std::string>& grps) throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getGroupsNames() const;
    std::vector<std::string> getFamiliesNames() const;
    void assignFamilyNameWithGroupName() throw(INTERP_KERNEL::Exception);
    void removeGroup(const char *name) throw(INTERP_KERNEL::Exception);
    void removeFamily(const char *name) throw(INTERP_KERNEL::Exception);
    void changeGroupName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception);
    void changeFamilyName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception);
    void changeFamilyId(int oldId, int newId) throw(INTERP_KERNEL::Exception);
    void setFamilyInfo(const std::map<std::string,int>& info);
    void setGroupInfo(const std::map<std::string, std::vector<std::string> >&info);
    int getFamilyId(const char *name) const throw(INTERP_KERNEL::Exception);
    int getMaxFamilyId() const throw(INTERP_KERNEL::Exception);
    int getMinFamilyId() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getAllFamiliesIdsReferenced() const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamiliesIds(const std::vector<std::string>& famNames) const throw(INTERP_KERNEL::Exception);
    std::string getFamilyNameGivenId(int id) const throw(INTERP_KERNEL::Exception);
    bool ensureDifferentFamIdsPerLevel() throw(INTERP_KERNEL::Exception);
    void normalizeFamIdsTrio() throw(INTERP_KERNEL::Exception);
    void normalizeFamIdsMEDFile() throw(INTERP_KERNEL::Exception);
    virtual int getMeshDimension() const throw(INTERP_KERNEL::Exception);
    virtual std::string simpleRepr() const;
    virtual std::string advancedRepr() const;
    //
    virtual MEDCouplingMesh *getGenMeshAtLevel(int meshDimRelToMax, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual void setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr) throw(INTERP_KERNEL::Exception);
    virtual void setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr) throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getGroupsArr(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getGroupArr(int meshDimRelToMaxExt, const char *grp, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getFamilyArr(int meshDimRelToMaxExt, const char *fam, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeGroupArr(const char *grp, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeGroupsArr(const std::vector<std::string>& grps, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeFamilyArr(const char *fam, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeFamiliesArr(const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception);
    %extend
       {
         MEDFileMesh(const char *fileName) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileMesh::New(fileName);
         }

         MEDFileMesh(const char *fileName, const char *mName, int dt=-1, int it=-1)
         {
           return MEDFileMesh::New(fileName,mName,dt,it);
         }
         
         std::string __str__() const
         {
           return self->simpleRepr();
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

         PyObject *getFamilyFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
         {
           const DataArrayInt *tmp=self->getFamilyFieldAtLevel(meshDimRelToMaxExt);
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

         PyObject *findOrCreateAndGiveFamilyWithId(int id, bool& created) throw(INTERP_KERNEL::Exception)
         {
           bool ret1;
           std::string ret0=self->findOrCreateAndGiveFamilyWithId(id,ret1);
           PyObject *ret=PyTuple_New(2);
           PyTuple_SetItem(ret,0,PyString_FromString(ret0.c_str()));
           PyTuple_SetItem(ret,1,SWIG_From_bool(ret1));
           return ret;
         }
       }
  };

  class MEDFileUMesh : public MEDFileMesh
  {
  public:
    static MEDFileUMesh *New(const char *fileName, const char *mName, int dt=-1, int it=-1) throw(INTERP_KERNEL::Exception);
    static MEDFileUMesh *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileUMesh *New();
    ~MEDFileUMesh();
    int getSpaceDimension() const throw(INTERP_KERNEL::Exception);
    //
    std::vector<int> getGrpNonEmptyLevels(const char *grp) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getGrpNonEmptyLevelsExt(const char *grp) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamNonEmptyLevels(const char *fam) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamNonEmptyLevelsExt(const char *fam) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getGrpsNonEmptyLevels(const std::vector<std::string>& grps) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getGrpsNonEmptyLevelsExt(const std::vector<std::string>& grps) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamsNonEmptyLevels(const std::vector<std::string>& fams) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamsNonEmptyLevelsExt(const std::vector<std::string>& fams) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getGroupsOnSpecifiedLev(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getGroup(int meshDimRelToMaxExt, const char *grp, bool renum=false) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getGroupArr(int meshDimRelToMaxExt, const char *grp, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getGroups(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum=false) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getGroupsArr(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getFamily(int meshDimRelToMaxExt, const char *fam, bool renum=false) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getFamilyArr(int meshDimRelToMaxExt, const char *fam, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getFamilies(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getNodeGroupArr(const char *grp, bool renum=false) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getNodeGroupsArr(const std::vector<std::string>& grps, bool renum=false) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getNodeFamilyArr(const char *fam, bool renum=false) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getNodeFamiliesArr(const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getMeshAtLevel(int meshDimRelToMaxExt, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getLevel0Mesh(bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getLevelM1Mesh(bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getLevelM2Mesh(bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getLevelM3Mesh(bool renum=false) const throw(INTERP_KERNEL::Exception);
    //
    void setFamilyNameAttachedOnId(int id, const std::string& newFamName) throw(INTERP_KERNEL::Exception);
    void setCoords(DataArrayDouble *coords) throw(INTERP_KERNEL::Exception);
    void eraseGroupsAtLevel(int meshDimRelToMaxExt) throw(INTERP_KERNEL::Exception);
    void addNodeGroup(const std::string& name, const std::vector<int>& ids) throw(INTERP_KERNEL::Exception);
    void removeMeshAtLevel(int meshDimRelToMax) throw(INTERP_KERNEL::Exception);
    void setMeshAtLevel(int meshDimRelToMax, MEDCouplingUMesh *m, bool newOrOld=false) throw(INTERP_KERNEL::Exception);
    void setMeshAtLevelGen(int meshDimRelToMax, MEDCouplingUMesh *m, bool newOrOld) throw(INTERP_KERNEL::Exception);
    void setGroupsFromScratch(int meshDimRelToMax, const std::vector<const MEDCouplingUMesh *>& ms) throw(INTERP_KERNEL::Exception);
    void setGroupsOnSetMesh(int meshDimRelToMax, const std::vector<const MEDCouplingUMesh *>& ms, bool renum) throw(INTERP_KERNEL::Exception);
    void optimizeFamilies() throw(INTERP_KERNEL::Exception);
    %extend
       { 
         MEDFileUMesh(const char *fileName, const char *mName, int dt=-1, int it=-1) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileUMesh::New(fileName,mName,dt,it);
         }

         MEDFileUMesh(const char *fileName) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileUMesh::New(fileName);
         }

         MEDFileUMesh()
         {
           return MEDFileUMesh::New();
         }
         
         PyObject *getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception)
         {
           const DataArrayInt *tmp=self->getRevNumberFieldAtLevel(meshDimRelToMaxExt);
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 );
         }
         
         void setGroupsAtLevel(int meshDimRelToMaxExt, PyObject *li, bool renum=false) throw(INTERP_KERNEL::Exception)
         {
           std::vector<const DataArrayInt *> grps;
           convertPyObjToVecDataArrayIntCst(li,grps);
           self->setGroupsAtLevel(meshDimRelToMaxExt,grps,renum);
         }

         DataArrayDouble *getCoords() const
         {
           DataArrayDouble *ret=self->getCoords();
           if(ret)
             ret->incrRef();
           return ret;
         }

         PyObject *duplicateNodesOnM1Group(const char *grpNameM1) throw(INTERP_KERNEL::Exception)
         {
           DataArrayInt *ret0=0,*ret1=0,*ret2=0;
           self->duplicateNodesOnM1Group(grpNameM1,ret0,ret1,ret2);
           PyObject *ret=PyTuple_New(3);
           PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           PyTuple_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(ret2),SWIGTYPE_p_ParaMEDMEM__DataArrayInt, SWIG_POINTER_OWN | 0 ));
           return ret;
         }
       }
  };

  class MEDFileCMesh : public MEDFileMesh
  {
  public:
    static MEDFileCMesh *New();
    static MEDFileCMesh *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileCMesh *New(const char *fileName, const char *mName, int dt=-1, int it=-1) throw(INTERP_KERNEL::Exception);
    void setMesh(MEDCouplingCMesh *m) throw(INTERP_KERNEL::Exception);
    %extend
       {
         MEDFileCMesh()
         {
           return MEDFileCMesh::New();
         }

         MEDFileCMesh(const char *fileName) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileCMesh::New(fileName);
         }

         MEDFileCMesh(const char *fileName, const char *mName, int dt=-1, int it=-1) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileCMesh::New(fileName,mName,dt,it);
         }
         
         PyObject *getMesh() const
         {
           const MEDCouplingCMesh *tmp=self->getMesh();
           if(tmp)
             tmp->incrRef();
           return SWIG_NewPointerObj(SWIG_as_voidptr(tmp),SWIGTYPE_p_ParaMEDMEM__MEDCouplingCMesh, SWIG_POINTER_OWN | 0 );
         }
       }
  };

  class MEDFileMeshMultiTS : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileMeshMultiTS *New();
    static MEDFileMeshMultiTS *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileMeshMultiTS *New(const char *fileName, const char *mName) throw(INTERP_KERNEL::Exception);
    const char *getName() const throw(INTERP_KERNEL::Exception);
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    void setOneTimeStep(MEDFileMesh *mesh1TimeStep) throw(INTERP_KERNEL::Exception);
    %extend
       { 
         MEDFileMeshMultiTS()
         {
           return MEDFileMeshMultiTS::New();
         }

         MEDFileMeshMultiTS(const char *fileName) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileMeshMultiTS::New(fileName);
         }

         MEDFileMeshMultiTS(const char *fileName, const char *mName) throw(INTERP_KERNEL::Exception)
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
    static MEDFileMeshes *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
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

         MEDFileMeshes(const char *fileName) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileMeshes::New(fileName);
         }

         std::string __str__() const
           {
             return self->simpleRepr();
           }

         MEDFileMesh *__getitem__(PyObject *obj) throw(INTERP_KERNEL::Exception)
         {
           if(PyInt_Check(obj))
             {
               MEDFileMesh *ret=self->getMeshAtPos((int)PyInt_AS_LONG(obj));
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
         
         MEDFileMesh *getMeshAtPos(int i) const throw(INTERP_KERNEL::Exception)
           {
             MEDFileMesh *ret=self->getMeshAtPos(i);
             if(ret)
               ret->incrRef();
             return ret;
           }
         MEDFileMesh *getMeshWithName(const char *mname) const throw(INTERP_KERNEL::Exception)
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
    const std::string& getName() const;
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
    void shallowCpyGlobs(const MEDFileFieldGlobsReal& other);
    std::vector<std::string> getPfls() const;
    std::vector<std::string> getLocs() const;
    virtual std::vector<std::string> getPflsReallyUsed() const;
    virtual std::vector<std::string> getLocsReallyUsed() const;
    virtual std::vector<std::string> getPflsReallyUsedMulti() const;
    virtual std::vector<std::string> getLocsReallyUsedMulti() const;
    void killProfileIds(const std::vector<int>& pflIds) throw(INTERP_KERNEL::Exception);
    void killLocalizationIds(const std::vector<int>& locIds) throw(INTERP_KERNEL::Exception);
    void changePflName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception);
    void changeLocName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception);
  %extend
     {
       PyObject *getProfile(const char *pflName) const throw(INTERP_KERNEL::Exception)
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
         loc->incrRef();
         return SWIG_NewPointerObj(SWIG_as_voidptr(loc),SWIGTYPE_p_ParaMEDMEM__MEDFileFieldLoc, SWIG_POINTER_OWN | 0 );
       }
       
       PyObject *getLocalization(const char *locName) const throw(INTERP_KERNEL::Exception)
       {
         const MEDFileFieldLoc *loc=&self->getLocalization(locName);
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
     }
  };

  class MEDFileField1TSWithoutDAS : public RefCountObject
  {
  public:
    void copyTinyInfoFrom(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception);
    //
    int getDimension() const;
    int getIteration() const;
    int getOrder() const;
    std::string getName();
    void setName(const char *name);
    std::string getMeshName();
    void setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception);
    int getNumberOfComponents() const;
    bool isDealingTS(int iteration, int order) const;
    const std::vector<std::string>& getInfo() const;
    void setTime(int iteration, int order, double val);
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

         PyObject *getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception)
         {
           std::vector<TypeOfField> ret=self->getTypesOfFieldAvailable();
           PyObject *ret2=PyList_New(ret.size());
           for(int i=0;i<(int)ret.size();i++)
             PyList_SetItem(ret2,i,SWIG_From_int(ret[i]));
           return ret2;
         }

         PyObject *getNonEmptyLevels(const char *mname=0) const throw(INTERP_KERNEL::Exception)
         {
           std::vector<int> ret1;
           int ret0=self->getNonEmptyLevels(mname,ret1);
           PyObject *elt=PyTuple_New(2);
           PyTuple_SetItem(elt,0,SWIG_From_int(ret0));
           PyTuple_SetItem(elt,1,convertIntArrToPyList2(ret1));
           return elt;
         }

         PyObject *getFieldSplitedByType(const char *mname=0) const throw(INTERP_KERNEL::Exception)
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

         PyObject *getFieldSplitedByType2(const char *mname=0) const throw(INTERP_KERNEL::Exception)
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
         
         bool changeMeshNames(PyObject *li) throw(INTERP_KERNEL::Exception)
         {
           std::vector< std::pair<std::string,std::string> > modifTab=convertVecPairStStFromPy(li);
           return self->changeMeshNames(modifTab);
         }
       }
  };

  class MEDFileField1TS : public MEDFileField1TSWithoutDAS, public MEDFileFieldGlobsReal, public MEDFileWritable
  {
  public:
    static MEDFileField1TS *New(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
    static MEDFileField1TS *New();
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtLevelOld(TypeOfField type, const char *mname, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    //
    void setFieldNoProfileSBT(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception);
    void setFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception);
    void setProfileNameOnLeaf(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const char *newPflName, bool forceRenameOnGlob=false) throw(INTERP_KERNEL::Exception);
    void setLocNameOnLeaf(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const char *newLocName, bool forceRenameOnGlob=false) throw(INTERP_KERNEL::Exception);
    %extend
       {
         MEDFileField1TS(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileField1TS::New(fileName,fieldName,iteration,order);
         }

         MEDFileField1TS()
         {
           return MEDFileField1TS::New();
         }
         
         std::string __str__() const
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

         void setProfileNameOnLeaf(INTERP_KERNEL::NormalizedCellType typ, int locId, const char *newPflName, bool forceRenameOnGlob=false) throw(INTERP_KERNEL::Exception)
           {
             self->setProfileNameOnLeaf(0,typ,locId,newPflName,forceRenameOnGlob);
           }
         
         void setLocNameOnLeaf(INTERP_KERNEL::NormalizedCellType typ, int locId, const char *newLocName, bool forceRenameOnGlob=false) throw(INTERP_KERNEL::Exception)
           {
             self->setLocNameOnLeaf(0,typ,locId,newLocName,forceRenameOnGlob);
           }
       }
  };

  class MEDFileFieldMultiTSWithoutDAS
  {
  public:
    int getNumberOfTS() const;
    void eraseEmptyTS() throw(INTERP_KERNEL::Exception);
    int getPosOfTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    int getPosGivenTime(double time, double eps=1e-8) const throw(INTERP_KERNEL::Exception);
    std::string getName() const;
    void setName(const char *name);
    std::string getMeshName() const throw(INTERP_KERNEL::Exception);
    void setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception);
    const std::vector<std::string>& getInfo() const;
    %extend
       {
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

         PyObject *getNonEmptyLevels(int iteration, int order, const char *mname=0) const throw(INTERP_KERNEL::Exception)
         {
           std::vector<int> ret1;
           int ret0=self->getNonEmptyLevels(iteration,order,mname,ret1);
           PyObject *elt=PyTuple_New(2);
           PyTuple_SetItem(elt,0,SWIG_From_int(ret0));
           PyTuple_SetItem(elt,1,convertIntArrToPyList2(ret1));
           return elt;
         }

         PyObject *getFieldSplitedByType(int iteration, int order, const char *mname=0) const throw(INTERP_KERNEL::Exception)
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

         PyObject *getFieldSplitedByType2(int iteration, int order, const char *mname=0) const throw(INTERP_KERNEL::Exception)
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
                     throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::__getitem__ : invalid input param ! input is a tuple of size 2 but two integers are expected in this tuple to request a time steps !");
                 }
               else
                 throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::__getitem__ : invalid input param ! input is a tuple of size != 2 ! two integers are expected in this tuple to request a time steps !");
             }
           else if(elt0 && PyFloat_Check(elt0))
             {
               double val=PyFloat_AS_DOUBLE(elt0);
               return self->getPosGivenTime(val);
             }
           else
             throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::__getitem__ : invalid input params ! expected fmts[int], fmts[int,int] or fmts[double] to request time step !");
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
                   ret[i]=ParaMEDMEM_MEDFileFieldMultiTSWithoutDAS_getTimeId(self,elt);
                 }
               return ret;
             }
           else
             {
               std::vector<int> ret(1);
               ret[0]=ParaMEDMEM_MEDFileFieldMultiTSWithoutDAS_getTimeId(self,elts);
               return ret;
             }
         }

         void __delitem__(PyObject *elts) throw(INTERP_KERNEL::Exception)
         {
           std::vector<int> idsToRemove=ParaMEDMEM_MEDFileFieldMultiTSWithoutDAS_getTimeIds(self,elts);
           if(!idsToRemove.empty())
             self->eraseTimeStepIds(&idsToRemove[0],&idsToRemove[0]+idsToRemove.size());
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
               throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutDAS::eraseTimeStepIds : unexpected input array type recognized !");
             }
         }

         bool changeMeshNames(PyObject *li) throw(INTERP_KERNEL::Exception)
         {
           std::vector< std::pair<std::string,std::string> > modifTab=convertVecPairStStFromPy(li);
           return self->changeMeshNames(modifTab);
         }
       }
  };

  class MEDFileFieldMultiTSIterator
  {
  public:
    %extend
    {
      PyObject *next() throw(INTERP_KERNEL::Exception)
      {
        MEDFileField1TS *ret=self->nextt();
        if(ret)
          return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__MEDFileField1TS,SWIG_POINTER_OWN | 0);
        else
          {
            PyErr_SetString(PyExc_StopIteration,"No more data.");
            return 0;
          }
      }
    }
  };

  class MEDFileFieldMultiTS : public MEDFileFieldMultiTSWithoutDAS, public MEDFileFieldGlobsReal, public MEDFileWritable
  {
  public:
    static MEDFileFieldMultiTS *New();
    static MEDFileFieldMultiTS *New(const char *fileName);
    static MEDFileFieldMultiTS *New(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception);
    //
    MEDFileField1TS *getTimeStepAtPos(int pos) const throw(INTERP_KERNEL::Exception);
    MEDFileField1TS *getTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    MEDFileField1TS *getTimeStepGivenTime(double time, double eps=1e-8) const throw(INTERP_KERNEL::Exception);
    //
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, int iteration, int order, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtLevelOld(TypeOfField type, const char *mname, int iteration, int order, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    //
    void appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception);
    void appendFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception);
    %extend
       {
         MEDFileFieldMultiTS()
         {
           return MEDFileFieldMultiTS::New();
         }

         MEDFileFieldMultiTS(const char *fileName) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileFieldMultiTS::New(fileName);
         }

         MEDFileFieldMultiTS(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileFieldMultiTS::New(fileName,fieldName);
         }
         
         std::string __str__() const
           {
             return self->simpleRepr();
           }

         MEDFileField1TS *__getitem__(PyObject *elt0) const throw(INTERP_KERNEL::Exception)
         {
           if(elt0 && PyInt_Check(elt0))
             {//fmts[3]
               int pos=PyInt_AS_LONG(elt0);
               return self->getTimeStepAtPos(pos);
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
                       return self->getTimeStep(iter,order);
                     }
                   else
                     throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::__getitem__ : invalid input param ! input is a tuple of size 2 but two integers are expected in this tuple to request a time steps !");
                 }
               else
                 throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::__getitem__ : invalid input param ! input is a tuple of size != 2 ! two integers are expected in this tuple to request a time steps !");
             }
           else if(elt0 && PyFloat_Check(elt0))
             {
               double val=PyFloat_AS_DOUBLE(elt0);
               return self->getTimeStepGivenTime(val);
             }
           else
             throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::__getitem__ : invalid input params ! expected fmts[int], fmts[int,int] or fmts[double] to request time step !");
         }

         MEDFileFieldMultiTSIterator *__iter__() throw(INTERP_KERNEL::Exception)
         {
           return self->iterator();
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
       }
  };

  class MEDFileFieldsIterator
  {
  public:
    %extend
    {
      PyObject *next() throw(INTERP_KERNEL::Exception)
      {
        MEDFileFieldMultiTS *ret=self->nextt();
        if(ret)
          return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_ParaMEDMEM__MEDFileFieldMultiTS,SWIG_POINTER_OWN | 0);
        else
          {
            PyErr_SetString(PyExc_StopIteration,"No more data.");
            return 0;
          }
      }
    }
  };

  class MEDFileFields : public RefCountObject, public MEDFileFieldGlobsReal, public MEDFileWritable
  {
  public:
    static MEDFileFields *New();
    static MEDFileFields *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    int getNumberOfFields() const;
    std::vector<std::string> getFieldsNames() const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getMeshesNames() const throw(INTERP_KERNEL::Exception);
    //
    void resize(int newSize) throw(INTERP_KERNEL::Exception);
    void pushField(MEDFileFieldMultiTS *field) throw(INTERP_KERNEL::Exception);
    void setFieldAtPos(int i, MEDFileFieldMultiTS *field) throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTS *getFieldAtPos(int i) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTS *getFieldWithName(const char *fieldName) const throw(INTERP_KERNEL::Exception);
    void destroyFieldAtPos(int i) throw(INTERP_KERNEL::Exception);
    %extend
       {
         MEDFileFields()
         {
           return MEDFileFields::New();
         }

         MEDFileFields(const char *fileName) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileFields::New(fileName);
         }
         
         std::string __str__() const
         {
           return self->simpleRepr();
         }
         
         MEDFileFieldMultiTS *__getitem__(PyObject *obj) throw(INTERP_KERNEL::Exception)
         {
           if(PyInt_Check(obj))
             {
               return self->getFieldAtPos((int)PyInt_AS_LONG(obj));
             }
           else if(PyString_Check(obj))
             {
               return self->getFieldWithName(PyString_AsString(obj));
             }
           else
             throw INTERP_KERNEL::Exception("MEDFileFields::__getitem__ : only integer or string with fieldname supported !");
         }

         MEDFileFields *__setitem__(int obj, MEDFileFieldMultiTS *field) throw(INTERP_KERNEL::Exception)
         {
           self->setFieldAtPos(obj,field);
           return self;
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
       }
  };

  class MEDFileData : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileData *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileData *New();
    void setFields(MEDFileFields *fields) throw(INTERP_KERNEL::Exception);
    void setMeshes(MEDFileMeshes *meshes) throw(INTERP_KERNEL::Exception);
    int getNumberOfFields() const throw(INTERP_KERNEL::Exception);
    int getNumberOfMeshes() const throw(INTERP_KERNEL::Exception);
    //
    bool changeMeshName(const char *oldMeshName, const char *newMeshName) throw(INTERP_KERNEL::Exception);
    //
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    %extend
       {
         MEDFileData(const char *fileName) throw(INTERP_KERNEL::Exception)
         {
           return MEDFileData::New(fileName);
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

         MEDFileFields *getFields() const
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
    static SauvReader* New(const char *fileName) throw(INTERP_KERNEL::Exception);
    MEDFileData * loadInMEDFileDS() throw(INTERP_KERNEL::Exception);
  };

  class SauvWriter : public RefCountObject
  {
  public:
    static SauvWriter * New();
    void setMEDFileDS(const MEDFileData* medData, unsigned meshIndex = 0);
    void write(const char* fileName);
  };

}
