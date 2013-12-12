// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (CEA/DEN)

#ifndef __MEDLOADER_HXX__
#define __MEDLOADER_HXX__

#include "MEDLoaderDefines.hxx"
#include "InterpKernelException.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "NormalizedUnstructuredMesh.hxx"

#include <list>
#include <vector>

namespace ParaMEDMEM
{
  class DataArrayInt;
  class MEDCouplingMesh;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class MEDFileWritable;
}

class MEDLOADER_EXPORT MEDLoader
{
 public:
  static void SetEpsilonForNodeComp(double val);
  static void SetCompPolicyForCell(int val);
  static void SetTooLongStrPolicy(int val);
  static bool HasXDR();
  static std::string MEDFileVersionStr();
  static void MEDFileVersion(int& major, int& minor, int& release);
  static void CheckFileForRead(const char *fileName);
  static std::vector<std::string> GetMeshNames(const char *fileName);
  static std::vector< std::vector< std::pair<INTERP_KERNEL::NormalizedCellType,int> > > GetUMeshGlobalInfo(const char *fileName, const char *meshName, int &meshDim, int& spaceDim, int& numberOfNodes);
  static std::vector< std::pair<std::string,std::string> > GetComponentsNamesOfField(const char *fileName, const char *fieldName);
  static std::vector<std::string> GetMeshNamesOnField(const char *fileName, const char *fieldName);
  static std::vector<std::string> GetMeshGroupsNames(const char *fileName, const char *meshName);
  static std::vector<std::string> GetMeshFamiliesNames(const char *fileName, const char *meshName);
  static std::vector<std::string> GetMeshFamiliesNamesOnGroup(const char *fileName, const char *meshName, const char *grpName);
  static std::vector<std::string> GetMeshGroupsNamesOnFamily(const char *fileName, const char *meshName, const char *famName);
  static std::vector<std::string> GetAllFieldNames(const char *fileName);
  static std::vector<std::string> GetAllFieldNamesOnMesh(const char *fileName, const char *meshName);
  static std::vector<ParaMEDMEM::TypeOfField> GetTypesOfField(const char *fileName, const char *meshName, const char *fieldName);
  static std::vector<std::string> GetFieldNamesOnMesh(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName);
  static std::vector<std::string> GetCellFieldNamesOnMesh(const char *fileName, const char *meshName);
  static std::vector<std::string> GetNodeFieldNamesOnMesh(const char *fileName, const char *meshName);
  static std::vector< std::pair<int,int> > GetFieldIterations(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, const char *fieldName);
  static std::vector< std::pair<int,int> > GetCellFieldIterations(const char *fileName, const char *meshName, const char *fieldName);
  static std::vector< std::pair<int,int> > GetNodeFieldIterations(const char *fileName, const char *meshName, const char *fieldName);
  static std::vector< std::pair< std::pair<int,int>, double> > GetAllFieldIterations(const char *fileName, const char *fieldName);
  static double GetTimeAttachedOnFieldIteration(const char *fileName, const char *fieldName, int iteration, int order);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFamilies(const char *fileName, const char *meshName, int meshDimRelToMax, const std::vector<std::string>& fams);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromGroups(const char *fileName, const char *meshName, int meshDimRelToMax, const std::vector<std::string>& grps);
  static ParaMEDMEM::MEDCouplingMesh *ReadMeshFromFile(const char *fileName, const char *meshName, int meshDimRelToMax=0);
  static ParaMEDMEM::MEDCouplingMesh *ReadMeshFromFile(const char *fileName, int meshDimRelToMax=0);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const char *fileName, const char *meshName, int meshDimRelToMax=0);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const char *fileName, int meshDimRelToMax=0);
  static int ReadUMeshDimFromFile(const char *fileName, const char *meshName);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadField(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order);
  static std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> ReadFieldsOnSameMesh(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName,
                                                                                const std::vector<std::pair<int,int> >& its) throw(INTERP_KERNEL::Exception);
  static std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> ReadFieldsCellOnSameMesh(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName,
                                                                                    const std::vector<std::pair<int,int> >& its) throw(INTERP_KERNEL::Exception);
  static std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> ReadFieldsNodeOnSameMesh(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName,
                                                                                    const std::vector<std::pair<int,int> >& its) throw(INTERP_KERNEL::Exception);
  static std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> ReadFieldsGaussOnSameMesh(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName,
                                                                                     const std::vector<std::pair<int,int> >& its) throw(INTERP_KERNEL::Exception);
  static std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> ReadFieldsGaussNEOnSameMesh(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName,
                                                                                       const std::vector<std::pair<int,int> >& its) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldCell(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldNode(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldGauss(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldGaussNE(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order);
  static void WriteMesh(const char *fileName, const ParaMEDMEM::MEDCouplingMesh *mesh, bool writeFromScratch);
  static void WriteUMesh(const char *fileName, const ParaMEDMEM::MEDCouplingUMesh *mesh, bool writeFromScratch);
  static void WriteUMeshDep(const char *fileName, const ParaMEDMEM::MEDCouplingUMesh *mesh, bool writeFromScratch);
  static void WriteUMeshesPartition(const char *fileName, const char *meshName, const std::vector<const ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool writeFromScratch);
  static void WriteUMeshesPartitionDep(const char *fileName, const char *meshName, const std::vector<const ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool writeFromScratch);
  static void WriteUMeshes(const char *fileName, const std::vector<const ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool writeFromScratch);
  static void WriteField(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch);
  static void WriteFieldDep(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch);
  static void WriteFieldUsingAlreadyWrittenMesh(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f);
 public:
  static void AssignStaticWritePropertiesTo(ParaMEDMEM::MEDFileWritable& obj);
 private:
  MEDLoader();
 public:
  static double _EPS_FOR_NODE_COMP;
  static int _COMP_FOR_CELL;
  static int _TOO_LONG_STR;
};

#endif
