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
  static void CheckFileForRead(const std::string& fileName);
  static std::vector<std::string> GetMeshNames(const std::string& fileName);
  static std::vector< std::vector< std::pair<INTERP_KERNEL::NormalizedCellType,int> > > GetUMeshGlobalInfo(const std::string& fileName, const std::string& meshName, int &meshDim, int& spaceDim, int& numberOfNodes);
  static std::vector< std::pair<std::string,std::string> > GetComponentsNamesOfField(const std::string& fileName, const std::string& fieldName);
  static std::vector<std::string> GetMeshNamesOnField(const std::string& fileName, const std::string& fieldName);
  static std::vector<std::string> GetMeshGroupsNames(const std::string& fileName, const std::string& meshName);
  static std::vector<std::string> GetMeshFamiliesNames(const std::string& fileName, const std::string& meshName);
  static std::vector<std::string> GetMeshFamiliesNamesOnGroup(const std::string& fileName, const std::string& meshName, const std::string& grpName);
  static std::vector<std::string> GetMeshGroupsNamesOnFamily(const std::string& fileName, const std::string& meshName, const std::string& famName);
  static std::vector<std::string> GetAllFieldNames(const std::string& fileName);
  static std::vector<std::string> GetAllFieldNamesOnMesh(const std::string& fileName, const std::string& meshName);
  static std::vector<ParaMEDMEM::TypeOfField> GetTypesOfField(const std::string& fileName, const std::string& meshName, const std::string& fieldName);
  static std::vector<std::string> GetFieldNamesOnMesh(ParaMEDMEM::TypeOfField type, const std::string& fileName, const std::string& meshName);
  static std::vector<std::string> GetCellFieldNamesOnMesh(const std::string& fileName, const std::string& meshName);
  static std::vector<std::string> GetNodeFieldNamesOnMesh(const std::string& fileName, const std::string& meshName);
  static std::vector< std::pair<int,int> > GetFieldIterations(ParaMEDMEM::TypeOfField type, const std::string& fileName, const std::string& meshName, const std::string& fieldName);
  static std::vector< std::pair<int,int> > GetCellFieldIterations(const std::string& fileName, const std::string& meshName, const std::string& fieldName);
  static std::vector< std::pair<int,int> > GetNodeFieldIterations(const std::string& fileName, const std::string& meshName, const std::string& fieldName);
  static std::vector< std::pair< std::pair<int,int>, double> > GetAllFieldIterations(const std::string& fileName, const std::string& fieldName);
  static double GetTimeAttachedOnFieldIteration(const std::string& fileName, const std::string& fieldName, int iteration, int order);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFamilies(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::vector<std::string>& fams);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromGroups(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::vector<std::string>& grps);
  static ParaMEDMEM::MEDCouplingMesh *ReadMeshFromFile(const std::string& fileName, const std::string& meshName, int meshDimRelToMax=0);
  static ParaMEDMEM::MEDCouplingMesh *ReadMeshFromFile(const std::string& fileName, int meshDimRelToMax=0);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const std::string& fileName, const std::string& meshName, int meshDimRelToMax=0);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const std::string& fileName, int meshDimRelToMax=0);
  static int ReadUMeshDimFromFile(const std::string& fileName, const std::string& meshName);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadField(ParaMEDMEM::TypeOfField type, const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order);
  static std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> ReadFieldsOnSameMesh(ParaMEDMEM::TypeOfField type, const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName,
      const std::vector<std::pair<int,int> >& its);
  static std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> ReadFieldsCellOnSameMesh(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName,
      const std::vector<std::pair<int,int> >& its);
  static std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> ReadFieldsNodeOnSameMesh(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName,
      const std::vector<std::pair<int,int> >& its);
  static std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> ReadFieldsGaussOnSameMesh(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName,
      const std::vector<std::pair<int,int> >& its);
  static std::vector<ParaMEDMEM::MEDCouplingFieldDouble *> ReadFieldsGaussNEOnSameMesh(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName,
      const std::vector<std::pair<int,int> >& its);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldCell(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldNode(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldGauss(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldGaussNE(const std::string& fileName, const std::string& meshName, int meshDimRelToMax, const std::string& fieldName, int iteration, int order);
  static void WriteMesh(const std::string& fileName, const ParaMEDMEM::MEDCouplingMesh *mesh, bool writeFromScratch);
  static void WriteUMesh(const std::string& fileName, const ParaMEDMEM::MEDCouplingUMesh *mesh, bool writeFromScratch);
  static void WriteUMeshDep(const std::string& fileName, const ParaMEDMEM::MEDCouplingUMesh *mesh, bool writeFromScratch);
  static void WriteUMeshesPartition(const std::string& fileName, const std::string& meshName, const std::vector<const ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool writeFromScratch);
  static void WriteUMeshesPartitionDep(const std::string& fileName, const std::string& meshName, const std::vector<const ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool writeFromScratch);
  static void WriteUMeshes(const std::string& fileName, const std::vector<const ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool writeFromScratch);
  static void WriteField(const std::string& fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch);
  static void WriteFieldDep(const std::string& fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch);
  static void WriteFieldUsingAlreadyWrittenMesh(const std::string& fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f);
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
