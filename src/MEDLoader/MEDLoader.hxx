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
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
}

class MEDLOADER_EXPORT MEDLoader
{
 public:
/// @cond INTERNAL
  class MEDConnOfOneElemType
  {
  public:
    MEDConnOfOneElemType(INTERP_KERNEL::NormalizedCellType type, int *conn, int *index, int *fam, int lgth, int connLgth);
    INTERP_KERNEL::NormalizedCellType getType() const { return _type; }
    int getLength() const { return _lgth; }
    int getConnLength() const { return _conn_lgth; }
    int *getArray() const { return _conn; }
    int *getIndex() const { return _index; }
    int *getFam() const { return _fam; }
    void setGlobal(int *global);
    const int *getGlobal() const { return _global; }
    void releaseArray();
  private:
    int _lgth;
    int *_fam;
    int *_conn;
    int *_index;
    int *_global;
    int _conn_lgth;
    INTERP_KERNEL::NormalizedCellType _type;
  };

  class MEDFieldDoublePerCellType
  {
  public:
    MEDFieldDoublePerCellType(INTERP_KERNEL::NormalizedCellType type, double *values, int ncomp, int ntuple, const int *cellIdPerType, const char *locName);
    INTERP_KERNEL::NormalizedCellType getType() const { return _type; }
    int getNbComp() const { return _ncomp; }
    int getNbOfTuple() const { return _ntuple; }
    int getNbOfValues() const { return _ncomp*_ntuple; }
    double *getArray() const { return _values; }
    const std::string& getLocName() const { return _loc_name; }
    const std::vector<int>& getCellIdPerType() const { return _cell_id_per_type; }
    void releaseArray();
  private:
    int _ntuple;
    int _ncomp;
    double *_values;
    std::string _loc_name;
    std::vector<int> _cell_id_per_type;
    INTERP_KERNEL::NormalizedCellType _type;
  };
/// @endcond
  static void setEpsilonForNodeComp(double val) throw(INTERP_KERNEL::Exception);
  static void setCompPolicyForCell(int val) throw(INTERP_KERNEL::Exception);
  static void setTooLongStrPolicy(int val) throw(INTERP_KERNEL::Exception);
  static void CheckFileForRead(const char *fileName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshNames(const char *fileName) throw(INTERP_KERNEL::Exception);
  static std::vector< std::vector< std::pair<INTERP_KERNEL::NormalizedCellType,int> > > GetUMeshGlobalInfo(const char *fileName, const char *meshName, int &meshDim, int& spaceDim, int& numberOfNodes) throw(INTERP_KERNEL::Exception);
  static std::vector< std::pair<std::string,std::string> > GetComponentsNamesOfField(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshNamesOnField(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshGroupsNames(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshFamiliesNames(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshFamiliesNamesOnGroup(const char *fileName, const char *meshName, const char *grpName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetMeshGroupsNamesOnFamily(const char *fileName, const char *meshName, const char *famName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetAllFieldNames(const char *fileName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetAllFieldNamesOnMesh(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception);
  static std::vector<ParaMEDMEM::TypeOfField> GetTypesOfField(const char *fileName, const char *meshName, const char *fieldName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetFieldNamesOnMesh(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetCellFieldNamesOnMesh(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception);
  static std::vector<std::string> GetNodeFieldNamesOnMesh(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception);
  static std::vector< std::pair<int,int> > GetFieldIterations(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, const char *fieldName) throw(INTERP_KERNEL::Exception);
  static std::vector< std::pair<int,int> > GetCellFieldIterations(const char *fileName, const char *meshName, const char *fieldName) throw(INTERP_KERNEL::Exception);
  static std::vector< std::pair<int,int> > GetNodeFieldIterations(const char *fileName, const char *meshName, const char *fieldName) throw(INTERP_KERNEL::Exception);
  static std::vector< std::pair< std::pair<int,int>, double> > GetAllFieldIterations(const char *fileName, const char *meshName, const char *fieldName) throw(INTERP_KERNEL::Exception);
  static double GetTimeAttachedOnFieldIteration(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFamilies(const char *fileName, const char *meshName, int meshDimRelToMax, const std::vector<std::string>& fams) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromGroups(const char *fileName, const char *meshName, int meshDimRelToMax, const std::vector<std::string>& grps) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const char *fileName, const char *meshName, int meshDimRelToMax=0) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const char *fileName, int meshDimRelToMax=0) throw(INTERP_KERNEL::Exception);
  static int ReadUMeshDimFromFile(const char *fileName, const char *meshName) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadField(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
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
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldCell(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldNode(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldGauss(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldGaussNE(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  static void WriteUMesh(const char *fileName, const ParaMEDMEM::MEDCouplingUMesh *mesh, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteUMeshDep(const char *fileName, const ParaMEDMEM::MEDCouplingUMesh *mesh, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteUMeshesPartition(const char *fileName, const char *meshName, const std::vector<const ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteUMeshesPartitionDep(const char *fileName, const char *meshName, const std::vector<const ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteUMeshes(const char *fileName, const std::vector<const ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteField(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteFieldDep(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch) throw(INTERP_KERNEL::Exception);
  static void WriteFieldUsingAlreadyWrittenMesh(const char *fileName, const ParaMEDMEM::MEDCouplingFieldDouble *f) throw(INTERP_KERNEL::Exception);
 private:
  MEDLoader();
 public:
  static double _EPS_FOR_NODE_COMP;
  static int _COMP_FOR_CELL;
  static int _TOO_LONG_STR;
};

#endif
