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
    MEDFieldDoublePerCellType(INTERP_KERNEL::NormalizedCellType type, double *values, int ncomp, int ntuple, const int *cellIdPerType);
    INTERP_KERNEL::NormalizedCellType getType() const { return _type; }
    int getNbComp() const { return _ncomp; }
    int getNbOfTuple() const { return _ntuple; }
    int getNbOfValues() const { return _ncomp*_ntuple; }
    double *getArray() const { return _values; }
    const std::vector<int>& getCellIdPerType() const { return _cell_id_per_type; }
    void releaseArray();
  private:
    int _ntuple;
    int _ncomp;
    double *_values;
    std::vector<int> _cell_id_per_type;
    INTERP_KERNEL::NormalizedCellType _type;
  };
  //
  static void setEpsilonForNodeComp(double val);
  static void setCompPolicyForCell(int val);
  static std::vector<std::string> GetMeshNames(const char *fileName);
  static std::vector<std::string> GetMeshGroupsNames(const char *fileName, const char *meshName);
  static std::vector<std::string> GetMeshFamilyNames(const char *fileName, const char *meshName);
  static std::vector<std::string> GetCellFieldNamesOnMesh(const char *fileName, const char *meshName);
  static std::vector<std::string> GetNodeFieldNamesOnMesh(const char *fileName, const char *meshName);
  static std::vector< std::pair<int,int> > GetFieldIterations(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, const char *fieldName);
  static std::vector< std::pair<int,int> > GetCellFieldIterations(const char *fileName, const char *meshName, const char *fieldName);
  static std::vector< std::pair<int,int> > GetNodeFieldIterations(const char *fileName, const char *meshName, const char *fieldName);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFamilies(const char *fileName, const char *meshName, int meshDimRelToMax, const std::vector<std::string>& fams);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromGroups(const char *fileName, const char *meshName, int meshDimRelToMax, const std::vector<std::string>& grps);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const char *fileName, const char *meshName, int meshDimRelToMax=0) throw(INTERP_KERNEL::Exception);
  static ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const char *fileName, int meshDimRelToMax=0) throw(INTERP_KERNEL::Exception);
  static int ReadUMeshDimFromFile(const char *fileName, const char *meshName);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldDouble(ParaMEDMEM::TypeOfField type, const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldDoubleCell(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order);
  static ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldDoubleNode(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order);
  static void WriteUMesh(const char *fileName, ParaMEDMEM::MEDCouplingUMesh *mesh, bool writeFromScratch);
  static void WriteUMeshes(const char *fileName, const char *meshName, const std::vector<ParaMEDMEM::MEDCouplingUMesh *>& meshes, bool writeFromScratch);
  static void WriteField(const char *fileName, ParaMEDMEM::MEDCouplingFieldDouble *f, bool writeFromScratch);
  static void WriteFieldUsingAlreadyWrittenMesh(const char *fileName, ParaMEDMEM::MEDCouplingFieldDouble *f);
 private:
  MEDLoader();
 public:
  static double _EPS_FOR_NODE_COMP;
  static int _COMP_FOR_CELL;
};

#endif
