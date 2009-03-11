//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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

#include "InterpKernelException.hxx"
#include "NormalizedUnstructuredMesh.hxx"

#include <list>
#include <vector>

namespace ParaMEDMEM
{
  class ParaMESH;
  class ParaFIELD;
  class DataArrayInt;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
}

namespace MEDLoader
{
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
    MEDFieldDoublePerCellType(INTERP_KERNEL::NormalizedCellType type, double *values, int ncomp, int nval);
    INTERP_KERNEL::NormalizedCellType getType() const { return _type; }
    int getNbComp() const { return _ncomp; }
    int getNbOfTuple() const { return _nval; }
    int getNbOfValues() const { return _ncomp*_nval; }
    double *getArray() const { return _values; }
    void releaseArray();
  private:
    int _nval;
    int _ncomp;
    double *_values;
    INTERP_KERNEL::NormalizedCellType _type;
  };
  //
  std::vector<std::string> GetMeshNames(const char *fileName);
  std::vector<std::string> GetMeshGroupsNames(const char *fileName, const char *meshName);
  std::vector<std::string> GetMeshFamilyNames(const char *fileName, const char *meshName);
  std::vector<std::string> GetCellFieldNamesOnMesh(const char *fileName, const char *meshName);
  std::vector<std::string> GetNodeFieldNamesOnMesh(const char *fileName, const char *meshName);
  std::vector< std::pair<int,int> > GetCellFieldIterations(const char *fileName, const char *fieldName);
  std::vector< std::pair<int,int> > GetNodeFieldIterations(const char *fileName, const char *fieldName);
  ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFamilies(const char *fileName, const char *meshName, int meshDimRelToMax, const std::vector<std::string>& fams);
  ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromGroups(const char *fileName, const char *meshName, int meshDimRelToMax, const std::vector<std::string>& grps);
  ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const char *fileName, const char *meshName=0, int meshDimRelToMax=0) throw(INTERP_KERNEL::Exception);
  ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldDoubleCell(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order);
  ParaMEDMEM::MEDCouplingFieldDouble *ReadFieldDoubleNode(const char *fileName, const char *meshName, int meshDimRelToMax, const char *fieldName, int iteration, int order);
  void writeUMesh(const char *fileName, ParaMEDMEM::MEDCouplingUMesh *mesh);
  void writeParaMesh(const char *fileName, ParaMEDMEM::ParaMESH *mesh);
  void writeParaField(const char *fileName, const char *meshName, ParaMEDMEM::ParaFIELD *f);
}

#endif
