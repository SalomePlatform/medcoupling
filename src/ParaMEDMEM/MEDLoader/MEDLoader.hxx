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
}

namespace MEDLoader
{
  class MEDConnOfOneElemType
  {
  public:
    MEDConnOfOneElemType(INTERP_KERNEL::NormalizedCellType type, int *conn, int lgth);
    INTERP_KERNEL::NormalizedCellType getType() const { return _type; }
    int getLength() const { return _lgth; }
    int *getArray() const { return _conn; }
    void setGlobal(int *global);
    void releaseArray();
  private:
    int _lgth;
    int *_conn;
    int *_global;
    INTERP_KERNEL::NormalizedCellType _type;
  };
  unsigned calculateHighestMeshDim(const std::list<MEDConnOfOneElemType>& conn);
  void keepSpecifiedMeshDim(std::list<MEDConnOfOneElemType>& conn, unsigned meshDim);
  void tradMEDFileCoreFrmt2MEDCouplingUMesh(const std::list<MEDConnOfOneElemType>& medConnFrmt,
                                            ParaMEDMEM::DataArrayInt* &conn,
                                            ParaMEDMEM::DataArrayInt* &connIndex);
  void releaseMEDFileCoreFrmt(std::list<MEDLoader::MEDConnOfOneElemType>& medConnFrmt);
  void writeMasterFile(const char *fileName, const std::vector<std::string>& fileNames, const char *meshName);
  int buildMEDSubConnectivityOfOneType(ParaMEDMEM::DataArrayInt *conn, ParaMEDMEM::DataArrayInt *connIndex,
                                       INTERP_KERNEL::NormalizedCellType type, std::vector<int>& conn4MEDFile);
  //
  ParaMEDMEM::MEDCouplingUMesh *ReadUMeshFromFile(const char *fileName, const char *meshName=0, int meshDimRelToMax=0) throw(INTERP_KERNEL::Exception);
  void writeUMesh(const char *fileName, ParaMEDMEM::MEDCouplingUMesh *mesh);
  void writeParaMesh(const char *fileName, ParaMEDMEM::ParaMESH *mesh);
  void writeParaField(const char *fileName, const char *meshName, ParaMEDMEM::ParaFIELD *f);
}

#endif
