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
// Author : Anthony Geay (EDF R&D)

#include "ParaMEDFileMesh.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "MEDFileMesh.hxx"
#include "MEDFileMeshLL.hxx"
#include "MEDLoader.hxx"

using namespace ParaMEDMEM;

MEDFileMesh *ParaMEDFileMesh::New(const MPI_Comm com, const MPI_Info nfo, const std::string& fileName, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  ParaMEDMEM::MEDCouplingMeshType meshType;
  MEDFileUtilities::AutoFid fid(MEDfileOpen(fileName.c_str(),MED_ACC_RDONLY));
  int dummy0,dummy1;
  std::string dummy2;
  MEDFileMeshL2::GetMeshIdFromName(fid,mName,meshType,dummy0,dummy1,dummy2);
  switch(meshType)
  {
    case UNSTRUCTURED:
      {
        return ParaMEDFileUMesh::New(com,nfo,fileName,mName,dt,it,mrs);
      }
    default:
      throw INTERP_KERNEL::Exception("ParaMEDFileMesh::New : only unstructured mesh supported for the moment !");
  }
}

MEDFileUMesh *ParaMEDFileUMesh::New(const MPI_Comm com, const MPI_Info nfo, const std::string& fileName, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileUMesh> ret(MEDFileUMesh::New());
  MEDFileUtilities::AutoFid fid(MEDparFileOpen(fileName.c_str(),MED_ACC_RDONLY,com,nfo));
  //
  int meshDim, spaceDim, numberOfNodes;
  std::vector< std::vector< std::pair<INTERP_KERNEL::NormalizedCellType,int> > > typesDistrib(MEDLoader::GetUMeshGlobalInfo(fileName,mName,meshDim,spaceDim,numberOfNodes));
  return ret.retn();
}

MEDFileMeshes *ParaMEDFileMeshes::New(const MPI_Comm com, const MPI_Info nfo, const std::string& fileName)
{
  std::vector<std::string> ms(MEDLoader::GetMeshNames(fileName));
  MEDCouplingAutoRefCountObjectPtr<MEDFileMeshes> ret(MEDFileMeshes::New());
  for(std::vector<std::string>::const_iterator it=ms.begin();it!=ms.end();it++)
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mesh(ParaMEDFileMesh::New(com,nfo,fileName,(*it)));
      ret->pushMesh(mesh);
    }
  return ret.retn();
}
