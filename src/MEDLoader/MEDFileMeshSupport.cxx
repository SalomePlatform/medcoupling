// Copyright (C) 2007-2017  CEA/DEN, EDF R&D
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

#include "MEDFileMeshSupport.hxx"

#include "MEDLoaderBase.hxx"
#include "MEDFileSafeCaller.txx"

#include "InterpKernelAutoPtr.hxx"

using namespace MEDCoupling;

MEDFileMeshSupport *MEDFileMeshSupport::New()
{
  return new MEDFileMeshSupport;
}

MEDFileMeshSupport *MEDFileMeshSupport::New(const std::string& fileName, int smid)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return New(fid,smid);
}

MEDFileMeshSupport *MEDFileMeshSupport::New(med_idt fid, int smid)
{
  return new MEDFileMeshSupport(fid,smid);
}

std::vector<const BigMemoryObject *> MEDFileMeshSupport::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

std::size_t MEDFileMeshSupport::getHeapMemorySizeWithoutChildren() const
{
  return 0;
}

int MEDFileMeshSupport::getSpaceDim() const
{
  return _space_dim;
}

void MEDFileMeshSupport::setSpaceDim(int dim)
{
  _space_dim=dim;
}

int MEDFileMeshSupport::getMeshDim() const
{
  return _mesh_dim;
}

void MEDFileMeshSupport::setMeshDim(int dim)
{
  _mesh_dim=dim;
}

void MEDFileMeshSupport::writeLL(med_idt fid) const
{
}

MEDFileMeshSupport::MEDFileMeshSupport(med_idt fid, int smid)
{
  INTERP_KERNEL::AutoPtr<char> msn(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> description(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  med_axis_type axType;
  int nAxis(MEDsupportMeshnAxis(fid,smid));
  INTERP_KERNEL::AutoPtr<char> axisName(new char[MED_SNAME_SIZE*nAxis+1]),axisUnit(new char[MED_SNAME_SIZE*nAxis+1]);
  MEDFILESAFECALLERRD0(MEDsupportMeshInfo,(fid,smid,msn,&_space_dim,&_mesh_dim,description,&axType,axisName,axisUnit));
}

MEDFileMeshSupport::MEDFileMeshSupport():_space_dim(-1),_mesh_dim(-1)
{
}
