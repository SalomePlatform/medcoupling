// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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

#include "MEDFileData.hxx"

using namespace ParaMEDMEM;

MEDFileData *MEDFileData::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileData(fileName);
}

MEDFileData *MEDFileData::New()
{
  return new MEDFileData;
}

MEDFileFields *MEDFileData::getFields() const
{
  return const_cast<MEDFileFields *>(static_cast<const MEDFileFields *>(_fields));
}

MEDFileMeshes *MEDFileData::getMeshes() const
{
  return const_cast<MEDFileMeshes *>(static_cast<const MEDFileMeshes *>(_meshes));
}

void MEDFileData::setFields(MEDFileFields *fields) throw(INTERP_KERNEL::Exception)
{
  if(!fields)
    throw INTERP_KERNEL::Exception("MEDFileData::setFields : input pointer is null !");
  fields->incrRef();
  _fields=fields;
}

void MEDFileData::setMeshes(MEDFileMeshes *meshes) throw(INTERP_KERNEL::Exception)
{
  if(!meshes)
    throw INTERP_KERNEL::Exception("MEDFileData::setMeshes : input pointer is null !");
  meshes->incrRef();
  _meshes=meshes;
}

int MEDFileData::getNumberOfFields() const throw(INTERP_KERNEL::Exception)
{
  const MEDFileFields *f=_fields;
  if(!f)
    throw INTERP_KERNEL::Exception("MEDFileData::getNumberOfFields : no fields set !");
  return f->getNumberOfFields();
}

int MEDFileData::getNumberOfMeshes() const throw(INTERP_KERNEL::Exception)
{
  const MEDFileMeshes *m=_meshes;
  if(!m)
    throw INTERP_KERNEL::Exception("MEDFileData::getNumberOfMeshes : no meshes set !");
  return m->getNumberOfMeshes();
}

std::string MEDFileData::simpleRepr() const
{
  std::ostringstream oss;
  oss << "(***************)\n(* MEDFileData *)\n(***************)\n\nFields part :\n*************\n\n";
  const MEDFileFields *tmp=_fields;
  if(tmp)
    {
      tmp->simpleReprWithoutHeader(oss);
      oss << std::endl;
    }
  else
    oss << "No fields set !!!\n\n";
  oss << "Meshes part :\n*************\n\n";
  const MEDFileMeshes *tmp2=_meshes;
  if(tmp2)
    {
      tmp2->simpleReprWithoutHeader(oss);
    }
  else
    oss << "No meshes set !!!\n";
  return oss.str();
}

MEDFileData::MEDFileData()
{
}

MEDFileData::MEDFileData(const char *fileName) throw(INTERP_KERNEL::Exception)
try
  {
    _fields=MEDFileFields::New(fileName);
    _meshes=MEDFileMeshes::New(fileName);
  }
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

void MEDFileData::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileMeshes *ms=_meshes;
  if(ms)
    ms->write(fileName,mode);
  int mode2=mode==2?0:mode;
  const MEDFileFields *fs=_fields;
  if(fs)
    fs->write(fileName,mode2);
}
