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
      tmp->simpleRepr(0,oss);
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

bool MEDFileData::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception)
{
  bool ret=false;
  MEDFileFields *fields=_fields;
  if(fields)
    ret=fields->changeMeshNames(modifTab) || ret;
  MEDFileMeshes *meshes=_meshes;
  if(meshes)
    ret=meshes->changeNames(modifTab) || ret;
  return ret;
}

bool MEDFileData::changeMeshName(const char *oldMeshName, const char *newMeshName) throw(INTERP_KERNEL::Exception)
{
  std::string oldName(oldMeshName);
  std::vector< std::pair<std::string,std::string> > v(1);
  v[0].first=oldName; v[0].second=newMeshName;
  return changeMeshNames(v);
}

/*!
 * This method performs unpolyzation in meshes in \a this and if it leads to a modification to one or more than one meshes in \a this 
 * the fields are automatically renumbered (for those impacted, that is to say here fields on cells and fields on gauss points on impacted fields)
 *
 * \return If true is returned it means that some meshes in \a this has been modified and eventually fields have been renumbered.
 *         \n If false \a this remains unchanged.
 */
bool MEDFileData::unPolyzeMeshes() throw(INTERP_KERNEL::Exception)
{
  MEDFileMeshes *ms=_meshes;
  if(!ms)
    return false;
  std::vector< MEDFileMesh * > meshesImpacted;
  std::vector< DataArrayInt * > renumParamsOfMeshImpacted;//same size as meshesImpacted
  std::vector< std::vector<int> > oldCodeOfMeshImpacted,newCodeOfMeshImpacted;//same size as meshesImpacted
  std::vector<MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > memSaverIfThrow;//same size as meshesImpacted
  for(int i=0;i<ms->getNumberOfMeshes();i++)
    {
      MEDFileMesh *m=ms->getMeshAtPos(i);
      if(m)
        {
          std::vector<int> oldCode,newCode;
          DataArrayInt *o2nRenumCell=0;
          bool modif=m->unPolyze(oldCode,newCode,o2nRenumCell);
          if(!modif)
            continue;
          renumParamsOfMeshImpacted.push_back(o2nRenumCell); memSaverIfThrow.push_back(o2nRenumCell);
          oldCodeOfMeshImpacted.push_back(oldCode);
          newCodeOfMeshImpacted.push_back(newCode);
          meshesImpacted.push_back(m);
        }
    }
  if(!meshesImpacted.empty())
    {
      MEDFileFields *fs=_fields;
      if(fs)
        for(std::size_t i=0;i<meshesImpacted.size();i++)
          fs->renumberEntitiesLyingOnMesh(meshesImpacted[i]->getName(),oldCodeOfMeshImpacted[i],newCodeOfMeshImpacted[i],renumParamsOfMeshImpacted[i]);
    }
  return !meshesImpacted.empty();
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
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,medmod);
  const MEDFileMeshes *ms=_meshes;
  if(ms)
    ms->write(fid);
  const MEDFileFields *fs=_fields;
  if(fs)
    fs->writeLL(fid);
}
