// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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

#include "MEDFileData.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileSafeCaller.txx"
#include "MEDFileBlowStrEltUp.hxx"

#include "InterpKernelAutoPtr.hxx"

using namespace MEDCoupling;

MEDFileData *MEDFileData::New(const std::string& fileName)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return New(fid);
}

MEDFileData *MEDFileData::New(med_idt fid)
{
  return new MEDFileData(fid);
}

MEDFileData *MEDFileData::New()
{
  return new MEDFileData;
}

MEDFileData *MEDFileData::deepCopy() const
{
  MCAuto<MEDFileFields> fields;
  if(_fields.isNotNull())
    fields=_fields->deepCopy();
  MCAuto<MEDFileMeshes> meshes;
  if(_meshes.isNotNull())
    meshes=_meshes->deepCopy();
  MCAuto<MEDFileParameters> params;
  if(_params.isNotNull())
    params=_params->deepCopy();
  MCAuto<MEDFileData> ret(MEDFileData::New());
  ret->_fields=fields; ret->_meshes=meshes; ret->_params=params;
  return ret.retn();
}

std::size_t MEDFileData::getHeapMemorySizeWithoutChildren() const
{
  return _header.capacity();
}

std::vector<const BigMemoryObject *> MEDFileData::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  ret.push_back((const MEDFileFields *)_fields);
  ret.push_back((const MEDFileMeshes *)_meshes);
  ret.push_back((const MEDFileParameters *)_params);
  ret.push_back((const MEDFileMeshSupports *)_mesh_supports);
  ret.push_back((const MEDFileStructureElements *)_struct_elems);
  return ret;

}

/** Return a borrowed reference (caller is not responsible for object destruction) */
MEDFileFields *MEDFileData::getFields() const
{
  return const_cast<MEDFileFields *>(static_cast<const MEDFileFields *>(_fields));
}

/** Return a borrowed reference (caller is not responsible for object destruction) */
MEDFileMeshes *MEDFileData::getMeshes() const
{
  return const_cast<MEDFileMeshes *>(static_cast<const MEDFileMeshes *>(_meshes));
}

/** Return a borrowed reference (caller is not responsible for object destruction) */
MEDFileParameters *MEDFileData::getParams() const
{
  return const_cast<MEDFileParameters *>(static_cast<const MEDFileParameters *>(_params));
}

void MEDFileData::setFields(MEDFileFields *fields)
{
  if(fields)
    fields->incrRef();
  _fields=fields;
}

void MEDFileData::setMeshes(MEDFileMeshes *meshes)
{
  if(meshes)
    meshes->incrRef();
  _meshes=meshes;
}

void MEDFileData::setParams(MEDFileParameters *params)
{
  if(params)
    params->incrRef();
  _params=params;
}

int MEDFileData::getNumberOfFields() const
{
  const MEDFileFields *f=_fields;
  if(!f)
    throw INTERP_KERNEL::Exception("MEDFileData::getNumberOfFields : no fields set !");
  return f->getNumberOfFields();
}

int MEDFileData::getNumberOfMeshes() const
{
  const MEDFileMeshes *m=_meshes;
  if(!m)
    throw INTERP_KERNEL::Exception("MEDFileData::getNumberOfMeshes : no meshes set !");
  return m->getNumberOfMeshes();
}

int MEDFileData::getNumberOfParams() const
{
  const MEDFileParameters *p=_params;
  if(!p)
    throw INTERP_KERNEL::Exception("MEDFileData::getNumberOfParams : no params set !");
  return p->getNumberOfParams();
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
    oss << "No meshes set !!!\n\n";
  oss << "Params part :\n*************\n\n";
  const MEDFileParameters *tmp3=_params;
  if(tmp3)
    {
      tmp3->simpleReprWithoutHeader(oss);
    }
  else
    oss << "No params set !!!\n";
  return oss.str();
}

bool MEDFileData::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab)
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

bool MEDFileData::changeMeshName(const std::string& oldMeshName, const std::string& newMeshName)
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
bool MEDFileData::unPolyzeMeshes()
{
  MEDFileMeshes *ms=_meshes;
  if(!ms)
    return false;
  std::vector< MEDFileMesh * > meshesImpacted;
  std::vector< DataArrayInt * > renumParamsOfMeshImpacted;//same size as meshesImpacted
  std::vector< std::vector<int> > oldCodeOfMeshImpacted,newCodeOfMeshImpacted;//same size as meshesImpacted
  std::vector<MCAuto<DataArrayInt> > memSaverIfThrow;//same size as meshesImpacted
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

void MEDFileData::dealWithStructureElements()
{
  if(_struct_elems.isNull())
    throw INTERP_KERNEL::Exception("MEDFileData::dealWithStructureElements : no structure elements in this !");
  if(_meshes.isNull() || _fields.isNull())
    throw INTERP_KERNEL::Exception("MEDFileData::dealWithStructureElements : meshes and fields must be not null !");
  MEDFileBlowStrEltUp::DealWithSE(_fields,_meshes,_struct_elems);
}

/*!
 * Precondition : all instances in \a mfds should have a single mesh with fields on it. If there is an instance with not exactly one mesh an exception will be thrown.
 * You can invoke MEDFileFields::partOfThisLyingOnSpecifiedMeshName method to make it work.
 */
MCAuto<MEDFileData> MEDFileData::Aggregate(const std::vector<const MEDFileData *>& mfds)
{
  if(mfds.empty())
    throw INTERP_KERNEL::Exception("MEDFileData::Aggregate : empty vector !");
  std::size_t sz(mfds.size()),i(0);
  MCAuto<MEDFileData> ret(MEDFileData::New());
  std::vector<const MEDFileUMesh *> ms(sz);
  std::vector< std::vector< std::pair<int,int> > > dts(sz);
  for(std::vector<const MEDFileData *>::const_iterator it=mfds.begin();it!=mfds.end();it++,i++)
    {
      const MEDFileData *elt(*it);
      if(!elt)
        throw INTERP_KERNEL::Exception("MEDFileData::Aggregate : presence of NULL pointer !");
      const MEDFileMeshes *meshes(elt->getMeshes());
      if(!meshes)
        throw INTERP_KERNEL::Exception("MEDFileData::Aggregate : presence of an instance with no meshes attached on it !");
      if(meshes->getNumberOfMeshes()!=1)
        throw INTERP_KERNEL::Exception("MEDFileData::Aggregate : all instances in input vector must lie on exactly one mesh ! To have it you can invoke partOfThisLyingOnSpecifiedMeshName method.");
      const MEDFileMesh *mesh(meshes->getMeshAtPos(0));
      if(!mesh)
        throw INTERP_KERNEL::Exception("MEDFileData::Aggregate : presence of null mesh in a MEDFileData instance among input vector !");
      const MEDFileUMesh *umesh(dynamic_cast<const MEDFileUMesh *>(mesh));
      if(!umesh)
        throw INTERP_KERNEL::Exception("MEDFileData::Aggregate : works only for unstructured meshes !");
      ms[i]=umesh;
      dts[i]=umesh->getAllDistributionOfTypes();
    }
  MCAuto<MEDFileUMesh> agg_m(MEDFileUMesh::Aggregate(ms));
  MCAuto<MEDFileMeshes> mss(MEDFileMeshes::New()); mss->pushMesh(agg_m);
  ret->setMeshes(mss);
  // fields
  std::vector<std::string> fieldNames(mfds[0]->getFields()->getFieldsNames());
  std::set<std::string> fieldNamess(fieldNames.begin(),fieldNames.end());
  if(fieldNames.size()!=fieldNamess.size())
    throw INTERP_KERNEL::Exception("MEDFileData::Aggregate : field names must be different each other !");
  std::vector< std::vector<const MEDFileAnyTypeFieldMultiTS *> > vectOfFields(fieldNames.size());
  std::vector< std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > > vectOfFields2(fieldNames.size());
  MCAuto<MEDFileFields> fss(MEDFileFields::New());
  for(std::vector<const MEDFileData *>::const_iterator it=mfds.begin();it!=mfds.end();it++)
    {
      std::vector<std::string> fieldNames0((*it)->getFields()->getFieldsNames());
      std::set<std::string> fieldNamess0(fieldNames0.begin(),fieldNames0.end());
      if(fieldNamess!=fieldNamess0)
        throw INTERP_KERNEL::Exception("MEDFileData::Aggregate : field names must be the same for all input instances !");
      i=0;
      for(std::vector<std::string>::const_iterator it1=fieldNames.begin();it1!=fieldNames.end();it1++,i++)
        {
          MCAuto<MEDFileAnyTypeFieldMultiTS> fmts((*it)->getFields()->getFieldWithName(*it1));
          if(fmts.isNull())
            throw INTERP_KERNEL::Exception("MEDFileData::Aggregate : internal error 1 !");
          vectOfFields2[i].push_back(fmts); vectOfFields[i].push_back(fmts);
        }
    }
  i=0;
  for(std::vector<std::string>::const_iterator it1=fieldNames.begin();it1!=fieldNames.end();it1++,i++)
    {
      MCAuto<MEDFileAnyTypeFieldMultiTS> fmts(MEDFileAnyTypeFieldMultiTS::Aggregate(vectOfFields[i],dts));
      fmts->setMeshName(agg_m->getName());
      fss->pushField(fmts);
    }
  ret->setFields(fss);
  return ret;
}

MEDFileData::MEDFileData()
{
}

MEDFileData::MEDFileData(med_idt fid)
try
{
  readHeader(fid);
  _mesh_supports=MEDFileMeshSupports::New(fid);
  _struct_elems=MEDFileStructureElements::New(fid,_mesh_supports);
  _fields=MEDFileFields::NewWithDynGT(fid,_struct_elems,true);
  _meshes=MEDFileMeshes::New(fid);
  _params=MEDFileParameters::New(fid);
}
catch(INTERP_KERNEL::Exception& e)
{
    throw e;
}

void MEDFileData::writeLL(med_idt fid) const
{
  writeHeader(fid);
  if(_meshes.isNotNull())
    _meshes->writeLL(fid);
  if(_fields.isNotNull())
    _fields->writeLL(fid);
  if(_params.isNotNull())
    _params->writeLL(fid);
  if(_mesh_supports.isNotNull())
    _mesh_supports->writeLL(fid);
  if(_struct_elems.isNotNull())
    _struct_elems->writeLL(fid);
}

std::string MEDFileData::getHeader() const
{
  return _header;
}


void MEDFileData::setHeader(const std::string& header)
{
  _header=header;
}

void MEDFileData::readHeader(med_idt fid)
{
  INTERP_KERNEL::AutoPtr<char> header(MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE));
  int ret(MEDfileCommentRd(fid,header));
  if(ret==0)
    _header=MEDLoaderBase::buildStringFromFortran(header,MED_COMMENT_SIZE);
}

void MEDFileData::writeHeader(med_idt fid) const
{
  INTERP_KERNEL::AutoPtr<char> header(MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE));
  MEDLoaderBase::safeStrCpy(_header.c_str(),MED_COMMENT_SIZE,header,_too_long_str);
  MEDFILESAFECALLERWR0(MEDfileCommentWr,(fid,header));
}
