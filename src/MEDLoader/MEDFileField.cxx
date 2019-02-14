// Copyright (C) 2007-2019  CEA/DEN, EDF R&D
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

#include "MEDFileField.hxx"
#include "MEDFileMesh.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDLoaderTraits.hxx"
#include "MEDFileSafeCaller.txx"
#include "MEDFileFieldOverView.hxx"
#include "MEDFileBlowStrEltUp.hxx"
#include "MEDFileFieldVisitor.hxx"

#include "MEDCouplingFieldDiscretization.hxx"

#include "InterpKernelAutoPtr.hxx"
#include "CellModel.hxx"

#include <algorithm>
#include <iterator>

extern INTERP_KERNEL::NormalizedCellType typmai2[MED_N_CELL_FIXED_GEO];
extern med_geometry_type typmainoeud[1];
extern med_geometry_type typmai3[34];

using namespace MEDCoupling;

//= MEDFileFields

MEDFileFields *MEDFileFields::New()
{
  return new MEDFileFields;
}

MEDFileFields *MEDFileFields::New(const std::string& fileName, bool loadAll)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return New(fid,loadAll);
}

MEDFileFields *MEDFileFields::NewAdv(const std::string& fileName, bool loadAll, const MEDFileEntities *entities)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return NewAdv(fid,loadAll,entities);
}

MEDFileFields *MEDFileFields::NewAdv(med_idt fid, bool loadAll, const MEDFileEntities *entities)
{
  return new MEDFileFields(fid,loadAll,0,entities);
}

MEDFileFields *MEDFileFields::NewWithDynGT(const std::string& fileName, const MEDFileStructureElements *se, bool loadAll)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return NewWithDynGT(fid,se,loadAll);
}

MEDFileFields *MEDFileFields::NewWithDynGT(med_idt fid, const MEDFileStructureElements *se, bool loadAll)
{
  if(!se)
    throw INTERP_KERNEL::Exception("MEDFileFields::NewWithDynGT : null struct element pointer !");
  INTERP_KERNEL::AutoCppPtr<MEDFileEntities> entities(MEDFileEntities::BuildFrom(*se));
  return new MEDFileFields(fid,loadAll,0,entities);
}

MEDFileFields *MEDFileFields::New(med_idt fid, bool loadAll)
{
  return new MEDFileFields(fid,loadAll,0,0);
}

MEDFileFields *MEDFileFields::LoadPartOf(const std::string& fileName, bool loadAll, const MEDFileMeshes *ms)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return new MEDFileFields(fid,loadAll,ms,0);
}

MEDFileFields *MEDFileFields::LoadSpecificEntities(const std::string& fileName, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> >& entities, bool loadAll)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  INTERP_KERNEL::AutoCppPtr<MEDFileEntities> ent(new MEDFileStaticEntities(entities));
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return new MEDFileFields(fid,loadAll,0,ent);
}

std::size_t MEDFileFields::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(MEDFileFieldGlobsReal::getHeapMemorySizeWithoutChildren());
  ret+=_fields.capacity()*sizeof(MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA>);
  return ret;
}

std::vector<const BigMemoryObject *> MEDFileFields::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    ret.push_back((const MEDFileAnyTypeFieldMultiTSWithoutSDA *)*it);
  return ret;
}

MEDFileFields *MEDFileFields::deepCopy() const
{
  MCAuto<MEDFileFields> ret(shallowCpy());
  std::size_t i(0);
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      if((const MEDFileAnyTypeFieldMultiTSWithoutSDA*)*it)
        ret->_fields[i]=(*it)->deepCopy();
    }
  ret->deepCpyGlobs(*this);
  return ret.retn();
}

MEDFileFields *MEDFileFields::shallowCpy() const
{
  return new MEDFileFields(*this);
}

/*!
 * This method scans for all fields in \a this which time steps ids are common. Time step are discriminated by the pair of integer (iteration,order) whatever
 * the double time value. If all returned time steps are \b exactly those for all fields in \a this output parameter \a areThereSomeForgottenTS will be set to false.
 * If \a areThereSomeForgottenTS is set to true, only the sorted intersection of time steps present for all fields in \a this will be returned.
 *
 * \param [out] areThereSomeForgottenTS - indicates to the caller if there is some time steps in \a this that are not present for all fields in \a this.
 * \return the sorted list of time steps (specified with a pair of integer iteration first and order second) present for all fields in \a this.
 * 
 * \sa MEDFileFields::partOfThisLyingOnSpecifiedTimeSteps, MEDFileFields::partOfThisNotLyingOnSpecifiedTimeSteps
 */
std::vector< std::pair<int,int> > MEDFileFields::getCommonIterations(bool& areThereSomeForgottenTS) const
{
  std::set< std::pair<int,int> > s;
  bool firstShot=true;
  areThereSomeForgottenTS=false;
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      if(!(const MEDFileAnyTypeFieldMultiTSWithoutSDA*)*it)
        continue;
      std::vector< std::pair<int,int> > v=(*it)->getIterations();
      std::set< std::pair<int,int> > s1; std::copy(v.begin(),v.end(),std::inserter(s1,s1.end()));
      if(firstShot)
        { s=s1; firstShot=false; }
      else
        {
          std::set< std::pair<int,int> > s2; std::set_intersection(s.begin(),s.end(),s1.begin(),s1.end(),std::inserter(s2,s2.end()));
          if(s!=s2)
            areThereSomeForgottenTS=true;
          s=s2;
        }
    }
  std::vector< std::pair<int,int> > ret;
  std::copy(s.begin(),s.end(),std::back_insert_iterator< std::vector< std::pair<int,int> > >(ret));
  return ret;
}

int MEDFileFields::getNumberOfFields() const
{
  return _fields.size();
}

std::vector<std::string> MEDFileFields::getFieldsNames() const
{
  std::vector<std::string> ret(_fields.size());
  int i(0);
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *f=(*it);
      if(f)
        {
          ret[i]=f->getName();
        }
      else
        {
          std::ostringstream oss; oss << "MEDFileFields::getFieldsNames : At rank #" << i << " field is not defined !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  return ret;
}

std::vector<std::string> MEDFileFields::getMeshesNames() const
{
  std::vector<std::string> ret;
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *cur(*it);
      if(cur)
        ret.push_back(cur->getMeshName());
    }
  return ret;
}

std::string MEDFileFields::simpleRepr() const
{
  std::ostringstream oss;
  oss << "(*****************)\n(* MEDFileFields *)\n(*****************)\n\n";
  simpleRepr(0,oss);
  return oss.str();
}

void MEDFileFields::simpleRepr(int bkOffset, std::ostream& oss) const
{
  int nbOfFields(getNumberOfFields());
  std::string startLine(bkOffset,' ');
  oss << startLine << "There are " << nbOfFields << " fields in this :" << std::endl;
  int i=0;
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *cur=(*it);
      if(cur)
        {
          oss << startLine << "  - # "<< i << " has the following name : \"" << cur->getName() << "\"." << std::endl;
        }
      else
        {
          oss << startLine << "  - not defined !" << std::endl;
        }
    }
  i=0;
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *cur=(*it);
      std::string chapter(17,'0'+i);
      oss << startLine << chapter << std::endl;
      if(cur)
        {
          cur->simpleRepr(bkOffset+2,oss,i);
        }
      else
        {
          oss << startLine << "  - not defined !" << std::endl;
        }
      oss << startLine << chapter << std::endl;
    }
  simpleReprGlobs(oss);
}

MEDFileFields::MEDFileFields()
{
}

MEDFileFields::MEDFileFields(med_idt fid, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
try:MEDFileFieldGlobsReal(fid)
{
  int nbFields(MEDnField(fid));
  _fields.resize(nbFields);
  med_field_type typcha;
  for(int i=0;i<nbFields;i++)
    {
      std::vector<std::string> infos;
      std::string fieldName,dtunit,meshName;
      int nbOfStep(MEDFileAnyTypeField1TS::LocateField2(fid,i,false,fieldName,typcha,infos,dtunit,meshName));
      switch(typcha)
      {
        case MED_FLOAT64:
          {
            _fields[i]=MEDFileFieldMultiTSWithoutSDA::New(fid,fieldName,meshName,typcha,infos,nbOfStep,dtunit,loadAll,ms,entities);
            break;
          }
        case MED_INT32:
          {
            _fields[i]=MEDFileIntFieldMultiTSWithoutSDA::New(fid,fieldName,meshName,typcha,infos,nbOfStep,dtunit,loadAll,ms,entities);
            break;
          }
        case MED_FLOAT32:
          {
            _fields[i]=MEDFileFloatFieldMultiTSWithoutSDA::New(fid,fieldName,meshName,typcha,infos,nbOfStep,dtunit,loadAll,ms,entities);
            break;
          }
        case MED_INT:
          {
            if(sizeof(med_int)==sizeof(int))
              {
                _fields[i]=MEDFileIntFieldMultiTSWithoutSDA::New(fid,fieldName,meshName,typcha,infos,nbOfStep,dtunit,loadAll,ms,entities);
                break;
              }
          }
        default:
          {
            std::ostringstream oss; oss << "constructor MEDFileFields(fileName) : file \'" << FileNameFromFID(fid) << "\' at pos #" << i << " field has name \'" << fieldName << "\' but the type of field is not in [MED_FLOAT64, MED_INT32, MED_FLOAT32] !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
    }
  loadAllGlobals(fid,entities);
}
catch(INTERP_KERNEL::Exception& e)
{
    throw e;
}

void MEDFileFields::writeLL(med_idt fid) const
{
  int i=0;
  writeGlobals(fid,*this);
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *elt=*it;
      if(!elt)
        {
          std::ostringstream oss; oss << "MEDFileFields::write : at rank #" << i << "/" << _fields.size() << " field is empty !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      elt->writeLL(fid,*this);
    }
}

/*!
 * This method alloc the arrays and load potentially huge arrays contained in this field.
 * This method should be called when a MEDFileAnyTypeFieldMultiTS::New constructor has been with false as the last parameter.
 * This method can be also called to refresh or reinit values from a file.
 * 
 * \throw If the fileName is not set or points to a non readable MED file.
 */
void MEDFileFields::loadArrays()
{
  if(getFileName().empty())
    throw INTERP_KERNEL::Exception("MEDFileFields::loadArrays : the structure does not come from a file !");
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(getFileName()));
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::iterator it=_fields.begin();it!=_fields.end();it++)
    {
      MEDFileAnyTypeFieldMultiTSWithoutSDA *elt(*it);
      if(elt)
        elt->loadBigArraysRecursively(fid,*elt);
    }
}

/*!
 * This method behaves as MEDFileFields::loadArrays does, the first call, if \a this was built using a file without loading big arrays.
 * But once data loaded once, this method does nothing.
 * 
 * \throw If the fileName is not set or points to a non readable MED file.
 * \sa MEDFileFields::loadArrays, MEDFileFields::unloadArrays
 */
void MEDFileFields::loadArraysIfNecessary()
{
  if(!getFileName().empty())
    {
      MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(getFileName()));
      for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::iterator it=_fields.begin();it!=_fields.end();it++)
        {
          MEDFileAnyTypeFieldMultiTSWithoutSDA *elt(*it);
          if(elt)
            elt->loadBigArraysRecursivelyIfNecessary(fid,*elt);
        }
    }
}

/*!
 * This method releases potentially big data arrays and so returns to the same heap memory than status loaded with 'loadAll' parameter set to false.
 * \b WARNING, this method does release arrays even if \a this does not come from a load of a MED file.
 * So this method can lead to a loss of data. If you want to unload arrays safely call MEDFileFields::unloadArraysWithoutDataLoss instead.
 * 
 * \sa MEDFileFields::loadArrays, MEDFileFields::loadArraysIfNecessary, MEDFileFields::unloadArraysWithoutDataLoss
 */
void MEDFileFields::unloadArrays()
{
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::iterator it=_fields.begin();it!=_fields.end();it++)
    {
      MEDFileAnyTypeFieldMultiTSWithoutSDA *elt(*it);
      if(elt)
        elt->unloadArrays();
    }
}

/*!
 * This method potentially releases big data arrays if \a this is coming from a file. If \a this has been built from scratch this method will have no effect.
 * This method is the symmetrical method of MEDFileFields::loadArraysIfNecessary.
 * This method is useful to reduce \b safely amount of heap memory necessary for \a this by using MED file as database.
 * 
 * \sa MEDFileFields::loadArraysIfNecessary
 */
void MEDFileFields::unloadArraysWithoutDataLoss()
{
  if(!getFileName().empty())
    unloadArrays();
}

std::vector<std::string> MEDFileFields::getPflsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MCAuto< MEDFileAnyTypeFieldMultiTSWithoutSDA > >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsed2();
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

std::vector<std::string> MEDFileFields::getLocsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MCAuto< MEDFileAnyTypeFieldMultiTSWithoutSDA > >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      std::vector<std::string> tmp((*it)->getLocsReallyUsed2());
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

std::vector<std::string> MEDFileFields::getPflsReallyUsedMulti() const
{
  std::vector<std::string> ret;
  for(std::vector< MCAuto< MEDFileAnyTypeFieldMultiTSWithoutSDA > >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      std::vector<std::string> tmp((*it)->getPflsReallyUsedMulti2());
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

std::vector<std::string> MEDFileFields::getLocsReallyUsedMulti() const
{
  std::vector<std::string> ret;
  for(std::vector< MCAuto< MEDFileAnyTypeFieldMultiTSWithoutSDA > >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      std::vector<std::string> tmp((*it)->getLocsReallyUsed2());
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

void MEDFileFields::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  for(std::vector< MCAuto< MEDFileAnyTypeFieldMultiTSWithoutSDA > >::iterator it=_fields.begin();it!=_fields.end();it++)
    (*it)->changePflsRefsNamesGen2(mapOfModif);
}

void MEDFileFields::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  for(std::vector< MCAuto< MEDFileAnyTypeFieldMultiTSWithoutSDA > >::iterator it=_fields.begin();it!=_fields.end();it++)
    (*it)->changeLocsRefsNamesGen2(mapOfModif);
}

void MEDFileFields::resize(int newSize)
{
  _fields.resize(newSize);
}

void MEDFileFields::pushFields(const std::vector<MEDFileAnyTypeFieldMultiTS *>& fields)
{
  for(std::vector<MEDFileAnyTypeFieldMultiTS *>::const_iterator it=fields.begin();it!=fields.end();it++)
    pushField(*it);
}

void MEDFileFields::pushField(MEDFileAnyTypeFieldMultiTS *field)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileFields::pushMesh : invalid input pointer ! should be different from 0 !");
  _fields.push_back(field->getContent());
  appendGlobs(*field,1e-12);
}

void MEDFileFields::setFieldAtPos(int i, MEDFileAnyTypeFieldMultiTS *field)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileFields::setFieldAtPos : invalid input pointer ! should be different from 0 !");
  if(i>=(int)_fields.size())
    _fields.resize(i+1);
  _fields[i]=field->getContent();
  appendGlobs(*field,1e-12);
}

void MEDFileFields::destroyFieldAtPos(int i)
{
  destroyFieldsAtPos(&i,&i+1);
}

void MEDFileFields::destroyFieldsAtPos(const int *startIds, const int *endIds)
{
  std::vector<bool> b(_fields.size(),true);
  for(const int *i=startIds;i!=endIds;i++)
    {
      if(*i<0 || *i>=(int)_fields.size())
        {
          std::ostringstream oss; oss << "MEDFileFields::destroyFieldsAtPos : Invalid given id in input (" << *i << ") should be in [0," << _fields.size() << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      b[*i]=false;
    }
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > fields(std::count(b.begin(),b.end(),true));
  std::size_t j=0;
  for(std::size_t i=0;i<_fields.size();i++)
    if(b[i])
      fields[j++]=_fields[i];
  _fields=fields;
}

void MEDFileFields::destroyFieldsAtPos2(int bg, int end, int step)
{
  static const char msg[]="MEDFileFields::destroyFieldsAtPos2";
  int nbOfEntriesToKill(DataArrayInt::GetNumberOfItemGivenBESRelative(bg,end,step,msg));
  std::vector<bool> b(_fields.size(),true);
  int k=bg;
  for(int i=0;i<nbOfEntriesToKill;i++,k+=step)
    {
      if(k<0 || k>=(int)_fields.size())
        {
          std::ostringstream oss; oss << "MEDFileFields::destroyFieldsAtPos2 : Invalid given id in input (" << k << ") should be in [0," << _fields.size() << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      b[k]=false;
    }
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > fields(std::count(b.begin(),b.end(),true));
  std::size_t j(0);
  for(std::size_t i=0;i<_fields.size();i++)
    if(b[i])
      fields[j++]=_fields[i];
  _fields=fields;
}

bool MEDFileFields::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab)
{
  bool ret(false);
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::iterator it=_fields.begin();it!=_fields.end();it++)
    {
      MEDFileAnyTypeFieldMultiTSWithoutSDA *cur(*it);
      if(cur)
        ret=cur->changeMeshNames(modifTab) || ret;
    }
  return ret;
}

/*!
 * \param [in] meshName the name of the mesh that will be renumbered.
 * \param [in] oldCode is of format returned by MEDCouplingUMesh::getDistributionOfTypes. And for each *i* oldCode[3*i+2] gives the position (MEDFileUMesh::PutInThirdComponentOfCodeOffset).
 *             This code corresponds to the distribution of types in the corresponding mesh.
 * \param [in] newCode idem to param \a oldCode except that here the new distribution is given.
 * \param [in] renumO2N the old to new renumber array.
 * \return If true a renumbering has been performed. The structure in \a this has been modified. If false, nothing has been done: it is typically the case if \a meshName is not referred by any 
 *         field in \a this.
 */
bool MEDFileFields::renumberEntitiesLyingOnMesh(const std::string& meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N)
{
  bool ret(false);
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::iterator it=_fields.begin();it!=_fields.end();it++)
    {
      MEDFileAnyTypeFieldMultiTSWithoutSDA *fmts(*it);
      if(fmts)
        {
          ret=fmts->renumberEntitiesLyingOnMesh(meshName,oldCode,newCode,renumO2N,*this) || ret;
        }
    }
  return ret;
}

/*!
 * Return an extraction of \a this using \a extractDef map to specify the extraction.
 * The keys of \a extractDef is level relative to max ext of \a mm mesh.
 *
 * \return A new object that the caller is responsible to deallocate.
 */
MEDFileFields *MEDFileFields::extractPart(const std::map<int, MCAuto<DataArrayInt> >& extractDef, MEDFileMesh *mm) const
{
  if(!mm)
    throw INTERP_KERNEL::Exception("MEDFileFields::extractPart : input mesh is NULL !");
  MCAuto<MEDFileFields> fsOut(MEDFileFields::New());
  int nbFields(getNumberOfFields());
  for(int i=0;i<nbFields;i++)
    {
      MCAuto<MEDFileAnyTypeFieldMultiTS> fmts(getFieldAtPos(i));
      if(!fmts)
        {
          std::ostringstream oss; oss << "MEDFileFields::extractPart : at pos #" << i << " field is null !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      MCAuto<MEDFileAnyTypeFieldMultiTS> fmtsOut(fmts->extractPart(extractDef,mm));
      fsOut->pushField(fmtsOut);
    }
  return fsOut.retn();
}

void MEDFileFields::accept(MEDFileFieldVisitor& visitor) const
{
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    if((*it).isNotNull())
      {
        visitor.newFieldEntry(*it);
        (*it)->accept(visitor);
        visitor.endFieldEntry(*it);
      }
}

class PFLData
{
public:
  PFLData():_add_pts_in_pfl(0) { }
  PFLData(const MCAuto<DataArrayInt>& mat, const MCAuto<DataArrayInt>& pfl, int nbOfNewPts):_matrix(mat),_pfl(pfl),_add_pts_in_pfl(nbOfNewPts) { }
  std::string getPflName() const { if(_pfl.isNull()) { return std::string(); } else { return _pfl->getName(); } }
  int getNbOfAddPtsInPfl() const { return _add_pts_in_pfl; }
  MCAuto<DataArrayInt> getProfile() const { return _pfl; }
  MCAuto<DataArrayInt> getMatrix() const { return _matrix; }
private:
  MCAuto<DataArrayInt> _matrix;
  MCAuto<DataArrayInt> _pfl;
  int _add_pts_in_pfl;
};

class MEDFileFieldLin2QuadVisitor : public MEDFileFieldVisitor
{
public:
  MEDFileFieldLin2QuadVisitor(const MEDFileUMesh *lin, const MEDFileUMesh *quad, const MEDFileFieldGlobsReal *linGlobs, MEDFileFields* outFs):_lin(lin),_quad(quad),_lin_globs(linGlobs),_out_fs(outFs),_gt(INTERP_KERNEL::NORM_ERROR),_1ts_update_requested(false) { }
  void newFieldEntry(const MEDFileAnyTypeFieldMultiTSWithoutSDA *field) { if(field->getMeshName()!=_lin->getName()) return; _cur_fmts=MEDFileFieldMultiTS::New(); }
  void endFieldEntry(const MEDFileAnyTypeFieldMultiTSWithoutSDA *field) { if(_cur_fmts.isNotNull()) { if(_cur_fmts->getNumberOfTS()>0) _out_fs->pushField(_cur_fmts); } }
  //
  void newTimeStepEntry(const MEDFileAnyTypeField1TSWithoutSDA *ts);
  void endTimeStepEntry(const MEDFileAnyTypeField1TSWithoutSDA *ts);
  //
  void newMeshEntry(const MEDFileFieldPerMesh *fpm);
  void endMeshEntry(const MEDFileFieldPerMesh *fpm) { }
  //
  void newPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt);
  void endPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt) { }
  //
  void newPerMeshPerTypePerDisc(const MEDFileFieldPerMeshPerTypePerDisc *pmptpd);
private:
  void updateData(MEDFileFieldPerMeshPerTypePerDisc *pmtd, int deltaNbNodes);
private:
  const MEDFileUMesh *_lin;
  const MEDFileUMesh *_quad;
  const MEDFileFieldGlobsReal *_lin_globs;
  MEDFileFields *_out_fs;
  MCAuto<MEDFileFieldMultiTS> _cur_fmts;
  MCAuto<MEDFileField1TS> _cur_f1ts;
  INTERP_KERNEL::NormalizedCellType _gt;
  // Info on 1TS modification
  bool _1ts_update_requested;
  // Cache of matrix to compute faster the values on newly created points
  std::map< std::string, PFLData > _cache;
  std::vector<std::string> _pfls_to_be_updated;
};

void MEDFileFieldLin2QuadVisitor::newPerMeshPerTypePerDisc(const MEDFileFieldPerMeshPerTypePerDisc *pmptpd)
{
  if(_cur_f1ts.isNull())
    return;
  if(pmptpd->getType()!=ON_NODES)
    throw INTERP_KERNEL::Exception("Not managed yet for ON_CELLS ON_GAUSS_NE and ON_GAUSS_PT");
  _1ts_update_requested=true;
  MEDFileAnyTypeField1TSWithoutSDA *ct(_cur_f1ts->contentNotNullBase());
  int locId(pmptpd->getFather()->locIdOfLeaf(pmptpd));
  MEDFileFieldPerMeshPerTypePerDisc *pmtdToModify(ct->getLeafGivenMeshAndTypeAndLocId(_lin->getName(),_gt,locId));
  std::string pflName(pmptpd->getProfile());
  _pfls_to_be_updated.push_back(pflName);
  std::map< std::string, PFLData >::iterator itCache(_cache.find(pflName));
  if(itCache!=_cache.end())
    {
      updateData(pmtdToModify,(*itCache).second.getNbOfAddPtsInPfl());
      return ;
    }
  MCAuto<DataArrayInt> pfl;
  if(pflName.empty())
    pfl=DataArrayInt::Range(0,pmptpd->getNumberOfVals(),1);
  else
    pfl=_lin_globs->getProfile(pflName)->deepCopy();
  //
  MCAuto<MEDCouplingUMesh> mesh3D(_lin->getMeshAtLevel(0)),mesh3DQuadratic(_quad->getMeshAtLevel(0));
  MCAuto<DataArrayInt> cellIds(mesh3D->getCellIdsLyingOnNodes(pfl->begin(),pfl->end(),true));
  MCAuto<MEDCouplingUMesh> mesh3DQuadraticRestricted(mesh3DQuadratic->buildPartOfMySelf(cellIds->begin(),cellIds->end(),true));
  MCAuto<DataArrayInt> mesh3DQuadraticRestrictedNodeIds(mesh3DQuadraticRestricted->computeFetchedNodeIds());
  mesh3DQuadraticRestrictedNodeIds->checkMonotonic(true);
  MCAuto<DataArrayInt> newPtsIds(mesh3DQuadraticRestrictedNodeIds->buildSubstraction(pfl));
  MCAuto<MEDCoupling1SGTUMesh> allSeg3;
  {
    MCAuto<DataArrayInt> a,b,c,d;
    MCAuto<MEDCouplingUMesh> seg3Tmp(mesh3DQuadraticRestricted->explodeIntoEdges(a,b,c,d));
    allSeg3=MEDCoupling1SGTUMesh::New(seg3Tmp);
  }
  if(allSeg3->getCellModelEnum()!=INTERP_KERNEL::NORM_SEG3)
    throw INTERP_KERNEL::Exception("MEDFileFieldLin2QuadVisitor::newPerMeshPerTypePerDisc : invalid situation where SEG3 expected !");
  MCAuto<DataArrayInt> midPts,cellSeg3Ids,matrix;
  {
    DataArrayInt *nodeConn(allSeg3->getNodalConnectivity());
    nodeConn->rearrange(3);
    {
      std::vector<int> v(1,2);
      midPts=nodeConn->keepSelectedComponents(v);
    }
    cellSeg3Ids=DataArrayInt::FindPermutationFromFirstToSecond(midPts,newPtsIds);
    {
      std::vector<int> v(2); v[0]=0; v[1]=1;
      MCAuto<DataArrayInt> tmp(nodeConn->keepSelectedComponents(v));
      matrix=tmp->selectByTupleId(cellSeg3Ids->begin(),cellSeg3Ids->end());
    }
    nodeConn->rearrange(1);
  }
  MCAuto<DataArrayInt> pflq;
  if(!pflName.empty())
    {
      std::vector<const DataArrayInt *> vs(2);
      vs[0]=pfl; vs[1]=newPtsIds;
      pflq=DataArrayInt::Aggregate(vs);
      pflq->setName(pflName);
    }
  PFLData pdata(matrix,pflq,newPtsIds->getNumberOfTuples());
  _cache[pflName]=pdata;
  updateData(pmtdToModify,pdata.getNbOfAddPtsInPfl());
}

void MEDFileFieldLin2QuadVisitor::updateData(MEDFileFieldPerMeshPerTypePerDisc *pmtd, int deltaNbNodes)
{
  pmtd->incrementNbOfVals(deltaNbNodes);
}

void MEDFileFieldLin2QuadVisitor::newPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt)
{
  const MEDFileFieldPerMeshPerType *pmpt2(dynamic_cast<const MEDFileFieldPerMeshPerType *>(pmpt));
  if(!pmpt2)
    throw INTERP_KERNEL::Exception("MEDFileFieldLin2QuadVisitor::newPerMeshPerTypeEntry : not managed for structure elements !");
  if(pmpt2->getNumberOfLoc()!=1)
    throw INTERP_KERNEL::Exception("MEDFileFieldLin2QuadVisitor::newPerMeshPerTypeEntry : not managed for multi discr per timestep !");
  _gt=pmpt->getGeoType();
}

void MEDFileFieldLin2QuadVisitor::newMeshEntry(const MEDFileFieldPerMesh *fpm)
{
  if(fpm->getMeshName()!=_lin->getName())
    throw INTERP_KERNEL::Exception("MEDFileFieldLin2QuadVisitor::newMeshEntry : mismatch into meshName !");
}

void MEDFileFieldLin2QuadVisitor::newTimeStepEntry(const MEDFileAnyTypeField1TSWithoutSDA *ts)
{
  _1ts_update_requested=false; _pfls_to_be_updated.clear();
  if(!ts)
    return ;
  const MEDFileField1TSWithoutSDA *tsd(dynamic_cast<const MEDFileField1TSWithoutSDA *>(ts));
  if(!tsd)
    return ;
  MCAuto<MEDFileAnyTypeField1TSWithoutSDA> contentCpy(ts->deepCopy());
  MCAuto<MEDFileField1TSWithoutSDA> contentCpy2(DynamicCastSafe<MEDFileAnyTypeField1TSWithoutSDA,MEDFileField1TSWithoutSDA>(contentCpy));
  if(contentCpy2.isNull())
    return;
  _cur_f1ts=MEDFileField1TS::New(*contentCpy2,true);
  _cur_f1ts->deepCpyGlobs(*_lin_globs);
}

void MEDFileFieldLin2QuadVisitor::endTimeStepEntry(const MEDFileAnyTypeField1TSWithoutSDA *ts)
{
  if(_cur_f1ts.isNull())
    return ;
  if(_1ts_update_requested)
    {
      MCAuto<DataArrayInt> matrix,oldPfl;
      for(std::vector<std::string>::const_iterator it=_pfls_to_be_updated.begin();it!=_pfls_to_be_updated.end();it++)
        {
          std::map< std::string, PFLData >::const_iterator it2(_cache.find(*it));
          if(it2==_cache.end())
            throw INTERP_KERNEL::Exception("MEDFileFieldLin2QuadVisitor::endTimeStepEntry : invalid situation !");
          matrix=(*it2).second.getMatrix();
          if((*it).empty())
            continue;
          int locId(_cur_f1ts->getProfileId(*it));
          oldPfl.takeRef(_cur_f1ts->getProfile(*it));
          {
            std::vector<int> locToKill(1,locId);
            _cur_f1ts->killProfileIds(locToKill);
          }
          _cur_f1ts->appendProfile((*it2).second.getProfile());
        }
      DataArrayDouble *arr(_cur_f1ts->getUndergroundDataArray());
      MCAuto<DataArrayDouble> res;
      {
        std::vector<int> v(1,0),v2(1,1);
        MCAuto<DataArrayInt> pts0(matrix->keepSelectedComponents(v));
        MCAuto<DataArrayInt> pts1(matrix->keepSelectedComponents(v2));
        if(oldPfl.isNotNull())
          {
            pts0=oldPfl->findIdForEach(pts0->begin(),pts0->end());
            pts1=oldPfl->findIdForEach(pts1->begin(),pts1->end());
          }
        MCAuto<DataArrayDouble> part0(arr->selectByTupleId(*pts0));
        MCAuto<DataArrayDouble> part1(arr->selectByTupleId(*pts1));
        res=DataArrayDouble::Add(part0,part1);
        res->applyLin(0.5,0.);
      }
      res=DataArrayDouble::Aggregate(arr,res);
      _cur_f1ts->setArray(res);
    }
  if(_cur_fmts.isNotNull())
    { _cur_fmts->pushBackTimeStep(_cur_f1ts); }
  _1ts_update_requested=false;
}

/*!
 * \a newQuad is expected to be the result of MEDFileUMesh::linearToQuadratic of \a oldLin
 */
MCAuto<MEDFileFields> MEDFileFields::linearToQuadratic(const MEDFileMeshes *oldLin, const MEDFileMeshes *newQuad) const
{
  if(!oldLin || !newQuad)
    throw INTERP_KERNEL::Exception("MEDFileFields::linearToQuadratic : input meshes must be non NULL !");
  MCAuto<MEDFileFields> ret(MEDFileFields::New());
  for(int i=0;i<oldLin->getNumberOfMeshes();i++)
    {
      MEDFileMesh *mm(oldLin->getMeshAtPos(i));
      if(!mm)
        continue;
      MEDFileUMesh *mmu(dynamic_cast<MEDFileUMesh *>(mm));
      if(!mmu)
        continue;
      MEDFileMesh *mmq(newQuad->getMeshWithName(mmu->getName()));
      MEDFileUMesh *mmqu(dynamic_cast<MEDFileUMesh *>(mmq));
      if(!mmqu)
        {
          std::ostringstream oss; oss << "MEDFileFields::linearToQuadratic : mismatch of name between input meshes for name \"" << mmu->getName() << "\"";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      MEDFileFieldLin2QuadVisitor vis(mmu,mmqu,this,ret);
      accept(vis);
    }
  return ret;
}

MEDFileAnyTypeFieldMultiTS *MEDFileFields::getFieldAtPos(int i) const
{
  if(i<0 || i>=(int)_fields.size())
    {
      std::ostringstream oss; oss << "MEDFileFields::getFieldAtPos : Invalid given id in input (" << i << ") should be in [0," << _fields.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *fmts=_fields[i];
  if(!fmts)
    return 0;
  MCAuto<MEDFileAnyTypeFieldMultiTS> ret;
  const MEDFileFieldMultiTSWithoutSDA *fmtsC(dynamic_cast<const MEDFileFieldMultiTSWithoutSDA *>(fmts));
  const MEDFileIntFieldMultiTSWithoutSDA *fmtsC2(dynamic_cast<const MEDFileIntFieldMultiTSWithoutSDA *>(fmts));
  const MEDFileFloatFieldMultiTSWithoutSDA *fmtsC3(dynamic_cast<const MEDFileFloatFieldMultiTSWithoutSDA *>(fmts));
  if(fmtsC)
    ret=MEDFileFieldMultiTS::New(*fmtsC,false);
  else if(fmtsC2)
    ret=MEDFileIntFieldMultiTS::New(*fmtsC2,false);
  else if(fmtsC3)
    ret=MEDFileFloatFieldMultiTS::New(*fmtsC3,false);
  else
    {
      std::ostringstream oss; oss << "MEDFileFields::getFieldAtPos : At pos #" << i << " field is neither double (FLOAT64) nor float (FLOAT32) nor integer (INT32) !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  ret->shallowCpyGlobs(*this);
  return ret.retn();
}

/*!
 * Return a shallow copy of \a this reduced to the fields ids defined in [ \a startIds , endIds ).
 * This method is accessible in python using __getitem__ with a list in input.
 * \return a new object that the caller should deal with.
 */
MEDFileFields *MEDFileFields::buildSubPart(const int *startIds, const int *endIds) const
{
  MCAuto<MEDFileFields> ret=shallowCpy();
  std::size_t sz=std::distance(startIds,endIds);
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > fields(sz);
  int j=0;
  for(const int *i=startIds;i!=endIds;i++,j++)
    {
      if(*i<0 || *i>=(int)_fields.size())
        {
          std::ostringstream oss; oss << "MEDFileFields::buildSubPart : Invalid given id in input (" << *i << ") should be in [0," << _fields.size() << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      fields[j]=_fields[*i];
    }
  ret->_fields=fields;
  return ret.retn();
}

MEDFileAnyTypeFieldMultiTS *MEDFileFields::getFieldWithName(const std::string& fieldName) const
{
  return getFieldAtPos(getPosFromFieldName(fieldName));
}

/*!
 * This method removes, if any, fields in \a this having no time steps.
 * If there is one or more than one such field in \a this true is returned and those fields will not be referenced anymore in \a this.
 * 
 * If false is returned \a this does not contain such fields. If false is returned this method can be considered as const.
 */
bool MEDFileFields::removeFieldsWithoutAnyTimeStep()
{
  std::vector<MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > newFields;
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *elt(*it);
      if(elt)
        {
          if(elt->getNumberOfTS()>0)
            newFields.push_back(*it);
        }
    }
  if(_fields.size()==newFields.size())
    return false;
  _fields=newFields;
  return true;
}

/*!
 * This method returns a new object containing part of \a this fields lying on mesh name specified by the input parameter \a meshName.
 * This method can be seen as a filter applied on \a this, that returns an object containing
 * reduced the list of fields compared to those in \a this. The returned object is a new object but the object on which it lies are only
 * shallow copied from \a this.
 * 
 * \param [in] meshName - the name of the mesh on w
 * \return a new object that the caller should deal with.
 */
MEDFileFields *MEDFileFields::partOfThisLyingOnSpecifiedMeshName(const std::string& meshName) const
{
  MCAuto<MEDFileFields> ret(MEDFileFields::New());
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *cur(*it);
      if(!cur)
        continue;
      if(cur->getMeshName()==meshName)
        {
          cur->incrRef();
          MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> cur2(const_cast<MEDFileAnyTypeFieldMultiTSWithoutSDA *>(cur));
          ret->_fields.push_back(cur2);
        }
    }
  ret->shallowCpyOnlyUsedGlobs(*this);
  return ret.retn();
}

/*!
 * This method returns a new object containing part of \a this fields lying ** exactly ** on the time steps specified by input parameter \a timeSteps.
 * Input time steps are specified using a pair of integer (iteration, order).
 * This method can be seen as a filter applied on \a this, that returns an object containing the same number of fields than those in \a this,
 * but for each multitimestep only the time steps in \a timeSteps are kept.
 * Typically the input parameter \a timeSteps comes from the call of MEDFileFields::getCommonIterations.
 * 
 * The returned object points to shallow copy of elements in \a this.
 * 
 * \param [in] timeSteps - the time steps given by a vector of pair of integers (iteration,order)
 * \throw If there is a field in \a this that is \b not defined on a time step in the input \a timeSteps.
 * \sa MEDFileFields::getCommonIterations, MEDFileFields::partOfThisNotLyingOnSpecifiedTimeSteps
 */
MEDFileFields *MEDFileFields::partOfThisLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const
{
  MCAuto<MEDFileFields> ret(MEDFileFields::New());
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *cur(*it);
      if(!cur)
        continue;
      MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> elt=cur->partOfThisLyingOnSpecifiedTimeSteps(timeSteps);
      ret->_fields.push_back(elt);
    }
  ret->shallowCpyOnlyUsedGlobs(*this);
  return ret.retn();
}

/*!
 * \sa MEDFileFields::getCommonIterations, MEDFileFields::partOfThisLyingOnSpecifiedTimeSteps
 */
MEDFileFields *MEDFileFields::partOfThisNotLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const
{
  MCAuto<MEDFileFields> ret=MEDFileFields::New();
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *cur(*it);
      if(!cur)
        continue;
      MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> elt=cur->partOfThisNotLyingOnSpecifiedTimeSteps(timeSteps);
      if(elt->getNumberOfTS()!=0)
        ret->_fields.push_back(elt);
    }
  ret->shallowCpyOnlyUsedGlobs(*this);
  return ret.retn();
}

bool MEDFileFields::presenceOfStructureElements() const
{
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    if((*it).isNotNull())
      if((*it)->presenceOfStructureElements())
        return true;
  return false;
}

void MEDFileFields::killStructureElements()
{
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > ret;
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::iterator it=_fields.begin();it!=_fields.end();it++)
    if((*it).isNotNull())
      {
        if((*it)->presenceOfStructureElements())
          {
            if(!(*it)->onlyStructureElements())
              {
                (*it)->killStructureElements();
                ret.push_back(*it);
              }
          }
        else
          {
            ret.push_back(*it);
          }
      }
  _fields=ret;
}

void MEDFileFields::keepOnlyStructureElements()
{
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > ret;
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::iterator it=_fields.begin();it!=_fields.end();it++)
    if((*it).isNotNull())
      {
        if((*it)->presenceOfStructureElements())
          {
            if(!(*it)->onlyStructureElements())
              (*it)->keepOnlyStructureElements();
            ret.push_back(*it);
          }
      }
  _fields=ret;
}

void MEDFileFields::keepOnlyOnMeshSE(const std::string& meshName, const std::string& seName)
{
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > ret;
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::iterator it=_fields.begin();it!=_fields.end();it++)
    if((*it).isNotNull())
      {
        if((*it)->getMeshName()!=meshName)
          continue;
        std::vector< std::pair<std::string,std::string> > ps;
        (*it)->getMeshSENames(ps);
        std::pair<std::string,std::string> p(meshName,seName);
        if(std::find(ps.begin(),ps.end(),p)!=ps.end())
          (*it)->keepOnlyOnSE(seName);
        ret.push_back(*it);
      }
  _fields=ret;
}

void MEDFileFields::getMeshSENames(std::vector< std::pair<std::string,std::string> >& ps) const
{
  for(std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    if((*it).isNotNull())
      (*it)->getMeshSENames(ps);
}

void MEDFileFields::blowUpSE(MEDFileMeshes *ms, const MEDFileStructureElements *ses)
{
  MEDFileBlowStrEltUp::DealWithSE(this,ms,ses);
}

MCAuto<MEDFileFields> MEDFileFields::partOfThisOnStructureElements() const
{
  MCAuto<MEDFileFields> ret(deepCopy());
  ret->keepOnlyStructureElements();
  return ret;
}

MCAuto<MEDFileFields> MEDFileFields::partOfThisLyingOnSpecifiedMeshSEName(const std::string& meshName, const std::string& seName) const
{
  MCAuto<MEDFileFields> ret(deepCopy());
  ret->keepOnlyOnMeshSE(meshName,seName);
  return ret;
}

void MEDFileFields::aggregate(const MEDFileFields& other)
{
  int nbFieldsToAdd(other.getNumberOfFields());
  std::vector<std::string> fsn(getFieldsNames());
  for(int i=0;i<nbFieldsToAdd;i++)
    {
      MCAuto<MEDFileAnyTypeFieldMultiTS> elt(other.getFieldAtPos(i));
      std::string name(elt->getName());
      if(std::find(fsn.begin(),fsn.end(),name)!=fsn.end())
        {
          std::ostringstream oss; oss << "MEDFileFields::aggregate : name \"" << name << "\" already appears !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      pushField(elt);
    }
}

MEDFileFieldsIterator *MEDFileFields::iterator()
{
  return new MEDFileFieldsIterator(this);
}

int MEDFileFields::getPosFromFieldName(const std::string& fieldName) const
{
  std::string tmp(fieldName);
  std::vector<std::string> poss;
  for(std::size_t i=0;i<_fields.size();i++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *f(_fields[i]);
      if(f)
        {
          std::string fname(f->getName());
          if(tmp==fname)
            return i;
          else
            poss.push_back(fname);
        }
    }
  std::ostringstream oss; oss << "MEDFileFields::getPosFromFieldName : impossible to find field '" << tmp << "' in this ! Possibilities are : ";
  std::copy(poss.begin(),poss.end(),std::ostream_iterator<std::string>(oss,", "));
  oss << " !";
  throw INTERP_KERNEL::Exception(oss.str());
}

MEDFileFieldsIterator::MEDFileFieldsIterator(MEDFileFields *fs):_fs(fs),_iter_id(0),_nb_iter(0)
{
  if(fs)
    {
      fs->incrRef();
      _nb_iter=fs->getNumberOfFields();
    }
}

MEDFileFieldsIterator::~MEDFileFieldsIterator() 
{
}

MEDFileAnyTypeFieldMultiTS *MEDFileFieldsIterator::nextt()
{
  if(_iter_id<_nb_iter)
    {
      MEDFileFields *fs(_fs);
      if(fs)
        return fs->getFieldAtPos(_iter_id++);
      else
        return 0;
    }
  else
    return 0;
}
