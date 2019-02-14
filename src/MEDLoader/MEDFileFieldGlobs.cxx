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

#include "MEDFileFieldGlobs.hxx"
#include "MEDFileField.txx"
#include "MEDFileMesh.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDLoaderTraits.hxx"
#include "MEDFileSafeCaller.txx"
#include "MEDFileFieldOverView.hxx"
#include "MEDFileBlowStrEltUp.hxx"
#include "MEDFileFieldVisitor.hxx"

#include "MEDCouplingFieldDiscretization.hxx"
#include "MCType.hxx"

#include "InterpKernelAutoPtr.hxx"
#include "CellModel.hxx"

#include <algorithm>
#include <iterator>

using namespace MEDCoupling;

void MEDFileFieldGlobs::loadProfileInFile(med_idt fid, int id, const std::string& pflName)
{
  if(id>=(int)_pfls.size())
    _pfls.resize(id+1);
  _pfls[id]=DataArrayInt::New();
  int lgth(MEDprofileSizeByName(fid,pflName.c_str()));
  _pfls[id]->setName(pflName);
  _pfls[id]->alloc(lgth,1);
  MEDFILESAFECALLERRD0(MEDprofileRd,(fid,pflName.c_str(),_pfls[id]->getPointer()));
  _pfls[id]->applyLin(1,-1,0);//Converting into C format
}

void MEDFileFieldGlobs::loadProfileInFile(med_idt fid, int i)
{
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  int sz;
  MEDFILESAFECALLERRD0(MEDprofileInfo,(fid,i+1,pflName,&sz));
  std::string pflCpp=MEDLoaderBase::buildStringFromFortran(pflName,MED_NAME_SIZE);
  if(i>=(int)_pfls.size())
    _pfls.resize(i+1);
  _pfls[i]=DataArrayInt::New();
  _pfls[i]->alloc(sz,1);
  _pfls[i]->setName(pflCpp.c_str());
  MEDFILESAFECALLERRD0(MEDprofileRd,(fid,pflName,_pfls[i]->getPointer()));
  _pfls[i]->applyLin(1,-1,0);//Converting into C format
}

void MEDFileFieldGlobs::writeGlobals(med_idt fid, const MEDFileWritable& opt) const
{
  int nbOfPfls=_pfls.size();
  for(int i=0;i<nbOfPfls;i++)
    {
      MCAuto<DataArrayInt> cpy=_pfls[i]->deepCopy();
      cpy->applyLin(1,1,0);
      INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      MEDLoaderBase::safeStrCpy(_pfls[i]->getName().c_str(),MED_NAME_SIZE,pflName,opt.getTooLongStrPolicy());
      MEDFILESAFECALLERWR0(MEDprofileWr,(fid,pflName,_pfls[i]->getNumberOfTuples(),cpy->getConstPointer()));
    }
  //
  int nbOfLocs=_locs.size();
  for(int i=0;i<nbOfLocs;i++)
    _locs[i]->writeLL(fid);
}

void MEDFileFieldGlobs::appendGlobs(const MEDFileFieldGlobs& other, double eps)
{
  std::vector<std::string> pfls=getPfls();
  for(std::vector< MCAuto<DataArrayInt> >::const_iterator it=other._pfls.begin();it!=other._pfls.end();it++)
    {
      std::vector<std::string>::iterator it2=std::find(pfls.begin(),pfls.end(),(*it)->getName());
      if(it2==pfls.end())
        {
          _pfls.push_back(*it);
        }
      else
        {
          int id=std::distance(pfls.begin(),it2);
          if(!(*it)->isEqual(*_pfls[id]))
            {
              std::ostringstream oss; oss << "MEDFileFieldGlobs::appendGlobs : Profile \"" << (*it)->getName() << "\" already exists and is different from those expecting to be append !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
        }
    }
  std::vector<std::string> locs=getLocs();
  for(std::vector< MCAuto<MEDFileFieldLoc> >::const_iterator it=other._locs.begin();it!=other._locs.end();it++)
    {
      std::vector<std::string>::iterator it2=std::find(locs.begin(),locs.end(),(*it)->getName());
      if(it2==locs.end())
        {
          _locs.push_back(*it);
        }
      else
        {
          int id=std::distance(locs.begin(),it2);
          if(!(*it)->isEqual(*_locs[id],eps))
            {
              std::ostringstream oss; oss << "MEDFileFieldGlobs::appendGlobs : Localization \"" << (*it)->getName() << "\" already exists and is different from those expecting to be append !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
        }
    }
}

void MEDFileFieldGlobs::checkGlobsPflsPartCoherency(const std::vector<std::string>& pflsUsed) const
{
  for(std::vector<std::string>::const_iterator it=pflsUsed.begin();it!=pflsUsed.end();it++)
    getProfile((*it).c_str());
}

void MEDFileFieldGlobs::checkGlobsLocsPartCoherency(const std::vector<std::string>& locsUsed) const
{
  for(std::vector<std::string>::const_iterator it=locsUsed.begin();it!=locsUsed.end();it++)
    getLocalization((*it).c_str());
}

void MEDFileFieldGlobs::loadGlobals(med_idt fid, const MEDFileFieldGlobsReal& real)
{
  std::vector<std::string> profiles=real.getPflsReallyUsed();
  int sz=profiles.size();
  _pfls.resize(sz);
  for(int i=0;i<sz;i++)
    loadProfileInFile(fid,i,profiles[i].c_str());
  //
  std::vector<std::string> locs=real.getLocsReallyUsed();
  sz=locs.size();
  _locs.resize(sz);
  for(int i=0;i<sz;i++)
    _locs[i]=MEDFileFieldLoc::New(fid,locs[i].c_str());
}

void MEDFileFieldGlobs::loadAllGlobals(med_idt fid, const MEDFileEntities *entities)
{
  int nProfil=MEDnProfile(fid);
  for(int i=0;i<nProfil;i++)
    loadProfileInFile(fid,i);
  int sz=MEDnLocalization(fid);
  _locs.resize(sz);
  for(int i=0;i<sz;i++)
    {
      _locs[i]=MEDFileFieldLoc::New(fid,i,entities);
    }
}

MEDFileFieldGlobs *MEDFileFieldGlobs::New(med_idt fid)
{
  return new MEDFileFieldGlobs(fid);
}

MEDFileFieldGlobs *MEDFileFieldGlobs::New()
{
  return new MEDFileFieldGlobs;
}

std::size_t MEDFileFieldGlobs::getHeapMemorySizeWithoutChildren() const
{
  return _file_name.capacity()+_pfls.capacity()*sizeof(MCAuto<DataArrayInt>)+_locs.capacity()*sizeof(MCAuto<MEDFileFieldLoc>);
}

std::vector<const BigMemoryObject *> MEDFileFieldGlobs::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MCAuto< DataArrayInt > >::const_iterator it=_pfls.begin();it!=_pfls.end();it++)
    ret.push_back((const DataArrayInt *)*it);
  for(std::vector< MCAuto<MEDFileFieldLoc> >::const_iterator it=_locs.begin();it!=_locs.end();it++)
    ret.push_back((const MEDFileFieldLoc *)*it);
  return ret;
}

MEDFileFieldGlobs *MEDFileFieldGlobs::deepCopy() const
{
  MCAuto<MEDFileFieldGlobs> ret=new MEDFileFieldGlobs(*this);
  std::size_t i=0;
  for(std::vector< MCAuto<DataArrayInt> >::const_iterator it=_pfls.begin();it!=_pfls.end();it++,i++)
    {
      if((const DataArrayInt *)*it)
        ret->_pfls[i]=(*it)->deepCopy();
    }
  i=0;
  for(std::vector< MCAuto<MEDFileFieldLoc> >::const_iterator it=_locs.begin();it!=_locs.end();it++,i++)
    {
      if((const MEDFileFieldLoc*)*it)
        ret->_locs[i]=(*it)->deepCopy();
    }
  return ret.retn();
}

/*!
 * \throw if a profile in \a pfls in not in \a this.
 * \throw if a localization in \a locs in not in \a this.
 * \sa MEDFileFieldGlobs::deepCpyPart
 */
MEDFileFieldGlobs *MEDFileFieldGlobs::shallowCpyPart(const std::vector<std::string>& pfls, const std::vector<std::string>& locs) const
{
  MCAuto<MEDFileFieldGlobs> ret=MEDFileFieldGlobs::New();
  for(std::vector<std::string>::const_iterator it1=pfls.begin();it1!=pfls.end();it1++)
    {
      DataArrayInt *pfl=const_cast<DataArrayInt *>(getProfile((*it1).c_str()));
      if(!pfl)
        throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::shallowCpyPart : internal error ! pfl null !");
      pfl->incrRef();
      MCAuto<DataArrayInt> pfl2(pfl);
      ret->_pfls.push_back(pfl2);
    }
  for(std::vector<std::string>::const_iterator it2=locs.begin();it2!=locs.end();it2++)
    {
      MEDFileFieldLoc *loc=const_cast<MEDFileFieldLoc *>(&getLocalization((*it2).c_str()));
      if(!loc)
        throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::shallowCpyPart : internal error ! loc null !");
      loc->incrRef();
      MCAuto<MEDFileFieldLoc> loc2(loc);
      ret->_locs.push_back(loc2);
    }
  ret->setFileName(getFileName());
  return ret.retn();
}

/*!
 * \throw if a profile in \a pfls in not in \a this.
 * \throw if a localization in \a locs in not in \a this.
 * \sa MEDFileFieldGlobs::shallowCpyPart
 */
MEDFileFieldGlobs *MEDFileFieldGlobs::deepCpyPart(const std::vector<std::string>& pfls, const std::vector<std::string>& locs) const
{
  MCAuto<MEDFileFieldGlobs> ret=MEDFileFieldGlobs::New();
  for(std::vector<std::string>::const_iterator it1=pfls.begin();it1!=pfls.end();it1++)
    {
      DataArrayInt *pfl=const_cast<DataArrayInt *>(getProfile((*it1).c_str()));
      if(!pfl)
        throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::deepCpyPart : internal error ! pfl null !");
      ret->_pfls.push_back(pfl->deepCopy());
    }
  for(std::vector<std::string>::const_iterator it2=locs.begin();it2!=locs.end();it2++)
    {
      MEDFileFieldLoc *loc=const_cast<MEDFileFieldLoc *>(&getLocalization((*it2).c_str()));
      if(!loc)
        throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::deepCpyPart : internal error ! loc null !");
      ret->_locs.push_back(loc->deepCopy());
    }
  ret->setFileName(getFileName());
  return ret.retn();
}

MEDFileFieldGlobs::MEDFileFieldGlobs(med_idt fid):_file_name(MEDFileWritable::FileNameFromFID(fid))
{
}

MEDFileFieldGlobs::MEDFileFieldGlobs()
{
}

MEDFileFieldGlobs::~MEDFileFieldGlobs()
{
}

void MEDFileFieldGlobs::simpleRepr(std::ostream& oss) const
{
  oss << "Profiles :\n";
  std::size_t n=_pfls.size();
  for(std::size_t i=0;i<n;i++)
    {
      oss << "  - #" << i << " ";
      const DataArrayInt *pfl=_pfls[i];
      if(pfl)
        oss << "\"" << pfl->getName() << "\"\n";
      else
        oss << "EMPTY !\n";
    }
  n=_locs.size();
  oss << "Localizations :\n";
  for(std::size_t i=0;i<n;i++)
    {
      oss << "  - #" << i << " ";
      const MEDFileFieldLoc *loc=_locs[i];
      if(loc)
        loc->simpleRepr(oss);
      else
        oss<< "EMPTY !\n";
    }
}

void MEDFileFieldGlobs::changePflsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  for(std::vector< MCAuto<DataArrayInt> >::iterator it=_pfls.begin();it!=_pfls.end();it++)
    {
      DataArrayInt *elt(*it);
      if(elt)
        {
          std::string name(elt->getName());
          for(std::vector< std::pair<std::vector<std::string>, std::string > >::const_iterator it2=mapOfModif.begin();it2!=mapOfModif.end();it2++)
            {
              if(std::find((*it2).first.begin(),(*it2).first.end(),name)!=(*it2).first.end())
                {
                  elt->setName((*it2).second.c_str());
                  return;
                }
            }
        }
    }
}

void MEDFileFieldGlobs::changeLocsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  for(std::vector< MCAuto<MEDFileFieldLoc> >::iterator it=_locs.begin();it!=_locs.end();it++)
    {
      MEDFileFieldLoc *elt(*it);
      if(elt)
        {
          std::string name(elt->getName());
          for(std::vector< std::pair<std::vector<std::string>, std::string > >::const_iterator it2=mapOfModif.begin();it2!=mapOfModif.end();it2++)
            {
              if(std::find((*it2).first.begin(),(*it2).first.end(),name)!=(*it2).first.end())
                {
                  elt->setName((*it2).second.c_str());
                  return;
                }
            }
        }
    }
}

int MEDFileFieldGlobs::getNbOfGaussPtPerCell(int locId) const
{
  if(locId<0 || locId>=(int)_locs.size())
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::getNbOfGaussPtPerCell : Invalid localization id !");
  return _locs[locId]->getNbOfGaussPtPerCell();
}

const MEDFileFieldLoc& MEDFileFieldGlobs::getLocalization(const std::string& locName) const
{
  return getLocalizationFromId(getLocalizationId(locName));
}

const MEDFileFieldLoc& MEDFileFieldGlobs::getLocalizationFromId(int locId) const
{
  if(locId<0 || locId>=(int)_locs.size())
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::getLocalizationFromId : Invalid localization id !");
  return *_locs[locId];
}

/// @cond INTERNAL
namespace MEDCouplingImpl
{
  class LocFinder
  {
  public:
    LocFinder(const std::string& loc):_loc(loc) { }
    bool operator() (const MCAuto<MEDFileFieldLoc>& loc) { return loc->isName(_loc); }
  private:
    const std::string &_loc;
  };

  class PflFinder
  {
  public:
    PflFinder(const std::string& pfl):_pfl(pfl) { }
    bool operator() (const MCAuto<DataArrayInt>& loc) { return loc->getName()==_pfl; }
  private:
    const std::string _pfl;
  };
}
/// @endcond

int MEDFileFieldGlobs::getLocalizationId(const std::string& loc) const
{
  std::vector< MCAuto<MEDFileFieldLoc> >::const_iterator it=std::find_if(_locs.begin(),_locs.end(),MEDCouplingImpl::LocFinder(loc));
  if(it==_locs.end())
    {
      std::ostringstream oss; oss << "MEDFileFieldGlobs::getLocalisationId : no such localisation name : \"" << loc << "\" Possible localizations are : ";
      for(it=_locs.begin();it!=_locs.end();it++)
        oss << "\"" << (*it)->getName() << "\", ";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  return std::distance(_locs.begin(),it);
}

int MEDFileFieldGlobs::getProfileId(const std::string& pfl) const
{
  std::vector< MCAuto<DataArrayInt> >::const_iterator it=std::find_if(_pfls.begin(),_pfls.end(),MEDCouplingImpl::PflFinder(pfl));
  if(it==_pfls.end())
    {
      std::ostringstream oss; oss << "MEDFileFieldGlobs::getProfileId : no such profile name : \"" << pfl << "\" Possible localizations are : ";
      for(it=_pfls.begin();it!=_pfls.end();it++)
        oss << "\"" << (*it)->getName() << "\", ";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  return std::distance(_pfls.begin(),it);
}

/*!
 * The returned value is never null.
 */
const DataArrayInt *MEDFileFieldGlobs::getProfile(const std::string& pflName) const
{
  return getProfileFromId(getProfileId(pflName));
}

const DataArrayInt *MEDFileFieldGlobs::getProfileFromId(int pflId) const
{
  if(pflId<0 || pflId>=(int)_pfls.size())
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::getProfileFromId : Invalid profile id !");
  return _pfls[pflId];
}

MEDFileFieldLoc& MEDFileFieldGlobs::getLocalizationFromId(int locId)
{
  if(locId<0 || locId>=(int)_locs.size())
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::getLocalizationFromId : Invalid localization id !");
  return *_locs[locId];
}

MEDFileFieldLoc& MEDFileFieldGlobs::getLocalization(const std::string& locName)
{
  return getLocalizationFromId(getLocalizationId(locName));
}

/*!
 * The returned value is never null. Borrowed reference returned.
 */
DataArrayInt *MEDFileFieldGlobs::getProfile(const std::string& pflName)
{
  std::string pflNameCpp(pflName);
  std::vector< MCAuto<DataArrayInt> >::iterator it=std::find_if(_pfls.begin(),_pfls.end(),MEDCouplingImpl::PflFinder(pflNameCpp));
  if(it==_pfls.end())
    {
      std::ostringstream oss; oss << "MEDFileFieldGlobs::getProfile: no such profile name : \"" << pflNameCpp << "\" Possible profiles are : ";
      for(it=_pfls.begin();it!=_pfls.end();it++)
        oss << "\"" << (*it)->getName() << "\", ";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  return *it;
}

DataArrayInt *MEDFileFieldGlobs::getProfileFromId(int pflId)
{
  if(pflId<0 || pflId>=(int)_pfls.size())
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::getProfileFromId : Invalid profile id !");
  return _pfls[pflId];
}

void MEDFileFieldGlobs::killProfileIds(const std::vector<int>& pflIds)
{
  std::vector< MCAuto<DataArrayInt> > newPfls;
  int i=0;
  for(std::vector< MCAuto<DataArrayInt> >::const_iterator it=_pfls.begin();it!=_pfls.end();it++,i++)
    {
      if(std::find(pflIds.begin(),pflIds.end(),i)==pflIds.end())
        newPfls.push_back(*it);
    }
  _pfls=newPfls;
}

void MEDFileFieldGlobs::killLocalizationIds(const std::vector<int>& locIds)
{
  std::vector< MCAuto<MEDFileFieldLoc> > newLocs;
  int i=0;
  for(std::vector< MCAuto<MEDFileFieldLoc> >::const_iterator it=_locs.begin();it!=_locs.end();it++,i++)
    {
      if(std::find(locIds.begin(),locIds.end(),i)==locIds.end())
        newLocs.push_back(*it);
    }
  _locs=newLocs;
}

void MEDFileFieldGlobs::killStructureElementsInGlobs()
{
  std::vector< MCAuto<MEDFileFieldLoc> > newLocs;
  for(std::vector< MCAuto<MEDFileFieldLoc> >::iterator it=_locs.begin();it!=_locs.end();it++)
    {
      if((*it).isNull())
        continue;
      if(!(*it)->isOnStructureElement())
        newLocs.push_back(*it);
    }
  _locs=newLocs;
}

std::vector<std::string> MEDFileFieldGlobs::getPfls() const
{
  int sz=_pfls.size();
  std::vector<std::string> ret(sz);
  for(int i=0;i<sz;i++)
    ret[i]=_pfls[i]->getName();
  return ret;
}

std::vector<std::string> MEDFileFieldGlobs::getLocs() const
{
  int sz=_locs.size();
  std::vector<std::string> ret(sz);
  for(int i=0;i<sz;i++)
    ret[i]=_locs[i]->getName();
  return ret;
}

bool MEDFileFieldGlobs::existsPfl(const std::string& pflName) const
{
  std::vector<std::string> v=getPfls();
  std::string s(pflName);
  return std::find(v.begin(),v.end(),s)!=v.end();
}

bool MEDFileFieldGlobs::existsLoc(const std::string& locName) const
{
  std::vector<std::string> v=getLocs();
  std::string s(locName);
  return std::find(v.begin(),v.end(),s)!=v.end();
}

std::vector< std::vector<int> > MEDFileFieldGlobs::whichAreEqualProfiles() const
{
  std::map<int,std::vector<int> > m;
  int i=0;
  for(std::vector< MCAuto<DataArrayInt> >::const_iterator it=_pfls.begin();it!=_pfls.end();it++,i++)
    {
      const DataArrayInt *tmp=(*it);
      if(tmp)
        {
          m[tmp->getHashCode()].push_back(i);
        }
    }
  std::vector< std::vector<int> > ret;
  for(std::map<int,std::vector<int> >::const_iterator it2=m.begin();it2!=m.end();it2++)
    {
      if((*it2).second.size()>1)
        {
          std::vector<int> ret0;
          bool equalityOrNot=false;
          for(std::vector<int>::const_iterator it3=(*it2).second.begin();it3!=(*it2).second.end();it3++)
            {
              std::vector<int>::const_iterator it4=it3; it4++;
              for(;it4!=(*it2).second.end();it4++)
                {
                  if(_pfls[*it3]->isEqualWithoutConsideringStr(*_pfls[*it4]))
                    {
                      if(!equalityOrNot)
                        ret0.push_back(*it3);
                      ret0.push_back(*it4);
                      equalityOrNot=true;
                    }
                }
            }
          if(!ret0.empty())
            ret.push_back(ret0);
        }
    }
  return ret;
}

std::vector< std::vector<int> > MEDFileFieldGlobs::whichAreEqualLocs(double eps) const
{
  throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::whichAreEqualLocs : no implemented yet ! Sorry !");
}

void MEDFileFieldGlobs::appendProfile(DataArrayInt *pfl)
{
  std::string name(pfl->getName());
  if(name.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::appendProfile : unsupported profiles with no name !");
  for(std::vector< MCAuto<DataArrayInt> >::const_iterator it=_pfls.begin();it!=_pfls.end();it++)
    if(name==(*it)->getName())
      {
        if(!pfl->isEqual(*(*it)))
          {
            std::ostringstream oss; oss << "MEDFileFieldGlobs::appendProfile : profile \"" << name << "\" already exists and is different from existing !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
  pfl->incrRef();
  _pfls.push_back(pfl);
}

void MEDFileFieldGlobs::appendLoc(const std::string& locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w)
{
  std::string name(locName);
  if(name.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::appendLoc : unsupported localizations with no name !");
  MCAuto<MEDFileFieldLoc> obj=MEDFileFieldLoc::New(locName,geoType,refCoo,gsCoo,w);
  for(std::vector< MCAuto<MEDFileFieldLoc> >::const_iterator it=_locs.begin();it!=_locs.end();it++)
    if((*it)->isName(locName))
      {
        if(!(*it)->isEqual(*obj,1e-12))
          {
            std::ostringstream oss; oss << "MEDFileFieldGlobs::appendLoc : localization \"" << name << "\" already exists and is different from existing !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
  _locs.push_back(obj);
}

std::string MEDFileFieldGlobs::createNewNameOfPfl() const
{
  std::vector<std::string> names=getPfls();
  return CreateNewNameNotIn("NewPfl_",names);
}

std::string MEDFileFieldGlobs::createNewNameOfLoc() const
{
  std::vector<std::string> names=getLocs();
  return CreateNewNameNotIn("NewLoc_",names);
}

std::string MEDFileFieldGlobs::CreateNewNameNotIn(const std::string& prefix, const std::vector<std::string>& namesToAvoid)
{
  for(std::size_t sz=0;sz<100000;sz++)
    {
      std::ostringstream tryName;
      tryName << prefix << sz;
      if(std::find(namesToAvoid.begin(),namesToAvoid.end(),tryName.str())==namesToAvoid.end())
        return tryName.str();
    }
  throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::CreateNewNameNotIn : impossible to create an additional profile limit of 100000 profiles reached !");
}

/*!
 * Creates a MEDFileFieldGlobsReal on a given file name. Nothing is read here.
 *  \param [in] fname - the file name.
 */
MEDFileFieldGlobsReal::MEDFileFieldGlobsReal(med_idt fid):_globals(MEDFileFieldGlobs::New(fid))
{
}

/*!
 * Creates an empty MEDFileFieldGlobsReal.
 */
MEDFileFieldGlobsReal::MEDFileFieldGlobsReal():_globals(MEDFileFieldGlobs::New())
{
}

std::size_t MEDFileFieldGlobsReal::getHeapMemorySizeWithoutChildren() const
{
  return 0;
}

std::vector<const BigMemoryObject *> MEDFileFieldGlobsReal::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  ret.push_back((const MEDFileFieldGlobs *)_globals);
  return ret;
}

/*!
 * Returns a string describing profiles and Gauss points held in \a this.
 *  \return std::string - the description string.
 */
void MEDFileFieldGlobsReal::simpleReprGlobs(std::ostream& oss) const
{
  const MEDFileFieldGlobs *glob=_globals;
  std::ostringstream oss2; oss2 << glob;
  std::string stars(oss2.str().length(),'*');
  oss << "Globals information on fields (at " << oss2.str() << "):" << "\n************************************" << stars  << "\n\n";
  if(glob)
    glob->simpleRepr(oss);
  else
    oss << "NO GLOBAL INFORMATION !\n";
}

void MEDFileFieldGlobsReal::resetContent()
{
  _globals=MEDFileFieldGlobs::New();
}

void MEDFileFieldGlobsReal::killStructureElementsInGlobs()
{
  contentNotNull()->killStructureElementsInGlobs();
}

MEDFileFieldGlobsReal::~MEDFileFieldGlobsReal()
{
}

/*!
 * Copies references to profiles and Gauss points from another MEDFileFieldGlobsReal.
 *  \param [in] other - the other MEDFileFieldGlobsReal to copy data from.
 */
void MEDFileFieldGlobsReal::shallowCpyGlobs(const MEDFileFieldGlobsReal& other)
{
  _globals=other._globals;
}

/*!
 * Copies references to ** only used ** by \a this, profiles and Gauss points from another MEDFileFieldGlobsReal.
 *  \param [in] other - the other MEDFileFieldGlobsReal to copy data from.
 */
void MEDFileFieldGlobsReal::shallowCpyOnlyUsedGlobs(const MEDFileFieldGlobsReal& other)
{
  const MEDFileFieldGlobs *otherg(other._globals);
  if(!otherg)
    return ;
  _globals=otherg->shallowCpyPart(getPflsReallyUsed(),getLocsReallyUsed());
}

/*!
 * Copies deeply to ** only used ** by \a this, profiles and Gauss points from another MEDFileFieldGlobsReal.
 *  \param [in] other - the other MEDFileFieldGlobsReal to copy data from.
 */
void MEDFileFieldGlobsReal::deepCpyOnlyUsedGlobs(const MEDFileFieldGlobsReal& other)
{
  const MEDFileFieldGlobs *otherg(other._globals);
  if(!otherg)
    return ;
  _globals=otherg->deepCpyPart(getPflsReallyUsed(),getLocsReallyUsed());
}

void MEDFileFieldGlobsReal::deepCpyGlobs(const MEDFileFieldGlobsReal& other)
{
  _globals=other._globals;
  if((const MEDFileFieldGlobs *)_globals)
    _globals=other._globals->deepCopy();
}

/*!
 * Adds profiles and Gauss points held by another MEDFileFieldGlobsReal to \a this one.
 *  \param [in] other - the MEDFileFieldGlobsReal to copy data from.
 *  \param [in] eps - a precision used to compare Gauss points with same name held by
 *         \a this and \a other MEDFileFieldGlobsReal.
 *  \throw If \a this and \a other hold profiles with equal names but different ids.
 *  \throw If  \a this and \a other hold different Gauss points with equal names.
 */
void MEDFileFieldGlobsReal::appendGlobs(const MEDFileFieldGlobsReal& other, double eps)
{
  const MEDFileFieldGlobs *thisGlobals(_globals),*otherGlobals(other._globals);
  if(thisGlobals==otherGlobals)
    return ;
  if(!thisGlobals)
    {
      _globals=other._globals;
      return ;
    }
  _globals->appendGlobs(*other._globals,eps);
}

void MEDFileFieldGlobsReal::checkGlobsCoherency() const
{
  checkGlobsPflsPartCoherency();
  checkGlobsLocsPartCoherency();
}

void MEDFileFieldGlobsReal::checkGlobsPflsPartCoherency() const
{
  contentNotNull()->checkGlobsPflsPartCoherency(getPflsReallyUsed());
}

void MEDFileFieldGlobsReal::checkGlobsLocsPartCoherency() const
{
  contentNotNull()->checkGlobsLocsPartCoherency(getLocsReallyUsed());
}

void MEDFileFieldGlobsReal::loadProfileInFile(med_idt fid, int id, const std::string& pflName)
{
  contentNotNull()->loadProfileInFile(fid,id,pflName);
}

void MEDFileFieldGlobsReal::loadProfileInFile(med_idt fid, int id)
{
  contentNotNull()->loadProfileInFile(fid,id);
}

void MEDFileFieldGlobsReal::loadGlobals(med_idt fid)
{
  contentNotNull()->loadGlobals(fid,*this);
}

void MEDFileFieldGlobsReal::loadAllGlobals(med_idt fid, const MEDFileEntities *entities)
{
  contentNotNull()->loadAllGlobals(fid,entities);
}

void MEDFileFieldGlobsReal::writeGlobals(med_idt fid, const MEDFileWritable& opt) const
{
  contentNotNull()->writeGlobals(fid,opt);
}

/*!
 * Returns names of all profiles. To get only used profiles call getPflsReallyUsed()
 * or getPflsReallyUsedMulti().
 *  \return std::vector<std::string> - a sequence of names of all profiles.
 */
std::vector<std::string> MEDFileFieldGlobsReal::getPfls() const
{
  return contentNotNull()->getPfls();
}

/*!
 * Returns names of all localizations. To get only used localizations call getLocsReallyUsed()
 * or getLocsReallyUsedMulti().
 *  \return std::vector<std::string> - a sequence of names of all localizations.
 */
std::vector<std::string> MEDFileFieldGlobsReal::getLocs() const
{
  return contentNotNull()->getLocs();
}

/*!
 * Checks if the profile with a given name exists.
 *  \param [in] pflName - the profile name of interest.
 *  \return bool - \c true if the profile named \a pflName exists.
 */
bool MEDFileFieldGlobsReal::existsPfl(const std::string& pflName) const
{
  return contentNotNull()->existsPfl(pflName);
}

/*!
 * Checks if the localization with a given name exists.
 *  \param [in] locName - the localization name of interest.
 *  \return bool - \c true if the localization named \a locName exists.
 */
bool MEDFileFieldGlobsReal::existsLoc(const std::string& locName) const
{
  return contentNotNull()->existsLoc(locName);
}

std::string MEDFileFieldGlobsReal::createNewNameOfPfl() const
{
  return contentNotNull()->createNewNameOfPfl();
}

std::string MEDFileFieldGlobsReal::createNewNameOfLoc() const
{
  return contentNotNull()->createNewNameOfLoc();
}

/*!
 * Sets the name of a MED file.
 *  \param [inout] fileName - the file name.
 */
void MEDFileFieldGlobsReal::setFileName(const std::string& fileName)
{
  contentNotNull()->setFileName(fileName);
}

/*!
 * Finds equal profiles. Two profiles are considered equal if they contain the same ids
 * in the same order.
 *  \return std::vector< std::vector<int> > - a sequence of groups of equal profiles.
 *          Each item of this sequence is a vector containing ids of equal profiles.
 */
std::vector< std::vector<int> > MEDFileFieldGlobsReal::whichAreEqualProfiles() const
{
  return contentNotNull()->whichAreEqualProfiles();
}

/*!
 * Finds equal localizations.
 *  \param [in] eps - a precision used to compare real values of the localizations.
 *  \return std::vector< std::vector<int> > - a sequence of groups of equal localizations.
 *          Each item of this sequence is a vector containing ids of equal localizations.
 */
std::vector< std::vector<int> > MEDFileFieldGlobsReal::whichAreEqualLocs(double eps) const
{
  return contentNotNull()->whichAreEqualLocs(eps);
}

/*!
 * Renames the profiles. References to profiles (a reference is a profile name) are not changed.
 * \param [in] mapOfModif - a sequence describing required renaming. Each element of
 *        this sequence is a pair whose 
 *        - the first item is a vector of profile names to replace by the second item,
 *        - the second item is a profile name to replace every profile name of the first item.
 */
void MEDFileFieldGlobsReal::changePflsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  contentNotNull()->changePflsNamesInStruct(mapOfModif);
}

/*!
 * Renames the localizations. References to localizations (a reference is a localization name) are not changed.
 * \param [in] mapOfModif - a sequence describing required renaming. Each element of
 *        this sequence is a pair whose 
 *        - the first item is a vector of localization names to replace by the second item,
 *        - the second item is a localization name to replace every localization name of the first item.
 */
void MEDFileFieldGlobsReal::changeLocsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  contentNotNull()->changeLocsNamesInStruct(mapOfModif);
}

/*!
 * Replaces references to some profiles (a reference is a profile name) by references
 * to other profiles and, contrary to changePflsRefsNamesGen(), renames the profiles
 * them-selves accordingly. <br>
 * This method is a generalization of changePflName().
 * \param [in] mapOfModif - a sequence describing required replacements. Each element of
 *        this sequence is a pair whose 
 *        - the first item is a vector of profile names to replace by the second item,
 *        - the second item is a profile name to replace every profile of the first item.
 * \sa changePflsRefsNamesGen()
 * \sa changePflName()
 */
void MEDFileFieldGlobsReal::changePflsNames(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  changePflsRefsNamesGen(mapOfModif);
  changePflsNamesInStruct(mapOfModif);
}

/*!
 * Replaces references to some localizations (a reference is a localization name) by references
 * to other localizations and, contrary to changeLocsRefsNamesGen(), renames the localizations
 * them-selves accordingly. <br>
 * This method is a generalization of changeLocName().
 * \param [in] mapOfModif - a sequence describing required replacements. Each element of
 *        this sequence is a pair whose 
 *        - the first item is a vector of localization names to replace by the second item,
 *        - the second item is a localization name to replace every localization of the first item.
 * \sa changeLocsRefsNamesGen()
 * \sa changeLocName()
 */
void MEDFileFieldGlobsReal::changeLocsNames(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  changeLocsRefsNamesGen(mapOfModif);
  changeLocsNamesInStruct(mapOfModif);
}

/*!
 * Renames the profile having a given name and updates references to this profile.
 *  \param [in] oldName - the name of the profile to rename.
 *  \param [in] newName - a new name of the profile.
 * \sa changePflsNames().
 */
void MEDFileFieldGlobsReal::changePflName(const std::string& oldName, const std::string& newName)
{
  std::vector< std::pair<std::vector<std::string>, std::string > > mapOfModif(1);
  std::pair<std::vector<std::string>, std::string > p(std::vector<std::string>(1,std::string(oldName)),std::string(newName));
  mapOfModif[0]=p;
  changePflsNames(mapOfModif);
}

/*!
 * Renames the localization having a given name and updates references to this localization.
 *  \param [in] oldName - the name of the localization to rename.
 *  \param [in] newName - a new name of the localization.
 * \sa changeLocsNames().
 */
void MEDFileFieldGlobsReal::changeLocName(const std::string& oldName, const std::string& newName)
{
  std::vector< std::pair<std::vector<std::string>, std::string > > mapOfModif(1);
  std::pair<std::vector<std::string>, std::string > p(std::vector<std::string>(1,std::string(oldName)),std::string(newName));
  mapOfModif[0]=p;
  changeLocsNames(mapOfModif);
}

/*!
 * Removes duplicated profiles. Returns a map used to update references to removed 
 * profiles via changePflsRefsNamesGen().
 * Equal profiles are found using whichAreEqualProfiles().
 *  \return std::vector< std::pair<std::vector<std::string>, std::string > > - 
 *          a sequence describing the performed replacements of profiles. Each element of
 *          this sequence is a pair whose
 *          - the first item is a vector of profile names replaced by the second item,
 *          - the second item is a profile name replacing every profile of the first item.
 */
std::vector< std::pair<std::vector<std::string>, std::string > > MEDFileFieldGlobsReal::zipPflsNames()
{
  std::vector< std::vector<int> > pseudoRet=whichAreEqualProfiles();
  std::vector< std::pair<std::vector<std::string>, std::string > > ret(pseudoRet.size());
  int i=0;
  for(std::vector< std::vector<int> >::const_iterator it=pseudoRet.begin();it!=pseudoRet.end();it++,i++)
    {
      std::vector< std::string > tmp((*it).size());
      int j=0;
      for(std::vector<int>::const_iterator it2=(*it).begin();it2!=(*it).end();it2++,j++)
        tmp[j]=std::string(getProfileFromId(*it2)->getName());
      std::pair<std::vector<std::string>, std::string > p(tmp,tmp.front());
      ret[i]=p;
      std::vector<int> tmp2((*it).begin()+1,(*it).end());
      killProfileIds(tmp2);
    }
  changePflsRefsNamesGen(ret);
  return ret;
}

/*!
 * Removes duplicated localizations. Returns a map used to update references to removed 
 * localizations via changeLocsRefsNamesGen().
 * Equal localizations are found using whichAreEqualLocs().
 *  \param [in] eps - a precision used to compare real values of the localizations.
 *  \return std::vector< std::pair<std::vector<std::string>, std::string > > - 
 *          a sequence describing the performed replacements of localizations. Each element of
 *          this sequence is a pair whose
 *          - the first item is a vector of localization names replaced by the second item,
 *          - the second item is a localization name replacing every localization of the first item.
 */
std::vector< std::pair<std::vector<std::string>, std::string > > MEDFileFieldGlobsReal::zipLocsNames(double eps)
{
  std::vector< std::vector<int> > pseudoRet=whichAreEqualLocs(eps);
  std::vector< std::pair<std::vector<std::string>, std::string > > ret(pseudoRet.size());
  int i=0;
  for(std::vector< std::vector<int> >::const_iterator it=pseudoRet.begin();it!=pseudoRet.end();it++,i++)
    {
      std::vector< std::string > tmp((*it).size());
      int j=0;
      for(std::vector<int>::const_iterator it2=(*it).begin();it2!=(*it).end();it2++,j++)
        tmp[j]=std::string(getLocalizationFromId(*it2).getName());
      std::pair<std::vector<std::string>, std::string > p(tmp,tmp.front());
      ret[i]=p;
      std::vector<int> tmp2((*it).begin()+1,(*it).end());
      killLocalizationIds(tmp2);
    }
  changeLocsRefsNamesGen(ret);
  return ret;
}

/*!
 * Returns number of Gauss points per cell in a given localization.
 *  \param [in] locId - an id of the localization of interest.
 *  \return int - the number of the Gauss points per cell.
 */
int MEDFileFieldGlobsReal::getNbOfGaussPtPerCell(int locId) const
{
  return contentNotNull()->getNbOfGaussPtPerCell(locId);
}

/*!
 * Returns an id of a localization by its name.
 *  \param [in] loc - the localization name of interest.
 *  \return int - the id of the localization.
 *  \throw If there is no a localization named \a loc.
 */
int MEDFileFieldGlobsReal::getLocalizationId(const std::string& loc) const
{
  return contentNotNull()->getLocalizationId(loc);
}

/*!
 * Returns an id of a profile by its name.
 *  \param [in] loc - the profile name of interest.
 *  \return int - the id of the profile.
 *  \throw If there is no a profile named \a loc.
 */
int MEDFileFieldGlobsReal::getProfileId(const std::string& pfl) const
{
  return contentNotNull()->getProfileId(pfl);
}

/*!
 * Returns the name of the MED file.
 *  \return const std::string&  - the MED file name.
 */
std::string MEDFileFieldGlobsReal::getFileName() const
{
  return contentNotNull()->getFileName();
}

/*!
 * Returns a localization object by its name.
 *  \param [in] locName - the name of the localization of interest.
 *  \return const MEDFileFieldLoc& - the localization object having the name \a locName.
 *  \throw If there is no a localization named \a locName.
 */
const MEDFileFieldLoc& MEDFileFieldGlobsReal::getLocalization(const std::string& locName) const
{
  return contentNotNull()->getLocalization(locName);
}

/*!
 * Returns a localization object by its id.
 *  \param [in] locId - the id of the localization of interest.
 *  \return const MEDFileFieldLoc& - the localization object having the id \a locId.
 *  \throw If there is no a localization with id \a locId.
 */
const MEDFileFieldLoc& MEDFileFieldGlobsReal::getLocalizationFromId(int locId) const
{
  return contentNotNull()->getLocalizationFromId(locId);
}

/*!
 * Returns a profile array by its name.
 *  \param [in] pflName - the name of the profile of interest.
 *  \return const DataArrayInt * - the profile array having the name \a pflName.
 *  \throw If there is no a profile named \a pflName.
 */
const DataArrayInt *MEDFileFieldGlobsReal::getProfile(const std::string& pflName) const
{
  return contentNotNull()->getProfile(pflName);
}

/*!
 * Returns a profile array by its id.
 *  \param [in] pflId - the id of the profile of interest.
 *  \return const DataArrayInt * - the profile array having the id \a pflId.
 *  \throw If there is no a profile with id \a pflId.
 */
const DataArrayInt *MEDFileFieldGlobsReal::getProfileFromId(int pflId) const
{
  return contentNotNull()->getProfileFromId(pflId);
}

/*!
 * Returns a localization object, apt for modification, by its id.
 *  \param [in] locId - the id of the localization of interest.
 *  \return MEDFileFieldLoc& - a non-const reference to the localization object
 *          having the id \a locId.
 *  \throw If there is no a localization with id \a locId.
 */
MEDFileFieldLoc& MEDFileFieldGlobsReal::getLocalizationFromId(int locId)
{
  return contentNotNull()->getLocalizationFromId(locId);
}

/*!
 * Returns a localization object, apt for modification, by its name.
 *  \param [in] locName - the name of the localization of interest.
 *  \return MEDFileFieldLoc& - a non-const reference to the localization object
 *          having the name \a locName.
 *  \throw If there is no a localization named \a locName.
 */
MEDFileFieldLoc& MEDFileFieldGlobsReal::getLocalization(const std::string& locName)
{
  return contentNotNull()->getLocalization(locName);
}

/*!
 * Returns a profile array, apt for modification, by its name.
 *  \param [in] pflName - the name of the profile of interest.
 *  \return DataArrayInt * - Borrowed reference - a non-const pointer to the profile array having the name \a pflName.
 *  \throw If there is no a profile named \a pflName.
 */
DataArrayInt *MEDFileFieldGlobsReal::getProfile(const std::string& pflName)
{
  return contentNotNull()->getProfile(pflName);
}

/*!
 * Returns a profile array, apt for modification, by its id.
 *  \param [in] pflId - the id of the profile of interest.
 *  \return DataArrayInt * - Borrowed reference - a non-const pointer to the profile array having the id \a pflId.
 *  \throw If there is no a profile with id \a pflId.
 */
DataArrayInt *MEDFileFieldGlobsReal::getProfileFromId(int pflId)
{
  return contentNotNull()->getProfileFromId(pflId);
}

/*!
 * Removes profiles given by their ids. No data is updated to track this removal.
 *  \param [in] pflIds - a sequence of ids of the profiles to remove.
 */
void MEDFileFieldGlobsReal::killProfileIds(const std::vector<int>& pflIds)
{
  contentNotNull()->killProfileIds(pflIds);
}

/*!
 * Removes localizations given by their ids. No data is updated to track this removal.
 *  \param [in] locIds - a sequence of ids of the localizations to remove.
 */
void MEDFileFieldGlobsReal::killLocalizationIds(const std::vector<int>& locIds)
{
  contentNotNull()->killLocalizationIds(locIds);
}

/*!
 * Stores a profile array.
 *  \param [in] pfl - the profile array to store.
 *  \throw If the name of \a pfl is empty.
 *  \throw If a profile with the same name as that of \a pfl already exists but contains
 *         different ids.
 */
void MEDFileFieldGlobsReal::appendProfile(DataArrayInt *pfl)
{
  contentNotNull()->appendProfile(pfl);
}

/*!
 * Adds a new localization of Gauss points.
 *  \param [in] locName - the name of the new localization.
 *  \param [in] geoType - a geometrical type of the reference cell.
 *  \param [in] refCoo - coordinates of points of the reference cell. Size of this vector
 *         must be \c nbOfNodesPerCell * \c dimOfType.
 *  \param [in] gsCoo - coordinates of Gauss points on the reference cell. Size of this vector
 *         must be  _wg_.size() * \c dimOfType.
 *  \param [in] w - the weights of Gauss points.
 *  \throw If \a locName is empty.
 *  \throw If a localization with the name \a locName already exists but is
 *         different form the new one.
 */
void MEDFileFieldGlobsReal::appendLoc(const std::string& locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w)
{
  contentNotNull()->appendLoc(locName,geoType,refCoo,gsCoo,w);
}

MEDFileFieldGlobs *MEDFileFieldGlobsReal::contentNotNull()
{
  MEDFileFieldGlobs *g(_globals);
  if(!g)
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobsReal::contentNotNull : no content in not const !");
  return g;
}

const MEDFileFieldGlobs *MEDFileFieldGlobsReal::contentNotNull() const
{
  const MEDFileFieldGlobs *g(_globals);
  if(!g)
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobsReal::contentNotNull : no content in const !");
  return g;
}

//= MEDFileFieldNameScope

MEDFileFieldNameScope::MEDFileFieldNameScope()
{
}

MEDFileFieldNameScope::MEDFileFieldNameScope(const std::string& fieldName, const std::string& meshName):_name(fieldName),_mesh_name(meshName)
{
}

/*!
 * Returns the name of \a this field.
 *  \return std::string - a string containing the field name.
 */
std::string MEDFileFieldNameScope::getName() const
{
  return _name;
}

/*!
 * Sets name of \a this field
 *  \param [in] name - the new field name.
 */
void MEDFileFieldNameScope::setName(const std::string& fieldName)
{
  _name=fieldName;
}

std::string MEDFileFieldNameScope::getDtUnit() const
{
  return _dt_unit;
}

void MEDFileFieldNameScope::setDtUnit(const std::string& dtUnit)
{
  _dt_unit=dtUnit;
}

void MEDFileFieldNameScope::copyNameScope(const MEDFileFieldNameScope& other)
{
  _name=other._name;
  _mesh_name=other._mesh_name;
  _dt_unit=other._dt_unit;
}

/*!
 * Returns the mesh name.
 *  \return std::string - a string holding the mesh name.
 *  \throw If \c _field_per_mesh.empty()
 */
std::string MEDFileFieldNameScope::getMeshName() const
{
  return _mesh_name;
}

void MEDFileFieldNameScope::setMeshName(const std::string& meshName)
{
  _mesh_name=meshName;
}

