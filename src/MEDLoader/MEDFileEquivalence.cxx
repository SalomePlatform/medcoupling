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
// Author : Anthony Geay (EDF R&D)

#include "MEDFileEquivalence.hxx"
#include "MEDFileSafeCaller.txx"
#include "MEDCouplingMemArray.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileMesh.hxx"
#include "InterpKernelAutoPtr.hxx"

extern med_geometry_type typmai[MED_N_CELL_FIXED_GEO];
extern INTERP_KERNEL::NormalizedCellType typmai2[MED_N_CELL_FIXED_GEO];
extern med_geometry_type typmai3[34];
extern med_geometry_type typmainoeud[1];

using namespace MEDCoupling;

MEDFileEquivalencePair *MEDFileEquivalencePair::Load(MEDFileEquivalences *father, med_idt fid, const std::string& name, const std::string &desc)
{
  if(!father)
    throw INTERP_KERNEL::Exception("MEDFileEquivalencePair::Load : father is NULL ! Should not !");
  MCAuto<MEDFileEquivalencePair> ret(new MEDFileEquivalencePair(father,name,desc));
  ret->load(fid);
  return ret.retn();
}

void MEDFileEquivalencePair::writeLL(med_idt fid) const
{
  std::string meshName(getFather()->getMeshName());
  INTERP_KERNEL::AutoPtr<char> meshName2(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> name(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> desc(MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE));
  MEDLoaderBase::safeStrCpy(meshName.c_str(),MED_NAME_SIZE,meshName2,getFather()->getMesh()->getTooLongStrPolicy());
  MEDLoaderBase::safeStrCpy(_name.c_str(),MED_NAME_SIZE,name,getFather()->getMesh()->getTooLongStrPolicy());
  MEDLoaderBase::safeStrCpy(_description.c_str(),MED_COMMENT_SIZE,desc,getFather()->getMesh()->getTooLongStrPolicy());
  MEDFILESAFECALLERWR0(MEDequivalenceCr,(fid,meshName2,name,desc));
  const MEDFileEquivalenceCell *cell(_cell);
  if(cell)
    cell->writeLL(fid);
  const MEDFileEquivalenceNode *node(_node);
  if(node)
    node->writeLL(fid);
}

const MEDFileMesh *MEDFileEquivalencePair::getMesh() const
{
  return getFather()->getMesh();
}

MEDFileMesh *MEDFileEquivalencePair::getMesh()
{
  return getFather()->getMesh();
}

MEDFileEquivalencePair *MEDFileEquivalencePair::deepCopy(MEDFileEquivalences *father) const
{
  MCAuto<MEDFileEquivalencePair> ret(new MEDFileEquivalencePair(father,_name,_description));
  const MEDFileEquivalenceCell *cell(_cell);
  if(cell)
    ret->_cell=cell->deepCopy(const_cast<MEDFileEquivalencePair *>(this));
  const MEDFileEquivalenceNode *node(_node);
  if(node)
    ret->_node=node->deepCopy(const_cast<MEDFileEquivalencePair *>(this));
  return ret.retn();
}

bool MEDFileEquivalencePair::isEqual(const MEDFileEquivalencePair *other, std::string& what) const
{
  if(_name!=other->_name)
    {
      std::ostringstream oss; oss << "Names differs : " << _name << " != " << other->_name << " !";
      what=oss.str();
      return false;
    }
  if(_description!=other->_description)
    {
       std::ostringstream oss; oss << "Description differs : " << _description << " != " << other->_description << " !";
       what=oss.str();
       return false;
    }
  const MEDFileEquivalenceCell *c1(_cell),*c2(other->_cell);
  if((c1 && !c2) || (!c1 && c2))
    {
      std::ostringstream oss; oss << "Cell def of Equiv " << _name << " are defined for this and not for other (or reversely) !";
      what=oss.str();
      return false;
    }
  if(c1 && c2)
    if(!c1->isEqual(c2,what))
      return false;
  const MEDFileEquivalenceNode *n1(_node),*n2(other->_node);
  if((n1 && !n2) || (!n1 && n2))
    {
      std::ostringstream oss; oss << "Node def of Equiv " << _name << " are defined for this and not for other (or reversely) !";
      what=oss.str();
      return false;
    }
  if(n1 && n2)
    if(!n1->isEqual(n2,what))
      return false;
  return true;
}

void MEDFileEquivalencePair::getRepr(std::ostream& oss) const
{
  const MEDFileEquivalenceNode *node(_node);
  const MEDFileEquivalenceCell *cell(_cell);
  oss << std::endl << "  name of equivalence : " << _name << std::endl;
  oss << "  description of equivalence : " << _description << std::endl;
  oss << "  Node : ";
  if(!node)
    oss << "None" << std::endl;
  else
    node->getRepr(oss);
  oss << "  Cell : ";
  if(!cell)
    oss << "None" << std::endl;
  else
    cell->getRepr(oss);
}

MEDFileEquivalencePair *MEDFileEquivalencePair::New(MEDFileEquivalences *father, const std::string& name)
{
  return new MEDFileEquivalencePair(father,name,std::string());
}

std::vector<const BigMemoryObject *> MEDFileEquivalencePair::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(2);
  ret[0]=_cell; ret[1]=_node;
  return ret;
}

void MEDFileEquivalencePair::setArray(int meshDimRelToMaxExt, DataArrayInt *da)
{
  if(meshDimRelToMaxExt>1)
    throw INTERP_KERNEL::Exception("MEDFileEquivalencePair::setArray : meshDimRelToMaxExt must be in [1,0,-1,-2,-3] at most !");
  if(meshDimRelToMaxExt==1)
    {
      MEDFileEquivalenceNode *node(_node);
      if(!node)
        {
          _node=new MEDFileEquivalenceNode(this,0);
          node=_node;
        }
      node->setArray(da);
    }
  else
    {
      MEDFileEquivalenceCell *cell(_cell);
      if(!cell)
        {
          _cell=new MEDFileEquivalenceCell(this);
          cell=_cell;
        }
      cell->setArray(meshDimRelToMaxExt,da);
    }
}

/*!
 * The returned pointer is a borrowed pointer.
 */
MEDFileEquivalenceCell *MEDFileEquivalencePair::initCell()
{
  _cell=new MEDFileEquivalenceCell(this);
  return _cell;
}

/*!
 * The returned pointer is a borrowed pointer.
 */
MEDFileEquivalenceNode *MEDFileEquivalencePair::initNode()
{
  _node=new MEDFileEquivalenceNode(this,0);
  return _node;
}

std::size_t MEDFileEquivalencePair::getHeapMemorySizeWithoutChildren() const
{
  return 0;
}

void MEDFileEquivalencePair::load(med_idt fid)
{
  std::string meshName(_father->getMeshName());
  int dt,it;
  _father->getDtIt(dt,it);
  med_int ncor;
  MEDFILESAFECALLERRD0(MEDequivalenceCorrespondenceSize,(fid,meshName.c_str(),_name.c_str(),dt,it,MED_NODE,MED_NONE,&ncor));
  if(ncor>0)
    {
      MCAuto<DataArrayInt> da(DataArrayInt::New());
      da->alloc(ncor*2);
      MEDFILESAFECALLERRD0(MEDequivalenceCorrespondenceRd,(fid,meshName.c_str(),_name.c_str(),dt,it,MED_NODE,MED_NONE,da->getPointer()));
      da->applyLin(1,-1);
      da->rearrange(2);
      MCAuto<MEDFileEquivalenceNode> node(new MEDFileEquivalenceNode(this,da));
      _node=node;
    }
  _cell=MEDFileEquivalenceCell::Load(fid,this);
}

std::vector<const BigMemoryObject *> MEDFileEquivalences::getDirectChildrenWithNull() const
{
  std::size_t sz(_equ.size());
  std::vector<const BigMemoryObject *> ret(sz);
  for(std::size_t i=0;i<sz;i++)
    ret[i]=_equ[i];
  return ret;
}

std::size_t MEDFileEquivalences::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(MEDFileEquivalences)+_equ.capacity()*sizeof(MEDFileEquivalencePair);
}

void MEDFileEquivalences::getDtIt(int &dt, int &it) const
{
  dt=_owner->getIteration(); it=_owner->getOrder();
}

std::string MEDFileEquivalences::getMeshName() const
{
  return _owner->getName();
}

void MEDFileEquivalences::pushEquivalence(MEDFileEquivalencePair *elt)
{
  MCAuto<MEDFileEquivalencePair> elta(elt);
  if(elt)
    elt->incrRef();
  _equ.push_back(elta);
}

MEDFileEquivalencePair *MEDFileEquivalences::getEquivalence(int i)
{
  int sz(size());
  if(i<0 || i>=sz)
    {
      std::ostringstream oss; oss << "MEDFileEquivalences::getEquivalence : invalid id ! Must be in [0," << sz << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return _equ[i];
}

MEDFileEquivalencePair *MEDFileEquivalences::getEquivalenceWithName(const std::string& name)
{
  for(std::vector< MCAuto<MEDFileEquivalencePair> >::iterator it=_equ.begin();it!=_equ.end();it++)
    {
      MEDFileEquivalencePair *elt(*it);
      if(elt)
        {
          if(elt->getName()==name)
            return elt;
        }
    }
  std::ostringstream oss; oss << "MEDFileEquivalences::getEquivalenceWithName : no equivalence with name \"" << name << "\" ! Must be in [ ";
  std::vector<std::string> eqs(getEquivalenceNames());
  std::copy(eqs.begin(),eqs.end(),std::ostream_iterator<std::string>(oss,", "));
  oss << "] !";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

int MEDFileEquivalences::size() const
{
  return _equ.size();
}

std::vector<std::string> MEDFileEquivalences::getEquivalenceNames() const
{
  std::vector<std::string> ret;
  for(std::vector< MCAuto<MEDFileEquivalencePair> >::const_iterator it=_equ.begin();it!=_equ.end();it++)
    {
      const MEDFileEquivalencePair *elt(*it);
      if(elt)
        {
          ret.push_back(elt->getName());
        }
    }
  return ret;
}

MEDFileEquivalencePair *MEDFileEquivalences::appendEmptyEquivalenceWithName(const std::string& name)
{
  MCAuto<MEDFileEquivalencePair> elt(MEDFileEquivalencePair::New(this,name));
  _equ.push_back(elt);
  return elt;
}

MEDFileEquivalences *MEDFileEquivalences::deepCopy(MEDFileMesh *owner) const
{
  MCAuto<MEDFileEquivalences> ret(new MEDFileEquivalences(owner));
  ret->deepCpyFrom(*this);
  return ret.retn();
}

bool MEDFileEquivalences::isEqual(const MEDFileEquivalences *other, std::string& what) const
{
  std::size_t sz(_equ.size());
  if(sz!=other->_equ.size())
    {
      what="Equivalences differs : not same number !";
      return false;
    }
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDFileEquivalencePair *thisp(_equ[i]),*otherp(other->_equ[i]);
      if(!thisp && !otherp)
        continue;
      if(thisp && otherp)
        {
          if(!thisp->isEqual(otherp,what))
            {
              std::ostringstream oss; oss << "At Eq #" << i << " there is a difference !";
              what=oss.str()+what;
              return false;
            }
        }
      else
        {
          std::ostringstream oss; oss << "At Eq #" << i << " defined in this not is other (or reversely) !";
          what=oss.str()+what;
          return false;
        }
    }
  return true;
}

void MEDFileEquivalences::getRepr(std::ostream& oss) const
{
  std::size_t ii(0);
  for(std::vector< MCAuto<MEDFileEquivalencePair> >::const_iterator it=_equ.begin();it!=_equ.end();it++,ii++)
    {
      const MEDFileEquivalencePair *elt(*it);
      oss << "Equivalence #" << ii << " : " ;
      if(elt)
        elt->getRepr(oss);
      else
        oss << "None" << std::endl;
    }
}

void MEDFileEquivalences::killEquivalenceWithName(const std::string& name)
{
  std::vector< MCAuto<MEDFileEquivalencePair> >::iterator it(_equ.begin());
  for(;it!=_equ.end();it++)
    {
      const MEDFileEquivalencePair *elt(*it);
      if(elt && elt->getName()==name)
        break;
    }
  if(it==_equ.end())
    {
      std::ostringstream oss; oss << "MEDFileEquivalences::killEquivalenceWithName : Equivalence with name \"" << name << "\" not found !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  _equ.erase(it);
}

void MEDFileEquivalences::killEquivalenceAt(int i)
{
  int sz(size());
  if(i<0 || i>=sz)
    {
      std::ostringstream oss; oss << "MEDFileEquivalences::killEquivalenceAt : Id must be in [0," << sz << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::vector< MCAuto<MEDFileEquivalencePair> >::iterator it(_equ.begin());
  for(int j=0;j<i;it++,j++);
  _equ.erase(it);
}

void MEDFileEquivalences::clear()
{
  _equ.clear();
}

void MEDFileEquivalences::writeLL(med_idt fid) const
{
  for(std::vector< MCAuto<MEDFileEquivalencePair> >::const_iterator it=_equ.begin();it!=_equ.end();it++)
    {
      const MEDFileEquivalencePair *elt(*it);
      if(elt)
        elt->writeLL(fid);
    }
}

int MEDFileEquivalences::PresenceOfEquivalences(med_idt fid, const std::string& meshName)
{
  med_int nequ(MEDnEquivalence(fid,meshName.c_str()));
  return nequ;
}

MEDFileEquivalences *MEDFileEquivalences::Load(med_idt fid, int nbOfEq, MEDFileMesh *owner)
{
  MCAuto<MEDFileEquivalences> ret(new MEDFileEquivalences(owner));
  if(!owner)
    throw INTERP_KERNEL::Exception("MEDFileEquivalences::Load : owner is NULL !");
  std::string meshName(owner->getName());
  for(int i=0;i<nbOfEq;i++)
    {
      INTERP_KERNEL::AutoPtr<char> equ(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
      INTERP_KERNEL::AutoPtr<char> desc(MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE));
      int nstep,nocstpncor;
      MEDFILESAFECALLERRD0(MEDequivalenceInfo,(fid,meshName.c_str(),i+1,equ,desc,&nstep,&nocstpncor));
      std::string eqName(MEDLoaderBase::buildStringFromFortran(equ,MED_NAME_SIZE)),eqDescName(MEDLoaderBase::buildStringFromFortran(desc,MED_COMMENT_SIZE));
      MCAuto<MEDFileEquivalencePair> eqv(MEDFileEquivalencePair::Load(ret,fid,eqName,eqDescName));
      ret->pushEquivalence(eqv);
    }
  return ret.retn();
}

void MEDFileEquivalences::CheckDataArray(const DataArrayInt *data)
{
  if(!data)
    return;
  data->checkAllocated();
  if(data->getNumberOfComponents()!=2)
    {
      std::ostringstream oss; oss << "MEDFileEquivalences::CheckDataArray : Input DataArray must have 2 components !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void MEDFileEquivalences::deepCpyFrom(const MEDFileEquivalences& other)
{
  for(std::vector< MCAuto<MEDFileEquivalencePair> >::const_iterator it=other._equ.begin();it!=other._equ.end();it++)
    {
      const MEDFileEquivalencePair *elt(*it);
      MCAuto<MEDFileEquivalencePair> eltCpy;
      if(elt)
        {
          eltCpy=elt->deepCopy(this);
        }
      _equ.push_back(eltCpy);
    }
}

MEDFileEquivalenceBase::MEDFileEquivalenceBase(MEDFileEquivalencePair *father):_father(father)
{
}

MEDFileEquivalenceData::MEDFileEquivalenceData(MEDFileEquivalencePair *owner, DataArrayInt *data):MEDFileEquivalenceBase(owner),_data(data)
{
  if(data)
    data->incrRef();
}

void MEDFileEquivalenceData::setArray(DataArrayInt *data)
{
  MEDFileEquivalences::CheckDataArray(data);
  _data=data;
  if(data)
    data->incrRef();
}

std::vector<const BigMemoryObject *> MEDFileEquivalenceData::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(1);
  ret[0]=_data;
  return ret;
}

bool MEDFileEquivalenceData::isEqual(const MEDFileEquivalenceData *other, std::string& what) const
{
  const DataArrayInt *d1(_data),*d2(other->_data);
  if((!d1 && d2) || (d1 && !d2))
    {
      what="Data array is defined in this not in other (or reversely) !";
      return false;
    }
  if(d1 && d2)
    {
      if(!d1->isEqualIfNotWhy(*d2,what))
        return false;
    }
  return true;
}

void MEDFileEquivalenceData::writeAdvanced(med_idt fid, med_entity_type medtype, med_geometry_type medgt) const
{
  
  const DataArrayInt *da(getArray());
  if(!da)
    return ;
  MEDFileEquivalences::CheckDataArray(da);
  const MEDFileMesh *mesh(getFather()->getMesh());
  int dt,it;
  mesh->getTime(dt,it);
  std::string meshName(mesh->getName());
  std::string equName(getFather()->getName());
  INTERP_KERNEL::AutoPtr<char> meshName2(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> name(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  MEDLoaderBase::safeStrCpy(meshName.c_str(),MED_NAME_SIZE,meshName2,getFather()->getMesh()->getTooLongStrPolicy());
  MEDLoaderBase::safeStrCpy(equName.c_str(),MED_NAME_SIZE,name,getFather()->getMesh()->getTooLongStrPolicy());
  MCAuto<DataArrayInt> da2(da->deepCopy()); da2->rearrange(1); da2->applyLin(1,1); da2->rearrange(2);
  MEDFILESAFECALLERWR0(MEDequivalenceCorrespondenceWr,(fid,meshName2,name,dt,it,medtype,medgt,da2->getNumberOfTuples(),da2->begin()));
}

std::size_t MEDFileEquivalenceCellType::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(MEDFileEquivalenceCellType);
}

MEDFileEquivalenceCellType *MEDFileEquivalenceCellType::deepCopy(MEDFileEquivalencePair *owner) const
{
  MCAuto<DataArrayInt> da;
  if(getArray())
    da=getArray()->deepCopy();
  return new MEDFileEquivalenceCellType(owner,_type,da);
}

bool MEDFileEquivalenceCellType::isEqual(const MEDFileEquivalenceCellType *other, std::string& what) const
{
  if(_type!=other->_type)
    {
      what="Geo types differs !";
      return false;
    }
  return MEDFileEquivalenceData::isEqual(other,what);
}

void MEDFileEquivalenceCellType::getRepr(std::ostream& oss) const
{
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(_type));
  const DataArrayInt *da(getArray());
  oss << cm.getRepr() << ":";
  if(da)
    oss << da->getNumberOfTuples() << " tuples";
  else
    oss << "no dataarray";
  oss << ",";
}

void MEDFileEquivalenceCellType::writeLL(med_idt fid) const
{
  writeAdvanced(fid,MED_CELL,typmai3[_type]);
}

std::vector<const BigMemoryObject *> MEDFileEquivalenceCell::getDirectChildrenWithNull() const
{
  std::size_t sz(_types.size());
  std::vector<const BigMemoryObject *> ret(sz);
  for(std::size_t i=0;i<sz;i++)
    ret[i]=_types[i];
  return ret;
}

std::size_t MEDFileEquivalenceCell::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(MEDFileEquivalenceCell)+_types.capacity()*sizeof(MEDFileEquivalenceCellType);
}

MEDFileEquivalenceCell *MEDFileEquivalenceCell::Load(med_idt fid, MEDFileEquivalencePair *owner)
{
  MCAuto<MEDFileEquivalenceCell> ret(new MEDFileEquivalenceCell(owner));
  ret->load(fid);
  if(ret->size()>0)
    return ret.retn();
  else
    return 0;
}

void MEDFileEquivalenceCell::writeLL(med_idt fid) const
{
  for(std::vector< MCAuto<MEDFileEquivalenceCellType> >::const_iterator it=_types.begin();it!=_types.end();it++)
    {
      const MEDFileEquivalenceCellType *ct(*it);
      if(ct)
        ct->writeLL(fid);
    }
}

MEDFileEquivalenceCell *MEDFileEquivalenceCell::deepCopy(MEDFileEquivalencePair *owner) const
{
  MCAuto<MEDFileEquivalenceCell> ret(new MEDFileEquivalenceCell(owner));
  for(std::vector< MCAuto<MEDFileEquivalenceCellType> >::const_iterator it=_types.begin();it!=_types.end();it++)
    {
      const MEDFileEquivalenceCellType *elt(*it);
      MCAuto<MEDFileEquivalenceCellType> eltCpy;
      if(elt)
        eltCpy=elt->deepCopy(owner);
      ret->_types.push_back(eltCpy);
    }
  return ret.retn();
}

bool MEDFileEquivalenceCell::isEqual(const MEDFileEquivalenceCell *other, std::string& what) const
{
  std::size_t sz(_types.size());
  if(sz!=other->_types.size())
    {
      std::ostringstream oss; oss << "Nb of geo types differs : " << sz << " != " << other->_types.size();
      what=oss.str();
      return false;
    }
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDFileEquivalenceCellType *ct1(_types[i]),*ct2(other->_types[i]);
      if((ct1 && !ct2) || (!ct1 && ct2))
        {
          std::ostringstream oss; oss << "At gt #" << i << " this is defined not other (or reversely !)";
          what=oss.str();
          return false;
        }
      if(ct1 && ct2)
        {
          if(!ct1->isEqual(ct2,what))
            {
              std::ostringstream oss; oss << "At gt #" << i << " of Eq " << getFather()->getName() << " it differs !";
              what=oss.str()+what;
              return false;
            }
        }
    }
  return true;
}

void MEDFileEquivalenceCell::getRepr(std::ostream& oss) const
{
  for(std::vector< MCAuto<MEDFileEquivalenceCellType> >::const_iterator it=_types.begin();it!=_types.end();it++)
    {
      const MEDFileEquivalenceCellType *elt(*it);
      if(elt)
        elt->getRepr(oss);
    }
}

DataArrayInt *MEDFileEquivalenceCell::getArray(INTERP_KERNEL::NormalizedCellType type)
{
  for(std::vector< MCAuto<MEDFileEquivalenceCellType> >::iterator it=_types.begin();it!=_types.end();it++)
    {
      MEDFileEquivalenceCellType *elt(*it);
      if(elt && elt->getType()==type)
        return elt->getArray();
    }
  std::ostringstream oss; oss << "MEDFileEquivalenceCell::getArray : In Equivalence \"" << getFather()->getName() << "\" the geotype " << type << " is not available !";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

void MEDFileEquivalenceCell::setArray(int meshDimRelToMax, DataArrayInt *da)
{
  if(!da)
    return ;
  MEDFileEquivalences::CheckDataArray(da);
  MEDFileMesh *mm(getMesh());
  int totalNbOfCells(mm->getNumberOfCellsAtLevel(meshDimRelToMax));
  //
  MCAuto<DataArrayInt> tmp(da->deepCopy()); tmp->rearrange(1);
  int maxv,minv;
  tmp->getMinMaxValues(minv,maxv);
  if((minv<0 || minv>=totalNbOfCells) || (maxv<0 || maxv>=totalNbOfCells))
    {
      std::ostringstream oss; oss << "MEDFileEquivalenceCell::setArray : Input 2 component DataArray has incorrect values ! all values must be in [0," << totalNbOfCells << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  //
  std::vector<INTERP_KERNEL::NormalizedCellType> gts(mm->getGeoTypesAtLevel(meshDimRelToMax));
  int startId(0),endId;
  std::vector<int> compS(1,0);
  for(std::vector<INTERP_KERNEL::NormalizedCellType>::const_iterator it=gts.begin();it!=gts.end();it++)
    {
      endId=startId+mm->getNumberOfCellsWithType(*it);
      MCAuto<DataArrayInt> da0(da->keepSelectedComponents(compS));
      MCAuto<DataArrayInt> ids(da0->findIdsInRange(startId,endId));
      MCAuto<DataArrayInt> da1(da->selectByTupleIdSafe(ids->begin(),ids->end()));
      da1->applyLin(1,-startId);
      setArrayForType(*it,da1);
      startId=endId;
    }
}

void MEDFileEquivalenceCell::setArrayForType(INTERP_KERNEL::NormalizedCellType type, DataArrayInt *da)
{
  for(std::vector< MCAuto<MEDFileEquivalenceCellType> >::iterator it=_types.begin();it!=_types.end();it++)
    {
      MEDFileEquivalenceCellType *elt(*it);
      if(elt && elt->getType()==type)
        {
          elt->setArray(da);
          return ;
        }
    }
  MCAuto<MEDFileEquivalenceCellType> newElt(new MEDFileEquivalenceCellType(getFather(),type,da));
  _types.push_back(newElt);
}

std::vector<INTERP_KERNEL::NormalizedCellType> MEDFileEquivalenceCell::getTypes() const
{
  std::vector<INTERP_KERNEL::NormalizedCellType> ret;
  for(std::vector< MCAuto<MEDFileEquivalenceCellType> >::const_iterator it=_types.begin();it!=_types.end();it++)
    {
      const MEDFileEquivalenceCellType *elt(*it);
      if(elt)
        ret.push_back(elt->getType());
    }
  return ret;
}

void MEDFileEquivalenceCell::load(med_idt fid)
{
  std::string meshName(getFather()->getFather()->getMeshName()),name(getName());
  int dt,it;
  getFather()->getFather()->getDtIt(dt,it);
  for(int i=0;i<MED_N_CELL_FIXED_GEO;i++)
    {
      med_int ncor;
      MEDFILESAFECALLERRD0(MEDequivalenceCorrespondenceSize,(fid,meshName.c_str(),name.c_str(),dt,it,MED_CELL,typmai[i],&ncor));
      if(ncor>0)
        {
          MCAuto<DataArrayInt> da(DataArrayInt::New());
          da->alloc(ncor*2);
          MEDFILESAFECALLERRD0(MEDequivalenceCorrespondenceRd,(fid,meshName.c_str(),name.c_str(),dt,it,MED_CELL,typmai[i],da->getPointer()));
          da->applyLin(1,-1);
          da->rearrange(2);
          MCAuto<MEDFileEquivalenceCellType> ct(new MEDFileEquivalenceCellType(getFather(),typmai2[i],da));
          _types.push_back(ct);
        }
    }
}

std::size_t MEDFileEquivalenceNode::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(MEDFileEquivalenceNode);
}

void MEDFileEquivalenceNode::writeLL(med_idt fid) const
{
  writeAdvanced(fid,MED_NODE,MED_NONE);
}

MEDFileEquivalenceNode *MEDFileEquivalenceNode::deepCopy(MEDFileEquivalencePair *owner) const
{
  MCAuto<DataArrayInt> da;
  if(getArray())
    da=getArray()->deepCopy();
  MCAuto<MEDFileEquivalenceNode> ret(new MEDFileEquivalenceNode(owner,da));
  return ret.retn();
}

bool MEDFileEquivalenceNode::isEqual(const MEDFileEquivalenceNode *other, std::string& what) const
{
  return MEDFileEquivalenceData::isEqual(other,what);
}

void MEDFileEquivalenceNode::getRepr(std::ostream& oss) const
{
  const DataArrayInt *da(getArray());
  if(!da)
    oss << " No dataarray defined !" << std::endl;
  else
    oss << da->getNumberOfTuples() << " tuples in node equivalence." << std::endl;
}
