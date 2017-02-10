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

#include "MEDFileStructureElement.hxx"
#include "MEDFileMeshSupport.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileMeshLL.hxx"
#include "MEDFileSafeCaller.txx"

#include "InterpKernelAutoPtr.hxx"

using namespace MEDCoupling;


std::string MEDFileSEHolder::getModelName() const
{
  return _father->getName();
}

std::string MEDFileSEHolder::getName() const
{
  return _name;
}

void MEDFileSEHolder::setName(const std::string& name)
{
  _name=name;
}

std::size_t MEDFileSEHolder::getHeapMemorySizeLoc() const
{
  return _name.capacity();
}

////////////////////

MEDFileSEConstAtt *MEDFileSEConstAtt::New(med_idt fid, MEDFileStructureElement *father, int idCstAtt, const MEDFileUMesh *mesh)
{
  return new MEDFileSEConstAtt(fid,father,idCstAtt,mesh);
}

MEDFileSEConstAtt::MEDFileSEConstAtt(med_idt fid, MEDFileStructureElement *father, int idCstAtt, const MEDFileUMesh *mesh):MEDFileSEHolder(father)
{
  std::string modelName(getModelName());
  INTERP_KERNEL::AutoPtr<char> constattname(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE)),profilename(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  med_attribute_type constatttype;
  int nbCompo;
  med_entity_type met;
  int pflSz;
  MEDFILESAFECALLERRD0(MEDstructElementConstAttInfo,(fid,modelName.c_str(),idCstAtt+1,constattname,&constatttype,&nbCompo,&met,profilename,&pflSz));
  std::string name(MEDLoaderBase::buildStringFromFortran(constattname,MED_NAME_SIZE));
  setName(name);
  setProfile(MEDLoaderBase::buildStringFromFortran(profilename,MED_NAME_SIZE));
  _tof=MEDFileMesh::ConvertFromMEDFileEntity(met);
  //
  _val=MEDFileStructureElement::BuildFrom(constatttype);
  nbCompo=MEDFileStructureElement::EffectiveNbCompo(constatttype,nbCompo);
  if(pflSz==0 && getProfile().empty())
    {
      switch(met)
        {
        case MED_CELL:
          {
            std::vector<INTERP_KERNEL::NormalizedCellType> gt(mesh->getAllGeoTypes());
            if(gt.size()!=1)
              throw INTERP_KERNEL::Exception("MEDFileSEConstAtt constr : only one cell type expected !");
            pflSz=mesh->getNumberOfCellsWithType(gt[0]);
            break;
          }
        case MED_NODE:
          {
            pflSz=mesh->getNumberOfNodes();
            break;
          }
        default:
          throw INTERP_KERNEL::Exception("MEDFileSEConstAtt cstr : not recognized entity type !");
        }
    }
  if(constatttype==MED_ATT_NAME)
    pflSz++;
  _val->alloc(pflSz,nbCompo);
  MEDFILESAFECALLERRD0(MEDstructElementConstAttRd,(fid,modelName.c_str(),name.c_str(),_val->getVoidStarPointer()));
  if(constatttype==MED_ATT_NAME)
    { pflSz--; _val->reAlloc(pflSz); }
}

std::vector<const BigMemoryObject *> MEDFileSEConstAtt::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  ret.push_back(_val);
  return ret;
}

std::size_t MEDFileSEConstAtt::getHeapMemorySizeWithoutChildren() const
{
  return getHeapMemorySizeLoc()+_pfl.capacity();
}

void MEDFileSEConstAtt::writeLL(med_idt fid) const
{
}

void MEDFileSEConstAtt::setProfile(const std::string& name)
{
  _pfl=name;
}

std::string MEDFileSEConstAtt::getProfile() const
{
  return _pfl;
}

////////////////////

MEDFileSEVarAtt *MEDFileSEVarAtt::New(med_idt fid, MEDFileStructureElement *father, int idVarAtt)
{
  return new MEDFileSEVarAtt(fid,father,idVarAtt);
}

MEDFileSEVarAtt::MEDFileSEVarAtt(med_idt fid, MEDFileStructureElement *father, int idVarAtt):MEDFileSEHolder(father)
{
  std::string modelName(getModelName());
  INTERP_KERNEL::AutoPtr<char> varattname(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE)),profilename(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  med_attribute_type varatttype;
  {
    int pflSz;
    MEDFILESAFECALLERRD0(MEDstructElementVarAttInfo,(fid,modelName.c_str(),idVarAtt+1,varattname,&varatttype,&_nb_compo));
  }
  setName(MEDLoaderBase::buildStringFromFortran(varattname,MED_NAME_SIZE));
  _gen=MEDFileStructureElement::BuildFrom(varatttype);
  _gen->alloc(0,1);
}

std::vector<const BigMemoryObject *> MEDFileSEVarAtt::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

std::size_t MEDFileSEVarAtt::getHeapMemorySizeWithoutChildren() const
{
  return getHeapMemorySizeLoc();
}

void MEDFileSEVarAtt::writeLL(med_idt fid) const
{
}

////////////////////

MEDFileStructureElement *MEDFileStructureElement::New(med_idt fid, int idSE, const MEDFileMeshSupports *ms)
{
  return new MEDFileStructureElement(fid,idSE,ms);
}

MEDFileStructureElement::MEDFileStructureElement(med_idt fid, int idSE, const MEDFileMeshSupports *ms)
{
  INTERP_KERNEL::AutoPtr<char> modelName(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE)),supportMeshName(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  med_geometry_type sgeoType;
  med_entity_type entiyType;
  int nConsAttr(0),nVarAttr(0);
  {
    med_bool anyPfl;
    int nnode(0),ncell(0);
    MEDFILESAFECALLERRD0(MEDstructElementInfo,(fid,idSE+1,modelName,&_id_type,&_dim,supportMeshName,&entiyType,&nnode,&ncell,&sgeoType,&nConsAttr,&anyPfl,&nVarAttr));
  }
  _name=MEDLoaderBase::buildStringFromFortran(modelName,MED_NAME_SIZE);
  _sup_mesh_name=MEDLoaderBase::buildStringFromFortran(supportMeshName,MED_NAME_SIZE);
  _geo_type=MEDFileMesh::ConvertFromMEDFileGeoType(sgeoType);
  _tof=MEDFileMesh::ConvertFromMEDFileEntity(entiyType);
  _cst_att.resize(nConsAttr);
  for(int i=0;i<nConsAttr;i++)
    _cst_att[i]=MEDFileSEConstAtt::New(fid,this,i,ms->getSupMeshWithName(_sup_mesh_name));
  _var_att.resize(nVarAttr);
  for(int i=0;i<nVarAttr;i++)
    _var_att[i]=MEDFileSEVarAtt::New(fid,this,i);
}

std::vector<const BigMemoryObject *> MEDFileStructureElement::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MCAuto<MEDFileSEConstAtt> >::const_iterator it=_cst_att.begin();it!=_cst_att.end();it++)
    ret.push_back(*it);
  for(std::vector< MCAuto<MEDFileSEVarAtt> >::const_iterator it=_var_att.begin();it!=_var_att.end();it++)
    ret.push_back(*it);
  return ret;
}

std::size_t MEDFileStructureElement::getHeapMemorySizeWithoutChildren() const
{
  return _name.capacity()+_cst_att.capacity()*sizeof(MCAuto<MEDFileSEConstAtt>)+_var_att.capacity()*sizeof(MCAuto<MEDFileSEVarAtt>);
}

void MEDFileStructureElement::writeLL(med_idt fid) const
{
}

std::string MEDFileStructureElement::getName() const
{
  return _name;
}

MCAuto<DataArray> MEDFileStructureElement::BuildFrom(med_attribute_type mat)
{
  MCAuto<DataArray> ret;
  switch(mat)
    {
    case MED_ATT_INT:
      {
        ret=DataArrayInt::New();
        break;
      }
    case MED_ATT_FLOAT64:
      {
        ret=DataArrayDouble::New();
        break;
      }
    case MED_ATT_NAME:
      {
        ret=DataArrayAsciiChar::New();
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileStructureElement::BuildFrom : not recognized type ! Only INT and FLOAT64 !");
    }
  return ret;
}

int MEDFileStructureElement::EffectiveNbCompo(med_attribute_type mat, int nbCompo)
{
  switch(mat)
    {
    case MED_ATT_INT:
    case MED_ATT_FLOAT64:
      return nbCompo;
    case MED_ATT_NAME:
      return nbCompo*MED_NAME_SIZE;
    default:
      throw INTERP_KERNEL::Exception("MEDFileStructureElement::BuildFrom : not recognized type ! Only INT and FLOAT64 !");
    }
}

int MEDFileStructureElement::getDynGT() const
{
  return _id_type;
}

TypeOfField MEDFileStructureElement::getEntity() const
{
  return _tof;
}

std::string MEDFileStructureElement::getMeshName() const
{
  return _sup_mesh_name;
}

std::vector<std::string> MEDFileStructureElement::getVarAtts() const
{
  std::vector<std::string> ret;
  for(std::vector< MCAuto<MEDFileSEVarAtt> >::const_iterator it=_var_att.begin();it!=_var_att.end();it++)
    if((*it).isNotNull())
      ret.push_back((*it)->getName());
  return ret;
}

const MEDFileSEVarAtt *MEDFileStructureElement::getVarAtt(const std::string& varName) const
{
  for(std::vector< MCAuto<MEDFileSEVarAtt> >::const_iterator it=_var_att.begin();it!=_var_att.end();it++)
    if((*it).isNotNull())
      if((*it)->getName()==varName)
        return *it;
  std::ostringstream oss; oss << "MEDFileStructureElement::getVarAtt : no var att with name \"" << varName << "\" !";
  throw INTERP_KERNEL::Exception(oss.str());
}

////////////////////

MEDFileStructureElements *MEDFileStructureElements::New(const std::string& fileName, const MEDFileMeshSupports *ms)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return New(fid,ms);
}

MEDFileStructureElements *MEDFileStructureElements::New(med_idt fid, const MEDFileMeshSupports *ms)
{
  return new MEDFileStructureElements(fid,ms);
}

MEDFileStructureElements *MEDFileStructureElements::New()
{
  return new MEDFileStructureElements;
}

std::vector<const BigMemoryObject *> MEDFileStructureElements::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MCAuto<MEDFileStructureElement> >::const_iterator it=_elems.begin();it!=_elems.end();it++)
    ret.push_back(*it);
  return ret;
}

std::size_t MEDFileStructureElements::getHeapMemorySizeWithoutChildren() const
{
  return _elems.capacity()*sizeof(MEDFileStructureElement);
}

void MEDFileStructureElements::writeLL(med_idt fid) const
{
}

MEDFileStructureElements::MEDFileStructureElements(med_idt fid, const MEDFileMeshSupports *ms)
{
  int nbSE(MEDnStructElement(fid));
  _elems.resize(nbSE);
  for(int i=0;i<nbSE;i++)
    _elems[i]=MEDFileStructureElement::New(fid,i,ms);
  _sup.takeRef(ms);
}

MEDFileStructureElements::MEDFileStructureElements()
{
}

MEDFileStructureElements::~MEDFileStructureElements()
{
}

int MEDFileStructureElements::getNumberOf() const
{
  return _elems.size();
}

std::vector<int> MEDFileStructureElements::getDynGTAvail() const
{
  std::vector<int> ret;
  for(std::vector< MCAuto<MEDFileStructureElement> >::const_iterator it=_elems.begin();it!=_elems.end();it++)
    {
      const MEDFileStructureElement *elt(*it);
      if(elt)
        ret.push_back(elt->getDynGT());
    }
  return ret;
}

const MEDFileStructureElement *MEDFileStructureElements::getWithGT(int idGT) const
{
  for(std::vector< MCAuto<MEDFileStructureElement> >::const_iterator it=_elems.begin();it!=_elems.end();it++)
    if((*it).isNotNull())
      {
        if((*it)->getDynGT()==idGT)
          return *it;
      }
  std::ostringstream oss; oss << "MEDFileStructureElements::getWithGT : no such geo type " << idGT << " !";
  throw INTERP_KERNEL::Exception(oss.str());
}

int MEDFileStructureElements::getNumberOfNodesPerSE(const std::string& seName) const
{
  if(seName=="MED_PARTICLE")
    return 1;
  const MEDFileStructureElement *se(getSEWithName(seName));
  std::string meshName(se->getMeshName());
  return _sup->getNumberOfNodesInConnOf(se->getEntity(),meshName);
}

const MEDFileStructureElement *MEDFileStructureElements::getSEWithName(const std::string& seName) const
{
  for(std::vector< MCAuto<MEDFileStructureElement> >::const_iterator it=_elems.begin();it!=_elems.end();it++)
    {
      if((*it).isNotNull())
        if((*it)->getName()==seName)
          return *it;
    }
  std::ostringstream oss; oss << "MEDFileStructureElements::getSEWithName : no such structure element with name " << seName << " !";
  throw INTERP_KERNEL::Exception(oss.str());
}

std::vector<std::string> MEDFileStructureElements::getVarAttsOf(const std::string& seName) const
{
  const MEDFileStructureElement *se(getSEWithName(seName));
  return se->getVarAtts();
}

const MEDFileSEVarAtt *MEDFileStructureElements::getVarAttOf(const std::string &seName, const std::string& varName) const
{
  const MEDFileStructureElement *se(getSEWithName(seName));
  return se->getVarAtt(varName);
}

const MEDFileUMesh *MEDFileStructureElements::getSupMeshWithName(const std::string& name) const
{
  return _sup->getSupMeshWithName(name);
}
