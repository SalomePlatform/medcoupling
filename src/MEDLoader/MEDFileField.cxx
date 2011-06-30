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

#include "MEDFileField.hxx"
#include "MEDFileMesh.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileUtilities.hxx"

#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldDiscretization.hxx"

#include "InterpKernelAutoPtr.hxx"
#include "CellModel.hxx"

#include <algorithm>

extern med_geometry_type typmai[MED_N_CELL_FIXED_GEO];
extern INTERP_KERNEL::NormalizedCellType typmai2[MED_N_CELL_FIXED_GEO];
extern med_geometry_type typmainoeud[1];
extern med_geometry_type typmai3[32];

using namespace ParaMEDMEM;

MEDFileFieldLoc *MEDFileFieldLoc::New(med_idt fid, const char *locName)
{
  return new MEDFileFieldLoc(fid,locName);
}

MEDFileFieldLoc *MEDFileFieldLoc::New(const char *locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w)
{
  return new MEDFileFieldLoc(locName,geoType,refCoo,gsCoo,w);
}

MEDFileFieldLoc::MEDFileFieldLoc(med_idt fid, const char *locName):_name(locName)
{
  med_geometry_type geotype;
  med_geometry_type sectiongeotype;
  int nsectionmeshcell;
  INTERP_KERNEL::AutoPtr<char> geointerpname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> sectionmeshname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  MEDlocalizationInfoByName(fid,locName,&geotype,&_dim,&_nb_gauss_pt,geointerpname,sectionmeshname,&nsectionmeshcell,&sectiongeotype);
  _geo_type=(INTERP_KERNEL::NormalizedCellType)(std::distance(typmai3,std::find(typmai3,typmai3+32,geotype)));
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_geo_type);
  _nb_node_per_cell=cm.getNumberOfNodes();
  _ref_coo.resize(_dim*_nb_node_per_cell);
  _gs_coo.resize(_dim*_nb_gauss_pt);
  _w.resize(_nb_gauss_pt);
  MEDlocalizationRd(fid,locName,MED_FULL_INTERLACE,&_ref_coo[0],&_gs_coo[0],&_w[0]);
}

MEDFileFieldLoc::MEDFileFieldLoc(const char *locName, INTERP_KERNEL::NormalizedCellType geoType,
                                 const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w):_name(locName),_geo_type(geoType),_ref_coo(refCoo),_gs_coo(gsCoo),
                                                                                                                                    _w(w)
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_geo_type);
  _dim=cm.getDimension();
  _nb_node_per_cell=cm.getNumberOfNodes();
  _nb_gauss_pt=_w.size();
}

bool MEDFileFieldLoc::isEqual(const MEDFileFieldLoc& other, double eps) const
{
  if(_name!=other._name)
    return false;
  if(_dim!=other._dim)
    return false;
  if(_nb_gauss_pt!=other._nb_gauss_pt)
    return false;
  if(_nb_node_per_cell!=other._nb_node_per_cell)
    return false;
  if(_geo_type!=other._geo_type)
    return false;
  if(MEDCouplingGaussLocalization::AreAlmostEqual(_ref_coo,other._ref_coo,eps))
    return false;
  if(MEDCouplingGaussLocalization::AreAlmostEqual(_gs_coo,other._gs_coo,eps))
    return false;
  if(MEDCouplingGaussLocalization::AreAlmostEqual(_w,other._w,eps))
    return false;
  
  return true;
}

void MEDFileFieldLoc::writeLL(med_idt fid) const
{
  MEDlocalizationWr(fid,_name.c_str(),typmai3[(int)_geo_type],_dim,&_ref_coo[0],MED_FULL_INTERLACE,_nb_gauss_pt,&_gs_coo[0],&_w[0],MED_NO_INTERPOLATION,MED_NO_MESH_SUPPORT);
}

void MEDFileFieldPerMeshPerTypePerDisc::assignFieldNoProfile(int offset, int nbOfCells, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobs& glob) throw(INTERP_KERNEL::Exception)
{
  _type=field->getTypeOfField();
  const DataArrayDouble *da=field->getArray();
  switch(_type)
    {
    case ON_CELLS:
      {
        _arr=da->selectByTupleId2(offset,offset+nbOfCells,1);
        _nval=nbOfCells;
        break;
      }
    case ON_GAUSS_NE:
      {
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=field->getDiscretization()->getOffsetArr(field->getMesh());
        const int *arrPtr=arr->getConstPointer();
        _nval=nbOfCells;
        _arr=da->selectByTupleId2(arrPtr[offset],arrPtr[offset+nbOfCells],1);
        break;
      }
    case ON_GAUSS_PT:
      {
        const MEDCouplingFieldDiscretization *disc=field->getDiscretization();
        const MEDCouplingGaussLocalization& gsLoc=field->getGaussLocalization(_loc_id);
        const MEDCouplingFieldDiscretizationGauss *disc2=dynamic_cast<const MEDCouplingFieldDiscretizationGauss *>(disc);
        if(!disc2)
          throw INTERP_KERNEL::Exception("assignFieldNoProfile : invalid call to this method ! Internal Error !");
        const DataArrayInt *dai=disc2->getArrayOfDiscIds();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> dai2=disc2->getOffsetArr(field->getMesh());
        const int *dai2Ptr=dai2->getConstPointer();
        int nbi=gsLoc.getWeights().size();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da2=dai->selectByTupleId2(offset,offset+nbOfCells,1);
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da3=da2->getIdsEqual(_loc_id);
        const int *da3Ptr=da3->getConstPointer();
        if(da3->getNumberOfTuples()!=nbOfCells)
          {//profile : for gauss even in NoProfile !!!
            std::ostringstream oss; oss << "Pfl_" << getName() << "_" << INTERP_KERNEL::CellModel::GetCellModel(getGeoType()).getRepr() << "_" << _loc_id;
            _profile=oss.str();
            da3->setName(_profile.c_str());
            glob.appendProfile(da3);
          }
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da4=DataArrayInt::New();
        _nval=da3->getNbOfElems();
        da4->alloc(_nval*nbi,1);
        int *da4Ptr=da4->getPointer();
        for(int i=0;i<_nval;i++)
          {
            int ref=dai2Ptr[offset+da3Ptr[i]];
            for(int j=0;j<nbi;j++)
              *da4Ptr++=ref+j;
          }
        std::ostringstream oss2; oss2 << "Loc_" << getName() << "_" << INTERP_KERNEL::CellModel::GetCellModel(getGeoType()).getRepr() << "_" << _loc_id;
        _localization=oss2.str();
        _arr=da->selectByTupleId(da4->getConstPointer(),da4->getConstPointer()+_nval*nbi);
        glob.appendLoc(_localization.c_str(),getGeoType(),gsLoc.getRefCoords(),gsLoc.getGaussCoords(),gsLoc.getWeights());
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::assignFieldNoProfile : not implemented yet for such discretization type of field !");
    }
}

void MEDFileFieldPerMeshPerTypePerDisc::assignFieldProfile(const char *pflName, const DataArrayInt *globIds, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobs& glob) throw(INTERP_KERNEL::Exception)
{
  if(pflName)
    _profile=pflName;
  else
    _profile.clear();
  _type=field->getTypeOfField();
  const DataArrayDouble *da=field->getArray();
  switch(_type)
    {
    case ON_CELLS:
      {
        _nval=globIds->getNumberOfTuples();
        _arr=da->selectByTupleId(globIds->getConstPointer(),globIds->getConstPointer()+_nval);
        break;
      }
    case ON_GAUSS_NE:
      {
        throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::assignFieldProfile : not implemented yet for profiles on gauss NE points !");
        /*MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=field->getDiscretization()->getOffsetArr(field->getMesh());
        const int *arrPtr=arr->getConstPointer();
        _nval=nbOfCells;
        _arr=da->selectByTupleId(arrPtr[offset],arrPtr[offset+nbOfCells],1);
        break;*/
      }
    case ON_GAUSS_PT:
      {
        throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::assignFieldProfile : not implemented yet for profiles on gauss points !");
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::assignFieldProfile : not implemented yet for such discretization type of field !");
    }
}

void MEDFileFieldPerMeshPerTypePerDisc::assignNodeFieldNoProfile(const MEDCouplingFieldDouble *field, MEDFieldFieldGlobs& glob) throw(INTERP_KERNEL::Exception)
{
  _arr=field->getArray()->deepCpy();
  _nval=field->getArray()->getNumberOfTuples();
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypePerDisc::New(MEDFileFieldPerMeshPerType *fath, med_idt fid, TypeOfField type, int profileIt) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMeshPerTypePerDisc(fath,fid,type,profileIt);
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypePerDisc::New(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int locId)
{
  return new MEDFileFieldPerMeshPerTypePerDisc(fath,type,locId);
}

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, med_idt fid, TypeOfField atype, int profileIt) throw(INTERP_KERNEL::Exception)
try:_type(atype),_father(fath),_profile_it(profileIt)
{
  INTERP_KERNEL::AutoPtr<char> locname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> pflname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  std::string fieldName=getName();
  std::string meshName=getMeshName();
  int iteration=getIteration();
  int order=getOrder();
  const std::vector<std::string>& infos=getInfo();
  TypeOfField type=getType();
  INTERP_KERNEL::NormalizedCellType geoType=getGeoType();
  int profilesize,nbi;
  med_geometry_type mgeoti;
  med_entity_type menti=MEDFileFieldPerMeshPerType::ConvertIntoMEDFileType(type,geoType,mgeoti);
  _nval=MEDfieldnValueWithProfile(fid,fieldName.c_str(),iteration,order,menti,mgeoti,profileIt,MED_COMPACT_PFLMODE,
                                  pflname,&profilesize,locname,&nbi);
  _arr=DataArrayDouble::New();
  _arr->alloc(_nval*nbi,infos.size());
  MEDfieldValueWithProfileRd(fid,fieldName.c_str(),iteration,order,menti,mgeoti,MED_COMPACT_PFLMODE,
                             pflname,MED_FULL_INTERLACE,MED_ALL_CONSTITUENT,reinterpret_cast<unsigned char*>(_arr->getPointer()));
  _profile=MEDLoaderBase::buildStringFromFortran(pflname,MED_NAME_SIZE);
  _localization=MEDLoaderBase::buildStringFromFortran(locname,MED_NAME_SIZE);
  if(type==ON_CELLS && !_localization.empty())
    setType(ON_GAUSS_PT);
}
catch(INTERP_KERNEL::Exception& e)
{
  throw e;
}

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int locId):_type(type),_father(fath),_loc_id(locId)
{
}

const MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerTypePerDisc::getFather() const
{
  return _father;
}

int MEDFileFieldPerMeshPerTypePerDisc::getIteration() const
{
  return _father->getIteration();
}

int MEDFileFieldPerMeshPerTypePerDisc::getOrder() const
{
  return _father->getOrder();
}

double MEDFileFieldPerMeshPerTypePerDisc::getTime() const
{
  return _father->getTime();
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getName() const
{
  return _father->getName();
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getMeshName() const
{
  return _father->getMeshName();
}

TypeOfField MEDFileFieldPerMeshPerTypePerDisc::getType() const
{
  return _type;
}

void MEDFileFieldPerMeshPerTypePerDisc::setType(TypeOfField newType)
{
  _type=newType;
}

INTERP_KERNEL::NormalizedCellType MEDFileFieldPerMeshPerTypePerDisc::getGeoType() const
{
  return _father->getGeoType();
}

int MEDFileFieldPerMeshPerTypePerDisc::getNumberOfComponents() const
{
  return _father->getNumberOfComponents();
}

int MEDFileFieldPerMeshPerTypePerDisc::getNumberOfTuples() const
{
  return _arr->getNumberOfTuples();
}

const std::vector<std::string>& MEDFileFieldPerMeshPerTypePerDisc::getInfo() const
{
  return _father->getInfo();
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getProfile() const
{
  return _profile;
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getLocalization() const
{
  return _localization;
}

void MEDFileFieldPerMeshPerTypePerDisc::getFieldAtLevel(TypeOfField type, const MEDFieldFieldGlobs *glob, std::vector<const DataArrayDouble *>& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs, std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const
{
  if(type!=_type)
    return ;
  dads.push_back(_arr);
  geoTypes.push_back(getGeoType());
  if(_profile.empty())
    pfls.push_back(0);
  else
    {
      pfls.push_back(glob->getProfile(_profile.c_str()));
    }
  if(_localization.empty())
    locs.push_back(-1);
  else
    {
      locs.push_back(glob->getLocalizationId(_localization.c_str()));
    }
}

void MEDFileFieldPerMeshPerTypePerDisc::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  TypeOfField type=getType();
  INTERP_KERNEL::NormalizedCellType geoType=getGeoType();
  med_geometry_type mgeoti;
  med_entity_type menti=MEDFileFieldPerMeshPerType::ConvertIntoMEDFileType(type,geoType,mgeoti);
  MEDfieldValueWithProfileWr(fid,getName().c_str(),getIteration(),getOrder(),getTime(),menti,mgeoti,
                             MED_COMPACT_PFLMODE,_profile.c_str(),_localization.c_str(),MED_FULL_INTERLACE,MED_ALL_CONSTITUENT,_nval,
                             reinterpret_cast<const unsigned char*>(_arr->getConstPointer()));
}

int MEDFileFieldPerMeshPerTypePerDisc::ConvertType(TypeOfField type, int locId) throw(INTERP_KERNEL::Exception)
{
  switch(type)
    {
    case ON_CELLS:
      return -2;
    case ON_GAUSS_NE:
      return -1;
    case ON_GAUSS_PT:
      return locId;
    default:
      throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::ConvertType : not managed type of field !");
    }
}

MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerType::New(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMeshPerType(fath,geoType);
}

void MEDFileFieldPerMeshPerType::assignFieldNoProfile(int offset, int nbOfCells, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobs& glob) throw(INTERP_KERNEL::Exception)
{
  std::vector<int> pos=addNewEntryIfNecessary(field,offset,nbOfCells);
  for(std::vector<int>::const_iterator it=pos.begin();it!=pos.end();it++)
    _field_pm_pt_pd[*it]->assignFieldNoProfile(offset,nbOfCells,field,glob);
}

void MEDFileFieldPerMeshPerType::assignFieldProfile(const DataArrayInt *globIds, DataArrayInt *locIds, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobs& glob) throw(INTERP_KERNEL::Exception)
{
  std::vector<int> pos=addNewEntryIfNecessary(field,globIds);
  if(locIds)
    {
      //
      std::string pflName(locIds->getName());
      if(pflName.empty())
        throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerType::assignFieldProfile : existing profile with empty name !");
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_geo_type);
      std::ostringstream oss; oss << locIds->getName() << "_" <<  cm.getRepr();
      locIds->setName(oss.str().c_str());
      glob.appendProfile(locIds);
      //
      for(std::vector<int>::const_iterator it=pos.begin();it!=pos.end();it++)
        _field_pm_pt_pd[*it]->assignFieldProfile(oss.str().c_str(),globIds,field,glob);
    }
  else
    {
      for(std::vector<int>::const_iterator it=pos.begin();it!=pos.end();it++)
        _field_pm_pt_pd[*it]->assignFieldProfile(0,globIds,field,glob);
    }
}

void MEDFileFieldPerMeshPerType::assignNodeFieldNoProfile(const MEDCouplingFieldDouble *field, MEDFieldFieldGlobs& glob) throw(INTERP_KERNEL::Exception)
{
  _field_pm_pt_pd.resize(1);
  _field_pm_pt_pd[0]=MEDFileFieldPerMeshPerTypePerDisc::New(this,ON_NODES,-3);
  _field_pm_pt_pd[0]->assignNodeFieldNoProfile(field,glob);
}

std::vector<int> MEDFileFieldPerMeshPerType::addNewEntryIfNecessary(const MEDCouplingFieldDouble *field, int offset, int nbOfCells) throw(INTERP_KERNEL::Exception)
{
  TypeOfField type=field->getTypeOfField();
  if(type!=ON_GAUSS_PT)
    {
      int locIdToFind=MEDFileFieldPerMeshPerTypePerDisc::ConvertType(type,0);
      int sz=_field_pm_pt_pd.size();
      bool found=false;
      for(int j=0;j<sz && !found;j++)
        {
          if(_field_pm_pt_pd[j]->getLocId()==locIdToFind)
            {
              _field_pm_pt_pd[j]=MEDFileFieldPerMeshPerTypePerDisc::New(this,type,locIdToFind);
              found=true;
            }
        }
      if(!found)
        {
          _field_pm_pt_pd.resize(sz+1);
          _field_pm_pt_pd[sz]=MEDFileFieldPerMeshPerTypePerDisc::New(this,type,locIdToFind);
        }
      std::vector<int> ret(1,0);
      return ret;
    }
  else
    {
      std::vector<int> ret2=addNewEntryIfNecessaryGauss(field,offset,nbOfCells);
      int sz2=ret2.size();
      std::vector<int> ret3(sz2);
      int k=0;
      for(int i=0;i<sz2;i++)
        {
          int sz=_field_pm_pt_pd.size();
          int locIdToFind=ret2[i];
          bool found=false;
          for(int j=0;j<sz && !found;j++)
            {
              if(_field_pm_pt_pd[j]->getLocId()==locIdToFind)
                {
                  _field_pm_pt_pd[j]=MEDFileFieldPerMeshPerTypePerDisc::New(this,type,locIdToFind);
                  ret3[k++]=j;
                  found=true;
                }
            }
          if(!found)
            {
              _field_pm_pt_pd.resize(sz+1);
              _field_pm_pt_pd[sz]=MEDFileFieldPerMeshPerTypePerDisc::New(this,type,locIdToFind);
              ret3[k++]=sz;
            }
        }
      return ret3;
    }
}

std::vector<int> MEDFileFieldPerMeshPerType::addNewEntryIfNecessaryGauss(const MEDCouplingFieldDouble *field, int offset, int nbOfCells) throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingFieldDiscretization *disc=field->getDiscretization();
  const MEDCouplingFieldDiscretizationGauss *disc2=dynamic_cast<const MEDCouplingFieldDiscretizationGauss *>(disc);
  if(!disc2)
    throw INTERP_KERNEL::Exception("addNewEntryIfNecessaryGauss : invalid call to this method ! Internal Error !");
  const DataArrayInt *da=disc2->getArrayOfDiscIds();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da2=da->selectByTupleId2(offset,offset+nbOfCells,1);
  std::set<int> retTmp=da2->getDifferentValues();
  if(retTmp.find(-1)!=retTmp.end())
    throw INTERP_KERNEL::Exception("addNewEntryIfNecessaryGauss : some cells have no dicretization description !");
  std::vector<int> ret(retTmp.begin(),retTmp.end());
  return ret;
}

std::vector<int> MEDFileFieldPerMeshPerType::addNewEntryIfNecessary(const MEDCouplingFieldDouble *field, const DataArrayInt *subCells) throw(INTERP_KERNEL::Exception)
{
  TypeOfField type=field->getTypeOfField();
  if(type!=ON_GAUSS_PT)
    {
      int locIdToFind=MEDFileFieldPerMeshPerTypePerDisc::ConvertType(type,0);
      int sz=_field_pm_pt_pd.size();
      bool found=false;
      for(int j=0;j<sz && !found;j++)
        {
          if(_field_pm_pt_pd[j]->getLocId()==locIdToFind)
            {
              _field_pm_pt_pd[j]=MEDFileFieldPerMeshPerTypePerDisc::New(this,type,locIdToFind);
              found=true;
            }
        }
      if(!found)
        {
          _field_pm_pt_pd.resize(sz+1);
          _field_pm_pt_pd[sz]=MEDFileFieldPerMeshPerTypePerDisc::New(this,type,locIdToFind);
        }
      std::vector<int> ret(1,0);
      return ret;
    }
  else
    {
      std::vector<int> ret2=addNewEntryIfNecessaryGauss(field,subCells);
      int sz2=ret2.size();
      std::vector<int> ret3(sz2);
      int k=0;
      for(int i=0;i<sz2;i++)
        {
          int sz=_field_pm_pt_pd.size();
          int locIdToFind=ret2[i];
          bool found=false;
          for(int j=0;j<sz && !found;j++)
            {
              if(_field_pm_pt_pd[j]->getLocId()==locIdToFind)
                {
                  _field_pm_pt_pd[j]=MEDFileFieldPerMeshPerTypePerDisc::New(this,type,locIdToFind);
                  ret3[k++]=j;
                  found=true;
                }
            }
          if(!found)
            {
              _field_pm_pt_pd.resize(sz+1);
              _field_pm_pt_pd[sz]=MEDFileFieldPerMeshPerTypePerDisc::New(this,type,locIdToFind);
              ret3[k++]=sz;
            }
        }
      return ret3;
    }
}

std::vector<int> MEDFileFieldPerMeshPerType::addNewEntryIfNecessaryGauss(const MEDCouplingFieldDouble *field, const DataArrayInt *subCells) throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingFieldDiscretization *disc=field->getDiscretization();
  const MEDCouplingFieldDiscretizationGauss *disc2=dynamic_cast<const MEDCouplingFieldDiscretizationGauss *>(disc);
  if(!disc2)
    throw INTERP_KERNEL::Exception("addNewEntryIfNecessaryGauss : invalid call to this method ! Internal Error !");
  const DataArrayInt *da=disc2->getArrayOfDiscIds();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da2=da->selectByTupleId(subCells->getConstPointer(),subCells->getConstPointer()+subCells->getNumberOfTuples());
  std::set<int> retTmp=da2->getDifferentValues();
  if(retTmp.find(-1)!=retTmp.end())
    throw INTERP_KERNEL::Exception("addNewEntryIfNecessaryGauss : some cells have no dicretization description !");
  std::vector<int> ret(retTmp.begin(),retTmp.end());
  return ret;
}

const MEDFileFieldPerMesh *MEDFileFieldPerMeshPerType::getFather() const
{
  return _father;
}

void MEDFileFieldPerMeshPerType::getDimension(int& dim) const
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_geo_type);
  int curDim=(int)cm.getDimension();
  dim=std::max(dim,curDim);
}

int MEDFileFieldPerMeshPerType::getIteration() const
{
  return _father->getIteration();
}

int MEDFileFieldPerMeshPerType::getOrder() const
{
  return _father->getOrder();
}

double MEDFileFieldPerMeshPerType::getTime() const
{
  return _father->getTime();
}

std::string MEDFileFieldPerMeshPerType::getName() const
{
  return _father->getName();
}

std::string MEDFileFieldPerMeshPerType::getMeshName() const
{
  return _father->getMeshName();
}

INTERP_KERNEL::NormalizedCellType MEDFileFieldPerMeshPerType::getGeoType() const
{
  return _geo_type;
}


int MEDFileFieldPerMeshPerType::getNumberOfComponents() const
{
  return _father->getNumberOfComponents();
}

const std::vector<std::string>& MEDFileFieldPerMeshPerType::getInfo() const
{
  return _father->getInfo();
}

std::vector<std::string> MEDFileFieldPerMeshPerType::getPflsReallyUsed() const
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it1=_field_pm_pt_pd.begin();it1!=_field_pm_pt_pd.end();it1++)
    {
      std::string tmp=(*it1)->getProfile();
      if(!tmp.empty())
        ret.push_back(tmp);
    }
  return ret;
}

std::vector<std::string> MEDFileFieldPerMeshPerType::getLocsReallyUsed() const
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it1=_field_pm_pt_pd.begin();it1!=_field_pm_pt_pd.end();it1++)
    {
      std::string tmp=(*it1)->getLocalization();
      if(!tmp.empty() && tmp!=MED_GAUSS_ELNO)
        ret.push_back(tmp);
    }
  return ret;
}

void MEDFileFieldPerMeshPerType::getFieldAtLevel(int meshDim, TypeOfField type, const MEDFieldFieldGlobs *glob, std::vector<const DataArrayDouble *>& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs, std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const
{
  if(_geo_type!=INTERP_KERNEL::NORM_ERROR)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_geo_type);
      if(meshDim!=(int)cm.getDimension())
        return ;
    }
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    (*it)->getFieldAtLevel(type,glob,dads,pfls,locs,geoTypes);
}

MEDFileFieldPerMeshPerType::MEDFileFieldPerMeshPerType(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception):_father(fath),_geo_type(geoType)
{
}

void MEDFileFieldPerMeshPerType::finishLoading(med_idt fid, TypeOfField type) throw(INTERP_KERNEL::Exception)
{
  INTERP_KERNEL::NormalizedCellType geoType=getGeoType();
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> locName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  med_geometry_type mgeoti;
  med_entity_type menti=ConvertIntoMEDFileType(type,geoType,mgeoti);
  int nbProfiles=MEDfieldnProfile(fid,getName().c_str(),getIteration(),getOrder(),menti,mgeoti,pflName,locName);
  _field_pm_pt_pd.resize(nbProfiles);
  for(int i=0;i<nbProfiles;i++)
    {
      _field_pm_pt_pd[i]=MEDFileFieldPerMeshPerTypePerDisc::New(this,fid,type,i+1);
    }
}

void MEDFileFieldPerMeshPerType::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    {
      (*it)->copyOptionsFrom(*this);
      (*it)->writeLL(fid);
    }
}

med_entity_type MEDFileFieldPerMeshPerType::ConvertIntoMEDFileType(TypeOfField ikType, INTERP_KERNEL::NormalizedCellType ikGeoType, med_geometry_type& medfGeoType)
{
  switch(ikType)
    {
    case ON_CELLS:
      medfGeoType=typmai3[(int)ikGeoType];
      return MED_CELL;
    case ON_NODES:
      medfGeoType=MED_NONE;
      return MED_NODE;
    case ON_GAUSS_NE:
      medfGeoType=typmai3[(int)ikGeoType];
      return MED_NODE_ELEMENT;
    case ON_GAUSS_PT:
      medfGeoType=typmai3[(int)ikGeoType];
      return MED_CELL;
    default:
      throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerType::ConvertIntoMEDFileType : unexpected entity type ! internal error");
    }
  return MED_UNDEF_ENTITY_TYPE;
}

MEDFileFieldPerMesh *MEDFileFieldPerMesh::New(MEDFileField1TSWithoutDAS *fath, int meshCsit, int meshIteration, int meshOrder)
{
  return new MEDFileFieldPerMesh(fath,meshCsit,meshIteration,meshOrder);
}

MEDFileFieldPerMesh *MEDFileFieldPerMesh::New(MEDFileField1TSWithoutDAS *fath, const MEDCouplingMesh *mesh)
{
  return new MEDFileFieldPerMesh(fath,mesh);
}

void MEDFileFieldPerMesh::copyTinyInfoFrom(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception)
{
  _mesh_name=mesh->getName();
  mesh->getTime(_mesh_iteration,_mesh_order);
}

void MEDFileFieldPerMesh::assignFieldProfile(const std::vector<int>& code, const std::vector<DataArrayInt *>& globIdsPerType, const std::vector<DataArrayInt *>& idsPerType, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobs& glob) throw(INTERP_KERNEL::Exception)
{
  int nbOfTypes=code.size()/3;
  bool isProfile=false;
  for(int i=0;i<nbOfTypes;i++)
    if(code[3*i+2]!=-1)
      isProfile=true;
  if(!isProfile)
    {
      if(globIdsPerType.empty())
        assignFieldNoProfileNoRenum(code,field,glob);
      else
        assignFieldProfileGeneral(code,globIdsPerType,idsPerType,field,glob);
    }
  else
    assignFieldProfileGeneral(code,globIdsPerType,idsPerType,field,glob);
}

void MEDFileFieldPerMesh::assignFieldNoProfileNoRenum(const std::vector<int>& code, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobs& glob) throw(INTERP_KERNEL::Exception)
{
  int nbOfTypes=code.size()/3;
  int offset=0;
  for(int i=0;i<nbOfTypes;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)code[3*i];
      int nbOfCells=code[3*i+1];
      int pos=addNewEntryIfNecessary(type);
      _field_pm_pt[pos]->assignFieldNoProfile(offset,nbOfCells,field,glob);
      offset+=nbOfCells;
    }
}

/*!
 * This method is the most general one. No optimization is done here.
 */
void MEDFileFieldPerMesh::assignFieldProfileGeneral(const std::vector<int>& code, const std::vector<DataArrayInt *>& globIdsPerType, const std::vector<DataArrayInt *>& idsPerType, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobs& glob) throw(INTERP_KERNEL::Exception)
{
  int nbOfTypes=code.size()/3;
  for(int i=0;i<nbOfTypes;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)code[3*i];
      int pos=addNewEntryIfNecessary(type);
      DataArrayInt *pfl=0;
      if(code[3*i+2]!=-1)
        pfl=idsPerType[code[3*i+2]];
      _field_pm_pt[pos]->assignFieldProfile(globIdsPerType[i],pfl,field,glob);
    }
}

void MEDFileFieldPerMesh::assignNodeFieldNoProfile(const MEDCouplingFieldDouble *field, MEDFieldFieldGlobs& glob) throw(INTERP_KERNEL::Exception)
{
  int pos=addNewEntryIfNecessary(INTERP_KERNEL::NORM_ERROR);
  _field_pm_pt[pos]->assignNodeFieldNoProfile(field,glob);
}

void MEDFileFieldPerMesh::finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  INTERP_KERNEL::AutoPtr<char> meshName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> locName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  _field_pm_pt.clear();
  for(int i=0;i<MED_N_CELL_FIXED_GEO;i++)
    {
      int nbProfile=MEDfield23nProfile(fid,getName().c_str(),getIteration(),getOrder(),MED_CELL,typmai[i],_mesh_csit,meshName,pflName,locName);
      if(nbProfile>0)
        {
          _field_pm_pt.resize(_field_pm_pt.size()+1);
          _field_pm_pt.back()=MEDFileFieldPerMeshPerType::New(this,typmai2[i]);
          _mesh_name=MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1);
          _field_pm_pt.back()->finishLoading(fid,ON_CELLS);
        }
      nbProfile=MEDfield23nProfile(fid,getName().c_str(),getIteration(),getOrder(),MED_NODE_ELEMENT,typmai[i],_mesh_csit,meshName,pflName,locName);
      if(nbProfile>0)
        {
          _field_pm_pt.resize(_field_pm_pt.size()+1);
          _field_pm_pt.back()=MEDFileFieldPerMeshPerType::New(this,typmai2[i]);
          _mesh_name=MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1);
          _field_pm_pt.back()->finishLoading(fid,ON_GAUSS_NE);
        }
    }
  int nbProfile=MEDfield23nProfile(fid,getName().c_str(),getIteration(),getOrder(),MED_NODE,MED_NONE,_mesh_csit,meshName,pflName,locName);
  if(nbProfile>0)
    {
      _field_pm_pt.resize(_field_pm_pt.size()+1);
      _field_pm_pt.back()=MEDFileFieldPerMeshPerType::New(this,INTERP_KERNEL::NORM_ERROR);
      _mesh_name=MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1);
      _field_pm_pt.back()->finishLoading(fid,ON_NODES);
    }
}

void MEDFileFieldPerMesh::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  int nbOfTypes=_field_pm_pt.size();
  for(int i=0;i<nbOfTypes;i++)
    {
      _field_pm_pt[i]->copyOptionsFrom(*this);
      _field_pm_pt[i]->writeLL(fid);
    }
}

void MEDFileFieldPerMesh::getDimension(int& dim) const
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->getDimension(dim);
}

double MEDFileFieldPerMesh::getTime() const
{
  return _father->getTime();
}

int MEDFileFieldPerMesh::getIteration() const
{
  return _father->getIteration();
}

const std::string& MEDFileFieldPerMesh::getDtUnit() const
{
  return _father->getDtUnit();
}

int MEDFileFieldPerMesh::getOrder() const
{
  return _father->getOrder();
}

std::string MEDFileFieldPerMesh::getName() const
{
  return _father->getName();
}

int MEDFileFieldPerMesh::getNumberOfComponents() const
{
  return _father->getNumberOfComponents();
}

const std::vector<std::string>& MEDFileFieldPerMesh::getInfo() const
{
  return _father->getInfo();
}

/*!
 * geoTypes,dads,pfls,locs are input parameters. They should have the same size.
 * Before the call of this method 'geoTypes','dads','pfls','locs' must be reorganized so that types in geoTypes are contiguous and ordered following typmai2 array.
 * It returns 2 output vectors :
 * - 'code' of size 3*sz where sz is the number of different values into 'geoTypes'
 * - 'notNullPfls' contains sz2 values that are extracted from 'pfls' in which null profiles have been removed.
 * 'code' and 'notNullPfls' are in MEDCouplingUMesh::checkTypeConsistencyAndContig format.
 */
void MEDFileFieldPerMesh::SortArraysPerType(const MEDFieldFieldGlobs *glob, const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector<const DataArrayDouble *>& dads, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& locs, std::vector<int>& code, std::vector<DataArrayInt *>& notNullPfls)
{
  int notNullPflsSz=0;
  int nbOfArrs=geoTypes.size();
  for(int i=0;i<nbOfArrs;i++)
    if(pfls[i])
      notNullPflsSz++;
  std::set<INTERP_KERNEL::NormalizedCellType> geoTypes3(geoTypes.begin(),geoTypes.end());
  int nbOfDiffGeoTypes=geoTypes3.size();
  code.resize(3*nbOfDiffGeoTypes);
  notNullPfls.resize(notNullPflsSz);
  notNullPflsSz=0;
  int j=0;
  for(int i=0;i<nbOfDiffGeoTypes;i++)
    {
      int startZone=j;
      INTERP_KERNEL::NormalizedCellType refType=geoTypes[j];
      std::vector<const DataArrayInt *> notNullTmp;
      if(pfls[j])
        notNullTmp.push_back(pfls[j]);
      j++;
      for(;j<nbOfArrs;j++)
        if(geoTypes[j]==refType)
          {
            if(pfls[j])
              notNullTmp.push_back(pfls[j]);
          }
        else
          break;
      std::vector<const DataArrayDouble *> tmpDads(dads.begin()+startZone,dads.begin()+j);
      std::vector<const DataArrayInt *> tmpPfls(pfls.begin()+startZone,pfls.begin()+j);
      std::vector<int> tmpLocs(locs.begin()+startZone,locs.begin()+j);
      code[3*i]=(int)refType;
      code[3*i+1]=ComputeNbOfElems(glob,tmpDads,tmpLocs);
      if(notNullTmp.empty())
        code[3*i+2]=-1;
      else
        {
          notNullPfls[notNullPflsSz]=DataArrayInt::Aggregate(notNullTmp);
          code[3*i+2]=notNullPflsSz++;
        }
    }
}

/*!
 * 'dads' and 'locs' are input parameters that should have same size sz. sz should be >=1.
 */
int MEDFileFieldPerMesh::ComputeNbOfElems(const MEDFieldFieldGlobs *glob, const std::vector<const DataArrayDouble *>& dads, const std::vector<int>& locs) throw(INTERP_KERNEL::Exception)
{
  int sz=dads.size();
  int ret=0;
  for(int i=0;i<sz;i++)
    {
      if(locs[i]==-1)
        ret+=dads[i]->getNumberOfTuples();
      else
        {
          int nbOfGaussPtPerCell=glob->getNbOfGaussPtPerCell(locs[i]);
          ret+=dads[i]->getNumberOfTuples()/nbOfGaussPtPerCell;
        }
    }
  return ret;
}

std::vector<std::string> MEDFileFieldPerMesh::getPflsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsed();
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

std::vector<std::string> MEDFileFieldPerMesh::getLocsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsed();
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

MEDCouplingFieldDouble *MEDFileFieldPerMesh::getFieldOnMeshAtLevel(TypeOfField type, int meshDim, const MEDFieldFieldGlobs *glob, const MEDCouplingMesh *mesh, bool& isPfl) const throw(INTERP_KERNEL::Exception)
{
  if(_field_pm_pt.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::getFieldOnMeshAtLevel : no types field set !");
  //
  std::vector<const DataArrayDouble *> dads;
  std::vector<const DataArrayInt *> pfls;
  std::vector<DataArrayInt *> notNullPflsPerGeoType;
  std::vector<int> locs,code;
  std::vector<INTERP_KERNEL::NormalizedCellType> geoTypes;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->getFieldAtLevel(meshDim,type,glob,dads,pfls,locs,geoTypes);
  // Sort by types
  SortArraysPerType(glob,geoTypes,dads,pfls,locs,code,notNullPflsPerGeoType);
  if(code.empty())
    {
      std::ostringstream oss; oss << "MEDFileFieldPerMesh::getFieldOnMeshAtLevel : " << "The field \"" << getName() << "\" exists but not with such spatial discretization or such dimension specified !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  //
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > notNullPflsPerGeoType2(notNullPflsPerGeoType.begin(),notNullPflsPerGeoType.end());
  std::vector< const DataArrayInt *> notNullPflsPerGeoType3(notNullPflsPerGeoType.begin(),notNullPflsPerGeoType.end());
  if(type!=ON_NODES)
    {
      DataArrayInt *arr=mesh->checkTypeConsistencyAndContig(code,notNullPflsPerGeoType3);
      if(!arr)
        return finishField(type,glob,dads,locs,mesh,isPfl);
      else
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr2(arr);
          return finishField2(type,glob,dads,locs,mesh,arr,isPfl);
        }
    }
  else
    {
      if(code.size()!=3)
        throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::getFieldOnMeshAtLevel : internal error #1 !");
      int nb=code[1];
      if(code[2]==-1)
        {
          if(nb!=mesh->getNumberOfNodes())
            {
              std::ostringstream oss; oss << "MEDFileFieldPerMesh::getFieldOnMeshAtLevel : There is a problem there is " << nb << " nodes in field whereas there is " << mesh->getNumberOfNodes();
              oss << " nodes in mesh !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          return finishField(type,glob,dads,locs,mesh,isPfl);
        }
      else
        return finishField3(glob,dads,locs,mesh,notNullPflsPerGeoType3[0],isPfl);
    }
}

int MEDFileFieldPerMesh::addNewEntryIfNecessary(INTERP_KERNEL::NormalizedCellType type)
{
  int i=0;
  int pos=std::distance(typmai2,std::find(typmai2,typmai2+MED_N_CELL_FIXED_GEO,type));
  std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::iterator it2=_field_pm_pt.begin();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++,i++)
    {
      INTERP_KERNEL::NormalizedCellType curType=(*it)->getGeoType();
      if(type==curType)
        return i;
      else
        {
          int pos2=std::distance(typmai2,std::find(typmai2,typmai2+MED_N_CELL_FIXED_GEO,curType));
          if(pos>pos2)
            it2=it+1;
        }
    }
  int ret=std::distance(_field_pm_pt.begin(),it2);
  _field_pm_pt.insert(it2,MEDFileFieldPerMeshPerType::New(this,type));
  return ret;
}

MEDCouplingFieldDouble *MEDFileFieldPerMesh::finishField(TypeOfField type, const MEDFieldFieldGlobs *glob,
                                                         const std::vector<const DataArrayDouble *>& dads, const std::vector<int>& locs,
                                                         const MEDCouplingMesh *mesh, bool& isPfl) const throw(INTERP_KERNEL::Exception)
{
  isPfl=false;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(type,ONE_TIME);
  ret->setMesh(mesh); ret->setName(getName().c_str()); ret->setTime(getTime(),getIteration(),getOrder()); ret->setTimeUnit(getDtUnit().c_str());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> da=DataArrayDouble::Aggregate(dads);
  const std::vector<std::string>& infos=getInfo();
  int nbOfComp=infos.size();
  for(int i=0;i<nbOfComp;i++)
    da->setInfoOnComponent(i,infos[i].c_str());
  ret->setArray(da);
  if(type==ON_GAUSS_PT)
    {
      int offset=0;
      int nbOfArrs=dads.size();
      for(int i=0;i<nbOfArrs;i++)
        {
          std::vector<const DataArrayDouble *> dads2(1,dads[i]); const std::vector<int> locs2(1,locs[i]);
          int nbOfElems=ComputeNbOfElems(glob,dads2,locs2);
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> di=DataArrayInt::New();
          di->alloc(nbOfElems,1);
          di->iota(offset);
          const MEDFileFieldLoc& fl=glob->getLocalizationFromId(locs[i]);
          ret->setGaussLocalizationOnCells(di->getConstPointer(),di->getConstPointer()+nbOfElems,fl.getRefCoords(),fl.getGaussCoords(),fl.getGaussWeights());
          offset+=nbOfElems;
        }
    }
  //
  ret->incrRef();
  return ret;
}

/*!
 * This method is an extension of MEDFileFieldPerMesh::finishField method. It deals with profiles. This method should be called when type is different from ON_NODES.
 * No check of this is performed. 'da' array contains an array in old2New style to be applyied to mesh to obtain the right support.
 * The order of cells in the returned field is those imposed by the profile.
 */
MEDCouplingFieldDouble *MEDFileFieldPerMesh::finishField2(TypeOfField type, const MEDFieldFieldGlobs *glob,
                                                          const std::vector<const DataArrayDouble *>& dads, const std::vector<int>& locs,
                                                          const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl) const throw(INTERP_KERNEL::Exception)
{
  if(da->isIdentity())
    {
      int nbOfTuples=da->getNumberOfTuples();
      if(nbOfTuples==ComputeNbOfElems(glob,dads,locs))
        return finishField(type,glob,dads,locs,mesh,isPfl);
    }
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=finishField(type,glob,dads,locs,mesh,isPfl);
  isPfl=true;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m2=mesh->buildPart(da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
  m2->setName(mesh->getName());
  ret->setMesh(m2);
  ret->incrRef();
  return ret;
}

/*!
 * This method is the complement of MEDFileFieldPerMesh::finishField2 method except that this method works for node profiles.
 */
MEDCouplingFieldDouble *MEDFileFieldPerMesh::finishField3(const MEDFieldFieldGlobs *glob,
                                                          const std::vector<const DataArrayDouble *>& dads, const std::vector<int>& locs,
                                                          const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl) const throw(INTERP_KERNEL::Exception)
{
  if(da->isIdentity())
    {
      int nbOfTuples=da->getNumberOfTuples();
      if(nbOfTuples==ComputeNbOfElems(glob,dads,locs))
        return finishField(ON_NODES,glob,dads,locs,mesh,isPfl);
    }
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=finishField(ON_NODES,glob,dads,locs,mesh,isPfl);
  isPfl=true;
  DataArrayInt *arr2=0;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cellIds=mesh->getCellIdsFullyIncludedInNodeIds(da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> mesh2=mesh->buildPartAndReduceNodes(cellIds->getConstPointer(),cellIds->getConstPointer()+cellIds->getNbOfElems(),arr2);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr3(arr2);
  int nnodes=mesh2->getNumberOfNodes();
  if(nnodes==da->getNbOfElems())
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da3=da->transformWithIndArrR(arr2->getConstPointer(),arr2->getConstPointer()+arr2->getNbOfElems());
      ret->getArray()->renumberInPlace(da3->getConstPointer());
      mesh2->setName(mesh->getName());
      ret->setMesh(mesh2);
      ret->incrRef();
      return ret;
    }
  else
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::finishField3 : not implemented yet !");
  return 0;
}

MEDFileFieldPerMesh::MEDFileFieldPerMesh(MEDFileField1TSWithoutDAS *fath, int meshCsit, int meshIteration, int meshOrder):_mesh_iteration(meshIteration),_mesh_order(meshOrder),
                                                                                                                          _mesh_csit(meshCsit),_father(fath)
{
}

MEDFileFieldPerMesh::MEDFileFieldPerMesh(MEDFileField1TSWithoutDAS *fath, const MEDCouplingMesh *mesh):_father(fath)
{
  copyTinyInfoFrom(mesh);
}

void MEDFieldFieldGlobs::loadProfileInFile(med_idt fid, int id, const char *pflName) throw(INTERP_KERNEL::Exception)
{
  _pfls[id]=DataArrayInt::New();
  int lgth=MEDprofileSizeByName(fid,pflName);
  _pfls[id]->setName(pflName);
  _pfls[id]->alloc(lgth,1);
  MEDprofileRd(fid,pflName,_pfls[id]->getPointer());
  _pfls[id]->applyLin(1,-1,0);//Converting into C format
}

void MEDFieldFieldGlobs::loadProfileInFile(med_idt fid, int i)
{
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  int sz;
  MEDprofileInfo(fid,i+1,pflName,&sz);
  std::string pflCpp=MEDLoaderBase::buildStringFromFortran(pflName,MED_NAME_SIZE);
  _pfls[i]=DataArrayInt::New();
  _pfls[i]->alloc(sz,1);
  _pfls[i]->setName(pflCpp.c_str());
  MEDprofileRd(fid,pflName,_pfls[i]->getPointer());
}

void MEDFieldFieldGlobs::writeGlobals(med_idt fid, const MEDFileWritable& opt) const throw(INTERP_KERNEL::Exception)
{
  int nbOfPfls=_pfls.size();
  for(int i=0;i<nbOfPfls;i++)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cpy=_pfls[i]->deepCpy();
      cpy->applyLin(1,1,0);
      INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      MEDLoaderBase::safeStrCpy(_pfls[i]->getName().c_str(),MED_NAME_SIZE,pflName,opt.getTooLongStrPolicy());
      MEDprofileWr(fid,pflName,_pfls[i]->getNumberOfTuples(),cpy->getConstPointer());
    }
  //
  int nbOfLocs=_locs.size();
  for(int i=0;i<nbOfLocs;i++)
    _locs[i]->writeLL(fid);
}

void MEDFieldFieldGlobs::loadGlobals(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> profiles=getPflsReallyUsed();
  int sz=profiles.size();
  _pfls.resize(sz);
  for(int i=0;i<sz;i++)
    loadProfileInFile(fid,i,profiles[i].c_str());
  //
  std::vector<std::string> locs=getLocsReallyUsed();
  sz=locs.size();
  _locs.resize(sz);
  for(int i=0;i<sz;i++)
    _locs[i]=MEDFileFieldLoc::New(fid,locs[i].c_str());
}

MEDFieldFieldGlobs::MEDFieldFieldGlobs(const char *fname):_file_name(fname)
{
}

MEDFieldFieldGlobs::MEDFieldFieldGlobs()
{
}

void MEDFieldFieldGlobs::setFileName(const char *fileName)
{
  _file_name=fileName;
}

int MEDFieldFieldGlobs::getNbOfGaussPtPerCell(int locId) const throw(INTERP_KERNEL::Exception)
{
  if(locId<0 || locId>=(int)_locs.size())
    throw INTERP_KERNEL::Exception("MEDFieldFieldGlobs::getNbOfGaussPtPerCell : Invalid localization id !");
  return _locs[locId]->getNbOfGaussPtPerCell();
}

const MEDFileFieldLoc& MEDFieldFieldGlobs::getLocalizationFromId(int locId) const throw(INTERP_KERNEL::Exception)
{
  if(locId<0 || locId>=(int)_locs.size())
    throw INTERP_KERNEL::Exception("MEDFieldFieldGlobs::getLocalizationFromId : Invalid localization id !");
  return *_locs[locId];
}

namespace ParaMEDMEMImpl
{
  class LocFinder
  {
  public:
    LocFinder(const char *loc):_loc(loc) { }
    bool operator() (const MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc>& loc) { return loc->isName(_loc); }
  private:
    const char *_loc;
  };

  class PflFinder
  {
  public:
    PflFinder(const std::string& pfl):_pfl(pfl) { }
    bool operator() (const MEDCouplingAutoRefCountObjectPtr<DataArrayInt>& pfl) { return _pfl==pfl->getName(); }
  private:
    const std::string& _pfl;
  };
}

int MEDFieldFieldGlobs::getLocalizationId(const char *loc) const throw(INTERP_KERNEL::Exception)
{
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> >::const_iterator it=std::find_if(_locs.begin(),_locs.end(),ParaMEDMEMImpl::LocFinder(loc));
  if(it==_locs.end())
    {
      std::ostringstream oss; oss << "MEDFieldFieldGlobs::getLocalisationId : no such localisation name : \"" << loc << "\" Possible localizations are : ";
      for(it=_locs.begin();it!=_locs.end();it++)
        oss << "\"" << (*it)->getName() << "\", ";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return std::distance(_locs.begin(),it);
}

const DataArrayInt *MEDFieldFieldGlobs::getProfile(const std::string& pflName) const throw(INTERP_KERNEL::Exception)
{
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> >::const_iterator it=std::find_if(_pfls.begin(),_pfls.end(),ParaMEDMEMImpl::PflFinder(pflName));
  if(it==_pfls.end())
    {
      std::ostringstream oss; oss << "MEDFieldFieldGlobs::getProfile: no such profile name : \"" << pflName << "\" Possible profiles are : ";
      for(it=_pfls.begin();it!=_pfls.end();it++)
        oss << "\"" << (*it)->getName() << "\", ";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return *it;
}

std::vector<std::string> MEDFieldFieldGlobs::getPfls() const
{
  int sz=_pfls.size();
  std::vector<std::string> ret(sz);
  for(int i=0;i<sz;i++)
    ret[i]=_pfls[i]->getName();
  return ret;
}

std::vector<std::string> MEDFieldFieldGlobs::getLocs() const
{
  int sz=_locs.size();
  std::vector<std::string> ret(sz);
  for(int i=0;i<sz;i++)
    ret[i]=_locs[i]->getName();
  return ret;
}

void MEDFieldFieldGlobs::appendProfile(DataArrayInt *pfl) throw(INTERP_KERNEL::Exception)
{
  std::string name(pfl->getName());
  if(name.empty())
    throw INTERP_KERNEL::Exception("MEDFieldFieldGlobs::appendProfile : unsupported profiles with no name !");
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> >::const_iterator it=_pfls.begin();it!=_pfls.end();it++)
    if(name==(*it)->getName())
      {
        if(!pfl->isEqual(*(*it)))
          {
            std::ostringstream oss; oss << "MEDFieldFieldGlobs::appendProfile : profile \"" << name << "\" already exists and is different from existing !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
  pfl->incrRef();
  _pfls.push_back(pfl);
}

void MEDFieldFieldGlobs::appendLoc(const char *locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w) throw(INTERP_KERNEL::Exception)
{
  std::string name(locName);
  if(name.empty())
    throw INTERP_KERNEL::Exception("MEDFieldFieldGlobs::appendLoc : unsupported localizations with no name !");
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> obj=MEDFileFieldLoc::New(locName,geoType,refCoo,gsCoo,w);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> >::const_iterator it=_locs.begin();it!=_locs.end();it++)
    if((*it)->isName(locName))
      {
        if(!(*it)->isEqual(*obj,1e-12))
          {
            std::ostringstream oss; oss << "MEDFieldFieldGlobs::appendLoc : localization \"" << name << "\" already exists and is different from existing !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
  _locs.push_back(obj);
}

/*!
 * This method returns the max dimension of 'this'.
 * This method returns -2 if 'this' is empty, -1 if only nodes are defined.
 */
int MEDFileField1TSWithoutDAS::getDimension() const
{
  int ret=-2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    (*it)->getDimension(ret);
  return ret;
}

void MEDFileField1TSWithoutDAS::CheckMeshDimRel(int meshDimRelToMax) throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMax>0)
    throw INTERP_KERNEL::Exception("CheckMeshDimRel : This is a meshDimRel not a meshDimRelExt ! So value should be <=0 !");
}

std::vector<int> MEDFileField1TSWithoutDAS::CheckSBTMesh(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception)
{
  //
  std::set<INTERP_KERNEL::NormalizedCellType> geoTypes=mesh->getAllGeoTypes();
  int nbOfTypes=geoTypes.size();
  std::vector<int> code(3*nbOfTypes);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr1=DataArrayInt::New();
  arr1->alloc(nbOfTypes,1);
  int *arrPtr=arr1->getPointer();
  std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator it=geoTypes.begin();
  for(int i=0;i<nbOfTypes;i++,it++)
    arrPtr[i]=std::distance(typmai2,std::find(typmai2,typmai2+MED_N_CELL_FIXED_GEO,*it));
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr2=arr1->checkAndPreparePermutation();
  const int *arrPtr2=arr2->getConstPointer();
  int i=0;
  for(it=geoTypes.begin();it!=geoTypes.end();it++,i++)
    {
      int pos=arrPtr2[i];
      int nbCells=mesh->getNumberOfCellsWithType(*it);
      code[3*pos]=(int)(*it);
      code[3*pos+1]=nbCells;
      code[3*pos+2]=-1;//no profiles
    }
  std::vector<const DataArrayInt *> idsPerType;//no profiles
  DataArrayInt *da=mesh->checkTypeConsistencyAndContig(code,idsPerType);
  if(da)
    {
      da->decrRef();
      throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::CheckSBTMesh : underlying mesh is not sorted by type as MED file expects !");
    }
  return code;
}

MEDFileField1TSWithoutDAS *MEDFileField1TSWithoutDAS::New(const char *fieldName, int csit, int iteration, int order, const std::vector<std::string>& infos)
{
  return new MEDFileField1TSWithoutDAS(fieldName,csit,iteration,order,infos);
}

void MEDFileField1TSWithoutDAS::copyTinyInfoFrom(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception)
{
  _name=field->getName();
  if(_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::copyTinyInfoFrom : unsupported fields with no name in MED file !");
  const DataArrayDouble *arr=field->getArray();
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::copyTinyInfoFrom : no array set !");
  _dt=field->getTime(_iteration,_order);
  _infos=arr->getInfoOnComponent();
}

std::string MEDFileField1TSWithoutDAS::getMeshName() const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::getMeshName : No field set !");
  return _field_per_mesh[0]->getMeshName();
}

int MEDFileField1TSWithoutDAS::getMeshIteration() const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::getMeshIteration : No field set !");
  return _field_per_mesh[0]->getMeshIteration();
}

int MEDFileField1TSWithoutDAS::getMeshOrder() const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::getMeshOrder : No field set !");
  return _field_per_mesh[0]->getMeshOrder();
}

bool MEDFileField1TSWithoutDAS::isDealingTS(int iteration, int order) const
{
  return iteration==_iteration && order==_order;
}

std::pair<int,int> MEDFileField1TSWithoutDAS::getDtIt() const
{
  std::pair<int,int> p;
  fillIteration(p);
  return p;
}

void MEDFileField1TSWithoutDAS::fillIteration(std::pair<int,int>& p) const
{
  p.first=_iteration;
  p.second=_order;
}

void MEDFileField1TSWithoutDAS::finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  med_int numdt,numit;
  med_float dt;
  med_int nmesh;
  med_bool localMesh;
  med_int meshnumdt,meshnumit;
  INTERP_KERNEL::AutoPtr<char> meshName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  MEDfieldComputingStepInfo(fid,_name.c_str(),_csit,&numdt,&numit,&_dt);
  MEDfield23ComputingStepMeshInfo(fid,_name.c_str(),_csit,&numdt,&numit,&dt,&nmesh,meshName,&localMesh,&meshnumdt,&meshnumit);
  if(_iteration!=numdt || _order!=numit)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::finishLoading : unexpected exception internal error !");
  _field_per_mesh.resize(nmesh);
  for(int i=0;i<nmesh;i++)
    _field_per_mesh[i]=MEDFileFieldPerMesh::New(this,i+1,meshnumdt,meshnumit);
  for(int i=0;i<nmesh;i++)
    _field_per_mesh[i]->finishLoading(fid);
}

std::vector<std::string> MEDFileField1TSWithoutDAS::getPflsReallyUsed2() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsed();
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

std::vector<std::string> MEDFileField1TSWithoutDAS::getLocsReallyUsed2() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsed();
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

void MEDFileField1TSWithoutDAS::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::writeLL : empty field !");
  if(_field_per_mesh.size()>1)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::writeLL : In MED3.0 mode in writting mode only ONE underlying mesh supported !");
  _field_per_mesh[0]->copyOptionsFrom(*this);
  _field_per_mesh[0]->writeLL(fid);
}

/*!
 * SBT means Sort By Type.
 * This method is the most basic method to assign field in this. Basic in sense that no renumbering is done. Underlying mesh in 'field' is globaly ignored except for type contiguous check.
 */
void MEDFileField1TSWithoutDAS::setFieldNoProfileSBT(const MEDCouplingFieldDouble *field, MEDFieldFieldGlobs& glob) throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingMesh *mesh=field->getMesh();
  //
  TypeOfField type=field->getTypeOfField();
  std::vector<DataArrayInt *> dummy;
  copyTinyInfoFrom(field);
  if(type!=ON_NODES)
    {
      std::vector<int> code=MEDFileField1TSWithoutDAS::CheckSBTMesh(mesh);
      //
      int pos=addNewEntryIfNecessary(mesh);
      _field_per_mesh[pos]->assignFieldProfile(code,dummy,dummy,field,glob);
    }
  else
    {
      int pos=addNewEntryIfNecessary(mesh);
      _field_per_mesh[pos]->assignNodeFieldNoProfile(field,glob);
    }
}

/*!
 * Generalization of MEDFileField1TSWithoutDAS::setFieldNoProfileSBT method.
 */
void MEDFileField1TSWithoutDAS::setFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFieldFieldGlobs& glob) throw(INTERP_KERNEL::Exception)
{
  TypeOfField type=field->getTypeOfField();
  copyTinyInfoFrom(field);
  if(type!=ON_NODES)
    {
      std::vector<int> code;
      std::vector<DataArrayInt *> globIdsPerType;
      std::vector<DataArrayInt *> idsPerType;
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=mesh->getGenMeshAtLevel(meshDimRelToMax);
      m->splitProfilePerType(profile,code,globIdsPerType,idsPerType);
      //
      std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > globIdsPerType2(globIdsPerType.size());
      for(std::size_t i=0;i<globIdsPerType.size();i++)
        globIdsPerType2[i]=globIdsPerType[i];
      std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > idsPerType2(idsPerType.size());
      for(std::size_t i=0;i<idsPerType.size();i++)
        idsPerType2[i]=idsPerType[i];
      //
      int pos=addNewEntryIfNecessary(m);
      _field_per_mesh[pos]->assignFieldProfile(code,globIdsPerType,idsPerType,field,glob);
    }
  else
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::setFieldProfile : not implemented yet !");
}

MEDCouplingFieldDouble *MEDFileField1TSWithoutDAS::getFieldAtLevel(TypeOfField type, int meshDimRelToMax, const char *mName, int renumPol, const MEDFieldFieldGlobs *glob) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm;
  if(mName==0)
    mm=MEDFileMesh::New(glob->getFileName(),getMeshName().c_str(),getMeshIteration(),getMeshOrder());
  else
    mm=MEDFileMesh::New(glob->getFileName(),mName,-1,-1);
  return MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,glob,mm);
}

MEDCouplingFieldDouble *MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol, const MEDFieldFieldGlobs *glob, const MEDFileMesh *mesh) const throw(INTERP_KERNEL::Exception)
{
  CheckMeshDimRel(meshDimRelToMax);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=mesh->getGenMeshAtLevel(meshDimRelToMax,false);
  const DataArrayInt *d=mesh->getNumberFieldAtLevel(meshDimRelToMax);
  const DataArrayInt *e=mesh->getNumberFieldAtLevel(1);
  return MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel(type,renumPol,glob,m,d,e);
}

MEDCouplingFieldDouble *MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel(TypeOfField type, int renumPol, const MEDFieldFieldGlobs *glob, const MEDCouplingMesh *mesh, const DataArrayInt *cellRenum, const DataArrayInt *nodeRenum) const throw(INTERP_KERNEL::Exception)
{
  static const char msg1[]="MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel : request for a renumbered field following mesh numbering whereas it is a profile field !";
  int dimRequested=mesh->getMeshDimension();
  int meshId=getMeshIdFromMeshName(mesh->getName());
  bool isPfl=false;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=_field_per_mesh[meshId]->getFieldOnMeshAtLevel(type,dimRequested,glob,mesh,isPfl);
  switch(renumPol)
    {
    case 0:
      {
        //no need to test _field_per_mesh.empty() because geMeshName has already done it
        ret->incrRef();
        return ret;
      }
    case 3:
    case 1:
      {
        if(isPfl)
          throw INTERP_KERNEL::Exception(msg1);
        //no need to test _field_per_mesh.empty() because geMeshName has already done it
        if(cellRenum)
          {
            if(cellRenum->getNbOfElems()!=mesh->getNumberOfCells())
              {
                std::ostringstream oss; oss << "MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel : Request of simple renumbering but it seems that underlying mesh \"" << mesh->getName() << "\" of requested field ";
                oss << "\"" << getName() << "\" has partial renumbering (some geotype has no renumber) !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
            ret->renumberCells(cellRenum->getConstPointer(),true);
          }
        if(renumPol==1)
          {
            ret->incrRef();
            return ret;
          }
      }
    case 2:
      {
        //no need to test _field_per_mesh.empty() because geMeshName has already done it
        if(isPfl)
          throw INTERP_KERNEL::Exception(msg1);
        if(nodeRenum)
          {
            if(nodeRenum->getNbOfElems()!=mesh->getNumberOfNodes())
              {
                std::ostringstream oss; oss << "MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel : Request of simple renumbering but it seems that underlying mesh \"" << mesh->getName() << "\" of requested field ";
                oss << "\"" << getName() << "\" not defined on all nodes !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
            MEDCouplingAutoRefCountObjectPtr<DataArrayInt> nodeRenumSafe=nodeRenum->checkAndPreparePermutation();
            ret->renumberNodes(nodeRenumSafe->getConstPointer());
          }
        ret->incrRef();
        return ret;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel : unsupported renum policy ! Dealing with policy 0 1 2 and 3 !");
    }
}

MEDFileField1TSWithoutDAS::MEDFileField1TSWithoutDAS(const char *fieldName, int csit, int iteration, int order,
                                                     const std::vector<std::string>& infos):_name(fieldName),_infos(infos),_csit(csit),_iteration(iteration),_order(order)
{
}

MEDFileField1TSWithoutDAS::MEDFileField1TSWithoutDAS()
{
}

int MEDFileField1TSWithoutDAS::addNewEntryIfNecessary(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception)
{
  std::string tmp(mesh->getName());
  if(tmp.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::addNewEntryIfNecessary : empty mesh name ! unsupported by MED file !");
  std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();
  int i=0;
  for(;it!=_field_per_mesh.end();it++,i++)
    {
      if((*it)->getMeshName()==tmp)
        return i;
    }
  int sz=_field_per_mesh.size();
  _field_per_mesh.resize(sz+1);
  _field_per_mesh[sz]=MEDFileFieldPerMesh::New(this,mesh);
  return sz;
}

int MEDFileField1TSWithoutDAS::getMeshIdFromMeshName(const char *mName) const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::getMeshIdFromMeshName : No field set !");
  std::string mName2(mName);
  int ret=0;
  std::vector<std::string> msg;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++,ret++)
    if(mName2==(*it)->getMeshName())
      return ret;
    else
      msg.push_back((*it)->getMeshName());
  std::ostringstream oss; oss << "MEDFileField1TSWithoutDAS::getMeshIdFromMeshName : No such mesh \"" << mName2 << "\" as underlying mesh of field \"" << getName() << "\" !\n";
  oss << "Possible meshes are : ";
  for(std::vector<std::string>::const_iterator it2=msg.begin();it2!=msg.end();it2++)
    oss << "\"" << (*it2) << "\" ";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

MEDFileField1TS *MEDFileField1TS::New(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileField1TS(fileName,fieldName,iteration,order);
}

MEDFileField1TS *MEDFileField1TS::New()
{
  return new MEDFileField1TS;
}

void MEDFileField1TS::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,medmod);
  int nbComp=_infos.size();
  INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  for(int i=0;i<nbComp;i++)
    {
      std::string info=_infos[i];
      std::string c,u;
      MEDLoaderBase::splitIntoNameAndUnit(info,c,u);
      MEDLoaderBase::safeStrCpy2(c.c_str(),MED_SNAME_SIZE-1,comp+i*MED_SNAME_SIZE,_too_long_str);
      MEDLoaderBase::safeStrCpy2(u.c_str(),MED_SNAME_SIZE-1,unit+i*MED_SNAME_SIZE,_too_long_str);
    }
  if(_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TS::write : MED file does not accept field with empty name !");
  MEDfieldCr(fid,_name.c_str(),MED_FLOAT64,nbComp,comp,unit,getDtUnit().c_str(),getMeshName().c_str());
  writeGlobals(fid,*this);
  writeLL(fid);
}

MEDFileField1TS::MEDFileField1TS(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
try:MEDFileField1TSWithoutDAS(fieldName,-1,iteration,order,std::vector<std::string>()),MEDFieldFieldGlobs(fileName)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  int nbFields=MEDnField(fid);
  med_field_type typcha;
  bool found=false;
  std::vector<std::string> fns(nbFields);
  int nbOfStep2=-1;
  for(int i=0;i<nbFields && !found;i++)
    {
      int ncomp=MEDfieldnComponent(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> dtunit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> nomMaa=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      med_bool localMesh;
      int nbOfStep;
      MEDfieldInfo(fid,i+1,nomcha,nomMaa,&localMesh,&typcha,comp,unit,dtunit,&nbOfStep);
      std::string tmp(nomcha);
      fns[i]=tmp;
      found=(tmp==fieldName);
      if(found)
        {
          nbOfStep2=nbOfStep;
          std::string mname=MEDLoaderBase::buildStringFromFortran(nomMaa,MED_NAME_SIZE);
          _infos.resize(ncomp);
          for(int j=0;j<ncomp;j++)
            _infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+j*MED_SNAME_SIZE,MED_SNAME_SIZE,(char *)unit+j*MED_SNAME_SIZE,MED_SNAME_SIZE);
        }
    }
  if(!found)
    {
      std::ostringstream oss; oss << "No such field '" << fieldName << "' in file '" << fileName << "' ! Available fields are : ";
      std::copy(fns.begin(),fns.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  found=false;
  std::vector< std::pair<int,int> > dtits(nbOfStep2);
  for(int i=0;i<nbOfStep2 && !found;i++)
    {
      med_int numdt,numit;
      med_float dt;
      MEDfieldComputingStepInfo(fid,fieldName,i+1,&numdt,&numit,&dt);
      if(numdt==iteration && numit==order)
        {
          found=true;
          _csit=i+1;
        }
      else
        dtits[i]=std::pair<int,int>(numdt,numit);
    }
  if(!found)
    {
      std::ostringstream oss; oss << "No such iteration (" << iteration << "," << order << ") in existing field '" << fieldName << "' in file '" << fileName << "' ! Available iterations are : ";
      for(std::vector< std::pair<int,int> >::const_iterator iter=dtits.begin();iter!=dtits.end();iter++)
        oss << "(" << (*iter).first << "," << (*iter).second << "), ";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  finishLoading(fid);
  //
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

MEDFileField1TS::MEDFileField1TS()
{
}

std::vector<std::string> MEDFileField1TS::getPflsReallyUsed() const
{
  return getPflsReallyUsed2();
}

std::vector<std::string> MEDFileField1TS::getLocsReallyUsed() const
{
  return getLocsReallyUsed2();
}

/*!
 * This method requests underlying file to perform the job, for mesh reading. If the current instance is not coming from a file and has been constructed from scratch
 * an exception will be thrown. In this case you should use MEDFileField1TS::getFieldOnMeshAtLevel method instead.
 * \b WARNING ! Parameter 'meshDimRelToMax' is relative from read mesh in file that can be different from the field in MED file !
 * It leads that the returned field of this method is always coherent.
 */
MEDCouplingFieldDouble *MEDFileField1TS::getFieldAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  if(_file_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TS::getFieldAtLevel : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtLevel method instead !");
  return MEDFileField1TSWithoutDAS::getFieldAtLevel(type,meshDimRelToMax,0,renumPol,this);
}

/*!
 * \b WARNING, there is a main difference with the two close methods (MEDFileField1TS::getFieldAtLevel and MEDFileField1TS::getFieldOnMeshAtLevel method) !
 * Here the mesh-dimension of 'mesh' is used by this to automatically request the right geoTypes regarding 'type'.
 * If no such element fufilled the deduced dimension and 'type' an exception will be thrown.
 * It leads that the returned field of this method is always coherent.
 */
MEDCouplingFieldDouble *MEDFileField1TS::getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  return MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel(type,renumPol,this,mesh,0,0);
}

/*!
 * This method can be called whatever the mode of instance feeding of this (MED file or from scratch).
 * But the parameter ''meshDimRelToMax' is applyied on 'mesh' (like MEDFileField1TS::getFieldAtLevel does). \b WARNING the dim of 'this' can be different from those in 'mesh' !
 * It leads that the returned field of this method is always coherent.
 */
MEDCouplingFieldDouble *MEDFileField1TS::getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  return MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,this,mesh);
}

/*!
 * This method is identical to MEDFileField1TS::getFieldAtLevel method except that meshName 'mname' should be specified.
 * This method is called "Old" because in MED3 norm a field has only one meshName attached. This method is only here for reader of MED2 files.
 * See MEDFileField1TS::getFieldAtLevel for more information.
 */
MEDCouplingFieldDouble *MEDFileField1TS::getFieldAtLevelOld(TypeOfField type, const char *mname, int meshDimRelToMax, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  if(_file_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TS::getFieldAtLevel : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtLevel method instead !");
  return MEDFileField1TSWithoutDAS::getFieldAtLevel(type,meshDimRelToMax,mname,renumPol,this);
}

/*!
 * SBT means Sort By Type.
 * This method is the most basic method to assign field in this. Basic in sense that no renumbering is done. Underlying mesh in 'field' is globaly ignored except for type contiguous check.
 * 
 */
void MEDFileField1TS::setFieldNoProfileSBT(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception)
{
  _file_name="";
  MEDFileField1TSWithoutDAS::setFieldNoProfileSBT(field,*this);
}

/*!
 * This method is a generalization of MEDFileField1TS::setFieldNoProfileSBT method. Here a profile array is given in input.
 * The support of field 'field' is \b not used by this method, so it can be null or incoherent with field.
 * This method uses input parameters 'mesh', 'meshDimRelToMax' and 'profile' to determine what is really the support of field 'field'. If field is incoherent regarding this deduced support,
 * an exception will be thrown.
 * This method is trying to reduce the size of MEDfile file so profile is created only if it is absolutely necessary. If it is necessary the name of 'profile' will be used to create it in 'this'.
 * In this case, if this profile name is empty an exception will be thrown.
 */
void MEDFileField1TS::setFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception)
{
  _file_name="";
  MEDFileField1TSWithoutDAS::setFieldProfile(field,mesh,meshDimRelToMax,profile,*this);
}

MEDFileFieldMultiTSWithoutDAS *MEDFileFieldMultiTSWithoutDAS::New(med_idt fid, const char *fieldName, int id, const std::vector<std::string>& infos, int nbOfStep) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTSWithoutDAS(fid,fieldName,id,infos,nbOfStep);
}

MEDFileFieldMultiTSWithoutDAS::MEDFileFieldMultiTSWithoutDAS()
{
}

MEDFileFieldMultiTSWithoutDAS::MEDFileFieldMultiTSWithoutDAS(const char *fieldName):_name(fieldName)
{
}

MEDFileFieldMultiTSWithoutDAS::MEDFileFieldMultiTSWithoutDAS(med_idt fid, const char *fieldName, int id, const std::vector<std::string>& infos, int nbOfStep) throw(INTERP_KERNEL::Exception)
try:_name(fieldName),_infos(infos)
{
  finishLoading(fid,nbOfStep);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

const std::vector<std::string>& MEDFileFieldMultiTSWithoutDAS::getInfo() const throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutDAS::getInfos : not time steps !");
  return _time_steps[0]->getInfo();
}

std::string MEDFileFieldMultiTSWithoutDAS::getMeshName() const throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutDAS::getMeshName : not time steps !");
  return _time_steps[0]->getMeshName();
}

std::string MEDFileFieldMultiTSWithoutDAS::getDtUnit() const throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutDAS::getMeshName : not time steps !");
  return _time_steps[0]->getDtUnit();
}

void MEDFileFieldMultiTSWithoutDAS::finishLoading(med_idt fid, int nbPdt) throw(INTERP_KERNEL::Exception)
{
  _time_steps.resize(nbPdt);
  for(int i=0;i<nbPdt;i++)
    {
      std::vector< std::pair<int,int> > ts;
      med_int numdt=0,numo=0;
      med_int meshIt=0,meshOrder=0;
      med_float dt=0.0;
      MEDfieldComputingStepMeshInfo(fid,_name.c_str(),i+1,&numdt,&numo,&dt,&meshIt,&meshOrder);
      _time_steps[i]=MEDFileField1TSWithoutDAS::New(_name.c_str(),i+1,numdt,numo,_infos);
      _time_steps[i]->finishLoading(fid);
    }
}

void MEDFileFieldMultiTSWithoutDAS::copyTinyInfoFrom(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception)
{
  _name=field->getName();
  if(_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutDAS::copyTinyInfoFrom : unsupported fields with no name in MED file !");
  const DataArrayDouble *arr=field->getArray();
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutDAS::copyTinyInfoFrom : no array set !");
  _infos=arr->getInfoOnComponent();
}

void MEDFileFieldMultiTSWithoutDAS::checkCoherencyOfTinyInfo(const MEDCouplingFieldDouble *field) const throw(INTERP_KERNEL::Exception)
{
  static const char MSG[]="MEDFileFieldMultiTSWithoutDAS::checkCoherencyOfTinyInfo : invalid ";
  if(_name!=field->getName())
    {
      std::ostringstream oss; oss << MSG << "name ! should be \"" << _name;
      oss << "\" and it is set in input field to \"" << field->getName() << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const DataArrayDouble *arr=field->getArray();
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutDAS::checkCoherencyOfTinyInfo : no array set !");
  if(_infos!=arr->getInfoOnComponent())
    {
      std::ostringstream oss; oss << MSG << "components ! should be \"";
      std::copy(_infos.begin(),_infos.end(),std::ostream_iterator<std::string>(oss,", "));
      oss << " But compo in input fields are : ";
      std::vector<std::string> tmp=arr->getInfoOnComponent();
      std::copy(tmp.begin(),tmp.end(),std::ostream_iterator<std::string>(oss,", "));
      oss << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void MEDFileFieldMultiTSWithoutDAS::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutDAS::writeLL : no time steps set !");
  std::vector<std::string> infos(getInfo());
  int nbComp=infos.size();
  INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  for(int i=0;i<nbComp;i++)
    {
      std::string info=infos[i];
      std::string c,u;
      MEDLoaderBase::splitIntoNameAndUnit(info,c,u);
      MEDLoaderBase::safeStrCpy2(c.c_str(),MED_SNAME_SIZE-1,comp+i*MED_SNAME_SIZE,_too_long_str);
      MEDLoaderBase::safeStrCpy2(u.c_str(),MED_SNAME_SIZE-1,unit+i*MED_SNAME_SIZE,_too_long_str);
    }
  if(_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutDAS::write : MED file does not accept field with empty name !");
  MEDfieldCr(fid,_name.c_str(),MED_FLOAT64,nbComp,comp,unit,getDtUnit().c_str(),getMeshName().c_str());
  int nbOfTS=_time_steps.size();
  for(int i=0;i<nbOfTS;i++)
    {
      _time_steps[i]->copyOptionsFrom(*this);
      _time_steps[i]->writeLL(fid);
    }
}

int MEDFileFieldMultiTSWithoutDAS::getNumberOfTS() const
{
  return _time_steps.size();
}

std::vector< std::pair<int,int> > MEDFileFieldMultiTSWithoutDAS::getIterations() const
{
  int lgth=_time_steps.size();
  std::vector< std::pair<int,int> > ret(lgth);
  for(int i=0;i<lgth;i++)
    _time_steps[i]->fillIteration(ret[i]);
  return ret;
}

const MEDFileField1TSWithoutDAS& MEDFileFieldMultiTSWithoutDAS::getTimeStepEntry(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutDAS>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it)->isDealingTS(iteration,order))
      return *(*it);
  std::ostringstream oss; oss << "MEDFileFieldMultiTS::getTimeStepEntry : Muli timestep field on time (" << iteration << "," << order << ") does not exist ! Available (iteration,order) are :\n";
  std::vector< std::pair<int,int> > vp=getIterations();
  for(std::vector< std::pair<int,int> >::const_iterator it2=vp.begin();it2!=vp.end();it2++)
    oss << "(" << (*it2).first << "," << (*it2).second << ") ";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

MEDFileField1TSWithoutDAS& MEDFileFieldMultiTSWithoutDAS::getTimeStepEntry(int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutDAS>  >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it)->isDealingTS(iteration,order))
      return *(*it);
  std::ostringstream oss; oss << "MEDFileFieldMultiTS::getTimeStepEntry : Muli timestep field on time (" << iteration << "," << order << ") does not exist ! Available (iteration,order) are :\n";
  std::vector< std::pair<int,int> > vp=getIterations();
  for(std::vector< std::pair<int,int> >::const_iterator it2=vp.begin();it2!=vp.end();it2++)
    oss << "(" << (*it2).first << "," << (*it2).second << ") ";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

std::vector<std::string> MEDFileFieldMultiTSWithoutDAS::getPflsReallyUsed2() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileField1TSWithoutDAS > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
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

std::vector<std::string> MEDFileFieldMultiTSWithoutDAS::getLocsReallyUsed2() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileField1TSWithoutDAS > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsed2();
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

MEDFileFieldMultiTS *MEDFileFieldMultiTS::New()
{
  return new MEDFileFieldMultiTS;
}

MEDFileFieldMultiTS *MEDFileFieldMultiTS::New(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTS(fileName,fieldName);
}

void MEDFileFieldMultiTS::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,medmod);
  writeGlobals(fid,*this);
  writeLL(fid);
}

/*!
 * Performs the job than MEDFileField1TS::getFieldAtLevel except that (iteration,order) couple should be specified !
 * If such couple does not exist an exception is thrown.
 */
MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileField1TSWithoutDAS& myF1TS=getTimeStepEntry(iteration,order);
  return myF1TS.getFieldAtLevel(type,meshDimRelToMax,0,renumPol,this);
}

/*!
 * Performs the job than MEDFileField1TS::getFieldOnMeshAtLevel except that (iteration,order) couple should be specified !
 * If such couple does not exist an exception is thrown.
 */
MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileField1TSWithoutDAS& myF1TS=getTimeStepEntry(iteration,order);
  return myF1TS.getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,this,mesh);
}

/*!
 * Performs the job than MEDFileField1TS::getFieldOnMeshAtLevel except that (iteration,order) couple should be specified !
 * If such couple does not exist an exception is thrown.
 */
MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileField1TSWithoutDAS& myF1TS=getTimeStepEntry(iteration,order);
  return myF1TS.getFieldOnMeshAtLevel(type,renumPol,this,mesh,0,0);
}

MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldAtLevelOld(TypeOfField type, const char *mname, int iteration, int order, int meshDimRelToMax, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileField1TSWithoutDAS& myF1TS=getTimeStepEntry(iteration,order);
  return myF1TS.getFieldAtLevel(type,meshDimRelToMax,mname,renumPol,this);
}

void MEDFileFieldMultiTS::appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutDAS> obj=new MEDFileField1TSWithoutDAS;
      obj->setFieldNoProfileSBT(field,*this);
      copyTinyInfoFrom(field);
      _time_steps.push_back(obj);
    }
  else
    {
      checkCoherencyOfTinyInfo(field);
      MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutDAS> obj=new MEDFileField1TSWithoutDAS;
      obj->setFieldNoProfileSBT(field,*this);
      _time_steps.push_back(obj);
    }
}

void MEDFileFieldMultiTS::appendFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutDAS> obj=new MEDFileField1TSWithoutDAS;
      obj->setFieldProfile(field,mesh,meshDimRelToMax,profile,*this);
      copyTinyInfoFrom(field);
      _time_steps.push_back(obj);
    }
  else
    {
      checkCoherencyOfTinyInfo(field);
      MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutDAS> obj=new MEDFileField1TSWithoutDAS;
      obj->setFieldProfile(field,mesh,meshDimRelToMax,profile,*this);
      _time_steps.push_back(obj);
    }
}

MEDFileFieldMultiTS::MEDFileFieldMultiTS()
{
}

MEDFileFieldMultiTS::MEDFileFieldMultiTS(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
try:MEDFileFieldMultiTSWithoutDAS(fieldName),MEDFieldFieldGlobs(fileName)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  int nbFields=MEDnField(fid);
  med_field_type typcha;
  bool found=false;
  std::vector<std::string> fns(nbFields);
  int nbstep2=-1;
  for(int i=0;i<nbFields && !found;i++)
    {
      int ncomp=MEDfieldnComponent(fid,i+1);
      INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> dtunit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      INTERP_KERNEL::AutoPtr<char> nomMaa=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
      med_bool localMesh;
      int nbOfStep;
      MEDfieldInfo(fid,i+1,nomcha,nomMaa,&localMesh,&typcha,comp,unit,dtunit,&nbOfStep);
      std::string tmp(nomcha);
      fns[i]=tmp;
      found=(tmp==fieldName);
      if(found)
        {
          nbstep2=nbOfStep;
          _infos.resize(ncomp);
          for(int j=0;j<ncomp;j++)
            _infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+j*MED_SNAME_SIZE,MED_SNAME_SIZE,(char *)unit+j*MED_SNAME_SIZE,MED_SNAME_SIZE);
        }
    }
  if(!found)
    {
      std::ostringstream oss; oss << "No such field '" << fieldName << "' in file '" << fileName << "' ! Available fields are : ";
      std::copy(fns.begin(),fns.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  finishLoading(fid,nbstep2);
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

std::vector<std::string> MEDFileFieldMultiTS::getPflsReallyUsed() const
{
  return getPflsReallyUsed2();
}

std::vector<std::string> MEDFileFieldMultiTS::getLocsReallyUsed() const
{
  return getLocsReallyUsed2();
}

MEDFileFields *MEDFileFields::New()
{
  return new MEDFileFields;
}

MEDFileFields *MEDFileFields::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFields(fileName);
}

int MEDFileFields::getNumberOfFields() const
{
  return _fields.size();
}

MEDFileFields::MEDFileFields()
{
}

MEDFileFields::MEDFileFields(const char *fileName) throw(INTERP_KERNEL::Exception)
try:MEDFieldFieldGlobs(fileName)
  {
    MEDFileUtilities::CheckFileForRead(fileName);
    MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
    int nbFields=MEDnField(fid);
    _fields.resize(nbFields);
    med_field_type typcha;
    for(int i=0;i<nbFields;i++)
      {
        int ncomp=MEDfieldnComponent(fid,i+1);
        INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
        INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
        INTERP_KERNEL::AutoPtr<char> dtunit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
        INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
        INTERP_KERNEL::AutoPtr<char> nomMaa=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
        med_bool localMesh;
        int nbOfStep;
        MEDfieldInfo(fid,i+1,nomcha,nomMaa,&localMesh,&typcha,comp,unit,dtunit,&nbOfStep);
        std::vector<std::string> infos(ncomp);
        for(int j=0;j<ncomp;j++)
          infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+j*MED_SNAME_SIZE,MED_SNAME_SIZE,(char *)unit+j*MED_SNAME_SIZE,MED_SNAME_SIZE);
        _fields[i]=MEDFileFieldMultiTSWithoutDAS::New(fid,nomcha,i+1,infos,nbOfStep);
      }
    int nProfil=MEDnProfile(fid);
    _pfls.resize(nProfil);
    for(int i=0;i<nProfil;i++)
      loadProfileInFile(fid,i);
  }
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

void MEDFileFields::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  int i=0;
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,medmod);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutDAS> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      const MEDFileFieldMultiTSWithoutDAS *elt=*it;
      if(!elt)
        {
          std::ostringstream oss; oss << "MEDFileFields::write : at rank #" << i << "/" << _fields.size() << " field is empty !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      elt->copyOptionsFrom(*this);
      elt->writeLL(fid);
    }
}

std::vector<std::string> MEDFileFields::getPflsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldMultiTSWithoutDAS > >::const_iterator it=_fields.begin();it!=_fields.end();it++)
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
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldMultiTSWithoutDAS > >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsed2();
      for(std::vector<std::string>::const_iterator it2=tmp.begin();it2!=tmp.end();it2++)
        if(ret2.find(*it2)==ret2.end())
          {
            ret.push_back(*it2);
            ret2.insert(*it2);
          }
    }
  return ret;
}

void MEDFileFields::resize(int newSize) throw(INTERP_KERNEL::Exception)
{
  _fields.resize(newSize);
}

void MEDFileFields::pushField(MEDFileFieldMultiTS *field) throw(INTERP_KERNEL::Exception)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileFields::pushMesh : invalid input pointer ! should be different from 0 !");
  field->incrRef();
  _fields.push_back(field);
}

void MEDFileFields::setFieldAtPos(int i, MEDFileFieldMultiTS *field) throw(INTERP_KERNEL::Exception)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileFields::setFieldAtPos : invalid input pointer ! should be different from 0 !");
  if(i>=(int)_fields.size())
    _fields.resize(i+1);
  field->incrRef();
  _fields[i]=field;
}

void MEDFileFields::destroyFieldAtPos(int i) throw(INTERP_KERNEL::Exception)
{
  if(i<0 || i>=(int)_fields.size())
    {
      std::ostringstream oss; oss << "MEDFileFields::destroyMeshAtPos : Invalid given id in input (" << i << ") should be in [0," << _fields.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  _fields.erase(_fields.begin()+i);
}
