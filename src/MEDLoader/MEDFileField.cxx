// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
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

#include "MEDFileField.hxx"
#include "MEDFileMesh.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileUtilities.hxx"

#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldDiscretization.hxx"

#include "InterpKernelAutoPtr.hxx"
#include "CellModel.hxx"

#include <algorithm>
#include <iterator>

extern med_geometry_type typmai[MED_N_CELL_FIXED_GEO];
extern INTERP_KERNEL::NormalizedCellType typmai2[MED_N_CELL_FIXED_GEO];
extern med_geometry_type typmainoeud[1];
extern med_geometry_type typmai3[32];

using namespace ParaMEDMEM;

const char MEDFileField1TSWithoutSDA::TYPE_STR[]="FLOAT64";
const char MEDFileIntField1TSWithoutSDA::TYPE_STR[]="INT32";

MEDFileFieldLoc *MEDFileFieldLoc::New(med_idt fid, const char *locName)
{
  return new MEDFileFieldLoc(fid,locName);
}

MEDFileFieldLoc *MEDFileFieldLoc::New(med_idt fid, int id)
{
  return new MEDFileFieldLoc(fid,id);
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

MEDFileFieldLoc::MEDFileFieldLoc(med_idt fid, int id)
{
  med_geometry_type geotype;
  med_geometry_type sectiongeotype;
  int nsectionmeshcell;
  INTERP_KERNEL::AutoPtr<char> locName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> geointerpname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> sectionmeshname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  MEDlocalizationInfo(fid,id+1,locName,&geotype,&_dim,&_nb_gauss_pt,geointerpname,sectionmeshname,&nsectionmeshcell,&sectiongeotype);
  _name=locName;
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

MEDFileFieldLoc *MEDFileFieldLoc::deepCpy() const
{
  return new MEDFileFieldLoc(*this);
}

std::size_t MEDFileFieldLoc::getHeapMemorySize() const
{
  return (_ref_coo.capacity()+_gs_coo.capacity()+_w.capacity())*sizeof(double)+_name.capacity();
}

void MEDFileFieldLoc::simpleRepr(std::ostream& oss) const
{
  static const char OFF7[]="\n    ";
  oss << "\"" << _name << "\"" << OFF7;
  oss << "GeoType=" << INTERP_KERNEL::CellModel::GetCellModel(_geo_type).getRepr() << OFF7;
  oss << "Dimension=" << _dim << OFF7;
  oss << "Number of Gauss points=" << _nb_gauss_pt << OFF7;
  oss << "Number of nodes per cell=" << _nb_node_per_cell << OFF7;
  oss << "RefCoords="; std::copy(_ref_coo.begin(),_ref_coo.end(),std::ostream_iterator<double>(oss," ")); oss << OFF7;
  oss << "Weights="; std::copy(_w.begin(),_w.end(),std::ostream_iterator<double>(oss," ")); oss << OFF7;
  oss << "GaussPtsCoords="; std::copy(_gs_coo.begin(),_gs_coo.end(),std::ostream_iterator<double>(oss," ")); oss << std::endl;
}

void MEDFileFieldLoc::setName(const char *name)
{
  _name=name;
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

std::string MEDFileFieldLoc::repr() const
{
  std::ostringstream oss; oss.precision(15);
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_geo_type);
  oss << "Localization \"" << _name << "\" :\n" << "  - Geometric Type : " << cm.getRepr();
  oss << "\n  - Dimension : " << _dim << "\n  - Number of gauss points : ";
  oss << _nb_gauss_pt << "\n  - Number of nodes in cell : " << _nb_node_per_cell;
  oss << "\n  - Ref coords are : ";
  int sz=_ref_coo.size();
  if(sz%_dim==0)
    {
      int nbOfTuples=sz/_dim;
      for(int i=0;i<nbOfTuples;i++)
        {
          oss << "(";
          for(int j=0;j<_dim;j++)
            { oss << _ref_coo[i*_dim+j]; if(j!=_dim-1) oss << ", "; }
          oss << ") ";
        }
    }
  else
    std::copy(_ref_coo.begin(),_ref_coo.end(),std::ostream_iterator<double>(oss," "));
  oss << "\n  - Gauss coords in reference element : ";
  sz=_gs_coo.size();
  if(sz%_dim==0)
    {
      int nbOfTuples=sz/_dim;
      for(int i=0;i<nbOfTuples;i++)
        {
          oss << "(";
          for(int j=0;j<_dim;j++)
            { oss << _gs_coo[i*_dim+j]; if(j!=_dim-1) oss << ", "; }
          oss << ") ";
        }
    }
  else
    std::copy(_gs_coo.begin(),_gs_coo.end(),std::ostream_iterator<double>(oss," "));
  oss << "\n  - Weights of Gauss coords are : "; std::copy(_w.begin(),_w.end(),std::ostream_iterator<double>(oss," "));
  return oss.str();
}

void MEDFileFieldPerMeshPerTypePerDisc::assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldDouble *field, const DataArray *arrr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  _type=field->getTypeOfField();
  _start=start;
  switch(_type)
    {
    case ON_CELLS:
      {
        getArray()->setContigPartOfSelectedValues2(_start,arrr,offset,offset+nbOfCells,1);
        _end=_start+nbOfCells;
        _nval=nbOfCells;
        break;
      }
    case ON_GAUSS_NE:
      {
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=field->getDiscretization()->getOffsetArr(field->getMesh());
        const int *arrPtr=arr->getConstPointer();
        getArray()->setContigPartOfSelectedValues2(_start,arrr,arrPtr[offset],arrPtr[offset+nbOfCells],1);
        _end=_start+(arrPtr[offset+nbOfCells]-arrPtr[offset]);
        _nval=nbOfCells;
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
            std::ostringstream oss; oss << "Pfl_" << nasc.getName() << "_" << INTERP_KERNEL::CellModel::GetCellModel(getGeoType()).getRepr() << "_" << _loc_id;
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
        std::ostringstream oss2; oss2 << "Loc_" << nasc.getName() << "_" << INTERP_KERNEL::CellModel::GetCellModel(getGeoType()).getRepr() << "_" << _loc_id;
        _localization=oss2.str();
        getArray()->setContigPartOfSelectedValues(_start,arrr,da4);
        _end=_start+_nval*nbi;
        glob.appendLoc(_localization.c_str(),getGeoType(),gsLoc.getRefCoords(),gsLoc.getGaussCoords(),gsLoc.getWeights());
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::assignFieldNoProfile : not implemented yet for such discretization type of field !");
    }
  start=_end;
}

/*!
 * Leaf method of field with profile assignement. This method is the most general one. No optimization is done here.
 * \param [in] pflName input containing name of profile if any. 0 if no profile (except for GAUSS_PT where a no profile can hide a profile when splitted by loc_id).
 * \param [in] multiTypePfl is the end user profile specified in high level API
 * \param [in] idsInPfl is the selection into the \a multiTypePfl whole profile that corresponds to the current geometric type.
 * \param [in] locIds is the profile needed to be created for MED file format. It can be null if all cells of current geometric type are fetched in \a multiTypePfl.
 *             \b WARNING if not null the MED file profile can be subdivided again in case of Gauss points.
 * \param [in] mesh is the mesh coming from the MEDFileMesh instance in correspondance with the MEDFileField. The mesh inside the \a field is simply ignored.
 */
void MEDFileFieldPerMeshPerTypePerDisc::assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, DataArrayInt *locIds, int nbOfEltsInWholeMesh, const MEDCouplingFieldDouble *field, const DataArray *arrr, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  _profile.clear();
  _type=field->getTypeOfField();
  std::string pflName(multiTypePfl->getName());
  std::ostringstream oss; oss << pflName;
  if(_type!=ON_NODES) { const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(getGeoType()); oss << "_" <<  cm.getRepr(); } else { oss << "_NODE"; }
  if(locIds)
    {
      if(pflName.empty())
        throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::assignFieldProfile : existing profile with empty name !");
      if(_type!=ON_GAUSS_PT)
        {
          locIds->setName(oss.str().c_str());
          glob.appendProfile(locIds);
          _profile=oss.str();
        }
    }
  _start=start;
  switch(_type)
    {
    case ON_NODES:
      {
         _nval=idsInPfl->getNumberOfTuples();
         getArray()->setContigPartOfSelectedValues2(_start,arrr,0,arrr->getNumberOfTuples(),1);
         _end=_start+_nval;
         break;
      }
    case ON_CELLS:
      {
        _nval=idsInPfl->getNumberOfTuples();
        getArray()->setContigPartOfSelectedValues(_start,arrr,idsInPfl);
        _end=_start+_nval;
        break;
      }
    case ON_GAUSS_NE:
      {
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=field->getDiscretization()->getOffsetArr(mesh);
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr2=arr->deltaShiftIndex();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr3=arr2->selectByTupleId(multiTypePfl->begin(),multiTypePfl->end());
        arr3->computeOffsets2();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp=idsInPfl->buildExplicitArrByRanges(arr3);
        int trueNval=tmp->getNumberOfTuples();
        _nval=idsInPfl->getNumberOfTuples();
        getArray()->setContigPartOfSelectedValues(_start,arrr,tmp);
        _end=_start+trueNval;
        break;
      }
    case ON_GAUSS_PT:
      {
        const MEDCouplingFieldDiscretizationGauss *disc2=dynamic_cast<const MEDCouplingFieldDiscretizationGauss *>(field->getDiscretization());
        if(!disc2)
          throw INTERP_KERNEL::Exception("addNewEntryIfNecessaryGauss : invalid call to this method ! Internal Error !");
        const DataArrayInt *da1=disc2->getArrayOfDiscIds();
        const MEDCouplingGaussLocalization& gsLoc=field->getGaussLocalization(_loc_id);
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da2=da1->selectByTupleId(idsInPfl->begin(),idsInPfl->end());
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da3=da2->getIdsEqual(_loc_id);
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da4=idsInPfl->selectByTupleId(da3->begin(),da3->end());
        //
        MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> mesh2=mesh->buildPart(multiTypePfl->begin(),multiTypePfl->end());
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=disc2->getOffsetArr(mesh2);
        //
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp=DataArrayInt::New();
        int trueNval=0;
        for(const int *pt=da4->begin();pt!=da4->end();pt++)
          trueNval+=arr->getIJ(*pt+1,0)-arr->getIJ(*pt,0);
        tmp->alloc(trueNval,1);
        int *tmpPtr=tmp->getPointer();
        for(const int *pt=da4->begin();pt!=da4->end();pt++)
          for(int j=arr->getIJ(*pt,0);j<arr->getIJ(*pt+1,0);j++)
            *tmpPtr++=j;
        //
        _nval=da4->getNumberOfTuples();
        getArray()->setContigPartOfSelectedValues(_start,arrr,tmp);
        _end=_start+trueNval;
        oss << "_loc_" << _loc_id;
        if(locIds)
          {
            MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da5=locIds->selectByTupleId(da3->begin(),da3->end());
            da5->setName(oss.str().c_str());
            glob.appendProfile(da5);
            _profile=oss.str();
          }
        else
          {
            if(da3->getNumberOfTuples()!=nbOfEltsInWholeMesh || !da3->isIdentity())
              {
                da3->setName(oss.str().c_str());
                glob.appendProfile(da3);
                _profile=oss.str();
              }
          }
        std::ostringstream oss2; oss2 << "Loc_" << nasc.getName() << "_" << INTERP_KERNEL::CellModel::GetCellModel(getGeoType()).getRepr() << "_" << _loc_id;
        _localization=oss2.str();
        glob.appendLoc(_localization.c_str(),getGeoType(),gsLoc.getRefCoords(),gsLoc.getGaussCoords(),gsLoc.getWeights());
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::assignFieldProfile : not implemented yet for such discretization type of field !");
    }
  start=_end;
}

void MEDFileFieldPerMeshPerTypePerDisc::assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, const DataArray *arrr, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  _start=start;
  _nval=arrr->getNumberOfTuples();
  getArray()->setContigPartOfSelectedValues2(_start,arrr,0,_nval,1);
  _end=_start+_nval;
  start=_end;
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypePerDisc::NewOnRead(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int profileIt) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMeshPerTypePerDisc(fath,type,profileIt);
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypePerDisc::New(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int locId)
{
  return new MEDFileFieldPerMeshPerTypePerDisc(fath,type,locId,std::string());
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypePerDisc::New(const MEDFileFieldPerMeshPerTypePerDisc& other)
{
  return new MEDFileFieldPerMeshPerTypePerDisc(other);
}

std::size_t MEDFileFieldPerMeshPerTypePerDisc::getHeapMemorySize() const
{
  return _profile.capacity()+_localization.capacity()+5*sizeof(int);
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypePerDisc::deepCpy(MEDFileFieldPerMeshPerType *father) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> ret=new MEDFileFieldPerMeshPerTypePerDisc(*this);
  ret->_father=father;
  return ret.retn();
}

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, TypeOfField atype, int profileIt) throw(INTERP_KERNEL::Exception)
try:_type(atype),_father(fath)
  {
  }
catch(INTERP_KERNEL::Exception& e)
{
  throw e;
}

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int locId, const std::string& dummy):_type(type),_father(fath),_loc_id(locId)
{
}

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc(const MEDFileFieldPerMeshPerTypePerDisc& other):_type(other._type),_father(0),_start(other._start),_end(other._end),_nval(other._nval),_profile(other._profile),_localization(other._localization),_loc_id(other._loc_id),_tmp_work1(other._tmp_work1)
{
}

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc():_type(ON_CELLS),_father(0),_start(-std::numeric_limits<int>::max()),_end(-std::numeric_limits<int>::max()),
                                                                       _nval(-std::numeric_limits<int>::max()),_loc_id(-std::numeric_limits<int>::max())
{
}

const MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerTypePerDisc::getFather() const
{
  return _father;
}

void MEDFileFieldPerMeshPerTypePerDisc::prepareLoading(med_idt fid, int profileIt, int& start, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  INTERP_KERNEL::AutoPtr<char> locname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> pflname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  std::string fieldName=nasc.getName();
  std::string meshName=getMeshName();
  int iteration=getIteration();
  int order=getOrder();
  TypeOfField type=getType();
  INTERP_KERNEL::NormalizedCellType geoType=getGeoType();
  int profilesize,nbi;
  med_geometry_type mgeoti;
  med_entity_type menti=MEDFileFieldPerMeshPerType::ConvertIntoMEDFileType(type,geoType,mgeoti);
  _nval=MEDfieldnValueWithProfile(fid,fieldName.c_str(),iteration,order,menti,mgeoti,profileIt,MED_COMPACT_PFLMODE,
                                  pflname,&profilesize,locname,&nbi);
  _profile=MEDLoaderBase::buildStringFromFortran(pflname,MED_NAME_SIZE);
  _localization=MEDLoaderBase::buildStringFromFortran(locname,MED_NAME_SIZE);
  _start=start;
  _end=start+_nval*nbi;
  start=_end;
  if(type==ON_CELLS && !_localization.empty())
    {
      if(_localization!="MED_GAUSS_ELNO")//For compatibily with MED2.3
        setType(ON_GAUSS_PT);
      else
        {
          setType(ON_GAUSS_NE);
          _localization.clear();
        }
    }
}

void MEDFileFieldPerMeshPerTypePerDisc::finishLoading(med_idt fid, int profileIt, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  std::string fieldName=nasc.getName();
  std::string meshName=getMeshName();
  int iteration=getIteration();
  int order=getOrder();
  TypeOfField type=getType();
  INTERP_KERNEL::NormalizedCellType geoType=getGeoType();
  med_geometry_type mgeoti;
  med_entity_type menti=MEDFileFieldPerMeshPerType::ConvertIntoMEDFileType(type,geoType,mgeoti);
  DataArray *arr=getArray();
  DataArrayDouble *arrD=dynamic_cast<DataArrayDouble *>(arr);
  if(arrD)
    {
      double *startFeeding=arrD->getPointer()+_start*arrD->getNumberOfComponents();
      MEDfieldValueWithProfileRd(fid,fieldName.c_str(),iteration,order,menti,mgeoti,MED_COMPACT_PFLMODE,
                                 _profile.c_str(),MED_FULL_INTERLACE,MED_ALL_CONSTITUENT,reinterpret_cast<unsigned char*>(startFeeding));
      return ;
    }
  DataArrayInt *arrI=dynamic_cast<DataArrayInt *>(arr);
  if(arrI)
    {
      int *startFeeding=arrI->getPointer()+_start*arrI->getNumberOfComponents();
      MEDfieldValueWithProfileRd(fid,fieldName.c_str(),iteration,order,menti,mgeoti,MED_COMPACT_PFLMODE,
                                 _profile.c_str(),MED_FULL_INTERLACE,MED_ALL_CONSTITUENT,reinterpret_cast<unsigned char*>(startFeeding));
      return ;
    }
  throw INTERP_KERNEL::Exception("Error on array reading ! Unrecognized type of field ! Should be in FLOAT64 or INT32 !");
}

/*!
 * Set a \c this->_start **and** \c this->_end keeping the same delta between the two.
 */
void MEDFileFieldPerMeshPerTypePerDisc::setNewStart(int newValueOfStart) throw(INTERP_KERNEL::Exception)
{
  int delta=_end-_start;
  _start=newValueOfStart;
  _end=_start+delta;
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

std::string MEDFileFieldPerMeshPerTypePerDisc::getMeshName() const
{
  return _father->getMeshName();
}

void MEDFileFieldPerMeshPerTypePerDisc::simpleRepr(int bkOffset, std::ostream& oss, int id) const
{
  const char startLine[]="    ## ";
  std::string startLine2(bkOffset,' ');
  startLine2+=startLine;
  MEDCouplingFieldDiscretization *tmp=MEDCouplingFieldDiscretization::New(_type);
  oss << startLine2 << "Localization #" << id << "." << std::endl;
  oss << startLine2 << "  Type=" << tmp->getRepr() << "." << std::endl;
  delete tmp;
  oss << startLine2 << "  This type discretization lies on profile : \"" << _profile << "\" and on the following localization : \"" << _localization << "\"." << std::endl;
  oss << startLine2 << "  This type discretization has " << _end-_start << " tuples (start=" << _start << ", end=" << _end << ")." << std::endl;
  oss << startLine2 << "  This type discretization has " << (_end-_start)/_nval << " integration points." << std::endl;
}

TypeOfField MEDFileFieldPerMeshPerTypePerDisc::getType() const
{
  return _type;
}

void MEDFileFieldPerMeshPerTypePerDisc::fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const throw(INTERP_KERNEL::Exception)
{
  types.insert(_type);
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
  return _end-_start;
}

DataArray *MEDFileFieldPerMeshPerTypePerDisc::getArray()
{
  return _father->getArray();
}

const DataArray *MEDFileFieldPerMeshPerTypePerDisc::getArray() const
{
  const MEDFileFieldPerMeshPerType *fath=_father;
  return fath->getArray();
}

DataArrayDouble *MEDFileFieldPerMeshPerTypePerDisc::getArrayDouble()
{
  return _father->getArrayDouble();
}

const DataArrayDouble *MEDFileFieldPerMeshPerTypePerDisc::getArrayDouble() const
{
  const MEDFileFieldPerMeshPerType *fath=_father;
  return fath->getArrayDouble();
}

const std::vector<std::string>& MEDFileFieldPerMeshPerTypePerDisc::getInfo() const
{
  return _father->getInfo();
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getProfile() const
{
  return _profile;
}

void MEDFileFieldPerMeshPerTypePerDisc::setProfile(const char *newPflName)
{
  _profile=newPflName;
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getLocalization() const
{
  return _localization;
}

void MEDFileFieldPerMeshPerTypePerDisc::setLocalization(const char *newLocName)
{
  _localization=newLocName;
}

void MEDFileFieldPerMeshPerTypePerDisc::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< std::pair<std::vector<std::string>, std::string > >::const_iterator it2=mapOfModif.begin();it2!=mapOfModif.end();it2++)
    {
      if(std::find((*it2).first.begin(),(*it2).first.end(),_profile)!=(*it2).first.end())
        {
          _profile=(*it2).second;
          return;
        }
    }
}

void MEDFileFieldPerMeshPerTypePerDisc::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< std::pair<std::vector<std::string>, std::string > >::const_iterator it2=mapOfModif.begin();it2!=mapOfModif.end();it2++)
    {
      if(std::find((*it2).first.begin(),(*it2).first.end(),_localization)!=(*it2).first.end())
        {
          _localization=(*it2).second;
          return;
        }
    }
}

void MEDFileFieldPerMeshPerTypePerDisc::getFieldAtLevel(TypeOfField type, const MEDFileFieldGlobsReal *glob, std::vector< std::pair<int,int> >& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs, std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const
{
  if(type!=_type)
    return ;
  dads.push_back(std::pair<int,int>(_start,_end));
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

void MEDFileFieldPerMeshPerTypePerDisc::fillValues(int discId, int& startEntryId, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const
{
  entries[startEntryId]=std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int> ,std::pair<int,int> >(std::pair<INTERP_KERNEL::NormalizedCellType,int>(getGeoType(),discId),std::pair<int,int>(_start,_end));
  startEntryId++;
}

void MEDFileFieldPerMeshPerTypePerDisc::writeLL(med_idt fid, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception)
{
  TypeOfField type=getType();
  INTERP_KERNEL::NormalizedCellType geoType=getGeoType();
  med_geometry_type mgeoti;
  med_entity_type menti=MEDFileFieldPerMeshPerType::ConvertIntoMEDFileType(type,geoType,mgeoti);
  const DataArray *arr=getArray();
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::writeLL : no array set !");
  const DataArrayDouble *arrD=dynamic_cast<const DataArrayDouble *>(arr);
  const DataArrayInt *arrI=dynamic_cast<const DataArrayInt *>(arr);
  const unsigned char *locToWrite=0;
  if(arrD)
    locToWrite=reinterpret_cast<const unsigned char *>(arrD->getConstPointer()+_start*arr->getNumberOfComponents());
  else if(arrI)
    locToWrite=reinterpret_cast<const unsigned char *>(arrI->getConstPointer()+_start*arr->getNumberOfComponents());
  else
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::writeLL : not recognized type of values ! Supported are FLOAT64 and INT32 !");
  MEDfieldValueWithProfileWr(fid,nasc.getName().c_str(),getIteration(),getOrder(),getTime(),menti,mgeoti,
                             MED_COMPACT_PFLMODE,_profile.c_str(),_localization.c_str(),MED_FULL_INTERLACE,MED_ALL_CONSTITUENT,_nval,
                             locToWrite);
}

void MEDFileFieldPerMeshPerTypePerDisc::getCoarseData(TypeOfField& type, std::pair<int,int>& dad, std::string& pfl, std::string& loc) const throw(INTERP_KERNEL::Exception)
{
  type=_type;
  pfl=_profile;
  loc=_localization;
  dad.first=_start; dad.second=_end;
}

/*!
 * \param [in] codeOfMesh is of format returned by MEDCouplingUMesh::getDistributionOfTypes. And for each *i* oldCode[3*i+2] gives the position (MEDFileUMesh::PutInThirdComponentOfCodeOffset).
 *             This code corresponds to the distribution of types in the corresponding mesh.
 * \param [out] ptToFill memory zone where the output will be stored.
 * \return the size of data pushed into output param \a ptToFill
 */
int MEDFileFieldPerMeshPerTypePerDisc::fillEltIdsFromCode(int offset, const std::vector<int>& codeOfMesh, const MEDFileFieldGlobsReal& glob, int *ptToFill) const throw(INTERP_KERNEL::Exception)
{
  _loc_id=offset;
  std::ostringstream oss;
  std::size_t nbOfType=codeOfMesh.size()/3;
  int found=-1;
  for(std::size_t i=0;i<nbOfType && found==-1;i++)
    if(getGeoType()==(INTERP_KERNEL::NormalizedCellType)codeOfMesh[3*i])
      found=(int)i;
  if(found==-1)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(getGeoType());
      oss << "MEDFileFieldPerMeshPerTypePerDisc::fillEltIdsFromCode : not found geometric type " << cm.getRepr() << " in the referenced mesh of field !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  int *work=ptToFill;
  if(_profile.empty())
    {
      if(_nval!=codeOfMesh[3*found+1])
        {
          const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(getGeoType());
          oss << "MEDFileFieldPerMeshPerTypePerDisc::fillEltIdsFromCode : for geometric type " << cm.getRepr() << " number of elt ids in mesh is equal to " << _nval;
          oss << " whereas mesh has " << codeOfMesh[3*found+1] << " for this geometric type !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      for(int ii=codeOfMesh[3*found+2];ii<codeOfMesh[3*found+2]+_nval;ii++)
        *work++=ii;
    }
  else
    {
      const DataArrayInt *pfl=glob.getProfile(_profile.c_str());
      if(pfl->getNumberOfTuples()!=_nval)
        {
          const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(getGeoType());
          oss << "MEDFileFieldPerMeshPerTypePerDisc::fillEltIdsFromCode : for geometric type " << cm.getRepr() << ", field is defined on profile \"" << _profile << "\" and size of profile is ";
          oss << _nval;
          oss << pfl->getNumberOfTuples() << " whereas the number of ids is set to " << _nval << " for this geometric type !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      int offset2=codeOfMesh[3*found+2];
      for(const int *pflId=pfl->begin();pflId!=pfl->end();pflId++)
        {
          if(*pflId<codeOfMesh[3*found+1])
            *work++=offset2+*pflId;
        }
    }
  return _nval;
}

int MEDFileFieldPerMeshPerTypePerDisc::fillTupleIds(int *ptToFill) const throw(INTERP_KERNEL::Exception)
{
  for(int i=_start;i<_end;i++)
    *ptToFill++=i;
  return _end-_start;
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

std::vector< std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> > MEDFileFieldPerMeshPerTypePerDisc::SplitPerDiscretization(const std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>& entries)
{
  int id=0;
  std::map<std::pair<std::string,TypeOfField>,int> m;
  std::vector< std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> > ret;
  for(std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>::const_iterator it=entries.begin();it!=entries.end();it++)
    if(m.find(std::pair<std::string,TypeOfField>((*it)->getLocalization(),(*it)->getType()))==m.end())
      m[std::pair<std::string,TypeOfField>((*it)->getLocalization(),(*it)->getType())]=id++;
  ret.resize(id);
  for(std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>::const_iterator it=entries.begin();it!=entries.end();it++)
    ret[m[std::pair<std::string,TypeOfField>((*it)->getLocalization(),(*it)->getType())]].push_back(*it);
  return ret;
}

/*!
 * - \c this->_loc_id mutable attribute is used for elt id in mesh offsets.
 * 
 * \param [in] offset the offset id used to take into account that \a result is not compulsary empty in input
 * \param [in] entriesOnSameDisc some entries **on same localization** if not the result can be invalid. The _start and _end on them are relative to \a arr parameter.
 * \param [in] explicitIdsInMesh ids in mesh of the considered chunk.
 * \param [in] newCode one of the input parameter to explicit the new geo type dispatch (in classical format same than those asked by MEDFileFields::renumberEntitiesLyingOnMesh)
 * \param [in,out] glob if necessary by the method, new profiles can be added to it
 * \param [in,out] arr after the call of this method \a arr is renumbered to be compliant with added entries to \a result.
 * \param [out] result All new entries will be appended on it.
 * \return false if the configuration of renumbering leads to an unnecessary resplit of input \a entriesOnSameDisc. If not true is returned (the most general case !)
 */
bool MEDFileFieldPerMeshPerTypePerDisc::RenumberChunks(int offset, const std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>& entriesOnSameDisc,
                                                       const DataArrayInt *explicitIdsInMesh,
                                                       const std::vector<int>& newCode,
                                                       MEDFileFieldGlobsReal& glob, DataArrayDouble *arr,
                                                       std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >& result)
{
  if(entriesOnSameDisc.empty())
    return false;
  TypeOfField type=entriesOnSameDisc[0]->getType();
  int szEntities=0,szTuples=0;
  for(std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>::const_iterator it=entriesOnSameDisc.begin();it!=entriesOnSameDisc.end();it++)
    { szEntities+=(*it)->_nval; szTuples+=(*it)->_end-(*it)->_start; }
  int nbi=szTuples/szEntities;
  if(szTuples%szEntities!=0)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::RenumberChunks : internal error the splitting into same dicretization failed !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> renumTuples=DataArrayInt::New(); renumTuples->alloc(szTuples,1);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ranges=MEDCouplingUMesh::ComputeRangesFromTypeDistribution(newCode);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > newGeoTypesPerChunk(entriesOnSameDisc.size());
  std::vector< const DataArrayInt * > newGeoTypesPerChunk2(entriesOnSameDisc.size());
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > newGeoTypesPerChunk_bis(entriesOnSameDisc.size());
  std::vector< const DataArrayInt * > newGeoTypesPerChunk3(entriesOnSameDisc.size());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newGeoTypesPerChunk4=DataArrayInt::New(); newGeoTypesPerChunk4->alloc(szEntities,nbi);
  int id=0;
  for(std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>::const_iterator it=entriesOnSameDisc.begin();it!=entriesOnSameDisc.end();it++,id++)
    {
      int startOfEltIdOfChunk=(*it)->_start;
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newEltIds=explicitIdsInMesh->substr(startOfEltIdOfChunk,startOfEltIdOfChunk+(*it)->_nval);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> rangeIdsForChunk=newEltIds->findRangeIdForEachTuple(ranges);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> idsInRrangeForChunk=newEltIds->findIdInRangeForEachTuple(ranges);
      //
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp=rangeIdsForChunk->duplicateEachTupleNTimes(nbi); rangeIdsForChunk->rearrange(nbi);
      newGeoTypesPerChunk4->setPartOfValues1(tmp,(*it)->_tmp_work1-offset,(*it)->_tmp_work1+(*it)->_nval*nbi-offset,1,0,nbi,1);
      //
      newGeoTypesPerChunk[id]=rangeIdsForChunk; newGeoTypesPerChunk2[id]=rangeIdsForChunk;
      newGeoTypesPerChunk_bis[id]=idsInRrangeForChunk; newGeoTypesPerChunk3[id]=idsInRrangeForChunk;
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newGeoTypesEltIdsAllGather=DataArrayInt::Aggregate(newGeoTypesPerChunk2); newGeoTypesPerChunk.clear(); newGeoTypesPerChunk2.clear();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> newGeoTypesEltIdsAllGather2=DataArrayInt::Aggregate(newGeoTypesPerChunk3); newGeoTypesPerChunk_bis.clear(); newGeoTypesPerChunk3.clear();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> diffVals=newGeoTypesEltIdsAllGather->getDifferentValues();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> renumEltIds=newGeoTypesEltIdsAllGather->buildPermArrPerLevel();
  //
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> renumTupleIds=newGeoTypesPerChunk4->buildPermArrPerLevel();
  //
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arrPart=arr->substr(offset,offset+szTuples);
  arrPart->renumberInPlace(renumTupleIds->begin());
  arr->setPartOfValues1(arrPart,offset,offset+szTuples,1,0,arrPart->getNumberOfComponents(),1);
  bool ret=false;
  const int *idIt=diffVals->begin();
  std::list<const MEDFileFieldPerMeshPerTypePerDisc *> li(entriesOnSameDisc.begin(),entriesOnSameDisc.end());
  int offset2=0;
  for(int i=0;i<diffVals->getNumberOfTuples();i++,idIt++)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ids=newGeoTypesEltIdsAllGather->getIdsEqual(*idIt);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> subIds=newGeoTypesEltIdsAllGather2->selectByTupleId(ids->begin(),ids->end());
      int nbEntityElts=subIds->getNumberOfTuples();
      bool ret2;
      MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> eltToAdd=MEDFileFieldPerMeshPerTypePerDisc::
        NewObjectOnSameDiscThanPool(type,(INTERP_KERNEL::NormalizedCellType)newCode[3*(*idIt)],subIds,!subIds->isIdentity() || nbEntityElts!=newCode[3*(*idIt)+1],nbi,
                                    offset+offset2,
                                    li,glob,ret2);
      ret=ret || ret2;
      result.push_back(eltToAdd);
      offset2+=nbEntityElts*nbi;
    }
  ret=ret || li.empty();
  return ret;
}

/*!
 * \param [in] typeF type of field of new chunk
 * \param [in] geoType the geometric type of the chunk
 * \param [in] idsOfMeshElt the entity ids of mesh (cells or nodes) of the new chunk.
 * \param [in] isPfl specifies if a profile is requested regarding size of \a idsOfMeshElt and the number of such entities regarding underlying mesh.
 * \param [in] nbi number of integration points
 * \param [in] offset The offset in the **global array of data**.
 * \param [in,out] entriesOnSameDisc the pool **on the same discretization** inside which it will be attempted to find an existing entry corresponding exactly
 *                 to the new chunk to create.
 * \param [in,out] glob the global shared info that will be requested for existing profiles or to append a new profile if needed.
 * \param [out] notInExisting If false the return newly allocated entry is not coming from \a entriesOnSameDisc. If true the output comes from copy of \a entriesOnSameDisc
 *              and corresponding entry erased from \a entriesOnSameDisc.
 * \return a newly allocated chunk
 */
MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypePerDisc::NewObjectOnSameDiscThanPool(TypeOfField typeF, INTERP_KERNEL::NormalizedCellType geoType, DataArrayInt *idsOfMeshElt,
                                                                                                  bool isPfl, int nbi, int offset,
                                                                                                  std::list< const MEDFileFieldPerMeshPerTypePerDisc *>& entriesOnSameDisc,
                                                                                                  MEDFileFieldGlobsReal& glob,
                                                                                                  bool &notInExisting) throw(INTERP_KERNEL::Exception)
{
  int nbMeshEntities=idsOfMeshElt->getNumberOfTuples();
  std::list< const MEDFileFieldPerMeshPerTypePerDisc *>::iterator it=entriesOnSameDisc.begin();
  for(;it!=entriesOnSameDisc.end();it++)
    {
      if(((INTERP_KERNEL::NormalizedCellType)(*it)->_loc_id)==geoType && (*it)->_nval==nbMeshEntities)
        {
          if(!isPfl)
            {
              if((*it)->_profile.empty())
                break;
              else
                if(!(*it)->_profile.empty())
                  {
                    const DataArrayInt *pfl=glob.getProfile((*it)->_profile.c_str());
                    if(pfl->isEqualWithoutConsideringStr(*idsOfMeshElt))
                      break;
                  }
            }
        }
    }
  if(it==entriesOnSameDisc.end())
    {
      notInExisting=true;
      MEDFileFieldPerMeshPerTypePerDisc *ret=new MEDFileFieldPerMeshPerTypePerDisc;
      ret->_type=typeF;
      ret->_loc_id=(int)geoType;
      ret->_nval=nbMeshEntities;
      ret->_start=offset;
      ret->_end=ret->_start+ret->_nval*nbi;
      if(isPfl)
        {
          idsOfMeshElt->setName(glob.createNewNameOfPfl().c_str());
          glob.appendProfile(idsOfMeshElt);
          ret->_profile=idsOfMeshElt->getName();
        }
      //tony treatment of localization
      return ret;
    }
  else
    {
      notInExisting=false;
      MEDFileFieldPerMeshPerTypePerDisc *ret=MEDFileFieldPerMeshPerTypePerDisc::New(*(*it));
      ret->_loc_id=(int)geoType;
      ret->setNewStart(offset);
      entriesOnSameDisc.erase(it);
      return ret;
    }
  
}

MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerType::NewOnRead(med_idt fid, MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMeshPerType(fid,fath,type,geoType,nasc);
}

MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerType::New(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMeshPerType(fath,geoType);
}

std::size_t MEDFileFieldPerMeshPerType::getHeapMemorySize() const
{
  std::size_t ret=_field_pm_pt_pd.capacity()*sizeof(MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc>);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    ret+=(*it)->getHeapMemorySize();
  return ret;
}

MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerType::deepCpy(MEDFileFieldPerMesh *father) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerType> ret=new MEDFileFieldPerMeshPerType(*this);
  ret->_father=father;
  std::size_t i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++,i++)
    {
      if((const MEDFileFieldPerMeshPerTypePerDisc *)*it)
        ret->_field_pm_pt_pd[i]=(*it)->deepCpy((MEDFileFieldPerMeshPerType *)ret);
    }
  return ret.retn();
}

void MEDFileFieldPerMeshPerType::assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  std::vector<int> pos=addNewEntryIfNecessary(field,offset,nbOfCells);
  for(std::vector<int>::const_iterator it=pos.begin();it!=pos.end();it++)
    _field_pm_pt_pd[*it]->assignFieldNoProfile(start,offset,nbOfCells,field,arr,glob,nasc);
}

/*!
 * This method is the most general one. No optimization is done here.
 * \param [in] multiTypePfl is the end user profile specified in high level API
 * \param [in] idsInPfl is the selection into the \a multiTypePfl whole profile that corresponds to the current geometric type.
 * \param [in] locIds is the profile needed to be created for MED file format. It can be null if all cells of current geometric type are fetched in \a multiTypePfl.
 *             \b WARNING if not null the MED file profile can be subdivided again in case of Gauss points.
 * \param [in] nbOfEltsInWholeMesh nb of elts of type \a this->_geo_type in \b WHOLE mesh
 * \param [in] mesh is the mesh coming from the MEDFileMesh instance in correspondance with the MEDFileField. The mesh inside the \a field is simply ignored.
 */
void MEDFileFieldPerMeshPerType::assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, DataArrayInt *locIds, int nbOfEltsInWholeMesh, const MEDCouplingFieldDouble *field, const DataArray *arr, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  std::vector<int> pos=addNewEntryIfNecessary(field,idsInPfl);
  for(std::vector<int>::const_iterator it=pos.begin();it!=pos.end();it++)
    _field_pm_pt_pd[*it]->assignFieldProfile(start,multiTypePfl,idsInPfl,locIds,nbOfEltsInWholeMesh,field,arr,mesh,glob,nasc);
}

void MEDFileFieldPerMeshPerType::assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  _field_pm_pt_pd.resize(1);
  _field_pm_pt_pd[0]=MEDFileFieldPerMeshPerTypePerDisc::New(this,ON_NODES,-3);
  _field_pm_pt_pd[0]->assignNodeFieldNoProfile(start,field,arr,glob);
}

void MEDFileFieldPerMeshPerType::assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> pfl2=pfl->deepCpy();
  //
  _field_pm_pt_pd.resize(1);
  _field_pm_pt_pd[0]=MEDFileFieldPerMeshPerTypePerDisc::New(this,ON_NODES,-3);
  _field_pm_pt_pd[0]->assignFieldProfile(start,pfl,pfl2,pfl2,-1,field,arr,0,glob,nasc);//mesh is not requested so 0 is send.
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
  if(!da)
    throw INTERP_KERNEL::Exception("addNewEntryIfNecessaryGauss (no profile) : no localization ids per cell array available ! The input Gauss node field is maybe invalid !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da2=da->selectByTupleId2(offset,offset+nbOfCells,1);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> retTmp=da2->getDifferentValues();
  if(retTmp->presenceOfValue(-1))
    throw INTERP_KERNEL::Exception("addNewEntryIfNecessaryGauss : some cells have no dicretization description !");
  std::vector<int> ret(retTmp->begin(),retTmp->end());
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
  if(!da)
    throw INTERP_KERNEL::Exception("addNewEntryIfNecessaryGauss : no localization ids per cell array available ! The input Gauss node field is maybe invalid !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da2=da->selectByTupleIdSafe(subCells->getConstPointer(),subCells->getConstPointer()+subCells->getNumberOfTuples());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> retTmp=da2->getDifferentValues();
  if(retTmp->presenceOfValue(-1))
    throw INTERP_KERNEL::Exception("addNewEntryIfNecessaryGauss : some cells have no dicretization description !");
  std::vector<int> ret(retTmp->begin(),retTmp->end());
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

void MEDFileFieldPerMeshPerType::fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    {
      (*it)->fillTypesOfFieldAvailable(types);
    }
}

void MEDFileFieldPerMeshPerType::fillFieldSplitedByType(std::vector< std::pair<int,int> >& dads, std::vector<TypeOfField>& types, std::vector<std::string>& pfls, std::vector<std::string>& locs) const throw(INTERP_KERNEL::Exception)
{
  int sz=_field_pm_pt_pd.size();
  dads.resize(sz); types.resize(sz); pfls.resize(sz); locs.resize(sz);
  for(int i=0;i<sz;i++)
    {
      _field_pm_pt_pd[i]->getCoarseData(types[i],dads[i],pfls[i],locs[i]);
    }
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

std::string MEDFileFieldPerMeshPerType::getMeshName() const
{
  return _father->getMeshName();
}

void MEDFileFieldPerMeshPerType::simpleRepr(int bkOffset, std::ostream& oss, int id) const
{
  const char startLine[]="  ## ";
  std::string startLine2(bkOffset,' ');
  std::string startLine3(startLine2);
  startLine3+=startLine;
  if(_geo_type!=INTERP_KERNEL::NORM_ERROR)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_geo_type);
      oss << startLine3 << "Entry geometry type #" << id << " is lying on geometry types " << cm.getRepr() << "." << std::endl;
    }
  else
    oss << startLine3 << "Entry geometry type #" << id << " is lying on NODES." << std::endl;
  oss << startLine3 << "Entry is defined on " <<  _field_pm_pt_pd.size() << " localizations." << std::endl;
  int i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++,i++)
    {
      const MEDFileFieldPerMeshPerTypePerDisc *cur=(*it);
      if(cur)
        cur->simpleRepr(bkOffset,oss,i);
      else
        {
          oss << startLine2 << "    ## " << "Localization #" << i << " is empty !" << std::endl;
        }
    }
}

void MEDFileFieldPerMeshPerType::getSizes(int& globalSz, int& nbOfEntries) const
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    {
      globalSz+=(*it)->getNumberOfTuples();
    }
  nbOfEntries+=(int)_field_pm_pt_pd.size();
}

INTERP_KERNEL::NormalizedCellType MEDFileFieldPerMeshPerType::getGeoType() const
{
  return _geo_type;
}


int MEDFileFieldPerMeshPerType::getNumberOfComponents() const
{
  return _father->getNumberOfComponents();
}

DataArray *MEDFileFieldPerMeshPerType::getArray()
{
  return _father->getArray();
}

const DataArray *MEDFileFieldPerMeshPerType::getArray() const
{
  const MEDFileFieldPerMesh *fath=_father;
  return fath->getArray();
}

DataArrayDouble *MEDFileFieldPerMeshPerType::getArrayDouble()
{
  return _father->getArrayDouble();
}

const DataArrayDouble *MEDFileFieldPerMeshPerType::getArrayDouble() const
{
  const MEDFileFieldPerMesh *fath=_father;
  return fath->getArrayDouble();
}

const std::vector<std::string>& MEDFileFieldPerMeshPerType::getInfo() const
{
  return _father->getInfo();
}

std::vector<std::string> MEDFileFieldPerMeshPerType::getPflsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it1=_field_pm_pt_pd.begin();it1!=_field_pm_pt_pd.end();it1++)
    {
      std::string tmp=(*it1)->getProfile();
      if(!tmp.empty())
        if(ret2.find(tmp)==ret2.end())
          {
            ret.push_back(tmp);
            ret2.insert(tmp);
          }
    }
  return ret;
}

std::vector<std::string> MEDFileFieldPerMeshPerType::getLocsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it1=_field_pm_pt_pd.begin();it1!=_field_pm_pt_pd.end();it1++)
    {
      std::string tmp=(*it1)->getLocalization();
      if(!tmp.empty() && tmp!=MED_GAUSS_ELNO)
        if(ret2.find(tmp)==ret2.end())
          {
            ret.push_back(tmp);
            ret2.insert(tmp);
          }
    }
  return ret;
}

std::vector<std::string> MEDFileFieldPerMeshPerType::getPflsReallyUsedMulti() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it1=_field_pm_pt_pd.begin();it1!=_field_pm_pt_pd.end();it1++)
    {
      std::string tmp=(*it1)->getProfile();
      if(!tmp.empty())
        ret.push_back(tmp);
    }
  return ret;
}

std::vector<std::string> MEDFileFieldPerMeshPerType::getLocsReallyUsedMulti() const
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

void MEDFileFieldPerMeshPerType::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::iterator it1=_field_pm_pt_pd.begin();it1!=_field_pm_pt_pd.end();it1++)
    (*it1)->changePflsRefsNamesGen(mapOfModif);
}

void MEDFileFieldPerMeshPerType::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::iterator it1=_field_pm_pt_pd.begin();it1!=_field_pm_pt_pd.end();it1++)
    (*it1)->changeLocsRefsNamesGen(mapOfModif);
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerType::getLeafGivenLocId(int locId) throw(INTERP_KERNEL::Exception)
{
  if(_field_pm_pt_pd.empty())
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_geo_type);
      std::ostringstream oss; oss << "MEDFileFieldPerMeshPerType::getLeafGivenLocId : no localizations for geotype \"" << cm.getRepr() << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(locId>=0 && locId<(int)_field_pm_pt_pd.size())
    return _field_pm_pt_pd[locId];
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_geo_type);
  std::ostringstream oss2; oss2 << "MEDFileFieldPerMeshPerType::getLeafGivenLocId : no such locId available (" << locId;
  oss2 << ") for geometric type \"" << cm.getRepr() << "\" It should be in [0," << _field_pm_pt_pd.size() << ") !";
  throw INTERP_KERNEL::Exception(oss2.str().c_str());
  return static_cast<MEDFileFieldPerMeshPerTypePerDisc*>(0);
}

const MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerType::getLeafGivenLocId(int locId) const throw(INTERP_KERNEL::Exception)
{
  if(_field_pm_pt_pd.empty())
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_geo_type);
      std::ostringstream oss; oss << "MEDFileFieldPerMeshPerType::getLeafGivenLocId : no localizations for geotype \"" << cm.getRepr() << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(locId>=0 && locId<(int)_field_pm_pt_pd.size())
    return _field_pm_pt_pd[locId];
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_geo_type);
  std::ostringstream oss2; oss2 << "MEDFileFieldPerMeshPerType::getLeafGivenLocId : no such locId available (" << locId;
  oss2 << ") for geometric type \"" << cm.getRepr() << "\" It should be in [0," << _field_pm_pt_pd.size() << ") !";
  throw INTERP_KERNEL::Exception(oss2.str().c_str());
  return static_cast<const MEDFileFieldPerMeshPerTypePerDisc*>(0);
}

void MEDFileFieldPerMeshPerType::getFieldAtLevel(int meshDim, TypeOfField type, const MEDFileFieldGlobsReal *glob, std::vector< std::pair<int,int> >& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs, std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const
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

void MEDFileFieldPerMeshPerType::fillValues(int& startEntryId, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const
{
  int i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++,i++)
    {
      (*it)->fillValues(i,startEntryId,entries);
    }
}

void MEDFileFieldPerMeshPerType::setLeaves(const std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerTypePerDisc > >& leaves) throw(INTERP_KERNEL::Exception)
{
  _field_pm_pt_pd=leaves;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    (*it)->setFather(this);
}

MEDFileFieldPerMeshPerType::MEDFileFieldPerMeshPerType(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception):_father(fath),_geo_type(geoType)
{
}

MEDFileFieldPerMeshPerType::MEDFileFieldPerMeshPerType(med_idt fid, MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception):_father(fath),_geo_type(geoType)
{
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> locName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  med_geometry_type mgeoti;
  med_entity_type menti=ConvertIntoMEDFileType(type,geoType,mgeoti);
  int nbProfiles=MEDfieldnProfile(fid,nasc.getName().c_str(),getIteration(),getOrder(),menti,mgeoti,pflName,locName);
  _field_pm_pt_pd.resize(nbProfiles);
  for(int i=0;i<nbProfiles;i++)
    {
      _field_pm_pt_pd[i]=MEDFileFieldPerMeshPerTypePerDisc::NewOnRead(this,type,i+1);
    }
}

void MEDFileFieldPerMeshPerType::prepareLoading(med_idt fid, int &start, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  int pflId=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++,pflId++)
    {
      (*it)->prepareLoading(fid,pflId+1,start,nasc);//tony
    }
}

void MEDFileFieldPerMeshPerType::finishLoading(med_idt fid, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  int pflId=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++,pflId++)
    {
      (*it)->finishLoading(fid,pflId+1,nasc);//tony
    }
}

void MEDFileFieldPerMeshPerType::writeLL(med_idt fid, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    {
      (*it)->copyOptionsFrom(*this);
      (*it)->writeLL(fid,nasc);
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

MEDFileFieldPerMesh *MEDFileFieldPerMesh::NewOnRead(med_idt fid, MEDFileAnyTypeField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMesh(fid,fath,meshCsit,meshIteration,meshOrder,nasc);
}

MEDFileFieldPerMesh *MEDFileFieldPerMesh::New(MEDFileAnyTypeField1TSWithoutSDA *fath, const MEDCouplingMesh *mesh)
{
  return new MEDFileFieldPerMesh(fath,mesh);
}

std::size_t MEDFileFieldPerMesh::getHeapMemorySize() const
{
  std::size_t ret=_mesh_name.capacity()+_field_pm_pt.capacity()*sizeof(MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType >);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    if((const MEDFileFieldPerMeshPerType *)*it)
      ret+=(*it)->getHeapMemorySize();
  return ret;
}

MEDFileFieldPerMesh *MEDFileFieldPerMesh::deepCpy(MEDFileAnyTypeField1TSWithoutSDA *father) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > ret=new MEDFileFieldPerMesh(*this);
  ret->_father=father;
  std::size_t i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++,i++)
    {
      if((const MEDFileFieldPerMeshPerType *)*it)
        ret->_field_pm_pt[i]=(*it)->deepCpy((MEDFileFieldPerMesh *)(ret));
    }
  return ret.retn();
}

void MEDFileFieldPerMesh::simpleRepr(int bkOffset, std::ostream& oss, int id) const
{
  std::string startLine(bkOffset,' ');
  oss << startLine << "## Field part (" << id << ") lying on mesh \"" << _mesh_name << "\", Mesh iteration=" << _mesh_iteration << ". Mesh order=" << _mesh_order << "." << std::endl;
  oss << startLine << "## Field is defined on " << _field_pm_pt.size() << " types." << std::endl;
  int i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++,i++)
    {
      const MEDFileFieldPerMeshPerType *cur=*it;
      if(cur)
        cur->simpleRepr(bkOffset,oss,i);
      else
        {
          oss << startLine << "  ## Entry geometry type #" << i << " is empty !" << std::endl;
        }
    }
}

void MEDFileFieldPerMesh::copyTinyInfoFrom(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception)
{
  _mesh_name=mesh->getName();
  mesh->getTime(_mesh_iteration,_mesh_order);
}

void MEDFileFieldPerMesh::assignFieldNoProfileNoRenum(int& start, const std::vector<int>& code, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  int nbOfTypes=code.size()/3;
  int offset=0;
  for(int i=0;i<nbOfTypes;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)code[3*i];
      int nbOfCells=code[3*i+1];
      int pos=addNewEntryIfNecessary(type);
      _field_pm_pt[pos]->assignFieldNoProfile(start,offset,nbOfCells,field,arr,glob,nasc);
      offset+=nbOfCells;
    }
}

/*!
 * This method is the most general one. No optimization is done here.
 * \param [in] multiTypePfl is the end user profile specified in high level API
 * \param [in] code is the code of \a mesh[multiTypePfl] mesh. It is of size of number of different geometric types into \a mesh[multiTypePfl].
 * \param [in] code2 is the code of the \b WHOLE mesh on the same level. So all types in \a code are in \a code2.
 * \param [in] idsInPflPerType is the selection into the \a multiTypePfl whole profile that corresponds to the given geometric type. This vector is always 3 times smaller than \a code.
 * \param [in] idsPerType is a vector containing the profiles needed to be created for MED file format. \b WARNING these processed MED file profiles can be subdivided again in case of Gauss points.
 * \param [in] mesh is the mesh coming from the MEDFileMesh instance in correspondance with the MEDFileField. The mesh inside the \a field is simply ignored.
 */
void MEDFileFieldPerMesh::assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const std::vector<int>& code, const std::vector<int>& code2, const std::vector<DataArrayInt *>& idsInPflPerType, const std::vector<DataArrayInt *>& idsPerType, const MEDCouplingFieldDouble *field, const DataArray *arr, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  int nbOfTypes=code.size()/3;
  for(int i=0;i<nbOfTypes;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)code[3*i];
      int pos=addNewEntryIfNecessary(type);
      DataArrayInt *pfl=0;
      if(code[3*i+2]!=-1)
        pfl=idsPerType[code[3*i+2]];
      int nbOfTupes2=code2.size()/3;
      int found=0;
      for(;found<nbOfTupes2;found++)
        if(code[3*i]==code2[3*found])
          break;
      if(found==nbOfTupes2)
        throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::assignFieldProfile : internal problem ! Should never happen ! Please report bug to anthony.geay@cea.fr !");
      _field_pm_pt[pos]->assignFieldProfile(start,multiTypePfl,idsInPflPerType[i],pfl,code2[3*found+1],field,arr,mesh,glob,nasc);
    }
}

void MEDFileFieldPerMesh::assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  int pos=addNewEntryIfNecessary(INTERP_KERNEL::NORM_ERROR);
  _field_pm_pt[pos]->assignNodeFieldNoProfile(start,field,arr,glob);
}

void MEDFileFieldPerMesh::assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  int pos=addNewEntryIfNecessary(INTERP_KERNEL::NORM_ERROR);
  _field_pm_pt[pos]->assignNodeFieldProfile(start,pfl,field,arr,glob,nasc);
}

void MEDFileFieldPerMesh::prepareLoading(med_idt fid, int& start, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->prepareLoading(fid,start,nasc);
}

void MEDFileFieldPerMesh::finishLoading(med_idt fid, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->finishLoading(fid,nasc);
}

void MEDFileFieldPerMesh::writeLL(med_idt fid, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception)
{
  int nbOfTypes=_field_pm_pt.size();
  for(int i=0;i<nbOfTypes;i++)
    {
      _field_pm_pt[i]->copyOptionsFrom(*this);
      _field_pm_pt[i]->writeLL(fid,nasc);
    }
}

void MEDFileFieldPerMesh::getDimension(int& dim) const
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->getDimension(dim);
}

void MEDFileFieldPerMesh::fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->fillTypesOfFieldAvailable(types);
}

std::vector< std::vector< std::pair<int,int> > > MEDFileFieldPerMesh::getFieldSplitedByType(std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> > & locs) const throw(INTERP_KERNEL::Exception)
{
  int sz=_field_pm_pt.size();
  std::vector< std::vector<std::pair<int,int> > > ret(sz);
  types.resize(sz); typesF.resize(sz); pfls.resize(sz); locs.resize(sz);
  for(int i=0;i<sz;i++)
    {
      types[i]=_field_pm_pt[i]->getGeoType();
      _field_pm_pt[i]->fillFieldSplitedByType(ret[i],typesF[i],pfls[i],locs[i]);
    }
  return ret;
}

double MEDFileFieldPerMesh::getTime() const
{
  int tmp1,tmp2;
  return _father->getTime(tmp1,tmp2);
}

int MEDFileFieldPerMesh::getIteration() const
{
  return _father->getIteration();
}

int MEDFileFieldPerMesh::getOrder() const
{
  return _father->getOrder();
}

int MEDFileFieldPerMesh::getNumberOfComponents() const
{
  return _father->getNumberOfComponents();
}

DataArray *MEDFileFieldPerMesh::getArray()
{
  if(!_father)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::getArray : no father ! internal error !");
  return _father->getOrCreateAndGetArray();
}

const DataArray *MEDFileFieldPerMesh::getArray() const
{
  if(!_father)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::getArray : no father ! internal error !");
  return _father->getOrCreateAndGetArray();
}

DataArrayDouble *MEDFileFieldPerMesh::getArrayDouble()
{
  MEDFileField1TSWithoutSDA *fatherC=dynamic_cast<MEDFileField1TSWithoutSDA *>(_father);
  if(!fatherC)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::getArrayDouble : Expected to be called on double array !");
  return fatherC->getOrCreateAndGetArrayDouble();
}

const DataArrayDouble *MEDFileFieldPerMesh::getArrayDouble() const
{
  const MEDFileField1TSWithoutSDA *fatherC=dynamic_cast<const MEDFileField1TSWithoutSDA *>(_father);
  if(!fatherC)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::getArrayDouble : Expected to be called on double array !");
  return fatherC->getOrCreateAndGetArrayDouble();
}

const std::vector<std::string>& MEDFileFieldPerMesh::getInfo() const
{
  return _father->getInfo();
}

/*!
 * type,geoTypes,dads,pfls,locs are input parameters. They should have the same size.
 * Before the call of this method 'geoTypes','dads','pfls','locs' must be reorganized so that types in geoTypes are contiguous and ordered following typmai2 array.
 * It returns 2 output vectors :
 * - 'code' of size 3*sz where sz is the number of different values into 'geoTypes'
 * - 'notNullPfls' contains sz2 values that are extracted from 'pfls' in which null profiles have been removed.
 * 'code' and 'notNullPfls' are in MEDCouplingUMesh::checkTypeConsistencyAndContig format.
 */
void MEDFileFieldPerMesh::SortArraysPerType(const MEDFileFieldGlobsReal *glob, TypeOfField type, const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector< std::pair<int,int> >& dads, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& locs, std::vector<int>& code, std::vector<DataArrayInt *>& notNullPfls)
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
      std::vector< std::pair<int,int> > tmpDads(dads.begin()+startZone,dads.begin()+j);
      std::vector<const DataArrayInt *> tmpPfls(pfls.begin()+startZone,pfls.begin()+j);
      std::vector<int> tmpLocs(locs.begin()+startZone,locs.begin()+j);
      code[3*i]=(int)refType;
      std::vector<INTERP_KERNEL::NormalizedCellType> refType2(1,refType);
      code[3*i+1]=ComputeNbOfElems(glob,type,refType2,tmpDads,tmpLocs);
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
 * 'dads' 'geoTypes' and 'locs' are input parameters that should have same size sz. sz should be >=1.
 */
int MEDFileFieldPerMesh::ComputeNbOfElems(const MEDFileFieldGlobsReal *glob, TypeOfField type, const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs) throw(INTERP_KERNEL::Exception)
{
  int sz=dads.size();
  int ret=0;
  for(int i=0;i<sz;i++)
    {
      if(locs[i]==-1)
        {
          if(type!=ON_GAUSS_NE)
            ret+=dads[i].second-dads[i].first;
          else
            {
              const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(geoTypes[i]);
              ret+=(dads[i].second-dads[i].first)/cm.getNumberOfNodes();
            }
        }
      else
        {
          int nbOfGaussPtPerCell=glob->getNbOfGaussPtPerCell(locs[i]);
          ret+=(dads[i].second-dads[i].first)/nbOfGaussPtPerCell;
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

std::vector<std::string> MEDFileFieldPerMesh::getPflsReallyUsedMulti() const
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsedMulti();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
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

std::vector<std::string> MEDFileFieldPerMesh::getLocsReallyUsedMulti() const
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsedMulti();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

bool MEDFileFieldPerMesh::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< std::pair<std::string,std::string> >::const_iterator it=modifTab.begin();it!=modifTab.end();it++)
    {
      if((*it).first==_mesh_name)
        {
          _mesh_name=(*it).second;
          return true;
        }
    }
  return false;
}

bool MEDFileFieldPerMesh::renumberEntitiesLyingOnMesh(const char *meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N,
                                                      MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  if(_mesh_name!=meshName)
    return false;
  std::set<INTERP_KERNEL::NormalizedCellType> typesToKeep;
  for(std::size_t i=0;i<oldCode.size()/3;i++) typesToKeep.insert((INTERP_KERNEL::NormalizedCellType)oldCode[3*i]);
  std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > > entries;
  std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> entriesKept;
  std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> otherEntries;
  getUndergroundDataArrayExt(entries);
  DataArrayDouble *arr=getArrayDouble();
  int sz=0;
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::renumberEntitiesLyingOnMesh : DataArrayDouble storing values of field is null !");
  for(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >::const_iterator it=entries.begin();it!=entries.end();it++)
    {
      if(typesToKeep.find((*it).first.first)!=typesToKeep.end())
        {
          entriesKept.push_back(getLeafGivenTypeAndLocId((*it).first.first,(*it).first.second));
          sz+=(*it).second.second-(*it).second.first;
        }
      else
        otherEntries.push_back(getLeafGivenTypeAndLocId((*it).first.first,(*it).first.second));
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> renumDefrag=DataArrayInt::New(); renumDefrag->alloc(arr->getNumberOfTuples(),1); renumDefrag->fillWithZero();
  ////////////////////
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> explicitIdsOldInMesh=DataArrayInt::New(); explicitIdsOldInMesh->alloc(sz,1);//sz is a majorant of the real size. A realloc will be done after
  int *workI2=explicitIdsOldInMesh->getPointer();
  int sz1=0,sz2=0,sid=1;
  std::vector< std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> > entriesKeptML=MEDFileFieldPerMeshPerTypePerDisc::SplitPerDiscretization(entriesKept);
  // std::vector<int> tupleIdOfStartOfNewChuncksV(entriesKeptML.size());
  for(std::vector< std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> >::const_iterator itL1=entriesKeptML.begin();itL1!=entriesKeptML.end();itL1++,sid++)
    {
      //  tupleIdOfStartOfNewChuncksV[sid-1]=sz2;
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> explicitIdsOldInArr=DataArrayInt::New(); explicitIdsOldInArr->alloc(sz,1);
      int *workI=explicitIdsOldInArr->getPointer();
      for(std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>::const_iterator itL2=(*itL1).begin();itL2!=(*itL1).end();itL2++)
        {
          int delta1=(*itL2)->fillTupleIds(workI); workI+=delta1; sz1+=delta1;
          (*itL2)->setLocId(sz2);
          (*itL2)->_tmp_work1=(*itL2)->getStart();
          int delta2=(*itL2)->fillEltIdsFromCode(sz2,oldCode,glob,workI2); workI2+=delta2; sz2+=delta2;
        }
      renumDefrag->setPartOfValuesSimple3(sid,explicitIdsOldInArr->begin(),explicitIdsOldInArr->end(),0,1,1);
    }
  explicitIdsOldInMesh->reAlloc(sz2);
  int tupleIdOfStartOfNewChuncks=arr->getNumberOfTuples()-sz2;
  ////////////////////
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> permArrDefrag=renumDefrag->buildPermArrPerLevel(); renumDefrag=0;
  // perform redispatching of non concerned MEDFileFieldPerMeshPerTypePerDisc
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> > otherEntriesNew;
  for(std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>::const_iterator it=otherEntries.begin();it!=otherEntries.end();it++)
    {
      otherEntriesNew.push_back(MEDFileFieldPerMeshPerTypePerDisc::New(*(*it)));
      otherEntriesNew.back()->setNewStart(permArrDefrag->getIJ((*it)->getStart(),0));
      otherEntriesNew.back()->setLocId((*it)->getGeoType());
    }
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> > entriesKeptNew;
  std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> entriesKeptNew2;
  for(std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>::const_iterator it=entriesKept.begin();it!=entriesKept.end();it++)
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> elt=MEDFileFieldPerMeshPerTypePerDisc::New(*(*it));
      int newStart=elt->getLocId();
      elt->setLocId((*it)->getGeoType());
      elt->setNewStart(newStart);
      elt->_tmp_work1=permArrDefrag->getIJ(elt->_tmp_work1,0);
      entriesKeptNew.push_back(elt);
      entriesKeptNew2.push_back(elt);
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2=arr->renumber(permArrDefrag->getConstPointer());
  // perform redispatching of concerned MEDFileFieldPerMeshPerTypePerDisc -> values are in arr2
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> explicitIdsNewInMesh=renumO2N->selectByTupleId(explicitIdsOldInMesh->begin(),explicitIdsOldInMesh->end());
  std::vector< std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> > entriesKeptPerDisc=MEDFileFieldPerMeshPerTypePerDisc::SplitPerDiscretization(entriesKeptNew2);
  bool ret=false;
  for(std::vector< std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> >::const_iterator it4=entriesKeptPerDisc.begin();it4!=entriesKeptPerDisc.end();it4++)
    {
      sid=0;
      /*for(std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>::const_iterator itL2=(*it4).begin();itL2!=(*it4).end();itL2++)
        {
          MEDFileFieldPerMeshPerTypePerDisc *curNC=const_cast<MEDFileFieldPerMeshPerTypePerDisc *>(*itL2);
          curNC->setNewStart(permArrDefrag->getIJ((*itL2)->getStart(),0)-tupleIdOfStartOfNewChuncks+tupleIdOfStartOfNewChuncksV[sid]);
          }*/
      ret=MEDFileFieldPerMeshPerTypePerDisc::RenumberChunks(tupleIdOfStartOfNewChuncks,*it4,explicitIdsNewInMesh,newCode,
                                                            glob,arr2,otherEntriesNew) || ret;
    }
  if(!ret)
    return false;
  // Assign new dispatching
  assignNewLeaves(otherEntriesNew);
  arr->cpyFrom(*arr2);
  return true;
}

void MEDFileFieldPerMesh::assignNewLeaves(const std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerTypePerDisc > >& leaves) throw(INTERP_KERNEL::Exception)
{
  std::map<INTERP_KERNEL::NormalizedCellType,std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerTypePerDisc> > > types;
  for( std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerTypePerDisc > >::const_iterator it=leaves.begin();it!=leaves.end();it++)
    types[(INTERP_KERNEL::NormalizedCellType)(*it)->getLocId()].push_back(*it);
  //
  std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > > fieldPmPt(types.size());
  std::map<INTERP_KERNEL::NormalizedCellType,std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerTypePerDisc> > >::const_iterator it1=types.begin();
  std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::iterator it2=fieldPmPt.begin();
  for(;it1!=types.end();it1++,it2++)
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerType> elt=MEDFileFieldPerMeshPerType::New(this,(INTERP_KERNEL::NormalizedCellType)((*it1).second[0]->getLocId()));
      elt->setLeaves((*it1).second);
      *it2=elt;
    }
  _field_pm_pt=fieldPmPt;
}

void MEDFileFieldPerMesh::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->changePflsRefsNamesGen(mapOfModif);
}

void MEDFileFieldPerMesh::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->changeLocsRefsNamesGen(mapOfModif);
}

/*!
 * \param [in] mesh is the whole mesh
 */
MEDCouplingFieldDouble *MEDFileFieldPerMesh::getFieldOnMeshAtLevel(TypeOfField type, const MEDFileFieldGlobsReal *glob, const MEDCouplingMesh *mesh, bool& isPfl, MEDCouplingAutoRefCountObjectPtr<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception)
{
  if(_field_pm_pt.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::getFieldOnMeshAtLevel : no types field set !");
  //
  std::vector< std::pair<int,int> > dads;
  std::vector<const DataArrayInt *> pfls;
  std::vector<DataArrayInt *> notNullPflsPerGeoType;
  std::vector<int> locs,code;
  std::vector<INTERP_KERNEL::NormalizedCellType> geoTypes;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->getFieldAtLevel(mesh->getMeshDimension(),type,glob,dads,pfls,locs,geoTypes);
  // Sort by types
  SortArraysPerType(glob,type,geoTypes,dads,pfls,locs,code,notNullPflsPerGeoType);
  if(code.empty())
    {
      std::ostringstream oss; oss << "MEDFileFieldPerMesh::getFieldOnMeshAtLevel : " << "The field \"" << nasc.getName() << "\" exists but not with such spatial discretization or such dimension specified !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  //
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > notNullPflsPerGeoType2(notNullPflsPerGeoType.begin(),notNullPflsPerGeoType.end());
  std::vector< const DataArrayInt *> notNullPflsPerGeoType3(notNullPflsPerGeoType.begin(),notNullPflsPerGeoType.end());
  if(type!=ON_NODES)
    {
      DataArrayInt *arr=mesh->checkTypeConsistencyAndContig(code,notNullPflsPerGeoType3);
      if(!arr)
        return finishField(type,glob,dads,locs,mesh,isPfl,arrOut,nasc);
      else
        {
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr2(arr);
          return finishField2(type,glob,dads,locs,geoTypes,mesh,arr,isPfl,arrOut,nasc);
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
          return finishField(type,glob,dads,locs,mesh,isPfl,arrOut,nasc);
        }
      else
        return finishFieldNode2(glob,dads,locs,mesh,notNullPflsPerGeoType3[0],isPfl,arrOut,nasc);
    }
}

DataArray *MEDFileFieldPerMesh::getFieldOnMeshAtLevelWithPfl(TypeOfField type, const MEDCouplingMesh *mesh, DataArrayInt *&pfl, const MEDFileFieldGlobsReal *glob, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception)
{
  if(_field_pm_pt.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::getFieldOnMeshAtLevel : no types field set !");
  //
  std::vector<std::pair<int,int> > dads;
  std::vector<const DataArrayInt *> pfls;
  std::vector<DataArrayInt *> notNullPflsPerGeoType;
  std::vector<int> locs,code;
  std::vector<INTERP_KERNEL::NormalizedCellType> geoTypes;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->getFieldAtLevel(mesh->getMeshDimension(),type,glob,dads,pfls,locs,geoTypes);
  // Sort by types
  SortArraysPerType(glob,type,geoTypes,dads,pfls,locs,code,notNullPflsPerGeoType);
  if(code.empty())
    {
      std::ostringstream oss; oss << "MEDFileFieldPerMesh::getFieldOnMeshAtLevelWithPfl : " << "The field \"" << nasc.getName() << "\" exists but not with such spatial discretization or such dimension specified !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > notNullPflsPerGeoType2(notNullPflsPerGeoType.begin(),notNullPflsPerGeoType.end());
  std::vector< const DataArrayInt *> notNullPflsPerGeoType3(notNullPflsPerGeoType.begin(),notNullPflsPerGeoType.end());
  if(type!=ON_NODES)
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=mesh->checkTypeConsistencyAndContig(code,notNullPflsPerGeoType3);
      return finishField4(dads,arr,mesh->getNumberOfCells(),pfl);
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
        }
      return finishField4(dads,code[2]==-1?0:notNullPflsPerGeoType3[0],mesh->getNumberOfNodes(),pfl);
    }
  //
  return 0;
}

void MEDFileFieldPerMesh::getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception)
{
  int globalSz=0;
  int nbOfEntries=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      (*it)->getSizes(globalSz,nbOfEntries);
    }
  entries.resize(nbOfEntries);
  nbOfEntries=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      (*it)->fillValues(nbOfEntries,entries);
    }
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMesh::getLeafGivenTypeAndLocId(INTERP_KERNEL::NormalizedCellType typ, int locId) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      if((*it)->getGeoType()==typ)
        return (*it)->getLeafGivenLocId(locId);
    }
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(typ);
  std::ostringstream oss; oss << "MEDFileFieldPerMesh::getLeafGivenTypeAndLocId : no such geometric type \"" << cm.getRepr() << "\" in this !" << std::endl;
  oss << "Possiblities are : ";
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      const INTERP_KERNEL::CellModel& cm2=INTERP_KERNEL::CellModel::GetCellModel((*it)->getGeoType());
      oss << "\"" << cm2.getRepr() << "\", ";
    }
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

const MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMesh::getLeafGivenTypeAndLocId(INTERP_KERNEL::NormalizedCellType typ, int locId) const throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      if((*it)->getGeoType()==typ)
        return (*it)->getLeafGivenLocId(locId);
    }
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(typ);
  std::ostringstream oss; oss << "MEDFileFieldPerMesh::getLeafGivenTypeAndLocId : no such geometric type \"" << cm.getRepr() << "\" in this !" << std::endl;
  oss << "Possiblities are : ";
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      const INTERP_KERNEL::CellModel& cm2=INTERP_KERNEL::CellModel::GetCellModel((*it)->getGeoType());
      oss << "\"" << cm2.getRepr() << "\", ";
    }
  throw INTERP_KERNEL::Exception(oss.str().c_str());
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

/*!
 * 'dads' and 'locs' input parameters have the same number of elements
 * \param [in] mesh is \b NOT the global mesh, but the possibly reduced mesh. \a mesh parameter will be directly aggregated in the returned field
 */
MEDCouplingFieldDouble *MEDFileFieldPerMesh::finishField(TypeOfField type, const MEDFileFieldGlobsReal *glob,
                                                         const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs,
                                                         const MEDCouplingMesh *mesh, bool& isPfl, MEDCouplingAutoRefCountObjectPtr<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception)
{
  isPfl=false;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(type,ONE_TIME);
  ret->setMesh(mesh); ret->setName(nasc.getName().c_str()); ret->setTime(getTime(),getIteration(),getOrder()); ret->setTimeUnit(nasc.getDtUnit().c_str());
  MEDCouplingAutoRefCountObjectPtr<DataArray> da=getArray()->selectByTupleRanges(dads);
  const std::vector<std::string>& infos=getInfo();
  da->setInfoOnComponents(infos);
  da->setName("");
  if(type==ON_GAUSS_PT)
    {
      int offset=0;
      int nbOfArrs=dads.size();
      for(int i=0;i<nbOfArrs;i++)
        {
          std::vector<std::pair<int,int> > dads2(1,dads[i]); const std::vector<int> locs2(1,locs[i]);
          const std::vector<INTERP_KERNEL::NormalizedCellType> geoTypes2(1,INTERP_KERNEL::NORM_ERROR);
          int nbOfElems=ComputeNbOfElems(glob,type,geoTypes2,dads2,locs2);
          MEDCouplingAutoRefCountObjectPtr<DataArrayInt> di=DataArrayInt::New();
          di->alloc(nbOfElems,1);
          di->iota(offset);
          const MEDFileFieldLoc& fl=glob->getLocalizationFromId(locs[i]);
          ret->setGaussLocalizationOnCells(di->getConstPointer(),di->getConstPointer()+nbOfElems,fl.getRefCoords(),fl.getGaussCoords(),fl.getGaussWeights());
          offset+=nbOfElems;
        }
    }
  arrOut=da;
  return ret.retn();
}

/*!
 * This method is an extension of MEDFileFieldPerMesh::finishField method. It deals with profiles. This method should be called when type is different from ON_NODES.
 * 'dads', 'locs' and 'geoTypes' input parameters have the same number of elements.
 * No check of this is performed. 'da' array contains an array in old2New style to be applyied to mesh to obtain the right support.
 * The order of cells in the returned field is those imposed by the profile.
 * \param [in] mesh is the global mesh.
 */
MEDCouplingFieldDouble *MEDFileFieldPerMesh::finishField2(TypeOfField type, const MEDFileFieldGlobsReal *glob,
                                                          const std::vector<std::pair<int,int> >& dads, const std::vector<int>& locs,
                                                          const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes,
                                                          const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl, MEDCouplingAutoRefCountObjectPtr<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception)
{
  if(da->isIdentity())
    {
      int nbOfTuples=da->getNumberOfTuples();
      if(nbOfTuples==mesh->getNumberOfCells())
        return finishField(type,glob,dads,locs,mesh,isPfl,arrOut,nasc);
    }
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m2=mesh->buildPart(da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
  m2->setName(mesh->getName());
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=finishField(type,glob,dads,locs,m2,isPfl,arrOut,nasc);
  isPfl=true;
  return ret.retn();
}

/*!
 * This method is the complement of MEDFileFieldPerMesh::finishField2 method except that this method works for node profiles.
 */
MEDCouplingFieldDouble *MEDFileFieldPerMesh::finishFieldNode2(const MEDFileFieldGlobsReal *glob,
                                                              const std::vector<std::pair<int,int> >& dads, const std::vector<int>& locs,
                                                              const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl, MEDCouplingAutoRefCountObjectPtr<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception)
{
  if(da->isIdentity())
    {
      int nbOfTuples=da->getNumberOfTuples();
      if(nbOfTuples==mesh->getNumberOfNodes())//No problem for NORM_ERROR because it is in context of node
        return finishField(ON_NODES,glob,dads,locs,mesh,isPfl,arrOut,nasc);
    }
  // Treatment of particular case where nodal field on pfl is requested with a meshDimRelToMax=1.
  const MEDCouplingUMesh *meshu=dynamic_cast<const MEDCouplingUMesh *>(mesh);
  if(meshu)
    {
      if(meshu->getNodalConnectivity()==0)
        {
          MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=finishField(ON_CELLS,glob,dads,locs,mesh,isPfl,arrOut,nasc);
          int nb=da->getNbOfElems();
          const int *ptr=da->getConstPointer();
          MEDCouplingUMesh *meshuc=const_cast<MEDCouplingUMesh *>(meshu);
          meshuc->allocateCells(nb);
          for(int i=0;i<nb;i++)
            meshuc->insertNextCell(INTERP_KERNEL::NORM_POINT1,1,ptr+i);
          meshuc->finishInsertingCells();
          ret->setMesh(meshuc);
          const MEDCouplingFieldDiscretization *disc=ret->getDiscretization();
          if(!disc) throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::finishFieldNode2 : internal error, no discretization on field !");
          disc->checkCoherencyBetween(meshuc,arrOut);
          return ret.retn();
        }
    }
  //
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=finishField(ON_NODES,glob,dads,locs,mesh,isPfl,arrOut,nasc);
  isPfl=true;
  DataArrayInt *arr2=0;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> cellIds=mesh->getCellIdsFullyIncludedInNodeIds(da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> mesh2=mesh->buildPartAndReduceNodes(cellIds->getConstPointer(),cellIds->getConstPointer()+cellIds->getNbOfElems(),arr2);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr3(arr2);
  int nnodes=mesh2->getNumberOfNodes();
  if(nnodes==(int)da->getNbOfElems())
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da3=da->transformWithIndArrR(arr2->begin(),arr2->end());
      arrOut->renumberInPlace(da3->getConstPointer());
      mesh2->setName(mesh->getName());
      ret->setMesh(mesh2);
      return ret.retn();
    }
  else
    {
      std::ostringstream oss; oss << "MEDFileFieldPerMesh::finishFieldNode2 : The field on nodes lies on a node profile so that it is impossible to find a submesh having exactly the same nodes of that profile !!!";
      oss << "So it is impossible to return a well definied MEDCouplingFieldDouble instance on specified mesh on a specified meshDim !" << std::endl;
      oss << "To retrieve correctly such a field you have 3 possibilities :" << std::endl;
      oss << " - use an another meshDim compatible with the field on nodes (MED file does not have such information)" << std::endl;
      oss << " - use an another a meshDimRelToMax equal to 1 -> it will return a mesh with artificial cell POINT1 containing the profile !" << std::endl;
      oss << " - if definitely the node profile has no link with mesh connectivity use MEDFileField1TS::getFieldWithProfile or MEDFileFieldMultiTS::getFieldWithProfile methods instead !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return 0;
}

/*!
 * This method is the most light method of field retrieving.
 */
DataArray *MEDFileFieldPerMesh::finishField4(const std::vector<std::pair<int,int> >& dads, const DataArrayInt *pflIn, int nbOfElems, DataArrayInt *&pflOut) const throw(INTERP_KERNEL::Exception)
{
  if(!pflIn)
    {
      pflOut=DataArrayInt::New();
      pflOut->alloc(nbOfElems,1);
      pflOut->iota(0);
    }
  else
    {
      pflOut=const_cast<DataArrayInt*>(pflIn);
      pflOut->incrRef();
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> safePfl(pflOut);
  MEDCouplingAutoRefCountObjectPtr<DataArray> da=getArray()->selectByTupleRanges(dads);
  const std::vector<std::string>& infos=getInfo();
  int nbOfComp=infos.size();
  for(int i=0;i<nbOfComp;i++)
    da->setInfoOnComponent(i,infos[i].c_str());
  safePfl->incrRef();
  return da.retn();
}

MEDFileFieldPerMesh::MEDFileFieldPerMesh(med_idt fid, MEDFileAnyTypeField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception):_mesh_iteration(meshIteration),_mesh_order(meshOrder),
                                                                                                                                                                                                                 _mesh_csit(meshCsit),_father(fath)
{
  INTERP_KERNEL::AutoPtr<char> meshName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> locName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  for(int i=0;i<MED_N_CELL_FIXED_GEO;i++)
    {
      int nbProfile=MEDfield23nProfile(fid,nasc.getName().c_str(),getIteration(),getOrder(),MED_CELL,typmai[i],_mesh_csit,meshName,pflName,locName);
      if(nbProfile>0)
        {
          _field_pm_pt.push_back(MEDFileFieldPerMeshPerType::NewOnRead(fid,this,ON_CELLS,typmai2[i],nasc));
          _mesh_name=MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1);
        }
      nbProfile=MEDfield23nProfile(fid,nasc.getName().c_str(),getIteration(),getOrder(),MED_NODE_ELEMENT,typmai[i],_mesh_csit,meshName,pflName,locName);
      if(nbProfile>0)
        {
          _field_pm_pt.push_back(MEDFileFieldPerMeshPerType::NewOnRead(fid,this,ON_GAUSS_NE,typmai2[i],nasc));
          _mesh_name=MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1);
        }
    }
  int nbProfile=MEDfield23nProfile(fid,nasc.getName().c_str(),getIteration(),getOrder(),MED_NODE,MED_NONE,_mesh_csit,meshName,pflName,locName);
  if(nbProfile>0)
    {
      _field_pm_pt.push_back(MEDFileFieldPerMeshPerType::NewOnRead(fid,this,ON_NODES,INTERP_KERNEL::NORM_ERROR,nasc));
      _mesh_name=MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1);
    }
}

MEDFileFieldPerMesh::MEDFileFieldPerMesh(MEDFileAnyTypeField1TSWithoutSDA *fath, const MEDCouplingMesh *mesh):_father(fath)
{
  copyTinyInfoFrom(mesh);
}

void MEDFileFieldGlobs::loadProfileInFile(med_idt fid, int id, const char *pflName) throw(INTERP_KERNEL::Exception)
{
  if(id>=(int)_pfls.size())
    _pfls.resize(id+1);
  _pfls[id]=DataArrayInt::New();
  int lgth=MEDprofileSizeByName(fid,pflName);
  _pfls[id]->setName(pflName);
  _pfls[id]->alloc(lgth,1);
  MEDprofileRd(fid,pflName,_pfls[id]->getPointer());
  _pfls[id]->applyLin(1,-1,0);//Converting into C format
}

void MEDFileFieldGlobs::loadProfileInFile(med_idt fid, int i)
{
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  int sz;
  MEDprofileInfo(fid,i+1,pflName,&sz);
  std::string pflCpp=MEDLoaderBase::buildStringFromFortran(pflName,MED_NAME_SIZE);
  if(i>=(int)_pfls.size())
    _pfls.resize(i+1);
  _pfls[i]=DataArrayInt::New();
  _pfls[i]->alloc(sz,1);
  _pfls[i]->setName(pflCpp.c_str());
  MEDprofileRd(fid,pflName,_pfls[i]->getPointer());
  _pfls[i]->applyLin(1,-1,0);//Converting into C format
}

void MEDFileFieldGlobs::writeGlobals(med_idt fid, const MEDFileWritable& opt) const throw(INTERP_KERNEL::Exception)
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

void MEDFileFieldGlobs::appendGlobs(const MEDFileFieldGlobs& other, double eps) throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> pfls=getPfls();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> >::const_iterator it=other._pfls.begin();it!=other._pfls.end();it++)
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
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
    }
  std::vector<std::string> locs=getLocs();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> >::const_iterator it=_locs.begin();it!=_locs.end();it++)
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
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
    }
}

void MEDFileFieldGlobs::checkGlobsPflsPartCoherency(const std::vector<std::string>& pflsUsed) const throw(INTERP_KERNEL::Exception)
{
  for(std::vector<std::string>::const_iterator it=pflsUsed.begin();it!=pflsUsed.end();it++)
    getProfile((*it).c_str());
}

void MEDFileFieldGlobs::checkGlobsLocsPartCoherency(const std::vector<std::string>& locsUsed) const throw(INTERP_KERNEL::Exception)
{
  for(std::vector<std::string>::const_iterator it=locsUsed.begin();it!=locsUsed.end();it++)
    getLocalization((*it).c_str());
}

void MEDFileFieldGlobs::loadGlobals(med_idt fid, const MEDFileFieldGlobsReal& real) throw(INTERP_KERNEL::Exception)
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

void MEDFileFieldGlobs::loadAllGlobals(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  int nProfil=MEDnProfile(fid);
  for(int i=0;i<nProfil;i++)
    loadProfileInFile(fid,i);
  int sz=MEDnLocalization(fid);
  _locs.resize(sz);
  for(int i=0;i<sz;i++)
    {
      _locs[i]=MEDFileFieldLoc::New(fid,i);
    }
}

MEDFileFieldGlobs *MEDFileFieldGlobs::New(const char *fname)
{
  return new MEDFileFieldGlobs(fname);
}

MEDFileFieldGlobs *MEDFileFieldGlobs::New()
{
  return new MEDFileFieldGlobs;
}

std::size_t MEDFileFieldGlobs::getHeapMemorySize() const
{
  std::size_t ret=_file_name.capacity()+_pfls.capacity()*sizeof(MEDCouplingAutoRefCountObjectPtr<DataArrayInt>)+_locs.capacity()*sizeof(MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc>);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> >::const_iterator it=_pfls.begin();it!=_pfls.end();it++)
    ret+=(*it)->getHeapMemorySize();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> >::const_iterator it=_locs.begin();it!=_locs.end();it++)
    ret+=(*it)->getHeapMemorySize();
  return ret;
}

MEDFileFieldGlobs *MEDFileFieldGlobs::deepCpy() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldGlobs> ret=new MEDFileFieldGlobs(*this);
  std::size_t i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> >::const_iterator it=_pfls.begin();it!=_pfls.end();it++,i++)
    {
      if((const DataArrayInt *)*it)
        ret->_pfls[i]=(*it)->deepCpy();
    }
  i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> >::const_iterator it=_locs.begin();it!=_locs.end();it++,i++)
    {
      if((const MEDFileFieldLoc*)*it)
        ret->_locs[i]=(*it)->deepCpy();
    }
  return ret.retn();
}

/*!
 * \throw if a profile in \a pfls in not in \a this.
 * \throw if a localization in \a locs in not in \a this.
 * \sa MEDFileFieldGlobs::deepCpyPart
 */
MEDFileFieldGlobs *MEDFileFieldGlobs::shallowCpyPart(const std::vector<std::string>& pfls, const std::vector<std::string>& locs) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldGlobs> ret=MEDFileFieldGlobs::New();
  for(std::vector<std::string>::const_iterator it1=pfls.begin();it1!=pfls.end();it1++)
    {
      DataArrayInt *pfl=const_cast<DataArrayInt *>(getProfile((*it1).c_str()));
      if(!pfl)
        throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::shallowCpyPart : internal error ! pfl null !");
      pfl->incrRef();
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> pfl2(pfl);
      ret->_pfls.push_back(pfl2);
    }
  for(std::vector<std::string>::const_iterator it2=locs.begin();it2!=locs.end();it2++)
    {
      MEDFileFieldLoc *loc=const_cast<MEDFileFieldLoc *>(&getLocalization((*it2).c_str()));
      if(!loc)
        throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::shallowCpyPart : internal error ! loc null !");
      loc->incrRef();
      MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> loc2(loc);
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
MEDFileFieldGlobs *MEDFileFieldGlobs::deepCpyPart(const std::vector<std::string>& pfls, const std::vector<std::string>& locs) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldGlobs> ret=MEDFileFieldGlobs::New();
  for(std::vector<std::string>::const_iterator it1=pfls.begin();it1!=pfls.end();it1++)
    {
      DataArrayInt *pfl=const_cast<DataArrayInt *>(getProfile((*it1).c_str()));
      if(!pfl)
        throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::deepCpyPart : internal error ! pfl null !");
      ret->_pfls.push_back(pfl->deepCpy());
    }
  for(std::vector<std::string>::const_iterator it2=locs.begin();it2!=locs.end();it2++)
    {
      MEDFileFieldLoc *loc=const_cast<MEDFileFieldLoc *>(&getLocalization((*it2).c_str()));
      if(!loc)
        throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::deepCpyPart : internal error ! loc null !");
      ret->_locs.push_back(loc->deepCpy());
    }
  ret->setFileName(getFileName());
  return ret.retn();
}

MEDFileFieldGlobs::MEDFileFieldGlobs(const char *fname):_file_name(fname)
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

void MEDFileFieldGlobs::setFileName(const char *fileName)
{
  _file_name=fileName;
}

void MEDFileFieldGlobs::changePflsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> >::iterator it=_pfls.begin();it!=_pfls.end();it++)
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

void MEDFileFieldGlobs::changeLocsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> >::iterator it=_locs.begin();it!=_locs.end();it++)
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

int MEDFileFieldGlobs::getNbOfGaussPtPerCell(int locId) const throw(INTERP_KERNEL::Exception)
{
  if(locId<0 || locId>=(int)_locs.size())
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::getNbOfGaussPtPerCell : Invalid localization id !");
  return _locs[locId]->getNbOfGaussPtPerCell();
}

const MEDFileFieldLoc& MEDFileFieldGlobs::getLocalization(const char *locName) const throw(INTERP_KERNEL::Exception)
{
  return getLocalizationFromId(getLocalizationId(locName));
}

const MEDFileFieldLoc& MEDFileFieldGlobs::getLocalizationFromId(int locId) const throw(INTERP_KERNEL::Exception)
{
  if(locId<0 || locId>=(int)_locs.size())
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::getLocalizationFromId : Invalid localization id !");
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

int MEDFileFieldGlobs::getLocalizationId(const char *loc) const throw(INTERP_KERNEL::Exception)
{
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> >::const_iterator it=std::find_if(_locs.begin(),_locs.end(),ParaMEDMEMImpl::LocFinder(loc));
  if(it==_locs.end())
    {
      std::ostringstream oss; oss << "MEDFileFieldGlobs::getLocalisationId : no such localisation name : \"" << loc << "\" Possible localizations are : ";
      for(it=_locs.begin();it!=_locs.end();it++)
        oss << "\"" << (*it)->getName() << "\", ";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return std::distance(_locs.begin(),it);
}

/*!
 * The returned value is never null.
 */
const DataArrayInt *MEDFileFieldGlobs::getProfile(const char *pflName) const throw(INTERP_KERNEL::Exception)
{
  std::string pflNameCpp(pflName);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> >::const_iterator it=std::find_if(_pfls.begin(),_pfls.end(),ParaMEDMEMImpl::PflFinder(pflNameCpp));
  if(it==_pfls.end())
    {
      std::ostringstream oss; oss << "MEDFileFieldGlobs::getProfile: no such profile name : \"" << pflNameCpp << "\" Possible profiles are : ";
      for(it=_pfls.begin();it!=_pfls.end();it++)
        oss << "\"" << (*it)->getName() << "\", ";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return *it;
}

const DataArrayInt *MEDFileFieldGlobs::getProfileFromId(int pflId) const throw(INTERP_KERNEL::Exception)
{
  if(pflId<0 || pflId>=(int)_pfls.size())
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::getProfileFromId : Invalid profile id !");
  return _pfls[pflId];
}

MEDFileFieldLoc& MEDFileFieldGlobs::getLocalizationFromId(int locId) throw(INTERP_KERNEL::Exception)
{
  if(locId<0 || locId>=(int)_locs.size())
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::getLocalizationFromId : Invalid localization id !");
  return *_locs[locId];
}

MEDFileFieldLoc& MEDFileFieldGlobs::getLocalization(const char *locName) throw(INTERP_KERNEL::Exception)
{
  return getLocalizationFromId(getLocalizationId(locName));
}

/*!
 * The returned value is never null.
 */
DataArrayInt *MEDFileFieldGlobs::getProfile(const char *pflName) throw(INTERP_KERNEL::Exception)
{
  std::string pflNameCpp(pflName);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> >::iterator it=std::find_if(_pfls.begin(),_pfls.end(),ParaMEDMEMImpl::PflFinder(pflNameCpp));
  if(it==_pfls.end())
    {
      std::ostringstream oss; oss << "MEDFileFieldGlobs::getProfile: no such profile name : \"" << pflNameCpp << "\" Possible profiles are : ";
      for(it=_pfls.begin();it!=_pfls.end();it++)
        oss << "\"" << (*it)->getName() << "\", ";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return *it;
}

DataArrayInt *MEDFileFieldGlobs::getProfileFromId(int pflId) throw(INTERP_KERNEL::Exception)
{
  if(pflId<0 || pflId>=(int)_pfls.size())
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::getProfileFromId : Invalid profile id !");
  return _pfls[pflId];
}

void MEDFileFieldGlobs::killProfileIds(const std::vector<int>& pflIds) throw(INTERP_KERNEL::Exception)
{
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > newPfls;
  int i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> >::const_iterator it=_pfls.begin();it!=_pfls.end();it++,i++)
    {
      if(std::find(pflIds.begin(),pflIds.end(),i)==pflIds.end())
        newPfls.push_back(*it);
    }
  _pfls=newPfls;
}

void MEDFileFieldGlobs::killLocalizationIds(const std::vector<int>& locIds) throw(INTERP_KERNEL::Exception)
{
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> > newLocs;
  int i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> >::const_iterator it=_locs.begin();it!=_locs.end();it++,i++)
    {
      if(std::find(locIds.begin(),locIds.end(),i)==locIds.end())
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

bool MEDFileFieldGlobs::existsPfl(const char *pflName) const
{
  std::vector<std::string> v=getPfls();
  std::string s(pflName);
  return std::find(v.begin(),v.end(),s)!=v.end();
}

bool MEDFileFieldGlobs::existsLoc(const char *locName) const
{
  std::vector<std::string> v=getLocs();
  std::string s(locName);
  return std::find(v.begin(),v.end(),s)!=v.end();
}

std::vector< std::vector<int> > MEDFileFieldGlobs::whichAreEqualProfiles() const
{
  std::map<int,std::vector<int> > m;
  int i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> >::const_iterator it=_pfls.begin();it!=_pfls.end();it++,i++)
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

void MEDFileFieldGlobs::appendProfile(DataArrayInt *pfl) throw(INTERP_KERNEL::Exception)
{
  std::string name(pfl->getName());
  if(name.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::appendProfile : unsupported profiles with no name !");
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> >::const_iterator it=_pfls.begin();it!=_pfls.end();it++)
    if(name==(*it)->getName())
      {
        if(!pfl->isEqual(*(*it)))
          {
            std::ostringstream oss; oss << "MEDFileFieldGlobs::appendProfile : profile \"" << name << "\" already exists and is different from existing !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
  pfl->incrRef();
  _pfls.push_back(pfl);
}

void MEDFileFieldGlobs::appendLoc(const char *locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w) throw(INTERP_KERNEL::Exception)
{
  std::string name(locName);
  if(name.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobs::appendLoc : unsupported localizations with no name !");
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> obj=MEDFileFieldLoc::New(locName,geoType,refCoo,gsCoo,w);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> >::const_iterator it=_locs.begin();it!=_locs.end();it++)
    if((*it)->isName(locName))
      {
        if(!(*it)->isEqual(*obj,1e-12))
          {
            std::ostringstream oss; oss << "MEDFileFieldGlobs::appendLoc : localization \"" << name << "\" already exists and is different from existing !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
  _locs.push_back(obj);
}

std::string MEDFileFieldGlobs::createNewNameOfPfl() const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> names=getPfls();
  return CreateNewNameNotIn("NewPfl_",names);
}

std::string MEDFileFieldGlobs::createNewNameOfLoc() const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> names=getLocs();
  return CreateNewNameNotIn("NewLoc_",names);
}

std::string MEDFileFieldGlobs::CreateNewNameNotIn(const char *prefix, const std::vector<std::string>& namesToAvoid) throw(INTERP_KERNEL::Exception)
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
MEDFileFieldGlobsReal::MEDFileFieldGlobsReal(const char *fname):_globals(MEDFileFieldGlobs::New(fname))
{
}

/*!
 * Creates an empty MEDFileFieldGlobsReal.
 */
MEDFileFieldGlobsReal::MEDFileFieldGlobsReal():_globals(MEDFileFieldGlobs::New())
{
}

std::size_t MEDFileFieldGlobsReal::getHeapMemorySize() const
{
  std::size_t ret=0;
  if((const MEDFileFieldGlobs *)_globals)
    ret+=_globals->getHeapMemorySize();
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
void MEDFileFieldGlobsReal::shallowCpyOnlyUsedGlobs(const MEDFileFieldGlobsReal& other) throw(INTERP_KERNEL::Exception)
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
void MEDFileFieldGlobsReal::deepCpyOnlyUsedGlobs(const MEDFileFieldGlobsReal& other) throw(INTERP_KERNEL::Exception)
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
    _globals=other._globals->deepCpy();
}

/*!
 * Adds profiles and Gauss points held by another MEDFileFieldGlobsReal to \a this one.
 *  \param [in] other - the MEDFileFieldGlobsReal to copy data from.
 *  \param [in] eps - a precision used to compare Gauss points with same name held by
 *         \a this and \a other MEDFileFieldGlobsReal.
 *  \throw If \a this and \a other hold profiles with equal names but different ids.
 *  \throw If  \a this and \a other hold different Gauss points with equal names.
 */
void MEDFileFieldGlobsReal::appendGlobs(const MEDFileFieldGlobsReal& other, double eps) throw(INTERP_KERNEL::Exception)
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

void MEDFileFieldGlobsReal::checkGlobsCoherency() const throw(INTERP_KERNEL::Exception)
{
  checkGlobsPflsPartCoherency();
  checkGlobsLocsPartCoherency();
}

void MEDFileFieldGlobsReal::checkGlobsPflsPartCoherency() const throw(INTERP_KERNEL::Exception)
{
  contentNotNull()->checkGlobsPflsPartCoherency(getPflsReallyUsed());
}

void MEDFileFieldGlobsReal::checkGlobsLocsPartCoherency() const throw(INTERP_KERNEL::Exception)
{
  contentNotNull()->checkGlobsLocsPartCoherency(getLocsReallyUsed());
}

void MEDFileFieldGlobsReal::loadProfileInFile(med_idt fid, int id, const char *pflName) throw(INTERP_KERNEL::Exception)
{
  contentNotNull()->loadProfileInFile(fid,id,pflName);
}

void MEDFileFieldGlobsReal::loadProfileInFile(med_idt fid, int id)
{
  contentNotNull()->loadProfileInFile(fid,id);
}

void MEDFileFieldGlobsReal::loadGlobals(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  contentNotNull()->loadGlobals(fid,*this);
}

void MEDFileFieldGlobsReal::loadAllGlobals(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  contentNotNull()->loadAllGlobals(fid);
}

void MEDFileFieldGlobsReal::writeGlobals(med_idt fid, const MEDFileWritable& opt) const throw(INTERP_KERNEL::Exception)
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
bool MEDFileFieldGlobsReal::existsPfl(const char *pflName) const
{
  return contentNotNull()->existsPfl(pflName);
}

/*!
 * Checks if the localization with a given name exists.
 *  \param [in] locName - the localization name of interest.
 *  \return bool - \c true if the localization named \a locName exists.
 */
bool MEDFileFieldGlobsReal::existsLoc(const char *locName) const
{
  return contentNotNull()->existsLoc(locName);
}

std::string MEDFileFieldGlobsReal::createNewNameOfPfl() const throw(INTERP_KERNEL::Exception)
{
  return contentNotNull()->createNewNameOfPfl();
}

std::string MEDFileFieldGlobsReal::createNewNameOfLoc() const throw(INTERP_KERNEL::Exception)
{
  return contentNotNull()->createNewNameOfLoc();
}

/*!
 * Sets the name of a MED file.
 *  \param [inout] fileName - the file name.
 */
void MEDFileFieldGlobsReal::setFileName(const char *fileName)
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
void MEDFileFieldGlobsReal::changePflsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
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
void MEDFileFieldGlobsReal::changeLocsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
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
void MEDFileFieldGlobsReal::changePflsNames(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
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
void MEDFileFieldGlobsReal::changeLocsNames(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
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
void MEDFileFieldGlobsReal::changePflName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception)
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
void MEDFileFieldGlobsReal::changeLocName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception)
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
std::vector< std::pair<std::vector<std::string>, std::string > > MEDFileFieldGlobsReal::zipPflsNames() throw(INTERP_KERNEL::Exception)
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
std::vector< std::pair<std::vector<std::string>, std::string > > MEDFileFieldGlobsReal::zipLocsNames(double eps) throw(INTERP_KERNEL::Exception)
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
int MEDFileFieldGlobsReal::getNbOfGaussPtPerCell(int locId) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNull()->getNbOfGaussPtPerCell(locId);
}

/*!
 * Returns an id of a localization by its name.
 *  \param [in] loc - the localization name of interest.
 *  \return int - the id of the localization.
 *  \throw If there is no a localization named \a loc.
 */
int MEDFileFieldGlobsReal::getLocalizationId(const char *loc) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNull()->getLocalizationId(loc);
}

/*!
 * Returns the name of the MED file.
 *  \return const char * - the MED file name.
 */
const char *MEDFileFieldGlobsReal::getFileName() const
{
  return contentNotNull()->getFileName();
}

std::string MEDFileFieldGlobsReal::getFileName2() const
{
  return contentNotNull()->getFileName2();
}

/*!
 * Returns a localization object by its name.
 *  \param [in] locName - the name of the localization of interest.
 *  \return const MEDFileFieldLoc& - the localization object having the name \a locName.
 *  \throw If there is no a localization named \a locName.
 */
const MEDFileFieldLoc& MEDFileFieldGlobsReal::getLocalization(const char *locName) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNull()->getLocalization(locName);
}

/*!
 * Returns a localization object by its id.
 *  \param [in] locId - the id of the localization of interest.
 *  \return const MEDFileFieldLoc& - the localization object having the id \a locId.
 *  \throw If there is no a localization with id \a locId.
 */
const MEDFileFieldLoc& MEDFileFieldGlobsReal::getLocalizationFromId(int locId) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNull()->getLocalizationFromId(locId);
}

/*!
 * Returns a profile array by its name.
 *  \param [in] pflName - the name of the profile of interest.
 *  \return const DataArrayInt * - the profile array having the name \a pflName.
 *  \throw If there is no a profile named \a pflName.
 */
const DataArrayInt *MEDFileFieldGlobsReal::getProfile(const char *pflName) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNull()->getProfile(pflName);
}

/*!
 * Returns a profile array by its id.
 *  \param [in] pflId - the id of the profile of interest.
 *  \return const DataArrayInt * - the profile array having the id \a pflId.
 *  \throw If there is no a profile with id \a pflId.
 */
const DataArrayInt *MEDFileFieldGlobsReal::getProfileFromId(int pflId) const throw(INTERP_KERNEL::Exception)
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
MEDFileFieldLoc& MEDFileFieldGlobsReal::getLocalizationFromId(int locId) throw(INTERP_KERNEL::Exception)
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
MEDFileFieldLoc& MEDFileFieldGlobsReal::getLocalization(const char *locName) throw(INTERP_KERNEL::Exception)
{
  return contentNotNull()->getLocalization(locName);
}

/*!
 * Returns a profile array, apt for modification, by its name.
 *  \param [in] pflName - the name of the profile of interest.
 *  \return DataArrayInt * - a non-const pointer to the profile array having the name \a pflName.
 *  \throw If there is no a profile named \a pflName.
 */
DataArrayInt *MEDFileFieldGlobsReal::getProfile(const char *pflName) throw(INTERP_KERNEL::Exception)
{
  return contentNotNull()->getProfile(pflName);
}

/*!
 * Returns a profile array, apt for modification, by its id.
 *  \param [in] pflId - the id of the profile of interest.
 *  \return DataArrayInt * - a non-const pointer to the profile array having the id \a pflId.
 *  \throw If there is no a profile with id \a pflId.
 */
DataArrayInt *MEDFileFieldGlobsReal::getProfileFromId(int pflId) throw(INTERP_KERNEL::Exception)
{
  return contentNotNull()->getProfileFromId(pflId);
}

/*!
 * Removes profiles given by their ids. No data is updated to track this removal.
 *  \param [in] pflIds - a sequence of ids of the profiles to remove.
 */
void MEDFileFieldGlobsReal::killProfileIds(const std::vector<int>& pflIds) throw(INTERP_KERNEL::Exception)
{
  contentNotNull()->killProfileIds(pflIds);
}

/*!
 * Removes localizations given by their ids. No data is updated to track this removal.
 *  \param [in] locIds - a sequence of ids of the localizations to remove.
 */
void MEDFileFieldGlobsReal::killLocalizationIds(const std::vector<int>& locIds) throw(INTERP_KERNEL::Exception)
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
void MEDFileFieldGlobsReal::appendProfile(DataArrayInt *pfl) throw(INTERP_KERNEL::Exception)
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
void MEDFileFieldGlobsReal::appendLoc(const char *locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w) throw(INTERP_KERNEL::Exception)
{
  contentNotNull()->appendLoc(locName,geoType,refCoo,gsCoo,w);
}

MEDFileFieldGlobs *MEDFileFieldGlobsReal::contentNotNull() throw(INTERP_KERNEL::Exception)
{
  MEDFileFieldGlobs *g(_globals);
  if(!g)
    throw INTERP_KERNEL::Exception("MEDFileFieldGlobsReal::contentNotNull : no content in not const !");
  return g;
}

const MEDFileFieldGlobs *MEDFileFieldGlobsReal::contentNotNull() const throw(INTERP_KERNEL::Exception)
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

MEDFileFieldNameScope::MEDFileFieldNameScope(const char *fieldName):_name(fieldName)
{
}

/*!
 * Returns the name of \a this field.
 *  \return std::string - a string containing the field name.
 */
std::string MEDFileFieldNameScope::getName() const throw(INTERP_KERNEL::Exception)
{
  return _name;
}

/*!
 * Sets name of \a this field
 *  \param [in] name - the new field name.
 */
void MEDFileFieldNameScope::setName(const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  _name=fieldName;
}

std::string MEDFileFieldNameScope::getDtUnit() const throw(INTERP_KERNEL::Exception)
{
  return _dt_unit;
}

void MEDFileFieldNameScope::setDtUnit(const char *dtUnit) throw(INTERP_KERNEL::Exception)
{
  _dt_unit=dtUnit;
}

void MEDFileFieldNameScope::copyNameScope(const MEDFileFieldNameScope& other)
{
  _name=other._name;
  _dt_unit=other._dt_unit;
}

//= MEDFileAnyTypeField1TSWithoutSDA

/*!
 * Prints a string describing \a this field into a stream. This string is outputted 
 * by \c print Python command.
 *  \param [in] bkOffset - number of white spaces printed at the beginning of each line.
 *  \param [in,out] oss - the out stream.
 *  \param [in] f1tsId - the field index within a MED file. If \a f1tsId < 0, the tiny
 *          info id printed, else, not.
 */
void MEDFileAnyTypeField1TSWithoutSDA::simpleRepr(int bkOffset, std::ostream& oss, int f1tsId) const
{
  std::string startOfLine(bkOffset,' ');
  oss << startOfLine << "Field ";
  if(bkOffset==0)
    oss << "[Type=" << getTypeStr() << "] ";
  oss << "on One time Step ";
  if(f1tsId>=0)
    oss << "(" << f1tsId << ") ";
  oss << "on iteration=" << _iteration << " order=" << _order << "." << std::endl;
  oss << startOfLine << "Time attached is : " << _dt << " [" << _dt_unit << "]." << std::endl;
  const DataArray *arr=getUndergroundDataArray();
  if(arr)
    {
      const std::vector<std::string> &comps=arr->getInfoOnComponents();
      if(f1tsId<0)
        {
          oss << startOfLine << "Field Name : \"" << arr->getName() << "\"." << std::endl;
          oss << startOfLine << "Field has " << comps.size() << " components with the following infos :" << std::endl;
          for(std::vector<std::string>::const_iterator it=comps.begin();it!=comps.end();it++)
            oss << startOfLine << "  -  \"" << (*it) << "\"" << std::endl;
        }
      if(arr->isAllocated())
        {
          oss << startOfLine << "Whole field contains " << arr->getNumberOfTuples() << " tuples." << std::endl;
        }
      else
        oss << startOfLine << "The array of the current field has not allocated yet !" << std::endl;
    }
  else
    {
      oss << startOfLine << "Field infos are empty ! Not defined yet !" << std::endl;
    }
  oss << startOfLine << "----------------------" << std::endl;
  if(!_field_per_mesh.empty())
    {
      int i=0;
      for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it2=_field_per_mesh.begin();it2!=_field_per_mesh.end();it2++,i++)
        {
          const MEDFileFieldPerMesh *cur=(*it2);
          if(cur)
            cur->simpleRepr(bkOffset,oss,i);
          else
            oss << startOfLine << "Field per mesh #" << i << " is not defined !" << std::endl;
        }
    }
  else
    {
      oss << startOfLine << "Field is not defined on any meshes !" << std::endl;
    }
  oss << startOfLine << "----------------------" << std::endl;
}

std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> > MEDFileAnyTypeField1TSWithoutSDA::splitComponents() const throw(INTERP_KERNEL::Exception)
{
  const DataArray *arr(getUndergroundDataArray());
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::splitComponents : no array defined !");
  int nbOfCompo=arr->getNumberOfComponents();
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> > ret(nbOfCompo);
  for(int i=0;i<nbOfCompo;i++)
    {
      ret[i]=deepCpy();
      std::vector<int> v(1,i);
      MEDCouplingAutoRefCountObjectPtr<DataArray> arr2=arr->keepSelectedComponents(v);
      ret[i]->setArray(arr2);
    }
  return ret;
}

MEDFileAnyTypeField1TSWithoutSDA::MEDFileAnyTypeField1TSWithoutSDA(const char *fieldName, int csit, int iteration, int order):MEDFileFieldNameScope(fieldName),_iteration(iteration),_order(order),_csit(csit)
{
}

MEDFileAnyTypeField1TSWithoutSDA::MEDFileAnyTypeField1TSWithoutSDA():_iteration(-1),_order(-1),_dt(0.),_csit(-1)
{
}

/*!
 * Returns the maximal dimension of supporting elements. Returns -2 if \a this is
 * empty. Returns -1 if this in on nodes.
 *  \return int - the dimension of \a this.
 */
int MEDFileAnyTypeField1TSWithoutSDA::getDimension() const
{
  int ret=-2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    (*it)->getDimension(ret);
  return ret;
}

/*!
 * Returns the mesh name.
 *  \return std::string - a string holding the mesh name.
 *  \throw If \c _field_per_mesh.empty()
 */
std::string MEDFileAnyTypeField1TSWithoutSDA::getMeshName() const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::getMeshName : No field set !");
  return _field_per_mesh[0]->getMeshName();
}

void MEDFileAnyTypeField1TSWithoutSDA::setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception)
{
  std::string oldName(getMeshName());
  std::vector< std::pair<std::string,std::string> > v(1);
  v[0].first=oldName; v[0].second=newMeshName;
  changeMeshNames(v);
}

bool MEDFileAnyTypeField1TSWithoutSDA::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception)
{
  bool ret=false;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    {
      MEDFileFieldPerMesh *cur(*it);
      if(cur)
        ret=cur->changeMeshNames(modifTab) || ret;
    }
  return ret;
}

/*!
 * Returns the number of iteration of the state of underlying mesh.
 *  \return int - the iteration number.
 *  \throw If \c _field_per_mesh.empty()
 */
int MEDFileAnyTypeField1TSWithoutSDA::getMeshIteration() const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::getMeshIteration : No field set !");
  return _field_per_mesh[0]->getMeshIteration();
}

/*!
 * Returns the order number of iteration of the state of underlying mesh.
 *  \return int - the order number.
 *  \throw If \c _field_per_mesh.empty()
 */
int MEDFileAnyTypeField1TSWithoutSDA::getMeshOrder() const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::getMeshOrder : No field set !");
  return _field_per_mesh[0]->getMeshOrder();
}

/*!
 * Checks if \a this field is tagged by a given iteration number and a given
 * iteration order number.
 *  \param [in] iteration - the iteration number of interest.
 *  \param [in] order - the iteration order number of interest.
 *  \return bool - \c true if \a this->getIteration() == \a iteration && 
 *          \a this->getOrder() == \a order.
 */
bool MEDFileAnyTypeField1TSWithoutSDA::isDealingTS(int iteration, int order) const
{
  return iteration==_iteration && order==_order;
}

/*!
 * Returns number of iteration and order number of iteration when
 * \a this field has been calculated.
 *  \return std::pair<int,int> - a pair of the iteration number and the iteration
 *          order number.
 */
std::pair<int,int> MEDFileAnyTypeField1TSWithoutSDA::getDtIt() const
{
  std::pair<int,int> p;
  fillIteration(p);
  return p;
}

/*!
 * Returns number of iteration and order number of iteration when
 * \a this field has been calculated.
 *  \param [in,out] p - a pair returning the iteration number and the iteration
 *          order number.
 */
void MEDFileAnyTypeField1TSWithoutSDA::fillIteration(std::pair<int,int>& p) const
{
  p.first=_iteration;
  p.second=_order;
}

/*!
 * Returns all types of spatial discretization of \a this field.
 *  \param [in,out] types - a sequence of types of \a this field.
 */
void MEDFileAnyTypeField1TSWithoutSDA::fillTypesOfFieldAvailable(std::vector<TypeOfField>& types) const throw(INTERP_KERNEL::Exception)
{
  std::set<TypeOfField> types2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    {
      (*it)->fillTypesOfFieldAvailable(types2);
    }
  std::back_insert_iterator< std::vector<TypeOfField> > bi(types);
  std::copy(types2.begin(),types2.end(),bi);
}

/*!
 * Returns all types of spatial discretization of \a this field.
 *  \return std::vector<TypeOfField> - a sequence of types of spatial discretization
 *          of \a this field.
 */
std::vector<TypeOfField> MEDFileAnyTypeField1TSWithoutSDA::getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception)
{
  std::vector<TypeOfField> ret;
  fillTypesOfFieldAvailable(ret);
  return ret;
}

std::vector<std::string> MEDFileAnyTypeField1TSWithoutSDA::getPflsReallyUsed2() const
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

std::vector<std::string> MEDFileAnyTypeField1TSWithoutSDA::getLocsReallyUsed2() const
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

std::vector<std::string> MEDFileAnyTypeField1TSWithoutSDA::getPflsReallyUsedMulti2() const
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsedMulti();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

std::vector<std::string> MEDFileAnyTypeField1TSWithoutSDA::getLocsReallyUsedMulti2() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsedMulti();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

void MEDFileAnyTypeField1TSWithoutSDA::changePflsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    (*it)->changePflsRefsNamesGen(mapOfModif);
}

void MEDFileAnyTypeField1TSWithoutSDA::changeLocsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    (*it)->changeLocsRefsNamesGen(mapOfModif);
}

/*!
 * Returns all attributes of parts of \a this field lying on a given mesh.
 * Each part differs from other ones by a type of supporting mesh entity. The _i_-th
 * item of every of returned sequences refers to the _i_-th part of \a this field.
 * Thus all sequences returned by this method are of the same length equal to number
 * of different types of supporting entities.<br>
 * A field part can include sub-parts with several different spatial discretizations,
 * \ref ParaMEDMEM::ON_CELLS "ON_CELLS" and \ref ParaMEDMEM::ON_GAUSS_PT "ON_GAUSS_PT"
 * for example. Hence, some of the returned sequences contains nested sequences, and an item
 * of a nested sequence corresponds to a type of spatial discretization.<br>
 * This method allows for iteration over MEDFile DataStructure without any overhead.
 *  \param [in] mname - a name of a mesh of interest. It can be \c NULL, which is valid
 *          for the case with only one underlying mesh. (Actually, the number of meshes is
 *          not checked if \a mname == \c NULL).
 *  \param [in,out] types - a sequence of types of underlying mesh entities. A type per
 *          a field part is returned. 
 *  \param [in,out] typesF - a sequence of sequences of types of spatial discretizations.
 *          This sequence is of the same length as \a types. 
 *  \param [in,out] pfls - a sequence returning a profile name per each type of spatial
 *          discretization. A profile name can be empty.
 *          Length of this and of nested sequences is the same as that of \a typesF.
 *  \param [in,out] locs - a sequence returning a localization name per each type of spatial
 *          discretization. A localization name can be empty.
 *          Length of this and of nested sequences is the same as that of \a typesF.
 *  \return std::vector< std::vector< std::pair<int,int> > > - a sequence holding a range
 *          of ids of tuples within the data array, per each type of spatial
 *          discretization within one mesh entity type. 
 *          Length of this and of nested sequences is the same as that of \a typesF.
 *  \throw If no field is lying on \a mname.
 */
std::vector< std::vector< std::pair<int,int> > > MEDFileAnyTypeField1TSWithoutSDA::getFieldSplitedByType(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  int meshId=0;
  if(mname)
    meshId=getMeshIdFromMeshName(mname);
  else
    if(_field_per_mesh.empty())
      throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getFieldSplitedByType : This is empty !");
  return _field_per_mesh[meshId]->getFieldSplitedByType(types,typesF,pfls,locs);
}

/*!
 * Returns dimensions of mesh elements \a this field lies on. The returned value is a
 * maximal absolute dimension and values returned via the out parameter \a levs are 
 * dimensions relative to the maximal absolute dimension. <br>
 * This method is designed for MEDFileField1TS instances that have a discretization
 * \ref ParaMEDMEM::ON_CELLS "ON_CELLS", 
 * \ref ParaMEDMEM::ON_GAUSS_PT "ON_GAUSS_PT", 
 * \ref ParaMEDMEM::ON_GAUSS_NE "ON_GAUSS_NE".
 * Only these 3 discretizations will be taken into account here. If \a this is
 * \ref ParaMEDMEM::ON_NODES "ON_NODES", -1 is returned and \a levs are empty.<br>
 * This method is useful to make the link between the dimension of the underlying mesh
 * and the levels of \a this, because it is possible that the highest dimension of \a this
 * field is not equal to the dimension of the underlying mesh.
 * 
 * Let's consider the following case:
 * - mesh \a m1 has a meshDimension 3 and has non empty levels [0,-1,-2] with elements
 * TETRA4, HEXA8, TRI3 and SEG2.
 * - field \a f1 lies on \a m1 and is defined on 3D and 1D elements TETRA4 and SEG2.
 * - field \a f2 lies on \a m1 and is defined on 2D and 1D elements TRI3 and SEG2.
 *
 * In this case \a f1->getNonEmptyLevels() returns (3,[0,-2]) and \a
 * f2->getNonEmptyLevels() returns (2,[0,-1]). <br>
 * The returned values can be used for example to retrieve a MEDCouplingFieldDouble lying
 * on elements of a certain relative level by calling getFieldAtLevel(). \a meshDimRelToMax
 * parameter of getFieldAtLevel() is computed basing on the returned values as this:
 * <em> meshDimRelToMax = absDim - meshDim + relativeLev </em>.
 * For example<br>
 * to retrieve the highest level of
 * \a f1: <em>f1->getFieldAtLevel( ON_CELLS, 3-3+0 ); // absDim - meshDim + relativeLev</em><br> 
 * to retrieve the lowest level of \a f1: <em>f1->getFieldAtLevel( ON_CELLS, 3-3+(-2) );</em><br>
 * to retrieve the highest level of \a f2: <em>f2->getFieldAtLevel( ON_CELLS, 2-3+0 );</em><br>
 * to retrieve the lowest level of \a f2: <em>f2->getFieldAtLevel( ON_CELLS, 2-3+(-1) )</em>.
 *  \param [in] mname - a name of a mesh of interest. It can be \c NULL, which is valid
 *          for the case with only one underlying mesh. (Actually, the number of meshes is
 *          not checked if \a mname == \c NULL).
 *  \param [in,out] levs - a sequence returning the dimensions relative to the maximal
 *          absolute one. They are in decreasing order. This sequence is cleared before
 *          filling it in.
 *  \return int - the maximal absolute dimension of elements \a this fields lies on.
 *  \throw If no field is lying on \a mname.
 */
int MEDFileAnyTypeField1TSWithoutSDA::getNonEmptyLevels(const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception)
{
  levs.clear();
  int meshId=getMeshIdFromMeshName(mname);
  std::vector<INTERP_KERNEL::NormalizedCellType> types;
  std::vector< std::vector<TypeOfField> > typesF;
  std::vector< std::vector<std::string> > pfls, locs;
  _field_per_mesh[meshId]->getFieldSplitedByType(types,typesF,pfls,locs);
  if(types.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getNonEmptyLevels : 'this' is empty !");
  std::set<INTERP_KERNEL::NormalizedCellType> st(types.begin(),types.end());
  if(st.size()==1 && (*st.begin())==INTERP_KERNEL::NORM_ERROR)
    return -1;
  st.erase(INTERP_KERNEL::NORM_ERROR);
  std::set<int> ret1;
  for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator it=st.begin();it!=st.end();it++)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(*it);
      ret1.insert((int)cm.getDimension());
    }
  int ret=*std::max_element(ret1.begin(),ret1.end());
  std::copy(ret1.rbegin(),ret1.rend(),std::back_insert_iterator<std::vector<int> >(levs));
  std::transform(levs.begin(),levs.end(),levs.begin(),std::bind2nd(std::plus<int>(),-ret));
  return ret;
}

/*!
 * \param [in] mName specifies the underlying mesh name. This value can be pointer 0 for users that do not deal with fields on multi mesh.
 * \param [in] typ is for the geometric cell type (or INTERP_KERNEL::NORM_ERROR for node field) entry to find the right MEDFileFieldPerMeshPerTypePerDisc instance to set.
 * \param [in] locId is the localization id to find the right MEDFileFieldPerMeshPerTypePerDisc instance to set. It corresponds to the position of 
 *             \c pfls[std::distance(types.begin(),std::find(types.begin(),typ)] vector in MEDFileField1TSWithoutSDA::getFieldSplitedByType. For non gausspoints field users, the value is 0.
 */
MEDFileFieldPerMeshPerTypePerDisc *MEDFileAnyTypeField1TSWithoutSDA::getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) throw(INTERP_KERNEL::Exception)
{
  int mid=getMeshIdFromMeshName(mName);
  return _field_per_mesh[mid]->getLeafGivenTypeAndLocId(typ,locId);
}

/*!
 * \param [in] mName specifies the underlying mesh name. This value can be pointer 0 for users that do not deal with fields on multi mesh.
 * \param [in] typ is for the geometric cell type (or INTERP_KERNEL::NORM_ERROR for node field) entry to find the right MEDFileFieldPerMeshPerTypePerDisc instance to set.
 * \param [in] locId is the localization id to find the right MEDFileFieldPerMeshPerTypePerDisc instance to set. It corresponds to the position of 
 *             \c pfls[std::distance(types.begin(),std::find(types.begin(),typ)] vector in MEDFileField1TSWithoutSDA::getFieldSplitedByType. For non gausspoints field users, the value is 0.
 */
const MEDFileFieldPerMeshPerTypePerDisc *MEDFileAnyTypeField1TSWithoutSDA::getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) const throw(INTERP_KERNEL::Exception)
{
  int mid=getMeshIdFromMeshName(mName);
  return _field_per_mesh[mid]->getLeafGivenTypeAndLocId(typ,locId);
}

/*!
 * \param [in] mName specifies the underlying mesh name. This value can be pointer 0 for users that do not deal with fields on multi mesh.
 */
int MEDFileAnyTypeField1TSWithoutSDA::getMeshIdFromMeshName(const char *mName) const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getMeshIdFromMeshName : No field set !");
  if(mName==0)
    return 0;
  std::string mName2(mName);
  int ret=0;
  std::vector<std::string> msg;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++,ret++)
    if(mName2==(*it)->getMeshName())
      return ret;
    else
      msg.push_back((*it)->getMeshName());
  std::ostringstream oss; oss << "MEDFileField1TSWithoutSDA::getMeshIdFromMeshName : No such mesh \"" << mName2 << "\" as underlying mesh of field \"" << getName() << "\" !\n";
  oss << "Possible meshes are : ";
  for(std::vector<std::string>::const_iterator it2=msg.begin();it2!=msg.end();it2++)
    oss << "\"" << (*it2) << "\" ";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

int MEDFileAnyTypeField1TSWithoutSDA::addNewEntryIfNecessary(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::addNewEntryIfNecessary : input mesh is NULL !");
  std::string tmp(mesh->getName());
  if(tmp.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::addNewEntryIfNecessary : empty mesh name ! unsupported by MED file !");
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

bool MEDFileAnyTypeField1TSWithoutSDA::renumberEntitiesLyingOnMesh(const char *meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N,
                                                            MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  bool ret=false;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    {
      MEDFileFieldPerMesh *fpm(*it);
      if(fpm)
        ret=fpm->renumberEntitiesLyingOnMesh(meshName,oldCode,newCode,renumO2N,glob) || ret;
    }
  return ret;
}

void MEDFileAnyTypeField1TSWithoutSDA::writeLL(med_idt fid, const MEDFileWritable& opts, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::writeLL : empty field !");
  if(_field_per_mesh.size()>1)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::writeLL : In MED3.0 mode in writting mode only ONE underlying mesh supported !");
  _field_per_mesh[0]->copyOptionsFrom(opts);
  _field_per_mesh[0]->writeLL(fid,nasc);
}

void MEDFileAnyTypeField1TSWithoutSDA::finishLoading(med_idt fid, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  med_int numdt,numit;
  med_float dt;
  med_int nmesh;
  med_bool localMesh;
  med_int meshnumdt,meshnumit;
  INTERP_KERNEL::AutoPtr<char> meshName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  MEDfieldComputingStepInfo(fid,nasc.getName().c_str(),_csit,&numdt,&numit,&_dt);
  MEDfield23ComputingStepMeshInfo(fid,nasc.getName().c_str(),_csit,&numdt,&numit,&dt,&nmesh,meshName,&localMesh,&meshnumdt,&meshnumit);
  if(_iteration!=numdt || _order!=numit)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::finishLoading : unexpected exception internal error !");
  _field_per_mesh.resize(nmesh);
  for(int i=0;i<nmesh;i++)
    _field_per_mesh[i]=MEDFileFieldPerMesh::NewOnRead(fid,this,i+1,meshnumdt,meshnumit,nasc);//tony
  int start=0;
  for(int i=0;i<nmesh;i++)
    {
      _field_per_mesh[i]->prepareLoading(fid,start,nasc);
    }
  getOrCreateAndGetArray()->alloc(start,getNumberOfComponents());
  for(int i=0;i<nmesh;i++)
    {
      _field_per_mesh[i]->finishLoading(fid,nasc);
    }
}

std::size_t MEDFileAnyTypeField1TSWithoutSDA::getHeapMemorySize() const
{
  std::size_t ret=_dt_unit.capacity()+_field_per_mesh.capacity()*sizeof(MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh >);
  if(getUndergroundDataArray())
    ret+=getUndergroundDataArray()->getHeapMemorySize();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    ret+=(*it)->getHeapMemorySize();
  return ret;
}

/*!
 * Adds a MEDCouplingFieldDouble to \a this. The underlying mesh of the given field is
 * checked if its elements are sorted suitable for writing to MED file ("STB" stands for
 * "Sort By Type"), if not, an exception is thrown. 
 *  \param [in] field - the field to add to \a this. The array of field \a field is ignored
 *  \param [in] arr - the array of values.
 *  \param [in,out] glob - the global data where profiles and localization present in
 *          \a field, if any, are added.
 *  \throw If the name of \a field is empty.
 *  \throw If the data array of \a field is not set.
 *  \throw If \a this->_arr is already allocated but has different number of components
 *         than \a field.
 *  \throw If the underlying mesh of \a field has no name.
 *  \throw If elements in the mesh are not in the order suitable for writing to the MED file.
 */
void MEDFileAnyTypeField1TSWithoutSDA::setFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingMesh *mesh=field->getMesh();
  //
  TypeOfField type=field->getTypeOfField();
  std::vector<DataArrayInt *> dummy;
  int start=copyTinyInfoFrom(field,arr);
  int pos=addNewEntryIfNecessary(mesh);
  if(type!=ON_NODES)
    {
      std::vector<int> code=MEDFileField1TSWithoutSDA::CheckSBTMesh(mesh);
      _field_per_mesh[pos]->assignFieldNoProfileNoRenum(start,code,field,arr,glob,nasc);
    }
  else
    _field_per_mesh[pos]->assignNodeFieldNoProfile(start,field,arr,glob);
}

/*!
 * Adds a MEDCouplingFieldDouble to \a this. Specified entities of a given dimension
 * of a given mesh are used as the support of the given field (a real support is not used). 
 * Elements of the given mesh must be sorted suitable for writing to MED file. 
 * Order of underlying mesh entities of the given field specified by \a profile parameter
 * is not prescribed; this method permutes field values to have them sorted by element
 * type as required for writing to MED file. A new profile is added only if no equal
 * profile is missing. 
 *  \param [in] field - the field to add to \a this. The field double values are ignored.
 *  \param [in] arrOfVals - the values of the field \a field used.
 *  \param [in] mesh - the supporting mesh of \a field.
 *  \param [in] meshDimRelToMax - a relative dimension of mesh entities \a field lies on.
 *  \param [in] profile - ids of mesh entities on which corresponding field values lie.
 *  \param [in,out] glob - the global data where profiles and localization present in
 *          \a field, if any, are added.
 *  \throw If either \a field or \a mesh or \a profile has an empty name.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in \a mesh.
 *  \throw If the data array of \a field is not set.
 *  \throw If \a this->_arr is already allocated but has different number of components
 *         than \a field.
 *  \throw If elements in \a mesh are not in the order suitable for writing to the MED file.
 *  \sa setFieldNoProfileSBT()
 */
void MEDFileAnyTypeField1TSWithoutSDA::setFieldProfile(const MEDCouplingFieldDouble *field, const DataArray *arrOfVals, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception)
{
  TypeOfField type=field->getTypeOfField();
  int start=copyTinyInfoFrom(field,arrOfVals);
  std::vector<DataArrayInt *> idsInPflPerType;
  std::vector<DataArrayInt *> idsPerType;
  std::vector<int> code,code2;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=mesh->getGenMeshAtLevel(meshDimRelToMax);
  if(type!=ON_NODES)
    {
      m->splitProfilePerType(profile,code,idsInPflPerType,idsPerType);
      code2=m->getDistributionOfTypes();
      //
      std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > idsInPflPerType2(idsInPflPerType.size());
      for(std::size_t i=0;i<idsInPflPerType.size();i++)
        idsInPflPerType2[i]=idsInPflPerType[i];
      std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > idsPerType2(idsPerType.size());
      for(std::size_t i=0;i<idsPerType.size();i++)
        idsPerType2[i]=idsPerType[i];
      //
      int pos=addNewEntryIfNecessary(m);
      _field_per_mesh[pos]->assignFieldProfile(start,profile,code,code2,idsInPflPerType,idsPerType,field,arrOfVals,m,glob,nasc);
    }
  else
    {
      int pos=addNewEntryIfNecessary(m);
      _field_per_mesh[pos]->assignNodeFieldProfile(start,profile,field,arrOfVals,glob,nasc);
    }
}

/*!
 * Copies tiny info and allocates \a this->_arr instance of DataArrayDouble to
 * append data of a given MEDCouplingFieldDouble. So that the size of \a this->_arr becomes
 * larger by the size of \a field. Returns an id of the first not filled
 * tuple of \a this->_arr.
 *  \param [in] field - the field to copy the info on components and the name from.
 *  \return int - the id of first not initialized tuple of \a this->_arr.
 *  \throw If the name of \a field is empty.
 *  \throw If the data array of \a field is not set.
 *  \throw If \a this->_arr is already allocated but has different number of components
 *         than \a field.
 */
int MEDFileAnyTypeField1TSWithoutSDA::copyTinyInfoFrom(const MEDCouplingFieldDouble *field, const DataArray *arr) throw(INTERP_KERNEL::Exception)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::copyTinyInfoFrom : input field is NULL !");
  std::string name(field->getName());
  setName(name.c_str());
  setDtUnit(field->getTimeUnit());
  if(name.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::copyTinyInfoFrom : unsupported fields with no name in MED file !");
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::copyTinyInfoFrom : no array set !");
  _dt=field->getTime(_iteration,_order);
  int nbOfComponents=arr->getNumberOfComponents();
  getOrCreateAndGetArray()->setInfoAndChangeNbOfCompo(arr->getInfoOnComponents());
  if(!getOrCreateAndGetArray()->isAllocated())
    {
      getOrCreateAndGetArray()->alloc(arr->getNumberOfTuples(),arr->getNumberOfComponents());
      return 0;
    }
  else
    {
      int oldNbOfTuples=getOrCreateAndGetArray()->getNumberOfTuples();
      int newNbOfTuples=oldNbOfTuples+arr->getNumberOfTuples();
      MEDCouplingAutoRefCountObjectPtr<DataArray> tmp=createNewEmptyDataArrayInstance();
      tmp->alloc(newNbOfTuples,nbOfComponents);
      tmp->copyStringInfoFrom(*getOrCreateAndGetArray());
      DataArray *arrr=getOrCreateAndGetArray();      
      tmp->setContigPartOfSelectedValues2(0,arrr,0,oldNbOfTuples,1);
      setArray(tmp);
      return oldNbOfTuples;
    }
}

/*!
 * Returns number of components in \a this field
 *  \return int - the number of components.
 */
int MEDFileAnyTypeField1TSWithoutSDA::getNumberOfComponents() const
{
  return getOrCreateAndGetArray()->getNumberOfComponents();
}

/*!
 * Change info on components in \a this.
 * \throw If size of \a infos is not equal to the number of components already in \a this.
 */
void MEDFileAnyTypeField1TSWithoutSDA::setInfo(const std::vector<std::string>& infos) throw(INTERP_KERNEL::Exception)
{
  DataArray *arr=getOrCreateAndGetArray();
  arr->setInfoOnComponents(infos);//will throw an exception if number of components mimatches
}

/*!
 * Returns info on components of \a this field.
 *  \return const std::vector<std::string>& - a sequence of strings each being an
 *          information on _i_-th component.
 */
const std::vector<std::string>& MEDFileAnyTypeField1TSWithoutSDA::getInfo() const
{
  const DataArray *arr=getOrCreateAndGetArray();
  return arr->getInfoOnComponents();
}

/*!
 * Returns a mutable info on components of \a this field.
 *  \return std::vector<std::string>& - a sequence of strings each being an
 *          information on _i_-th component.
 */
std::vector<std::string>& MEDFileAnyTypeField1TSWithoutSDA::getInfo()
{
  DataArray *arr=getOrCreateAndGetArray();
  return arr->getInfoOnComponents();
}

/*!
 * Returns a new MEDCouplingFieldDouble of given type lying on a given support.
 *  \param [in] type - a spatial discretization of the new field.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [in] mName - a name of the supporting mesh.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \param [in] glob - the global data storing profiles and localization.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If the MED file is not readable.
 *  \throw If there is no mesh named \a mName in the MED file.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in the mesh.
 *  \throw If no field of \a this is lying on the mesh \a mName.
 *  \throw If no field values of the given \a type or given \a meshDimRelToMax are available.
 */
MEDCouplingFieldDouble *MEDFileAnyTypeField1TSWithoutSDA::getFieldAtLevel(TypeOfField type, int meshDimRelToMax, const char *mName, int renumPol, const MEDFileFieldGlobsReal *glob, MEDCouplingAutoRefCountObjectPtr<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm;
  if(mName==0)
    mm=MEDFileMesh::New(glob->getFileName(),getMeshName().c_str(),getMeshIteration(),getMeshOrder());
  else
    mm=MEDFileMesh::New(glob->getFileName(),mName,getMeshIteration(),getMeshOrder());
  return MEDFileAnyTypeField1TSWithoutSDA::getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,glob,mm,arrOut,nasc);
}

/*!
 * Returns a new MEDCouplingFieldDouble of given type lying on a given support.
 *  \param [in] type - a spatial discretization of the new field.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \param [in] glob - the global data storing profiles and localization.
 *  \param [in] mesh - the supporting mesh.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If the MED file is not readable.
 *  \throw If no field of \a this is lying on \a mesh.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in the mesh.
 *  \throw If no field values of the given \a type or given \a meshDimRelToMax are available.
 */
MEDCouplingFieldDouble *MEDFileAnyTypeField1TSWithoutSDA::getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol, const MEDFileFieldGlobsReal *glob, const MEDFileMesh *mesh, MEDCouplingAutoRefCountObjectPtr<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=mesh->getGenMeshAtLevel(meshDimRelToMax,false);
  const DataArrayInt *d=mesh->getNumberFieldAtLevel(meshDimRelToMax);
  const DataArrayInt *e=mesh->getNumberFieldAtLevel(1);
  if(meshDimRelToMax==1)
    (static_cast<MEDCouplingUMesh *>((MEDCouplingMesh *)m))->setMeshDimension(0);
  return MEDFileAnyTypeField1TSWithoutSDA::getFieldOnMeshAtLevel(type,renumPol,glob,m,d,e,arrOut,nasc);
}

/*!
 * Returns a new MEDCouplingFieldDouble of a given type lying on the top level cells of a
 * given mesh. 
 *  \param [in] type - a spatial discretization of the new field.
 *  \param [in] mName - a name of the supporting mesh.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \param [in] glob - the global data storing profiles and localization.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If the MED file is not readable.
 *  \throw If there is no mesh named \a mName in the MED file.
 *  \throw If there are no mesh entities in the mesh.
 *  \throw If no field values of the given \a type are available.
 */
MEDCouplingFieldDouble *MEDFileAnyTypeField1TSWithoutSDA::getFieldAtTopLevel(TypeOfField type, const char *mName, int renumPol, const MEDFileFieldGlobsReal *glob, MEDCouplingAutoRefCountObjectPtr<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception)
{
   MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm;
  if(mName==0)
    mm=MEDFileMesh::New(glob->getFileName(),getMeshName().c_str(),getMeshIteration(),getMeshOrder());
  else
    mm=MEDFileMesh::New(glob->getFileName(),mName,getMeshIteration(),getMeshOrder());
  int absDim=getDimension();
  int meshDimRelToMax=absDim-mm->getMeshDimension();
  return MEDFileAnyTypeField1TSWithoutSDA::getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,glob,mm,arrOut,nasc);
}

/*!
 * Returns a new MEDCouplingFieldDouble of given type lying on a given support.
 *  \param [in] type - a spatial discretization of the new field.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \param [in] glob - the global data storing profiles and localization.
 *  \param [in] mesh - the supporting mesh.
 *  \param [in] cellRenum - the cell numbers array used for permutation of the result
 *         field according to \a renumPol.
 *  \param [in] nodeRenum - the node numbers array used for permutation of the result
 *         field according to \a renumPol.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in the mesh.
 *  \throw If no field of \a this is lying on \a mesh.
 *  \throw If no field values of the given \a type or given \a meshDimRelToMax are available.
 */
MEDCouplingFieldDouble *MEDFileAnyTypeField1TSWithoutSDA::getFieldOnMeshAtLevel(TypeOfField type, int renumPol, const MEDFileFieldGlobsReal *glob, const MEDCouplingMesh *mesh, const DataArrayInt *cellRenum, const DataArrayInt *nodeRenum, MEDCouplingAutoRefCountObjectPtr<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception)
{
  static const char msg1[]="MEDFileField1TSWithoutSDA::getFieldOnMeshAtLevel : request for a renumbered field following mesh numbering whereas it is a profile field !";
  int meshId=getMeshIdFromMeshName(mesh->getName());
  bool isPfl=false;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=_field_per_mesh[meshId]->getFieldOnMeshAtLevel(type,glob,mesh,isPfl,arrOut,nasc);
  switch(renumPol)
    {
    case 0:
      {
        //no need to test _field_per_mesh.empty() because geMeshName has already done it
        return ret.retn();
      }
    case 3:
    case 1:
      {
        if(isPfl)
          throw INTERP_KERNEL::Exception(msg1);
        //no need to test _field_per_mesh.empty() because geMeshName has already done it
        if(cellRenum)
          {
            if((int)cellRenum->getNbOfElems()!=mesh->getNumberOfCells())
              {
                std::ostringstream oss; oss << "MEDFileField1TSWithoutSDA::getFieldOnMeshAtLevel : Request of simple renumbering but it seems that underlying mesh \"" << mesh->getName() << "\" of requested field ";
                oss << "\"" << getName() << "\" has partial renumbering (some geotype has no renumber) !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
            MEDCouplingFieldDiscretization *disc=ret->getDiscretization();
            if(!disc) throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::getFieldOnMeshAtLevel : internal error, no discretization on field !");
            std::vector<DataArray *> arrOut2(1,arrOut);
            // 2 following lines replace ret->renumberCells(cellRenum->getConstPointer()) if not DataArrayDouble
            disc->renumberArraysForCell(ret->getMesh(),arrOut2,cellRenum->getConstPointer(),true);
            (const_cast<MEDCouplingMesh*>(ret->getMesh()))->renumberCells(cellRenum->getConstPointer(),true);
          }
        if(renumPol==1)
          return ret.retn();
      }
    case 2:
      {
        //no need to test _field_per_mesh.empty() because geMeshName has already done it
        if(isPfl)
          throw INTERP_KERNEL::Exception(msg1);
        if(nodeRenum)
          {
            if((int)nodeRenum->getNbOfElems()!=mesh->getNumberOfNodes())
              {
                std::ostringstream oss; oss << "MEDFileField1TSWithoutSDA::getFieldOnMeshAtLevel : Request of simple renumbering but it seems that underlying mesh \"" << mesh->getName() << "\" of requested field ";
                oss << "\"" << nasc.getName() << "\" not defined on all nodes !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
            MEDCouplingAutoRefCountObjectPtr<DataArrayInt> nodeRenumSafe=nodeRenum->checkAndPreparePermutation();
            if(!dynamic_cast<DataArrayDouble *>((DataArray *)arrOut))
              throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getFieldOnMeshAtLevel : node renumbering not implemented for not double DataArrays !");
            ret->renumberNodes(nodeRenumSafe->getConstPointer());
          }
        return ret.retn();
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getFieldOnMeshAtLevel : unsupported renum policy ! Dealing with policy 0 1 2 and 3 !");
    }
}

/*!
 * Returns values and a profile of the field of a given type lying on a given support.
 *  \param [in] type - a spatial discretization of the field.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [in] mesh - the supporting mesh.
 *  \param [out] pfl - a new instance of DataArrayInt holding ids of mesh entities the
 *          field of interest lies on. If the field lies on all entities of the given
 *          dimension, all ids in \a pfl are zero. The caller is to delete this array
 *          using decrRef() as it is no more needed.  
 *  \param [in] glob - the global data storing profiles and localization.
 *  \return DataArrayDouble * - a new instance of DataArrayDouble holding values of the
 *          field. The caller is to delete this array using decrRef() as it is no more needed.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in \a mesh.
 *  \throw If no field of \a this is lying on \a mesh.
 *  \throw If no field values of the given \a type are available.
 */
DataArray *MEDFileAnyTypeField1TSWithoutSDA::getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl, const MEDFileFieldGlobsReal *glob, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=mesh->getGenMeshAtLevel(meshDimRelToMax);
  int meshId=getMeshIdFromMeshName(mesh->getName());
  MEDCouplingAutoRefCountObjectPtr<DataArray> ret=_field_per_mesh[meshId]->getFieldOnMeshAtLevelWithPfl(type,m,pfl,glob,nasc);
  ret->setName(nasc.getName().c_str());
  return ret.retn();
}

//= MEDFileField1TSWithoutSDA

/*!
 * Throws if a given value is not a valid (non-extended) relative dimension.
 *  \param [in] meshDimRelToMax - the relative dimension value.
 *  \throw If \a meshDimRelToMax > 0.
 */
void MEDFileField1TSWithoutSDA::CheckMeshDimRel(int meshDimRelToMax) throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMax>0)
    throw INTERP_KERNEL::Exception("CheckMeshDimRel : This is a meshDimRel not a meshDimRelExt ! So value should be <=0 !");
}

/*!
 * Checks if elements of a given mesh are in the order suitable for writing 
 * to the MED file. If this is not so, an exception is thrown. In a case of success, returns a
 * vector describing types of elements and their number.
 *  \param [in] mesh - the mesh to check.
 *  \return std::vector<int> - a vector holding for each element type (1) item of
 *          INTERP_KERNEL::NormalizedCellType, (2) number of elements, (3) -1. 
 *          These values are in full-interlace mode.
 *  \throw If elements in \a mesh are not in the order suitable for writing to the MED file.
 */
std::vector<int> MEDFileField1TSWithoutSDA::CheckSBTMesh(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::CheckSBTMesh : input mesh is NULL !");
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
      throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::CheckSBTMesh : underlying mesh is not sorted by type as MED file expects !");
    }
  return code;
}

MEDFileField1TSWithoutSDA *MEDFileField1TSWithoutSDA::New(const char *fieldName, int csit, int iteration, int order, const std::vector<std::string>& infos)
{
  return new MEDFileField1TSWithoutSDA(fieldName,csit,iteration,order,infos);
}

/*!
 * Returns all attributes and values of parts of \a this field lying on a given mesh.
 * Each part differs from other ones by a type of supporting mesh entity. The _i_-th
 * item of every of returned sequences refers to the _i_-th part of \a this field.
 * Thus all sequences returned by this method are of the same length equal to number
 * of different types of supporting entities.<br>
 * A field part can include sub-parts with several different spatial discretizations,
 * \ref ParaMEDMEM::ON_CELLS "ON_CELLS" and \ref ParaMEDMEM::ON_GAUSS_PT "ON_GAUSS_PT"
 * for example. Hence, some of the returned sequences contains nested sequences, and an item
 * of a nested sequence corresponds to a type of spatial discretization.<br>
 * This method allows for iteration over MEDFile DataStructure with a reduced overhead.
 * The overhead is due to selecting values into new instances of DataArrayDouble.
 *  \param [in] mname - a name of a mesh of interest. It can be \c NULL, which is valid
 *          for the case with only one underlying mesh. (Actually, the number of meshes is
 *          not checked if \a mname == \c NULL).
 *  \param [in,out] types - a sequence of types of underlying mesh entities. A type per
 *          a field part is returned. 
 *  \param [in,out] typesF - a sequence of sequences of types of spatial discretizations.
 *          A field part can include sub-parts with several different spatial discretizations,
 *          \ref ParaMEDMEM::ON_CELLS "ON_CELLS" and 
 *          \ref ParaMEDMEM::ON_GAUSS_PT "ON_GAUSS_PT" for example.
 *          This sequence is of the same length as \a types. 
 *  \param [in,out] pfls - a sequence returning a profile name per each type of spatial
 *          discretization. A profile name can be empty.
 *          Length of this and of nested sequences is the same as that of \a typesF.
 *  \param [in,out] locs - a sequence returning a localization name per each type of spatial
 *          discretization. A localization name can be empty.
 *          Length of this and of nested sequences is the same as that of \a typesF.
 *  \return std::vector< std::vector<DataArrayDouble *> > - a sequence holding arrays of values
 *          per each type of spatial discretization within one mesh entity type.
 *          The caller is to delete each DataArrayDouble using decrRef() as it is no more needed.
 *          Length of this and of nested sequences is the same as that of \a typesF.
 *  \throw If no field is lying on \a mname.
 */
std::vector< std::vector<DataArrayDouble *> > MEDFileField1TSWithoutSDA::getFieldSplitedByType2(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  int meshId=0;
  if(mname)
    meshId=getMeshIdFromMeshName(mname);
  else
    if(_field_per_mesh.empty())
      throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getFieldSplitedByType : This is empty !");
  std::vector< std::vector< std::pair<int,int> > > ret0=_field_per_mesh[meshId]->getFieldSplitedByType(types,typesF,pfls,locs);
  int nbOfRet=ret0.size();
  std::vector< std::vector<DataArrayDouble *> > ret(nbOfRet);
  for(int i=0;i<nbOfRet;i++)
    {
      const std::vector< std::pair<int,int> >& p=ret0[i];
      int nbOfRet1=p.size();
      ret[i].resize(nbOfRet1);
      for(int j=0;j<nbOfRet1;j++)
        {
          DataArrayDouble *tmp=_arr->selectByTupleId2(p[j].first,p[j].second,1);
          ret[i][j]=tmp;
        }
    }
  return ret;
}

/*!
 * Returns a pointer to the underground DataArrayDouble instance. So the
 * caller should not decrRef() it. This method allows for a direct access to the field
 * values. This method is quite unusable if there is more than a nodal field or a cell
 * field on single geometric cell type. 
 *  \return DataArrayDouble * - the pointer to the field values array.
 */
DataArrayDouble *MEDFileField1TSWithoutSDA::getUndergroundDataArrayDouble() const throw(INTERP_KERNEL::Exception)
{
  const DataArrayDouble *ret=_arr;
  if(ret)
    return const_cast<DataArrayDouble *>(ret);
  else
    return 0;
}

const char *MEDFileField1TSWithoutSDA::getTypeStr() const throw(INTERP_KERNEL::Exception)
{
  return TYPE_STR;
}

/*!
 * Returns a pointer to the underground DataArrayDouble instance. So the
 * caller should not decrRef() it. This method allows for a direct access to the field
 * values. This method is quite unusable if there is more than a nodal field or a cell
 * field on single geometric cell type. 
 *  \return DataArrayDouble * - the pointer to the field values array.
 */
DataArray *MEDFileField1TSWithoutSDA::getUndergroundDataArray() const throw(INTERP_KERNEL::Exception)
{
  return getUndergroundDataArrayDouble();
}

/*!
 * Returns a pointer to the underground DataArrayDouble instance and a
 * sequence describing parameters of a support of each part of \a this field. The
 * caller should not decrRef() the returned DataArrayDouble. This method allows for a
 * direct access to the field values. This method is intended for the field lying on one
 * mesh only.
 *  \param [in,out] entries - the sequence describing parameters of a support of each
 *         part of \a this field. Each item of this sequence consists of two parts. The
 *         first part describes a type of mesh entity and an id of discretization of a
 *         current field part. The second part describes a range of values [begin,end)
 *         within the returned array relating to the current field part.
 *  \return DataArrayDouble * - the pointer to the field values array.
 *  \throw If the number of underlying meshes is not equal to 1.
 *  \throw If no field values are available.
 *  \sa getUndergroundDataArray()
 */
DataArrayDouble *MEDFileField1TSWithoutSDA::getUndergroundDataArrayDoubleExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.size()!=1)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt : field lies on several meshes, this method has no sense !");
  if(_field_per_mesh[0]==0)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt : no field specified !");
  _field_per_mesh[0]->getUndergroundDataArrayExt(entries);
  return getUndergroundDataArrayDouble();
}

/*!
 * Returns a pointer to the underground DataArrayDouble instance and a
 * sequence describing parameters of a support of each part of \a this field. The
 * caller should not decrRef() the returned DataArrayDouble. This method allows for a
 * direct access to the field values. This method is intended for the field lying on one
 * mesh only.
 *  \param [in,out] entries - the sequence describing parameters of a support of each
 *         part of \a this field. Each item of this sequence consists of two parts. The
 *         first part describes a type of mesh entity and an id of discretization of a
 *         current field part. The second part describes a range of values [begin,end)
 *         within the returned array relating to the current field part.
 *  \return DataArrayDouble * - the pointer to the field values array.
 *  \throw If the number of underlying meshes is not equal to 1.
 *  \throw If no field values are available.
 *  \sa getUndergroundDataArray()
 */
DataArray *MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception)
{
  return getUndergroundDataArrayDoubleExt(entries);
}

MEDFileField1TSWithoutSDA::MEDFileField1TSWithoutSDA(const char *fieldName, int csit, int iteration, int order,
                                                     const std::vector<std::string>& infos):MEDFileAnyTypeField1TSWithoutSDA(fieldName,csit,iteration,order)
{
  DataArrayDouble *arr=getOrCreateAndGetArrayDouble();
  arr->setInfoAndChangeNbOfCompo(infos);
}

MEDFileField1TSWithoutSDA::MEDFileField1TSWithoutSDA():MEDFileAnyTypeField1TSWithoutSDA()
{
}

MEDFileAnyTypeField1TSWithoutSDA *MEDFileField1TSWithoutSDA::shallowCpy() const throw(INTERP_KERNEL::Exception)
{
  return new MEDFileField1TSWithoutSDA(*this);
}

MEDFileAnyTypeField1TSWithoutSDA *MEDFileField1TSWithoutSDA::deepCpy() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA> ret=static_cast<MEDFileField1TSWithoutSDA *>(shallowCpy());
  if((const DataArrayDouble *)_arr)
    ret->_arr=_arr->deepCpy();
  std::size_t i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++,i++)
    {
      if((const MEDFileFieldPerMesh *)*it)
        ret->_field_per_mesh[i]=(*it)->deepCpy((MEDFileField1TSWithoutSDA *)ret);
    }
  return ret.retn();
}

void MEDFileField1TSWithoutSDA::setArray(DataArray *arr) throw(INTERP_KERNEL::Exception)
{
  if(!arr)
    _arr=0;
  DataArrayDouble *arrC=dynamic_cast<DataArrayDouble *>(arr);
  if(!arrC)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::setArray : the input not null array is not of type DataArrayDouble !");
  arrC->incrRef();
  _arr=arrC;
}

DataArray *MEDFileField1TSWithoutSDA::createNewEmptyDataArrayInstance() const
{
  return DataArrayDouble::New();
}

DataArrayDouble *MEDFileField1TSWithoutSDA::getOrCreateAndGetArrayDouble()
{
  DataArrayDouble *ret=_arr;
  if(ret)
    return ret;
  _arr=DataArrayDouble::New();
  return _arr;
}

DataArray *MEDFileField1TSWithoutSDA::getOrCreateAndGetArray()
{
  return getOrCreateAndGetArrayDouble();
}

const DataArrayDouble *MEDFileField1TSWithoutSDA::getOrCreateAndGetArrayDouble() const
{
  const DataArrayDouble *ret=_arr;
  if(ret)
    return ret;
  DataArrayDouble *ret2=DataArrayDouble::New();
  const_cast<MEDFileField1TSWithoutSDA *>(this)->_arr=DataArrayDouble::New();
  return ret2;
}

const DataArray *MEDFileField1TSWithoutSDA::getOrCreateAndGetArray() const
{
  return getOrCreateAndGetArrayDouble();
}

//= MEDFileIntField1TSWithoutSDA

MEDFileIntField1TSWithoutSDA *MEDFileIntField1TSWithoutSDA::New(const char *fieldName, int csit, int iteration, int order,
                                                                const std::vector<std::string>& infos)
{
  return new MEDFileIntField1TSWithoutSDA(fieldName,csit,iteration,order,infos);
}

MEDFileIntField1TSWithoutSDA::MEDFileIntField1TSWithoutSDA():MEDFileAnyTypeField1TSWithoutSDA()
{
}

MEDFileIntField1TSWithoutSDA::MEDFileIntField1TSWithoutSDA(const char *fieldName, int csit, int iteration, int order,
                                                           const std::vector<std::string>& infos):MEDFileAnyTypeField1TSWithoutSDA(fieldName,csit,iteration,order)
{
  DataArrayInt *arr=getOrCreateAndGetArrayInt();
  arr->setInfoAndChangeNbOfCompo(infos);
}

const char *MEDFileIntField1TSWithoutSDA::getTypeStr() const throw(INTERP_KERNEL::Exception)
{
  return TYPE_STR;
}

/*!
 * Returns a pointer to the underground DataArrayInt instance. So the
 * caller should not decrRef() it. This method allows for a direct access to the field
 * values. This method is quite unusable if there is more than a nodal field or a cell
 * field on single geometric cell type. 
 *  \return DataArrayInt * - the pointer to the field values array.
 */
DataArray *MEDFileIntField1TSWithoutSDA::getUndergroundDataArray() const throw(INTERP_KERNEL::Exception)
{
  return getUndergroundDataArrayInt();
}

/*!
 * Returns a pointer to the underground DataArrayInt instance. So the
 * caller should not decrRef() it. This method allows for a direct access to the field
 * values. This method is quite unusable if there is more than a nodal field or a cell
 * field on single geometric cell type. 
 *  \return DataArrayInt * - the pointer to the field values array.
 */
DataArrayInt *MEDFileIntField1TSWithoutSDA::getUndergroundDataArrayInt() const throw(INTERP_KERNEL::Exception)
{
  const DataArrayInt *ret=_arr;
  if(ret)
    return const_cast<DataArrayInt *>(ret);
  else
    return 0;
}

/*!
 * Returns a pointer to the underground DataArrayInt instance and a
 * sequence describing parameters of a support of each part of \a this field. The
 * caller should not decrRef() the returned DataArrayInt. This method allows for a
 * direct access to the field values. This method is intended for the field lying on one
 * mesh only.
 *  \param [in,out] entries - the sequence describing parameters of a support of each
 *         part of \a this field. Each item of this sequence consists of two parts. The
 *         first part describes a type of mesh entity and an id of discretization of a
 *         current field part. The second part describes a range of values [begin,end)
 *         within the returned array relating to the current field part.
 *  \return DataArrayInt * - the pointer to the field values array.
 *  \throw If the number of underlying meshes is not equal to 1.
 *  \throw If no field values are available.
 *  \sa getUndergroundDataArray()
 */
DataArray *MEDFileIntField1TSWithoutSDA::getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception)
{
  return getUndergroundDataArrayIntExt(entries);
}

/*!
 * Returns a pointer to the underground DataArrayInt instance and a
 * sequence describing parameters of a support of each part of \a this field. The
 * caller should not decrRef() the returned DataArrayInt. This method allows for a
 * direct access to the field values. This method is intended for the field lying on one
 * mesh only.
 *  \param [in,out] entries - the sequence describing parameters of a support of each
 *         part of \a this field. Each item of this sequence consists of two parts. The
 *         first part describes a type of mesh entity and an id of discretization of a
 *         current field part. The second part describes a range of values [begin,end)
 *         within the returned array relating to the current field part.
 *  \return DataArrayInt * - the pointer to the field values array.
 *  \throw If the number of underlying meshes is not equal to 1.
 *  \throw If no field values are available.
 *  \sa getUndergroundDataArray()
 */
DataArrayInt *MEDFileIntField1TSWithoutSDA::getUndergroundDataArrayIntExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.size()!=1)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt : field lies on several meshes, this method has no sense !");
  if(_field_per_mesh[0]==0)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt : no field specified !");
  _field_per_mesh[0]->getUndergroundDataArrayExt(entries);
  return getUndergroundDataArrayInt();
}

MEDFileAnyTypeField1TSWithoutSDA *MEDFileIntField1TSWithoutSDA::shallowCpy() const throw(INTERP_KERNEL::Exception)
{
  return new MEDFileIntField1TSWithoutSDA(*this);
}

MEDFileAnyTypeField1TSWithoutSDA *MEDFileIntField1TSWithoutSDA::deepCpy() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileIntField1TSWithoutSDA> ret=static_cast<MEDFileIntField1TSWithoutSDA *>(shallowCpy());
  if((const DataArrayInt *)_arr)
    ret->_arr=_arr->deepCpy();
  std::size_t i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++,i++)
    {
      if((const MEDFileFieldPerMesh *)*it)
        ret->_field_per_mesh[i]=(*it)->deepCpy((MEDFileIntField1TSWithoutSDA *)ret);
    }
  return ret.retn();
}

void MEDFileIntField1TSWithoutSDA::setArray(DataArray *arr) throw(INTERP_KERNEL::Exception)
{
  if(!arr)
    _arr=0;
  DataArrayInt *arrC=dynamic_cast<DataArrayInt *>(arr);
  if(!arrC)
    throw INTERP_KERNEL::Exception("MEDFileIntField1TSWithoutSDA::setArray : the input not null array is not of type DataArrayInt !");
  _arr=arrC;
}

DataArray *MEDFileIntField1TSWithoutSDA::createNewEmptyDataArrayInstance() const
{
  return DataArrayInt::New();
}

DataArrayInt *MEDFileIntField1TSWithoutSDA::getOrCreateAndGetArrayInt()
{
  DataArrayInt *ret=_arr;
  if(ret)
    return ret;
  _arr=DataArrayInt::New();
  return _arr;
}

DataArray *MEDFileIntField1TSWithoutSDA::getOrCreateAndGetArray()
{
  return getOrCreateAndGetArrayInt();
}

const DataArrayInt *MEDFileIntField1TSWithoutSDA::getOrCreateAndGetArrayInt() const
{
  const DataArrayInt *ret=_arr;
  if(ret)
    return ret;
  DataArrayInt *ret2=DataArrayInt::New();
  const_cast<MEDFileIntField1TSWithoutSDA *>(this)->_arr=DataArrayInt::New();
  return ret2;
}

const DataArray *MEDFileIntField1TSWithoutSDA::getOrCreateAndGetArray() const
{
  return getOrCreateAndGetArrayInt();
}

MEDFileAnyTypeField1TS::MEDFileAnyTypeField1TS()
{
}

//= MEDFileAnyTypeField1TS

MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeField1TS::BuildContentFrom(med_idt fid, const char *fileName) throw(INTERP_KERNEL::Exception)
{
  med_field_type typcha;
  //
  std::vector<std::string> infos;
  std::string dtunit,fieldName;
  LocateField2(fid,fileName,0,true,fieldName,typcha,infos,dtunit);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> ret;
  switch(typcha)
    {
    case MED_FLOAT64:
      {
        ret=MEDFileField1TSWithoutSDA::New(fieldName.c_str(),-1,-1/*iteration*/,-1/*order*/,std::vector<std::string>());
        break;
      }
    case MED_INT32:
      {
        ret=MEDFileIntField1TSWithoutSDA::New(fieldName.c_str(),-1,-1/*iteration*/,-1/*order*/,std::vector<std::string>());
        break;
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileAnyTypeField1TS::BuildContentFrom(fileName) : file \'" << fileName << "\' contains field with name \'" << fieldName << "\' but the type of the first field is not in [MED_FLOAT64, MED_INT32] !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    }
  ret->setDtUnit(dtunit.c_str());
  ret->getOrCreateAndGetArray()->setInfoAndChangeNbOfCompo(infos);
  //
  med_int numdt,numit;
  med_float dt;
  MEDfieldComputingStepInfo(fid,fieldName.c_str(),1,&numdt,&numit,&dt);
  ret->setTime(numdt,numit,dt);
  ret->_csit=1;
  ret->finishLoading(fid,*((const MEDFileAnyTypeField1TSWithoutSDA*)ret));
  return ret.retn();
}

MEDFileAnyTypeField1TS::MEDFileAnyTypeField1TS(const char *fileName) throw(INTERP_KERNEL::Exception)
try:MEDFileFieldGlobsReal(fileName)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  _content=BuildContentFrom(fid,fileName);
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeField1TS::BuildContentFrom(med_idt fid, const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  med_field_type typcha;
  std::vector<std::string> infos;
  std::string dtunit;
  int iii=-1;
  int nbSteps=LocateField(fid,fileName,fieldName,iii,typcha,infos,dtunit);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> ret;
  switch(typcha)
    {
    case MED_FLOAT64:
      {
        ret=MEDFileField1TSWithoutSDA::New(fieldName,-1,-1/*iteration*/,-1/*order*/,std::vector<std::string>());
        break;
      }
    case MED_INT32:
      {
        ret=MEDFileIntField1TSWithoutSDA::New(fieldName,-1,-1/*iteration*/,-1/*order*/,std::vector<std::string>());
        break;
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileAnyTypeField1TS::BuildContentFrom(fileName,fieldName) : file \'" << fileName << "\' contains field with name \'" << fieldName << "\' but the type of field is not in [MED_FLOAT64, MED_INT32] !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    }
  ret->setDtUnit(dtunit.c_str());
  ret->getOrCreateAndGetArray()->setInfoAndChangeNbOfCompo(infos);
  //
  if(nbSteps<1)
    {
      std::ostringstream oss; oss << "MEDFileField1TS(fileName,fieldName) : file \'" << fileName << "\' contains field with name \'" << fieldName << "\' but there is no time steps on it !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  //
  med_int numdt,numit;
  med_float dt;
  MEDfieldComputingStepInfo(fid,fieldName,1,&numdt,&numit,&dt);
  ret->setTime(numdt,numit,dt);
  ret->_csit=1;
  ret->finishLoading(fid,*((const MEDFileAnyTypeField1TSWithoutSDA*)ret));
  return ret.retn();
}

MEDFileAnyTypeField1TS::MEDFileAnyTypeField1TS(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
try:MEDFileFieldGlobsReal(fileName)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  _content=BuildContentFrom(fid,fileName,fieldName);
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::BuildNewInstanceFromContent(MEDFileAnyTypeField1TSWithoutSDA *c, const char *fileName) throw(INTERP_KERNEL::Exception)
{
  if(!c)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TS::BuildNewInstanceFromContent : empty content in input : unable to build a new instance !");
  if(dynamic_cast<const MEDFileField1TSWithoutSDA *>(c))
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> ret=MEDFileField1TS::New();
      ret->setFileName(fileName);
      ret->_content=c; c->incrRef();
      return ret.retn();
    }
  if(dynamic_cast<const MEDFileIntField1TSWithoutSDA *>(c))
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileIntField1TS> ret=MEDFileIntField1TS::New();
      ret->setFileName(fileName);
      ret->_content=c; c->incrRef();
      return ret.retn();
    }
  throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TS::BuildNewInstanceFromContent : internal error ! a content of type different from FLOAT64 and INT32 has been built but not intercepted !");
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> c=BuildContentFrom(fid,fileName);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TS> ret=BuildNewInstanceFromContent(c,fileName);
  ret->loadGlobals(fid);
  return ret.retn();
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::New(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> c=BuildContentFrom(fid,fileName,fieldName);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TS> ret=BuildNewInstanceFromContent(c,fileName);
  ret->loadGlobals(fid);
  return ret.retn();
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::New(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> c=BuildContentFrom(fid,fileName,fieldName,iteration,order);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TS> ret=BuildNewInstanceFromContent(c,fileName);
  ret->loadGlobals(fid);
  return ret.retn();
}

MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeField1TS::BuildContentFrom(med_idt fid, const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  med_field_type typcha;
  std::vector<std::string> infos;
  std::string dtunit;
  int iii=-1;
  int nbOfStep2=LocateField(fid,fileName,fieldName,iii,typcha,infos,dtunit);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> ret;
  switch(typcha)
    {
    case MED_FLOAT64:
      {
        ret=MEDFileField1TSWithoutSDA::New(fieldName,-1,iteration,order,std::vector<std::string>());
        break;
      }
    case MED_INT32:
      {
        ret=MEDFileIntField1TSWithoutSDA::New(fieldName,-1,iteration,order,std::vector<std::string>());
        break;
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileAnyTypeField1TS::BuildContentFrom(fileName,fieldName,iteration,order) : file \'" << fileName << "\' contains field with name \'" << fieldName << "\' but the type of field is not in [MED_FLOAT64, MED_INT32] !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    }
  ret->setDtUnit(dtunit.c_str());
  ret->getOrCreateAndGetArray()->setInfoAndChangeNbOfCompo(infos);
  //
  bool found=false;
  std::vector< std::pair<int,int> > dtits(nbOfStep2);
  for(int i=0;i<nbOfStep2 && !found;i++)
    {
      med_int numdt,numit;
      med_float dt;
      MEDfieldComputingStepInfo(fid,fieldName,i+1,&numdt,&numit,&dt);
      if(numdt==iteration && numit==order)
        {
          found=true;
          ret->_csit=i+1;
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
  ret->finishLoading(fid,*((const MEDFileAnyTypeField1TSWithoutSDA*)ret));
  return ret.retn();
}

MEDFileAnyTypeField1TS::MEDFileAnyTypeField1TS(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
try:MEDFileFieldGlobsReal(fileName)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  _content=BuildContentFrom(fid,fileName,fieldName,iteration,order);
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

/*!
 * This constructor is a shallow copy constructor. If \a shallowCopyOfContent is true the content of \a other is shallow copied.
 * If \a shallowCopyOfContent is false, \a other is taken to be the content of \a this.
 *
 * \warning this is a shallow copy constructor
 */
MEDFileAnyTypeField1TS::MEDFileAnyTypeField1TS(const MEDFileAnyTypeField1TSWithoutSDA& other, bool shallowCopyOfContent)
{
  if(!shallowCopyOfContent)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *otherPtr(&other);
      otherPtr->incrRef();
      _content=const_cast<MEDFileAnyTypeField1TSWithoutSDA *>(otherPtr);
    }
  else
    {
      _content=other.shallowCpy();
    }
}

int MEDFileAnyTypeField1TS::LocateField2(med_idt fid, const char *fileName, int fieldIdCFormat, bool checkFieldId, std::string& fieldName, med_field_type& typcha, std::vector<std::string>& infos, std::string& dtunitOut) throw(INTERP_KERNEL::Exception)
{
  if(checkFieldId)
    {
      int nbFields=MEDnField(fid);
      if(fieldIdCFormat>=nbFields)
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeField1TS::LocateField2(fileName) : in file \'" << fileName << "\' number of fields is " << nbFields << " ! Trying to request for id " << fieldIdCFormat << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  int ncomp=MEDfieldnComponent(fid,fieldIdCFormat+1);
  INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> dtunit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> nomMaa=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  med_bool localMesh;
  int nbOfStep;
  MEDfieldInfo(fid,fieldIdCFormat+1,nomcha,nomMaa,&localMesh,&typcha,comp,unit,dtunit,&nbOfStep);
  fieldName=MEDLoaderBase::buildStringFromFortran(nomcha,MED_NAME_SIZE);
  dtunitOut=MEDLoaderBase::buildStringFromFortran(dtunit,MED_LNAME_SIZE);
  infos.clear(); infos.resize(ncomp);
  for(int j=0;j<ncomp;j++)
    infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+j*MED_SNAME_SIZE,MED_SNAME_SIZE,(char *)unit+j*MED_SNAME_SIZE,MED_SNAME_SIZE);
  return nbOfStep;
}

/*!
 * This method throws an INTERP_KERNEL::Exception if \a fieldName field is not in file pointed by \a fid and with name \a fileName.
 * 
 * \param [out]
 * \return in case of success the number of time steps available for the field with name \a fieldName.
 */
int MEDFileAnyTypeField1TS::LocateField(med_idt fid, const char *fileName, const char *fieldName, int& posCFormat, med_field_type& typcha, std::vector<std::string>& infos, std::string& dtunitOut) throw(INTERP_KERNEL::Exception)
{
  int nbFields=MEDnField(fid);
  bool found=false;
  std::vector<std::string> fns(nbFields);
  int nbOfStep2=-1;
  for(int i=0;i<nbFields && !found;i++)
    {
      std::string tmp;
      nbOfStep2=LocateField2(fid,fileName,i,false,tmp,typcha,infos,dtunitOut);
      fns[i]=tmp;
      found=(tmp==fieldName);
      if(found)
        posCFormat=i;
    }
  if(!found)
    {
      std::ostringstream oss; oss << "No such field '" << fieldName << "' in file '" << fileName << "' ! Available fields are : ";
      for(std::vector<std::string>::const_iterator it=fns.begin();it!=fns.end();it++)
        oss << "\"" << *it << "\" ";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return nbOfStep2;
}

/*!
 * This method as MEDFileField1TSW::setLocNameOnLeaf, is dedicated for advanced user that a want a very fine control on their data structure
 * without overhead. This method can be called only regarding information returned by MEDFileField1TSWithoutSDA::getFieldSplitedByType or MEDFileField1TSWithoutSDA::getFieldSplitedByType2.
 * This method changes the attribute (here it's profile name) of the leaf datastructure (MEDFileFieldPerMeshPerTypePerDisc instance).
 * It is the responsability of the caller to invoke MEDFileFieldGlobs::appendProfile or MEDFileFieldGlobs::getProfile
 * to keep a valid instance.
 * If \b this do not have any leaf that correspond to the request of the input parameter (\b mName, \b typ, \b locId) an INTERP_KERNEL::Exception will be thrown.
 * If \b newPflName profile name does not already exist the profile with old name will be renamed with name \b newPflName.
 * If \b newPflName already exists and that \b forceRenameOnGlob is false (the default) an INTERP_KERNEL::Exception will be thrown to avoid big confusion. In this case the called should rename before the profile name with name \b newPflName.
 *
 * \param [in] mName specifies the underlying mesh name. This value can be pointer 0 for users that do not deal with fields on multi mesh.
 * \param [in] typ is for the geometric cell type (or INTERP_KERNEL::NORM_ERROR for node field) entry to find the right MEDFileFieldPerMeshPerTypePerDisc instance to set.
 * \param [in] locId is the localization id to find the right MEDFileFieldPerMeshPerTypePerDisc instance to set. It corresponds to the position of 
 *             \c pfls[std::distance(types.begin(),std::find(types.begin(),typ)] vector in MEDFileField1TSWithoutSDA::getFieldSplitedByType. For non gausspoints field users, the value is 0.
 * \param [in] newLocName is the new localization name.
 * \param [in] forceRenameOnGlob specifies the behaviour in case of profile \b newPflName already exists. If true, the renaming is done without check. It can lead to major bug.
 *             If false, an exception will be thrown to force user to change previously the name of the profile with name \b newPflName
 */
void MEDFileAnyTypeField1TS::setProfileNameOnLeaf(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const char *newPflName, bool forceRenameOnGlob) throw(INTERP_KERNEL::Exception)
{
  MEDFileFieldPerMeshPerTypePerDisc *disc=getLeafGivenMeshAndTypeAndLocId(mName,typ,locId);
  std::string oldPflName=disc->getProfile();
  std::vector<std::string> vv=getPflsReallyUsedMulti();
  int nbOfOcc=std::count(vv.begin(),vv.end(),oldPflName);
  if(forceRenameOnGlob || (!existsPfl(newPflName) && nbOfOcc==1))
    {
      disc->setProfile(newPflName);
      DataArrayInt *pfl=getProfile(oldPflName.c_str());
      pfl->setName(newPflName);
    }
  else
    {
      std::ostringstream oss; oss << "MEDFileField1TS::setProfileNameOnLeaf : Profile \"" << newPflName << "\" already exists or referenced more than one !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * This method as MEDFileField1TSW::setProfileNameOnLeaf, is dedicated for advanced user that a want a very fine control on their data structure
 * without overhead. This method can be called only regarding information returned by MEDFileField1TSWithoutSDA::getFieldSplitedByType or MEDFileField1TSWithoutSDA::getFieldSplitedByType2.
 * This method changes the attribute (here it's localization name) of the leaf datastructure (MEDFileFieldPerMeshPerTypePerDisc instance).
 * It is the responsability of the caller to invoke MEDFileFieldGlobs::appendProfile or MEDFileFieldGlobs::getProfile
 * to keep a valid instance.
 * If \b this do not have any leaf that correspond to the request of the input parameter (\b mName, \b typ, \b locId) an INTERP_KERNEL::Exception will be thrown.
 * This method is an extension of MEDFileField1TSWithoutSDA::setProfileNameOnLeafExt method because it performs a modification of global info.
 * If \b newLocName profile name does not already exist the localization with old name will be renamed with name \b newLocName.
 * If \b newLocName already exists an INTERP_KERNEL::Exception will be thrown to avoid big confusion. In this case the called should rename before the profile name with name \b newLocName.
 *
 * \param [in] mName specifies the underlying mesh name. This value can be pointer 0 for users that do not deal with fields on multi mesh.
 * \param [in] typ is for the geometric cell type (or INTERP_KERNEL::NORM_ERROR for node field) entry to find the right MEDFileFieldPerMeshPerTypePerDisc instance to set.
 * \param [in] locId is the localization id to find the right MEDFileFieldPerMeshPerTypePerDisc instance to set. It corresponds to the position of 
 *             \c pfls[std::distance(types.begin(),std::find(types.begin(),typ)] vector in MEDFileField1TSWithoutSDA::getFieldSplitedByType. For non gausspoints field users, the value is 0.
 * \param [in] newLocName is the new localization name.
 * \param [in] forceRenameOnGlob specifies the behaviour in case of profile \b newLocName already exists. If true, the renaming is done without check. It can lead to major bug.
 *             If false, an exception will be thrown to force user to change previously the name of the profile with name \b newLocName
 */
void MEDFileAnyTypeField1TS::setLocNameOnLeaf(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const char *newLocName, bool forceRenameOnGlob) throw(INTERP_KERNEL::Exception)
{
  MEDFileFieldPerMeshPerTypePerDisc *disc=getLeafGivenMeshAndTypeAndLocId(mName,typ,locId);
  std::string oldLocName=disc->getLocalization();
  std::vector<std::string> vv=getLocsReallyUsedMulti();
  int nbOfOcc=std::count(vv.begin(),vv.end(),oldLocName);
  if(forceRenameOnGlob || (!existsLoc(newLocName) && nbOfOcc==1))
    {
      disc->setLocalization(newLocName);
      MEDFileFieldLoc& loc=getLocalization(oldLocName.c_str());
      loc.setName(newLocName);
    }
  else
    {
      std::ostringstream oss; oss << "MEDFileField1TS::setLocNameOnLeaf : Localization \"" << newLocName << "\" already exists or referenced more than one !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeField1TS::contentNotNullBase() throw(INTERP_KERNEL::Exception)
{
  MEDFileAnyTypeField1TSWithoutSDA *ret=_content;
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TS : content is expected to be not null !");
  return ret;
}

const MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeField1TS::contentNotNullBase() const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA *ret=_content;
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TS : const content is expected to be not null !");
  return ret;
}

/*!
 * Writes \a this field into a MED file specified by its name.
 *  \param [in] fileName - the MED file name.
 *  \param [in] mode - the writing mode. For more on \a mode, see \ref AdvMEDLoaderBasics.
 * - 2 - erase; an existing file is removed.
 * - 1 - append; same data should not be present in an existing file.
 * - 0 - overwrite; same data present in an existing file is overwritten.
 *  \throw If the field name is not set.
 *  \throw If no field data is set.
 *  \throw If \a mode == 1 and the same data is present in an existing file.
 */
void MEDFileAnyTypeField1TS::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,medmod);
  writeLL(fid);
}

void MEDFileAnyTypeField1TS::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  int nbComp=getNumberOfComponents();
  INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  for(int i=0;i<nbComp;i++)
    {
      std::string info=getInfo()[i];
      std::string c,u;
      MEDLoaderBase::splitIntoNameAndUnit(info,c,u);
      MEDLoaderBase::safeStrCpy2(c.c_str(),MED_SNAME_SIZE,comp+i*MED_SNAME_SIZE,_too_long_str);
      MEDLoaderBase::safeStrCpy2(u.c_str(),MED_SNAME_SIZE,unit+i*MED_SNAME_SIZE,_too_long_str);
    }
  if(getName().empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TS::write : MED file does not accept field with empty name !");
  MEDfieldCr(fid,getName().c_str(),getMEDFileFieldType(),nbComp,comp,unit,getDtUnit().c_str(),getMeshName().c_str());
  writeGlobals(fid,*this);
  contentNotNullBase()->writeLL(fid,*this,*contentNotNullBase());
}

std::size_t MEDFileAnyTypeField1TS::getHeapMemorySize() const
{
  std::size_t ret=0;
  if((const MEDFileAnyTypeField1TSWithoutSDA *)_content)
    ret+=_content->getHeapMemorySize();
  return ret+MEDFileFieldGlobsReal::getHeapMemorySize();
}

/*!
 * Returns a string describing \a this field. This string is outputted 
 * by \c print Python command.
 */
std::string MEDFileAnyTypeField1TS::simpleRepr() const
{
  std::ostringstream oss;
  contentNotNullBase()->simpleRepr(0,oss,-1);
  simpleReprGlobs(oss);
  return oss.str();
}

/*!
 * This method returns all profiles whose name is non empty used.
 * \b WARNING If profile is used several times it will be reported \b only \b once.
 * To get non empty name profiles as time as they appear in \b this call MEDFileField1TS::getPflsReallyUsedMulti instead.
 */
std::vector<std::string> MEDFileAnyTypeField1TS::getPflsReallyUsed() const
{
  return contentNotNullBase()->getPflsReallyUsed2();
}

/*!
 * This method returns all localizations whose name is non empty used.
 * \b WARNING If localization is used several times it will be reported \b only \b once.
 */
std::vector<std::string> MEDFileAnyTypeField1TS::getLocsReallyUsed() const
{
  return contentNotNullBase()->getLocsReallyUsed2();
}

/*!
 * This method returns all profiles whose name is non empty used.
 * \b WARNING contrary to MEDFileField1TS::getPflsReallyUsed, if profile is used several times it will be reported as time as it appears.
 */
std::vector<std::string> MEDFileAnyTypeField1TS::getPflsReallyUsedMulti() const
{
  return contentNotNullBase()->getPflsReallyUsedMulti2();
}

/*!
 * This method returns all localizations whose name is non empty used.
 * \b WARNING contrary to MEDFileField1TS::getLocsReallyUsed if localization is used several times it will be reported as time as it appears.
 */
std::vector<std::string> MEDFileAnyTypeField1TS::getLocsReallyUsedMulti() const
{
  return contentNotNullBase()->getLocsReallyUsedMulti2();
}

void MEDFileAnyTypeField1TS::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  contentNotNullBase()->changePflsRefsNamesGen2(mapOfModif);
}

void MEDFileAnyTypeField1TS::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  contentNotNullBase()->changeLocsRefsNamesGen2(mapOfModif);
}

int MEDFileAnyTypeField1TS::getDimension() const
{
  return contentNotNullBase()->getDimension();
}

int MEDFileAnyTypeField1TS::getIteration() const
{
  return contentNotNullBase()->getIteration();
}

int MEDFileAnyTypeField1TS::getOrder() const
{
  return contentNotNullBase()->getOrder();
}

double MEDFileAnyTypeField1TS::getTime(int& iteration, int& order) const
{
  return contentNotNullBase()->getTime(iteration,order);
}

void MEDFileAnyTypeField1TS::setTime(int iteration, int order, double val)
{
  contentNotNullBase()->setTime(iteration,order,val);
}

std::string MEDFileAnyTypeField1TS::getName() const
{
  return contentNotNullBase()->getName();
}

void MEDFileAnyTypeField1TS::setName(const char *name)
{
  contentNotNullBase()->setName(name);
}

void MEDFileAnyTypeField1TS::simpleRepr(int bkOffset, std::ostream& oss, int f1tsId) const
{
  contentNotNullBase()->simpleRepr(bkOffset,oss,f1tsId);
}

std::string MEDFileAnyTypeField1TS::getDtUnit() const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getDtUnit();
}

void MEDFileAnyTypeField1TS::setDtUnit(const char *dtUnit) throw(INTERP_KERNEL::Exception)
{
  contentNotNullBase()->setDtUnit(dtUnit);
}

std::string MEDFileAnyTypeField1TS::getMeshName() const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getMeshName();
}

void MEDFileAnyTypeField1TS::setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception)
{
  contentNotNullBase()->setMeshName(newMeshName);
}

bool MEDFileAnyTypeField1TS::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->changeMeshNames(modifTab);
}

int MEDFileAnyTypeField1TS::getMeshIteration() const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getMeshIteration();
}

int MEDFileAnyTypeField1TS::getMeshOrder() const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getMeshOrder();
}

int MEDFileAnyTypeField1TS::getNumberOfComponents() const
{
  return contentNotNullBase()->getNumberOfComponents();
}

bool MEDFileAnyTypeField1TS::isDealingTS(int iteration, int order) const
{
  return contentNotNullBase()->isDealingTS(iteration,order);
}

std::pair<int,int> MEDFileAnyTypeField1TS::getDtIt() const
{
  return contentNotNullBase()->getDtIt();
}

void MEDFileAnyTypeField1TS::fillIteration(std::pair<int,int>& p) const
{
  contentNotNullBase()->fillIteration(p);
}

void MEDFileAnyTypeField1TS::fillTypesOfFieldAvailable(std::vector<TypeOfField>& types) const throw(INTERP_KERNEL::Exception)
{
  contentNotNullBase()->fillTypesOfFieldAvailable(types);
}

void MEDFileAnyTypeField1TS::setInfo(const std::vector<std::string>& infos) throw(INTERP_KERNEL::Exception)
{
  contentNotNullBase()->setInfo(infos);
}

const std::vector<std::string>& MEDFileAnyTypeField1TS::getInfo() const
{
  return contentNotNullBase()->getInfo();
}
std::vector<std::string>& MEDFileAnyTypeField1TS::getInfo()
{
  return contentNotNullBase()->getInfo();
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileAnyTypeField1TS::getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getLeafGivenMeshAndTypeAndLocId(mName,typ,locId);
}

const MEDFileFieldPerMeshPerTypePerDisc *MEDFileAnyTypeField1TS::getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getLeafGivenMeshAndTypeAndLocId(mName,typ,locId);
}

int MEDFileAnyTypeField1TS::getNonEmptyLevels(const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getNonEmptyLevels(mname,levs);
}

std::vector<TypeOfField> MEDFileAnyTypeField1TS::getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getTypesOfFieldAvailable();
}

std::vector< std::vector<std::pair<int,int> > > MEDFileAnyTypeField1TS::getFieldSplitedByType(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF,
                                                                                       std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getFieldSplitedByType(mname,types,typesF,pfls,locs);
}

/*!
 * This method returns as MEDFileAnyTypeField1TS new instances as number of components in \a this.
 * The returned instances are deep copy of \a this except that for globals that are share with those contained in \a this.
 * ** WARNING ** do no forget to rename the ouput instances to avoid to write n-times in the same MED file field !
 */
std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TS > > MEDFileAnyTypeField1TS::splitComponents() const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA *content(_content);
  if(!content)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TS::splitComponents : no content in this ! Unable to split components !");
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> > contentsSplit=content->splitComponents();
  std::size_t sz(contentsSplit.size());
  std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TS > > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=shallowCpy();
      ret[i]->_content=contentsSplit[i];
    }
  return ret;
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::deepCpy() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TS> ret=shallowCpy();
  if((const MEDFileAnyTypeField1TSWithoutSDA *)_content)
    ret->_content=_content->deepCpy();
  ret->deepCpyGlobs(*this);
  return ret.retn();
}

int MEDFileAnyTypeField1TS::copyTinyInfoFrom(const MEDCouplingFieldDouble *field, const DataArray *arr) throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->copyTinyInfoFrom(field,arr);
}

//= MEDFileField1TS

/*!
 * Returns a new instance of MEDFileField1TS holding data of the first time step of 
 * the first field that has been read from a specified MED file.
 *  \param [in] fileName - the name of the MED file to read.
 *  \return MEDFileField1TS * - a new instance of MEDFileFieldMultiTS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 *  \throw If reading the file fails.
 */
MEDFileField1TS *MEDFileField1TS::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> ret=new MEDFileField1TS(fileName);
  ret->contentNotNull();
  return ret.retn();
}

/*!
 * Returns a new instance of MEDFileField1TS holding data of the first time step of 
 * a given field that has been read from a specified MED file.
 *  \param [in] fileName - the name of the MED file to read.
 *  \param [in] fieldName - the name of the field to read.
 *  \return MEDFileField1TS * - a new instance of MEDFileFieldMultiTS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 *  \throw If reading the file fails.
 *  \throw If there is no field named \a fieldName in the file.
 */
MEDFileField1TS *MEDFileField1TS::New(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> ret=new MEDFileField1TS(fileName,fieldName);
  ret->contentNotNull();
  return ret.retn();
}

/*!
 * Returns a new instance of MEDFileField1TS holding data of a given time step of 
 * a given field that has been read from a specified MED file.
 *  \param [in] fileName - the name of the MED file to read.
 *  \param [in] fieldName - the name of the field to read.
 *  \param [in] iteration - the iteration number of a required time step.
 *  \param [in] order - the iteration order number of required time step.
 *  \return MEDFileField1TS * - a new instance of MEDFileFieldMultiTS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 *  \throw If reading the file fails.
 *  \throw If there is no field named \a fieldName in the file.
 *  \throw If the required time step is missing from the file.
 */
MEDFileField1TS *MEDFileField1TS::New(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> ret=new MEDFileField1TS(fileName,fieldName,iteration,order);
  ret->contentNotNull();
  return ret.retn();
}

/*!
 * Returns a new instance of MEDFileField1TS. If \a shallowCopyOfContent is true the content of \a other is shallow copied.
 * If \a shallowCopyOfContent is false, \a other is taken to be the content of \a this.
 *
 * Returns a new instance of MEDFileField1TS holding either a shallow copy
 * of a given MEDFileField1TSWithoutSDA ( \a other ) or \a other itself.
 * \warning this is a shallow copy constructor
 *  \param [in] other - a MEDFileField1TSWithoutSDA to copy.
 *  \param [in] shallowCopyOfContent - if \c true, a shallow copy of \a other is created.
 *  \return MEDFileField1TS * - a new instance of MEDFileField1TS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 */
MEDFileField1TS *MEDFileField1TS::New(const MEDFileField1TSWithoutSDA& other, bool shallowCopyOfContent)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> ret=new MEDFileField1TS(other,shallowCopyOfContent);
  ret->contentNotNull();
  return ret.retn();
}

/*!
 * Returns a new empty instance of MEDFileField1TS.
 *  \return MEDFileField1TS * - a new instance of MEDFileField1TS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 */
MEDFileField1TS *MEDFileField1TS::New()
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> ret=new MEDFileField1TS;
  ret->contentNotNull();
  return ret.retn();
}

const MEDFileField1TSWithoutSDA *MEDFileField1TS::contentNotNull() const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA *pt(_content);
  if(!pt)
    throw INTERP_KERNEL::Exception("MEDFileField1TS::contentNotNull : the content pointer is null !");
  const MEDFileField1TSWithoutSDA *ret=dynamic_cast<const MEDFileField1TSWithoutSDA *>(pt);
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileField1TS::contentNotNull : the content pointer is not null but it is not of type double ! Reason is maybe that the read field has not the type FLOAT64 !");
  return ret;
}

MEDFileField1TSWithoutSDA *MEDFileField1TS::contentNotNull() throw(INTERP_KERNEL::Exception)
{
  MEDFileAnyTypeField1TSWithoutSDA *pt(_content);
  if(!pt)
    throw INTERP_KERNEL::Exception("MEDFileField1TS::contentNotNull : the non const content pointer is null !");
  MEDFileField1TSWithoutSDA *ret=dynamic_cast<MEDFileField1TSWithoutSDA *>(pt);
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileField1TS::contentNotNull : the non const content pointer is not null but it is not of type double ! Reason is maybe that the read field has not the type FLOAT64 !");
  return ret;
}

void MEDFileField1TS::SetDataArrayDoubleInField(MEDCouplingFieldDouble *f, MEDCouplingAutoRefCountObjectPtr<DataArray>& arr) throw(INTERP_KERNEL::Exception)
{
  if(!f)
    throw INTERP_KERNEL::Exception("MEDFileField1TS::SetDataArrayDoubleInField : input field is NULL !");
  if(!((DataArray*)arr))
    throw INTERP_KERNEL::Exception("MEDFileField1TS::SetDataArrayDoubleInField : no array !");
  DataArrayDouble *arrOutC=dynamic_cast<DataArrayDouble *>((DataArray*)arr);
  if(!arrOutC)
    throw INTERP_KERNEL::Exception("MEDFileField1TS::SetDataArrayDoubleInField : mismatch between dataArrays type and MEDFileField1TS ! Expected double !");
  f->setArray(arrOutC);
}

DataArrayDouble *MEDFileField1TS::ReturnSafelyDataArrayDouble(MEDCouplingAutoRefCountObjectPtr<DataArray>& arr) throw(INTERP_KERNEL::Exception)
{
  if(!((DataArray*)arr))
    throw INTERP_KERNEL::Exception("MEDFileField1TS::ReturnSafelyDataArrayDouble : no array !");
  DataArrayDouble *arrOutC=dynamic_cast<DataArrayDouble *>((DataArray*)arr);
  if(!arrOutC)
    throw INTERP_KERNEL::Exception("MEDFileField1TS::ReturnSafelyDataArrayDouble : mismatch between dataArrays type and MEDFileField1TS ! Expected double !");
  arrOutC->incrRef();
  return arrOutC;
}

MEDFileField1TS::MEDFileField1TS(const char *fileName) throw(INTERP_KERNEL::Exception)
try:MEDFileAnyTypeField1TS(fileName)
{
}
catch(INTERP_KERNEL::Exception& e)
  { throw e; }

MEDFileField1TS::MEDFileField1TS(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
try:MEDFileAnyTypeField1TS(fileName,fieldName)
{
}
catch(INTERP_KERNEL::Exception& e)
  { throw e; }

MEDFileField1TS::MEDFileField1TS(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
try:MEDFileAnyTypeField1TS(fileName,fieldName,iteration,order)
{
}
catch(INTERP_KERNEL::Exception& e)
  { throw e; }

/*!
 * This constructor is a shallow copy constructor. If \a shallowCopyOfContent is true the content of \a other is shallow copied.
 * If \a shallowCopyOfContent is false, \a other is taken to be the content of \a this.
 *
 * \warning this is a shallow copy constructor
 */
MEDFileField1TS::MEDFileField1TS(const MEDFileField1TSWithoutSDA& other, bool shallowCopyOfContent)
try:MEDFileAnyTypeField1TS(other,shallowCopyOfContent)
{
}
catch(INTERP_KERNEL::Exception& e)
  { throw e; }

MEDFileField1TS::MEDFileField1TS()
{
  _content=new MEDFileField1TSWithoutSDA;
}

/*!
 * Returns a new MEDCouplingFieldDouble of a given type lying on
 * mesh entities of a given dimension of the first mesh in MED file. If \a this field 
 * has not been constructed via file reading, an exception is thrown.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of interest.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If \a this field has not been constructed via file reading.
 *  \throw If the MED file is not readable.
 *  \throw If there is no mesh in the MED file.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in the mesh.
 *  \throw If no field values of the given \a type or given \a meshDimRelToMax are available.
 *  \sa getFieldOnMeshAtLevel()
 */
MEDCouplingFieldDouble *MEDFileField1TS::getFieldAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  if(getFileName2().empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TS::getFieldAtLevel : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtLevel method instead !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arrOut;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=contentNotNull()->getFieldAtLevel(type,meshDimRelToMax,0,renumPol,this,arrOut,*contentNotNull());
  MEDFileField1TS::SetDataArrayDoubleInField(ret,arrOut);
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble of a given type lying on
 * the top level cells of the first mesh in MED file. If \a this field 
 * has not been constructed via file reading, an exception is thrown.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of interest.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If \a this field has not been constructed via file reading.
 *  \throw If the MED file is not readable.
 *  \throw If there is no mesh in the MED file.
 *  \throw If no field values of the given \a type.
 *  \throw If no field values lying on the top level support.
 *  \sa getFieldAtLevel()
 */
MEDCouplingFieldDouble *MEDFileField1TS::getFieldAtTopLevel(TypeOfField type, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  if(getFileName2().empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TS::getFieldAtTopLevel : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtTopLevel method instead !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arrOut;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=contentNotNull()->getFieldAtTopLevel(type,0,renumPol,this,arrOut,*contentNotNull());
  MEDFileField1TS::SetDataArrayDoubleInField(ret,arrOut);
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble of given type lying on a given mesh.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of the new field.
 *  \param [in] mesh - the supporting mesh.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If no field of \a this is lying on \a mesh.
 *  \throw If the mesh is empty.
 *  \throw If no field values of the given \a type are available.
 *  \sa getFieldAtLevel()
 *  \sa getFieldOnMeshAtLevel() 
 */
MEDCouplingFieldDouble *MEDFileField1TS::getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArray> arrOut;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=contentNotNull()->getFieldOnMeshAtLevel(type,renumPol,this,mesh,0,0,arrOut,*contentNotNull());
  MEDFileField1TS::SetDataArrayDoubleInField(ret,arrOut);
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble of a given type lying on a given support.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of interest.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [in] mesh - the supporting mesh.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in the mesh.
 *  \throw If no field of \a this is lying on \a mesh.
 *  \throw If no field values of the given \a type or given \a meshDimRelToMax are available.
 *  \sa getFieldAtLevel()
 *  \sa getFieldOnMeshAtLevel() 
 */
MEDCouplingFieldDouble *MEDFileField1TS::getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArray> arrOut;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=contentNotNull()->getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,this,mesh,arrOut,*contentNotNull());
  MEDFileField1TS::SetDataArrayDoubleInField(ret,arrOut);
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble of a given type lying on a given support.
 * This method is called "Old" because in MED3 norm a field has only one meshName
 * attached, so this method is for readers of MED2 files. If \a this field 
 * has not been constructed via file reading, an exception is thrown.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of interest.
 *  \param [in] mName - a name of the supporting mesh.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If the MED file is not readable.
 *  \throw If there is no mesh named \a mName in the MED file.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in the mesh.
 *  \throw If \a this field has not been constructed via file reading.
 *  \throw If no field of \a this is lying on the mesh named \a mName.
 *  \throw If no field values of the given \a type or given \a meshDimRelToMax are available.
 *  \sa getFieldAtLevel()
 */
MEDCouplingFieldDouble *MEDFileField1TS::getFieldAtLevelOld(TypeOfField type, const char *mname, int meshDimRelToMax, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  if(getFileName2().empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TS::getFieldAtLevelOld : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtLevel method instead !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arrOut;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=contentNotNull()->getFieldAtLevel(type,meshDimRelToMax,mname,renumPol,this,arrOut,*contentNotNull());
  MEDFileField1TS::SetDataArrayDoubleInField(ret,arrOut);
  return ret.retn();
}

/*!
 * Returns values and a profile of the field of a given type lying on a given support.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of the field.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [in] mesh - the supporting mesh.
 *  \param [out] pfl - a new instance of DataArrayInt holding ids of mesh entities the
 *          field of interest lies on. If the field lies on all entities of the given
 *          dimension, all ids in \a pfl are zero. The caller is to delete this array
 *          using decrRef() as it is no more needed.  
 *  \return DataArrayDouble * - a new instance of DataArrayDouble holding values of the
 *          field. The caller is to delete this array using decrRef() as it is no more needed.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in \a mesh.
 *  \throw If no field of \a this is lying on \a mesh.
 *  \throw If no field values of the given \a type or given \a meshDimRelToMax are available.
 */
DataArrayDouble *MEDFileField1TS::getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArray> ret=contentNotNull()->getFieldWithProfile(type,meshDimRelToMax,mesh,pfl,this,*contentNotNull());
  return MEDFileField1TS::ReturnSafelyDataArrayDouble(ret);
}

/*!
 * Adds a MEDCouplingFieldDouble to \a this. The underlying mesh of the given field is
 * checked if its elements are sorted suitable for writing to MED file ("STB" stands for
 * "Sort By Type"), if not, an exception is thrown. 
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] field - the field to add to \a this.
 *  \throw If the name of \a field is empty.
 *  \throw If the data array of \a field is not set.
 *  \throw If the data array is already allocated but has different number of components
 *         than \a field.
 *  \throw If the underlying mesh of \a field has no name.
 *  \throw If elements in the mesh are not in the order suitable for writing to the MED file.
 */
void MEDFileField1TS::setFieldNoProfileSBT(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception)
{
  setFileName("");
  contentNotNull()->setFieldNoProfileSBT(field,field->getArray(),*this,*contentNotNull());
}

/*!
 * Adds a MEDCouplingFieldDouble to \a this. Specified entities of a given dimension
 * of a given mesh are used as the support of the given field (a real support is not used). 
 * Elements of the given mesh must be sorted suitable for writing to MED file.
 * Order of underlying mesh entities of the given field specified by \a profile parameter
 * is not prescribed; this method permutes field values to have them sorted by element
 * type as required for writing to MED file. A new profile is added only if no equal
 * profile is missing.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] field - the field to add to \a this.
 *  \param [in] mesh - the supporting mesh of \a field.
 *  \param [in] meshDimRelToMax - a relative dimension of mesh entities \a field lies on.
 *  \param [in] profile - ids of mesh entities on which corresponding field values lie.
 *  \throw If either \a field or \a mesh or \a profile has an empty name.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in \a mesh.
 *  \throw If the data array of \a field is not set.
 *  \throw If the data array of \a this is already allocated but has different number of
 *         components than \a field.
 *  \throw If elements in \a mesh are not in the order suitable for writing to the MED file.
 *  \sa setFieldNoProfileSBT()
 */
void MEDFileField1TS::setFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception)
{
  setFileName("");
  contentNotNull()->setFieldProfile(field,field->getArray(),mesh,meshDimRelToMax,profile,*this,*contentNotNull());
}

MEDFileAnyTypeField1TS *MEDFileField1TS::shallowCpy() const throw(INTERP_KERNEL::Exception)
{
  return new MEDFileField1TS(*this);
}

DataArrayDouble *MEDFileField1TS::getUndergroundDataArray() const throw(INTERP_KERNEL::Exception)
{
  return contentNotNull()->getUndergroundDataArrayDouble();
}

DataArrayDouble *MEDFileField1TS::getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNull()->getUndergroundDataArrayDoubleExt(entries);
}

std::vector< std::vector<DataArrayDouble *> > MEDFileField1TS::getFieldSplitedByType2(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF,
                                                                                      std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNull()->getFieldSplitedByType2(mname,types,typesF,pfls,locs);
}

//= MEDFileIntField1TS

MEDFileIntField1TS *MEDFileIntField1TS::New()
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileIntField1TS> ret=new MEDFileIntField1TS;
  ret->contentNotNull();
  return ret.retn();
}

MEDFileIntField1TS *MEDFileIntField1TS::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileIntField1TS> ret=new MEDFileIntField1TS(fileName);
  ret->contentNotNull();
  return ret.retn();
}

MEDFileIntField1TS *MEDFileIntField1TS::New(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileIntField1TS> ret=new MEDFileIntField1TS(fileName,fieldName);
  ret->contentNotNull();
  return ret.retn();
}

MEDFileIntField1TS *MEDFileIntField1TS::New(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileIntField1TS> ret=new MEDFileIntField1TS(fileName,fieldName,iteration,order);
  ret->contentNotNull();
  return ret.retn();
}

MEDFileIntField1TS *MEDFileIntField1TS::New(const MEDFileIntField1TSWithoutSDA& other, bool shallowCopyOfContent)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileIntField1TS> ret=new MEDFileIntField1TS(other,shallowCopyOfContent);
  ret->contentNotNull();
  return ret.retn();
}

MEDFileIntField1TS::MEDFileIntField1TS()
{
  _content=new MEDFileIntField1TSWithoutSDA;
}

MEDFileIntField1TS::MEDFileIntField1TS(const char *fileName) throw(INTERP_KERNEL::Exception)
try:MEDFileAnyTypeField1TS(fileName)
{
}
catch(INTERP_KERNEL::Exception& e)
  { throw e; }

MEDFileIntField1TS::MEDFileIntField1TS(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
try:MEDFileAnyTypeField1TS(fileName,fieldName)
{
}
catch(INTERP_KERNEL::Exception& e)
  { throw e; }

MEDFileIntField1TS::MEDFileIntField1TS(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
try:MEDFileAnyTypeField1TS(fileName,fieldName,iteration,order)
{
}
catch(INTERP_KERNEL::Exception& e)
  { throw e; }

/*!
 * This constructor is a shallow copy constructor. If \a shallowCopyOfContent is true the content of \a other is shallow copied.
 * If \a shallowCopyOfContent is false, \a other is taken to be the content of \a this.
 *
 * \warning this is a shallow copy constructor
 */
MEDFileIntField1TS::MEDFileIntField1TS(const MEDFileIntField1TSWithoutSDA& other, bool shallowCopyOfContent):MEDFileAnyTypeField1TS(other,shallowCopyOfContent)
{
}

MEDFileAnyTypeField1TS *MEDFileIntField1TS::shallowCpy() const throw(INTERP_KERNEL::Exception)
{
  return new MEDFileIntField1TS(*this);
}

/*!
 * Adds a MEDCouplingFieldDouble to \a this. The underlying mesh of the given field is
 * checked if its elements are sorted suitable for writing to MED file ("STB" stands for
 * "Sort By Type"), if not, an exception is thrown. 
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] field - the field to add to \a this. The field double values are ignored.
 *  \param [in] arrOfVals - the values of the field \a field used.
 *  \throw If the name of \a field is empty.
 *  \throw If the data array of \a field is not set.
 *  \throw If the data array is already allocated but has different number of components
 *         than \a field.
 *  \throw If the underlying mesh of \a field has no name.
 *  \throw If elements in the mesh are not in the order suitable for writing to the MED file.
 */
void MEDFileIntField1TS::setFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals) throw(INTERP_KERNEL::Exception)
{
  setFileName("");
  contentNotNull()->setFieldNoProfileSBT(field,arrOfVals,*this,*contentNotNull());
}

/*!
 * Adds a MEDCouplingFieldDouble to \a this. Specified entities of a given dimension
 * of a given mesh are used as the support of the given field (a real support is not used). 
 * Elements of the given mesh must be sorted suitable for writing to MED file.
 * Order of underlying mesh entities of the given field specified by \a profile parameter
 * is not prescribed; this method permutes field values to have them sorted by element
 * type as required for writing to MED file. A new profile is added only if no equal
 * profile is missing.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] field - the field to add to \a this. The field double values are ignored.
 *  \param [in] arrOfVals - the values of the field \a field used.
 *  \param [in] mesh - the supporting mesh of \a field.
 *  \param [in] meshDimRelToMax - a relative dimension of mesh entities \a field lies on.
 *  \param [in] profile - ids of mesh entities on which corresponding field values lie.
 *  \throw If either \a field or \a mesh or \a profile has an empty name.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in \a mesh.
 *  \throw If the data array of \a field is not set.
 *  \throw If the data array of \a this is already allocated but has different number of
 *         components than \a field.
 *  \throw If elements in \a mesh are not in the order suitable for writing to the MED file.
 *  \sa setFieldNoProfileSBT()
 */
void MEDFileIntField1TS::setFieldProfile(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception)
{
  setFileName("");
  contentNotNull()->setFieldProfile(field,arrOfVals,mesh,meshDimRelToMax,profile,*this,*contentNotNull());
}

const MEDFileIntField1TSWithoutSDA *MEDFileIntField1TS::contentNotNull() const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA *pt(_content);
  if(!pt)
    throw INTERP_KERNEL::Exception("MEDFileIntField1TS::contentNotNull : the content pointer is null !");
  const MEDFileIntField1TSWithoutSDA *ret=dynamic_cast<const MEDFileIntField1TSWithoutSDA *>(pt);
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileIntField1TS::contentNotNull : the content pointer is not null but it is not of type int32 ! Reason is maybe that the read field has not the type INT32 !");
  return ret;
}

MEDCouplingFieldDouble *MEDFileIntField1TS::getFieldAtLevel(TypeOfField type, int meshDimRelToMax, DataArrayInt* &arrOut, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  if(getFileName2().empty())
    throw INTERP_KERNEL::Exception("MEDFileIntField1TS::getFieldAtLevel : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtLevel method instead !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arrOut2;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=contentNotNull()->getFieldAtLevel(type,meshDimRelToMax,0,renumPol,this,arrOut2,*contentNotNull());
  DataArrayInt *arrOutC=dynamic_cast<DataArrayInt *>((DataArray *)arrOut2);
  if(!arrOutC)
    throw INTERP_KERNEL::Exception("MEDFileIntField1TS::getFieldAtLevelOld : mismatch between dataArrays type and MEDFileIntField1TS ! Expected int32 !");
  arrOut=arrOutC;
  return ret.retn();
}

DataArrayInt *MEDFileIntField1TS::ReturnSafelyDataArrayInt(MEDCouplingAutoRefCountObjectPtr<DataArray>& arr) throw(INTERP_KERNEL::Exception)
{
  if(!((DataArray *)arr))
    throw INTERP_KERNEL::Exception("MEDFileIntField1TS::ReturnSafelyDataArrayInt : input DataArray is NULL !");
  DataArrayInt *arrC=dynamic_cast<DataArrayInt *>((DataArray *)arr);
  if(!arrC)
    throw INTERP_KERNEL::Exception("MEDFileIntField1TS::ReturnSafelyDataArrayInt : input DataArray is not of type INT32 !");
  arrC->incrRef();
  return arrC;
}

/*!
 * Returns a new MEDCouplingFieldDouble of a given type lying on
 * the top level cells of the first mesh in MED file. If \a this field 
 * has not been constructed via file reading, an exception is thrown.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of interest.
 *  \param [out] arrOut - the DataArrayInt containing values of field.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If \a this field has not been constructed via file reading.
 *  \throw If the MED file is not readable.
 *  \throw If there is no mesh in the MED file.
 *  \throw If no field values of the given \a type.
 *  \throw If no field values lying on the top level support.
 *  \sa getFieldAtLevel()
 */
MEDCouplingFieldDouble *MEDFileIntField1TS::getFieldAtTopLevel(TypeOfField type, DataArrayInt* &arrOut, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  if(getFileName2().empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TS::getFieldAtTopLevel : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtTopLevel method instead !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arr;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=contentNotNull()->getFieldAtTopLevel(type,0,renumPol,this,arr,*contentNotNull());
  arrOut=MEDFileIntField1TS::ReturnSafelyDataArrayInt(arr);
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble of given type lying on a given mesh.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of the new field.
 *  \param [in] mesh - the supporting mesh.
 *  \param [out] arrOut - the DataArrayInt containing values of field.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If no field of \a this is lying on \a mesh.
 *  \throw If the mesh is empty.
 *  \throw If no field values of the given \a type are available.
 *  \sa getFieldAtLevel()
 *  \sa getFieldOnMeshAtLevel() 
 */
MEDCouplingFieldDouble *MEDFileIntField1TS::getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, DataArrayInt* &arrOut, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArray> arr;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=contentNotNull()->getFieldOnMeshAtLevel(type,renumPol,this,mesh,0,0,arr,*contentNotNull());
  arrOut=MEDFileIntField1TS::ReturnSafelyDataArrayInt(arr);
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble of a given type lying on a given support.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of interest.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [out] arrOut - the DataArrayInt containing values of field.
 *  \param [in] mesh - the supporting mesh.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in the mesh.
 *  \throw If no field of \a this is lying on \a mesh.
 *  \throw If no field values of the given \a type or given \a meshDimRelToMax are available.
 *  \sa getFieldAtLevel()
 *  \sa getFieldOnMeshAtLevel() 
 */
MEDCouplingFieldDouble *MEDFileIntField1TS::getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt* &arrOut, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArray> arr;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=contentNotNull()->getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,this,mesh,arr,*contentNotNull());
  arrOut=MEDFileIntField1TS::ReturnSafelyDataArrayInt(arr);
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble of a given type lying on a given support.
 * This method is called "Old" because in MED3 norm a field has only one meshName
 * attached, so this method is for readers of MED2 files. If \a this field 
 * has not been constructed via file reading, an exception is thrown.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of interest.
 *  \param [in] mName - a name of the supporting mesh.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [out] arrOut - the DataArrayInt containing values of field.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If the MED file is not readable.
 *  \throw If there is no mesh named \a mName in the MED file.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in the mesh.
 *  \throw If \a this field has not been constructed via file reading.
 *  \throw If no field of \a this is lying on the mesh named \a mName.
 *  \throw If no field values of the given \a type or given \a meshDimRelToMax are available.
 *  \sa getFieldAtLevel()
 */
MEDCouplingFieldDouble *MEDFileIntField1TS::getFieldAtLevelOld(TypeOfField type, const char *mname, int meshDimRelToMax, DataArrayInt* &arrOut, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  if(getFileName2().empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TS::getFieldAtLevelOld : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtLevel method instead !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arr;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=contentNotNull()->getFieldAtLevel(type,meshDimRelToMax,mname,renumPol,this,arr,*contentNotNull());
  arrOut=MEDFileIntField1TS::ReturnSafelyDataArrayInt(arr);
  return ret.retn();
}

/*!
 * Returns values and a profile of the field of a given type lying on a given support.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of the field.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [in] mesh - the supporting mesh.
 *  \param [out] pfl - a new instance of DataArrayInt holding ids of mesh entities the
 *          field of interest lies on. If the field lies on all entities of the given
 *          dimension, all ids in \a pfl are zero. The caller is to delete this array
 *          using decrRef() as it is no more needed.  
 *  \return DataArrayInt * - a new instance of DataArrayInt holding values of the
 *          field. The caller is to delete this array using decrRef() as it is no more needed.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in \a mesh.
 *  \throw If no field of \a this is lying on \a mesh.
 *  \throw If no field values of the given \a type or given \a meshDimRelToMax are available.
 */
DataArrayInt *MEDFileIntField1TS::getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArray> arr=contentNotNull()->getFieldWithProfile(type,meshDimRelToMax,mesh,pfl,this,*contentNotNull());
  return MEDFileIntField1TS::ReturnSafelyDataArrayInt(arr);
}

MEDFileIntField1TSWithoutSDA *MEDFileIntField1TS::contentNotNull() throw(INTERP_KERNEL::Exception)
{
  MEDFileAnyTypeField1TSWithoutSDA *pt(_content);
  if(!pt)
    throw INTERP_KERNEL::Exception("MEDFileIntField1TS::contentNotNull : the non const content pointer is null !");
  MEDFileIntField1TSWithoutSDA *ret=dynamic_cast<MEDFileIntField1TSWithoutSDA *>(pt);
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileIntField1TS::contentNotNull : the non const content pointer is not null but it is not of type int32 ! Reason is maybe that the read field has not the type INT32 !");
  return ret;
}

DataArrayInt *MEDFileIntField1TS::getUndergroundDataArray() const throw(INTERP_KERNEL::Exception)
{
  return contentNotNull()->getUndergroundDataArrayInt();
}

//= MEDFileAnyTypeFieldMultiTSWithoutSDA

MEDFileAnyTypeFieldMultiTSWithoutSDA::MEDFileAnyTypeFieldMultiTSWithoutSDA()
{
}

MEDFileAnyTypeFieldMultiTSWithoutSDA::MEDFileAnyTypeFieldMultiTSWithoutSDA(const char *fieldName):MEDFileFieldNameScope(fieldName)
{
}

/*!
 * \param [in] fieldId field id in C mode
 */
MEDFileAnyTypeFieldMultiTSWithoutSDA::MEDFileAnyTypeFieldMultiTSWithoutSDA(med_idt fid, int fieldId) throw(INTERP_KERNEL::Exception)
{
  med_field_type typcha;
  std::string dtunitOut;
  int nbOfStep=MEDFileAnyTypeField1TS::LocateField2(fid,"",fieldId,false,_name,typcha,_infos,dtunitOut);
  setDtUnit(dtunitOut.c_str());
  finishLoading(fid,nbOfStep,typcha);
}

MEDFileAnyTypeFieldMultiTSWithoutSDA::MEDFileAnyTypeFieldMultiTSWithoutSDA(med_idt fid, const char *fieldName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit) throw(INTERP_KERNEL::Exception)
try:MEDFileFieldNameScope(fieldName),_infos(infos)
{
  setDtUnit(dtunit.c_str());
  finishLoading(fid,nbOfStep,fieldTyp);
}
catch(INTERP_KERNEL::Exception& e)
{
  throw e;
}

std::size_t MEDFileAnyTypeFieldMultiTSWithoutSDA::getHeapMemorySize() const
{
  std::size_t ret=_name.capacity()+_infos.capacity()*sizeof(std::string)+_time_steps.capacity()*sizeof(MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA>);
  for(std::vector<std::string>::const_iterator it=_infos.begin();it!=_infos.end();it++)
    ret+=(*it).capacity();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((const MEDFileAnyTypeField1TSWithoutSDA *)(*it))
      ret+=(*it)->getHeapMemorySize();
  return ret;
}

/*!
 * If one of the id in [ \a startIds , \a endIds ) points to a null element, there is not throw. Simply, this empty element is added as if it were not
 * NULL.
 */
MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::buildFromTimeStepIds(const int *startIds, const int *endIds) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> ret=createNew();
  ret->setInfo(_infos);
  int sz=(int)_time_steps.size();
  for(const int *id=startIds;id!=endIds;id++)
    {
      if(*id>=0 && *id<sz)
        {
          const MEDFileAnyTypeField1TSWithoutSDA *tse=_time_steps[*id];
          MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> tse2;
          if(tse)
            {
              tse->incrRef();
              tse2=(const_cast<MEDFileAnyTypeField1TSWithoutSDA *>(tse));
            }
          ret->pushBackTimeStep(tse2);
        }
      else
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::buildFromTimeStepIds : At pos #" << std::distance(startIds,id) << " value is " << *id;
          oss << " ! Should be in [0," << sz << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  if(ret->getNumberOfTS()>0)
    ret->synchronizeNameScope();
  ret->copyNameScope(*this);
  return ret.retn();
}

/*!
 * If one of the id in the input range points to a null element, there is not throw. Simply, this empty element is added as if it were not
 * NULL.
 */
MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::buildFromTimeStepIds2(int bg, int end, int step) const throw(INTERP_KERNEL::Exception)
{
  static const char msg[]="MEDFileAnyTypeFieldMultiTSWithoutSDA::buildFromTimeStepIds2";
  int nbOfEntriesToKeep=DataArrayInt::GetNumberOfItemGivenBESRelative(bg,end,step,msg);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> ret=createNew();
  ret->setInfo(_infos);
  int sz=(int)_time_steps.size();
  int j=bg;
  for(int i=0;i<nbOfEntriesToKeep;i++,j+=step)
    {
      if(j>=0 && j<sz)
        {
          const MEDFileAnyTypeField1TSWithoutSDA *tse=_time_steps[j];
          MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> tse2;
          if(tse)
            {
              tse->incrRef();
              tse2=(const_cast<MEDFileAnyTypeField1TSWithoutSDA *>(tse));
            }
          ret->pushBackTimeStep(tse2);
        }
      else
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::buildFromTimeStepIds : At pos #" << i << " value is " << j;
          oss << " ! Should be in [0," << sz << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  if(ret->getNumberOfTS()>0)
    ret->synchronizeNameScope();
  ret->copyNameScope(*this);
  return ret.retn();
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::partOfThisLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const throw(INTERP_KERNEL::Exception)
{
  int id=0;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ids=DataArrayInt::New(); ids->alloc(0,1);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,id++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *cur(*it);
      if(!cur)
        continue;
      std::pair<int,int> p(cur->getIteration(),cur->getOrder());
      if(std::find(timeSteps.begin(),timeSteps.end(),p)!=timeSteps.end())
        ids->pushBackSilent(id);
    }
  return buildFromTimeStepIds(ids->begin(),ids->end());
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::partOfThisNotLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const throw(INTERP_KERNEL::Exception)
{
  int id=0;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ids=DataArrayInt::New(); ids->alloc(0,1);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,id++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *cur(*it);
      if(!cur)
        continue;
      std::pair<int,int> p(cur->getIteration(),cur->getOrder());
      if(std::find(timeSteps.begin(),timeSteps.end(),p)==timeSteps.end())
        ids->pushBackSilent(id);
    }
  return buildFromTimeStepIds(ids->begin(),ids->end());
}

const std::vector<std::string>& MEDFileAnyTypeFieldMultiTSWithoutSDA::getInfo() const throw(INTERP_KERNEL::Exception)
{
  return _infos;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::setInfo(const std::vector<std::string>& info) throw(INTERP_KERNEL::Exception)
{
  _infos=info;
}

int MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepPos(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  int ret=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,ret++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *pt(*it);
      if(pt->isDealingTS(iteration,order))
        return ret;
    }
  std::ostringstream oss; oss << "MEDFileFieldMultiTS::getTimeStepPos : Muli timestep field on time (" << iteration << "," << order << ") does not exist ! Available (iteration,order) are :\n";
  std::vector< std::pair<int,int> > vp=getIterations();
  for(std::vector< std::pair<int,int> >::const_iterator it2=vp.begin();it2!=vp.end();it2++)
    oss << "(" << (*it2).first << "," << (*it2).second << ") ";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

const MEDFileAnyTypeField1TSWithoutSDA& MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepEntry(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  return *_time_steps[getTimeStepPos(iteration,order)];
}

MEDFileAnyTypeField1TSWithoutSDA& MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepEntry(int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  return *_time_steps[getTimeStepPos(iteration,order)];
}

std::string MEDFileAnyTypeFieldMultiTSWithoutSDA::getMeshName() const throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::getMeshName : not time steps !");
  return _time_steps[0]->getMeshName();
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception)
{
  std::string oldName(getMeshName());
  std::vector< std::pair<std::string,std::string> > v(1);
  v[0].first=oldName; v[0].second=newMeshName;
  changeMeshNames(v);
}

bool MEDFileAnyTypeFieldMultiTSWithoutSDA::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception)
{
  bool ret=false;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileAnyTypeField1TSWithoutSDA *cur(*it);
      if(cur)
        ret=cur->changeMeshNames(modifTab) || ret;
    }
  return ret;
}

/*!
 * See doc at MEDFileField1TSWithoutSDA::getUndergroundDataArray
 */
DataArray *MEDFileAnyTypeFieldMultiTSWithoutSDA::getUndergroundDataArray(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  return getTimeStepEntry(iteration,order).getUndergroundDataArray();
}

/*!
 * See doc at MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt
 */
DataArray *MEDFileAnyTypeFieldMultiTSWithoutSDA::getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception)
{
  return getTimeStepEntry(iteration,order).getUndergroundDataArrayExt(entries);
}

bool MEDFileAnyTypeFieldMultiTSWithoutSDA::renumberEntitiesLyingOnMesh(const char *meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N,
                                                                MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  bool ret=false;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileAnyTypeField1TSWithoutSDA *f1ts(*it);
      if(f1ts)
        ret=f1ts->renumberEntitiesLyingOnMesh(meshName,oldCode,newCode,renumO2N,glob) || ret;
    }
  return ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const
{
  std::string startLine(bkOffset,' ');
  oss << startLine << "Field multi time steps [Type=" << getTypeStr() << "]";
  if(fmtsId>=0)
    oss << " (" << fmtsId << ")";
  oss << " has the following name: \"" << _name << "\"." << std::endl;
  oss << startLine << "Field multi time steps has " << _infos.size() << " components with the following infos :" << std::endl;
  for(std::vector<std::string>::const_iterator it=_infos.begin();it!=_infos.end();it++)
    {
      oss << startLine << "  -  \"" << *it << "\"" << std::endl;
    }
  int i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,i++)
    {
      std::string chapter(17,'0'+i);
      oss << startLine << chapter << std::endl;
      const MEDFileAnyTypeField1TSWithoutSDA *cur=(*it);
      if(cur)
        cur->simpleRepr(bkOffset+2,oss,i);
      else
        oss << startLine << "  Field on one time step #" << i << " is not defined !" << std::endl;
      oss << startLine << chapter << std::endl;
    }
}

std::vector< std::pair<int,int> > MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeSteps(std::vector<double>& ret1) const throw(INTERP_KERNEL::Exception)
{
  std::size_t sz=_time_steps.size();
  std::vector< std::pair<int,int> > ret(sz);
  ret1.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *f1ts=_time_steps[i];
      if(f1ts)
        {
          ret1[i]=f1ts->getTime(ret[i].first,ret[i].second);
        }
      else
        {
          std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::getTimeSteps : At rank #" << i << " time step is not defined. Invoke eraseEmptyTS method !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::pushBackTimeStep(MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA>& tse) throw(INTERP_KERNEL::Exception)
{
  MEDFileAnyTypeField1TSWithoutSDA *tse2(tse);
  if(!tse2)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::pushBackTimeStep : input content object is null !");
  checkCoherencyOfType(tse2);
  if(_time_steps.empty())
    {
      setName(tse2->getName().c_str());
      setInfo(tse2->getInfo());
    }
  checkThatComponentsMatch(tse2->getInfo());
  _time_steps.push_back(tse);
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::synchronizeNameScope() throw(INTERP_KERNEL::Exception)
{
  std::size_t nbOfCompo=_infos.size();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileAnyTypeField1TSWithoutSDA *cur=(*it);
      if(cur)
        {
          if((cur->getInfo()).size()!=nbOfCompo)
            {
              std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::synchronizeNameScope : Mismatch in the number of components of parts ! Should be " << nbOfCompo;
              oss << " ! but the field at iteration=" << cur->getIteration() << " order=" << cur->getOrder() << " has " << (cur->getInfo()).size() << " components !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          cur->copyNameScope(*this);
        }
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::finishLoading(med_idt fid, int nbPdt, med_field_type fieldTyp) throw(INTERP_KERNEL::Exception)
{
  _time_steps.resize(nbPdt);
  for(int i=0;i<nbPdt;i++)
    {
      std::vector< std::pair<int,int> > ts;
      med_int numdt=0,numo=0;
      med_int meshIt=0,meshOrder=0;
      med_float dt=0.0;
      MEDfieldComputingStepMeshInfo(fid,_name.c_str(),i+1,&numdt,&numo,&dt,&meshIt,&meshOrder);
      switch(fieldTyp)
        {
        case MED_FLOAT64:
          {
            _time_steps[i]=MEDFileField1TSWithoutSDA::New(_name.c_str(),i+1,numdt,numo,_infos);
            break;
          }
        case MED_INT32:
          {
            _time_steps[i]=MEDFileIntField1TSWithoutSDA::New(_name.c_str(),i+1,numdt,numo,_infos);
            break;
          }
        default:
          throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::finishLoading : managed field type are : FLOAT64, INT32 !");
        }
      _time_steps[i]->finishLoading(fid,*this);
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::writeLL(med_idt fid, const MEDFileWritable& opts) const throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::writeLL : no time steps set !");
  checkThatNbOfCompoOfTSMatchThis();
  std::vector<std::string> infos(getInfo());
  int nbComp=infos.size();
  INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  for(int i=0;i<nbComp;i++)
    {
      std::string info=infos[i];
      std::string c,u;
      MEDLoaderBase::splitIntoNameAndUnit(info,c,u);
      MEDLoaderBase::safeStrCpy2(c.c_str(),MED_SNAME_SIZE,comp+i*MED_SNAME_SIZE,opts.getTooLongStrPolicy());
      MEDLoaderBase::safeStrCpy2(u.c_str(),MED_SNAME_SIZE,unit+i*MED_SNAME_SIZE,opts.getTooLongStrPolicy());
    }
  if(_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::write : MED file does not accept field with empty name !");
  MEDfieldCr(fid,_name.c_str(),getMEDFileFieldType(),nbComp,comp,unit,getDtUnit().c_str(),getMeshName().c_str());
  int nbOfTS=_time_steps.size();
  for(int i=0;i<nbOfTS;i++)
    _time_steps[i]->writeLL(fid,opts,*this);
}

int MEDFileAnyTypeFieldMultiTSWithoutSDA::getNumberOfTS() const
{
  return _time_steps.size();
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::eraseEmptyTS() throw(INTERP_KERNEL::Exception)
{
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA>  > newTS;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *tmp=(*it);
      if(tmp)
        newTS.push_back(*it);
    }
  _time_steps=newTS;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::eraseTimeStepIds(const int *startIds, const int *endIds) throw(INTERP_KERNEL::Exception)
{
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> > newTS;
  int maxId=(int)_time_steps.size();
  int ii=0;
  std::set<int> idsToDel;
  for(const int *id=startIds;id!=endIds;id++,ii++)
    {
      if(*id>=0 && *id<maxId)
        {
          idsToDel.insert(*id);
        }
      else
        {
          std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::eraseTimeStepIds : At pos #" << ii << " request for id=" << *id << " not in [0," << maxId << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  for(int iii=0;iii<maxId;iii++)
    if(idsToDel.find(iii)==idsToDel.end())
      newTS.push_back(_time_steps[iii]);
  _time_steps=newTS;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::eraseTimeStepIds2(int bg, int end, int step) throw(INTERP_KERNEL::Exception)
{
  static const char msg[]="MEDFileAnyTypeFieldMultiTSWithoutSDA::eraseTimeStepIds2";
  int nbOfEntriesToKill=DataArrayInt::GetNumberOfItemGivenBESRelative(bg,end,step,msg);
  if(nbOfEntriesToKill==0)
    return ;
  std::size_t sz=_time_steps.size();
  std::vector<bool> b(sz,true);
  int j=bg;
  for(int i=0;i<nbOfEntriesToKill;i++,j+=step)
    b[j]=false;
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> > newTS;
  for(std::size_t i=0;i<sz;i++)
    if(b[i])
      newTS.push_back(_time_steps[i]);
  _time_steps=newTS;
}

int MEDFileAnyTypeFieldMultiTSWithoutSDA::getPosOfTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  int ret=0;
  std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::getPosOfTimeStep : No such time step (" << iteration << "," << order << ") !\nPossibilities are : "; 
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,ret++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *tmp(*it);
      if(tmp)
        {
          int it2,ord;
          tmp->getTime(it2,ord);
          if(it2==iteration && order==ord)
            return ret;
          else
            oss << "(" << it2 << ","  << ord << "), ";
        }
    }
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

int MEDFileAnyTypeFieldMultiTSWithoutSDA::getPosGivenTime(double time, double eps) const throw(INTERP_KERNEL::Exception)
{
  int ret=0;
  std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::getPosGivenTime : No such time step " << time << "! \nPossibilities are : ";
  oss.precision(15);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,ret++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *tmp(*it);
      if(tmp)
        {
          int it2,ord;
          double ti=tmp->getTime(it2,ord);
          if(fabs(time-ti)<eps)
            return ret;
          else
            oss << ti << ", ";
        }
    }
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

std::vector< std::pair<int,int> > MEDFileAnyTypeFieldMultiTSWithoutSDA::getIterations() const
{
  int lgth=_time_steps.size();
  std::vector< std::pair<int,int> > ret(lgth);
  for(int i=0;i<lgth;i++)
    _time_steps[i]->fillIteration(ret[i]);
  return ret;
}

/*!
 * This method has 3 inputs 'iteration' 'order' 'mname'. 'mname' can be null if the user is the general case where there is only one meshName lying on 'this'
 * This method returns two things.
 * - The absolute dimension of 'this' in first parameter. 
 * - The available ext levels relative to the absolute dimension returned in first parameter. These relative levels are relative
 *   to the first output parameter. The values in 'levs' will be returned in decreasing order.
 *
 * This method is designed for MEDFileFieldMultiTS instances that have a discritization ON_CELLS, ON_GAUSS_NE and ON_GAUSS.
 * Only these 3 discretizations will be taken into account here.
 *
 * If 'this' is empty this method will throw an INTERP_KERNEL::Exception.
 * If there is \b only node fields defined in 'this' -1 is returned and 'levs' output parameter will be empty. In this
 * case the caller has to know the underlying mesh it refers to. By defaut it is the level 0 of the corresponding mesh.
 *
 * This method is usefull to make the link between meshDimension of the underlying mesh in 'this' and the levels on 'this'.
 * It is possible (even if it is not common) that the highest level in 'this' were not equal to the meshDimension of the underlying mesh in 'this'.
 * 
 * Let's consider the typical following case :
 * - a mesh 'm1' has a meshDimension 3 and has the following non empty levels
 * [0,-1,-2] for example 'm1' lies on TETRA4, HEXA8 TRI3 and SEG2
 * - 'f1' lies on 'm1' and is defined on 3D and 1D cells for example
 *   TETRA4 and SEG2
 * - 'f2' lies on 'm1' too and is defined on 2D and 1D cells for example TRI3 and SEG2
 *
 * In this case f1->getNonEmptyLevelsExt will return (3,[0,-2]) and f2->getNonEmptyLevelsExt will return (2,[0,-1])
 * 
 * To retrieve the highest level of f1 it should be done, f1->getFieldAtLevel(ON_CELLS,3-3+0);//absDim-meshDim+relativeLev
 * To retrieve the lowest level of f1 it should be done, f1->getFieldAtLevel(ON_CELLS,3-3+(-2));//absDim-meshDim+relativeLev
 * To retrieve the highest level of f2 it should be done, f1->getFieldAtLevel(ON_CELLS,2-3+0);//absDim-meshDim+relativeLev
 * To retrieve the lowest level of f2 it should be done, f1->getFieldAtLevel(ON_CELLS,2-3+(-1));//absDim-meshDim+relativeLev
 */
int MEDFileAnyTypeFieldMultiTSWithoutSDA::getNonEmptyLevels(int iteration, int order, const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception)
{
  return getTimeStepEntry(iteration,order).getNonEmptyLevels(mname,levs);
}

const MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2(int pos) const throw(INTERP_KERNEL::Exception)
{
  if(pos<0 || pos>=(int)_time_steps.size())
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << " whereas should be in [0," << _time_steps.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const MEDFileAnyTypeField1TSWithoutSDA *item=_time_steps[pos];
  if(item==0)
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << ", this pos id exists but the underlying Field1TS is null !";
      oss << "\nTry to use following method eraseEmptyTS !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return item;
}

MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2(int pos) throw(INTERP_KERNEL::Exception)
{
  if(pos<0 || pos>=(int)_time_steps.size())
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << " whereas should be in [0," << _time_steps.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDFileAnyTypeField1TSWithoutSDA *item=_time_steps[pos];
  if(item==0)
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << ", this pos id exists but the underlying Field1TS is null !";
      oss << "\nTry to use following method eraseEmptyTS !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return item;
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTSWithoutSDA::getPflsReallyUsed2() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
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

std::vector<std::string> MEDFileAnyTypeFieldMultiTSWithoutSDA::getLocsReallyUsed2() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
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

std::vector<std::string> MEDFileAnyTypeFieldMultiTSWithoutSDA::getPflsReallyUsedMulti2() const
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsedMulti2();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTSWithoutSDA::getLocsReallyUsedMulti2() const
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsedMulti2();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::changePflsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TSWithoutSDA > >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    (*it)->changePflsRefsNamesGen2(mapOfModif);
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::changeLocsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TSWithoutSDA > >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    (*it)->changeLocsRefsNamesGen2(mapOfModif);
}

std::vector< std::vector<TypeOfField> > MEDFileAnyTypeFieldMultiTSWithoutSDA::getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception)
{
  int lgth=_time_steps.size();
  std::vector< std::vector<TypeOfField> > ret(lgth);
  for(int i=0;i<lgth;i++)
    _time_steps[i]->fillTypesOfFieldAvailable(ret[i]);
  return ret;
}

/*!
 * entry point for users that want to iterate into MEDFile DataStructure without any overhead.
 */
std::vector< std::vector< std::pair<int,int> > > MEDFileAnyTypeFieldMultiTSWithoutSDA::getFieldSplitedByType(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  return getTimeStepEntry(iteration,order).getFieldSplitedByType(mname,types,typesF,pfls,locs);
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::deepCpy() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> ret=shallowCpy();
  std::size_t i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,i++)
    {
      if((const MEDFileAnyTypeField1TSWithoutSDA *)*it)
        ret->_time_steps[i]=(*it)->deepCpy();
    }
  return ret.retn();
}

std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> > MEDFileAnyTypeFieldMultiTSWithoutSDA::splitComponents() const throw(INTERP_KERNEL::Exception)
{
  std::size_t sz(_infos.size()),sz2(_time_steps.size());
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> > ret(sz);
  std::vector< std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> > > ts(sz2);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=shallowCpy();
      ret[i]->_infos.resize(1); ret[i]->_infos[0]=_infos[i];
    }
  for(std::size_t i=0;i<sz2;i++)
    {
      std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> > ret1=_time_steps[i]->splitComponents();
      if(ret1.size()!=sz)
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::splitComponents : At rank #" << i << " number of components is " << ret1.size() << " whereas it should be for all time steps " << sz << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      ts[i]=ret1;
    }
  for(std::size_t i=0;i<sz;i++)
    for(std::size_t j=0;j<sz2;j++)
      ret[i]->_time_steps[j]=ts[j][i];
  return ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::copyTinyInfoFrom(const MEDCouplingFieldDouble *field, const DataArray *arr) throw(INTERP_KERNEL::Exception)
{
  _name=field->getName();
  if(_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::copyTinyInfoFrom : unsupported fields with no name in MED file !");
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::copyTinyInfoFrom : no array set !");
  _infos=arr->getInfoOnComponents();
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::checkCoherencyOfTinyInfo(const MEDCouplingFieldDouble *field, const DataArray *arr) const throw(INTERP_KERNEL::Exception)
{
  static const char MSG[]="MEDFileFieldMultiTSWithoutSDA::checkCoherencyOfTinyInfo : invalid ";
  if(_name!=field->getName())
    {
      std::ostringstream oss; oss << MSG << "name ! should be \"" << _name;
      oss << "\" and it is set in input field to \"" << field->getName() << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::checkCoherencyOfTinyInfo : no array set !");
  checkThatComponentsMatch(arr->getInfoOnComponents());
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::checkThatComponentsMatch(const std::vector<std::string>& compos) const throw(INTERP_KERNEL::Exception)
{
  static const char MSG[]="MEDFileFieldMultiTSWithoutSDA::checkThatComponentsMatch : ";
  if(getInfo().size()!=compos.size())
    {
      std::ostringstream oss; oss << MSG << "mismatch of number of components between this (" << getInfo().size() << ") and ";
      oss << " number of components of element to append (" << compos.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(_infos!=compos)
    {
      std::ostringstream oss; oss << MSG << "components have same size but are different ! should be \"";
      std::copy(_infos.begin(),_infos.end(),std::ostream_iterator<std::string>(oss,", "));
      oss << " But compo in input fields are : ";
      std::copy(compos.begin(),compos.end(),std::ostream_iterator<std::string>(oss,", "));
      oss << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::checkThatNbOfCompoOfTSMatchThis() const throw(INTERP_KERNEL::Exception)
{
  std::size_t sz=_infos.size();
  int j=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,j++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *elt(*it);
      if(elt)
        if(elt->getInfo().size()!=sz)
          {
            std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::checkThatNbOfCompoOfTSMatchThis : At pos #" << j << " the number of components is equal to ";
            oss << elt->getInfo().size() << " whereas it is expected to be equal to " << sz << " !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::appendFieldNoProfileSBT : input field is NULL !");
  if(!_time_steps.empty())
    checkCoherencyOfTinyInfo(field,arr);
  MEDFileAnyTypeField1TSWithoutSDA *objC=createNew1TSWithoutSDAEmptyInstance();
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> obj(objC);
  objC->setFieldNoProfileSBT(field,arr,glob,*this);
  copyTinyInfoFrom(field,arr);
  _time_steps.push_back(obj);
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::appendFieldProfile(const MEDCouplingFieldDouble *field, const DataArray *arr, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTSWithoutSDA::appendFieldNoProfileSBT : input field is NULL !");
  if(!_time_steps.empty())
    checkCoherencyOfTinyInfo(field,arr);
  MEDFileField1TSWithoutSDA *objC=new MEDFileField1TSWithoutSDA;
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> obj(objC);
  objC->setFieldProfile(field,arr,mesh,meshDimRelToMax,profile,glob,*this);
  copyTinyInfoFrom(field,arr);
  _time_steps.push_back(obj);
}

//= MEDFileFieldMultiTSWithoutSDA

MEDFileFieldMultiTSWithoutSDA *MEDFileFieldMultiTSWithoutSDA::New(med_idt fid, const char *fieldName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTSWithoutSDA(fid,fieldName,fieldTyp,infos,nbOfStep,dtunit);
}

MEDFileFieldMultiTSWithoutSDA::MEDFileFieldMultiTSWithoutSDA()
{
}

MEDFileFieldMultiTSWithoutSDA::MEDFileFieldMultiTSWithoutSDA(const char *fieldName):MEDFileAnyTypeFieldMultiTSWithoutSDA(fieldName)
{
}

/*!
 * \param [in] fieldId field id in C mode
 */
MEDFileFieldMultiTSWithoutSDA::MEDFileFieldMultiTSWithoutSDA(med_idt fid, int fieldId) throw(INTERP_KERNEL::Exception)
try:MEDFileAnyTypeFieldMultiTSWithoutSDA(fid,fieldId)
{
}
catch(INTERP_KERNEL::Exception& e)
  { throw e; }

MEDFileFieldMultiTSWithoutSDA::MEDFileFieldMultiTSWithoutSDA(med_idt fid, const char *fieldName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit) throw(INTERP_KERNEL::Exception)
try:MEDFileAnyTypeFieldMultiTSWithoutSDA(fid,fieldName,fieldTyp,infos,nbOfStep,dtunit)
{
}
catch(INTERP_KERNEL::Exception& e)
{ throw e; }

MEDFileAnyTypeField1TSWithoutSDA *MEDFileFieldMultiTSWithoutSDA::createNew1TSWithoutSDAEmptyInstance() const throw(INTERP_KERNEL::Exception)
{
  return new MEDFileField1TSWithoutSDA;
}

void MEDFileFieldMultiTSWithoutSDA::checkCoherencyOfType(const MEDFileAnyTypeField1TSWithoutSDA *f1ts) const throw(INTERP_KERNEL::Exception)
{
  if(!f1ts)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::checkCoherencyOfType : input field1TS is NULL ! Impossible to check !");
  const MEDFileField1TSWithoutSDA *f1tsC=dynamic_cast<const MEDFileField1TSWithoutSDA *>(f1ts);
  if(!f1tsC)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::checkCoherencyOfType : the input field1TS is not a FLOAT64 type !");
}

const char *MEDFileFieldMultiTSWithoutSDA::getTypeStr() const throw(INTERP_KERNEL::Exception)
{
  return MEDFileField1TSWithoutSDA::TYPE_STR;
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileFieldMultiTSWithoutSDA::shallowCpy() const throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTSWithoutSDA(*this);
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileFieldMultiTSWithoutSDA::createNew() const throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTSWithoutSDA;
}

/*!
 * entry point for users that want to iterate into MEDFile DataStructure with a reduced overhead because output arrays are extracted (created) specially
 * for the call of this method. That's why the DataArrayDouble instance in returned vector of vector should be dealed by the caller.
 */
std::vector< std::vector<DataArrayDouble *> > MEDFileFieldMultiTSWithoutSDA::getFieldSplitedByType2(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA& myF1TS=getTimeStepEntry(iteration,order);
  const MEDFileField1TSWithoutSDA *myF1TSC=dynamic_cast<const MEDFileField1TSWithoutSDA *>(&myF1TS);
  if(!myF1TSC)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::getFieldSplitedByType2 : mismatch of type of field expecting FLOAT64 !");
  return myF1TSC->getFieldSplitedByType2(mname,types,typesF,pfls,locs);
}

MEDFileAnyTypeFieldMultiTS::MEDFileAnyTypeFieldMultiTS()
{
}

MEDFileAnyTypeFieldMultiTS::MEDFileAnyTypeFieldMultiTS(const char *fileName) throw(INTERP_KERNEL::Exception)
try:MEDFileFieldGlobsReal(fileName)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  _content=BuildContentFrom(fid,fileName);
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTS::BuildContentFrom(med_idt fid, const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  med_field_type typcha;
  std::vector<std::string> infos;
  std::string dtunit;
  int i=-1;
  MEDFileAnyTypeField1TS::LocateField(fid,fileName,fieldName,i,typcha,infos,dtunit);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> ret;
  switch(typcha)
    {
    case MED_FLOAT64:
      {
        ret=new MEDFileFieldMultiTSWithoutSDA(fid,i);
        break;
      }
    case MED_INT32:
      {
        ret=new MEDFileIntFieldMultiTSWithoutSDA(fid,i);
        break;
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::BuildContentFrom(fileName,fieldName) : file \'" << fileName << "\' contains field with name \'" << fieldName << "\' but the type of field is not in [MED_FLOAT64, MED_INT32] !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    }
  ret->setDtUnit(dtunit.c_str());
  return ret.retn();
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTS::BuildContentFrom(med_idt fid, const char *fileName) throw(INTERP_KERNEL::Exception)
{
  med_field_type typcha;
  //
  std::vector<std::string> infos;
  std::string dtunit,fieldName;
  MEDFileAnyTypeField1TS::LocateField2(fid,fileName,0,true,fieldName,typcha,infos,dtunit);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> ret;
  switch(typcha)
    {
    case MED_FLOAT64:
      {
        ret=new MEDFileFieldMultiTSWithoutSDA(fid,0);
        break;
      }
    case MED_INT32:
      {
        ret=new MEDFileIntFieldMultiTSWithoutSDA(fid,0);
        break;
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::BuildContentFrom(fileName) : file \'" << fileName << "\' contains field with name \'" << fieldName << "\' but the type of the first field is not in [MED_FLOAT64, MED_INT32] !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    }
  ret->setDtUnit(dtunit.c_str());
  return ret.retn();
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::BuildNewInstanceFromContent(MEDFileAnyTypeFieldMultiTSWithoutSDA *c, const char *fileName) throw(INTERP_KERNEL::Exception)
{
  if(!c)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::BuildNewInstanceFromContent : empty content in input : unable to build a new instance !");
  if(dynamic_cast<const MEDFileFieldMultiTSWithoutSDA *>(c))
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTS> ret=MEDFileFieldMultiTS::New();
      ret->setFileName(fileName);
      ret->_content=c;  c->incrRef();
      return ret.retn();
    }
  if(dynamic_cast<const MEDFileIntFieldMultiTSWithoutSDA *>(c))
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileIntFieldMultiTS> ret=MEDFileIntFieldMultiTS::New();
      ret->setFileName(fileName);
      ret->_content=c;  c->incrRef();
      return ret.retn();
    }
  throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::BuildNewInstanceFromContent : internal error ! a content of type different from FLOAT64 and INT32 has been built but not intercepted !");
}

MEDFileAnyTypeFieldMultiTS::MEDFileAnyTypeFieldMultiTS(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
try:MEDFileFieldGlobsReal(fileName)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  _content=BuildContentFrom(fid,fileName,fieldName);
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

//= MEDFileIntFieldMultiTSWithoutSDA

MEDFileIntFieldMultiTSWithoutSDA *MEDFileIntFieldMultiTSWithoutSDA::New(med_idt fid, const char *fieldName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileIntFieldMultiTSWithoutSDA(fid,fieldName,fieldTyp,infos,nbOfStep,dtunit);
}

MEDFileIntFieldMultiTSWithoutSDA::MEDFileIntFieldMultiTSWithoutSDA()
{
}

MEDFileIntFieldMultiTSWithoutSDA::MEDFileIntFieldMultiTSWithoutSDA(const char *fieldName):MEDFileAnyTypeFieldMultiTSWithoutSDA(fieldName)
{
}

MEDFileIntFieldMultiTSWithoutSDA::MEDFileIntFieldMultiTSWithoutSDA(med_idt fid, const char *fieldName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit) throw(INTERP_KERNEL::Exception)
try:MEDFileAnyTypeFieldMultiTSWithoutSDA(fid,fieldName,fieldTyp,infos,nbOfStep,dtunit)
{
}
catch(INTERP_KERNEL::Exception& e)
{ throw e; }

/*!
 * \param [in] fieldId field id in C mode
 */
MEDFileIntFieldMultiTSWithoutSDA::MEDFileIntFieldMultiTSWithoutSDA(med_idt fid, int fieldId) throw(INTERP_KERNEL::Exception)
try:MEDFileAnyTypeFieldMultiTSWithoutSDA(fid,fieldId)
{
}
catch(INTERP_KERNEL::Exception& e)
  { throw e; }

MEDFileAnyTypeField1TSWithoutSDA *MEDFileIntFieldMultiTSWithoutSDA::createNew1TSWithoutSDAEmptyInstance() const throw(INTERP_KERNEL::Exception)
{
  return new MEDFileIntField1TSWithoutSDA;
}

void MEDFileIntFieldMultiTSWithoutSDA::checkCoherencyOfType(const MEDFileAnyTypeField1TSWithoutSDA *f1ts) const throw(INTERP_KERNEL::Exception)
{
  if(!f1ts)
    throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTSWithoutSDA::checkCoherencyOfType : input field1TS is NULL ! Impossible to check !");
  const MEDFileIntField1TSWithoutSDA *f1tsC=dynamic_cast<const MEDFileIntField1TSWithoutSDA *>(f1ts);
  if(!f1tsC)
    throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTSWithoutSDA::checkCoherencyOfType : the input field1TS is not a INT32 type !");
}

const char *MEDFileIntFieldMultiTSWithoutSDA::getTypeStr() const throw(INTERP_KERNEL::Exception)
{
  return MEDFileIntField1TSWithoutSDA::TYPE_STR;
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileIntFieldMultiTSWithoutSDA::shallowCpy() const throw(INTERP_KERNEL::Exception)
{
  return new MEDFileIntFieldMultiTSWithoutSDA(*this);
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileIntFieldMultiTSWithoutSDA::createNew() const throw(INTERP_KERNEL::Exception)
{
  return new MEDFileIntFieldMultiTSWithoutSDA;
}

//= MEDFileAnyTypeFieldMultiTS

/*!
 * Returns a new instance of MEDFileFieldMultiTS or MEDFileIntFieldMultiTS holding data of the first field
 * that has been read from a specified MED file.
 *  \param [in] fileName - the name of the MED file to read.
 *  \return MEDFileFieldMultiTS * - a new instance of MEDFileFieldMultiTS or MEDFileIntFieldMultiTS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 *  \throw If reading the file fails.
 */
MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> c=BuildContentFrom(fid,fileName);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTS> ret=BuildNewInstanceFromContent(c,fileName);
  ret->loadGlobals(fid);
  return ret.retn();
}

/*!
 * Returns a new instance of MEDFileFieldMultiTS or MEDFileIntFieldMultiTS holding data of a given field
 * that has been read from a specified MED file.
 *  \param [in] fileName - the name of the MED file to read.
 *  \param [in] fieldName - the name of the field to read.
 *  \return MEDFileFieldMultiTS * - a new instance of MEDFileFieldMultiTS or MEDFileIntFieldMultiTS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 *  \throw If reading the file fails.
 *  \throw If there is no field named \a fieldName in the file.
 */
MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::New(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> c=BuildContentFrom(fid,fileName,fieldName);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTS> ret=BuildNewInstanceFromContent(c,fileName);
  ret->loadGlobals(fid);
  return ret.retn();
}

/*!
 * This constructor is a shallow copy constructor. If \a shallowCopyOfContent is true the content of \a other is shallow copied.
 * If \a shallowCopyOfContent is false, \a other is taken to be the content of \a this.
 *
 * \warning this is a shallow copy constructor
 */
MEDFileAnyTypeFieldMultiTS::MEDFileAnyTypeFieldMultiTS(const MEDFileAnyTypeFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent)
{
  if(!shallowCopyOfContent)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *otherPtr(&other);
      otherPtr->incrRef();
      _content=const_cast<MEDFileAnyTypeFieldMultiTSWithoutSDA *>(otherPtr);
    }
  else
    {
      _content=other.shallowCpy();
    }
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTS::contentNotNullBase() throw(INTERP_KERNEL::Exception)
{
  MEDFileAnyTypeFieldMultiTSWithoutSDA *ret=_content;
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS : content is expected to be not null !");
  return ret;
}

const MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTS::contentNotNullBase() const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *ret=_content;
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS : const content is expected to be not null !");
  return ret;
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTS::getPflsReallyUsed() const
{
  return contentNotNullBase()->getPflsReallyUsed2();
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTS::getLocsReallyUsed() const
{
  return contentNotNullBase()->getLocsReallyUsed2();
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTS::getPflsReallyUsedMulti() const
{
  return contentNotNullBase()->getPflsReallyUsedMulti2();
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTS::getLocsReallyUsedMulti() const
{
  return contentNotNullBase()->getLocsReallyUsedMulti2();
}

void MEDFileAnyTypeFieldMultiTS::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  contentNotNullBase()->changePflsRefsNamesGen2(mapOfModif);
}

void MEDFileAnyTypeFieldMultiTS::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  contentNotNullBase()->changeLocsRefsNamesGen2(mapOfModif);
}

int MEDFileAnyTypeFieldMultiTS::getNumberOfTS() const
{
  return contentNotNullBase()->getNumberOfTS();
}

void MEDFileAnyTypeFieldMultiTS::eraseEmptyTS() throw(INTERP_KERNEL::Exception)
{
  contentNotNullBase()->eraseEmptyTS();
}

void MEDFileAnyTypeFieldMultiTS::eraseTimeStepIds(const int *startIds, const int *endIds) throw(INTERP_KERNEL::Exception)
{
  contentNotNullBase()->eraseTimeStepIds(startIds,endIds);
}

void MEDFileAnyTypeFieldMultiTS::eraseTimeStepIds2(int bg, int end, int step) throw(INTERP_KERNEL::Exception)
{
  contentNotNullBase()->eraseTimeStepIds2(bg,end,step);
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::buildSubPart(const int *startIds, const int *endIds) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> c=contentNotNullBase()->buildFromTimeStepIds(startIds,endIds);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTS> ret=shallowCpy();
  ret->_content=c;
  return ret.retn();
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::buildSubPartSlice(int bg, int end, int step) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> c=contentNotNullBase()->buildFromTimeStepIds2(bg,end,step);
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTS> ret=shallowCpy();
  ret->_content=c;
  return ret.retn();
}

std::vector< std::pair<int,int> > MEDFileAnyTypeFieldMultiTS::getIterations() const
{
  return contentNotNullBase()->getIterations();
}

void MEDFileAnyTypeFieldMultiTS::pushBackTimeSteps(const std::vector<MEDFileAnyTypeField1TS *>& f1ts) throw(INTERP_KERNEL::Exception)
{
  for(std::vector<MEDFileAnyTypeField1TS *>::const_iterator it=f1ts.begin();it!=f1ts.end();it++)
    pushBackTimeStep(*it);
}

void MEDFileAnyTypeFieldMultiTS::pushBackTimeStep(MEDFileAnyTypeField1TS *f1ts) throw(INTERP_KERNEL::Exception)
{
  if(!f1ts)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::pushBackTimeStep : input pointer is NULL !");
  checkCoherencyOfType(f1ts);
  f1ts->incrRef();
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TS> f1tsSafe(f1ts);
  MEDFileAnyTypeField1TSWithoutSDA *c=f1ts->contentNotNullBase();
  c->incrRef();
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> cSafe(c);
  if(!((MEDFileAnyTypeFieldMultiTSWithoutSDA *)_content))
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::pushBackTimeStep : no content in this !");
  _content->pushBackTimeStep(cSafe);
  appendGlobs(*f1ts,1e-12);
}

void MEDFileAnyTypeFieldMultiTS::synchronizeNameScope() throw(INTERP_KERNEL::Exception)
{
  contentNotNullBase()->synchronizeNameScope();
}

int MEDFileAnyTypeFieldMultiTS::getPosOfTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getPosOfTimeStep(iteration,order);
}

int MEDFileAnyTypeFieldMultiTS::getPosGivenTime(double time, double eps) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getPosGivenTime(time,eps);
}

int MEDFileAnyTypeFieldMultiTS::getNonEmptyLevels(int iteration, int order, const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getNonEmptyLevels(iteration,order,mname,levs);
}

std::vector< std::vector<TypeOfField> > MEDFileAnyTypeFieldMultiTS::getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getTypesOfFieldAvailable();
}

std::vector< std::vector< std::pair<int,int> > > MEDFileAnyTypeFieldMultiTS::getFieldSplitedByType(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getFieldSplitedByType(iteration,order,mname,types,typesF,pfls,locs);
}

std::string MEDFileAnyTypeFieldMultiTS::getName() const
{
  return contentNotNullBase()->getName();
}

void MEDFileAnyTypeFieldMultiTS::setName(const char *name)
{
  contentNotNullBase()->setName(name);
}

std::string MEDFileAnyTypeFieldMultiTS::getDtUnit() const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getDtUnit();
}

void MEDFileAnyTypeFieldMultiTS::setDtUnit(const char *dtUnit) throw(INTERP_KERNEL::Exception)
{
  contentNotNullBase()->setDtUnit(dtUnit);
}

void MEDFileAnyTypeFieldMultiTS::simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const
{
  contentNotNullBase()->simpleRepr(bkOffset,oss,fmtsId);
}

std::vector< std::pair<int,int> > MEDFileAnyTypeFieldMultiTS::getTimeSteps(std::vector<double>& ret1) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getTimeSteps(ret1);
}

std::string MEDFileAnyTypeFieldMultiTS::getMeshName() const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getMeshName();
}

void MEDFileAnyTypeFieldMultiTS::setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception)
{
  contentNotNullBase()->setMeshName(newMeshName);
}

bool MEDFileAnyTypeFieldMultiTS::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->changeMeshNames(modifTab);
}

const std::vector<std::string>& MEDFileAnyTypeFieldMultiTS::getInfo() const throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->getInfo();
}

void MEDFileAnyTypeFieldMultiTS::setInfo(const std::vector<std::string>& info) throw(INTERP_KERNEL::Exception)
{
  return contentNotNullBase()->setInfo(info);
}

int MEDFileAnyTypeFieldMultiTS::getNumberOfComponents() const throw(INTERP_KERNEL::Exception)
{
  const std::vector<std::string> ret=getInfo();
  return (int)ret.size();
}

void MEDFileAnyTypeFieldMultiTS::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  writeGlobals(fid,*this);
  contentNotNullBase()->writeLL(fid,*this);
}

/*!
 * Writes \a this field into a MED file specified by its name.
 *  \param [in] fileName - the MED file name.
 *  \param [in] mode - the writing mode. For more on \a mode, see \ref AdvMEDLoaderBasics.
 * - 2 - erase; an existing file is removed.
 * - 1 - append; same data should not be present in an existing file.
 * - 0 - overwrite; same data present in an existing file is overwritten.
 *  \throw If the field name is not set.
 *  \throw If no field data is set.
 *  \throw If \a mode == 1 and the same data is present in an existing file.
 */
void MEDFileAnyTypeFieldMultiTS::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,medmod);
  writeLL(fid);
}

std::string MEDFileAnyTypeFieldMultiTS::simpleRepr() const
{
  std::ostringstream oss;
  contentNotNullBase()->simpleRepr(0,oss,-1);
  simpleReprGlobs(oss);
  return oss.str();
}

std::size_t MEDFileAnyTypeFieldMultiTS::getHeapMemorySize() const
{
  std::size_t ret=0;
  if((const MEDFileAnyTypeFieldMultiTSWithoutSDA*)_content)
    ret+=_content->getHeapMemorySize();
  return ret+MEDFileFieldGlobsReal::getHeapMemorySize();
}

/*!
 * This method returns as MEDFileAnyTypeField1TS new instances as number of components in \a this.
 * The returned instances are deep copy of \a this except that for globals that are share with those contained in \a this.
 * ** WARNING ** do no forget to rename the ouput instances to avoid to write n-times in the same MED file field !
 */
std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTS > > MEDFileAnyTypeFieldMultiTS::splitComponents() const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *content(_content);
  if(!content)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::splitComponents : no content in this ! Unable to split components !");
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> > contentsSplit=content->splitComponents();
  std::size_t sz(contentsSplit.size());
  std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTS > > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=shallowCpy();
      ret[i]->_content=contentsSplit[i];
    }
  return ret;
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::deepCpy() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTS> ret=shallowCpy();
  if((const MEDFileAnyTypeFieldMultiTSWithoutSDA *)_content)
    ret->_content=_content->deepCpy();
  ret->deepCpyGlobs(*this);
  return ret.retn();
}

MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> MEDFileAnyTypeFieldMultiTS::getContent()
{
  return _content;
}

/*!
 * Returns a new MEDFileField1TS or MEDFileIntField1TS holding data of a given time step of \a this field.
 *  \param [in] iteration - the iteration number of a required time step.
 *  \param [in] order - the iteration order number of required time step.
 *  \return MEDFileField1TS * or MEDFileIntField1TS *- a new instance of MEDFileField1TS or MEDFileIntField1TS. The caller is to
 *          delete this field using decrRef() as it is no more needed.
 *  \throw If there is no required time step in \a this field.
 */
MEDFileAnyTypeField1TS *MEDFileAnyTypeFieldMultiTS::getTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  int pos=getPosOfTimeStep(iteration,order);
  return getTimeStepAtPos(pos);
}

/*!
 * Returns a new MEDFileField1TS or MEDFileIntField1TS holding data of a given time step of \a this field.
 *  \param [in] time - the time of the time step of interest.
 *  \param [in] eps - a precision used to compare time values.
 *  \return MEDFileField1TS * - a new instance of MEDFileField1TS. The caller is to
 *          delete this field using decrRef() as it is no more needed.
 *  \throw If there is no required time step in \a this field.
 */
MEDFileAnyTypeField1TS *MEDFileAnyTypeFieldMultiTS::getTimeStepGivenTime(double time, double eps) const throw(INTERP_KERNEL::Exception)
{
  int pos=getPosGivenTime(time,eps);
  return getTimeStepAtPos(pos);
}

MEDFileAnyTypeFieldMultiTSIterator *MEDFileAnyTypeFieldMultiTS::iterator() throw(INTERP_KERNEL::Exception)
{
  return new MEDFileAnyTypeFieldMultiTSIterator(this);
}

//= MEDFileFieldMultiTS

/*!
 * Returns a new empty instance of MEDFileFieldMultiTS.
 *  \return MEDFileFieldMultiTS * - a new instance of MEDFileFieldMultiTS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 */
MEDFileFieldMultiTS *MEDFileFieldMultiTS::New()
{
  return new MEDFileFieldMultiTS;
}

/*!
 * Returns a new instance of MEDFileFieldMultiTS holding data of the first field
 * that has been read from a specified MED file.
 *  \param [in] fileName - the name of the MED file to read.
 *  \return MEDFileFieldMultiTS * - a new instance of MEDFileFieldMultiTS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 *  \throw If reading the file fails.
 */
MEDFileFieldMultiTS *MEDFileFieldMultiTS::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTS> ret=new MEDFileFieldMultiTS(fileName);
  ret->contentNotNull();//to check that content type matches with \a this type.
  return ret.retn();
}

/*!
 * Returns a new instance of MEDFileFieldMultiTS holding data of a given field
 * that has been read from a specified MED file.
 *  \param [in] fileName - the name of the MED file to read.
 *  \param [in] fieldName - the name of the field to read.
 *  \return MEDFileFieldMultiTS * - a new instance of MEDFileFieldMultiTS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 *  \throw If reading the file fails.
 *  \throw If there is no field named \a fieldName in the file.
 */
MEDFileFieldMultiTS *MEDFileFieldMultiTS::New(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTS> ret=new MEDFileFieldMultiTS(fileName,fieldName);
  ret->contentNotNull();//to check that content type matches with \a this type.
  return ret.retn();
}

/*!
 * Returns a new instance of MEDFileFieldMultiTS. If \a shallowCopyOfContent is true the content of \a other is shallow copied.
 * If \a shallowCopyOfContent is false, \a other is taken to be the content of \a this.
 *
 * Returns a new instance of MEDFileFieldMultiTS holding either a shallow copy
 * of a given MEDFileFieldMultiTSWithoutSDA ( \a other ) or \a other itself.
 * \warning this is a shallow copy constructor
 *  \param [in] other - a MEDFileField1TSWithoutSDA to copy.
 *  \param [in] shallowCopyOfContent - if \c true, a shallow copy of \a other is created.
 *  \return MEDFileFieldMultiTS * - a new instance of MEDFileFieldMultiTS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 */
MEDFileFieldMultiTS *MEDFileFieldMultiTS::New(const MEDFileFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent)
{
  return new MEDFileFieldMultiTS(other,shallowCopyOfContent);
}

MEDFileAnyTypeFieldMultiTS *MEDFileFieldMultiTS::shallowCpy() const throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTS(*this);
}

void MEDFileFieldMultiTS::checkCoherencyOfType(const MEDFileAnyTypeField1TS *f1ts) const throw(INTERP_KERNEL::Exception)
{
  if(!f1ts)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::checkCoherencyOfType : input field1TS is NULL ! Impossible to check !");
  const MEDFileField1TS *f1tsC=dynamic_cast<const MEDFileField1TS *>(f1ts);
  if(!f1tsC)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::checkCoherencyOfType : the input field1TS is not a FLOAT64 type !");
}

/*!
 * Returns a new MEDFileField1TS holding data of a given time step of \a this field.
 *  \param [in] pos - a time step id.
 *  \return MEDFileField1TS * - a new instance of MEDFileField1TS. The caller is to
 *          delete this field using decrRef() as it is no more needed.
 *  \throw If \a pos is not a valid time step id.
 */
MEDFileAnyTypeField1TS *MEDFileFieldMultiTS::getTimeStepAtPos(int pos) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA *item=contentNotNullBase()->getTimeStepAtPos2(pos);
  if(!item)
    {
      std::ostringstream oss; oss << "MEDFileFieldMultiTS::getTimeStepAtPos : field at pos #" << pos << " is null !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const MEDFileField1TSWithoutSDA *itemC=dynamic_cast<const MEDFileField1TSWithoutSDA *>(item);
  if(itemC)
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> ret=MEDFileField1TS::New(*itemC,false);
      ret->shallowCpyGlobs(*this);
      return ret.retn();
    }
  std::ostringstream oss; oss << "MEDFileFieldMultiTS::getTimeStepAtPos : type of field at pos #" << pos << " is not FLOAT64 !";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

/*!
 * Returns a new MEDCouplingFieldDouble of a given type, of a given time step, lying on
 * mesh entities of a given dimension of the first mesh in MED file.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of interest.
 *  \param [in] iteration - the iteration number of a required time step.
 *  \param [in] order - the iteration order number of required time step.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If the MED file is not readable.
 *  \throw If there is no mesh in the MED file.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in the mesh.
 *  \throw If no field values of the required parameters are available.
 */
MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA& myF1TS=contentNotNullBase()->getTimeStepEntry(iteration,order);
  const MEDFileField1TSWithoutSDA *myF1TSC=dynamic_cast<const MEDFileField1TSWithoutSDA *>(&myF1TS);
  if(!myF1TSC)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::getFieldAtLevel : mismatch of type of field expecting FLOAT64 !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arrOut;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=myF1TSC->getFieldAtLevel(type,meshDimRelToMax,0,renumPol,this,arrOut,*contentNotNullBase());
  MEDFileField1TS::SetDataArrayDoubleInField(ret,arrOut);
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble of a given type, of a given time step, lying on
 * the top level cells of the first mesh in MED file.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of interest.
 *  \param [in] iteration - the iteration number of a required time step.
 *  \param [in] order - the iteration order number of required time step.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If the MED file is not readable.
 *  \throw If there is no mesh in the MED file.
 *  \throw If no field values of the required parameters are available.
 */
MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldAtTopLevel(TypeOfField type, int iteration, int order, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA& myF1TS=contentNotNullBase()->getTimeStepEntry(iteration,order);
  const MEDFileField1TSWithoutSDA *myF1TSC=dynamic_cast<const MEDFileField1TSWithoutSDA *>(&myF1TS);
  if(!myF1TSC)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::getFieldAtTopLevel : mismatch of type of field !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arrOut;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=myF1TSC->getFieldAtTopLevel(type,0,renumPol,this,arrOut,*contentNotNullBase());
  MEDFileField1TS::SetDataArrayDoubleInField(ret,arrOut);
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble of a given type, of a given time step, lying on
 * a given support.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of interest.
 *  \param [in] iteration - the iteration number of a required time step.
 *  \param [in] order - the iteration order number of required time step.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [in] mesh - the supporting mesh.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in the mesh.
 *  \throw If no field of \a this is lying on \a mesh.
 *  \throw If no field values of the required parameters are available.
 */
MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA& myF1TS=contentNotNullBase()->getTimeStepEntry(iteration,order);
  const MEDFileField1TSWithoutSDA *myF1TSC=dynamic_cast<const MEDFileField1TSWithoutSDA *>(&myF1TS);
  if(!myF1TSC)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::getFieldOnMeshAtLevel : mismatch of type of field !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arrOut;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=myF1TSC->getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,this,mesh,arrOut,*contentNotNullBase());
  MEDFileField1TS::SetDataArrayDoubleInField(ret,arrOut);
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble of given type, of a given time step, lying on a
 * given support. 
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of the new field.
 *  \param [in] iteration - the iteration number of a required time step.
 *  \param [in] order - the iteration order number of required time step.
 *  \param [in] mesh - the supporting mesh.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If no field of \a this is lying on \a mesh.
 *  \throw If no field values of the required parameters are available.
 */
MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA& myF1TS=contentNotNullBase()->getTimeStepEntry(iteration,order);
  const MEDFileField1TSWithoutSDA *myF1TSC=dynamic_cast<const MEDFileField1TSWithoutSDA *>(&myF1TS);
  if(!myF1TSC)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::getFieldOnMeshAtLevel : mismatch of type of field !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arrOut;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=myF1TSC->getFieldOnMeshAtLevel(type,renumPol,this,mesh,0,0,arrOut,*contentNotNullBase());
  MEDFileField1TS::SetDataArrayDoubleInField(ret,arrOut);
  return ret.retn();
}

/*!
 * This method has a close behaviour than MEDFileFieldMultiTS::getFieldAtLevel.
 * This method is called 'old' because the user should give the mesh name he wants to use for it's field.
 * This method is useful for MED2 file format when field on different mesh was autorized.
 */
MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldAtLevelOld(TypeOfField type, const char *mname, int iteration, int order, int meshDimRelToMax, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA& myF1TS=contentNotNullBase()->getTimeStepEntry(iteration,order);
  const MEDFileField1TSWithoutSDA *myF1TSC=dynamic_cast<const MEDFileField1TSWithoutSDA *>(&myF1TS);
  if(!myF1TSC)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::getFieldAtLevelOld : mismatch of type of field !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arrOut;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=myF1TSC->getFieldAtLevel(type,meshDimRelToMax,mname,renumPol,this,arrOut,*contentNotNullBase());
  MEDFileField1TS::SetDataArrayDoubleInField(ret,arrOut);
  return ret.retn();
}

/*!
 * Returns values and a profile of the field of a given type, of a given time step,
 * lying on a given support.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of the field.
 *  \param [in] iteration - the iteration number of a required time step.
 *  \param [in] order - the iteration order number of required time step.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [in] mesh - the supporting mesh.
 *  \param [out] pfl - a new instance of DataArrayInt holding ids of mesh entities the
 *          field of interest lies on. If the field lies on all entities of the given
 *          dimension, all ids in \a pfl are zero. The caller is to delete this array
 *          using decrRef() as it is no more needed.  
 *  \param [in] glob - the global data storing profiles and localization.
 *  \return DataArrayDouble * - a new instance of DataArrayDouble holding values of the
 *          field. The caller is to delete this array using decrRef() as it is no more needed.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in \a mesh.
 *  \throw If no field of \a this is lying on \a mesh.
 *  \throw If no field values of the required parameters are available.
 */
DataArrayDouble *MEDFileFieldMultiTS::getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA& myF1TS=contentNotNullBase()->getTimeStepEntry(iteration,order);
  const MEDFileField1TSWithoutSDA *myF1TSC=dynamic_cast<const MEDFileField1TSWithoutSDA *>(&myF1TS);
  if(!myF1TSC)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::getFieldWithProfile : mismatch of type of field !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> ret=myF1TSC->getFieldWithProfile(type,meshDimRelToMax,mesh,pfl,this,*contentNotNullBase());
  return MEDFileField1TS::ReturnSafelyDataArrayDouble(ret);
}

const MEDFileFieldMultiTSWithoutSDA *MEDFileFieldMultiTS::contentNotNull() const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *pt(_content);
  if(!pt)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::contentNotNull : the content pointer is null !");
  const MEDFileFieldMultiTSWithoutSDA *ret=dynamic_cast<const MEDFileFieldMultiTSWithoutSDA *>(pt);
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::contentNotNull : the content pointer is not null but it is not of type double ! Reason is maybe that the read field has not the type FLOAT64 !");
  return ret;
}

 MEDFileFieldMultiTSWithoutSDA *MEDFileFieldMultiTS::contentNotNull() throw(INTERP_KERNEL::Exception)
{
  MEDFileAnyTypeFieldMultiTSWithoutSDA *pt(_content);
  if(!pt)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::contentNotNull : the non const content pointer is null !");
  MEDFileFieldMultiTSWithoutSDA *ret=dynamic_cast<MEDFileFieldMultiTSWithoutSDA *>(pt);
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::contentNotNull : the non const content pointer is not null but it is not of type double ! Reason is maybe that the read field has not the type FLOAT64 !");
  return ret;
}

/*!
 * Adds a MEDCouplingFieldDouble to \a this as another time step. The underlying mesh of
 * the given field is checked if its elements are sorted suitable for writing to MED file
 * ("STB" stands for "Sort By Type"), if not, an exception is thrown. 
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] field - the field to add to \a this.
 *  \throw If the name of \a field is empty.
 *  \throw If the data array of \a field is not set.
 *  \throw If existing time steps have different name or number of components than \a field.
 *  \throw If the underlying mesh of \a field has no name.
 *  \throw If elements in the mesh are not in the order suitable for writing to the MED file.
 */
void MEDFileFieldMultiTS::appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception)
{
  const DataArrayDouble *arr=0;
  if(field)
    arr=field->getArray();
  contentNotNull()->appendFieldNoProfileSBT(field,arr,*this);
}

/*!
 * Adds a MEDCouplingFieldDouble to \a this as another time step. Specified entities of
 * a given dimension of a given mesh are used as the support of the given field.
 * Elements of the given mesh must be sorted suitable for writing to MED file. 
 * Order of underlying mesh entities of the given field specified by \a profile parameter
 * is not prescribed; this method permutes field values to have them sorted by element
 * type as required for writing to MED file.  
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] field - the field to add to \a this.
 *  \param [in] mesh - the supporting mesh of \a field.
 *  \param [in] meshDimRelToMax - a relative dimension of mesh entities \a field lies on.
 *  \param [in] profile - ids of mesh entities on which corresponding field values lie.
 *  \throw If either \a field or \a mesh or \a profile has an empty name.
 *  \throw If existing time steps have different name or number of components than \a field.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in \a mesh.
 *  \throw If the data array of \a field is not set.
 *  \throw If elements in \a mesh are not in the order suitable for writing to the MED file.
 */
void MEDFileFieldMultiTS::appendFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception)
{
  const DataArrayDouble *arr=0;
  if(field)
    arr=field->getArray();
  contentNotNull()->appendFieldProfile(field,arr,mesh,meshDimRelToMax,profile,*this);
}

MEDFileFieldMultiTS::MEDFileFieldMultiTS()
{
  _content=new MEDFileFieldMultiTSWithoutSDA;
}

MEDFileFieldMultiTS::MEDFileFieldMultiTS(const char *fileName) throw(INTERP_KERNEL::Exception)
try:MEDFileAnyTypeFieldMultiTS(fileName)
{
}
catch(INTERP_KERNEL::Exception& e)
  { throw e; }

MEDFileFieldMultiTS::MEDFileFieldMultiTS(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
try:MEDFileAnyTypeFieldMultiTS(fileName,fieldName)
{
}
catch(INTERP_KERNEL::Exception& e)
  { throw e; }

MEDFileFieldMultiTS::MEDFileFieldMultiTS(const MEDFileFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent):MEDFileAnyTypeFieldMultiTS(other,shallowCopyOfContent)
{
}

std::vector< std::vector<DataArrayDouble *> > MEDFileFieldMultiTS::getFieldSplitedByType2(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  return contentNotNull()->getFieldSplitedByType2(iteration,order,mname,types,typesF,pfls,locs);
}

DataArrayDouble *MEDFileFieldMultiTS::getUndergroundDataArray(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  return static_cast<DataArrayDouble *>(contentNotNull()->getUndergroundDataArray(iteration,order));
}

DataArrayDouble *MEDFileFieldMultiTS::getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception)
{
  return static_cast<DataArrayDouble *>(contentNotNull()->getUndergroundDataArrayExt(iteration,order,entries));
}

//= MEDFileAnyTypeFieldMultiTSIterator

MEDFileAnyTypeFieldMultiTSIterator::MEDFileAnyTypeFieldMultiTSIterator(MEDFileAnyTypeFieldMultiTS *fmts):_fmts(fmts),_iter_id(0),_nb_iter(0)
{
  if(fmts)
    {
      fmts->incrRef();
      _nb_iter=fmts->getNumberOfTS();
    }
}

MEDFileAnyTypeFieldMultiTSIterator::~MEDFileAnyTypeFieldMultiTSIterator() 
{
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeFieldMultiTSIterator::nextt() throw(INTERP_KERNEL::Exception)
{
  if(_iter_id<_nb_iter)
    {
      MEDFileAnyTypeFieldMultiTS *fmts(_fmts);
      if(fmts)
        return fmts->getTimeStepAtPos(_iter_id++);
      else
        return 0;
    }
  else
    return 0;
}

//= MEDFileIntFieldMultiTS

/*!
 * Returns a new empty instance of MEDFileFieldMultiTS.
 *  \return MEDFileIntFieldMultiTS * - a new instance of MEDFileIntFieldMultiTS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 */
MEDFileIntFieldMultiTS *MEDFileIntFieldMultiTS::New()
{
  return new MEDFileIntFieldMultiTS;
}

/*!
 * Returns a new instance of MEDFileIntFieldMultiTS holding data of the first field
 * that has been read from a specified MED file.
 *  \param [in] fileName - the name of the MED file to read.
 *  \return MEDFileFieldMultiTS * - a new instance of MEDFileIntFieldMultiTS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 *  \throw If reading the file fails.
 */
MEDFileIntFieldMultiTS *MEDFileIntFieldMultiTS::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileIntFieldMultiTS> ret=new MEDFileIntFieldMultiTS(fileName);
  ret->contentNotNull();//to check that content type matches with \a this type.
  return ret.retn();
}

/*!
 * Returns a new instance of MEDFileIntFieldMultiTS holding data of a given field
 * that has been read from a specified MED file.
 *  \param [in] fileName - the name of the MED file to read.
 *  \param [in] fieldName - the name of the field to read.
 *  \return MEDFileFieldMultiTS * - a new instance of MEDFileIntFieldMultiTS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 *  \throw If reading the file fails.
 *  \throw If there is no field named \a fieldName in the file.
 */
MEDFileIntFieldMultiTS *MEDFileIntFieldMultiTS::New(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileIntFieldMultiTS> ret=new MEDFileIntFieldMultiTS(fileName,fieldName);
  ret->contentNotNull();//to check that content type matches with \a this type.
  return ret.retn();
}

/*!
 * Returns a new instance of MEDFileIntFieldMultiTS. If \a shallowCopyOfContent is true the content of \a other is shallow copied.
 * If \a shallowCopyOfContent is false, \a other is taken to be the content of \a this.
 *
 * Returns a new instance of MEDFileIntFieldMultiTS holding either a shallow copy
 * of a given MEDFileIntFieldMultiTSWithoutSDA ( \a other ) or \a other itself.
 * \warning this is a shallow copy constructor
 *  \param [in] other - a MEDFileIntField1TSWithoutSDA to copy.
 *  \param [in] shallowCopyOfContent - if \c true, a shallow copy of \a other is created.
 *  \return MEDFileIntFieldMultiTS * - a new instance of MEDFileIntFieldMultiTS. The caller
 *          is to delete this field using decrRef() as it is no more needed.
 */
MEDFileIntFieldMultiTS *MEDFileIntFieldMultiTS::New(const MEDFileIntFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent)
{
  return new MEDFileIntFieldMultiTS(other,shallowCopyOfContent);
}

MEDFileAnyTypeFieldMultiTS *MEDFileIntFieldMultiTS::shallowCpy() const throw(INTERP_KERNEL::Exception)
{
  return new MEDFileIntFieldMultiTS(*this);
}

void MEDFileIntFieldMultiTS::checkCoherencyOfType(const MEDFileAnyTypeField1TS *f1ts) const throw(INTERP_KERNEL::Exception)
{
  if(!f1ts)
    throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTS::checkCoherencyOfType : input field1TS is NULL ! Impossible to check !");
  const MEDFileIntField1TS *f1tsC=dynamic_cast<const MEDFileIntField1TS *>(f1ts);
  if(!f1tsC)
    throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTS::checkCoherencyOfType : the input field1TS is not a INT32 type !");
}

/*!
 * Returns a new MEDCouplingFieldDouble of a given type, of a given time step, lying on
 * mesh entities of a given dimension of the first mesh in MED file.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of interest.
 *  \param [in] iteration - the iteration number of a required time step.
 *  \param [in] order - the iteration order number of required time step.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [out] arrOut - the DataArrayInt containing values of field.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If the MED file is not readable.
 *  \throw If there is no mesh in the MED file.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in the mesh.
 *  \throw If no field values of the required parameters are available.
 */
MEDCouplingFieldDouble *MEDFileIntFieldMultiTS::getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, DataArrayInt* &arrOut, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA& myF1TS=contentNotNullBase()->getTimeStepEntry(iteration,order);
  const MEDFileIntField1TSWithoutSDA *myF1TSC=dynamic_cast<const MEDFileIntField1TSWithoutSDA *>(&myF1TS);
  if(!myF1TSC)
    throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTS::getFieldAtLevel : mismatch of type of field expecting INT32 !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arr;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=myF1TSC->getFieldAtLevel(type,meshDimRelToMax,0,renumPol,this,arr,*contentNotNullBase());
  arrOut=MEDFileIntField1TS::ReturnSafelyDataArrayInt(arr);
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble of a given type, of a given time step, lying on
 * the top level cells of the first mesh in MED file.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of interest.
 *  \param [in] iteration - the iteration number of a required time step.
 *  \param [in] order - the iteration order number of required time step.
 *  \param [out] arrOut - the DataArrayInt containing values of field.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If the MED file is not readable.
 *  \throw If there is no mesh in the MED file.
 *  \throw If no field values of the required parameters are available.
 */
MEDCouplingFieldDouble *MEDFileIntFieldMultiTS::getFieldAtTopLevel(TypeOfField type, int iteration, int order, DataArrayInt* &arrOut, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA& myF1TS=contentNotNullBase()->getTimeStepEntry(iteration,order);
  const MEDFileIntField1TSWithoutSDA *myF1TSC=dynamic_cast<const MEDFileIntField1TSWithoutSDA *>(&myF1TS);
  if(!myF1TSC)
    throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTS::getFieldAtTopLevel : mismatch of type of field ! INT32 expected !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arr;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=myF1TSC->getFieldAtTopLevel(type,0,renumPol,this,arr,*contentNotNullBase());
  arrOut=MEDFileIntField1TS::ReturnSafelyDataArrayInt(arr);
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble of a given type, of a given time step, lying on
 * a given support.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of interest.
 *  \param [in] iteration - the iteration number of a required time step.
 *  \param [in] order - the iteration order number of required time step.
 *  \param [out] arrOut - the DataArrayInt containing values of field.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [in] mesh - the supporting mesh.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in the mesh.
 *  \throw If no field of \a this is lying on \a mesh.
 *  \throw If no field values of the required parameters are available.
 */
MEDCouplingFieldDouble *MEDFileIntFieldMultiTS::getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt* &arrOut, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA& myF1TS=contentNotNullBase()->getTimeStepEntry(iteration,order);
  const MEDFileIntField1TSWithoutSDA *myF1TSC=dynamic_cast<const MEDFileIntField1TSWithoutSDA *>(&myF1TS);
  if(!myF1TSC)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::getFieldOnMeshAtLevel : mismatch of type of field ! INT32 expected !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arr;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=myF1TSC->getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,this,mesh,arr,*contentNotNullBase());
  arrOut=MEDFileIntField1TS::ReturnSafelyDataArrayInt(arr);
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble of given type, of a given time step, lying on a
 * given support. 
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of the new field.
 *  \param [in] iteration - the iteration number of a required time step.
 *  \param [in] order - the iteration order number of required time step.
 *  \param [in] mesh - the supporting mesh.
 *  \param [out] arrOut - the DataArrayInt containing values of field.
 *  \param [in] renumPol - specifies how to permute values of the result field according to
 *          the optional numbers of cells and nodes, if any. The valid values are
 *          - 0 - do not permute.
 *          - 1 - permute cells.
 *          - 2 - permute nodes.
 *          - 3 - permute cells and nodes.
 *
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *          caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If no field of \a this is lying on \a mesh.
 *  \throw If no field values of the required parameters are available.
 */
MEDCouplingFieldDouble *MEDFileIntFieldMultiTS::getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, DataArrayInt* &arrOut, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA& myF1TS=contentNotNullBase()->getTimeStepEntry(iteration,order);
  const MEDFileIntField1TSWithoutSDA *myF1TSC=dynamic_cast<const MEDFileIntField1TSWithoutSDA *>(&myF1TS);
  if(!myF1TSC)
    throw INTERP_KERNEL::Exception("MEDFileFieldIntMultiTS::getFieldOnMeshAtLevel : mismatch of type of field ! INT32 expected !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arr;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=myF1TSC->getFieldOnMeshAtLevel(type,renumPol,this,mesh,0,0,arr,*contentNotNullBase());
  arrOut=MEDFileIntField1TS::ReturnSafelyDataArrayInt(arr);
  return ret.retn();
}

/*!
 * This method has a close behaviour than MEDFileIntFieldMultiTS::getFieldAtLevel.
 * This method is called 'old' because the user should give the mesh name he wants to use for it's field.
 * This method is useful for MED2 file format when field on different mesh was autorized.
 */
MEDCouplingFieldDouble *MEDFileIntFieldMultiTS::getFieldAtLevelOld(TypeOfField type, int iteration, int order, const char *mname, int meshDimRelToMax, DataArrayInt* &arrOut, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA& myF1TS=contentNotNullBase()->getTimeStepEntry(iteration,order);
  const MEDFileIntField1TSWithoutSDA *myF1TSC=dynamic_cast<const MEDFileIntField1TSWithoutSDA *>(&myF1TS);
  if(!myF1TSC)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::getFieldOnMeshAtLevel : mismatch of type of field ! INT32 expected !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> arr;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=myF1TSC->getFieldAtLevel(type,meshDimRelToMax,mname,renumPol,this,arr,*contentNotNullBase());
  arrOut=MEDFileIntField1TS::ReturnSafelyDataArrayInt(arr);
  return ret.retn();
}

/*!
 * Returns values and a profile of the field of a given type, of a given time step,
 * lying on a given support.
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] type - a spatial discretization of the field.
 *  \param [in] iteration - the iteration number of a required time step.
 *  \param [in] order - the iteration order number of required time step.
 *  \param [in] meshDimRelToMax - a relative dimension of the supporting mesh entities.
 *  \param [in] mesh - the supporting mesh.
 *  \param [out] pfl - a new instance of DataArrayInt holding ids of mesh entities the
 *          field of interest lies on. If the field lies on all entities of the given
 *          dimension, all ids in \a pfl are zero. The caller is to delete this array
 *          using decrRef() as it is no more needed.  
 *  \param [in] glob - the global data storing profiles and localization.
 *  \return DataArrayInt * - a new instance of DataArrayInt holding values of the
 *          field. The caller is to delete this array using decrRef() as it is no more needed.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in \a mesh.
 *  \throw If no field of \a this is lying on \a mesh.
 *  \throw If no field values of the required parameters are available.
 */
DataArrayInt *MEDFileIntFieldMultiTS::getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA& myF1TS=contentNotNullBase()->getTimeStepEntry(iteration,order);
  const MEDFileIntField1TSWithoutSDA *myF1TSC=dynamic_cast<const MEDFileIntField1TSWithoutSDA *>(&myF1TS);
  if(!myF1TSC)
    throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTS::getFieldWithProfile : mismatch of type of field ! INT32 expected !");
  MEDCouplingAutoRefCountObjectPtr<DataArray> ret=myF1TSC->getFieldWithProfile(type,meshDimRelToMax,mesh,pfl,this,*contentNotNullBase());
  return MEDFileIntField1TS::ReturnSafelyDataArrayInt(ret);
}

/*!
 * Returns a new MEDFileIntField1TS holding data of a given time step of \a this field.
 *  \param [in] pos - a time step id.
 *  \return MEDFileIntField1TS * - a new instance of MEDFileIntField1TS. The caller is to
 *          delete this field using decrRef() as it is no more needed.
 *  \throw If \a pos is not a valid time step id.
 */
MEDFileAnyTypeField1TS *MEDFileIntFieldMultiTS::getTimeStepAtPos(int pos) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeField1TSWithoutSDA *item=contentNotNullBase()->getTimeStepAtPos2(pos);
  if(!item)
    {
      std::ostringstream oss; oss << "MEDFileIntFieldMultiTS::getTimeStepAtPos : field at pos #" << pos << " is null !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const MEDFileIntField1TSWithoutSDA *itemC=dynamic_cast<const MEDFileIntField1TSWithoutSDA *>(item);
  if(itemC)
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileIntField1TS> ret=MEDFileIntField1TS::New(*itemC,false);
      ret->shallowCpyGlobs(*this);
      return ret.retn();
    }
  std::ostringstream oss; oss << "MEDFileIntFieldMultiTS::getTimeStepAtPos : type of field at pos #" << pos << " is not INT32 !";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

/*!
 * Adds a MEDCouplingFieldDouble to \a this as another time step. The underlying mesh of
 * the given field is checked if its elements are sorted suitable for writing to MED file
 * ("STB" stands for "Sort By Type"), if not, an exception is thrown. 
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] field - the field to add to \a this.
 *  \throw If the name of \a field is empty.
 *  \throw If the data array of \a field is not set.
 *  \throw If existing time steps have different name or number of components than \a field.
 *  \throw If the underlying mesh of \a field has no name.
 *  \throw If elements in the mesh are not in the order suitable for writing to the MED file.
 */
void MEDFileIntFieldMultiTS::appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals) throw(INTERP_KERNEL::Exception)
{
  contentNotNull()->appendFieldNoProfileSBT(field,arrOfVals,*this);
}

/*!
 * Adds a MEDCouplingFieldDouble to \a this as another time step. Specified entities of
 * a given dimension of a given mesh are used as the support of the given field.
 * Elements of the given mesh must be sorted suitable for writing to MED file. 
 * Order of underlying mesh entities of the given field specified by \a profile parameter
 * is not prescribed; this method permutes field values to have them sorted by element
 * type as required for writing to MED file.  
 * For more info, see \ref AdvMEDLoaderAPIFieldRW
 *  \param [in] field - the field to add to \a this.
 *  \param [in] mesh - the supporting mesh of \a field.
 *  \param [in] meshDimRelToMax - a relative dimension of mesh entities \a field lies on.
 *  \param [in] profile - ids of mesh entities on which corresponding field values lie.
 *  \throw If either \a field or \a mesh or \a profile has an empty name.
 *  \throw If existing time steps have different name or number of components than \a field.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in \a mesh.
 *  \throw If the data array of \a field is not set.
 *  \throw If elements in \a mesh are not in the order suitable for writing to the MED file.
 */
void MEDFileIntFieldMultiTS::appendFieldProfile(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception)
{
  contentNotNull()->appendFieldProfile(field,arrOfVals,mesh,meshDimRelToMax,profile,*this);
}

const MEDFileIntFieldMultiTSWithoutSDA *MEDFileIntFieldMultiTS::contentNotNull() const throw(INTERP_KERNEL::Exception)
{
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *pt(_content);
  if(!pt)
    throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTS::contentNotNull : the content pointer is null !");
  const MEDFileIntFieldMultiTSWithoutSDA *ret=dynamic_cast<const MEDFileIntFieldMultiTSWithoutSDA *>(pt);
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTS::contentNotNull : the content pointer is not null but it is not of type int ! Reason is maybe that the read field has not the type INT32 !");
  return ret;
}

 MEDFileIntFieldMultiTSWithoutSDA *MEDFileIntFieldMultiTS::contentNotNull() throw(INTERP_KERNEL::Exception)
{
  MEDFileAnyTypeFieldMultiTSWithoutSDA *pt(_content);
  if(!pt)
    throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTS::contentNotNull : the non const content pointer is null !");
  MEDFileIntFieldMultiTSWithoutSDA *ret=dynamic_cast<MEDFileIntFieldMultiTSWithoutSDA *>(pt);
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTS::contentNotNull : the non const content pointer is not null but it is not of type int ! Reason is maybe that the read field has not the type INT32 !");
  return ret;
}

MEDFileIntFieldMultiTS::MEDFileIntFieldMultiTS()
{
  _content=new MEDFileIntFieldMultiTSWithoutSDA;
}

MEDFileIntFieldMultiTS::MEDFileIntFieldMultiTS(const MEDFileIntFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent):MEDFileAnyTypeFieldMultiTS(other,shallowCopyOfContent)
{
}

MEDFileIntFieldMultiTS::MEDFileIntFieldMultiTS(const char *fileName) throw(INTERP_KERNEL::Exception)
try:MEDFileAnyTypeFieldMultiTS(fileName)
{
}
catch(INTERP_KERNEL::Exception& e)
  { throw e; }

MEDFileIntFieldMultiTS::MEDFileIntFieldMultiTS(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
try:MEDFileAnyTypeFieldMultiTS(fileName,fieldName)
{
}
catch(INTERP_KERNEL::Exception& e)
  { throw e; }

DataArrayInt *MEDFileIntFieldMultiTS::getUndergroundDataArray(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  return static_cast<DataArrayInt *>(contentNotNull()->getUndergroundDataArray(iteration,order));
}

//= MEDFileFields

MEDFileFields *MEDFileFields::New()
{
  return new MEDFileFields;
}

MEDFileFields *MEDFileFields::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFields(fileName);
}

std::size_t MEDFileFields::getHeapMemorySize() const
{
  std::size_t ret=_fields.capacity()*sizeof(MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA>);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    if((const MEDFileAnyTypeFieldMultiTSWithoutSDA *)*it)
      ret+=(*it)->getHeapMemorySize();
  return ret+MEDFileFieldGlobsReal::getHeapMemorySize();
}

MEDFileFields *MEDFileFields::deepCpy() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileFields> ret=shallowCpy();
  std::size_t i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      if((const MEDFileAnyTypeFieldMultiTSWithoutSDA*)*it)
        ret->_fields[i]=(*it)->deepCpy();
    }
  ret->deepCpyGlobs(*this);
  return ret.retn();
}

MEDFileFields *MEDFileFields::shallowCpy() const throw(INTERP_KERNEL::Exception)
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
std::vector< std::pair<int,int> > MEDFileFields::getCommonIterations(bool& areThereSomeForgottenTS) const throw(INTERP_KERNEL::Exception)
{
  std::set< std::pair<int,int> > s;
  bool firstShot=true;
  areThereSomeForgottenTS=false;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
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

std::vector<std::string> MEDFileFields::getFieldsNames() const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> ret(_fields.size());
  int i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *f=(*it);
      if(f)
        {
          ret[i]=f->getName();
        }
      else
        {
          std::ostringstream oss; oss << "MEDFileFields::getFieldsNames : At rank #" << i << " field is not defined !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return ret;
}

std::vector<std::string> MEDFileFields::getMeshesNames() const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
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
  int nbOfFields=getNumberOfFields();
  std::string startLine(bkOffset,' ');
  oss << startLine << "There are " << nbOfFields << " fields in this :" << std::endl;
  int i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
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
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
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

MEDFileFields::MEDFileFields(const char *fileName) throw(INTERP_KERNEL::Exception)
try:MEDFileFieldGlobsReal(fileName)
  {
    MEDFileUtilities::CheckFileForRead(fileName);
    MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
    int nbFields=MEDnField(fid);
    _fields.resize(nbFields);
    med_field_type typcha;
    for(int i=0;i<nbFields;i++)
      {
        std::vector<std::string> infos;
        std::string fieldName,dtunit;
        int nbOfStep=MEDFileAnyTypeField1TS::LocateField2(fid,fileName,i,false,fieldName,typcha,infos,dtunit);
        switch(typcha)
          {
          case MED_FLOAT64:
            {
              _fields[i]=MEDFileFieldMultiTSWithoutSDA::New(fid,fieldName.c_str(),typcha,infos,nbOfStep,dtunit);
              break;
            }
          case MED_INT32:
            {
              _fields[i]=MEDFileIntFieldMultiTSWithoutSDA::New(fid,fieldName.c_str(),typcha,infos,nbOfStep,dtunit);
              break;
            }
          default:
            {
              std::ostringstream oss; oss << "constructor MEDFileFields(fileName) : file \'" << fileName << "\' at pos #" << i << " field has name \'" << fieldName << "\' but the type of field is not in [MED_FLOAT64, MED_INT32] !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          }
      }
    loadAllGlobals(fid);
  }
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

void MEDFileFields::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  int i=0;
  writeGlobals(fid,*this);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *elt=*it;
      if(!elt)
        {
          std::ostringstream oss; oss << "MEDFileFields::write : at rank #" << i << "/" << _fields.size() << " field is empty !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      elt->writeLL(fid,*this);
    }
}

void MEDFileFields::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,medmod);
  writeLL(fid);
}

std::vector<std::string> MEDFileFields::getPflsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTSWithoutSDA > >::const_iterator it=_fields.begin();it!=_fields.end();it++)
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
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTSWithoutSDA > >::const_iterator it=_fields.begin();it!=_fields.end();it++)
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

std::vector<std::string> MEDFileFields::getPflsReallyUsedMulti() const
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTSWithoutSDA > >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsedMulti2();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

std::vector<std::string> MEDFileFields::getLocsReallyUsedMulti() const
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTSWithoutSDA > >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsed2();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

void MEDFileFields::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTSWithoutSDA > >::iterator it=_fields.begin();it!=_fields.end();it++)
    (*it)->changePflsRefsNamesGen2(mapOfModif);
}

void MEDFileFields::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTSWithoutSDA > >::iterator it=_fields.begin();it!=_fields.end();it++)
    (*it)->changeLocsRefsNamesGen2(mapOfModif);
}

void MEDFileFields::resize(int newSize) throw(INTERP_KERNEL::Exception)
{
  _fields.resize(newSize);
}

void MEDFileFields::pushFields(const std::vector<MEDFileAnyTypeFieldMultiTS *>& fields) throw(INTERP_KERNEL::Exception)
{
  for(std::vector<MEDFileAnyTypeFieldMultiTS *>::const_iterator it=fields.begin();it!=fields.end();it++)
    pushField(*it);
}

void MEDFileFields::pushField(MEDFileAnyTypeFieldMultiTS *field) throw(INTERP_KERNEL::Exception)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileFields::pushMesh : invalid input pointer ! should be different from 0 !");
  _fields.push_back(field->getContent());
  appendGlobs(*field,1e-12);
}

void MEDFileFields::setFieldAtPos(int i, MEDFileAnyTypeFieldMultiTS *field) throw(INTERP_KERNEL::Exception)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileFields::setFieldAtPos : invalid input pointer ! should be different from 0 !");
  if(i>=(int)_fields.size())
    _fields.resize(i+1);
  _fields[i]=field->getContent();
  appendGlobs(*field,1e-12);
}

void MEDFileFields::destroyFieldAtPos(int i) throw(INTERP_KERNEL::Exception)
{
  destroyFieldsAtPos(&i,&i+1);
}

void MEDFileFields::destroyFieldsAtPos(const int *startIds, const int *endIds) throw(INTERP_KERNEL::Exception)
{
  std::vector<bool> b(_fields.size(),true);
  for(const int *i=startIds;i!=endIds;i++)
    {
      if(*i<0 || *i>=(int)_fields.size())
        {
          std::ostringstream oss; oss << "MEDFileFields::destroyFieldsAtPos : Invalid given id in input (" << *i << ") should be in [0," << _fields.size() << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      b[*i]=false;
    }
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> > fields(std::count(b.begin(),b.end(),true));
  std::size_t j=0;
  for(std::size_t i=0;i<_fields.size();i++)
    if(b[i])
      fields[j++]=_fields[i];
  _fields=fields;
}

void MEDFileFields::destroyFieldsAtPos2(int bg, int end, int step) throw(INTERP_KERNEL::Exception)
{
  static const char msg[]="MEDFileFields::destroyFieldsAtPos2";
  int nbOfEntriesToKill=DataArrayInt::GetNumberOfItemGivenBESRelative(bg,end,step,msg);
  std::vector<bool> b(_fields.size(),true);
  int k=bg;
  for(int i=0;i<nbOfEntriesToKill;i++,k+=step)
    {
      if(k<0 || k>=(int)_fields.size())
        {
          std::ostringstream oss; oss << "MEDFileFields::destroyFieldsAtPos2 : Invalid given id in input (" << k << ") should be in [0," << _fields.size() << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      b[k]=false;
    }
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> > fields(std::count(b.begin(),b.end(),true));
  std::size_t j=0;
  for(std::size_t i=0;i<_fields.size();i++)
    if(b[i])
      fields[j++]=_fields[i];
  _fields=fields;
}

bool MEDFileFields::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception)
{
  bool ret=false;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::iterator it=_fields.begin();it!=_fields.end();it++)
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
 * \return If true a renumbering has been performed. The structure in \a this has been modified. If false, nothing has been done: it is typically the case if \a meshName is not refered by any 
 *         field in \a this.
 */
bool MEDFileFields::renumberEntitiesLyingOnMesh(const char *meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N) throw(INTERP_KERNEL::Exception)
{
  bool ret=false;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::iterator it=_fields.begin();it!=_fields.end();it++)
    {
      MEDFileAnyTypeFieldMultiTSWithoutSDA *fmts(*it);
      if(fmts)
        {
          ret=fmts->renumberEntitiesLyingOnMesh(meshName,oldCode,newCode,renumO2N,*this) || ret;
        }
    }
  return ret;
}

MEDFileAnyTypeFieldMultiTS *MEDFileFields::getFieldAtPos(int i) const throw(INTERP_KERNEL::Exception)
{
  if(i<0 || i>=(int)_fields.size())
    {
      std::ostringstream oss; oss << "MEDFileFields::getFieldAtPos : Invalid given id in input (" << i << ") should be in [0," << _fields.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *fmts=_fields[i];
  if(!fmts)
    return 0;
  MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTS> ret;
  const MEDFileFieldMultiTSWithoutSDA *fmtsC=dynamic_cast<const MEDFileFieldMultiTSWithoutSDA *>(fmts);
  const MEDFileIntFieldMultiTSWithoutSDA *fmtsC2=dynamic_cast<const MEDFileIntFieldMultiTSWithoutSDA *>(fmts);
  if(fmtsC)
    ret=MEDFileFieldMultiTS::New(*fmtsC,false);
  else if(fmtsC2)
    ret=MEDFileIntFieldMultiTS::New(*fmtsC2,false);
  else
    {
      std::ostringstream oss; oss << "MEDFileFields::getFieldAtPos : At pos #" << i << " field is neither double (FLOAT64) nor integer (INT32) !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  ret->shallowCpyGlobs(*this);
  return ret.retn();
}

/*!
 * Return a shallow copy of \a this reduced to the fields ids defined in [ \a startIds , endIds ).
 * This method is accessible in python using __getitem__ with a list in input.
 * \return a new object that the caller should deal with.
 */
MEDFileFields *MEDFileFields::buildSubPart(const int *startIds, const int *endIds) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileFields> ret=shallowCpy();
  std::size_t sz=std::distance(startIds,endIds);
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> > fields(sz);
  int j=0;
  for(const int *i=startIds;i!=endIds;i++,j++)
    {
      if(*i<0 || *i>=(int)_fields.size())
        {
          std::ostringstream oss; oss << "MEDFileFields::buildSubPart : Invalid given id in input (" << *i << ") should be in [0," << _fields.size() << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      fields[j]=_fields[*i];
    }
  ret->_fields=fields;
  return ret.retn();
}

MEDFileAnyTypeFieldMultiTS *MEDFileFields::getFieldWithName(const char *fieldName) const throw(INTERP_KERNEL::Exception)
{
  return getFieldAtPos(getPosFromFieldName(fieldName));
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
MEDFileFields *MEDFileFields::partOfThisLyingOnSpecifiedMeshName(const char *meshName) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileFields> ret=MEDFileFields::New();
  ret->shallowCpyOnlyUsedGlobs(*this);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *cur=(*it);
      if(!cur)
        continue;
      if(cur->getMeshName()==meshName)
        {
          cur->incrRef();
          MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> cur2(const_cast<MEDFileAnyTypeFieldMultiTSWithoutSDA *>(cur));
          ret->_fields.push_back(cur2);
        }
    }
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
MEDFileFields *MEDFileFields::partOfThisLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileFields> ret=MEDFileFields::New();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *cur=(*it);
      if(!cur)
        continue;
      MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> elt=cur->partOfThisLyingOnSpecifiedTimeSteps(timeSteps);
      ret->_fields.push_back(elt);
    }
  ret->shallowCpyOnlyUsedGlobs(*this);
  return ret.retn();
}

/*!
 * \sa MEDFileFields::getCommonIterations, MEDFileFields::partOfThisLyingOnSpecifiedTimeSteps
 */
MEDFileFields *MEDFileFields::partOfThisNotLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileFields> ret=MEDFileFields::New();
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *cur=(*it);
      if(!cur)
        continue;
      MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> elt=cur->partOfThisNotLyingOnSpecifiedTimeSteps(timeSteps);
      if(elt->getNumberOfTS()!=0)
        ret->_fields.push_back(elt);
    }
  ret->shallowCpyOnlyUsedGlobs(*this);
  return ret.retn();
}

MEDFileFieldsIterator *MEDFileFields::iterator() throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldsIterator(this);
}

int MEDFileFields::getPosFromFieldName(const char *fieldName) const throw(INTERP_KERNEL::Exception)
{
  std::string tmp(fieldName);
  std::vector<std::string> poss;
  for(std::size_t i=0;i<_fields.size();i++)
    {
      const MEDFileAnyTypeFieldMultiTSWithoutSDA *f=_fields[i];
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
  throw INTERP_KERNEL::Exception(oss.str().c_str());
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
