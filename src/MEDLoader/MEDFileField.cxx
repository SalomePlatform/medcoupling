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

void MEDFileFieldPerMeshPerTypePerDisc::assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  _type=field->getTypeOfField();
  const DataArrayDouble *da=field->getArray();
  _start=start;
  switch(_type)
    {
    case ON_CELLS:
      {
        getArray()->setContigPartOfSelectedValues2(_start,da,offset,offset+nbOfCells,1);
        _end=_start+nbOfCells;
        _nval=nbOfCells;
        break;
      }
    case ON_GAUSS_NE:
      {
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=field->getDiscretization()->getOffsetArr(field->getMesh());
        const int *arrPtr=arr->getConstPointer();
        getArray()->setContigPartOfSelectedValues2(_start,da,arrPtr[offset],arrPtr[offset+nbOfCells],1);
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
        getArray()->setContigPartOfSelectedValues(_start,da,da4);
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
 * Leaf method of field with profile assignement.
 * @param pflName input containing name of profile if any. 0 if no profile.
 * @param multiTypePfl input containing the profile array \b including \b all \b types. This array is usefull only for GAUSS_NE.
 * @param idsInPfl input containing the ids in the profile 'multiTypePfl' concerning the current geo type.
 */
void MEDFileFieldPerMeshPerTypePerDisc::assignFieldProfile(int& start, const char *pflName, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  if(pflName)
    _profile=pflName;
  else
    _profile.clear();
  _type=field->getTypeOfField();
  const DataArrayDouble *da=field->getArray();
  _start=start;
  switch(_type)
    {
    case ON_NODES:
      {
         _nval=idsInPfl->getNumberOfTuples();
         getArray()->setContigPartOfSelectedValues2(_start,da,0,da->getNumberOfTuples(),1);
         _end=_start+_nval;
         break;
      }
    case ON_CELLS:
      {
        _nval=idsInPfl->getNumberOfTuples();
        getArray()->setContigPartOfSelectedValues(_start,da,idsInPfl);
        _end=_start+_nval;
        break;
      }
    case ON_GAUSS_NE:
      {
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr=field->getDiscretization()->getOffsetArr(mesh);
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr2=arr->deltaShiftIndex();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> arr3=arr2->selectByTupleId(multiTypePfl->getConstPointer(),multiTypePfl->getConstPointer()+multiTypePfl->getNumberOfTuples());
        arr3->computeOffsets2();
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp=idsInPfl->buildExplicitArrByRanges(arr3);
        int trueNval=tmp->getNumberOfTuples();
        _nval=idsInPfl->getNumberOfTuples();
        getArray()->setContigPartOfSelectedValues(_start,da,tmp);
        _end=_start+trueNval;
        break;
      }
    case ON_GAUSS_PT:
      {
        throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::assignFieldProfile : not implemented yet for profiles on gauss points !");
      }
    default:
      throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::assignFieldProfile : not implemented yet for such discretization type of field !");
    }
  start=_end;
}

void MEDFileFieldPerMeshPerTypePerDisc::assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  _start=start;
  _nval=field->getArray()->getNumberOfTuples();
  getArray()->setContigPartOfSelectedValues2(_start,field->getArray(),0,_nval,1);
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

void MEDFileFieldPerMeshPerTypePerDisc::prepareLoading(med_idt fid, int profileIt, int& start) throw(INTERP_KERNEL::Exception)
{
  INTERP_KERNEL::AutoPtr<char> locname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> pflname=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  std::string fieldName=getName();
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

void MEDFileFieldPerMeshPerTypePerDisc::finishLoading(med_idt fid, int profileIt, int ft) throw(INTERP_KERNEL::Exception)
{
  std::string fieldName=getName();
  std::string meshName=getMeshName();
  int iteration=getIteration();
  int order=getOrder();
  TypeOfField type=getType();
  INTERP_KERNEL::NormalizedCellType geoType=getGeoType();
  med_geometry_type mgeoti;
  med_entity_type menti=MEDFileFieldPerMeshPerType::ConvertIntoMEDFileType(type,geoType,mgeoti);
  DataArrayDouble *arr=getArray();
  double *startFeeding=arr->getPointer()+_start*arr->getNumberOfComponents();
  switch(ft)
    {
    case 0:
      {
        MEDfieldValueWithProfileRd(fid,fieldName.c_str(),iteration,order,menti,mgeoti,MED_COMPACT_PFLMODE,
                                   _profile.c_str(),MED_FULL_INTERLACE,MED_ALL_CONSTITUENT,reinterpret_cast<unsigned char*>(startFeeding));
        break;
      }
    case 1:
      {
        INTERP_KERNEL::AutoPtr<int> tmpp=new int[(_end-_start)*arr->getNumberOfComponents()];
        MEDfieldValueWithProfileRd(fid,fieldName.c_str(),iteration,order,menti,mgeoti,MED_COMPACT_PFLMODE,
                                   _profile.c_str(),MED_FULL_INTERLACE,MED_ALL_CONSTITUENT,reinterpret_cast<unsigned char*>((int *)tmpp));
        std::copy((const int *)tmpp,(const int *)tmpp+(_end-_start)*arr->getNumberOfComponents(),startFeeding);
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("Error on array reading ! Unrecognized type of field ! Should be in FLOAT64 or INT32 !");
    }
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

std::string MEDFileFieldPerMeshPerTypePerDisc::getName() const
{
  return _father->getName();
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

DataArrayDouble *MEDFileFieldPerMeshPerTypePerDisc::getArray()
{
  return _father->getArray();
}

const DataArrayDouble *MEDFileFieldPerMeshPerTypePerDisc::getArray() const
{
  const MEDFileFieldPerMeshPerType *fath=_father;
  return fath->getArray();
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

void MEDFileFieldPerMeshPerTypePerDisc::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  TypeOfField type=getType();
  INTERP_KERNEL::NormalizedCellType geoType=getGeoType();
  med_geometry_type mgeoti;
  med_entity_type menti=MEDFileFieldPerMeshPerType::ConvertIntoMEDFileType(type,geoType,mgeoti);
  const DataArrayDouble *arr=getArray();
  const double *locToWrite=arr->getConstPointer()+_start*arr->getNumberOfComponents();
  MEDfieldValueWithProfileWr(fid,getName().c_str(),getIteration(),getOrder(),getTime(),menti,mgeoti,
                             MED_COMPACT_PFLMODE,_profile.c_str(),_localization.c_str(),MED_FULL_INTERLACE,MED_ALL_CONSTITUENT,_nval,
                             reinterpret_cast<const unsigned char*>(locToWrite));
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
  std::size_t found=-1;
  for(std::size_t i=0;i<nbOfType && found==-1;i++)
    if(getGeoType()==(INTERP_KERNEL::NormalizedCellType)codeOfMesh[3*i])
      found=i;
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
      int offset=codeOfMesh[3*found+2];
      for(const int *pflId=pfl->begin();pflId!=pfl->end();pflId++)
        {
          if(*pflId<codeOfMesh[3*found+1])
            *work++=offset+*pflId;
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
  std::set<int> diffVals=newGeoTypesEltIdsAllGather->getDifferentValues();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> renumEltIds=newGeoTypesEltIdsAllGather->buildPermArrPerLevel();
  //
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> renumTupleIds=newGeoTypesPerChunk4->buildPermArrPerLevel();
  //
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arrPart=arr->substr(offset,offset+szTuples);
  arrPart->renumberInPlace(renumTupleIds->begin());
  arr->setPartOfValues1(arrPart,offset,offset+szTuples,1,0,arrPart->getNumberOfComponents(),1);
  bool ret=false;
  std::set<int>::const_iterator idIt=diffVals.begin();
  std::list<const MEDFileFieldPerMeshPerTypePerDisc *> li(entriesOnSameDisc.begin(),entriesOnSameDisc.end());
  int offset2=0;
  for(std::size_t i=0;i<diffVals.size();i++,idIt++)
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

MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerType::NewOnRead(med_idt fid, MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMeshPerType(fid,fath,type,geoType);
}

MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerType::New(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMeshPerType(fath,geoType);
}

void MEDFileFieldPerMeshPerType::assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  std::vector<int> pos=addNewEntryIfNecessary(field,offset,nbOfCells);
  for(std::vector<int>::const_iterator it=pos.begin();it!=pos.end();it++)
    _field_pm_pt_pd[*it]->assignFieldNoProfile(start,offset,nbOfCells,field,glob);
}

void MEDFileFieldPerMeshPerType::assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, DataArrayInt *locIds, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  std::vector<int> pos=addNewEntryIfNecessary(field,idsInPfl);
  if(locIds)
    {
      //
      std::string pflName(locIds->getName());
      if(pflName.empty())
        throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerType::assignFieldProfile : existing profile with empty name !");
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_geo_type);
      std::ostringstream oss; oss << pflName << "_" <<  cm.getRepr();
      locIds->setName(oss.str().c_str());
      glob.appendProfile(locIds);
      //
      for(std::vector<int>::const_iterator it=pos.begin();it!=pos.end();it++)
        _field_pm_pt_pd[*it]->assignFieldProfile(start,oss.str().c_str(),multiTypePfl,idsInPfl,field,mesh,glob);
    }
  else
    {
      for(std::vector<int>::const_iterator it=pos.begin();it!=pos.end();it++)
        _field_pm_pt_pd[*it]->assignFieldProfile(start,0,multiTypePfl,idsInPfl,field,mesh,glob);
    }
}

void MEDFileFieldPerMeshPerType::assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  _field_pm_pt_pd.resize(1);
  _field_pm_pt_pd[0]=MEDFileFieldPerMeshPerTypePerDisc::New(this,ON_NODES,-3);
  _field_pm_pt_pd[0]->assignNodeFieldNoProfile(start,field,glob);
}

void MEDFileFieldPerMeshPerType::assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  std::string pflName(pfl->getName());
  if(pflName.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerType::assignNodeFieldProfile : existing profile with empty name !");
  std::ostringstream oss; oss << pflName << "_NODE";
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> pfl2=pfl->deepCpy();
  pfl2->setName(oss.str().c_str());
  glob.appendProfile(pfl2);
  //
  _field_pm_pt_pd.resize(1);
  _field_pm_pt_pd[0]=MEDFileFieldPerMeshPerTypePerDisc::New(this,ON_NODES,-3);
  _field_pm_pt_pd[0]->assignFieldProfile(start,oss.str().c_str(),pfl,pfl2,field,0,glob);//mesh is not requested so 0 is send.
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

std::string MEDFileFieldPerMeshPerType::getName() const
{
  return _father->getName();
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

DataArrayDouble *MEDFileFieldPerMeshPerType::getArray()
{
  return _father->getArray();
}

const DataArrayDouble *MEDFileFieldPerMeshPerType::getArray() const
{
  const MEDFileFieldPerMesh *fath=_father;
  return fath->getArray();
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

MEDFileFieldPerMeshPerType::MEDFileFieldPerMeshPerType(med_idt fid, MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception):_father(fath),_geo_type(geoType)
{
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> locName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  med_geometry_type mgeoti;
  med_entity_type menti=ConvertIntoMEDFileType(type,geoType,mgeoti);
  int nbProfiles=MEDfieldnProfile(fid,getName().c_str(),getIteration(),getOrder(),menti,mgeoti,pflName,locName);
  _field_pm_pt_pd.resize(nbProfiles);
  for(int i=0;i<nbProfiles;i++)
    {
      _field_pm_pt_pd[i]=MEDFileFieldPerMeshPerTypePerDisc::NewOnRead(this,type,i+1);
    }
}

void MEDFileFieldPerMeshPerType::prepareLoading(med_idt fid, int &start) throw(INTERP_KERNEL::Exception)
{
  int pflId=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++,pflId++)
    {
      (*it)->prepareLoading(fid,pflId+1,start);
    }
}

void MEDFileFieldPerMeshPerType::finishLoading(med_idt fid, int ft) throw(INTERP_KERNEL::Exception)
{
  int pflId=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >::iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++,pflId++)
    {
      (*it)->finishLoading(fid,pflId+1,ft);
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

MEDFileFieldPerMesh *MEDFileFieldPerMesh::NewOnRead(med_idt fid, MEDFileField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMesh(fid,fath,meshCsit,meshIteration,meshOrder);
}

MEDFileFieldPerMesh *MEDFileFieldPerMesh::New(MEDFileField1TSWithoutSDA *fath, const MEDCouplingMesh *mesh)
{
  return new MEDFileFieldPerMesh(fath,mesh);
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

void MEDFileFieldPerMesh::assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const std::vector<int>& code, const std::vector<DataArrayInt *>& idsInPflPerType, const std::vector<DataArrayInt *>& idsPerType, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  int nbOfTypes=code.size()/3;
  bool isProfile=false;
  for(int i=0;i<nbOfTypes;i++)
    if(code[3*i+2]!=-1)
      isProfile=true;
  if(!isProfile)
    {
      if(idsInPflPerType.empty())
        assignFieldNoProfileNoRenum(start,code,field,glob);
      else
        assignFieldProfileGeneral(start,multiTypePfl,code,idsInPflPerType,idsPerType,field,mesh,glob);
    }
  else
    assignFieldProfileGeneral(start,multiTypePfl,code,idsInPflPerType,idsPerType,field,mesh,glob);
}

void MEDFileFieldPerMesh::assignFieldNoProfileNoRenum(int& start, const std::vector<int>& code, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  int nbOfTypes=code.size()/3;
  int offset=0;
  for(int i=0;i<nbOfTypes;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)code[3*i];
      int nbOfCells=code[3*i+1];
      int pos=addNewEntryIfNecessary(type);
      _field_pm_pt[pos]->assignFieldNoProfile(start,offset,nbOfCells,field,glob);
      offset+=nbOfCells;
    }
}

/*!
 * This method is the most general one. No optimization is done here.
 */
void MEDFileFieldPerMesh::assignFieldProfileGeneral(int& start, const DataArrayInt *multiTypePfl, const std::vector<int>& code, const std::vector<DataArrayInt *>& idsInPflPerType, const std::vector<DataArrayInt *>& idsPerType, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  int nbOfTypes=code.size()/3;
  for(int i=0;i<nbOfTypes;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)code[3*i];
      int pos=addNewEntryIfNecessary(type);
      DataArrayInt *pfl=0;
      if(code[3*i+2]!=-1)
        pfl=idsPerType[code[3*i+2]];
      _field_pm_pt[pos]->assignFieldProfile(start,multiTypePfl,idsInPflPerType[i],pfl,field,mesh,glob);
    }
}

void MEDFileFieldPerMesh::assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  int pos=addNewEntryIfNecessary(INTERP_KERNEL::NORM_ERROR);
  _field_pm_pt[pos]->assignNodeFieldNoProfile(start,field,glob);
}

void MEDFileFieldPerMesh::assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  int pos=addNewEntryIfNecessary(INTERP_KERNEL::NORM_ERROR);
  _field_pm_pt[pos]->assignNodeFieldProfile(start,pfl,field,glob);
}

void MEDFileFieldPerMesh::prepareLoading(med_idt fid, int& start) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->prepareLoading(fid,start);
}

void MEDFileFieldPerMesh::finishLoading(med_idt fid, int ft) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->finishLoading(fid,ft);
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

DataArrayDouble *MEDFileFieldPerMesh::getArray()
{
  return _father->getOrCreateAndGetArray();
}

const DataArrayDouble *MEDFileFieldPerMesh::getArray() const
{
  const MEDFileField1TSWithoutSDA *fath=_father;
  return fath->getOrCreateAndGetArray();
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
  DataArrayDouble *arr=getUndergroundDataArrayExt(entries);
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

MEDCouplingFieldDouble *MEDFileFieldPerMesh::getFieldOnMeshAtLevel(TypeOfField type, const MEDFileFieldGlobsReal *glob, const MEDCouplingMesh *mesh, bool& isPfl) const throw(INTERP_KERNEL::Exception)
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
          return finishField2(type,glob,dads,locs,geoTypes,mesh,arr,isPfl);
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

DataArrayDouble *MEDFileFieldPerMesh::getFieldOnMeshAtLevelWithPfl(TypeOfField type, const MEDCouplingMesh *mesh, DataArrayInt *&pfl, const MEDFileFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception)
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
      std::ostringstream oss; oss << "MEDFileFieldPerMesh::getFieldOnMeshAtLevelWithPfl : " << "The field \"" << getName() << "\" exists but not with such spatial discretization or such dimension specified !";
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

DataArrayDouble *MEDFileFieldPerMesh::getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception)
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
  return _father->getUndergroundDataArray();
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
 * 'dads' and 'locs' input parameters have the same number of elements.
 */
MEDCouplingFieldDouble *MEDFileFieldPerMesh::finishField(TypeOfField type, const MEDFileFieldGlobsReal *glob,
                                                         const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs,
                                                         const MEDCouplingMesh *mesh, bool& isPfl) const throw(INTERP_KERNEL::Exception)
{
  isPfl=false;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(type,ONE_TIME);
  ret->setMesh(mesh); ret->setName(getName().c_str()); ret->setTime(getTime(),getIteration(),getOrder()); ret->setTimeUnit(getDtUnit().c_str());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> da=getArray()->selectByTupleRanges(dads);
  const std::vector<std::string>& infos=getInfo();
  da->setInfoOnComponents(infos);
  da->setName("");
  ret->setArray(da);
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
  //
  ret->incrRef();
  return ret;
}

/*!
 * This method is an extension of MEDFileFieldPerMesh::finishField method. It deals with profiles. This method should be called when type is different from ON_NODES.
 * 'dads', 'locs' and 'geoTypes' input parameters have the same number of elements.
 * No check of this is performed. 'da' array contains an array in old2New style to be applyied to mesh to obtain the right support.
 * The order of cells in the returned field is those imposed by the profile.
 */
MEDCouplingFieldDouble *MEDFileFieldPerMesh::finishField2(TypeOfField type, const MEDFileFieldGlobsReal *glob,
                                                          const std::vector<std::pair<int,int> >& dads, const std::vector<int>& locs,
                                                          const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes,
                                                          const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl) const throw(INTERP_KERNEL::Exception)
{
  if(da->isIdentity())
    {
      int nbOfTuples=da->getNumberOfTuples();
      if(nbOfTuples==mesh->getNumberOfCells())
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
MEDCouplingFieldDouble *MEDFileFieldPerMesh::finishField3(const MEDFileFieldGlobsReal *glob,
                                                          const std::vector<std::pair<int,int> >& dads, const std::vector<int>& locs,
                                                          const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl) const throw(INTERP_KERNEL::Exception)
{
  if(da->isIdentity())
    {
      int nbOfTuples=da->getNumberOfTuples();
      const std::vector<INTERP_KERNEL::NormalizedCellType> geoTypes2(1,INTERP_KERNEL::NORM_ERROR);
      if(nbOfTuples==ComputeNbOfElems(glob,ON_NODES,geoTypes2,dads,locs))//No problem for NORM_ERROR because it is in context of node
        return finishField(ON_NODES,glob,dads,locs,mesh,isPfl);
    }
  // Treatment of particular case where nodal field on pfl is requested with a meshDimRelToMax=1.
  const MEDCouplingUMesh *meshu=dynamic_cast<const MEDCouplingUMesh *>(mesh);
  if(meshu)
    {
      if(meshu->getNodalConnectivity()==0)
        {
          MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=finishField(ON_CELLS,glob,dads,locs,mesh,isPfl);
          int nb=da->getNbOfElems();
          const int *ptr=da->getConstPointer();
          MEDCouplingUMesh *meshuc=const_cast<MEDCouplingUMesh *>(meshu);
          meshuc->allocateCells(nb);
          for(int i=0;i<nb;i++)
            meshuc->insertNextCell(INTERP_KERNEL::NORM_POINT1,1,ptr+i);
          meshuc->finishInsertingCells();
          ret->setMesh(meshuc);
          ret->checkCoherency();
          ret->incrRef();
          return ret;
        }
    }
  //
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
    {
      std::ostringstream oss; oss << "MEDFileFieldPerMesh::finishField3 : The field on nodes lies on a node profile so that it is impossible to find a submesh having exactly the same nodes of that profile !!!";
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
DataArrayDouble *MEDFileFieldPerMesh::finishField4(const std::vector<std::pair<int,int> >& dads, const DataArrayInt *pflIn, int nbOfElems, DataArrayInt *&pflOut) const throw(INTERP_KERNEL::Exception)
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
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> da=getArray()->selectByTupleRanges(dads);
  const std::vector<std::string>& infos=getInfo();
  int nbOfComp=infos.size();
  for(int i=0;i<nbOfComp;i++)
    da->setInfoOnComponent(i,infos[i].c_str());
  safePfl->incrRef();
  da->incrRef();
  return da;
}

MEDFileFieldPerMesh::MEDFileFieldPerMesh(med_idt fid, MEDFileField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder) throw(INTERP_KERNEL::Exception):_mesh_iteration(meshIteration),_mesh_order(meshOrder),
                                                                                                                                                                       _mesh_csit(meshCsit),_father(fath)
{
  INTERP_KERNEL::AutoPtr<char> meshName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> locName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  for(int i=0;i<MED_N_CELL_FIXED_GEO;i++)
    {
      int nbProfile=MEDfield23nProfile(fid,getName().c_str(),getIteration(),getOrder(),MED_CELL,typmai[i],_mesh_csit,meshName,pflName,locName);
      if(nbProfile>0)
        {
          _field_pm_pt.push_back(MEDFileFieldPerMeshPerType::NewOnRead(fid,this,ON_CELLS,typmai2[i]));
          _mesh_name=MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1);
        }
      nbProfile=MEDfield23nProfile(fid,getName().c_str(),getIteration(),getOrder(),MED_NODE_ELEMENT,typmai[i],_mesh_csit,meshName,pflName,locName);
      if(nbProfile>0)
        {
          _field_pm_pt.push_back(MEDFileFieldPerMeshPerType::NewOnRead(fid,this,ON_GAUSS_NE,typmai2[i]));
          _mesh_name=MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1);
        }
    }
  int nbProfile=MEDfield23nProfile(fid,getName().c_str(),getIteration(),getOrder(),MED_NODE,MED_NONE,_mesh_csit,meshName,pflName,locName);
  if(nbProfile>0)
    {
      _field_pm_pt.push_back(MEDFileFieldPerMeshPerType::NewOnRead(fid,this,ON_NODES,INTERP_KERNEL::NORM_ERROR));
      _mesh_name=MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1);
    }
}

MEDFileFieldPerMesh::MEDFileFieldPerMesh(MEDFileField1TSWithoutSDA *fath, const MEDCouplingMesh *mesh):_father(fath)
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

MEDFileFieldGlobsReal::MEDFileFieldGlobsReal(const char *fname):_globals(MEDFileFieldGlobs::New(fname))
{
}

MEDFileFieldGlobsReal::MEDFileFieldGlobsReal():_globals(MEDFileFieldGlobs::New())
{
}

void MEDFileFieldGlobsReal::simpleRepr(std::ostream& oss) const
{
  oss << "Globals information on fields :" << "\n*******************************\n\n";
  const MEDFileFieldGlobs *glob=_globals;
  if(glob)
    glob->simpleRepr(oss);
  else
    oss << "NO GLOBAL INFORMATION !\n";
}

MEDFileFieldGlobsReal::~MEDFileFieldGlobsReal()
{
}

void MEDFileFieldGlobsReal::shallowCpyGlobs(const MEDFileFieldGlobsReal& other)
{
  _globals=other._globals;
}

void MEDFileFieldGlobsReal::appendGlobs(const MEDFileFieldGlobsReal& other, double eps) throw(INTERP_KERNEL::Exception)
{
  _globals->appendGlobs(*other._globals,eps);
}

void MEDFileFieldGlobsReal::loadProfileInFile(med_idt fid, int id, const char *pflName) throw(INTERP_KERNEL::Exception)
{
  _globals->loadProfileInFile(fid,id,pflName);
}

void MEDFileFieldGlobsReal::loadProfileInFile(med_idt fid, int id)
{
  _globals->loadProfileInFile(fid,id);
}

void MEDFileFieldGlobsReal::loadGlobals(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  _globals->loadGlobals(fid,*this);
}

void MEDFileFieldGlobsReal::loadAllGlobals(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  _globals->loadAllGlobals(fid);
}

void MEDFileFieldGlobsReal::writeGlobals(med_idt fid, const MEDFileWritable& opt) const throw(INTERP_KERNEL::Exception)
{
  _globals->writeGlobals(fid,opt);
}

std::vector<std::string> MEDFileFieldGlobsReal::getPfls() const
{
  return _globals->getPfls();
}

std::vector<std::string> MEDFileFieldGlobsReal::getLocs() const
{
  return _globals->getLocs();
}

bool MEDFileFieldGlobsReal::existsPfl(const char *pflName) const
{
  return _globals->existsPfl(pflName);
}

bool MEDFileFieldGlobsReal::existsLoc(const char *locName) const
{
  return _globals->existsLoc(locName);
}

std::string MEDFileFieldGlobsReal::createNewNameOfPfl() const throw(INTERP_KERNEL::Exception)
{
  return _globals->createNewNameOfPfl();
}

std::string MEDFileFieldGlobsReal::createNewNameOfLoc() const throw(INTERP_KERNEL::Exception)
{
  return _globals->createNewNameOfLoc();
}

void MEDFileFieldGlobsReal::setFileName(const char *fileName)
{
  _globals->setFileName(fileName);
}

std::vector< std::vector<int> > MEDFileFieldGlobsReal::whichAreEqualProfiles() const
{
  return _globals->whichAreEqualProfiles();
}

std::vector< std::vector<int> > MEDFileFieldGlobsReal::whichAreEqualLocs(double eps) const
{
  return _globals->whichAreEqualLocs(eps);
}

void MEDFileFieldGlobsReal::changePflsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  _globals->changePflsNamesInStruct(mapOfModif);
}

void MEDFileFieldGlobsReal::changeLocsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  _globals->changeLocsNamesInStruct(mapOfModif);
}

/*!
 * This method is a generalization of MEDFileFieldGlobsReal::changePflName.
 * This method contrary to abstract method MEDFileFieldGlobsReal::changePflsRefsNamesGen updates in addition of MEDFileFieldGlobsReal::changePflsRefsNamesGen,
 * the profiles themselves and not only leaves of field.
 */
void MEDFileFieldGlobsReal::changePflsNames(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  changePflsRefsNamesGen(mapOfModif);
  changePflsNamesInStruct(mapOfModif);
}

/*!
 * This method is a generalization of MEDFileFieldGlobsReal::changePflName.
 * This method contrary to abstract method MEDFileFieldGlobsReal::changeLocsRefsNamesGen updates in addition of MEDFileFieldGlobsReal::changeLocsRefsNamesGen,
 * the localizations themselves and not only leaves of field.
 */
void MEDFileFieldGlobsReal::changeLocsNames(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  changeLocsRefsNamesGen(mapOfModif);
  changeLocsNamesInStruct(mapOfModif);
}

/*!
 * This method is a more friendly API but less general method than MEDFileFieldGlobsReal::changePflsNames.
 */
void MEDFileFieldGlobsReal::changePflName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception)
{
  std::vector< std::pair<std::vector<std::string>, std::string > > mapOfModif(1);
  std::pair<std::vector<std::string>, std::string > p(std::vector<std::string>(1,std::string(oldName)),std::string(newName));
  mapOfModif[0]=p;
  changePflsNames(mapOfModif);
}

/*!
 * This method is a more friendly API but less general method than MEDFileFieldGlobsReal::changeLocsNames.
 */
void MEDFileFieldGlobsReal::changeLocName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception)
{
  std::vector< std::pair<std::vector<std::string>, std::string > > mapOfModif(1);
  std::pair<std::vector<std::string>, std::string > p(std::vector<std::string>(1,std::string(oldName)),std::string(newName));
  mapOfModif[0]=p;
  changeLocsNames(mapOfModif);
}

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

int MEDFileFieldGlobsReal::getNbOfGaussPtPerCell(int locId) const throw(INTERP_KERNEL::Exception)
{
  return _globals->getNbOfGaussPtPerCell(locId);
}

int MEDFileFieldGlobsReal::getLocalizationId(const char *loc) const throw(INTERP_KERNEL::Exception)
{
  return _globals->getLocalizationId(loc);
}

const char *MEDFileFieldGlobsReal::getFileName() const
{
  return _globals->getFileName();
}

std::string MEDFileFieldGlobsReal::getFileName2() const
{
  return _globals->getFileName2();
}

const MEDFileFieldLoc& MEDFileFieldGlobsReal::getLocalization(const char *locName) const throw(INTERP_KERNEL::Exception)
{
  return _globals->getLocalization(locName);
}

const MEDFileFieldLoc& MEDFileFieldGlobsReal::getLocalizationFromId(int locId) const throw(INTERP_KERNEL::Exception)
{
  return _globals->getLocalizationFromId(locId);
}

const DataArrayInt *MEDFileFieldGlobsReal::getProfile(const char *pflName) const throw(INTERP_KERNEL::Exception)
{
  return _globals->getProfile(pflName);
}

const DataArrayInt *MEDFileFieldGlobsReal::getProfileFromId(int pflId) const throw(INTERP_KERNEL::Exception)
{
  return _globals->getProfileFromId(pflId);
}

MEDFileFieldLoc& MEDFileFieldGlobsReal::getLocalizationFromId(int locId) throw(INTERP_KERNEL::Exception)
{
  return _globals->getLocalizationFromId(locId);
}

MEDFileFieldLoc& MEDFileFieldGlobsReal::getLocalization(const char *locName) throw(INTERP_KERNEL::Exception)
{
  return _globals->getLocalization(locName);
}

DataArrayInt *MEDFileFieldGlobsReal::getProfile(const char *pflName) throw(INTERP_KERNEL::Exception)
{
  return _globals->getProfile(pflName);
}

DataArrayInt *MEDFileFieldGlobsReal::getProfileFromId(int pflId) throw(INTERP_KERNEL::Exception)
{
  return _globals->getProfileFromId(pflId);
}

void MEDFileFieldGlobsReal::killProfileIds(const std::vector<int>& pflIds) throw(INTERP_KERNEL::Exception)
{
  _globals->killProfileIds(pflIds);
}

void MEDFileFieldGlobsReal::killLocalizationIds(const std::vector<int>& locIds) throw(INTERP_KERNEL::Exception)
{
  _globals->killLocalizationIds(locIds);
}

void MEDFileFieldGlobsReal::appendProfile(DataArrayInt *pfl) throw(INTERP_KERNEL::Exception)
{
  _globals->appendProfile(pfl);
}

void MEDFileFieldGlobsReal::appendLoc(const char *locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w) throw(INTERP_KERNEL::Exception)
{
  _globals->appendLoc(locName,geoType,refCoo,gsCoo,w);
}

/*!
 * This method returns the max dimension of 'this'.
 * This method returns -2 if 'this' is empty, -1 if only nodes are defined.
 */
int MEDFileField1TSWithoutSDA::getDimension() const
{
  int ret=-2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    (*it)->getDimension(ret);
  return ret;
}

void MEDFileField1TSWithoutSDA::CheckMeshDimRel(int meshDimRelToMax) throw(INTERP_KERNEL::Exception)
{
  if(meshDimRelToMax>0)
    throw INTERP_KERNEL::Exception("CheckMeshDimRel : This is a meshDimRel not a meshDimRelExt ! So value should be <=0 !");
}

std::vector<int> MEDFileField1TSWithoutSDA::CheckSBTMesh(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception)
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
      throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::CheckSBTMesh : underlying mesh is not sorted by type as MED file expects !");
    }
  return code;
}

MEDFileField1TSWithoutSDA *MEDFileField1TSWithoutSDA::New(const char *fieldName, int csit, int fieldtype, int iteration, int order, const std::vector<std::string>& infos)
{
  return new MEDFileField1TSWithoutSDA(fieldName,csit,fieldtype,iteration,order,infos);
}

/*!
 * This method copyies tiny info but also preallocated the DataArrayDouble instance in this->_arr.
 * This not allocated it allocates to the size of 'field' array. If already allocated it grows the array to
 * the previous size + the size of the array of the input 'field'.
 * This method returns the position (in tuple id) where to start to feed 'this->_arr'
 */
int MEDFileField1TSWithoutSDA::copyTinyInfoFrom(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception)
{
  std::string name(field->getName());
  getOrCreateAndGetArray()->setName(name.c_str());
  if(name.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::copyTinyInfoFrom : unsupported fields with no name in MED file !");
  const DataArrayDouble *arr=field->getArray();
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::copyTinyInfoFrom : no array set !");
  _dt=field->getTime(_iteration,_order);
  int nbOfComponents=arr->getNumberOfComponents();
  getOrCreateAndGetArray()->setInfoAndChangeNbOfCompo(arr->getInfoOnComponents());
  if(!getOrCreateAndGetArray()->isAllocated())
    {
      _arr->alloc(arr->getNumberOfTuples(),arr->getNumberOfComponents());
      return 0;
    }
  else
    {
      int oldNbOfTuples=getOrCreateAndGetArray()->getNumberOfTuples();
      int newNbOfTuples=oldNbOfTuples+arr->getNumberOfTuples();
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> tmp=DataArrayDouble::New();
      tmp->alloc(newNbOfTuples,nbOfComponents);
      tmp->copyStringInfoFrom(*_arr);
      std::copy(_arr->begin(),_arr->end(),tmp->getPointer());
      _arr=tmp;
      return oldNbOfTuples;
    }
}

std::string MEDFileField1TSWithoutSDA::getName() const
{
  const DataArrayDouble *arr=getOrCreateAndGetArray();
  return arr->getName();
}

void MEDFileField1TSWithoutSDA::setName(const char *name)
{
  DataArrayDouble *arr=getOrCreateAndGetArray();
  arr->setName(name);
}

void MEDFileField1TSWithoutSDA::simpleRepr(int bkOffset, std::ostream& oss, int f1tsId) const
{
  std::string startOfLine(bkOffset,' ');
  oss << startOfLine << "Field on One time Step ";
  if(f1tsId>=0)
    oss << "(" << f1tsId << ") ";
  oss << "on iteration=" << _iteration << " order=" << _order << "." << std::endl;
  oss << startOfLine << "Time attached is : " << _dt << " [" << _dt_unit << "]." << std::endl;
  const DataArrayDouble *arr=_arr;
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

std::string MEDFileField1TSWithoutSDA::getMeshName() const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::getMeshName : No field set !");
  return _field_per_mesh[0]->getMeshName();
}

void MEDFileField1TSWithoutSDA::setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception)
{
  std::string oldName(getMeshName());
  std::vector< std::pair<std::string,std::string> > v(1);
  v[0].first=oldName; v[0].second=newMeshName;
  changeMeshNames(v);
}

bool MEDFileField1TSWithoutSDA::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception)
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

bool MEDFileField1TSWithoutSDA::renumberEntitiesLyingOnMesh(const char *meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N,
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

int MEDFileField1TSWithoutSDA::getMeshIteration() const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::getMeshIteration : No field set !");
  return _field_per_mesh[0]->getMeshIteration();
}

int MEDFileField1TSWithoutSDA::getMeshOrder() const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::getMeshOrder : No field set !");
  return _field_per_mesh[0]->getMeshOrder();
}

int MEDFileField1TSWithoutSDA::getNumberOfComponents() const
{
  return getOrCreateAndGetArray()->getNumberOfComponents();
}

bool MEDFileField1TSWithoutSDA::isDealingTS(int iteration, int order) const
{
  return iteration==_iteration && order==_order;
}

std::pair<int,int> MEDFileField1TSWithoutSDA::getDtIt() const
{
  std::pair<int,int> p;
  fillIteration(p);
  return p;
}

void MEDFileField1TSWithoutSDA::fillIteration(std::pair<int,int>& p) const
{
  p.first=_iteration;
  p.second=_order;
}

void MEDFileField1TSWithoutSDA::fillTypesOfFieldAvailable(std::vector<TypeOfField>& types) const throw(INTERP_KERNEL::Exception)
{
  std::set<TypeOfField> types2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    {
      (*it)->fillTypesOfFieldAvailable(types2);
    }
  std::back_insert_iterator< std::vector<TypeOfField> > bi(types);
  std::copy(types2.begin(),types2.end(),bi);
}

const std::vector<std::string>& MEDFileField1TSWithoutSDA::getInfo() const
{
  const DataArrayDouble *arr=getOrCreateAndGetArray();
  return arr->getInfoOnComponents();
}

std::vector<std::string>& MEDFileField1TSWithoutSDA::getInfo()
{
  DataArrayDouble *arr=getOrCreateAndGetArray();
  return arr->getInfoOnComponents();
}

/*!
 * This method has one input 'mname'. It can be null if the user is the general case where there is only one meshName lying on 'this'
 * This method returns two things.
 * - The absolute dimension of 'this' in first parameter. 
 * - The available ext levels relative to the absolute dimension returned in first parameter. These relative levels are relative
 *   to the first output parameter. The values in 'levs' will be returned in decreasing order.
 *
 * This method is designed for MEDFileField1TS instances that have a discritization ON_CELLS, ON_GAUSS_NE and ON_GAUSS.
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
int MEDFileField1TSWithoutSDA::getNonEmptyLevels(const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception)
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

std::vector<TypeOfField> MEDFileField1TSWithoutSDA::getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception)
{
  std::vector<TypeOfField> ret;
  fillTypesOfFieldAvailable(ret);
  return ret;
}

/*!
 * entry point for users that want to iterate into MEDFile DataStructure without any overhead.
 */
std::vector< std::vector< std::pair<int,int> > > MEDFileField1TSWithoutSDA::getFieldSplitedByType(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
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
 * entry point for users that want to iterate into MEDFile DataStructure with a reduced overhead because output arrays are extracted (created) specially
 * for the call of this method. That's why the DataArrayDouble instance in returned vector of vector should be dealed by the caller.
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

void MEDFileField1TSWithoutSDA::finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  med_int numdt,numit;
  med_float dt;
  med_int nmesh;
  med_bool localMesh;
  med_int meshnumdt,meshnumit;
  INTERP_KERNEL::AutoPtr<char> meshName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  MEDfieldComputingStepInfo(fid,getName().c_str(),_csit,&numdt,&numit,&_dt);
  MEDfield23ComputingStepMeshInfo(fid,getName().c_str(),_csit,&numdt,&numit,&dt,&nmesh,meshName,&localMesh,&meshnumdt,&meshnumit);
  if(_iteration!=numdt || _order!=numit)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::finishLoading : unexpected exception internal error !");
  _field_per_mesh.resize(nmesh);
  for(int i=0;i<nmesh;i++)
    _field_per_mesh[i]=MEDFileFieldPerMesh::NewOnRead(fid,this,i+1,meshnumdt,meshnumit);
  int start=0;
  for(int i=0;i<nmesh;i++)
    {
      _field_per_mesh[i]->prepareLoading(fid,start);
    }
  getOrCreateAndGetArray()->alloc(start,getNumberOfComponents());
  for(int i=0;i<nmesh;i++)
    {
      _field_per_mesh[i]->finishLoading(fid,_field_type);
    }
}

std::vector<std::string> MEDFileField1TSWithoutSDA::getPflsReallyUsed2() const
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

std::vector<std::string> MEDFileField1TSWithoutSDA::getLocsReallyUsed2() const
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

std::vector<std::string> MEDFileField1TSWithoutSDA::getPflsReallyUsedMulti2() const
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsedMulti();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

std::vector<std::string> MEDFileField1TSWithoutSDA::getLocsReallyUsedMulti2() const
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

void MEDFileField1TSWithoutSDA::changePflsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    (*it)->changePflsRefsNamesGen(mapOfModif);
}

void MEDFileField1TSWithoutSDA::changeLocsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    (*it)->changeLocsRefsNamesGen(mapOfModif);
}

void MEDFileField1TSWithoutSDA::writeLL(med_idt fid, const MEDFileWritable& opts) const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::writeLL : empty field !");
  if(_field_per_mesh.size()>1)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::writeLL : In MED3.0 mode in writting mode only ONE underlying mesh supported !");
  _field_per_mesh[0]->copyOptionsFrom(opts);
  _field_per_mesh[0]->writeLL(fid);
}

/*!
 * SBT means Sort By Type.
 * This method is the most basic method to assign field in this. Basic in sense that no renumbering is done. Underlying mesh in 'field' is globaly ignored except for type contiguous check.
 */
void MEDFileField1TSWithoutSDA::setFieldNoProfileSBT(const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingMesh *mesh=field->getMesh();
  //
  TypeOfField type=field->getTypeOfField();
  std::vector<DataArrayInt *> dummy;
  int start=copyTinyInfoFrom(field);
  if(type!=ON_NODES)
    {
      std::vector<int> code=MEDFileField1TSWithoutSDA::CheckSBTMesh(mesh);
      //
      int pos=addNewEntryIfNecessary(mesh);
      _field_per_mesh[pos]->assignFieldProfile(start,0,code,dummy,dummy,field,0,glob);//mesh is set to 0 because no external mesh is needed to be sent because no profile.
    }
  else
    {
      int pos=addNewEntryIfNecessary(mesh);
      _field_per_mesh[pos]->assignNodeFieldNoProfile(start,field,glob);
    }
}

/*!
 * Generalization of MEDFileField1TSWithoutSDA::setFieldNoProfileSBT method.
 */
void MEDFileField1TSWithoutSDA::setFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  TypeOfField type=field->getTypeOfField();
  int start=copyTinyInfoFrom(field);
  std::vector<DataArrayInt *> idsInPflPerType;
  std::vector<DataArrayInt *> idsPerType;
  std::vector<int> code;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=mesh->getGenMeshAtLevel(meshDimRelToMax);
  if(type!=ON_NODES)
    {
      m->splitProfilePerType(profile,code,idsInPflPerType,idsPerType);
      //
      std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > idsInPflPerType2(idsInPflPerType.size());
      for(std::size_t i=0;i<idsInPflPerType.size();i++)
        idsInPflPerType2[i]=idsInPflPerType[i];
      std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > idsPerType2(idsPerType.size());
      for(std::size_t i=0;i<idsPerType.size();i++)
        idsPerType2[i]=idsPerType[i];
      //
      int pos=addNewEntryIfNecessary(m);
      _field_per_mesh[pos]->assignFieldProfile(start,profile,code,idsInPflPerType,idsPerType,field,m,glob);
    }
  else
    {
      int pos=addNewEntryIfNecessary(m);
      _field_per_mesh[pos]->assignNodeFieldProfile(start,profile,field,glob);
    }
}

MEDCouplingFieldDouble *MEDFileField1TSWithoutSDA::getFieldAtLevel(TypeOfField type, int meshDimRelToMax, const char *mName, int renumPol, const MEDFileFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm;
  if(mName==0)
    mm=MEDFileMesh::New(glob->getFileName(),getMeshName().c_str(),getMeshIteration(),getMeshOrder());
  else
    mm=MEDFileMesh::New(glob->getFileName(),mName,-1,-1);
  return MEDFileField1TSWithoutSDA::getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,glob,mm);
}

MEDCouplingFieldDouble *MEDFileField1TSWithoutSDA::getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol, const MEDFileFieldGlobsReal *glob, const MEDFileMesh *mesh) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=mesh->getGenMeshAtLevel(meshDimRelToMax,false);
  const DataArrayInt *d=mesh->getNumberFieldAtLevel(meshDimRelToMax);
  const DataArrayInt *e=mesh->getNumberFieldAtLevel(1);
  if(meshDimRelToMax==1)
    (static_cast<MEDCouplingUMesh *>((MEDCouplingMesh *)m))->setMeshDimension(0);
  return MEDFileField1TSWithoutSDA::getFieldOnMeshAtLevel(type,renumPol,glob,m,d,e);
}

MEDCouplingFieldDouble *MEDFileField1TSWithoutSDA::getFieldAtTopLevel(TypeOfField type, const char *mName, int renumPol, const MEDFileFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception)
{
   MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm;
  if(mName==0)
    mm=MEDFileMesh::New(glob->getFileName(),getMeshName().c_str(),getMeshIteration(),getMeshOrder());
  else
    mm=MEDFileMesh::New(glob->getFileName(),mName,-1,-1);
  int absDim=getDimension();
  int meshDimRelToMax=absDim-mm->getMeshDimension();
  return MEDFileField1TSWithoutSDA::getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,glob,mm);
}

MEDCouplingFieldDouble *MEDFileField1TSWithoutSDA::getFieldOnMeshAtLevel(TypeOfField type, int renumPol, const MEDFileFieldGlobsReal *glob, const MEDCouplingMesh *mesh, const DataArrayInt *cellRenum, const DataArrayInt *nodeRenum) const throw(INTERP_KERNEL::Exception)
{
  static const char msg1[]="MEDFileField1TSWithoutSDA::getFieldOnMeshAtLevel : request for a renumbered field following mesh numbering whereas it is a profile field !";
  int meshId=getMeshIdFromMeshName(mesh->getName());
  bool isPfl=false;
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=_field_per_mesh[meshId]->getFieldOnMeshAtLevel(type,glob,mesh,isPfl);
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
                std::ostringstream oss; oss << "MEDFileField1TSWithoutSDA::getFieldOnMeshAtLevel : Request of simple renumbering but it seems that underlying mesh \"" << mesh->getName() << "\" of requested field ";
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
                std::ostringstream oss; oss << "MEDFileField1TSWithoutSDA::getFieldOnMeshAtLevel : Request of simple renumbering but it seems that underlying mesh \"" << mesh->getName() << "\" of requested field ";
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
      throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getFieldOnMeshAtLevel : unsupported renum policy ! Dealing with policy 0 1 2 and 3 !");
    }
}

DataArrayDouble *MEDFileField1TSWithoutSDA::getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl, const MEDFileFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=mesh->getGenMeshAtLevel(meshDimRelToMax);
  int meshId=getMeshIdFromMeshName(mesh->getName());
  return _field_per_mesh[meshId]->getFieldOnMeshAtLevelWithPfl(type,m,pfl,glob);
}

/*!
 * This method retrieves direct access to the underground ParaMEDMEM::DataArrayDouble instance. The returned array is not a newly
 * created array so it should \b not be dealed by the caller.
 * This method allows to the user a direct access to the values.
 * This method is quite unusable if there is more than a nodal field or a cell field on single geometric cell type.
 */
DataArrayDouble *MEDFileField1TSWithoutSDA::getUndergroundDataArray() const throw(INTERP_KERNEL::Exception)
{
  const DataArrayDouble *ret=_arr;
  if(ret)
    return const_cast<DataArrayDouble *>(ret);
  else
    return 0;
}

/*!
 * This method returns an array that the caller must deal with (contrary to those returned by MEDFileField1TSWithoutSDA::getUndergroundDataArray method).
 * The returned array is the result of the aggregation of all sub arrays stored in the MED file. So to allow the caller to select the output param
 * 'entries' is returned. This output param is a vector of a pair of 2 pairs. The first pair of pair informs of the geometric type it refers to and the discretization
 * id attached to it. The second pair of pair precise the range [begin,end) into the returned array.
 * This method makes the hypothesis that the field lies only on one mesh. If it is not the case an exception will be thrown.
 */
DataArrayDouble *MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.size()!=1)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt : field lies on several meshes, this method has no sense !");
  if(_field_per_mesh[0]==0)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt : no field specified !");
  return _field_per_mesh[0]->getUndergroundDataArrayExt(entries);
}

MEDFileField1TSWithoutSDA::MEDFileField1TSWithoutSDA(const char *fieldName, int csit, int fieldtype, int iteration, int order,
                                                     const std::vector<std::string>& infos):_csit(csit),_field_type(fieldtype),_iteration(iteration),_order(order)
{
  DataArrayDouble *arr=getOrCreateAndGetArray();
  arr->setName(fieldName);
  arr->setInfoAndChangeNbOfCompo(infos);
}

MEDFileField1TSWithoutSDA::MEDFileField1TSWithoutSDA():_csit(-1),_field_type(-1)
{
}

int MEDFileField1TSWithoutSDA::addNewEntryIfNecessary(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception)
{
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

/*!
 * \param [in] mName specifies the underlying mesh name. This value can be pointer 0 for users that do not deal with fields on multi mesh.
 */
int MEDFileField1TSWithoutSDA::getMeshIdFromMeshName(const char *mName) const throw(INTERP_KERNEL::Exception)
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

/*!
 * \param [in] mName specifies the underlying mesh name. This value can be pointer 0 for users that do not deal with fields on multi mesh.
 * \param [in] typ is for the geometric cell type (or INTERP_KERNEL::NORM_ERROR for node field) entry to find the right MEDFileFieldPerMeshPerTypePerDisc instance to set.
 * \param [in] locId is the localization id to find the right MEDFileFieldPerMeshPerTypePerDisc instance to set. It corresponds to the position of 
 *             \c pfls[std::distance(types.begin(),std::find(types.begin(),typ)] vector in MEDFileField1TSWithoutSDA::getFieldSplitedByType. For non gausspoints field users, the value is 0.
 */
MEDFileFieldPerMeshPerTypePerDisc *MEDFileField1TSWithoutSDA::getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) throw(INTERP_KERNEL::Exception)
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
const MEDFileFieldPerMeshPerTypePerDisc *MEDFileField1TSWithoutSDA::getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) const throw(INTERP_KERNEL::Exception)
{
  int mid=getMeshIdFromMeshName(mName);
  return _field_per_mesh[mid]->getLeafGivenTypeAndLocId(typ,locId);
}

DataArrayDouble *MEDFileField1TSWithoutSDA::getOrCreateAndGetArray()
{
  DataArrayDouble *ret=_arr;
  if(ret)
    return ret;
  _arr=DataArrayDouble::New();
  return _arr;
}

const DataArrayDouble *MEDFileField1TSWithoutSDA::getOrCreateAndGetArray() const
{
  const DataArrayDouble *ret=_arr;
  if(ret)
    return ret;
  DataArrayDouble *ret2=DataArrayDouble::New();
  const_cast<MEDFileField1TSWithoutSDA *>(this)->_arr=DataArrayDouble::New();
  return ret2;
}

MEDFileField1TS *MEDFileField1TS::New(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileField1TS(fileName,fieldName,iteration,order);
}

/*!
 * \warning this is a shallow copy constructor
 */
MEDFileField1TS *MEDFileField1TS::New(const MEDFileField1TSWithoutSDA& other, bool deepCpy)
{
  return new MEDFileField1TS(other,deepCpy);
}

MEDFileField1TS *MEDFileField1TS::New()
{
  return new MEDFileField1TS;
}

std::string MEDFileField1TS::simpleRepr() const
{
  std::ostringstream oss;
  _content->simpleRepr(0,oss,-1);
  MEDFileFieldGlobsReal::simpleRepr(oss);
  return oss.str();
}

void MEDFileField1TS::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  int nbComp=getNumberOfComponents();
  INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  for(int i=0;i<nbComp;i++)
    {
      std::string info=getInfo()[i];
      std::string c,u;
      MEDLoaderBase::splitIntoNameAndUnit(info,c,u);
      MEDLoaderBase::safeStrCpy2(c.c_str(),MED_SNAME_SIZE-1,comp+i*MED_SNAME_SIZE,_too_long_str);
      MEDLoaderBase::safeStrCpy2(u.c_str(),MED_SNAME_SIZE-1,unit+i*MED_SNAME_SIZE,_too_long_str);
    }
  if(getName().empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TS::write : MED file does not accept field with empty name !");
  MEDfieldCr(fid,getName().c_str(),MED_FLOAT64,nbComp,comp,unit,getDtUnit().c_str(),getMeshName().c_str());
  writeGlobals(fid,*this);
  _content->writeLL(fid,*this);
}

void MEDFileField1TS::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,medmod);
  writeLL(fid);
}

MEDFileField1TS::MEDFileField1TS(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
try:_content(MEDFileField1TSWithoutSDA::New(fieldName,-1,-1,iteration,order,std::vector<std::string>())),MEDFileFieldGlobsReal(fileName)
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
          std::vector<std::string> infos(ncomp);
          for(int j=0;j<ncomp;j++)
            infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+j*MED_SNAME_SIZE,MED_SNAME_SIZE,(char *)unit+j*MED_SNAME_SIZE,MED_SNAME_SIZE);
          _content->getOrCreateAndGetArray()->setInfoAndChangeNbOfCompo(infos);
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
          _content->_csit=i+1;
          _content->_field_type=MEDFileUtilities::TraduceFieldType(typcha);
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
  _content->finishLoading(fid);
  //
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

/*!
 * \warning this is a shallow copy constructor
 */
MEDFileField1TS::MEDFileField1TS(const MEDFileField1TSWithoutSDA& other, bool deepCpy)
{
  if(!deepCpy)
    {
      const MEDFileField1TSWithoutSDA *otherPtr(&other);
      otherPtr->incrRef();
      _content=const_cast<MEDFileField1TSWithoutSDA *>(otherPtr);
    }
  else
    {
      _content=new MEDFileField1TSWithoutSDA(other);
    }
}

MEDFileField1TS::MEDFileField1TS():_content(new MEDFileField1TSWithoutSDA)
{
}

/*!
 * This method returns all profiles whose name is non empty used.
 * \b WARNING If profile is used several times it will be reported \b only \b once.
 * To get non empty name profiles as time as they appear in \b this call MEDFileField1TS::getPflsReallyUsedMulti instead.
 */
std::vector<std::string> MEDFileField1TS::getPflsReallyUsed() const
{
  return _content->getPflsReallyUsed2();
}

/*!
 * This method returns all localizations whose name is non empty used.
 * \b WARNING If localization is used several times it will be reported \b only \b once.
 */
std::vector<std::string> MEDFileField1TS::getLocsReallyUsed() const
{
  return _content->getLocsReallyUsed2();
}

/*!
 * This method returns all profiles whose name is non empty used.
 * \b WARNING contrary to MEDFileField1TS::getPflsReallyUsed, if profile is used several times it will be reported as time as it appears.
 */
std::vector<std::string> MEDFileField1TS::getPflsReallyUsedMulti() const
{
  return _content->getPflsReallyUsedMulti2();
}

/*!
 * This method returns all localizations whose name is non empty used.
 * \b WARNING contrary to MEDFileField1TS::getLocsReallyUsed if localization is used several times it will be reported as time as it appears.
 */
std::vector<std::string> MEDFileField1TS::getLocsReallyUsedMulti() const
{
  return _content->getLocsReallyUsedMulti2();
}

void MEDFileField1TS::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  _content->changePflsRefsNamesGen2(mapOfModif);
}

void MEDFileField1TS::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  _content->changeLocsRefsNamesGen2(mapOfModif);
}

/*!
 * This method requests underlying file to perform the job, for mesh reading. If the current instance is not coming from a file and has been constructed from scratch
 * an exception will be thrown. In this case you should use MEDFileField1TS::getFieldOnMeshAtLevel method instead.
 * \b WARNING ! Parameter 'meshDimRelToMax' is relative from read mesh in file that can be different from the field in MED file !
 * It leads that the returned field of this method is always coherent.
 */
MEDCouplingFieldDouble *MEDFileField1TS::getFieldAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  if(getFileName2().empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TS::getFieldAtLevel : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtLevel method instead !");
  return _content->getFieldAtLevel(type,meshDimRelToMax,0,renumPol,this);
}

/*!
 * This method is close to MEDFileField1TS::getFieldAtLevel except that here the 'meshDimRelToMax' param is ignored and the maximal dimension is taken
 * automatically. If the field lies on different level and that an another level than the maximal is requested MEDFileField1TS::getFieldAtLevel
 * should be called instead.
 */
MEDCouplingFieldDouble *MEDFileField1TS::getFieldAtTopLevel(TypeOfField type, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  if(getFileName2().empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TS::getFieldAtTopLevel : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtTopLevel method instead !");
  return _content->getFieldAtTopLevel(type,0,renumPol,this);
}

/*!
 * \b WARNING, there is a main difference with the two close methods (MEDFileField1TS::getFieldAtLevel and MEDFileField1TS::getFieldOnMeshAtLevel method) !
 * Here the mesh-dimension of 'mesh' is used by this to automatically request the right geoTypes regarding 'type'.
 * If no such element fufilled the deduced dimension and 'type' an exception will be thrown.
 * It leads that the returned field of this method is always coherent.
 */
MEDCouplingFieldDouble *MEDFileField1TS::getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  return _content->getFieldOnMeshAtLevel(type,renumPol,this,mesh,0,0);
}

/*!
 * This method can be called whatever the mode of instance feeding of this (MED file or from scratch).
 * But the parameter ''meshDimRelToMax' is applyied on 'mesh' (like MEDFileField1TS::getFieldAtLevel does). \b WARNING the dim of 'this' can be different from those in 'mesh' !
 * It leads that the returned field of this method is always coherent.
 */
MEDCouplingFieldDouble *MEDFileField1TS::getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  return _content->getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,this,mesh);
}

/*!
 * This method is identical to MEDFileField1TS::getFieldAtLevel method except that meshName 'mname' should be specified.
 * This method is called "Old" because in MED3 norm a field has only one meshName attached. This method is only here for reader of MED2 files.
 * See MEDFileField1TS::getFieldAtLevel for more information.
 */
MEDCouplingFieldDouble *MEDFileField1TS::getFieldAtLevelOld(TypeOfField type, const char *mname, int meshDimRelToMax, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  if(getFileName2().empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TS::getFieldAtLevel : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtLevel method instead !");
  return _content->getFieldAtLevel(type,meshDimRelToMax,mname,renumPol,this);
}

DataArrayDouble *MEDFileField1TS::getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const throw(INTERP_KERNEL::Exception)
{
  return _content->getFieldWithProfile(type,meshDimRelToMax,mesh,pfl,this);
}

/*!
 * SBT means Sort By Type.
 * This method is the most basic method to assign field in this. Basic in sense that no renumbering is done. Underlying mesh in 'field' is globaly ignored except for type contiguous check.
 * 
 */
void MEDFileField1TS::setFieldNoProfileSBT(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception)
{
  setFileName("");
  _content->setFieldNoProfileSBT(field,*this);
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
  setFileName("");
  _content->setFieldProfile(field,mesh,meshDimRelToMax,profile,*this);
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
void MEDFileField1TS::setProfileNameOnLeaf(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const char *newPflName, bool forceRenameOnGlob) throw(INTERP_KERNEL::Exception)
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
void MEDFileField1TS::setLocNameOnLeaf(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const char *newLocName, bool forceRenameOnGlob) throw(INTERP_KERNEL::Exception)
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

int MEDFileField1TS::copyTinyInfoFrom(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception)
{
  return _content->copyTinyInfoFrom(field);
}

int MEDFileField1TS::getDimension() const
{
  return _content->getDimension();
}

int MEDFileField1TS::getIteration() const
{
  return _content->getIteration();
}

int MEDFileField1TS::getOrder() const
{
  return _content->getOrder();
}

double MEDFileField1TS::getTime(int& iteration, int& order) const
{
  return _content->getTime(iteration,order);
}

void MEDFileField1TS::setTime(int iteration, int order, double val)
{
  _content->setTime(iteration,order,val);
}

std::string MEDFileField1TS::getName() const
{
  return _content->getName();
}

void MEDFileField1TS::setName(const char *name)
{
  _content->setName(name);
}

void MEDFileField1TS::simpleRepr(int bkOffset, std::ostream& oss, int f1tsId) const
{
  _content->simpleRepr(bkOffset,oss,f1tsId);
}

const std::string& MEDFileField1TS::getDtUnit() const
{
  return _content->getDtUnit();
}

std::string MEDFileField1TS::getMeshName() const throw(INTERP_KERNEL::Exception)
{
  return _content->getMeshName();
}

void MEDFileField1TS::setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception)
{
  _content->setMeshName(newMeshName);
}

bool MEDFileField1TS::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception)
{
  return _content->changeMeshNames(modifTab);
}

int MEDFileField1TS::getMeshIteration() const throw(INTERP_KERNEL::Exception)
{
  return _content->getMeshIteration();
}

int MEDFileField1TS::getMeshOrder() const throw(INTERP_KERNEL::Exception)
{
  return _content->getMeshOrder();
}

int MEDFileField1TS::getNumberOfComponents() const
{
  return _content->getNumberOfComponents();
}

bool MEDFileField1TS::isDealingTS(int iteration, int order) const
{
  return _content->isDealingTS(iteration,order);
}

std::pair<int,int> MEDFileField1TS::getDtIt() const
{
  return _content->getDtIt();
}

void MEDFileField1TS::fillIteration(std::pair<int,int>& p) const
{
  _content->fillIteration(p);
}
void MEDFileField1TS::fillTypesOfFieldAvailable(std::vector<TypeOfField>& types) const throw(INTERP_KERNEL::Exception)
{
  _content->fillTypesOfFieldAvailable(types);
}

const std::vector<std::string>& MEDFileField1TS::getInfo() const
{
  return _content->getInfo();
}
std::vector<std::string>& MEDFileField1TS::getInfo()
{
  return _content->getInfo();
}

DataArrayDouble *MEDFileField1TS::getUndergroundDataArray() const throw(INTERP_KERNEL::Exception)
{
  return _content->getUndergroundDataArray();
}

DataArrayDouble *MEDFileField1TS::getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception)
{
  return _content->getUndergroundDataArrayExt(entries);
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileField1TS::getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) throw(INTERP_KERNEL::Exception)
{
  return _content->getLeafGivenMeshAndTypeAndLocId(mName,typ,locId);
}

const MEDFileFieldPerMeshPerTypePerDisc *MEDFileField1TS::getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) const throw(INTERP_KERNEL::Exception)
{
  return _content->getLeafGivenMeshAndTypeAndLocId(mName,typ,locId);
}

int MEDFileField1TS::getNonEmptyLevels(const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception)
{
  return _content->getNonEmptyLevels(mname,levs);
}

std::vector<TypeOfField> MEDFileField1TS::getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception)
{
  return _content->getTypesOfFieldAvailable();
}

std::vector< std::vector<std::pair<int,int> > > MEDFileField1TS::getFieldSplitedByType(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF,
                                                                                       std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  return _content->getFieldSplitedByType(mname,types,typesF,pfls,locs);
}
std::vector< std::vector<DataArrayDouble *> > MEDFileField1TS::getFieldSplitedByType2(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF,
                                                                                      std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  return _content->getFieldSplitedByType2(mname,types,typesF,pfls,locs);
}

MEDFileFieldMultiTSWithoutSDA *MEDFileFieldMultiTSWithoutSDA::New(med_idt fid, const char *fieldName, int id, int ft, const std::vector<std::string>& infos, int nbOfStep) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTSWithoutSDA(fid,fieldName,id,ft,infos,nbOfStep);
}

MEDFileFieldMultiTSWithoutSDA::MEDFileFieldMultiTSWithoutSDA():_field_type(-1)
{
}

MEDFileFieldMultiTSWithoutSDA::MEDFileFieldMultiTSWithoutSDA(const char *fieldName):_name(fieldName),_field_type(-1)
{
}

/*!
 * \param [in] fieldId field id in C mode
 */
MEDFileFieldMultiTSWithoutSDA::MEDFileFieldMultiTSWithoutSDA(med_idt fid, int fieldId) throw(INTERP_KERNEL::Exception)
try:_name("")
{
  med_field_type typcha;
  int nbstep2=-1;
  //
  int ncomp=MEDfieldnComponent(fid,fieldId+1);
  INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> dtunit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> nomMaa=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  med_bool localMesh;
  int nbOfStep;
  MEDfieldInfo(fid,fieldId+1,nomcha,nomMaa,&localMesh,&typcha,comp,unit,dtunit,&nbOfStep);
  _name=MEDLoaderBase::buildStringFromFortran(nomcha,MED_NAME_SIZE);
  _field_type=MEDFileUtilities::TraduceFieldType(typcha);
  _infos.resize(ncomp);
  for(int j=0;j<ncomp;j++)
    _infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+j*MED_SNAME_SIZE,MED_SNAME_SIZE,(char *)unit+j*MED_SNAME_SIZE,MED_SNAME_SIZE);
  //
  finishLoading(fid,nbOfStep);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

MEDFileFieldMultiTSWithoutSDA::MEDFileFieldMultiTSWithoutSDA(med_idt fid, const char *fieldName, int id, int ft, const std::vector<std::string>& infos, int nbOfStep) throw(INTERP_KERNEL::Exception)
try:_name(fieldName),_infos(infos),_field_type(ft)
{
  finishLoading(fid,nbOfStep);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

const std::vector<std::string>& MEDFileFieldMultiTSWithoutSDA::getInfo() const throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::getInfos : not time steps !");
  return _time_steps[0]->getInfo();
}

/*!
 * See doc at MEDFileField1TSWithoutSDA::getUndergroundDataArray
 */
DataArrayDouble *MEDFileFieldMultiTSWithoutSDA::getUndergroundDataArray(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  return getTimeStepEntry(iteration,order).getUndergroundDataArray();
}

/*!
 * See doc at MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt
 */
DataArrayDouble *MEDFileFieldMultiTSWithoutSDA::getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception)
{
  return getTimeStepEntry(iteration,order).getUndergroundDataArrayExt(entries);
}

std::string MEDFileFieldMultiTSWithoutSDA::getMeshName() const throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::getMeshName : not time steps !");
  return _time_steps[0]->getMeshName();
}

void MEDFileFieldMultiTSWithoutSDA::setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception)
{
  std::string oldName(getMeshName());
  std::vector< std::pair<std::string,std::string> > v(1);
  v[0].first=oldName; v[0].second=newMeshName;
  changeMeshNames(v);
}

bool MEDFileFieldMultiTSWithoutSDA::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception)
{
  bool ret=false;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileField1TSWithoutSDA *cur(*it);
      if(cur)
        ret=cur->changeMeshNames(modifTab) || ret;
    }
  return ret;
}

bool MEDFileFieldMultiTSWithoutSDA::renumberEntitiesLyingOnMesh(const char *meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N,
                                                                MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  bool ret=false;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileField1TSWithoutSDA *f1ts(*it);
      if(f1ts)
        ret=f1ts->renumberEntitiesLyingOnMesh(meshName,oldCode,newCode,renumO2N,glob) || ret;
    }
  return ret;
}

void MEDFileFieldMultiTSWithoutSDA::appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA> obj=new MEDFileField1TSWithoutSDA;
      obj->setFieldNoProfileSBT(field,glob);
      copyTinyInfoFrom(field);
      _time_steps.push_back(obj);
    }
  else
    {
      checkCoherencyOfTinyInfo(field);
      MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA> obj=new MEDFileField1TSWithoutSDA;
      obj->setFieldNoProfileSBT(field,glob);
      _time_steps.push_back(obj);
    }
}

void MEDFileFieldMultiTSWithoutSDA::appendFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    {
      MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA> obj=new MEDFileField1TSWithoutSDA;
      obj->setFieldProfile(field,mesh,meshDimRelToMax,profile,glob);
      copyTinyInfoFrom(field);
      _time_steps.push_back(obj);
    }
  else
    {
      checkCoherencyOfTinyInfo(field);
      MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA> obj=new MEDFileField1TSWithoutSDA;
      obj->setFieldProfile(field,mesh,meshDimRelToMax,profile,glob);
      _time_steps.push_back(obj);
    }
}

std::string MEDFileFieldMultiTSWithoutSDA::getDtUnit() const throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::getMeshName : not time steps !");
  return _time_steps[0]->getDtUnit();
}

std::string MEDFileFieldMultiTSWithoutSDA::getName() const
{
  return _name;
}

void MEDFileFieldMultiTSWithoutSDA::setName(const char *name)
{
  _name=name;
}

void MEDFileFieldMultiTSWithoutSDA::simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const
{
  std::string startLine(bkOffset,' ');
  oss << startLine << "Field multi time steps";
  if(fmtsId>=0)
    oss << " (" << fmtsId << ")";
  oss << " has the following name: \"" << _name << "\"." << std::endl;
  oss << startLine << "Field multi time steps has " << _infos.size() << " components with the following infos :" << std::endl;
  for(std::vector<std::string>::const_iterator it=_infos.begin();it!=_infos.end();it++)
    {
      oss << startLine << "  -  \"" << *it << "\"" << std::endl;
    }
  int i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,i++)
    {
      std::string chapter(17,'0'+i);
      oss << startLine << chapter << std::endl;
      const MEDFileField1TSWithoutSDA *cur=(*it);
      if(cur)
        cur->simpleRepr(bkOffset+2,oss,i);
      else
        oss << startLine << "  Field on one time step #" << i << " is not defined !" << std::endl;
      oss << startLine << chapter << std::endl;
    }
}

std::vector< std::pair<int,int> > MEDFileFieldMultiTSWithoutSDA::getTimeSteps(std::vector<double>& ret1) const throw(INTERP_KERNEL::Exception)
{
  std::size_t sz=_time_steps.size();
  std::vector< std::pair<int,int> > ret(sz);
  ret1.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDFileField1TSWithoutSDA *f1ts=_time_steps[i];
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

void MEDFileFieldMultiTSWithoutSDA::finishLoading(med_idt fid, int nbPdt) throw(INTERP_KERNEL::Exception)
{
  _time_steps.resize(nbPdt);
  for(int i=0;i<nbPdt;i++)
    {
      std::vector< std::pair<int,int> > ts;
      med_int numdt=0,numo=0;
      med_int meshIt=0,meshOrder=0;
      med_float dt=0.0;
      MEDfieldComputingStepMeshInfo(fid,_name.c_str(),i+1,&numdt,&numo,&dt,&meshIt,&meshOrder);
      _time_steps[i]=MEDFileField1TSWithoutSDA::New(_name.c_str(),i+1,_field_type,numdt,numo,_infos);
      _time_steps[i]->finishLoading(fid);
    }
}

void MEDFileFieldMultiTSWithoutSDA::copyTinyInfoFrom(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception)
{
  _name=field->getName();
  if(_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::copyTinyInfoFrom : unsupported fields with no name in MED file !");
  const DataArrayDouble *arr=field->getArray();
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::copyTinyInfoFrom : no array set !");
  _infos=arr->getInfoOnComponents();
}

void MEDFileFieldMultiTSWithoutSDA::checkCoherencyOfTinyInfo(const MEDCouplingFieldDouble *field) const throw(INTERP_KERNEL::Exception)
{
  static const char MSG[]="MEDFileFieldMultiTSWithoutSDA::checkCoherencyOfTinyInfo : invalid ";
  if(_name!=field->getName())
    {
      std::ostringstream oss; oss << MSG << "name ! should be \"" << _name;
      oss << "\" and it is set in input field to \"" << field->getName() << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const DataArrayDouble *arr=field->getArray();
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::checkCoherencyOfTinyInfo : no array set !");
  if(_infos!=arr->getInfoOnComponents())
    {
      std::ostringstream oss; oss << MSG << "components ! should be \"";
      std::copy(_infos.begin(),_infos.end(),std::ostream_iterator<std::string>(oss,", "));
      oss << " But compo in input fields are : ";
      std::vector<std::string> tmp=arr->getInfoOnComponents();
      std::copy(tmp.begin(),tmp.end(),std::ostream_iterator<std::string>(oss,", "));
      oss << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void MEDFileFieldMultiTSWithoutSDA::writeLL(med_idt fid, const MEDFileWritable& opts) const throw(INTERP_KERNEL::Exception)
{
  if(_time_steps.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::writeLL : no time steps set !");
  std::vector<std::string> infos(getInfo());
  int nbComp=infos.size();
  INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE);
  for(int i=0;i<nbComp;i++)
    {
      std::string info=infos[i];
      std::string c,u;
      MEDLoaderBase::splitIntoNameAndUnit(info,c,u);
      MEDLoaderBase::safeStrCpy2(c.c_str(),MED_SNAME_SIZE-1,comp+i*MED_SNAME_SIZE,opts.getTooLongStrPolicy());
      MEDLoaderBase::safeStrCpy2(u.c_str(),MED_SNAME_SIZE-1,unit+i*MED_SNAME_SIZE,opts.getTooLongStrPolicy());
    }
  if(_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::write : MED file does not accept field with empty name !");
  MEDfieldCr(fid,_name.c_str(),MED_FLOAT64,nbComp,comp,unit,getDtUnit().c_str(),getMeshName().c_str());
  int nbOfTS=_time_steps.size();
  for(int i=0;i<nbOfTS;i++)
    _time_steps[i]->writeLL(fid,opts);
}

int MEDFileFieldMultiTSWithoutSDA::getNumberOfTS() const
{
  return _time_steps.size();
}

void MEDFileFieldMultiTSWithoutSDA::eraseEmptyTS() throw(INTERP_KERNEL::Exception)
{
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA>  > newTS;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      const MEDFileField1TSWithoutSDA *tmp=(*it);
      if(tmp)
        newTS.push_back(*it);
    }
  _time_steps=newTS;
}

void MEDFileFieldMultiTSWithoutSDA::eraseTimeStepIds(const int *startIds, const int *endIds) throw(INTERP_KERNEL::Exception)
{
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA>  > newTS;
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

int MEDFileFieldMultiTSWithoutSDA::getPosOfTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  int ret=0;
  std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::getPosOfTimeStep : No such time step (" << iteration << "," << order << ") !\nPossibilities are : "; 
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,ret++)
    {
      const MEDFileField1TSWithoutSDA *tmp(*it);
      if(tmp)
        {
          int it,ord;
          tmp->getTime(it,ord);
          if(it==iteration && order==ord)
            return ret;
          else
            oss << "(" << it << ","  << ord << "), ";
        }
    }
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

int MEDFileFieldMultiTSWithoutSDA::getPosGivenTime(double time, double eps) const throw(INTERP_KERNEL::Exception)
{
  int ret=0;
  std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::getPosGivenTime : No such time step " << time << "! \nPossibilities are : ";
  oss.precision(15);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,ret++)
    {
      const MEDFileField1TSWithoutSDA *tmp(*it);
      if(tmp)
        {
          int it,ord;
          double ti=tmp->getTime(it,ord);
          if(fabs(time-ti)<eps)
            return ret;
          else
            oss << ti << ", ";
        }
    }
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

std::vector< std::pair<int,int> > MEDFileFieldMultiTSWithoutSDA::getIterations() const
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
int MEDFileFieldMultiTSWithoutSDA::getNonEmptyLevels(int iteration, int order, const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception)
{
  return getTimeStepEntry(iteration,order).getNonEmptyLevels(mname,levs);
}

std::vector< std::vector<TypeOfField> > MEDFileFieldMultiTSWithoutSDA::getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception)
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
std::vector< std::vector< std::pair<int,int> > > MEDFileFieldMultiTSWithoutSDA::getFieldSplitedByType(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  return getTimeStepEntry(iteration,order).getFieldSplitedByType(mname,types,typesF,pfls,locs);
}

/*!
 * entry point for users that want to iterate into MEDFile DataStructure with a reduced overhead because output arrays are extracted (created) specially
 * for the call of this method. That's why the DataArrayDouble instance in returned vector of vector should be dealed by the caller.
 */
std::vector< std::vector<DataArrayDouble *> > MEDFileFieldMultiTSWithoutSDA::getFieldSplitedByType2(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  return getTimeStepEntry(iteration,order).getFieldSplitedByType2(mname,types,typesF,pfls,locs);
}

const MEDFileField1TSWithoutSDA& MEDFileFieldMultiTSWithoutSDA::getTimeStepEntry(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it)->isDealingTS(iteration,order))
      return *(*it);
  std::ostringstream oss; oss << "MEDFileFieldMultiTS::getTimeStepEntry : Muli timestep field on time (" << iteration << "," << order << ") does not exist ! Available (iteration,order) are :\n";
  std::vector< std::pair<int,int> > vp=getIterations();
  for(std::vector< std::pair<int,int> >::const_iterator it2=vp.begin();it2!=vp.end();it2++)
    oss << "(" << (*it2).first << "," << (*it2).second << ") ";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

MEDFileField1TSWithoutSDA& MEDFileFieldMultiTSWithoutSDA::getTimeStepEntry(int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA>  >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it)->isDealingTS(iteration,order))
      return *(*it);
  std::ostringstream oss; oss << "MEDFileFieldMultiTS::getTimeStepEntry : Muli timestep field on time (" << iteration << "," << order << ") does not exist ! Available (iteration,order) are :\n";
  std::vector< std::pair<int,int> > vp=getIterations();
  for(std::vector< std::pair<int,int> >::const_iterator it2=vp.begin();it2!=vp.end();it2++)
    oss << "(" << (*it2).first << "," << (*it2).second << ") ";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

const MEDFileField1TSWithoutSDA *MEDFileFieldMultiTSWithoutSDA::getTimeStepAtPos2(int pos) const throw(INTERP_KERNEL::Exception)
{
  if(pos<0 || pos>=(int)_time_steps.size())
    {
      std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << " whereas should be in [0," << _time_steps.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const MEDFileField1TSWithoutSDA *item=_time_steps[pos];
  if(item==0)
    {
      std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << ", this pos id exists but the underlying Field1TS is null !";
      oss << "\nTry to use following method eraseEmptyTS !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return item;
}

MEDFileField1TSWithoutSDA *MEDFileFieldMultiTSWithoutSDA::getTimeStepAtPos2(int pos) throw(INTERP_KERNEL::Exception)
{
  if(pos<0 || pos>=(int)_time_steps.size())
    {
      std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << " whereas should be in [0," << _time_steps.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDFileField1TSWithoutSDA *item=_time_steps[pos];
  if(item==0)
    {
      std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << ", this pos id exists but the underlying Field1TS is null !";
      oss << "\nTry to use following method eraseEmptyTS !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return item;
}

std::vector<std::string> MEDFileFieldMultiTSWithoutSDA::getPflsReallyUsed2() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
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

std::vector<std::string> MEDFileFieldMultiTSWithoutSDA::getLocsReallyUsed2() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
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

std::vector<std::string> MEDFileFieldMultiTSWithoutSDA::getPflsReallyUsedMulti2() const
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsedMulti2();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

std::vector<std::string> MEDFileFieldMultiTSWithoutSDA::getLocsReallyUsedMulti2() const
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsedMulti2();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

void MEDFileFieldMultiTSWithoutSDA::changePflsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileField1TSWithoutSDA > >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    (*it)->changePflsRefsNamesGen2(mapOfModif);
}

void MEDFileFieldMultiTSWithoutSDA::changeLocsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileField1TSWithoutSDA > >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    (*it)->changeLocsRefsNamesGen2(mapOfModif);
}

MEDFileFieldMultiTS *MEDFileFieldMultiTS::New()
{
  return new MEDFileFieldMultiTS;
}

MEDFileFieldMultiTS *MEDFileFieldMultiTS::New(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTS(fileName);
}

MEDFileFieldMultiTS *MEDFileFieldMultiTS::New(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTS(fileName,fieldName);
}

MEDFileFieldMultiTS *MEDFileFieldMultiTS::New(const MEDFileFieldMultiTSWithoutSDA& other, bool deepCpy)
{
  return new MEDFileFieldMultiTS(other,deepCpy);
}

MEDFileField1TS *MEDFileFieldMultiTS::getTimeStepAtPos(int pos) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileField1TSWithoutSDA *item=_content->getTimeStepAtPos2(pos);
  MEDCouplingAutoRefCountObjectPtr<MEDFileField1TS> ret=MEDFileField1TS::New(*item,false);
  ret->shallowCpyGlobs(*this);
  ret->incrRef();
  return ret;
}

MEDFileField1TS *MEDFileFieldMultiTS::getTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  int pos=getPosOfTimeStep(iteration,order);
  return getTimeStepAtPos(pos);
}

MEDFileField1TS *MEDFileFieldMultiTS::getTimeStepGivenTime(double time, double eps) const throw(INTERP_KERNEL::Exception)
{
  int pos=getPosGivenTime(time,eps);
  return getTimeStepAtPos(pos);
}

MEDFileFieldMultiTSIterator *MEDFileFieldMultiTS::iterator() throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTSIterator(this);
}

std::string MEDFileFieldMultiTS::simpleRepr() const
{
  std::ostringstream oss;
  _content->simpleRepr(0,oss,-1);
  MEDFileFieldGlobsReal::simpleRepr(oss);
  return oss.str();
}

void MEDFileFieldMultiTS::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  writeGlobals(fid,*this);
  _content->writeLL(fid,*this);
}

void MEDFileFieldMultiTS::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,medmod);
  writeLL(fid);
}

/*!
 * Performs the job than MEDFileField1TS::getFieldAtLevel except that (iteration,order) couple should be specified !
 * If such couple does not exist an exception is thrown.
 */
MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileField1TSWithoutSDA& myF1TS=_content->getTimeStepEntry(iteration,order);
  return myF1TS.getFieldAtLevel(type,meshDimRelToMax,0,renumPol,this);
}

MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldAtTopLevel(TypeOfField type, int iteration, int order, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileField1TSWithoutSDA& myF1TS=_content->getTimeStepEntry(iteration,order);
  return myF1TS.getFieldAtTopLevel(type,0,renumPol,this);
}

/*!
 * Performs the job than MEDFileField1TS::getFieldOnMeshAtLevel except that (iteration,order) couple should be specified !
 * If such couple does not exist an exception is thrown.
 */
MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileField1TSWithoutSDA& myF1TS=_content->getTimeStepEntry(iteration,order);
  return myF1TS.getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,this,mesh);
}

/*!
 * Performs the job than MEDFileField1TS::getFieldOnMeshAtLevel except that (iteration,order) couple should be specified !
 * If such couple does not exist an exception is thrown.
 */
MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileField1TSWithoutSDA& myF1TS=_content->getTimeStepEntry(iteration,order);
  return myF1TS.getFieldOnMeshAtLevel(type,renumPol,this,mesh,0,0);
}

/*!
 * This method has a close behaviour than MEDFileFieldMultiTS::getFieldAtLevel.
 * This method is called 'old' because the user should give the mesh name he wants to use for it's field.
 * This method is useful for MED2 file format when field on different mesh was autorized.
 */
MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldAtLevelOld(TypeOfField type, const char *mname, int iteration, int order, int meshDimRelToMax, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileField1TSWithoutSDA& myF1TS=_content->getTimeStepEntry(iteration,order);
  return myF1TS.getFieldAtLevel(type,meshDimRelToMax,mname,renumPol,this);
}

DataArrayDouble *MEDFileFieldMultiTS::getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileField1TSWithoutSDA& myF1TS=_content->getTimeStepEntry(iteration,order);
  return myF1TS.getFieldWithProfile(type,meshDimRelToMax,mesh,pfl,this);
}

void MEDFileFieldMultiTS::appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception)
{
  _content->appendFieldNoProfileSBT(field,*this);
}

void MEDFileFieldMultiTS::appendFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception)
{
  _content->appendFieldProfile(field,mesh,meshDimRelToMax,profile,*this);
}

MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutSDA> MEDFileFieldMultiTS::getContent()
{
  return _content;
}

MEDFileFieldMultiTS::MEDFileFieldMultiTS():_content(new MEDFileFieldMultiTSWithoutSDA)
{
}



MEDFileFieldMultiTS::MEDFileFieldMultiTS(const char *fileName) throw(INTERP_KERNEL::Exception)
try:MEDFileFieldGlobsReal(fileName)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  int nbFields=MEDnField(fid);
  if(nbFields<1)
    {
      std::ostringstream oss; oss << "MEDFileFieldMultiTS(const char *fileName) constructor : no fields in file \"" << fileName << "\" !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  _content=new MEDFileFieldMultiTSWithoutSDA(fid,0);
  //
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

MEDFileFieldMultiTS::MEDFileFieldMultiTS(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception)
try:MEDFileFieldGlobsReal(fileName)
{
  MEDFileUtilities::CheckFileForRead(fileName);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,MED_ACC_RDONLY);
  int nbFields=MEDnField(fid);
  med_field_type typcha;
  bool found=false;
  std::vector<std::string> fns(nbFields);
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
        _content=new MEDFileFieldMultiTSWithoutSDA(fid,i);
    }
  if(!found)
    {
      std::ostringstream oss; oss << "No such field '" << fieldName << "' in file '" << fileName << "' ! Available fields are : ";
      std::copy(fns.begin(),fns.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  //
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
  {
    throw e;
  }

MEDFileFieldMultiTS::MEDFileFieldMultiTS(const MEDFileFieldMultiTSWithoutSDA& other, bool deepCpy)
{
  if(!deepCpy)
    {
      const MEDFileFieldMultiTSWithoutSDA *otherPtr(&other);
      otherPtr->incrRef();
      _content=const_cast<MEDFileFieldMultiTSWithoutSDA *>(otherPtr);
    }
  else
    {
      _content=new MEDFileFieldMultiTSWithoutSDA(other);
    }
}

std::vector<std::string> MEDFileFieldMultiTS::getPflsReallyUsed() const
{
  return _content->getPflsReallyUsed2();
}

std::vector<std::string> MEDFileFieldMultiTS::getLocsReallyUsed() const
{
  return _content->getLocsReallyUsed2();
}

std::vector<std::string> MEDFileFieldMultiTS::getPflsReallyUsedMulti() const
{
  return _content->getPflsReallyUsedMulti2();
}

std::vector<std::string> MEDFileFieldMultiTS::getLocsReallyUsedMulti() const
{
  return _content->getLocsReallyUsedMulti2();
}

void MEDFileFieldMultiTS::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  _content->changePflsRefsNamesGen2(mapOfModif);
}

void MEDFileFieldMultiTS::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  _content->changeLocsRefsNamesGen2(mapOfModif);
}

int MEDFileFieldMultiTS::getNumberOfTS() const
{
  return _content->getNumberOfTS();
}

void MEDFileFieldMultiTS::eraseEmptyTS() throw(INTERP_KERNEL::Exception)
{
  _content->eraseEmptyTS();
}

void MEDFileFieldMultiTS::eraseTimeStepIds(const int *startIds, const int *endIds) throw(INTERP_KERNEL::Exception)
{
  _content->eraseTimeStepIds(startIds,endIds);
}

std::vector< std::pair<int,int> > MEDFileFieldMultiTS::getIterations() const
{
  return _content->getIterations();
}

int MEDFileFieldMultiTS::getPosOfTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  return _content->getPosOfTimeStep(iteration,order);
}

int MEDFileFieldMultiTS::getPosGivenTime(double time, double eps) const throw(INTERP_KERNEL::Exception)
{
  return _content->getPosGivenTime(time,eps);
}

int MEDFileFieldMultiTS::getNonEmptyLevels(int iteration, int order, const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception)
{
  return _content->getNonEmptyLevels(iteration,order,mname,levs);
}

std::vector< std::vector<TypeOfField> > MEDFileFieldMultiTS::getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception)
{
  return _content->getTypesOfFieldAvailable();
}

std::vector< std::vector< std::pair<int,int> > > MEDFileFieldMultiTS::getFieldSplitedByType(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  return _content->getFieldSplitedByType(iteration,order,mname,types,typesF,pfls,locs);
}

std::vector< std::vector<DataArrayDouble *> > MEDFileFieldMultiTS::getFieldSplitedByType2(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  return _content->getFieldSplitedByType2(iteration,order,mname,types,typesF,pfls,locs);
}

std::string MEDFileFieldMultiTS::getName() const
{
  return _content->getName();
}

void MEDFileFieldMultiTS::setName(const char *name)
{
  _content->setName(name);
}

void MEDFileFieldMultiTS::simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const
{
  _content->simpleRepr(bkOffset,oss,fmtsId);
}

std::vector< std::pair<int,int> > MEDFileFieldMultiTS::getTimeSteps(std::vector<double>& ret1) const throw(INTERP_KERNEL::Exception)
{
  return _content->getTimeSteps(ret1);
}

std::string MEDFileFieldMultiTS::getMeshName() const throw(INTERP_KERNEL::Exception)
{
  return _content->getMeshName();
}

void MEDFileFieldMultiTS::setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception)
{
  _content->setMeshName(newMeshName);
}

bool MEDFileFieldMultiTS::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception)
{
  return _content->changeMeshNames(modifTab);
}

const std::vector<std::string>& MEDFileFieldMultiTS::getInfo() const throw(INTERP_KERNEL::Exception)
{
  return _content->getInfo();
}

DataArrayDouble *MEDFileFieldMultiTS::getUndergroundDataArray(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  return _content->getUndergroundDataArray(iteration,order);
}

DataArrayDouble *MEDFileFieldMultiTS::getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception)
{
  return _content->getUndergroundDataArrayExt(iteration,order,entries);
}

MEDFileFieldMultiTSIterator::MEDFileFieldMultiTSIterator(MEDFileFieldMultiTS *fmts):_fmts(fmts),_iter_id(0),_nb_iter(0)
{
  if(fmts)
    {
      fmts->incrRef();
      _nb_iter=fmts->getNumberOfTS();
    }
}

MEDFileFieldMultiTSIterator::~MEDFileFieldMultiTSIterator() 
{
}

MEDFileField1TS *MEDFileFieldMultiTSIterator::nextt()
{
  if(_iter_id<_nb_iter)
    {
      MEDFileFieldMultiTS *fmts(_fmts);
      if(fmts)
        return fmts->getTimeStepAtPos(_iter_id++);
      else
        return 0;
    }
  else
    return 0;
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

std::vector<std::string> MEDFileFields::getFieldsNames() const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> ret(_fields.size());
  int i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      const MEDFileFieldMultiTSWithoutSDA *f=(*it);
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
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      const MEDFileFieldMultiTSWithoutSDA *cur(*it);
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
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      const MEDFileFieldMultiTSWithoutSDA *cur=(*it);
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
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      const MEDFileFieldMultiTSWithoutSDA *cur=(*it);
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
  MEDFileFieldGlobsReal::simpleRepr(oss);
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
        int ncomp=MEDfieldnComponent(fid,i+1);
        INTERP_KERNEL::AutoPtr<char> comp=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
        INTERP_KERNEL::AutoPtr<char> unit=MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE);
        INTERP_KERNEL::AutoPtr<char> dtunit=MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE);
        INTERP_KERNEL::AutoPtr<char> nomcha=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
        INTERP_KERNEL::AutoPtr<char> nomMaa=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
        med_bool localMesh;
        int nbOfStep;
        MEDfieldInfo(fid,i+1,nomcha,nomMaa,&localMesh,&typcha,comp,unit,dtunit,&nbOfStep);
        int ft=MEDFileUtilities::TraduceFieldType(typcha);
        std::vector<std::string> infos(ncomp);
        for(int j=0;j<ncomp;j++)
          infos[j]=MEDLoaderBase::buildUnionUnit((char *)comp+j*MED_SNAME_SIZE,MED_SNAME_SIZE,(char *)unit+j*MED_SNAME_SIZE,MED_SNAME_SIZE);
        _fields[i]=MEDFileFieldMultiTSWithoutSDA::New(fid,nomcha,i+1,ft,infos,nbOfStep);
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
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutSDA> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      const MEDFileFieldMultiTSWithoutSDA *elt=*it;
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
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldMultiTSWithoutSDA > >::const_iterator it=_fields.begin();it!=_fields.end();it++)
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
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldMultiTSWithoutSDA > >::const_iterator it=_fields.begin();it!=_fields.end();it++)
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
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldMultiTSWithoutSDA > >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsedMulti2();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

std::vector<std::string> MEDFileFields::getLocsReallyUsedMulti() const
{
  std::vector<std::string> ret;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldMultiTSWithoutSDA > >::const_iterator it=_fields.begin();it!=_fields.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsed2();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

void MEDFileFields::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldMultiTSWithoutSDA > >::iterator it=_fields.begin();it!=_fields.end();it++)
    (*it)->changePflsRefsNamesGen2(mapOfModif);
}

void MEDFileFields::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception)
{
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldMultiTSWithoutSDA > >::iterator it=_fields.begin();it!=_fields.end();it++)
    (*it)->changeLocsRefsNamesGen2(mapOfModif);
}

void MEDFileFields::resize(int newSize) throw(INTERP_KERNEL::Exception)
{
  _fields.resize(newSize);
}

void MEDFileFields::pushField(MEDFileFieldMultiTS *field) throw(INTERP_KERNEL::Exception)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileFields::pushMesh : invalid input pointer ! should be different from 0 !");
  _fields.push_back(field->getContent());
  appendGlobs(*field,1e-12);
}

void MEDFileFields::setFieldAtPos(int i, MEDFileFieldMultiTS *field) throw(INTERP_KERNEL::Exception)
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
  if(i<0 || i>=(int)_fields.size())
    {
      std::ostringstream oss; oss << "MEDFileFields::destroyMeshAtPos : Invalid given id in input (" << i << ") should be in [0," << _fields.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  _fields.erase(_fields.begin()+i);
}

bool MEDFileFields::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception)
{
  bool ret=false;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutSDA> >::iterator it=_fields.begin();it!=_fields.end();it++)
    {
      MEDFileFieldMultiTSWithoutSDA *cur(*it);
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
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutSDA> >::iterator it=_fields.begin();it!=_fields.end();it++)
    {
      MEDFileFieldMultiTSWithoutSDA *fmts(*it);
      if(fmts)
        {
          ret=fmts->renumberEntitiesLyingOnMesh(meshName,oldCode,newCode,renumO2N,*this) || ret;
        }
    }
  return ret;
}

MEDFileFieldMultiTS *MEDFileFields::getFieldAtPos(int i) const throw(INTERP_KERNEL::Exception)
{
  if(i<0 || i>=(int)_fields.size())
    {
      std::ostringstream oss; oss << "MEDFileFields::getFieldAtPos : Invalid given id in input (" << i << ") should be in [0," << _fields.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const MEDFileFieldMultiTSWithoutSDA *fmts=_fields[i];
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTS> ret=MEDFileFieldMultiTS::New(*fmts,false);
  ret->shallowCpyGlobs(*this);
  ret->incrRef();
  return ret;
}

MEDFileFieldMultiTS *MEDFileFields::getFieldWithName(const char *fieldName) const throw(INTERP_KERNEL::Exception)
{
  return getFieldAtPos(getPosFromFieldName(fieldName));
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
      const MEDFileFieldMultiTSWithoutSDA *f=_fields[i];
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

MEDFileFieldMultiTS *MEDFileFieldsIterator::nextt()
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
