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

void MEDFileFieldPerMeshPerTypePerDisc::assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
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
void MEDFileFieldPerMeshPerTypePerDisc::assignFieldProfile(int& start, const char *pflName, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
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

void MEDFileFieldPerMeshPerTypePerDisc::assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
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

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, TypeOfField atype, int profileIt) throw(INTERP_KERNEL::Exception)
try:_type(atype),_father(fath),_profile_it(profileIt)
  {
  }
catch(INTERP_KERNEL::Exception& e)
{
  throw e;
}

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int locId, const std::string& dummy):_type(type),_father(fath),_loc_id(locId)
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

std::string MEDFileFieldPerMeshPerTypePerDisc::getLocalization() const
{
  return _localization;
}

void MEDFileFieldPerMeshPerTypePerDisc::getFieldAtLevel(TypeOfField type, const MEDFieldFieldGlobsReal *glob, std::vector< std::pair<int,int> >& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs, std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const
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

MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerType::NewOnRead(med_idt fid, MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMeshPerType(fid,fath,type,geoType);
}

MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerType::New(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMeshPerType(fath,geoType);
}

void MEDFileFieldPerMeshPerType::assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  std::vector<int> pos=addNewEntryIfNecessary(field,offset,nbOfCells);
  for(std::vector<int>::const_iterator it=pos.begin();it!=pos.end();it++)
    _field_pm_pt_pd[*it]->assignFieldNoProfile(start,offset,nbOfCells,field,glob);
}

void MEDFileFieldPerMeshPerType::assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, DataArrayInt *locIds, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
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

void MEDFileFieldPerMeshPerType::assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  _field_pm_pt_pd.resize(1);
  _field_pm_pt_pd[0]=MEDFileFieldPerMeshPerTypePerDisc::New(this,ON_NODES,-3);
  _field_pm_pt_pd[0]->assignNodeFieldNoProfile(start,field,glob);
}

void MEDFileFieldPerMeshPerType::assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
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

void MEDFileFieldPerMeshPerType::getFieldAtLevel(int meshDim, TypeOfField type, const MEDFieldFieldGlobsReal *glob, std::vector< std::pair<int,int> >& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs, std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const
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

MEDFileFieldPerMesh *MEDFileFieldPerMesh::NewOnRead(med_idt fid, MEDFileField1TSWithoutDAS *fath, int meshCsit, int meshIteration, int meshOrder) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldPerMesh(fid,fath,meshCsit,meshIteration,meshOrder);
}

MEDFileFieldPerMesh *MEDFileFieldPerMesh::New(MEDFileField1TSWithoutDAS *fath, const MEDCouplingMesh *mesh)
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

void MEDFileFieldPerMesh::assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const std::vector<int>& code, const std::vector<DataArrayInt *>& idsInPflPerType, const std::vector<DataArrayInt *>& idsPerType, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
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

void MEDFileFieldPerMesh::assignFieldNoProfileNoRenum(int& start, const std::vector<int>& code, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
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
void MEDFileFieldPerMesh::assignFieldProfileGeneral(int& start, const DataArrayInt *multiTypePfl, const std::vector<int>& code, const std::vector<DataArrayInt *>& idsInPflPerType, const std::vector<DataArrayInt *>& idsPerType, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
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

void MEDFileFieldPerMesh::assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  int pos=addNewEntryIfNecessary(INTERP_KERNEL::NORM_ERROR);
  _field_pm_pt[pos]->assignNodeFieldNoProfile(start,field,glob);
}

void MEDFileFieldPerMesh::assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
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
  const MEDFileField1TSWithoutDAS *fath=_father;
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
void MEDFileFieldPerMesh::SortArraysPerType(const MEDFieldFieldGlobsReal *glob, TypeOfField type, const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector< std::pair<int,int> >& dads, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& locs, std::vector<int>& code, std::vector<DataArrayInt *>& notNullPfls)
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
int MEDFileFieldPerMesh::ComputeNbOfElems(const MEDFieldFieldGlobsReal *glob, TypeOfField type, const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs) throw(INTERP_KERNEL::Exception)
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

MEDCouplingFieldDouble *MEDFileFieldPerMesh::getFieldOnMeshAtLevel(TypeOfField type, const MEDFieldFieldGlobsReal *glob, const MEDCouplingMesh *mesh, bool& isPfl) const throw(INTERP_KERNEL::Exception)
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

DataArrayDouble *MEDFileFieldPerMesh::getFieldOnMeshAtLevelWithPfl(TypeOfField type, const MEDCouplingMesh *mesh, DataArrayInt *&pfl, const MEDFieldFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception)
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
MEDCouplingFieldDouble *MEDFileFieldPerMesh::finishField(TypeOfField type, const MEDFieldFieldGlobsReal *glob,
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
MEDCouplingFieldDouble *MEDFileFieldPerMesh::finishField2(TypeOfField type, const MEDFieldFieldGlobsReal *glob,
                                                          const std::vector<std::pair<int,int> >& dads, const std::vector<int>& locs,
                                                          const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes,
                                                          const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl) const throw(INTERP_KERNEL::Exception)
{
  if(da->isIdentity())
    {
      int nbOfTuples=da->getNumberOfTuples();
      if(nbOfTuples==ComputeNbOfElems(glob,type,geoTypes,dads,locs))
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
MEDCouplingFieldDouble *MEDFileFieldPerMesh::finishField3(const MEDFieldFieldGlobsReal *glob,
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

MEDFileFieldPerMesh::MEDFileFieldPerMesh(med_idt fid, MEDFileField1TSWithoutDAS *fath, int meshCsit, int meshIteration, int meshOrder) throw(INTERP_KERNEL::Exception):_mesh_iteration(meshIteration),_mesh_order(meshOrder),
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

MEDFileFieldPerMesh::MEDFileFieldPerMesh(MEDFileField1TSWithoutDAS *fath, const MEDCouplingMesh *mesh):_father(fath)
{
  copyTinyInfoFrom(mesh);
}

void MEDFieldFieldGlobs::loadProfileInFile(med_idt fid, int id, const char *pflName) throw(INTERP_KERNEL::Exception)
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

void MEDFieldFieldGlobs::loadProfileInFile(med_idt fid, int i)
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

void MEDFieldFieldGlobs::appendGlobs(const MEDFieldFieldGlobs& other, double eps) throw(INTERP_KERNEL::Exception)
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
              std::ostringstream oss; oss << "MEDFieldFieldGlobs::appendGlobs : Profile \"" << (*it)->getName() << "\" already exists and is different from those expecting to be append !";
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
              std::ostringstream oss; oss << "MEDFieldFieldGlobs::appendGlobs : Localization \"" << (*it)->getName() << "\" already exists and is different from those expecting to be append !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
    }
}

void MEDFieldFieldGlobs::loadGlobals(med_idt fid, const MEDFieldFieldGlobsReal& real) throw(INTERP_KERNEL::Exception)
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

void MEDFieldFieldGlobs::loadAllGlobals(med_idt fid) throw(INTERP_KERNEL::Exception)
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

MEDFieldFieldGlobs *MEDFieldFieldGlobs::New(const char *fname)
{
  return new MEDFieldFieldGlobs(fname);
}

MEDFieldFieldGlobs *MEDFieldFieldGlobs::New()
{
  return new MEDFieldFieldGlobs;
}

MEDFieldFieldGlobs::MEDFieldFieldGlobs(const char *fname):_file_name(fname)
{
}

MEDFieldFieldGlobs::MEDFieldFieldGlobs()
{
}

MEDFieldFieldGlobs::~MEDFieldFieldGlobs()
{
}

void MEDFieldFieldGlobs::simpleRepr(std::ostream& oss) const
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

const MEDFileFieldLoc& MEDFieldFieldGlobs::getLocalization(const char *pflName) const throw(INTERP_KERNEL::Exception)
{
  return getLocalizationFromId(getLocalizationId(pflName));
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

const DataArrayInt *MEDFieldFieldGlobs::getProfile(const char *pflName) const throw(INTERP_KERNEL::Exception)
{
  std::string pflNameCpp(pflName);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> >::const_iterator it=std::find_if(_pfls.begin(),_pfls.end(),ParaMEDMEMImpl::PflFinder(pflNameCpp));
  if(it==_pfls.end())
    {
      std::ostringstream oss; oss << "MEDFieldFieldGlobs::getProfile: no such profile name : \"" << pflNameCpp << "\" Possible profiles are : ";
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

MEDFieldFieldGlobsReal::MEDFieldFieldGlobsReal(const char *fname):_globals(MEDFieldFieldGlobs::New(fname))
{
}

MEDFieldFieldGlobsReal::MEDFieldFieldGlobsReal():_globals(MEDFieldFieldGlobs::New())
{
}

void MEDFieldFieldGlobsReal::simpleRepr(std::ostream& oss) const
{
  oss << "Globals information on fields :" << "\n*******************************\n\n";
  const MEDFieldFieldGlobs *glob=_globals;
  if(glob)
    glob->simpleRepr(oss);
  else
    oss << "NO GLOBAL INFORMATION !\n";
}

MEDFieldFieldGlobsReal::~MEDFieldFieldGlobsReal()
{
}

void MEDFieldFieldGlobsReal::shallowCpyGlobs(const MEDFieldFieldGlobsReal& other)
{
  _globals=other._globals;
}

void MEDFieldFieldGlobsReal::appendGlobs(const MEDFieldFieldGlobsReal& other, double eps) throw(INTERP_KERNEL::Exception)
{
  _globals->appendGlobs(*other._globals,eps);
}

void MEDFieldFieldGlobsReal::loadProfileInFile(med_idt fid, int id, const char *pflName) throw(INTERP_KERNEL::Exception)
{
  _globals->loadProfileInFile(fid,id,pflName);
}

void MEDFieldFieldGlobsReal::loadProfileInFile(med_idt fid, int id)
{
  _globals->loadProfileInFile(fid,id);
}

void MEDFieldFieldGlobsReal::loadGlobals(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  _globals->loadGlobals(fid,*this);
}

void MEDFieldFieldGlobsReal::loadAllGlobals(med_idt fid) throw(INTERP_KERNEL::Exception)
{
  _globals->loadAllGlobals(fid);
}

void MEDFieldFieldGlobsReal::writeGlobals(med_idt fid, const MEDFileWritable& opt) const throw(INTERP_KERNEL::Exception)
{
  _globals->writeGlobals(fid,opt);
}

std::vector<std::string> MEDFieldFieldGlobsReal::getPfls() const
{
  return _globals->getPfls();
}

std::vector<std::string> MEDFieldFieldGlobsReal::getLocs() const
{
  return _globals->getLocs();
}

void MEDFieldFieldGlobsReal::setFileName(const char *fileName)
{
  _globals->setFileName(fileName);
}

int MEDFieldFieldGlobsReal::getNbOfGaussPtPerCell(int locId) const throw(INTERP_KERNEL::Exception)
{
  return _globals->getNbOfGaussPtPerCell(locId);
}

int MEDFieldFieldGlobsReal::getLocalizationId(const char *loc) const throw(INTERP_KERNEL::Exception)
{
  return _globals->getLocalizationId(loc);
}

const char *MEDFieldFieldGlobsReal::getFileName() const
{
  return _globals->getFileName();
}

std::string MEDFieldFieldGlobsReal::getFileName2() const
{
  return _globals->getFileName2();
}

const MEDFileFieldLoc& MEDFieldFieldGlobsReal::getLocalization(const char *pflName) const throw(INTERP_KERNEL::Exception)
{
  return _globals->getLocalization(pflName);
}

const MEDFileFieldLoc& MEDFieldFieldGlobsReal::getLocalizationFromId(int locId) const throw(INTERP_KERNEL::Exception)
{
  return _globals->getLocalizationFromId(locId);
}

const DataArrayInt *MEDFieldFieldGlobsReal::getProfile(const char *pflName) const throw(INTERP_KERNEL::Exception)
{
  return _globals->getProfile(pflName);
}

void MEDFieldFieldGlobsReal::appendProfile(DataArrayInt *pfl) throw(INTERP_KERNEL::Exception)
{
  _globals->appendProfile(pfl);
}

void MEDFieldFieldGlobsReal::appendLoc(const char *locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w) throw(INTERP_KERNEL::Exception)
{
  _globals->appendLoc(locName,geoType,refCoo,gsCoo,w);
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

MEDFileField1TSWithoutDAS *MEDFileField1TSWithoutDAS::New(const char *fieldName, int csit, int fieldtype, int iteration, int order, const std::vector<std::string>& infos)
{
  return new MEDFileField1TSWithoutDAS(fieldName,csit,fieldtype,iteration,order,infos);
}

/*!
 * This method copyies tiny info but also preallocated the DataArrayDouble instance in this->_arr.
 * This not allocated it allocates to the size of 'field' array. If already allocated it grows the array to
 * the previous size + the size of the array of the input 'field'.
 * This method returns the position (in tuple id) where to start to feed 'this->_arr'
 */
int MEDFileField1TSWithoutDAS::copyTinyInfoFrom(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception)
{
  std::string name(field->getName());
  getOrCreateAndGetArray()->setName(name.c_str());
  if(name.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::copyTinyInfoFrom : unsupported fields with no name in MED file !");
  const DataArrayDouble *arr=field->getArray();
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::copyTinyInfoFrom : no array set !");
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

std::string MEDFileField1TSWithoutDAS::getName() const
{
  const DataArrayDouble *arr=getOrCreateAndGetArray();
  return arr->getName();
}

void MEDFileField1TSWithoutDAS::simpleRepr(int bkOffset, std::ostream& oss, int f1tsId) const
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

int MEDFileField1TSWithoutDAS::getNumberOfComponents() const
{
  return getOrCreateAndGetArray()->getNumberOfComponents();
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

void MEDFileField1TSWithoutDAS::fillTypesOfFieldAvailable(std::vector<TypeOfField>& types) const throw(INTERP_KERNEL::Exception)
{
  std::set<TypeOfField> types2;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    {
      (*it)->fillTypesOfFieldAvailable(types2);
    }
  std::back_insert_iterator< std::vector<TypeOfField> > bi(types);
  std::copy(types2.begin(),types2.end(),bi);
}

const std::vector<std::string>& MEDFileField1TSWithoutDAS::getInfo() const
{
  const DataArrayDouble *arr=getOrCreateAndGetArray();
  return arr->getInfoOnComponents();
}

std::vector<std::string>& MEDFileField1TSWithoutDAS::getInfo()
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
int MEDFileField1TSWithoutDAS::getNonEmptyLevels(const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception)
{
  levs.clear();
  int meshId=getMeshIdFromMeshName(mname);
  std::vector<INTERP_KERNEL::NormalizedCellType> types;
  std::vector< std::vector<TypeOfField> > typesF;
  std::vector< std::vector<std::string> > pfls, locs;
  _field_per_mesh[meshId]->getFieldSplitedByType(types,typesF,pfls,locs);
  if(types.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::getNonEmptyLevels : 'this' is empty !");
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

std::vector<TypeOfField> MEDFileField1TSWithoutDAS::getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception)
{
  std::vector<TypeOfField> ret;
  fillTypesOfFieldAvailable(ret);
  return ret;
}

/*!
 * entry point for users that want to iterate into MEDFile DataStructure without any overhead.
 */
std::vector< std::vector< std::pair<int,int> > > MEDFileField1TSWithoutDAS::getFieldSplitedByType(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  int meshId=0;
  if(mname)
    meshId=getMeshIdFromMeshName(mname);
  else
    if(_field_per_mesh.empty())
      throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::getFieldSplitedByType : This is empty !");
  return _field_per_mesh[meshId]->getFieldSplitedByType(types,typesF,pfls,locs);
}

/*!
 * entry point for users that want to iterate into MEDFile DataStructure with a reduced overhead because output arrays are extracted (created) specially
 * for the call of this method. That's why the DataArrayDouble instance in returned vector of vector should be dealed by the caller.
 */
std::vector< std::vector<DataArrayDouble *> > MEDFileField1TSWithoutDAS::getFieldSplitedByType2(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  int meshId=0;
  if(mname)
    meshId=getMeshIdFromMeshName(mname);
  else
    if(_field_per_mesh.empty())
      throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::getFieldSplitedByType : This is empty !");
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

void MEDFileField1TSWithoutDAS::finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception)
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
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::finishLoading : unexpected exception internal error !");
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
void MEDFileField1TSWithoutDAS::setFieldNoProfileSBT(const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
{
  const MEDCouplingMesh *mesh=field->getMesh();
  //
  TypeOfField type=field->getTypeOfField();
  std::vector<DataArrayInt *> dummy;
  int start=copyTinyInfoFrom(field);
  if(type!=ON_NODES)
    {
      std::vector<int> code=MEDFileField1TSWithoutDAS::CheckSBTMesh(mesh);
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
 * Generalization of MEDFileField1TSWithoutDAS::setFieldNoProfileSBT method.
 */
void MEDFileField1TSWithoutDAS::setFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception)
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

MEDCouplingFieldDouble *MEDFileField1TSWithoutDAS::getFieldAtLevel(TypeOfField type, int meshDimRelToMax, const char *mName, int renumPol, const MEDFieldFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm;
  if(mName==0)
    mm=MEDFileMesh::New(glob->getFileName(),getMeshName().c_str(),getMeshIteration(),getMeshOrder());
  else
    mm=MEDFileMesh::New(glob->getFileName(),mName,-1,-1);
  return MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,glob,mm);
}

MEDCouplingFieldDouble *MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol, const MEDFieldFieldGlobsReal *glob, const MEDFileMesh *mesh) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> m=mesh->getGenMeshAtLevel(meshDimRelToMax,false);
  const DataArrayInt *d=mesh->getNumberFieldAtLevel(meshDimRelToMax);
  const DataArrayInt *e=mesh->getNumberFieldAtLevel(1);
  if(meshDimRelToMax==1)
    (static_cast<MEDCouplingUMesh *>((MEDCouplingMesh *)m))->setMeshDimension(0);
  return MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel(type,renumPol,glob,m,d,e);
}

MEDCouplingFieldDouble *MEDFileField1TSWithoutDAS::getFieldAtTopLevel(TypeOfField type, const char *mName, int renumPol, const MEDFieldFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception)
{
   MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> mm;
  if(mName==0)
    mm=MEDFileMesh::New(glob->getFileName(),getMeshName().c_str(),getMeshIteration(),getMeshOrder());
  else
    mm=MEDFileMesh::New(glob->getFileName(),mName,-1,-1);
  int absDim=getDimension();
  int meshDimRelToMax=absDim-mm->getMeshDimension();
  return MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,glob,mm);
}

MEDCouplingFieldDouble *MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel(TypeOfField type, int renumPol, const MEDFieldFieldGlobsReal *glob, const MEDCouplingMesh *mesh, const DataArrayInt *cellRenum, const DataArrayInt *nodeRenum) const throw(INTERP_KERNEL::Exception)
{
  static const char msg1[]="MEDFileField1TSWithoutDAS::getFieldOnMeshAtLevel : request for a renumbered field following mesh numbering whereas it is a profile field !";
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

DataArrayDouble *MEDFileField1TSWithoutDAS::getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl, const MEDFieldFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception)
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
DataArrayDouble *MEDFileField1TSWithoutDAS::getUndergroundDataArray() const throw(INTERP_KERNEL::Exception)
{
  const DataArrayDouble *ret=_arr;
  if(ret)
    return const_cast<DataArrayDouble *>(ret);
  else
    return 0;
}

/*!
 * This method returns an array that the caller must deal with (contrary to those returned by MEDFileField1TSWithoutDAS::getUndergroundDataArray method).
 * The returned array is the result of the aggregation of all sub arrays stored in the MED file. So to allow the caller to select the output param
 * 'entries' is returned. This output param is a vector of a pair of 2 pairs. The first pair of pair informs of the geometric type it refers to and the discretization
 * id attached to it. The second pair of pair precise the range [begin,end) into the returned array.
 * This method makes the hypothesis that the field lies only on one mesh. If it is not the case an exception will be thrown.
 */
DataArrayDouble *MEDFileField1TSWithoutDAS::getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception)
{
  if(_field_per_mesh.size()!=1)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::getUndergroundDataArrayExt : field lies on several meshes, this method has no sense !");
  if(_field_per_mesh[0]==0)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutDAS::getUndergroundDataArrayExt : no field specified !");
  return _field_per_mesh[0]->getUndergroundDataArrayExt(entries);
}

MEDFileField1TSWithoutDAS::MEDFileField1TSWithoutDAS(const char *fieldName, int csit, int fieldtype, int iteration, int order,
                                                     const std::vector<std::string>& infos):_csit(csit),_field_type(fieldtype),_iteration(iteration),_order(order)
{
  DataArrayDouble *arr=getOrCreateAndGetArray();
  arr->setName(fieldName);
  arr->setInfoAndChangeNbOfCompo(infos);
}

MEDFileField1TSWithoutDAS::MEDFileField1TSWithoutDAS():_csit(-1),_field_type(-1)
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
  std::ostringstream oss; oss << "MEDFileField1TSWithoutDAS::getMeshIdFromMeshName : No such mesh \"" << mName2 << "\" as underlying mesh of field \"" << getName() << "\" !\n";
  oss << "Possible meshes are : ";
  for(std::vector<std::string>::const_iterator it2=msg.begin();it2!=msg.end();it2++)
    oss << "\"" << (*it2) << "\" ";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

DataArrayDouble *MEDFileField1TSWithoutDAS::getOrCreateAndGetArray()
{
  DataArrayDouble *ret=_arr;
  if(ret)
    return ret;
  _arr=DataArrayDouble::New();
  return _arr;
}

const DataArrayDouble *MEDFileField1TSWithoutDAS::getOrCreateAndGetArray() const
{
  const DataArrayDouble *ret=_arr;
  if(ret)
    return ret;
  DataArrayDouble *ret2=DataArrayDouble::New();
  const_cast<MEDFileField1TSWithoutDAS *>(this)->_arr=DataArrayDouble::New();
  return ret2;
}

MEDFileField1TS *MEDFileField1TS::New(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileField1TS(fileName,fieldName,iteration,order);
}

MEDFileField1TS *MEDFileField1TS::New()
{
  return new MEDFileField1TS;
}

std::string MEDFileField1TS::simpleRepr() const
{
  std::ostringstream oss;
  MEDFileField1TSWithoutDAS::simpleRepr(0,oss,-1);
  MEDFieldFieldGlobsReal::simpleRepr(oss);
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
  MEDFileField1TSWithoutDAS::writeLL(fid);
}

void MEDFileField1TS::write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception)
{
  med_access_mode medmod=MEDFileUtilities::TraduceWriteMode(mode);
  MEDFileUtilities::AutoFid fid=MEDfileOpen(fileName,medmod);
  writeLL(fid);
}

MEDFileField1TS::MEDFileField1TS(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception)
try:MEDFileField1TSWithoutDAS(fieldName,-1,-1,iteration,order,std::vector<std::string>()),MEDFieldFieldGlobsReal(fileName)
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
          getOrCreateAndGetArray()->setInfoAndChangeNbOfCompo(infos);
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
          _field_type=MEDFileUtilities::TraduceFieldType(typcha);
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
  if(getFileName2().empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TS::getFieldAtLevel : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtLevel method instead !");
  return MEDFileField1TSWithoutDAS::getFieldAtLevel(type,meshDimRelToMax,0,renumPol,this);
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
  return MEDFileField1TSWithoutDAS::getFieldAtTopLevel(type,0,renumPol,this);
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
  if(getFileName2().empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TS::getFieldAtLevel : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtLevel method instead !");
  return MEDFileField1TSWithoutDAS::getFieldAtLevel(type,meshDimRelToMax,mname,renumPol,this);
}

DataArrayDouble *MEDFileField1TS::getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const throw(INTERP_KERNEL::Exception)
{
  return MEDFileField1TSWithoutDAS::getFieldWithProfile(type,meshDimRelToMax,mesh,pfl,this);
}

/*!
 * SBT means Sort By Type.
 * This method is the most basic method to assign field in this. Basic in sense that no renumbering is done. Underlying mesh in 'field' is globaly ignored except for type contiguous check.
 * 
 */
void MEDFileField1TS::setFieldNoProfileSBT(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception)
{
  setFileName("");
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
  setFileName("");
  MEDFileField1TSWithoutDAS::setFieldProfile(field,mesh,meshDimRelToMax,profile,*this);
}

MEDFileFieldMultiTSWithoutDAS *MEDFileFieldMultiTSWithoutDAS::New(med_idt fid, const char *fieldName, int id, int ft, const std::vector<std::string>& infos, int nbOfStep) throw(INTERP_KERNEL::Exception)
{
  return new MEDFileFieldMultiTSWithoutDAS(fid,fieldName,id,ft,infos,nbOfStep);
}

MEDFileFieldMultiTSWithoutDAS::MEDFileFieldMultiTSWithoutDAS():_field_type(-1)
{
}

MEDFileFieldMultiTSWithoutDAS::MEDFileFieldMultiTSWithoutDAS(const char *fieldName):_name(fieldName),_field_type(-1)
{
}

MEDFileFieldMultiTSWithoutDAS::MEDFileFieldMultiTSWithoutDAS(med_idt fid, const char *fieldName, int id, int ft, const std::vector<std::string>& infos, int nbOfStep) throw(INTERP_KERNEL::Exception)
try:_name(fieldName),_infos(infos),_field_type(ft)
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

/*!
 * See doc at MEDFileField1TSWithoutDAS::getUndergroundDataArray
 */
DataArrayDouble *MEDFileFieldMultiTSWithoutDAS::getUndergroundDataArray(int iteration, int order) const throw(INTERP_KERNEL::Exception)
{
  return getTimeStepEntry(iteration,order).getUndergroundDataArray();
}

/*!
 * See doc at MEDFileField1TSWithoutDAS::getUndergroundDataArrayExt
 */
DataArrayDouble *MEDFileFieldMultiTSWithoutDAS::getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception)
{
  return getTimeStepEntry(iteration,order).getUndergroundDataArrayExt(entries);
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

std::string MEDFileFieldMultiTSWithoutDAS::getName() const
{
  return _name;
}

void MEDFileFieldMultiTSWithoutDAS::simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const
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
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutDAS>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,i++)
    {
      std::string chapter(17,'0'+i);
      oss << startLine << chapter << std::endl;
      const MEDFileField1TSWithoutDAS *cur=(*it);
      if(cur)
        cur->simpleRepr(bkOffset+2,oss,i);
      else
        oss << startLine << "  Field on one time step #" << i << " is not defined !" << std::endl;
      oss << startLine << chapter << std::endl;
    }
}

std::vector< std::pair<int,int> > MEDFileFieldMultiTSWithoutDAS::getTimeSteps(std::vector<double>& ret1) const throw(INTERP_KERNEL::Exception)
{
  std::size_t sz=_time_steps.size();
  std::vector< std::pair<int,int> > ret(sz);
  ret1.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDFileField1TSWithoutDAS *f1ts=_time_steps[i];
      if(f1ts)
        {
          ret1[i]=f1ts->getTime(ret[i].first,ret[i].second);
        }
      else
        {
          std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutDAS::getTimeSteps : At rank #" << i << " time step is not defined !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return ret;
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
      _time_steps[i]=MEDFileField1TSWithoutDAS::New(_name.c_str(),i+1,_field_type,numdt,numo,_infos);
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
  _infos=arr->getInfoOnComponents();
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
int MEDFileFieldMultiTSWithoutDAS::getNonEmptyLevels(int iteration, int order, const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception)
{
  return getTimeStepEntry(iteration,order).getNonEmptyLevels(mname,levs);
}

std::vector< std::vector<TypeOfField> > MEDFileFieldMultiTSWithoutDAS::getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception)
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
std::vector< std::vector< std::pair<int,int> > > MEDFileFieldMultiTSWithoutDAS::getFieldSplitedByType(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  return getTimeStepEntry(iteration,order).getFieldSplitedByType(mname,types,typesF,pfls,locs);
}

/*!
 * entry point for users that want to iterate into MEDFile DataStructure with a reduced overhead because output arrays are extracted (created) specially
 * for the call of this method. That's why the DataArrayDouble instance in returned vector of vector should be dealed by the caller.
 */
std::vector< std::vector<DataArrayDouble *> > MEDFileFieldMultiTSWithoutDAS::getFieldSplitedByType2(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception)
{
  return getTimeStepEntry(iteration,order).getFieldSplitedByType2(mname,types,typesF,pfls,locs);
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

MEDFileFieldMultiTS *MEDFileFieldMultiTS::New(const MEDFileFieldMultiTSWithoutDAS& other)
{
  return new MEDFileFieldMultiTS(other);
}

std::string MEDFileFieldMultiTS::simpleRepr() const
{
  std::ostringstream oss;
  MEDFileFieldMultiTSWithoutDAS::simpleRepr(0,oss,-1);
  MEDFieldFieldGlobsReal::simpleRepr(oss);
  return oss.str();
}

void MEDFileFieldMultiTS::writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception)
{
  writeGlobals(fid,*this);
  MEDFileFieldMultiTSWithoutDAS::writeLL(fid);
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
  const MEDFileField1TSWithoutDAS& myF1TS=getTimeStepEntry(iteration,order);
  return myF1TS.getFieldAtLevel(type,meshDimRelToMax,0,renumPol,this);
}

MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldAtTopLevel(TypeOfField type, int iteration, int order, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileField1TSWithoutDAS& myF1TS=getTimeStepEntry(iteration,order);
  return myF1TS.getFieldAtTopLevel(type,0,renumPol,this);
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

/*!
 * This method has a close behaviour than MEDFileFieldMultiTS::getFieldAtLevel.
 * This method is called 'old' because the user should give the mesh name he wants to use for it's field.
 * This method is useful for MED2 file format when field on different mesh was autorized.
 */
MEDCouplingFieldDouble *MEDFileFieldMultiTS::getFieldAtLevelOld(TypeOfField type, const char *mname, int iteration, int order, int meshDimRelToMax, int renumPol) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileField1TSWithoutDAS& myF1TS=getTimeStepEntry(iteration,order);
  return myF1TS.getFieldAtLevel(type,meshDimRelToMax,mname,renumPol,this);
}

DataArrayDouble *MEDFileFieldMultiTS::getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const throw(INTERP_KERNEL::Exception)
{
  const MEDFileField1TSWithoutDAS& myF1TS=getTimeStepEntry(iteration,order);
  return myF1TS.getFieldWithProfile(type,meshDimRelToMax,mesh,pfl,this);
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
try:MEDFileFieldMultiTSWithoutDAS(fieldName),MEDFieldFieldGlobsReal(fileName)
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
          _field_type=MEDFileUtilities::TraduceFieldType(typcha);
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

MEDFileFieldMultiTS::MEDFileFieldMultiTS(const MEDFileFieldMultiTSWithoutDAS& other):MEDFileFieldMultiTSWithoutDAS(other)
{
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

std::vector<std::string> MEDFileFields::getFieldsNames() const throw(INTERP_KERNEL::Exception)
{
  std::vector<std::string> ret(_fields.size());
  int i=0;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutDAS> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      const MEDFileFieldMultiTSWithoutDAS *f=(*it);
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
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutDAS> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      const MEDFileFieldMultiTSWithoutDAS *cur=(*it);
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
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutDAS> >::const_iterator it=_fields.begin();it!=_fields.end();it++,i++)
    {
      const MEDFileFieldMultiTSWithoutDAS *cur=(*it);
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
  MEDFieldFieldGlobsReal::simpleRepr(oss);
}

MEDFileFields::MEDFileFields()
{
}

MEDFileFields::MEDFileFields(const char *fileName) throw(INTERP_KERNEL::Exception)
try:MEDFieldFieldGlobsReal(fileName)
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
        _fields[i]=MEDFileFieldMultiTSWithoutDAS::New(fid,nomcha,i+1,ft,infos,nbOfStep);
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
  appendGlobs(*field,1e-12);
}

void MEDFileFields::setFieldAtPos(int i, MEDFileFieldMultiTS *field) throw(INTERP_KERNEL::Exception)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileFields::setFieldAtPos : invalid input pointer ! should be different from 0 !");
  if(i>=(int)_fields.size())
    _fields.resize(i+1);
  field->incrRef();
  _fields[i]=field;
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

MEDFileFieldMultiTS *MEDFileFields::getFieldAtPos(int i) const throw(INTERP_KERNEL::Exception)
{
  if(i<0 || i>=(int)_fields.size())
    {
      std::ostringstream oss; oss << "MEDFileFields::getFieldAtPos : Invalid given id in input (" << i << ") should be in [0," << _fields.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const MEDFileFieldMultiTSWithoutDAS *fmts=_fields[i];
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTS> ret=MEDFileFieldMultiTS::New(*fmts);
  ret->shallowCpyGlobs(*this);
  ret->incrRef();
  return ret;
}

MEDFileFieldMultiTS *MEDFileFields::getFieldWithName(const char *fieldName) const throw(INTERP_KERNEL::Exception)
{
  return getFieldAtPos(getPosFromFieldName(fieldName));
}

int MEDFileFields::getPosFromFieldName(const char *fieldName) const throw(INTERP_KERNEL::Exception)
{
  std::string tmp(fieldName);
  std::vector<std::string> poss;
  for(std::size_t i=0;i<_fields.size();i++)
    {
      const MEDFileFieldMultiTSWithoutDAS *f=_fields[i];
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
