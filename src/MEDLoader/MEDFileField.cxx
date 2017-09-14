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

#include "MEDFileField.txx"
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

extern med_geometry_type typmai[MED_N_CELL_FIXED_GEO];
extern INTERP_KERNEL::NormalizedCellType typmai2[MED_N_CELL_FIXED_GEO];
extern med_geometry_type typmainoeud[1];
extern med_geometry_type typmai3[34];

using namespace MEDCoupling;

template class MEDCoupling::MEDFileField1TSTemplateWithoutSDA<int>;
template class MEDCoupling::MEDFileField1TSTemplateWithoutSDA<float>;
template class MEDCoupling::MEDFileField1TSTemplateWithoutSDA<double>;
template class MEDCoupling::MEDFileField1TSNDTemplateWithoutSDA<int>;
template class MEDCoupling::MEDFileField1TSNDTemplateWithoutSDA<float>;
template class MEDCoupling::MEDFileNDTemplateFieldMultiTSWithoutSDA<int>;
template class MEDCoupling::MEDFileNDTemplateFieldMultiTSWithoutSDA<float>;
template class MEDCoupling::MEDFileTemplateField1TS<int>;
template class MEDCoupling::MEDFileTemplateField1TS<float>;
template class MEDCoupling::MEDFileTemplateField1TS<double>;
template class MEDCoupling::MEDFileNDTemplateField1TS<int>;
template class MEDCoupling::MEDFileNDTemplateField1TS<float>;
template class MEDCoupling::MEDFileTemplateFieldMultiTS<int>;
template class MEDCoupling::MEDFileTemplateFieldMultiTS<float>;
template class MEDCoupling::MEDFileTemplateFieldMultiTS<double>;
template class MEDCoupling::MEDFileNDTemplateFieldMultiTS<int>;
template class MEDCoupling::MEDFileNDTemplateFieldMultiTS<float>;

const char MEDFileField1TSWithoutSDA::TYPE_STR[]="FLOAT64";
const char MEDFileIntField1TSWithoutSDA::TYPE_STR[]="INT32";
const char MEDFileFloatField1TSWithoutSDA::TYPE_STR[]="FLOAT32";

MEDFileGTKeeper::~MEDFileGTKeeper()
{
}

MEDFileGTKeeper *MEDFileGTKeeperSta::deepCopy() const
{
  return new MEDFileGTKeeperSta(_geo_type);
}

INTERP_KERNEL::NormalizedCellType MEDFileGTKeeperSta::getGeoType() const
{
  return _geo_type;
}

std::string MEDFileGTKeeperSta::getRepr() const
{
  return INTERP_KERNEL::CellModel::GetCellModel(_geo_type).getRepr();
}

bool MEDFileGTKeeperSta::isEqual(const MEDFileGTKeeper *other) const
{
  const MEDFileGTKeeperSta *otherC(dynamic_cast<const MEDFileGTKeeperSta *>(other));
  if(!otherC)
    return false;
  return _geo_type==otherC->_geo_type;
}

MEDFileGTKeeperDyn::MEDFileGTKeeperDyn(const MEDFileUMesh *mesh, const MEDFileUMesh *section, const MEDFileStructureElement *se):_mesh(mesh),_section(section),_se(se)
{
  if(mesh)
    mesh->incrRef();
  if(section)
    section->incrRef();
  if(se)
    se->incrRef();
  if(_mesh.isNull() || _section.isNull() || _se.isNull())
    throw INTERP_KERNEL::Exception("MEDFileGTKeeperDyn constructor : null pointer not allowed !");
}

MEDFileGTKeeper *MEDFileGTKeeperDyn::deepCopy() const
{
  return new MEDFileGTKeeperDyn(_mesh,_section,_se);
}

INTERP_KERNEL::NormalizedCellType MEDFileGTKeeperDyn::getGeoType() const
{
  throw INTERP_KERNEL::Exception("MEDFileGTKeeperDyn::getGeoType : not valid !");
}

std::string MEDFileGTKeeperDyn::getRepr() const
{
  std::ostringstream oss;
  oss << _se->getDynGT();
  return oss.str();
}

bool MEDFileGTKeeperDyn::isEqual(const MEDFileGTKeeper *other) const
{
  const MEDFileGTKeeperDyn *otherC(dynamic_cast<const MEDFileGTKeeperDyn *>(other));
  if(!otherC)
    return false;
  return this==otherC;
}

MEDFileFieldLoc *MEDFileFieldLoc::New(med_idt fid, const std::string& locName)
{
  return new MEDFileFieldLoc(fid,locName);
}

MEDFileFieldLoc *MEDFileFieldLoc::New(med_idt fid, int id, const MEDFileEntities *entities)
{
  return new MEDFileFieldLoc(fid,id,entities);
}

MEDFileFieldLoc *MEDFileFieldLoc::New(const std::string& locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w)
{
  return new MEDFileFieldLoc(locName,geoType,refCoo,gsCoo,w);
}

MEDFileFieldLoc::MEDFileFieldLoc(med_idt fid, const std::string& locName):_name(locName)
{
  med_geometry_type geotype;
  med_geometry_type sectiongeotype;
  int nsectionmeshcell;
  INTERP_KERNEL::AutoPtr<char> geointerpname(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> sectionmeshname(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  MEDlocalizationInfoByName(fid,locName.c_str(),&geotype,&_dim,&_nb_gauss_pt,geointerpname,sectionmeshname,&nsectionmeshcell,&sectiongeotype);
  _gt=new MEDFileGTKeeperSta((INTERP_KERNEL::NormalizedCellType)(std::distance(typmai3,std::find(typmai3,typmai3+34,geotype))));
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(getGeoType()));
  _nb_node_per_cell=cm.getNumberOfNodes();
  _ref_coo.resize(_dim*_nb_node_per_cell);
  _gs_coo.resize(_dim*_nb_gauss_pt);
  _w.resize(_nb_gauss_pt);
  MEDlocalizationRd(fid,locName.c_str(),MED_FULL_INTERLACE,&_ref_coo[0],&_gs_coo[0],&_w[0]);
}

MEDFileFieldLoc::MEDFileFieldLoc(med_idt fid, int id, const MEDFileEntities *entities)
{
  med_geometry_type geotype;
  med_geometry_type sectiongeotype;
  int nsectionmeshcell;
  INTERP_KERNEL::AutoPtr<char> locName(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> geointerpname(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> sectionmeshname(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  MEDlocalizationInfo(fid,id+1,locName,&geotype,&_dim,&_nb_gauss_pt,geointerpname,sectionmeshname,&nsectionmeshcell,&sectiongeotype);
  _name=locName;
  std::string sectionName(MEDLoaderBase::buildStringFromFortran(sectionmeshname,MED_NAME_SIZE));
  if(sectionName.empty())
    {
      _gt=new MEDFileGTKeeperSta((INTERP_KERNEL::NormalizedCellType)(std::distance(typmai3,std::find(typmai3,typmai3+34,geotype))));
      const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(getGeoType()));
      _nb_node_per_cell=cm.getNumberOfNodes();
    }
  else
    {
      const MEDFileAllStaticEntitiesPlusDyn *entities2(dynamic_cast<const MEDFileAllStaticEntitiesPlusDyn *>(entities));
      if(!entities2)
        {
          std::ostringstream oss; oss << "MEDFileFieldLoc cstr : for loc \"" << _name << "\" presence of non static type ! Expect entities !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      const MEDFileStructureElement *se(entities2->getWithGT(geotype));
      const MEDFileUMesh *um(entities2->getSupMeshWithName(se->getMeshName()));
      const MEDFileUMesh *section(entities2->getSupMeshWithName(sectionName));
      _gt=new MEDFileGTKeeperDyn(um,section,se);
      {
        int dummy;
        MEDFILESAFECALLERWR0(MEDmeshGeotypeParameter,(fid,geotype,&dummy,&_nb_node_per_cell));
      }
    }
  _ref_coo.resize(_dim*_nb_node_per_cell);
  _gs_coo.resize(_dim*_nb_gauss_pt);
  _w.resize(_nb_gauss_pt);
  MEDlocalizationRd(fid,locName,MED_FULL_INTERLACE,&_ref_coo[0],&_gs_coo[0],&_w[0]);
}

MEDFileFieldLoc::MEDFileFieldLoc(const std::string& locName, INTERP_KERNEL::NormalizedCellType geoType,
                                 const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w):_name(locName),_gt(new MEDFileGTKeeperSta(geoType)),_ref_coo(refCoo),_gs_coo(gsCoo),_w(w)
{
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(getGeoType()));
  _dim=cm.getDimension();
  _nb_node_per_cell=cm.getNumberOfNodes();
  _nb_gauss_pt=_w.size();
}


MEDFileFieldLoc::MEDFileFieldLoc(const MEDFileFieldLoc& other):_dim(other._dim),_nb_gauss_pt(other._nb_gauss_pt),_gt(other._gt->deepCopy()),_nb_node_per_cell(other._nb_node_per_cell),_name(other._name),_ref_coo(other._ref_coo),_gs_coo(other._gs_coo),_w(other._w)
{
}

MEDFileFieldLoc *MEDFileFieldLoc::deepCopy() const
{
  return new MEDFileFieldLoc(*this);
}

bool MEDFileFieldLoc::isOnStructureElement() const
{
  const MEDFileGTKeeper *gt(_gt);
  if(!gt)
    throw INTERP_KERNEL::Exception("MEDFileFieldLoc::isOnStructureElement : null pointer !");
  const MEDFileGTKeeperDyn *gt2(dynamic_cast<const MEDFileGTKeeperDyn *>(gt));
  return gt2!=NULL;
}

std::size_t MEDFileFieldLoc::getHeapMemorySizeWithoutChildren() const
{
  return (_ref_coo.capacity()+_gs_coo.capacity()+_w.capacity())*sizeof(double)+_name.capacity();
}

std::vector<const BigMemoryObject *> MEDFileFieldLoc::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

void MEDFileFieldLoc::simpleRepr(std::ostream& oss) const
{
  static const char OFF7[]="\n    ";
  oss << "\"" << _name << "\"" << OFF7;
  oss << "GeoType=" << _gt->getRepr() << OFF7;
  oss << "Dimension=" << _dim << OFF7;
  oss << "Number of Gauss points=" << _nb_gauss_pt << OFF7;
  oss << "Number of nodes per cell=" << _nb_node_per_cell << OFF7;
  oss << "RefCoords="; std::copy(_ref_coo.begin(),_ref_coo.end(),std::ostream_iterator<double>(oss," ")); oss << OFF7;
  oss << "Weights="; std::copy(_w.begin(),_w.end(),std::ostream_iterator<double>(oss," ")); oss << OFF7;
  oss << "GaussPtsCoords="; std::copy(_gs_coo.begin(),_gs_coo.end(),std::ostream_iterator<double>(oss," ")); oss << std::endl;
}

void MEDFileFieldLoc::setName(const std::string& name)
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
  if(!_gt->isEqual(other._gt))
    return false;
  if(!MEDCouplingGaussLocalization::AreAlmostEqual(_ref_coo,other._ref_coo,eps))
    return false;
  if(!MEDCouplingGaussLocalization::AreAlmostEqual(_gs_coo,other._gs_coo,eps))
    return false;
  if(!MEDCouplingGaussLocalization::AreAlmostEqual(_w,other._w,eps))
    return false;

  return true;
}

void MEDFileFieldLoc::writeLL(med_idt fid) const
{
  MEDlocalizationWr(fid,_name.c_str(),typmai3[(int)getGeoType()],_dim,&_ref_coo[0],MED_FULL_INTERLACE,_nb_gauss_pt,&_gs_coo[0],&_w[0],MED_NO_INTERPOLATION,MED_NO_MESH_SUPPORT);
}

std::string MEDFileFieldLoc::repr() const
{
  std::ostringstream oss; oss.precision(15);
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(getGeoType()));
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

void MEDFileFieldPerMeshPerTypePerDisc::assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldTemplate *field, const DataArray *arrr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc)
{
  _type=field->getTypeOfField();
  _start=start;
  switch(_type)
  {
    case ON_CELLS:
      {
        getOrCreateAndGetArray()->setContigPartOfSelectedValuesSlice(_start,arrr,offset,offset+nbOfCells,1);
        _end=_start+nbOfCells;
        _nval=nbOfCells;
        break;
      }
    case ON_GAUSS_NE:
      {
        MCAuto<DataArrayInt> arr=field->getDiscretization()->getOffsetArr(field->getMesh());
        const int *arrPtr=arr->getConstPointer();
        getOrCreateAndGetArray()->setContigPartOfSelectedValuesSlice(_start,arrr,arrPtr[offset],arrPtr[offset+nbOfCells],1);
        _end=_start+(arrPtr[offset+nbOfCells]-arrPtr[offset]);
        _nval=nbOfCells;
        break;
      }
    case ON_GAUSS_PT:
      {
        const MEDCouplingFieldDiscretization *disc(field->getDiscretization());
        const MEDCouplingGaussLocalization& gsLoc(field->getGaussLocalization(_loc_id));
        const MEDCouplingFieldDiscretizationGauss *disc2(dynamic_cast<const MEDCouplingFieldDiscretizationGauss *>(disc));
        if(!disc2)
          throw INTERP_KERNEL::Exception("assignFieldNoProfile : invalid call to this method ! Internal Error !");
        const DataArrayInt *dai(disc2->getArrayOfDiscIds());
        MCAuto<DataArrayInt> dai2(disc2->getOffsetArr(field->getMesh()));
        const int *dai2Ptr(dai2->getConstPointer());
        int nbi(gsLoc.getWeights().size());
        MCAuto<DataArrayInt> da2(dai->selectByTupleIdSafeSlice(offset,offset+nbOfCells,1));
        MCAuto<DataArrayInt> da3(da2->findIdsEqual(_loc_id));
        const int *da3Ptr(da3->getConstPointer());
        if(da3->getNumberOfTuples()!=nbOfCells)
          {//profile : for gauss even in NoProfile !!!
            std::ostringstream oss; oss << "Pfl_" << nasc.getName() << "_" << INTERP_KERNEL::CellModel::GetCellModel(getGeoType()).getRepr() << "_" << _loc_id;
            _profile=oss.str();
            da3->setName(_profile.c_str());
            glob.appendProfile(da3);
          }
        MCAuto<DataArrayInt> da4(DataArrayInt::New());
        _nval=da3->getNbOfElems();
        da4->alloc(_nval*nbi,1);
        int *da4Ptr(da4->getPointer());
        for(int i=0;i<_nval;i++)
          {
            int ref=dai2Ptr[offset+da3Ptr[i]];
            for(int j=0;j<nbi;j++)
              *da4Ptr++=ref+j;
          }
        std::ostringstream oss2; oss2 << "Loc_" << nasc.getName() << "_" << INTERP_KERNEL::CellModel::GetCellModel(getGeoType()).getRepr() << "_" << _loc_id;
        _localization=oss2.str();
        getOrCreateAndGetArray()->setContigPartOfSelectedValues(_start,arrr,da4);
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
void MEDFileFieldPerMeshPerTypePerDisc::assignFieldProfile(bool isPflAlone, int& start, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, DataArrayInt *locIds, int nbOfEltsInWholeMesh, const MEDCouplingFieldTemplate *field, const DataArray *arrr, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc)
{
  _profile.clear();
  _type=field->getTypeOfField();
  std::string pflName(multiTypePfl->getName());
  std::ostringstream oss; oss << pflName;
  if(_type!=ON_NODES)
    {
      if(!isPflAlone)
        { const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(getGeoType()); oss << "_" <<  cm.getRepr(); }
    }
  else
    { oss << "_NODE"; }
  if(locIds)
    {
      if(pflName.empty())
        throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::assignFieldProfile : existing profile with empty name !");
      if(_type!=ON_GAUSS_PT)
        {
          locIds->setName(oss.str());
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
        getOrCreateAndGetArray()->setContigPartOfSelectedValuesSlice(_start,arrr,0,arrr->getNumberOfTuples(),1);
        _end=_start+_nval;
        break;
      }
    case ON_CELLS:
      {
        _nval=idsInPfl->getNumberOfTuples();
        getOrCreateAndGetArray()->setContigPartOfSelectedValues(_start,arrr,idsInPfl);
        _end=_start+_nval;
        break;
      }
    case ON_GAUSS_NE:
      {
        MCAuto<DataArrayInt> arr=field->getDiscretization()->getOffsetArr(mesh);
        MCAuto<DataArrayInt> arr2=arr->deltaShiftIndex();
        MCAuto<DataArrayInt> arr3=arr2->selectByTupleId(multiTypePfl->begin(),multiTypePfl->end());
        arr3->computeOffsetsFull();
        MCAuto<DataArrayInt> tmp=idsInPfl->buildExplicitArrByRanges(arr3);
        int trueNval=tmp->getNumberOfTuples();
        _nval=idsInPfl->getNumberOfTuples();
        getOrCreateAndGetArray()->setContigPartOfSelectedValues(_start,arrr,tmp);
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
        MCAuto<DataArrayInt> da2=da1->selectByTupleId(idsInPfl->begin(),idsInPfl->end());
        MCAuto<DataArrayInt> da3=da2->findIdsEqual(_loc_id);
        MCAuto<DataArrayInt> da4=idsInPfl->selectByTupleId(da3->begin(),da3->end());
        //
        MCAuto<MEDCouplingMesh> mesh2=mesh->buildPart(multiTypePfl->begin(),multiTypePfl->end());
        MCAuto<DataArrayInt> arr=disc2->getOffsetArr(mesh2);
        //
        MCAuto<DataArrayInt> tmp=DataArrayInt::New();
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
        getOrCreateAndGetArray()->setContigPartOfSelectedValues(_start,arrr,tmp);
        _end=_start+trueNval;
        oss << "_loc_" << _loc_id;
        if(locIds)
          {
            MCAuto<DataArrayInt> da5=locIds->selectByTupleId(da3->begin(),da3->end());
            da5->setName(oss.str());
            glob.appendProfile(da5);
            _profile=oss.str();
          }
        else
          {
            if(!da3->isIota(nbOfEltsInWholeMesh))
              {
                da3->setName(oss.str());
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

void MEDFileFieldPerMeshPerTypePerDisc::assignNodeFieldNoProfile(int& start, const MEDCouplingFieldTemplate *field, const DataArray *arrr, MEDFileFieldGlobsReal& glob)
{
  _start=start;
  _nval=arrr->getNumberOfTuples();
  getOrCreateAndGetArray()->setContigPartOfSelectedValuesSlice(_start,arrr,0,_nval,1);
  _end=_start+_nval;
  start=_end;
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypePerDisc::NewOnRead(MEDFileFieldPerMeshPerTypeCommon *fath, TypeOfField type, int profileIt, const PartDefinition *pd)
{
  return new MEDFileFieldPerMeshPerTypePerDisc(fath,type,profileIt,pd);
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypePerDisc::New(MEDFileFieldPerMeshPerTypeCommon *fath, TypeOfField type, int locId)
{
  return new MEDFileFieldPerMeshPerTypePerDisc(fath,type,locId,std::string());
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypePerDisc::New(const MEDFileFieldPerMeshPerTypePerDisc& other)
{
  return new MEDFileFieldPerMeshPerTypePerDisc(other);
}

std::size_t MEDFileFieldPerMeshPerTypePerDisc::getHeapMemorySizeWithoutChildren() const
{
  return _profile.capacity()+_localization.capacity()+sizeof(MEDFileFieldPerMeshPerTypePerDisc);
}

std::vector<const BigMemoryObject *> MEDFileFieldPerMeshPerTypePerDisc::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(1);
  ret[0]=(const PartDefinition*)_pd;
  return ret;
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypePerDisc::deepCopy(MEDFileFieldPerMeshPerTypeCommon *father) const
{
  MCAuto<MEDFileFieldPerMeshPerTypePerDisc> ret(new MEDFileFieldPerMeshPerTypePerDisc(*this));
  ret->_father=father;
  return ret.retn();
}

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerTypeCommon *fath, TypeOfField atype, int profileIt, const PartDefinition *pd)
try:_type(atype),_father(fath),_profile_it(profileIt),_pd(const_cast<PartDefinition *>(pd))
{
  if(pd)
    pd->incrRef();
}
catch(INTERP_KERNEL::Exception& e)
{
    throw e;
}

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerTypeCommon *fath, TypeOfField type, int locId, const std::string& dummy):_type(type),_father(fath),_loc_id(locId)
{
}

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc(const MEDFileFieldPerMeshPerTypePerDisc& other):RefCountObject(other),_type(other._type),_father(0),_start(other._start),_end(other._end),_nval(other._nval),_profile(other._profile),_localization(other._localization),_loc_id(other._loc_id),_profile_it(other._profile_it),_pd(other._pd),_tmp_work1(other._tmp_work1)
{
}

MEDFileFieldPerMeshPerTypePerDisc::MEDFileFieldPerMeshPerTypePerDisc():_type(ON_CELLS),_father(0),_start(-std::numeric_limits<int>::max()),_end(-std::numeric_limits<int>::max()),
    _nval(-std::numeric_limits<int>::max()),_loc_id(-std::numeric_limits<int>::max())
{
}

void MEDFileFieldPerMeshPerTypePerDisc::goReadZeValuesInFile(med_idt fid, const std::string& fieldName, int nbOfCompo, int iteration, int order, med_entity_type menti, med_geometry_type mgeoti, unsigned char *startFeedingPtr)
{
  const PartDefinition *pd(_pd);
  if(!pd)
    {
      med_entity_type mentiCpy(menti);
      INTERP_KERNEL::AutoPtr<char> locname(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
      int nbi,tmp1;
      med_int nbValsInFile(MEDfieldnValueWithProfileByName(fid,fieldName.c_str(),iteration,order,menti,mgeoti,_profile.c_str(),MED_COMPACT_PFLMODE,&tmp1,locname,&nbi));
      if(nbValsInFile==0 && menti==MED_CELL)
        {//
          nbValsInFile=MEDfieldnValueWithProfileByName(fid,fieldName.c_str(),iteration,order,MED_DESCENDING_FACE,mgeoti,_profile.c_str(),MED_COMPACT_PFLMODE,&tmp1,locname,&nbi);
          if(nbValsInFile==0)
            {
              nbValsInFile=MEDfieldnValueWithProfileByName(fid,fieldName.c_str(),iteration,order,MED_DESCENDING_EDGE,mgeoti,_profile.c_str(),MED_COMPACT_PFLMODE,&tmp1,locname,&nbi);
              if(nbValsInFile!=0)
                { mentiCpy=MED_DESCENDING_EDGE; }
            }
          else
            { mentiCpy=MED_DESCENDING_FACE; }
        }
      if(_end-_start!=nbValsInFile*nbi)
        {
          std::ostringstream oss; oss << "MEDFileFieldPerMeshPerTypePerDisc::goReadZeValuesInFile : The number of tuples to read is " << nbValsInFile << "*" << nbi <<  " (nb integration points) ! But in data structure it values " << _end-_start << " is expected !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      MEDFILESAFECALLERRD0(MEDfieldValueWithProfileRd,(fid,fieldName.c_str(),iteration,order,mentiCpy,mgeoti,MED_COMPACT_PFLMODE,_profile.c_str(),MED_FULL_INTERLACE,MED_ALL_CONSTITUENT,startFeedingPtr));
    }
  else
    {
      if(!_profile.empty())
        throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::goReadZeValuesInFile : not implemented !");
      INTERP_KERNEL::AutoPtr<char> pflname(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE)),locname(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
      int profilesize,nbi;
      int overallNval(MEDfieldnValueWithProfile(fid,fieldName.c_str(),iteration,order,menti,mgeoti,_profile_it+1,MED_COMPACT_PFLMODE,pflname,&profilesize,locname,&nbi));
      const SlicePartDefinition *spd(dynamic_cast<const SlicePartDefinition *>(pd));
      if(spd)
        {
          int start,stop,step;
          spd->getSlice(start,stop,step);
          int nbOfEltsToLoad(DataArray::GetNumberOfItemGivenBES(start,stop,step,"MEDFileFieldPerMeshPerTypePerDisc::goReadZeValuesInFile"));
          med_filter filter=MED_FILTER_INIT;
          MEDfilterBlockOfEntityCr(fid,/*nentity*/overallNval,/*nvaluesperentity*/nbi,/*nconstituentpervalue*/nbOfCompo,
                                   MED_ALL_CONSTITUENT,MED_FULL_INTERLACE,MED_COMPACT_STMODE,MED_NO_PROFILE,
                                   /*start*/start+1,/*stride*/step,/*count*/1,/*blocksize*/nbOfEltsToLoad,
                                   /*lastblocksize=useless because count=1*/0,&filter);
          MEDFILESAFECALLERRD0(MEDfieldValueAdvancedRd,(fid,fieldName.c_str(),iteration,order,menti,mgeoti,&filter,startFeedingPtr));
          MEDfilterClose(&filter);
          return ;
        }
      const DataArrayPartDefinition *dpd(dynamic_cast<const DataArrayPartDefinition *>(pd));
      if(dpd)
        {
          dpd->checkConsistencyLight();
          MCAuto<DataArrayInt> myIds(dpd->toDAI());
          int a(myIds->getMinValueInArray()),b(myIds->getMaxValueInArray());
          myIds=myIds->deepCopy();// WARNING deep copy here because _pd is modified by applyLin !!!
          myIds->applyLin(1,-a);
          int nbOfEltsToLoad(b-a+1);
          med_filter filter=MED_FILTER_INIT;
          {//TODO : manage int32 !
            MCAuto<DataArrayDouble> tmp(DataArrayDouble::New());
            tmp->alloc(nbOfEltsToLoad,nbOfCompo);
            MEDfilterBlockOfEntityCr(fid,/*nentity*/overallNval,/*nvaluesperentity*/nbi,/*nconstituentpervalue*/nbOfCompo,
                                     MED_ALL_CONSTITUENT,MED_FULL_INTERLACE,MED_COMPACT_STMODE,MED_NO_PROFILE,
                                     /*start*/a+1,/*stride*/1,/*count*/1,/*blocksize*/nbOfEltsToLoad,
                                     /*lastblocksize=useless because count=1*/0,&filter);
            MEDFILESAFECALLERRD0(MEDfieldValueAdvancedRd,(fid,fieldName.c_str(),iteration,order,menti,mgeoti,&filter,reinterpret_cast<unsigned char *>(tmp->getPointer())));
            MCAuto<DataArrayDouble> feeder(DataArrayDouble::New());
            feeder->useExternalArrayWithRWAccess(reinterpret_cast<double *>(startFeedingPtr),_nval,nbOfCompo);
            feeder->setContigPartOfSelectedValues(0,tmp,myIds);
          }
          MEDfilterClose(&filter);
        }
      else
        throw INTERP_KERNEL::Exception("Not implemented yet for not slices!");
    }
}

const MEDFileFieldPerMeshPerTypeCommon *MEDFileFieldPerMeshPerTypePerDisc::getFather() const
{
  return _father;
}

void MEDFileFieldPerMeshPerTypePerDisc::loadOnlyStructureOfDataRecursively(med_idt fid, int& start, const MEDFileFieldNameScope& nasc)
{
  INTERP_KERNEL::AutoPtr<char> locname(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> pflname(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  std::string fieldName(nasc.getName()),meshName(getMeshName());
  int iteration(getIteration()),order(getOrder()),profilesize,nbi;
  TypeOfField type(getType());
  med_geometry_type mgeoti;
  med_entity_type menti;
  _father->entriesForMEDfile(type,mgeoti,menti);
  int zeNVal(MEDfieldnValueWithProfile(fid,fieldName.c_str(),iteration,order,menti,mgeoti,_profile_it+1,MED_COMPACT_PFLMODE,pflname,&profilesize,locname,&nbi));
  if(zeNVal==0 && type==ON_CELLS)
    {//eheh maybe there's a surprise :)
      int zeNVal1(MEDfieldnValueWithProfile(fid,fieldName.c_str(),iteration,order,MED_DESCENDING_FACE,mgeoti,_profile_it+1,MED_COMPACT_PFLMODE,pflname,&profilesize,locname,&nbi));
      if(zeNVal1==0)
        {
          int zeNVal2(MEDfieldnValueWithProfile(fid,fieldName.c_str(),iteration,order,MED_DESCENDING_EDGE,mgeoti,_profile_it+1,MED_COMPACT_PFLMODE,pflname,&profilesize,locname,&nbi));
          if(zeNVal2!=0)
            zeNVal=zeNVal2;
        }
      else
        {
          zeNVal=zeNVal1;
        }
    }
  _profile=MEDLoaderBase::buildStringFromFortran(pflname,MED_NAME_SIZE);
  _localization=MEDLoaderBase::buildStringFromFortran(locname,MED_NAME_SIZE);
  const PartDefinition *pd(_pd);
  if(!pd)
    {
      _nval=zeNVal;
    }
  else
    {
      if(!_profile.empty())
        throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::loadOnlyStructureOfDataRecursively : profiles are not managed yet with part of def !");
      _nval=pd->getNumberOfElems();
    }
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

void MEDFileFieldPerMeshPerTypePerDisc::loadBigArray(med_idt fid, const MEDFileFieldNameScope& nasc)
{
  std::string fieldName(nasc.getName()),meshName(getMeshName());
  int iteration(getIteration()),order(getOrder());
  TypeOfField type(getType());
  med_geometry_type mgeoti;
  med_entity_type menti;
  _father->entriesForMEDfile(type,mgeoti,menti);
  if(_start>_end)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::loadBigArray : internal error in range !");
  if(_start==_end)
    return ;
  DataArray *arr(getOrCreateAndGetArray());//arr is not null due to the spec of getOrCreateAndGetArray
  if(_start<0 || _start>=arr->getNumberOfTuples())
    {
      std::ostringstream oss; oss << "MEDFileFieldPerMeshPerTypePerDisc::loadBigArray : Invalid start ("<< _start << ") regarding admissible range of allocated array [0," << arr->getNumberOfTuples() << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  if(_end<0 || _end>arr->getNumberOfTuples())
    {
      std::ostringstream oss; oss << "MEDFileFieldPerMeshPerTypePerDisc::loadBigArray : Invalid start ("<< _start << ") regarding admissible range of allocated array [0," << arr->getNumberOfTuples() << "] !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  int nbOfCompo(arr->getNumberOfComponents());
  DataArrayDouble *arrD(dynamic_cast<DataArrayDouble *>(arr));
  if(arrD)
    {
      double *startFeeding(arrD->getPointer()+_start*nbOfCompo);
      goReadZeValuesInFile(fid,fieldName,nbOfCompo,iteration,order,menti,mgeoti,reinterpret_cast<unsigned char*>(startFeeding));
      return ;
    }
  DataArrayInt *arrI(dynamic_cast<DataArrayInt *>(arr));
  if(arrI)
    {
      int *startFeeding(arrI->getPointer()+_start*nbOfCompo);
      goReadZeValuesInFile(fid,fieldName,nbOfCompo,iteration,order,menti,mgeoti,reinterpret_cast<unsigned char*>(startFeeding));
      return ;
    }
  DataArrayFloat *arrF(dynamic_cast<DataArrayFloat *>(arr));
  if(arrF)
    {
      float *startFeeding(arrF->getPointer()+_start*nbOfCompo);
      goReadZeValuesInFile(fid,fieldName,nbOfCompo,iteration,order,menti,mgeoti,reinterpret_cast<unsigned char*>(startFeeding));
      return ;
    }
  throw INTERP_KERNEL::Exception("Error on array reading ! Unrecognized type of field ! Should be in FLOAT64 FLOAT32 or INT32 !");
}

/*!
 * Set a \c this->_start **and** \c this->_end keeping the same delta between the two.
 */
void MEDFileFieldPerMeshPerTypePerDisc::setNewStart(int newValueOfStart)
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
  INTERP_KERNEL::AutoCppPtr<MEDCouplingFieldDiscretization> tmp(MEDCouplingFieldDiscretization::New(_type));
  oss << startLine2 << "Localization #" << id << "." << std::endl;
  oss << startLine2 << "  Type=" << tmp->getRepr() << "." << std::endl;
  oss << startLine2 << "  This type discretization lies on profile : \"" << _profile << "\" and on the following localization : \"" << _localization << "\"." << std::endl;
  oss << startLine2 << "  This type discretization has " << _end-_start << " tuples (start=" << _start << ", end=" << _end << ")." << std::endl;
  oss << startLine2 << "  This type discretization has " << (_end-_start)/_nval << " integration points." << std::endl;
}

TypeOfField MEDFileFieldPerMeshPerTypePerDisc::getType() const
{
  return _type;
}

INTERP_KERNEL::NormalizedCellType MEDFileFieldPerMeshPerTypePerDisc::getGeoType() const
{
  return _father->getGeoType();
}

void MEDFileFieldPerMeshPerTypePerDisc::fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const
{
  types.insert(_type);
}

void MEDFileFieldPerMeshPerTypePerDisc::setType(TypeOfField newType)
{
  _type=newType;
}

int MEDFileFieldPerMeshPerTypePerDisc::getNumberOfComponents() const
{
  return _father->getNumberOfComponents();
}

int MEDFileFieldPerMeshPerTypePerDisc::getNumberOfTuples() const
{
  return _end-_start;
}

DataArray *MEDFileFieldPerMeshPerTypePerDisc::getOrCreateAndGetArray()
{
  return _father->getOrCreateAndGetArray();
}

const DataArray *MEDFileFieldPerMeshPerTypePerDisc::getOrCreateAndGetArray() const
{
  const MEDFileFieldPerMeshPerTypeCommon *fath=_father;
  return fath->getOrCreateAndGetArray();
}

const std::vector<std::string>& MEDFileFieldPerMeshPerTypePerDisc::getInfo() const
{
  return _father->getInfo();
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getProfile() const
{
  return _profile;
}

void MEDFileFieldPerMeshPerTypePerDisc::setProfile(const std::string& newPflName)
{
  _profile=newPflName;
}

std::string MEDFileFieldPerMeshPerTypePerDisc::getLocalization() const
{
  return _localization;
}

void MEDFileFieldPerMeshPerTypePerDisc::setLocalization(const std::string& newLocName)
{
  _localization=newLocName;
}

void MEDFileFieldPerMeshPerTypePerDisc::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
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

void MEDFileFieldPerMeshPerTypePerDisc::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
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

void MEDFileFieldPerMeshPerTypePerDisc::writeLL(med_idt fid, const MEDFileFieldNameScope& nasc) const
{
  TypeOfField type=getType();
  INTERP_KERNEL::NormalizedCellType geoType(getGeoType());
  med_geometry_type mgeoti;
  med_entity_type menti;
  _father->entriesForMEDfile(getType(),mgeoti,menti);
  const DataArray *arr(getOrCreateAndGetArray());
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::writeLL : no array set !");
  if(!arr->isAllocated())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::writeLL : the array to be written is not allocated !");
  const DataArrayDouble *arrD(dynamic_cast<const DataArrayDouble *>(arr));
  const DataArrayInt *arrI(dynamic_cast<const DataArrayInt *>(arr));
  const DataArrayFloat *arrF(dynamic_cast<const DataArrayFloat *>(arr));
  const unsigned char *locToWrite=0;
  if(arrD)
    locToWrite=reinterpret_cast<const unsigned char *>(arrD->getConstPointer()+_start*arr->getNumberOfComponents());
  else if(arrI)
    locToWrite=reinterpret_cast<const unsigned char *>(arrI->getConstPointer()+_start*arr->getNumberOfComponents());
  else if(arrF)
    locToWrite=reinterpret_cast<const unsigned char *>(arrF->getConstPointer()+_start*arr->getNumberOfComponents());
  else
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::writeLL : not recognized type of values ! Supported are FLOAT64 FLOAT32 and INT32 !");
  MEDFILESAFECALLERWR0(MEDfieldValueWithProfileWr,(fid,nasc.getName().c_str(),getIteration(),getOrder(),getTime(),menti,mgeoti,
                                                   MED_COMPACT_PFLMODE,_profile.c_str(),_localization.c_str(),MED_FULL_INTERLACE,MED_ALL_CONSTITUENT,_nval,
                                                   locToWrite));
}

void MEDFileFieldPerMeshPerTypePerDisc::getCoarseData(TypeOfField& type, std::pair<int,int>& dad, std::string& pfl, std::string& loc) const
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
int MEDFileFieldPerMeshPerTypePerDisc::fillEltIdsFromCode(int offset, const std::vector<int>& codeOfMesh, const MEDFileFieldGlobsReal& glob, int *ptToFill) const
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
      throw INTERP_KERNEL::Exception(oss.str());
    }
  int *work=ptToFill;
  if(_profile.empty())
    {
      if(_nval!=codeOfMesh[3*found+1])
        {
          const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(getGeoType());
          oss << "MEDFileFieldPerMeshPerTypePerDisc::fillEltIdsFromCode : for geometric type " << cm.getRepr() << " number of elt ids in mesh is equal to " << _nval;
          oss << " whereas mesh has " << codeOfMesh[3*found+1] << " for this geometric type !";
          throw INTERP_KERNEL::Exception(oss.str());
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
          throw INTERP_KERNEL::Exception(oss.str());
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

int MEDFileFieldPerMeshPerTypePerDisc::fillTupleIds(int *ptToFill) const
{
  for(int i=_start;i<_end;i++)
    *ptToFill++=i;
  return _end-_start;
}

int MEDFileFieldPerMeshPerTypePerDisc::ConvertType(TypeOfField type, int locId)
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
                                                       std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >& result)
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
  MCAuto<DataArrayInt> renumTuples=DataArrayInt::New(); renumTuples->alloc(szTuples,1);
  MCAuto<DataArrayInt> ranges=MEDCouplingUMesh::ComputeRangesFromTypeDistribution(newCode);
  std::vector< MCAuto<DataArrayInt> > newGeoTypesPerChunk(entriesOnSameDisc.size());
  std::vector< const DataArrayInt * > newGeoTypesPerChunk2(entriesOnSameDisc.size());
  std::vector< MCAuto<DataArrayInt> > newGeoTypesPerChunk_bis(entriesOnSameDisc.size());
  std::vector< const DataArrayInt * > newGeoTypesPerChunk3(entriesOnSameDisc.size());
  MCAuto<DataArrayInt> newGeoTypesPerChunk4=DataArrayInt::New(); newGeoTypesPerChunk4->alloc(szEntities,nbi);
  int id=0;
  for(std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>::const_iterator it=entriesOnSameDisc.begin();it!=entriesOnSameDisc.end();it++,id++)
    {
      int startOfEltIdOfChunk=(*it)->_start;
      MCAuto<DataArrayInt> newEltIds=explicitIdsInMesh->subArray(startOfEltIdOfChunk,startOfEltIdOfChunk+(*it)->_nval);
      MCAuto<DataArrayInt> rangeIdsForChunk=newEltIds->findRangeIdForEachTuple(ranges);
      MCAuto<DataArrayInt> idsInRrangeForChunk=newEltIds->findIdInRangeForEachTuple(ranges);
      //
      MCAuto<DataArrayInt> tmp=rangeIdsForChunk->duplicateEachTupleNTimes(nbi); rangeIdsForChunk->rearrange(nbi);
      newGeoTypesPerChunk4->setPartOfValues1(tmp,(*it)->_tmp_work1-offset,(*it)->_tmp_work1+(*it)->_nval*nbi-offset,1,0,nbi,1);
      //
      newGeoTypesPerChunk[id]=rangeIdsForChunk; newGeoTypesPerChunk2[id]=rangeIdsForChunk;
      newGeoTypesPerChunk_bis[id]=idsInRrangeForChunk; newGeoTypesPerChunk3[id]=idsInRrangeForChunk;
    }
  MCAuto<DataArrayInt> newGeoTypesEltIdsAllGather=DataArrayInt::Aggregate(newGeoTypesPerChunk2); newGeoTypesPerChunk.clear(); newGeoTypesPerChunk2.clear();
  MCAuto<DataArrayInt> newGeoTypesEltIdsAllGather2=DataArrayInt::Aggregate(newGeoTypesPerChunk3); newGeoTypesPerChunk_bis.clear(); newGeoTypesPerChunk3.clear();
  MCAuto<DataArrayInt> diffVals=newGeoTypesEltIdsAllGather->getDifferentValues();
  MCAuto<DataArrayInt> renumEltIds=newGeoTypesEltIdsAllGather->buildPermArrPerLevel();
  //
  MCAuto<DataArrayInt> renumTupleIds=newGeoTypesPerChunk4->buildPermArrPerLevel();
  //
  MCAuto<DataArrayDouble> arrPart=arr->subArray(offset,offset+szTuples);
  arrPart->renumberInPlace(renumTupleIds->begin());
  arr->setPartOfValues1(arrPart,offset,offset+szTuples,1,0,arrPart->getNumberOfComponents(),1);
  bool ret=false;
  const int *idIt=diffVals->begin();
  std::list<const MEDFileFieldPerMeshPerTypePerDisc *> li(entriesOnSameDisc.begin(),entriesOnSameDisc.end());
  int offset2=0;
  for(int i=0;i<diffVals->getNumberOfTuples();i++,idIt++)
    {
      MCAuto<DataArrayInt> ids=newGeoTypesEltIdsAllGather->findIdsEqual(*idIt);
      MCAuto<DataArrayInt> subIds=newGeoTypesEltIdsAllGather2->selectByTupleId(ids->begin(),ids->end());
      int nbEntityElts=subIds->getNumberOfTuples();
      bool ret2;
      MCAuto<MEDFileFieldPerMeshPerTypePerDisc> eltToAdd=MEDFileFieldPerMeshPerTypePerDisc::
          NewObjectOnSameDiscThanPool(type,(INTERP_KERNEL::NormalizedCellType)newCode[3*(*idIt)],subIds,!subIds->isIota(newCode[3*(*idIt)+1]),nbi,
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
                                                                                                  bool &notInExisting)
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

////////////////////////////////////

MEDFileFieldPerMeshPerTypeCommon::~MEDFileFieldPerMeshPerTypeCommon()
{
}

void MEDFileFieldPerMeshPerTypeCommon::setFather(MEDFileFieldPerMesh *father)
{
  _father=father;
}

void MEDFileFieldPerMeshPerTypeCommon::accept(MEDFileFieldVisitor& visitor) const
{
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    if((*it).isNotNull())
      {
        visitor.newPerMeshPerTypePerDisc(*it);
      }
}

void MEDFileFieldPerMeshPerTypeCommon::deepCopyElements()
{
  std::size_t i=0;
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++,i++)
    {
      if((const MEDFileFieldPerMeshPerTypePerDisc *)*it)
        _field_pm_pt_pd[i]=(*it)->deepCopy(this);
    }
}

std::size_t MEDFileFieldPerMeshPerTypeCommon::getHeapMemorySizeWithoutChildren() const
{
  return _field_pm_pt_pd.capacity()*sizeof(MCAuto<MEDFileFieldPerMeshPerTypePerDisc>);
}

std::vector<const BigMemoryObject *> MEDFileFieldPerMeshPerTypeCommon::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    ret.push_back((const MEDFileFieldPerMeshPerTypePerDisc *)*it);
  return ret;
}

void MEDFileFieldPerMeshPerTypeCommon::assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldTemplate *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc)
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
void MEDFileFieldPerMeshPerTypeCommon::assignFieldProfile(bool isPflAlone, int& start, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, DataArrayInt *locIds, int nbOfEltsInWholeMesh, const MEDCouplingFieldTemplate *field, const DataArray *arr, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc)
{
  std::vector<int> pos=addNewEntryIfNecessary(field,idsInPfl);
  for(std::vector<int>::const_iterator it=pos.begin();it!=pos.end();it++)
    _field_pm_pt_pd[*it]->assignFieldProfile(isPflAlone,start,multiTypePfl,idsInPfl,locIds,nbOfEltsInWholeMesh,field,arr,mesh,glob,nasc);
}

void MEDFileFieldPerMeshPerTypeCommon::assignNodeFieldNoProfile(int& start, const MEDCouplingFieldTemplate *field, const DataArray *arr, MEDFileFieldGlobsReal& glob)
{
  _field_pm_pt_pd.resize(1);
  _field_pm_pt_pd[0]=MEDFileFieldPerMeshPerTypePerDisc::New(this,ON_NODES,-3);
  _field_pm_pt_pd[0]->assignNodeFieldNoProfile(start,field,arr,glob);
}

void MEDFileFieldPerMeshPerTypeCommon::assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldTemplate *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc)
{
  MCAuto<DataArrayInt> pfl2=pfl->deepCopy();
  if(!arr || !arr->isAllocated())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypeCommon::assignNodeFieldProfile : input array is null, or not allocated !");
  _field_pm_pt_pd.resize(1);
  _field_pm_pt_pd[0]=MEDFileFieldPerMeshPerTypePerDisc::New(this,ON_NODES,-3);
  _field_pm_pt_pd[0]->assignFieldProfile(true,start,pfl,pfl2,pfl2,-1,field,arr,0,glob,nasc);//mesh is not requested so 0 is send.
}

std::vector<int> MEDFileFieldPerMeshPerTypeCommon::addNewEntryIfNecessary(const MEDCouplingFieldTemplate *field, int offset, int nbOfCells)
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
      std::vector<int> ret(1,(int)sz);
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

std::vector<int> MEDFileFieldPerMeshPerTypeCommon::addNewEntryIfNecessaryGauss(const MEDCouplingFieldTemplate *field, int offset, int nbOfCells)
{
  const MEDCouplingFieldDiscretization *disc=field->getDiscretization();
  const MEDCouplingFieldDiscretizationGauss *disc2=dynamic_cast<const MEDCouplingFieldDiscretizationGauss *>(disc);
  if(!disc2)
    throw INTERP_KERNEL::Exception("addNewEntryIfNecessaryGauss : invalid call to this method ! Internal Error !");
  const DataArrayInt *da=disc2->getArrayOfDiscIds();
  if(!da)
    throw INTERP_KERNEL::Exception("addNewEntryIfNecessaryGauss (no profile) : no localization ids per cell array available ! The input Gauss node field is maybe invalid !");
  MCAuto<DataArrayInt> da2=da->selectByTupleIdSafeSlice(offset,offset+nbOfCells,1);
  MCAuto<DataArrayInt> retTmp=da2->getDifferentValues();
  if(retTmp->presenceOfValue(-1))
    throw INTERP_KERNEL::Exception("addNewEntryIfNecessaryGauss : some cells have no dicretization description !");
  std::vector<int> ret(retTmp->begin(),retTmp->end());
  return ret;
}

std::vector<int> MEDFileFieldPerMeshPerTypeCommon::addNewEntryIfNecessary(const MEDCouplingFieldTemplate *field, const DataArrayInt *subCells)
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

std::vector<int> MEDFileFieldPerMeshPerTypeCommon::addNewEntryIfNecessaryGauss(const MEDCouplingFieldTemplate *field, const DataArrayInt *subCells)
{
  const MEDCouplingFieldDiscretization *disc=field->getDiscretization();
  const MEDCouplingFieldDiscretizationGauss *disc2=dynamic_cast<const MEDCouplingFieldDiscretizationGauss *>(disc);
  if(!disc2)
    throw INTERP_KERNEL::Exception("addNewEntryIfNecessaryGauss : invalid call to this method ! Internal Error !");
  const DataArrayInt *da=disc2->getArrayOfDiscIds();
  if(!da)
    throw INTERP_KERNEL::Exception("addNewEntryIfNecessaryGauss : no localization ids per cell array available ! The input Gauss node field is maybe invalid !");
  MCAuto<DataArrayInt> da2=da->selectByTupleIdSafe(subCells->getConstPointer(),subCells->getConstPointer()+subCells->getNumberOfTuples());
  MCAuto<DataArrayInt> retTmp=da2->getDifferentValues();
  if(retTmp->presenceOfValue(-1))
    throw INTERP_KERNEL::Exception("addNewEntryIfNecessaryGauss : some cells have no dicretization description !");
  std::vector<int> ret(retTmp->begin(),retTmp->end());
  return ret;
}

const MEDFileFieldPerMesh *MEDFileFieldPerMeshPerTypeCommon::getFather() const
{
  return _father;
}

bool MEDFileFieldPerMeshPerTypeCommon::isUniqueLevel(int& dim) const
{
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(getGeoType()));
  int curDim((int)cm.getDimension());
  if(dim!=std::numeric_limits<int>::max())
    {
      if(dim!=curDim)
        return false;
    }
  else
    dim=curDim;
  return true;
}

void MEDFileFieldPerMeshPerTypeCommon::fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const
{
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    {
      (*it)->fillTypesOfFieldAvailable(types);
    }
}

void MEDFileFieldPerMeshPerTypeCommon::fillFieldSplitedByType(std::vector< std::pair<int,int> >& dads, std::vector<TypeOfField>& types, std::vector<std::string>& pfls, std::vector<std::string>& locs) const
{
  int sz=_field_pm_pt_pd.size();
  dads.resize(sz); types.resize(sz); pfls.resize(sz); locs.resize(sz);
  for(int i=0;i<sz;i++)
    {
      _field_pm_pt_pd[i]->getCoarseData(types[i],dads[i],pfls[i],locs[i]);
    }
}

int MEDFileFieldPerMeshPerTypeCommon::getIteration() const
{
  return _father->getIteration();
}

int MEDFileFieldPerMeshPerTypeCommon::getOrder() const
{
  return _father->getOrder();
}

double MEDFileFieldPerMeshPerTypeCommon::getTime() const
{
  return _father->getTime();
}

std::string MEDFileFieldPerMeshPerTypeCommon::getMeshName() const
{
  return _father->getMeshName();
}

void MEDFileFieldPerMeshPerTypeCommon::getSizes(int& globalSz, int& nbOfEntries) const
{
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    {
      globalSz+=(*it)->getNumberOfTuples();
    }
  nbOfEntries+=(int)_field_pm_pt_pd.size();
}

int MEDFileFieldPerMeshPerTypeCommon::getNumberOfComponents() const
{
  return _father->getNumberOfComponents();
}

bool MEDFileFieldPerMeshPerTypeCommon::presenceOfMultiDiscPerGeoType() const
{
  std::size_t nb(0);
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    {
      const MEDFileFieldPerMeshPerTypePerDisc *fmtd(*it);
      if(fmtd)
        nb++;
    }
  return nb>1;
}

void MEDFileFieldPerMeshPerTypeCommon::pushDiscretization(MEDFileFieldPerMeshPerTypePerDisc *disc)
{
  MCAuto<MEDFileFieldPerMeshPerTypePerDisc> elt;
  elt.takeRef(disc);
  _field_pm_pt_pd.push_back(elt);
}

DataArray *MEDFileFieldPerMeshPerTypeCommon::getOrCreateAndGetArray()
{
  return _father->getOrCreateAndGetArray();
}

const DataArray *MEDFileFieldPerMeshPerTypeCommon::getOrCreateAndGetArray() const
{
  const MEDFileFieldPerMesh *fath=_father;
  return fath->getOrCreateAndGetArray();
}

const std::vector<std::string>& MEDFileFieldPerMeshPerTypeCommon::getInfo() const
{
  return _father->getInfo();
}

std::vector<std::string> MEDFileFieldPerMeshPerTypeCommon::getPflsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it1=_field_pm_pt_pd.begin();it1!=_field_pm_pt_pd.end();it1++)
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

std::vector<std::string> MEDFileFieldPerMeshPerTypeCommon::getLocsReallyUsed() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it1=_field_pm_pt_pd.begin();it1!=_field_pm_pt_pd.end();it1++)
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

std::vector<std::string> MEDFileFieldPerMeshPerTypeCommon::getPflsReallyUsedMulti() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it1=_field_pm_pt_pd.begin();it1!=_field_pm_pt_pd.end();it1++)
    {
      std::string tmp=(*it1)->getProfile();
      if(!tmp.empty())
        ret.push_back(tmp);
    }
  return ret;
}

std::vector<std::string> MEDFileFieldPerMeshPerTypeCommon::getLocsReallyUsedMulti() const
{
  std::vector<std::string> ret;
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it1=_field_pm_pt_pd.begin();it1!=_field_pm_pt_pd.end();it1++)
    {
      std::string tmp=(*it1)->getLocalization();
      if(!tmp.empty() && tmp!=MED_GAUSS_ELNO)
        ret.push_back(tmp);
    }
  return ret;
}

void MEDFileFieldPerMeshPerTypeCommon::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::iterator it1=_field_pm_pt_pd.begin();it1!=_field_pm_pt_pd.end();it1++)
    (*it1)->changePflsRefsNamesGen(mapOfModif);
}

void MEDFileFieldPerMeshPerTypeCommon::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::iterator it1=_field_pm_pt_pd.begin();it1!=_field_pm_pt_pd.end();it1++)
    (*it1)->changeLocsRefsNamesGen(mapOfModif);
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypeCommon::getLeafGivenLocId(int locId)
{
  if(_field_pm_pt_pd.empty())
    {
      std::ostringstream oss; oss << "MEDFileFieldPerMeshPerTypeCommon::getLeafGivenLocId : no localizations for geotype \"" << getGeoTypeRepr() << "\" !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  if(locId>=0 && locId<(int)_field_pm_pt_pd.size())
    return _field_pm_pt_pd[locId];
  std::ostringstream oss2; oss2 << "MEDFileFieldPerMeshPerTypeCommon::getLeafGivenLocId : no such locId available (" << locId;
  oss2 << ") for geometric type \"" << getGeoTypeRepr() << "\" It should be in [0," << _field_pm_pt_pd.size() << ") !";
  throw INTERP_KERNEL::Exception(oss2.str().c_str());
  return static_cast<MEDFileFieldPerMeshPerTypePerDisc*>(0);
}

const MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMeshPerTypeCommon::getLeafGivenLocId(int locId) const
{
  if(_field_pm_pt_pd.empty())
    {
      std::ostringstream oss; oss << "MEDFileFieldPerMeshPerTypeCommon::getLeafGivenLocId : no localizations for geotype \"" << getGeoTypeRepr() << "\" !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  if(locId>=0 && locId<(int)_field_pm_pt_pd.size())
    return _field_pm_pt_pd[locId];
  std::ostringstream oss2; oss2 << "MEDFileFieldPerMeshPerTypeCommon::getLeafGivenLocId : no such locId available (" << locId;
  oss2 << ") for geometric type \"" << getGeoTypeRepr() << "\" It should be in [0," << _field_pm_pt_pd.size() << ") !";
  throw INTERP_KERNEL::Exception(oss2.str().c_str());
  return static_cast<const MEDFileFieldPerMeshPerTypePerDisc*>(0);
}

void MEDFileFieldPerMeshPerTypeCommon::fillValues(int& startEntryId, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const
{
  int i=0;
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++,i++)
    {
      (*it)->fillValues(i,startEntryId,entries);
    }
}

void MEDFileFieldPerMeshPerTypeCommon::setLeaves(const std::vector< MCAuto< MEDFileFieldPerMeshPerTypePerDisc > >& leaves)
{
  _field_pm_pt_pd=leaves;
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    (*it)->setFather(this);
}

/*!
 *  \param [in,out] globalNum a global numbering counter for the renumbering. 
 *  \param [out] its - list of pair (start,stop) kept
 *  \return bool - false if the type of field \a tof is not contained in \a this.
 */
bool MEDFileFieldPerMeshPerTypeCommon::keepOnlySpatialDiscretization(TypeOfField tof, int &globalNum, std::vector< std::pair<int,int> >& its)
{
  bool ret(false);
  std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> > newPmPtPd;
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    if((*it)->getType()==tof)
      {
        newPmPtPd.push_back(*it);
        std::pair<int,int> bgEnd; bgEnd.first=(*it)->getStart(); bgEnd.second=(*it)->getEnd();
        (*it)->setNewStart(globalNum);
        globalNum=(*it)->getEnd();
        its.push_back(bgEnd);
        ret=true;
      }
  if(ret)
    _field_pm_pt_pd=newPmPtPd;
  return ret;
}

/*!
 *  \param [in,out] globalNum a global numbering counter for the renumbering.
 *  \param [out] its - list of pair (start,stop) kept
 *  \return bool - false if the type of field \a tof is not contained in \a this.
 */
bool MEDFileFieldPerMeshPerTypeCommon::keepOnlyGaussDiscretization(std::size_t idOfDisc, int &globalNum, std::vector< std::pair<int,int> >& its)
{
  if(_field_pm_pt_pd.size()<=idOfDisc)
    return false;
  MCAuto<MEDFileFieldPerMeshPerTypePerDisc> elt(_field_pm_pt_pd[idOfDisc]);
  std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> > newPmPtPd(1,elt);
  std::pair<int,int> bgEnd; bgEnd.first=_field_pm_pt_pd[idOfDisc]->getStart(); bgEnd.second=_field_pm_pt_pd[idOfDisc]->getEnd();
  elt->setNewStart(globalNum);
  globalNum=elt->getEnd();
  its.push_back(bgEnd);
  _field_pm_pt_pd=newPmPtPd;
  return true;
}

void MEDFileFieldPerMeshPerTypeCommon::loadOnlyStructureOfDataRecursively(med_idt fid, int &start, const MEDFileFieldNameScope& nasc)
{
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    (*it)->loadOnlyStructureOfDataRecursively(fid,start,nasc);
}

void MEDFileFieldPerMeshPerTypeCommon::loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc)
{
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    (*it)->loadBigArray(fid,nasc);
}

void MEDFileFieldPerMeshPerTypeCommon::writeLL(med_idt fid, const MEDFileFieldNameScope& nasc) const
{
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    {
      (*it)->copyOptionsFrom(*this);
      (*it)->writeLL(fid,nasc);
    }
}

med_entity_type MEDFileFieldPerMeshPerTypeCommon::ConvertIntoMEDFileType(TypeOfField ikType, INTERP_KERNEL::NormalizedCellType ikGeoType, med_geometry_type& medfGeoType)
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
      throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypeCommon::ConvertIntoMEDFileType : unexpected entity type ! internal error");
  }
  return MED_UNDEF_ENTITY_TYPE;
}

//////////////////////////////////////////////////

MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerType::NewOnRead(med_idt fid, MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, const MEDFileFieldNameScope& nasc, const PartDefinition *pd)
{
  return new MEDFileFieldPerMeshPerType(fid,fath,type,geoType,nasc,pd);
}

MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerType::New(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType)
{
  return new MEDFileFieldPerMeshPerType(fath,geoType);
}

MEDFileFieldPerMeshPerType *MEDFileFieldPerMeshPerType::deepCopy(MEDFileFieldPerMesh *father) const
{
  MCAuto<MEDFileFieldPerMeshPerType> ret=new MEDFileFieldPerMeshPerType(*this);
  ret->setFather(father);
  ret->deepCopyElements();
  return ret.retn();
}

void MEDFileFieldPerMeshPerType::getFieldAtLevel(int meshDim, TypeOfField type, const MEDFileFieldGlobsReal *glob, std::vector< std::pair<int,int> >& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs, std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const
{
  if(_geo_type!=INTERP_KERNEL::NORM_ERROR)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_geo_type);
      if(meshDim!=(int)cm.getDimension())
        return ;
    }
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++)
    (*it)->getFieldAtLevel(type,glob,dads,pfls,locs,geoTypes);
}

INTERP_KERNEL::NormalizedCellType MEDFileFieldPerMeshPerType::getGeoType() const
{
  return _geo_type;
}

void MEDFileFieldPerMeshPerType::entriesForMEDfile(TypeOfField mct, med_geometry_type& gt, med_entity_type& ent) const
{
  ent=MEDFileFieldPerMeshPerTypeCommon::ConvertIntoMEDFileType(mct,_geo_type,gt);
}

void MEDFileFieldPerMeshPerType::getDimension(int& dim) const
{
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(_geo_type));
  int curDim((int)cm.getDimension());
  dim=std::max(dim,curDim);
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
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++,i++)
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

std::string MEDFileFieldPerMeshPerType::getGeoTypeRepr() const
{
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(_geo_type));
  return std::string(cm.getRepr());
}

MEDFileFieldPerMeshPerType::MEDFileFieldPerMeshPerType(MEDFileFieldPerMesh *father, INTERP_KERNEL::NormalizedCellType gt):MEDFileFieldPerMeshPerTypeCommon(father),_geo_type(gt)
{
}

MEDFileFieldPerMeshPerType::MEDFileFieldPerMeshPerType(med_idt fid, MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, const MEDFileFieldNameScope& nasc, const PartDefinition *pd):MEDFileFieldPerMeshPerTypeCommon(fath),_geo_type(geoType)
{
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> locName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  med_geometry_type mgeoti;
  med_entity_type menti(ConvertIntoMEDFileType(type,geoType,mgeoti));
  int nbProfiles(MEDfieldnProfile(fid,nasc.getName().c_str(),getIteration(),getOrder(),menti,mgeoti,pflName,locName));
  _field_pm_pt_pd.resize(nbProfiles);
  for(int i=0;i<nbProfiles;i++)
    {
      _field_pm_pt_pd[i]=MEDFileFieldPerMeshPerTypePerDisc::NewOnRead(this,type,i,pd);
    }
  if(type==ON_CELLS)
    {
      int nbProfiles2(MEDfieldnProfile(fid,nasc.getName().c_str(),getIteration(),getOrder(),MED_NODE_ELEMENT,mgeoti,pflName,locName));
      for(int i=0;i<nbProfiles2;i++)
        _field_pm_pt_pd.push_back(MEDFileFieldPerMeshPerTypePerDisc::NewOnRead(this,ON_GAUSS_NE,i,pd));
    }
  if(!_field_pm_pt_pd.empty() || type!=ON_CELLS)
    return ;
  // dark side of the force.
  {
    int nbProfiles1(MEDfieldnProfile(fid,nasc.getName().c_str(),getIteration(),getOrder(),MED_DESCENDING_FACE,mgeoti,pflName,locName));
    int nbProfiles2(MEDfieldnProfile(fid,nasc.getName().c_str(),getIteration(),getOrder(),MED_DESCENDING_EDGE,mgeoti,pflName,locName));
    if(nbProfiles1==0 && nbProfiles2==0)
      return ;// OK definitely nothing in field
    menti=nbProfiles1>=nbProfiles2?MED_DESCENDING_FACE:MED_DESCENDING_EDGE;//not enough words to describe the beauty
    nbProfiles=std::max(nbProfiles1,nbProfiles2);
    _field_pm_pt_pd.resize(nbProfiles);
    for(int i=0;i<nbProfiles;i++)
      _field_pm_pt_pd[i]=MEDFileFieldPerMeshPerTypePerDisc::NewOnRead(this,ON_CELLS,i,pd);
  }
}

MCAuto<MEDFileFieldPerMeshPerType> MEDFileFieldPerMeshPerType::Aggregate(int &start, const std::vector<std::pair<int,const MEDFileFieldPerMeshPerType *> >& pms, const std::vector< std::vector< std::pair<int,int> > >& dts, INTERP_KERNEL::NormalizedCellType gt, MEDFileFieldPerMesh *father, std::vector<std::pair< int, std::pair<int,int> > >& extractInfo)
{
  MCAuto<MEDFileFieldPerMeshPerType> ret(MEDFileFieldPerMeshPerType::New(father,gt));
  std::map<TypeOfField, std::vector< std::pair<int,const MEDFileFieldPerMeshPerTypePerDisc * > > > m;
  for(std::vector<std::pair<int,const MEDFileFieldPerMeshPerType *> >::const_iterator it=pms.begin();it!=pms.end();it++)
    {
      for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it2=(*it).second->_field_pm_pt_pd.begin();it2!=(*it).second->_field_pm_pt_pd.end();it2++)
        m[(*it2)->getType()].push_back(std::pair<int,const MEDFileFieldPerMeshPerTypePerDisc * >((*it).first,*it2));
    }
  for(std::map<TypeOfField, std::vector< std::pair<int,const MEDFileFieldPerMeshPerTypePerDisc * > > >::const_iterator it=m.begin();it!=m.end();it++)
    {
      MCAuto<MEDFileFieldPerMeshPerTypePerDisc> agg(MEDFileFieldPerMeshPerTypePerDisc::Aggregate(start,(*it).second,dts,(*it).first,ret,extractInfo));
      ret->_field_pm_pt_pd.push_back(agg);
    }
  return ret;
}

//////////////////////////////////////////////////

MEDFileFieldPerMeshPerTypeDyn *MEDFileFieldPerMeshPerTypeDyn::NewOnRead(med_idt fid, MEDFileFieldPerMesh *fath, const MEDFileEntities *entities, int idGT, const MEDFileFieldNameScope& nasc)
{
  if(!entities)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypeDyn::NewOnRead : null pointer !");
  const MEDFileAllStaticEntitiesPlusDyn *entities2(dynamic_cast<const MEDFileAllStaticEntitiesPlusDyn *>(entities));
  if(!entities2)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypeDyn::NewOnRead : invalid type of entities !");
  const MEDFileStructureElement *se(entities2->getWithGT(idGT));
  return new MEDFileFieldPerMeshPerTypeDyn(fid,fath,se,nasc);
}

MEDFileFieldPerMeshPerTypeDyn::MEDFileFieldPerMeshPerTypeDyn(med_idt fid, MEDFileFieldPerMesh *fath, const MEDFileStructureElement *se, const MEDFileFieldNameScope& nasc):MEDFileFieldPerMeshPerTypeCommon(fath)
{
  _se.takeRef(se);
  INTERP_KERNEL::AutoPtr<char> pflName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  INTERP_KERNEL::AutoPtr<char> locName=MEDLoaderBase::buildEmptyString(MED_NAME_SIZE);
  int nbProfiles(MEDfieldnProfile(fid,nasc.getName().c_str(),getIteration(),getOrder(),MED_STRUCT_ELEMENT,_se->getDynGT(),pflName,locName));
  _field_pm_pt_pd.resize(nbProfiles);
  for(int i=0;i<nbProfiles;i++)
    {
      _field_pm_pt_pd[i]=MEDFileFieldPerMeshPerTypePerDisc::NewOnRead(this,_se->getEntity(),i,NULL);
    }
}

int MEDFileFieldPerMeshPerTypeDyn::getDynGT() const
{
  return _se->getDynGT();
}

std::string MEDFileFieldPerMeshPerTypeDyn::getModelName() const
{
  return _se->getName();
}

void MEDFileFieldPerMeshPerTypeDyn::getDimension(int& dim) const
{
  throw INTERP_KERNEL::Exception("not implemented yet !");
}

void MEDFileFieldPerMeshPerTypeDyn::entriesForMEDfile(TypeOfField mct, med_geometry_type& gt, med_entity_type& ent) const
{
  gt=getDynGT();
  ent=MED_STRUCT_ELEMENT;
}

INTERP_KERNEL::NormalizedCellType MEDFileFieldPerMeshPerTypeDyn::getGeoType() const
{
  throw INTERP_KERNEL::Exception("not implemented yet !");
}

void MEDFileFieldPerMeshPerTypeDyn::simpleRepr(int bkOffset, std::ostream& oss, int id) const
{
  const char startLine[]="  ## ";
  std::string startLine2(bkOffset,' ');
  std::string startLine3(startLine2);
  startLine3+=startLine;
  oss << startLine3 << "Entry geometry type #" << id << " is lying on geometry STRUCTURE_ELEMENT type " << getDynGT() << "." << std::endl;
  oss << startLine3 << "Entry is defined on " <<  _field_pm_pt_pd.size() << " localizations." << std::endl;
  int i=0;
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++,i++)
    {
      if((*it).isNotNull())
        (*it)->simpleRepr(bkOffset,oss,i);
      else
        {
          oss << startLine2 << "    ## " << "Localization #" << i << " is empty !" << std::endl;
        }
    }
}

std::string MEDFileFieldPerMeshPerTypeDyn::getGeoTypeRepr() const
{
  throw INTERP_KERNEL::Exception("not implemented yet !");
}

MEDFileFieldPerMeshPerTypeDyn *MEDFileFieldPerMeshPerTypeDyn::deepCopy(MEDFileFieldPerMesh *father) const
{
  MCAuto<MEDFileFieldPerMeshPerTypeDyn> ret(new MEDFileFieldPerMeshPerTypeDyn(*this));
  ret->setFather(father);
  ret->deepCopyElements();
  return ret.retn();
}

void MEDFileFieldPerMeshPerTypeDyn::getFieldAtLevel(int meshDim, TypeOfField type, const MEDFileFieldGlobsReal *glob, std::vector< std::pair<int,int> >& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs, std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const
{
  throw INTERP_KERNEL::Exception("not implemented yet !");
}

//////////////////////////////////////////////////

MEDFileFieldPerMesh *MEDFileFieldPerMesh::NewOnRead(med_idt fid, MEDFileAnyTypeField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder, const MEDFileFieldNameScope& nasc, const MEDFileMesh *mm, const MEDFileEntities *entities)
{
  return new MEDFileFieldPerMesh(fid,fath,meshCsit,meshIteration,meshOrder,nasc,mm,entities);
}

MEDFileFieldPerMesh *MEDFileFieldPerMesh::New(MEDFileAnyTypeField1TSWithoutSDA *fath, const MEDCouplingMesh *mesh)
{
  return new MEDFileFieldPerMesh(fath,mesh);
}

std::size_t MEDFileFieldPerMesh::getHeapMemorySizeWithoutChildren() const
{
  return _field_pm_pt.capacity()*sizeof(MCAuto< MEDFileFieldPerMeshPerType >);
}

std::vector<const BigMemoryObject *> MEDFileFieldPerMesh::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    ret.push_back(*it);
  return ret;
}

MEDFileFieldPerMesh *MEDFileFieldPerMesh::deepCopy(MEDFileAnyTypeField1TSWithoutSDA *father) const
{
  MCAuto< MEDFileFieldPerMesh > ret=new MEDFileFieldPerMesh(*this);
  ret->_father=father;
  std::size_t i=0;
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++,i++)
    {
      if((*it).isNotNull())
        ret->_field_pm_pt[i]=(*it)->deepCopy((MEDFileFieldPerMesh *)(ret));
    }
  return ret.retn();
}

void MEDFileFieldPerMesh::simpleRepr(int bkOffset, std::ostream& oss, int id) const
{
  std::string startLine(bkOffset,' ');
  oss << startLine << "## Field part (" << id << ") lying on mesh \"" << getMeshName() << "\", Mesh iteration=" << _mesh_iteration << ". Mesh order=" << _mesh_order << "." << std::endl;
  oss << startLine << "## Field is defined on " << _field_pm_pt.size() << " types." << std::endl;
  int i=0;
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++,i++)
    {
      if((*it).isNotNull())
        (*it)->simpleRepr(bkOffset,oss,i);
      else
        {
          oss << startLine << "  ## Entry geometry type #" << i << " is empty !" << std::endl;
        }
    }
}

void MEDFileFieldPerMesh::copyTinyInfoFrom(const MEDCouplingMesh *mesh)
{
  mesh->getTime(_mesh_iteration,_mesh_order);
}

void MEDFileFieldPerMesh::assignFieldNoProfileNoRenum(int& start, const std::vector<int>& code, const MEDCouplingFieldTemplate *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc)
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
void MEDFileFieldPerMesh::assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const std::vector<int>& code, const std::vector<int>& code2, const std::vector<DataArrayInt *>& idsInPflPerType, const std::vector<DataArrayInt *>& idsPerType, const MEDCouplingFieldTemplate *field, const DataArray *arr, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc)
{
  int nbOfTypes(code.size()/3);
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
      _field_pm_pt[pos]->assignFieldProfile(nbOfTypes==1,start,multiTypePfl,idsInPflPerType[i],pfl,code2[3*found+1],field,arr,mesh,glob,nasc);
    }
}

void MEDFileFieldPerMesh::assignNodeFieldNoProfile(int& start, const MEDCouplingFieldTemplate *field, const DataArray *arr, MEDFileFieldGlobsReal& glob)
{
  int pos=addNewEntryIfNecessary(INTERP_KERNEL::NORM_ERROR);
  _field_pm_pt[pos]->assignNodeFieldNoProfile(start,field,arr,glob);
}

void MEDFileFieldPerMesh::assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldTemplate *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc)
{
  int pos=addNewEntryIfNecessary(INTERP_KERNEL::NORM_ERROR);
  _field_pm_pt[pos]->assignNodeFieldProfile(start,pfl,field,arr,glob,nasc);
}

void MEDFileFieldPerMesh::loadOnlyStructureOfDataRecursively(med_idt fid, int& start, const MEDFileFieldNameScope& nasc)
{
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->loadOnlyStructureOfDataRecursively(fid,start,nasc);
}

void MEDFileFieldPerMesh::loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc)
{
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->loadBigArraysRecursively(fid,nasc);
}

void MEDFileFieldPerMesh::writeLL(med_idt fid, const MEDFileFieldNameScope& nasc) const
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
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->getDimension(dim);
}

bool MEDFileFieldPerMesh::isUniqueLevel(int& dim) const
{
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    if(!(*it)->isUniqueLevel(dim))
      return false;
  return true;
}

void MEDFileFieldPerMesh::fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const
{
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->fillTypesOfFieldAvailable(types);
}

std::vector< std::vector< std::pair<int,int> > > MEDFileFieldPerMesh::getFieldSplitedByType(std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> > & locs) const
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

std::string MEDFileFieldPerMesh::getMeshName() const
{
  return _father->getMeshName();
}

void MEDFileFieldPerMesh::setMeshName(const std::string& meshName)
{
  _father->setMeshName(meshName);
}

bool MEDFileFieldPerMesh::presenceOfMultiDiscPerGeoType() const
{
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      if((*it).isNull())
        continue;
      if((*it)->presenceOfMultiDiscPerGeoType())
        return true;
    }
  return false;
}

bool MEDFileFieldPerMesh::presenceOfStructureElements() const
{
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    if((*it).isNotNull())
      {
        const MEDFileFieldPerMeshPerTypeDyn *pt(dynamic_cast<const MEDFileFieldPerMeshPerTypeDyn *>((const MEDFileFieldPerMeshPerTypeCommon *)*it));
        if(pt)
          return true;
      }
  return false;
}

bool MEDFileFieldPerMesh::onlyStructureElements() const
{
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    if((*it).isNotNull())
      {
        const MEDFileFieldPerMeshPerTypeDyn *pt(dynamic_cast<const MEDFileFieldPerMeshPerTypeDyn *>((const MEDFileFieldPerMeshPerTypeCommon *)*it));
        if(!pt)
          return false;
      }
  return true;
}

void MEDFileFieldPerMesh::killStructureElements()
{
  std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > > res;
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      if((*it).isNotNull())
        {
          const MEDFileFieldPerMeshPerTypeDyn *pt(dynamic_cast<const MEDFileFieldPerMeshPerTypeDyn *>((const MEDFileFieldPerMeshPerTypeCommon *)*it));
          if(!pt)
            res.push_back(*it);
        }
    }
  _field_pm_pt=res;
}

void MEDFileFieldPerMesh::keepOnlyStructureElements()
{
  std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > > res;
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      if((*it).isNotNull())
        {
          const MEDFileFieldPerMeshPerTypeDyn *pt(dynamic_cast<const MEDFileFieldPerMeshPerTypeDyn *>((const MEDFileFieldPerMeshPerTypeCommon *)*it));
          if(pt)
            res.push_back(*it);
        }
    }
  _field_pm_pt=res;
}

void MEDFileFieldPerMesh::keepOnlyOnSE(const std::string& seName)
{
  std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > > res;
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      if((*it).isNotNull())
        {
          const MEDFileFieldPerMeshPerTypeDyn *pt(dynamic_cast<const MEDFileFieldPerMeshPerTypeDyn *>((const MEDFileFieldPerMeshPerTypeCommon *)*it));
          if(!pt)
            throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::keepOnlyOnSE : presence of non SE !");
          if(pt->getModelName()==seName)
            res.push_back(*it);
        }
    }
  _field_pm_pt=res;
}

void MEDFileFieldPerMesh::getMeshSENames(std::vector< std::pair<std::string,std::string> >& ps) const
{
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      if((*it).isNotNull())
        {
          const MEDFileFieldPerMeshPerTypeDyn *pt(dynamic_cast<const MEDFileFieldPerMeshPerTypeDyn *>((const MEDFileFieldPerMeshPerTypeCommon *)*it));
          if(pt)
            {
              ps.push_back(std::pair<std::string,std::string>(getMeshName(),pt->getModelName()));
            }
          else
            throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::getMeshSENames : presence of a non structure element part !");
        }
    }
}

DataArray *MEDFileFieldPerMesh::getOrCreateAndGetArray()
{
  if(!_father)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::getOrCreateAndGetArray : no father ! internal error !");
  return _father->getOrCreateAndGetArray();
}

const DataArray *MEDFileFieldPerMesh::getOrCreateAndGetArray() const
{
  if(!_father)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::getOrCreateAndGetArray : no father ! internal error !");
  return _father->getOrCreateAndGetArray();
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
int MEDFileFieldPerMesh::ComputeNbOfElems(const MEDFileFieldGlobsReal *glob, TypeOfField type, const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs)
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
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
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
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
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
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
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
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsedMulti();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

bool MEDFileFieldPerMesh::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab)
{
  for(std::vector< std::pair<std::string,std::string> >::const_iterator it=modifTab.begin();it!=modifTab.end();it++)
    {
      if((*it).first==getMeshName())
        {
          setMeshName((*it).second);
          return true;
        }
    }
  return false;
}

void MEDFileFieldPerMesh::convertMedBallIntoClassic()
{
  if(_field_pm_pt.size()!=1)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::convertMedBallIntoClassic : Only managed for single mesh !");
  if(_field_pm_pt[0].isNull())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::convertMedBallIntoClassic : null pointer !");
  MEDFileFieldPerMeshPerTypeDyn *pt(dynamic_cast<MEDFileFieldPerMeshPerTypeDyn *>((MEDFileFieldPerMeshPerTypeCommon *)_field_pm_pt[0]));
  if(!pt)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::convertMedBallIntoClassic : this is expected to be marked as structure element !");
  if(pt->getNumberOfLoc()!=1)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::convertMedBallIntoClassic : only one loc managed !");
  const MEDFileFieldPerMeshPerTypePerDisc *disc(pt->getLeafGivenLocId(0));
  if(!disc)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::convertMedBallIntoClassic : internal error !");
  MCAuto<MEDFileFieldPerMeshPerTypePerDisc> disc2(MEDFileFieldPerMeshPerTypePerDisc::New(*disc));
  disc2->setType(ON_NODES);
  MCAuto<MEDFileFieldPerMeshPerType> pt2(MEDFileFieldPerMeshPerType::New(this,INTERP_KERNEL::NORM_ERROR));
  disc2->setFather(pt2);
  pt2->setFather(this);
  pt2->pushDiscretization(disc2);
  _field_pm_pt[0]=DynamicCast<MEDFileFieldPerMeshPerType,MEDFileFieldPerMeshPerTypeCommon>(pt2);
}

bool MEDFileFieldPerMesh::renumberEntitiesLyingOnMesh(const std::string& meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N,
                                                      MEDFileFieldGlobsReal& glob)
{
  if(getMeshName()!=meshName)
    return false;
  std::set<INTERP_KERNEL::NormalizedCellType> typesToKeep;
  for(std::size_t i=0;i<oldCode.size()/3;i++) typesToKeep.insert((INTERP_KERNEL::NormalizedCellType)oldCode[3*i]);
  std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > > entries;
  std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> entriesKept;
  std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> otherEntries;
  getUndergroundDataArrayExt(entries);
  DataArray *arr0(getOrCreateAndGetArray());//tony
  if(!arr0)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::renumberEntitiesLyingOnMesh : DataArray storing values of field is null !");
  DataArrayDouble *arr(dynamic_cast<DataArrayDouble *>(arr0));//tony
  if(!arr0)
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::renumberEntitiesLyingOnMesh : DataArray storing values is double ! Not managed for the moment !");
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
  MCAuto<DataArrayInt> renumDefrag=DataArrayInt::New(); renumDefrag->alloc(arr->getNumberOfTuples(),1); renumDefrag->fillWithZero();
  ////////////////////
  MCAuto<DataArrayInt> explicitIdsOldInMesh=DataArrayInt::New(); explicitIdsOldInMesh->alloc(sz,1);//sz is a majorant of the real size. A realloc will be done after
  int *workI2=explicitIdsOldInMesh->getPointer();
  int sz1=0,sz2=0,sid=1;
  std::vector< std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> > entriesKeptML=MEDFileFieldPerMeshPerTypePerDisc::SplitPerDiscretization(entriesKept);
  // std::vector<int> tupleIdOfStartOfNewChuncksV(entriesKeptML.size());
  for(std::vector< std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> >::const_iterator itL1=entriesKeptML.begin();itL1!=entriesKeptML.end();itL1++,sid++)
    {
      //  tupleIdOfStartOfNewChuncksV[sid-1]=sz2;
      MCAuto<DataArrayInt> explicitIdsOldInArr=DataArrayInt::New(); explicitIdsOldInArr->alloc(sz,1);
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
  MCAuto<DataArrayInt> permArrDefrag=renumDefrag->buildPermArrPerLevel(); renumDefrag=0;
  // perform redispatching of non concerned MEDFileFieldPerMeshPerTypePerDisc
  std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> > otherEntriesNew;
  for(std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>::const_iterator it=otherEntries.begin();it!=otherEntries.end();it++)
    {
      otherEntriesNew.push_back(MEDFileFieldPerMeshPerTypePerDisc::New(*(*it)));
      otherEntriesNew.back()->setNewStart(permArrDefrag->getIJ((*it)->getStart(),0));
      otherEntriesNew.back()->setLocId((*it)->getGeoType());
    }
  std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> > entriesKeptNew;
  std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> entriesKeptNew2;
  for(std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>::const_iterator it=entriesKept.begin();it!=entriesKept.end();it++)
    {
      MCAuto<MEDFileFieldPerMeshPerTypePerDisc> elt=MEDFileFieldPerMeshPerTypePerDisc::New(*(*it));
      int newStart=elt->getLocId();
      elt->setLocId((*it)->getGeoType());
      elt->setNewStart(newStart);
      elt->_tmp_work1=permArrDefrag->getIJ(elt->_tmp_work1,0);
      entriesKeptNew.push_back(elt);
      entriesKeptNew2.push_back(elt);
    }
  MCAuto<DataArrayDouble> arr2=arr->renumber(permArrDefrag->getConstPointer());
  // perform redispatching of concerned MEDFileFieldPerMeshPerTypePerDisc -> values are in arr2
  MCAuto<DataArrayInt> explicitIdsNewInMesh=renumO2N->selectByTupleId(explicitIdsOldInMesh->begin(),explicitIdsOldInMesh->end());
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
  arr->deepCopyFrom(*arr2);
  return true;
}

/*!
 * \param [in,out] globalNum a global numbering counter for the renumbering.
 * \param [out] its - list of pair (start,stop) kept
 */
void MEDFileFieldPerMesh::keepOnlySpatialDiscretization(TypeOfField tof, int &globalNum, std::vector< std::pair<int,int> >& its)
{
  std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > > ret;
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      std::vector< std::pair<int,int> > its2;
      if((*it)->keepOnlySpatialDiscretization(tof,globalNum,its2))
        {
          ret.push_back(*it);
          its.insert(its.end(),its2.begin(),its2.end());
        }
    }
  _field_pm_pt=ret;
}

/*!
 * \param [in,out] globalNum a global numbering counter for the renumbering.
 * \param [out] its - list of pair (start,stop) kept
 */
void MEDFileFieldPerMesh::keepOnlyGaussDiscretization(std::size_t idOfDisc, int &globalNum, std::vector< std::pair<int,int> >& its)
{
  std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > > ret;
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      std::vector< std::pair<int,int> > its2;
      if((*it)->keepOnlyGaussDiscretization(idOfDisc,globalNum,its2))
        {
          ret.push_back(*it);
          its.insert(its.end(),its2.begin(),its2.end());
        }
    }
  _field_pm_pt=ret;
}

void MEDFileFieldPerMesh::assignNewLeaves(const std::vector< MCAuto< MEDFileFieldPerMeshPerTypePerDisc > >& leaves)
{
  std::map<INTERP_KERNEL::NormalizedCellType,std::vector< MCAuto< MEDFileFieldPerMeshPerTypePerDisc> > > types;
  for( std::vector< MCAuto< MEDFileFieldPerMeshPerTypePerDisc > >::const_iterator it=leaves.begin();it!=leaves.end();it++)
    types[(INTERP_KERNEL::NormalizedCellType)(*it)->getLocId()].push_back(*it);
  //
  std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > > fieldPmPt(types.size());
  std::map<INTERP_KERNEL::NormalizedCellType,std::vector< MCAuto< MEDFileFieldPerMeshPerTypePerDisc> > >::const_iterator it1=types.begin();
  std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::iterator it2=fieldPmPt.begin();
  for(;it1!=types.end();it1++,it2++)
    {
      MCAuto<MEDFileFieldPerMeshPerType> elt=MEDFileFieldPerMeshPerType::New(this,(INTERP_KERNEL::NormalizedCellType)((*it1).second[0]->getLocId()));
      elt->setLeaves((*it1).second);
      MCAuto<MEDFileFieldPerMeshPerTypeCommon> elt2(DynamicCast<MEDFileFieldPerMeshPerType,MEDFileFieldPerMeshPerTypeCommon>(elt));
      *it2=elt2;
    }
  _field_pm_pt=fieldPmPt;
}

void MEDFileFieldPerMesh::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->changePflsRefsNamesGen(mapOfModif);
}

void MEDFileFieldPerMesh::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->changeLocsRefsNamesGen(mapOfModif);
}

/*!
 * \param [in] mesh is the whole mesh
 */
MEDCouplingFieldDouble *MEDFileFieldPerMesh::getFieldOnMeshAtLevel(TypeOfField type, const MEDFileFieldGlobsReal *glob, const MEDCouplingMesh *mesh, bool& isPfl, MCAuto<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const
{
  if(_field_pm_pt.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::getFieldOnMeshAtLevel : no types field set !");
  //
  std::vector< std::pair<int,int> > dads;
  std::vector<const DataArrayInt *> pfls;
  std::vector<DataArrayInt *> notNullPflsPerGeoType;
  std::vector<int> locs,code;
  std::vector<INTERP_KERNEL::NormalizedCellType> geoTypes;
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->getFieldAtLevel(mesh->getMeshDimension(),type,glob,dads,pfls,locs,geoTypes);
  // Sort by types
  SortArraysPerType(glob,type,geoTypes,dads,pfls,locs,code,notNullPflsPerGeoType);
  if(code.empty())
    {
      std::ostringstream oss; oss << "MEDFileFieldPerMesh::getFieldOnMeshAtLevel : " << "The field \"" << nasc.getName() << "\" exists but not with such spatial discretization or such dimension specified !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  //
  std::vector< MCAuto<DataArrayInt> > notNullPflsPerGeoType2(notNullPflsPerGeoType.begin(),notNullPflsPerGeoType.end());
  std::vector< const DataArrayInt *> notNullPflsPerGeoType3(notNullPflsPerGeoType.begin(),notNullPflsPerGeoType.end());
  if(type!=ON_NODES)
    {
      DataArrayInt *arr=mesh->checkTypeConsistencyAndContig(code,notNullPflsPerGeoType3);
      if(!arr)
        return finishField(type,glob,dads,locs,mesh,isPfl,arrOut,nasc);
      else
        {
          MCAuto<DataArrayInt> arr2(arr);
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
              throw INTERP_KERNEL::Exception(oss.str());
            }
          return finishField(type,glob,dads,locs,mesh,isPfl,arrOut,nasc);
        }
      else
        return finishFieldNode2(glob,dads,locs,mesh,notNullPflsPerGeoType3[0],isPfl,arrOut,nasc);
    }
}

DataArray *MEDFileFieldPerMesh::getFieldOnMeshAtLevelWithPfl(TypeOfField type, const MEDCouplingMesh *mesh, DataArrayInt *&pfl, const MEDFileFieldGlobsReal *glob, const MEDFileFieldNameScope& nasc) const
{
  if(_field_pm_pt.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::getFieldOnMeshAtLevel : no types field set !");
  //
  std::vector<std::pair<int,int> > dads;
  std::vector<const DataArrayInt *> pfls;
  std::vector<DataArrayInt *> notNullPflsPerGeoType;
  std::vector<int> locs,code;
  std::vector<INTERP_KERNEL::NormalizedCellType> geoTypes;
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    (*it)->getFieldAtLevel(mesh->getMeshDimension(),type,glob,dads,pfls,locs,geoTypes);
  // Sort by types
  SortArraysPerType(glob,type,geoTypes,dads,pfls,locs,code,notNullPflsPerGeoType);
  if(code.empty())
    {
      std::ostringstream oss; oss << "MEDFileFieldPerMesh::getFieldOnMeshAtLevelWithPfl : " << "The field \"" << nasc.getName() << "\" exists but not with such spatial discretization or such dimension specified !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  std::vector< MCAuto<DataArrayInt> > notNullPflsPerGeoType2(notNullPflsPerGeoType.begin(),notNullPflsPerGeoType.end());
  std::vector< const DataArrayInt *> notNullPflsPerGeoType3(notNullPflsPerGeoType.begin(),notNullPflsPerGeoType.end());
  if(type!=ON_NODES)
    {
      MCAuto<DataArrayInt> arr=mesh->checkTypeConsistencyAndContig(code,notNullPflsPerGeoType3);
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
              throw INTERP_KERNEL::Exception(oss.str());
            }
        }
      return finishField4(dads,code[2]==-1?0:notNullPflsPerGeoType3[0],mesh->getNumberOfNodes(),pfl);
    }
  //
  return 0;
}

void MEDFileFieldPerMesh::accept(MEDFileFieldVisitor& visitor) const
{
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    if((*it).isNotNull())
      {
        visitor.newPerMeshPerTypeEntry(*it);
        (*it)->accept(visitor);
        visitor.endPerMeshPerTypeEntry(*it);
      }
}

void MEDFileFieldPerMesh::getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const
{
  int globalSz=0;
  int nbOfEntries=0;
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      (*it)->getSizes(globalSz,nbOfEntries);
    }
  entries.resize(nbOfEntries);
  nbOfEntries=0;
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      (*it)->fillValues(nbOfEntries,entries);
    }
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMesh::getLeafGivenTypeAndLocId(INTERP_KERNEL::NormalizedCellType typ, int locId)
{
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      if((*it)->getGeoType()==typ)
        return (*it)->getLeafGivenLocId(locId);
    }
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(typ);
  std::ostringstream oss; oss << "MEDFileFieldPerMesh::getLeafGivenTypeAndLocId : no such geometric type \"" << cm.getRepr() << "\" in this !" << std::endl;
  oss << "Possiblities are : ";
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      const INTERP_KERNEL::CellModel& cm2=INTERP_KERNEL::CellModel::GetCellModel((*it)->getGeoType());
      oss << "\"" << cm2.getRepr() << "\", ";
    }
  throw INTERP_KERNEL::Exception(oss.str());
}

const MEDFileFieldPerMeshPerTypePerDisc *MEDFileFieldPerMesh::getLeafGivenTypeAndLocId(INTERP_KERNEL::NormalizedCellType typ, int locId) const
{
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      if((*it)->getGeoType()==typ)
        return (*it)->getLeafGivenLocId(locId);
    }
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(typ);
  std::ostringstream oss; oss << "MEDFileFieldPerMesh::getLeafGivenTypeAndLocId : no such geometric type \"" << cm.getRepr() << "\" in this !" << std::endl;
  oss << "Possiblities are : ";
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++)
    {
      const INTERP_KERNEL::CellModel& cm2=INTERP_KERNEL::CellModel::GetCellModel((*it)->getGeoType());
      oss << "\"" << cm2.getRepr() << "\", ";
    }
  throw INTERP_KERNEL::Exception(oss.str());
}

/*!
 * \param [in,out] start - Integer that gives the current position in the final aggregated array
 * \param [in] pms - list of elements to aggregate. integer gives the mesh id 
 * \param [in] dts - (Distribution of types) = level 1 : meshes to aggregate. Level 2 : all geo type. Level 3 pair specifying geo type and number of elem in geotype.
 * \param [out] extractInfo - Gives information about the where the data comes from. It is a vector of triplet. First element in the triplet the mesh pos. The 2nd one the start pos. The 3rd the end pos.
 */
MCAuto<MEDFileFieldPerMeshPerTypePerDisc> MEDFileFieldPerMeshPerTypePerDisc::Aggregate(int &start, const std::vector< std::pair<int,const MEDFileFieldPerMeshPerTypePerDisc *> >& pms, const std::vector< std::vector< std::pair<int,int> > >& dts, TypeOfField tof, MEDFileFieldPerMeshPerType *father, std::vector<std::pair< int, std::pair<int,int> > >& extractInfo)
{
  MCAuto<MEDFileFieldPerMeshPerTypePerDisc> ret(new MEDFileFieldPerMeshPerTypePerDisc(father,tof));
  if(pms.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::Aggregate : empty input vector !");
  for(std::vector<std::pair<int,const MEDFileFieldPerMeshPerTypePerDisc *> >::const_iterator it=pms.begin();it!=pms.end();it++)
    {
      if(!(*it).second)
        throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::Aggregate : presence of null pointer !");
      if(!(*it).second->getProfile().empty())
        throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::Aggregate : not implemented yet for profiles !");
      if(!(*it).second->getLocalization().empty())
        throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::Aggregate : not implemented yet for gauss pts !");
    }
  INTERP_KERNEL::NormalizedCellType gt(pms[0].second->getGeoType());
  std::size_t i(0);
  std::vector< std::pair<int,int> > filteredDTS;
  for(std::vector< std::vector< std::pair<int,int> > >::const_iterator it=dts.begin();it!=dts.end();it++,i++)
    for(std::vector< std::pair<int,int> >::const_iterator it2=(*it).begin();it2!=(*it).end();it2++)
      if((*it2).first==gt)
        filteredDTS.push_back(std::pair<int,int>(i,(*it2).second));
  if(pms.size()!=filteredDTS.size())
    throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::Aggregate : not implemented yet for generated profiles !");
  std::vector<std::pair<int,const MEDFileFieldPerMeshPerTypePerDisc *> >::const_iterator it1(pms.begin());
  std::vector< std::pair<int,int> >::const_iterator it2(filteredDTS.begin());
  int zeStart(start),nval(0);
  for(;it1!=pms.end();it1++,it2++)
    {
      if((*it1).first!=(*it2).first)
        throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypePerDisc::Aggregate : not implemented yet for generated profiles 2 !");
      int s1((*it1).second->getStart()),e1((*it1).second->getEnd());
      extractInfo.push_back(std::pair<int, std::pair<int,int> >((*it1).first,std::pair<int,int>(s1,e1)));
      start+=e1-s1;
      nval+=((*it1).second)->getNumberOfVals();
    }
  ret->_start=zeStart; ret->_end=start; ret->_nval=nval;
  return ret;
}

MCAuto<MEDFileFieldPerMesh> MEDFileFieldPerMesh::Aggregate(int &start, const std::vector<const MEDFileFieldPerMesh *>& pms, const std::vector< std::vector< std::pair<int,int> > >& dts, MEDFileAnyTypeField1TSWithoutSDA *father, std::vector<std::pair< int, std::pair<int,int> > >& extractInfo)
{
  MCAuto<MEDFileFieldPerMesh> ret(new MEDFileFieldPerMesh(father,pms[0]->getMeshName(),pms[0]->getMeshIteration(),pms[0]->getMeshOrder()));
  std::map<INTERP_KERNEL::NormalizedCellType, std::vector< std::pair<int,const MEDFileFieldPerMeshPerType *> > > m;
  std::size_t i(0);
  for(std::vector<const MEDFileFieldPerMesh *>::const_iterator it=pms.begin();it!=pms.end();it++,i++)
    {
      const std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >& v((*it)->_field_pm_pt);
      for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::const_iterator it2=v.begin();it2!=v.end();it2++)
        {
          INTERP_KERNEL::NormalizedCellType gt((*it2)->getGeoType());
          const MEDFileFieldPerMeshPerType *elt(dynamic_cast<const MEDFileFieldPerMeshPerType *>((const MEDFileFieldPerMeshPerTypeCommon *)(*it2)));
          if(!elt)
            throw INTERP_KERNEL::Exception("MEDFileFieldPerMesh::Aggregate : not managed for structelement !");
          m[gt].push_back(std::pair<int,const MEDFileFieldPerMeshPerType *>(i,elt));
        }
    }
  for(std::map<INTERP_KERNEL::NormalizedCellType, std::vector< std::pair<int,const MEDFileFieldPerMeshPerType *> > >::const_iterator it=m.begin();it!=m.end();it++)
    {
      MCAuto<MEDFileFieldPerMeshPerType> agg(MEDFileFieldPerMeshPerType::Aggregate(start,(*it).second,dts,(*it).first,ret,extractInfo));
      MCAuto<MEDFileFieldPerMeshPerTypeCommon> agg2(DynamicCast<MEDFileFieldPerMeshPerType,MEDFileFieldPerMeshPerTypeCommon>(agg));
      ret->_field_pm_pt.push_back(agg2);
    }
  return ret;
}

int MEDFileFieldPerMesh::addNewEntryIfNecessary(INTERP_KERNEL::NormalizedCellType type)
{
  int i=0;
  int pos=std::distance(typmai2,std::find(typmai2,typmai2+MED_N_CELL_FIXED_GEO,type));
  std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::iterator it2=_field_pm_pt.begin();
  for(std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > >::iterator it=_field_pm_pt.begin();it!=_field_pm_pt.end();it++,i++)
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
                                                         const MEDCouplingMesh *mesh, bool& isPfl, MCAuto<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const
{
  isPfl=false;
  MCAuto<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(type,ONE_TIME);
  ret->setMesh(mesh); ret->setName(nasc.getName().c_str()); ret->setTime(getTime(),getIteration(),getOrder()); ret->setTimeUnit(nasc.getDtUnit().c_str());
  MCAuto<DataArray> da=getOrCreateAndGetArray()->selectByTupleRanges(dads);
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
          MCAuto<DataArrayInt> di=DataArrayInt::New();
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
                                                          const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl, MCAuto<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const
{
  if(da->isIota(mesh->getNumberOfCells()))
    return finishField(type,glob,dads,locs,mesh,isPfl,arrOut,nasc);
  MCAuto<MEDCouplingMesh> m2=mesh->buildPart(da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
  m2->setName(mesh->getName().c_str());
  MCAuto<MEDCouplingFieldDouble> ret=finishField(type,glob,dads,locs,m2,isPfl,arrOut,nasc);
  isPfl=true;
  return ret.retn();
}

/*!
 * This method is the complement of MEDFileFieldPerMesh::finishField2 method except that this method works for node profiles.
 */
MEDCouplingFieldDouble *MEDFileFieldPerMesh::finishFieldNode2(const MEDFileFieldGlobsReal *glob,
                                                              const std::vector<std::pair<int,int> >& dads, const std::vector<int>& locs,
                                                              const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl, MCAuto<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const
{
  if(da->isIota(mesh->getNumberOfNodes()))
    return finishField(ON_NODES,glob,dads,locs,mesh,isPfl,arrOut,nasc);
  // Treatment of particular case where nodal field on pfl is requested with a meshDimRelToMax=1.
  const MEDCouplingUMesh *meshu=dynamic_cast<const MEDCouplingUMesh *>(mesh);
  if(meshu)
    {
      if(meshu->getNodalConnectivity()==0)
        {
          MCAuto<MEDCouplingFieldDouble> ret=finishField(ON_CELLS,glob,dads,locs,mesh,isPfl,arrOut,nasc);
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
  MCAuto<MEDCouplingFieldDouble> ret=finishField(ON_NODES,glob,dads,locs,mesh,isPfl,arrOut,nasc);
  isPfl=true;
  DataArrayInt *arr2=0;
  MCAuto<DataArrayInt> cellIds=mesh->getCellIdsFullyIncludedInNodeIds(da->getConstPointer(),da->getConstPointer()+da->getNbOfElems());
  MCAuto<MEDCouplingMesh> mesh2=mesh->buildPartAndReduceNodes(cellIds->getConstPointer(),cellIds->getConstPointer()+cellIds->getNbOfElems(),arr2);
  MCAuto<DataArrayInt> arr3(arr2);
  int nnodes=mesh2->getNumberOfNodes();
  if(nnodes==(int)da->getNbOfElems())
    {
      MCAuto<DataArrayInt> da3=da->transformWithIndArrR(arr2->begin(),arr2->end());
      arrOut->renumberInPlace(da3->getConstPointer());
      mesh2->setName(mesh->getName().c_str());
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
      throw INTERP_KERNEL::Exception(oss.str());
    }
  return 0;
}

/*!
 * This method is the most light method of field retrieving.
 */
DataArray *MEDFileFieldPerMesh::finishField4(const std::vector<std::pair<int,int> >& dads, const DataArrayInt *pflIn, int nbOfElems, DataArrayInt *&pflOut) const
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
  MCAuto<DataArrayInt> safePfl(pflOut);
  MCAuto<DataArray> da=getOrCreateAndGetArray()->selectByTupleRanges(dads);
  const std::vector<std::string>& infos=getInfo();
  int nbOfComp=infos.size();
  for(int i=0;i<nbOfComp;i++)
    da->setInfoOnComponent(i,infos[i].c_str());
  safePfl->incrRef();
  return da.retn();
}


/// @cond INTERNAL

class MFFPMIter
{
public:
  static MFFPMIter *NewCell(const MEDFileEntities *entities);
  static bool IsPresenceOfNode(const MEDFileEntities *entities);
  virtual ~MFFPMIter() { }
  virtual void begin() = 0;
  virtual bool finished() const = 0;
  virtual void next() = 0;
  virtual int current() const = 0;
};

class MFFPMIterSimple : public MFFPMIter
{
public:
  MFFPMIterSimple():_pos(0) { }
  void begin() { _pos=0; }
  bool finished() const { return _pos>=MED_N_CELL_FIXED_GEO; }
  void next() { _pos++; }
  int current() const { return _pos; }
private:
  int _pos;
};

class MFFPMIter2 : public MFFPMIter
{
public:
  MFFPMIter2(const std::vector<INTERP_KERNEL::NormalizedCellType>& cts);
  void begin() { _it=_ids.begin(); }
  bool finished() const { return _it==_ids.end(); }
  void next() { _it++; }
  int current() const { return *_it; }
private:
  std::vector<int> _ids;
  std::vector<int>::const_iterator _it;
};

MFFPMIter *MFFPMIter::NewCell(const MEDFileEntities *entities)
{
  if(!entities)
    return new MFFPMIterSimple;
  else
    {
      const MEDFileStaticEntities *entities2(dynamic_cast<const MEDFileStaticEntities *>(entities));
      if(entities2)
        {
          std::vector<INTERP_KERNEL::NormalizedCellType> tmp;
          const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> >& myEnt(entities2->getEntries());
          for(std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> >::const_iterator it=myEnt.begin();it!=myEnt.end();it++)
            {
              if((*it).first==ON_CELLS || (*it).first==ON_GAUSS_NE || (*it).first==ON_GAUSS_PT)
                tmp.push_back((*it).second);
            }
          return new MFFPMIter2(tmp);
        }
      return new MFFPMIterSimple;// for MEDFileAllStaticEntites and MEDFileAllStaticEntitiesPlusDyn cells are in
    }
}

bool MFFPMIter::IsPresenceOfNode(const MEDFileEntities *entities)
{
  if(!entities)
    return true;
  else
    {
      const MEDFileStaticEntities *entities2(dynamic_cast<const MEDFileStaticEntities *>(entities));
      if(entities2)
        {
          const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> >& myEnt(entities2->getEntries());
          for(std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> >::const_iterator it=myEnt.begin();it!=myEnt.end();it++)
            if((*it).first==ON_NODES)
              return true;
          return false;
        }
      return true;// for MEDFileAllStaticEntites and MEDFileAllStaticEntitiesPlusDyn nodes are in
    }
}

MFFPMIter2::MFFPMIter2(const std::vector<INTERP_KERNEL::NormalizedCellType>& cts)
{
  std::size_t sz(cts.size());
  _ids.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      INTERP_KERNEL::NormalizedCellType *loc(std::find(typmai2,typmai2+MED_N_CELL_FIXED_GEO,cts[i]));
      if(loc!=typmai2+MED_N_CELL_FIXED_GEO)
        _ids[i]=(int)std::distance(typmai2,loc);
      else
        throw INTERP_KERNEL::Exception("MFFPMIter2 : The specified geo type does not exists !");
    }
}

/// @endcond

MEDFileFieldPerMesh::MEDFileFieldPerMesh(med_idt fid, MEDFileAnyTypeField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder, const MEDFileFieldNameScope& nasc, const MEDFileMesh *mm, const MEDFileEntities *entities):_mesh_iteration(meshIteration),_mesh_order(meshOrder),
    _father(fath)
{
  INTERP_KERNEL::AutoPtr<char> meshName(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> pflName(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> locName(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  const MEDFileUMesh *mmu(dynamic_cast<const MEDFileUMesh *>(mm));
  INTERP_KERNEL::AutoCppPtr<MFFPMIter> iter0(MFFPMIter::NewCell(entities));
  for(iter0->begin();!iter0->finished();iter0->next())
    {
      int nbProfile (MEDfield23nProfile(fid,nasc.getName().c_str(),getIteration(),getOrder(),MED_CELL        ,typmai[iter0->current()],meshCsit+1,meshName,pflName,locName));
      std::string name0(MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1));
      int nbProfile2(MEDfield23nProfile(fid,nasc.getName().c_str(),getIteration(),getOrder(),MED_NODE_ELEMENT,typmai[iter0->current()],meshCsit+1,meshName,pflName,locName));
      std::string name1(MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1));
      if(nbProfile>0 || nbProfile2>0)
        {
          const PartDefinition *pd(0);
          if(mmu)
            pd=mmu->getPartDefAtLevel(mmu->getRelativeLevOnGeoType(typmai2[iter0->current()]),typmai2[iter0->current()]);
          _field_pm_pt.push_back(MEDFileFieldPerMeshPerType::NewOnRead(fid,this,ON_CELLS,typmai2[iter0->current()],nasc,pd));
          if(nbProfile>0)
            setMeshName(name0);
          else
            setMeshName(name1);
        }
    }
  if(MFFPMIter::IsPresenceOfNode(entities))
    {
      int nbProfile(MEDfield23nProfile(fid,nasc.getName().c_str(),getIteration(),getOrder(),MED_NODE,MED_NONE,meshCsit+1,meshName,pflName,locName));
      if(nbProfile>0)
        {
          const PartDefinition *pd(0);
          if(mmu)
            pd=mmu->getPartDefAtLevel(1,INTERP_KERNEL::NORM_ERROR);
          _field_pm_pt.push_back(MEDFileFieldPerMeshPerType::NewOnRead(fid,this,ON_NODES,INTERP_KERNEL::NORM_ERROR,nasc,pd));
          setMeshName(MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE));
        }
    }
  if(!entities)
    return ;
  std::vector<int> dynGT(entities->getDynGTAvail());
  for(std::vector<int>::const_iterator it=dynGT.begin();it!=dynGT.end();it++)
    {
      int nbPfl(MEDfieldnProfile(fid,nasc.getName().c_str(),getIteration(),getOrder(),MED_STRUCT_ELEMENT,*it,pflName,locName));
      if(nbPfl>0)
        {
          _field_pm_pt.push_back(MEDFileFieldPerMeshPerTypeDyn::NewOnRead(fid,this,entities,*it,nasc));
          setMeshName(MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE));
        }
    }
  if(!_field_pm_pt.empty())
    return;
  //for vicious users using MED_ARETE MED_FACE in fields. the last try. For Others not overhead to pay.
  iter0=MFFPMIter::NewCell(entities);
  for(iter0->begin();!iter0->finished();iter0->next())
    {
      int nbProfile (MEDfield23nProfile(fid,nasc.getName().c_str(),getIteration(),getOrder(),MED_DESCENDING_FACE,typmai[iter0->current()],meshCsit+1,meshName,pflName,locName));
      std::string name0(MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1));
      int nbProfile2(MEDfield23nProfile(fid,nasc.getName().c_str(),getIteration(),getOrder(),MED_DESCENDING_EDGE,typmai[iter0->current()],meshCsit+1,meshName,pflName,locName));
      std::string name1(MEDLoaderBase::buildStringFromFortran(meshName,MED_NAME_SIZE+1));
      if(nbProfile>0 || nbProfile2>0)
        {
          _field_pm_pt.push_back(MEDFileFieldPerMeshPerType::NewOnRead(fid,this,ON_CELLS,typmai2[iter0->current()],nasc,NULL));
          if(nbProfile>0)
            setMeshName(name0);
          else
            setMeshName(name1);
        }
    }
}

MEDFileFieldPerMesh::MEDFileFieldPerMesh(MEDFileAnyTypeField1TSWithoutSDA *fath, const MEDCouplingMesh *mesh):_father(fath)
{
  copyTinyInfoFrom(mesh);
}

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
    bool operator() (const MCAuto<DataArrayInt>& pfl) { return _pfl==pfl->getName(); }
  private:
    const std::string& _pfl;
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

/*!
 * The returned value is never null.
 */
const DataArrayInt *MEDFileFieldGlobs::getProfile(const std::string& pflName) const
{
  std::string pflNameCpp(pflName);
  std::vector< MCAuto<DataArrayInt> >::const_iterator it=std::find_if(_pfls.begin(),_pfls.end(),MEDCouplingImpl::PflFinder(pflNameCpp));
  if(it==_pfls.end())
    {
      std::ostringstream oss; oss << "MEDFileFieldGlobs::getProfile: no such profile name : \"" << pflNameCpp << "\" Possible profiles are : ";
      for(it=_pfls.begin();it!=_pfls.end();it++)
        oss << "\"" << (*it)->getName() << "\", ";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  return *it;
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
 * The returned value is never null.
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
 *  \return DataArrayInt * - a non-const pointer to the profile array having the name \a pflName.
 *  \throw If there is no a profile named \a pflName.
 */
DataArrayInt *MEDFileFieldGlobsReal::getProfile(const std::string& pflName)
{
  return contentNotNull()->getProfile(pflName);
}

/*!
 * Returns a profile array, apt for modification, by its id.
 *  \param [in] pflId - the id of the profile of interest.
 *  \return DataArrayInt * - a non-const pointer to the profile array having the id \a pflId.
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

//= MEDFileAnyTypeField1TSWithoutSDA

void MEDFileAnyTypeField1TSWithoutSDA::deepCpyLeavesFrom(const MEDFileAnyTypeField1TSWithoutSDA& other)
{
  _field_per_mesh.resize(other._field_per_mesh.size());
  std::size_t i=0;
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it=other._field_per_mesh.begin();it!=other._field_per_mesh.end();it++,i++)
    {
      if((const MEDFileFieldPerMesh *)*it)
        _field_per_mesh[i]=(*it)->deepCopy(this);
    }
}

void MEDFileAnyTypeField1TSWithoutSDA::accept(MEDFileFieldVisitor& visitor) const
{
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    if((*it).isNotNull())
      {
        visitor.newMeshEntry(*it);
        (*it)->accept(visitor);
        visitor.endMeshEntry(*it);
      }
}

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
    oss << "[Type=" << getTypeStr() << "] with name \"" << getName() << "\" ";
  oss << "on one time Step ";
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
      for(std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it2=_field_per_mesh.begin();it2!=_field_per_mesh.end();it2++,i++)
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

std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > MEDFileAnyTypeField1TSWithoutSDA::splitComponents() const
{
  const DataArray *arr(getUndergroundDataArray());
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::splitComponents : no array defined !");
  int nbOfCompo=arr->getNumberOfComponents();
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > ret(nbOfCompo);
  for(int i=0;i<nbOfCompo;i++)
    {
      ret[i]=deepCopy();
      std::vector<int> v(1,i);
      MCAuto<DataArray> arr2=arr->keepSelectedComponents(v);
      ret[i]->setArray(arr2);
    }
  return ret;
}

MEDFileAnyTypeField1TSWithoutSDA::MEDFileAnyTypeField1TSWithoutSDA(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order):MEDFileFieldNameScope(fieldName,meshName),_iteration(iteration),_order(order),_csit(csit),_nb_of_tuples_to_be_allocated(-2)
{
}

MEDFileAnyTypeField1TSWithoutSDA::MEDFileAnyTypeField1TSWithoutSDA():_iteration(-1),_order(-1),_dt(0.),_csit(-1),_nb_of_tuples_to_be_allocated(-1)
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
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    (*it)->getDimension(ret);
  return ret;
}

bool MEDFileAnyTypeField1TSWithoutSDA::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab)
{
  bool ret=false;
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
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
int MEDFileAnyTypeField1TSWithoutSDA::getMeshIteration() const
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
int MEDFileAnyTypeField1TSWithoutSDA::getMeshOrder() const
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
void MEDFileAnyTypeField1TSWithoutSDA::fillTypesOfFieldAvailable(std::vector<TypeOfField>& types) const
{
  std::set<TypeOfField> types2;
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
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
std::vector<TypeOfField> MEDFileAnyTypeField1TSWithoutSDA::getTypesOfFieldAvailable() const
{
  std::vector<TypeOfField> ret;
  fillTypesOfFieldAvailable(ret);
  return ret;
}

std::vector<std::string> MEDFileAnyTypeField1TSWithoutSDA::getPflsReallyUsed2() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
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
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
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
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
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
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsedMulti();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

void MEDFileAnyTypeField1TSWithoutSDA::changePflsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    (*it)->changePflsRefsNamesGen(mapOfModif);
}

void MEDFileAnyTypeField1TSWithoutSDA::changeLocsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    (*it)->changeLocsRefsNamesGen(mapOfModif);
}

/*!
 * Returns all attributes of parts of \a this field lying on a given mesh.
 * Each part differs from other ones by a type of supporting mesh entity. The _i_-th
 * item of every of returned sequences refers to the _i_-th part of \a this field.
 * Thus all sequences returned by this method are of the same length equal to number
 * of different types of supporting entities.<br>
 * A field part can include sub-parts with several different spatial discretizations,
 * \ref MEDCoupling::ON_CELLS "ON_CELLS" and \ref MEDCoupling::ON_GAUSS_PT "ON_GAUSS_PT"
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
std::vector< std::vector< std::pair<int,int> > > MEDFileAnyTypeField1TSWithoutSDA::getFieldSplitedByType(const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getFieldSplitedByType : This is empty !");
  return _field_per_mesh[0]->getFieldSplitedByType(types,typesF,pfls,locs);
}

/*!
 * Returns dimensions of mesh elements \a this field lies on. The returned value is a
 * maximal absolute dimension and values returned via the out parameter \a levs are 
 * dimensions relative to the maximal absolute dimension. <br>
 * This method is designed for MEDFileField1TS instances that have a discretization
 * \ref MEDCoupling::ON_CELLS "ON_CELLS", 
 * \ref MEDCoupling::ON_GAUSS_PT "ON_GAUSS_PT", 
 * \ref MEDCoupling::ON_GAUSS_NE "ON_GAUSS_NE".
 * Only these 3 discretizations will be taken into account here. If \a this is
 * \ref MEDCoupling::ON_NODES "ON_NODES", -1 is returned and \a levs are empty.<br>
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
int MEDFileAnyTypeField1TSWithoutSDA::getNonEmptyLevels(const std::string& mname, std::vector<int>& levs) const
{
  levs.clear();
  std::vector<INTERP_KERNEL::NormalizedCellType> types;
  std::vector< std::vector<TypeOfField> > typesF;
  std::vector< std::vector<std::string> > pfls, locs;
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::getNonEmptyLevels : This is empty !");
  _field_per_mesh[0]->getFieldSplitedByType(types,typesF,pfls,locs);
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

void MEDFileAnyTypeField1TSWithoutSDA::convertMedBallIntoClassic()
{
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it<_field_per_mesh.end();it++)
    if((*it).isNotNull())
      (*it)->convertMedBallIntoClassic();
}

void MEDFileAnyTypeField1TSWithoutSDA::makeReduction(INTERP_KERNEL::NormalizedCellType ct, TypeOfField tof, const DataArrayInt *pfl)
{
  if(!pfl)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::makeReduction : null pfl !");
  std::string name(pfl->getName());
  pfl->checkAllocated();
  if(pfl->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::makeReduction : non mono compo array !");
  if(name.empty())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::makeReduction : empty pfl name !");
  if(_field_per_mesh.size()!=1)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::makeReduction : only single mesh supported !");
  MCAuto<MEDFileFieldPerMesh> fpm(_field_per_mesh[0]);
  if(fpm.isNull())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::makeReduction : only single not null mesh supported !");
  MEDFileFieldPerMeshPerTypePerDisc *disc(fpm->getLeafGivenTypeAndLocId(ct,0));
  if(disc->getType()!=tof)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::makeReduction : error !");
  int s(disc->getStart()),e(disc->getEnd()),nt(pfl->getNumberOfTuples());
  DataArray *arr(getUndergroundDataArray());
  int nt2(arr->getNumberOfTuples()),delta((e-s)-nt);
  if(delta<0)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::makeReduction : internal error !");
  MCAuto<DataArray> arr0(arr->selectByTupleIdSafeSlice(0,s,1)),arr1(arr->selectByTupleIdSafeSlice(s,e,1)),arr2(arr->selectByTupleIdSafeSlice(e,nt2,1));
  MCAuto<DataArray> arr11(arr1->selectByTupleIdSafe(pfl->begin(),pfl->end()));
  MCAuto<DataArray> arrOut(arr->buildNewEmptyInstance());
  arrOut->alloc(nt2-delta,arr->getNumberOfComponents());
  arrOut->copyStringInfoFrom(*arr);
  arrOut->setContigPartOfSelectedValuesSlice(0,arr0,0,s,1);
  arrOut->setContigPartOfSelectedValuesSlice(s,arr11,0,nt,1);
  arrOut->setContigPartOfSelectedValuesSlice(e-delta,arr2,0,nt2-e,1);
  setArray(arrOut);
  disc->setEnd(e-delta);
  disc->setProfile(name);
}

/*!
 * \param [in] mName specifies the underlying mesh name. This value can be pointer 0 for users that do not deal with fields on multi mesh.
 * \param [in] typ is for the geometric cell type (or INTERP_KERNEL::NORM_ERROR for node field) entry to find the right MEDFileFieldPerMeshPerTypePerDisc instance to set.
 * \param [in] locId is the localization id to find the right MEDFileFieldPerMeshPerTypePerDisc instance to set. It corresponds to the position of 
 *             \c pfls[std::distance(types.begin(),std::find(types.begin(),typ)] vector in MEDFileField1TSWithoutSDA::getFieldSplitedByType. For non gausspoints field users, the value is 0.
 */
MEDFileFieldPerMeshPerTypePerDisc *MEDFileAnyTypeField1TSWithoutSDA::getLeafGivenMeshAndTypeAndLocId(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId)
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::getLeafGivenMeshAndTypeAndLocId : This is empty !");
  return _field_per_mesh[0]->getLeafGivenTypeAndLocId(typ,locId);
}

/*!
 * \param [in] mName specifies the underlying mesh name. This value can be pointer 0 for users that do not deal with fields on multi mesh.
 * \param [in] typ is for the geometric cell type (or INTERP_KERNEL::NORM_ERROR for node field) entry to find the right MEDFileFieldPerMeshPerTypePerDisc instance to set.
 * \param [in] locId is the localization id to find the right MEDFileFieldPerMeshPerTypePerDisc instance to set. It corresponds to the position of 
 *             \c pfls[std::distance(types.begin(),std::find(types.begin(),typ)] vector in MEDFileField1TSWithoutSDA::getFieldSplitedByType. For non gausspoints field users, the value is 0.
 */
const MEDFileFieldPerMeshPerTypePerDisc *MEDFileAnyTypeField1TSWithoutSDA::getLeafGivenMeshAndTypeAndLocId(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId) const
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::getLeafGivenMeshAndTypeAndLocId : This is empty !");
  return _field_per_mesh[0]->getLeafGivenTypeAndLocId(typ,locId);
}

/*!
 * \param [in] mName specifies the underlying mesh name. This value can be pointer 0 for users that do not deal with fields on multi mesh.
 */
int MEDFileAnyTypeField1TSWithoutSDA::getMeshIdFromMeshName(const std::string& mName) const
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getMeshIdFromMeshName : No field set !");
  if(mName.empty())
    return 0;
  std::string mName2(mName);
  int ret=0;
  std::vector<std::string> msg;
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++,ret++)
    if(mName2==(*it)->getMeshName())
      return ret;
    else
      msg.push_back((*it)->getMeshName());
  std::ostringstream oss; oss << "MEDFileField1TSWithoutSDA::getMeshIdFromMeshName : No such mesh \"" << mName2 << "\" as underlying mesh of field \"" << getName() << "\" !\n";
  oss << "Possible meshes are : ";
  for(std::vector<std::string>::const_iterator it2=msg.begin();it2!=msg.end();it2++)
    oss << "\"" << (*it2) << "\" ";
  throw INTERP_KERNEL::Exception(oss.str());
}

int MEDFileAnyTypeField1TSWithoutSDA::addNewEntryIfNecessary(const MEDCouplingMesh *mesh)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::addNewEntryIfNecessary : input mesh is NULL !");
  std::string tmp(mesh->getName());
  if(tmp.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::addNewEntryIfNecessary : empty mesh name ! unsupported by MED file !");
  setMeshName(tmp);
  std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();
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

bool MEDFileAnyTypeField1TSWithoutSDA::renumberEntitiesLyingOnMesh(const std::string& meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N,
                                                                   MEDFileFieldGlobsReal& glob)
{
  bool ret=false;
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    {
      MEDFileFieldPerMesh *fpm(*it);
      if(fpm)
        ret=fpm->renumberEntitiesLyingOnMesh(meshName,oldCode,newCode,renumO2N,glob) || ret;
    }
  return ret;
}

/*!
 * This method splits \a this into several sub-parts so that each sub parts have exactly one spatial discretization. This method implements the minimal
 * splitting that leads to single spatial discretization of this.
 *
 * \sa splitMultiDiscrPerGeoTypes
 */
std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > MEDFileAnyTypeField1TSWithoutSDA::splitDiscretizations() const
{
  std::vector<INTERP_KERNEL::NormalizedCellType> types;
  std::vector< std::vector<TypeOfField> > typesF;
  std::vector< std::vector<std::string> > pfls,locs;
  std::vector< std::vector<std::pair<int,int> > > bgEnd(getFieldSplitedByType(getMeshName().c_str(),types,typesF,pfls,locs));
  std::set<TypeOfField> allEnt;
  for(std::vector< std::vector<TypeOfField> >::const_iterator it1=typesF.begin();it1!=typesF.end();it1++)
    for(std::vector<TypeOfField>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
      allEnt.insert(*it2);
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > ret(allEnt.size());
  std::set<TypeOfField>::const_iterator it3(allEnt.begin());
  for(std::size_t i=0;i<allEnt.size();i++,it3++)
    {
      std::vector< std::pair<int,int> > its;
      ret[i]=shallowCpy();
      int newLgth(ret[i]->keepOnlySpatialDiscretization(*it3,its));
      ret[i]->updateData(newLgth,its);
    }
  return ret;
}

/*!
 * This method performs a sub splitting as splitDiscretizations does but finer. This is the finest spliting level that can be done.
 * This method implements the minimal splitting so that each returned elements are mono Gauss discretization per geometric type.
 *
 * \sa splitDiscretizations
 */
std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > MEDFileAnyTypeField1TSWithoutSDA::splitMultiDiscrPerGeoTypes() const
{
  std::vector<INTERP_KERNEL::NormalizedCellType> types;
  std::vector< std::vector<TypeOfField> > typesF;
  std::vector< std::vector<std::string> > pfls,locs;
  std::vector< std::vector<std::pair<int,int> > > bgEnd(getFieldSplitedByType(getMeshName().c_str(),types,typesF,pfls,locs));
  std::set<TypeOfField> allEnt;
  std::size_t nbOfMDPGT(0),ii(0);
  for(std::vector< std::vector<TypeOfField> >::const_iterator it1=typesF.begin();it1!=typesF.end();it1++,ii++)
    {
      nbOfMDPGT=std::max(nbOfMDPGT,locs[ii].size());
      for(std::vector<TypeOfField>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
        allEnt.insert(*it2);
    }
  if(allEnt.size()!=1)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::splitMultiDiscrPerGeoTypes : this field is expected to be defined only on one spatial discretization !");
  if(nbOfMDPGT==0)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::splitMultiDiscrPerGeoTypes : empty field !");
  if(nbOfMDPGT==1)
    {
      std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > ret0(1);
      ret0[0]=const_cast<MEDFileAnyTypeField1TSWithoutSDA *>(this); this->incrRef();
      return ret0;
    }
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > ret(nbOfMDPGT);
  for(std::size_t i=0;i<nbOfMDPGT;i++)
    {
      std::vector< std::pair<int,int> > its;
      ret[i]=shallowCpy();
      int newLgth(ret[i]->keepOnlyGaussDiscretization(i,its));
      ret[i]->updateData(newLgth,its);
    }
  return ret;
}

int MEDFileAnyTypeField1TSWithoutSDA::keepOnlySpatialDiscretization(TypeOfField tof, std::vector< std::pair<int,int> >& its)
{
  int globalCounter(0);
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    (*it)->keepOnlySpatialDiscretization(tof,globalCounter,its);
  return globalCounter;
}

int MEDFileAnyTypeField1TSWithoutSDA::keepOnlyGaussDiscretization(std::size_t idOfDisc, std::vector< std::pair<int,int> >& its)
{
  int globalCounter(0);
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    (*it)->keepOnlyGaussDiscretization(idOfDisc,globalCounter,its);
  return globalCounter;
}

void MEDFileAnyTypeField1TSWithoutSDA::updateData(int newLgth, const std::vector< std::pair<int,int> >& oldStartStops)
{
  if(_nb_of_tuples_to_be_allocated>=0)
    {
      _nb_of_tuples_to_be_allocated=newLgth;
      const DataArray *oldArr(getUndergroundDataArray());
      if(oldArr)
        {
          MCAuto<DataArray> newArr(createNewEmptyDataArrayInstance());
          newArr->setInfoAndChangeNbOfCompo(oldArr->getInfoOnComponents());
          setArray(newArr);
          _nb_of_tuples_to_be_allocated=newLgth;//force the _nb_of_tuples_to_be_allocated because setArray has been used specialy
        }
      return ;
    }
  if(_nb_of_tuples_to_be_allocated==-1)
    return ;
  if(_nb_of_tuples_to_be_allocated==-2 || _nb_of_tuples_to_be_allocated==-3)
    {
      const DataArray *oldArr(getUndergroundDataArray());
      if(!oldArr || !oldArr->isAllocated())
        throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::updateData : internal error 1 !");
      MCAuto<DataArray> newArr(createNewEmptyDataArrayInstance());
      newArr->alloc(newLgth,getNumberOfComponents());
      if(oldArr)
        newArr->copyStringInfoFrom(*oldArr);
      int pos=0;
      for(std::vector< std::pair<int,int> >::const_iterator it=oldStartStops.begin();it!=oldStartStops.end();it++)
        {
          if((*it).second<(*it).first)
            throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::updateData : the range in the leaves was invalid !");
          newArr->setContigPartOfSelectedValuesSlice(pos,oldArr,(*it).first,(*it).second,1);
          pos+=(*it).second-(*it).first;
        }
      setArray(newArr);
      return ;
    }
  throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::updateData : internal error 2 !");
}

void MEDFileAnyTypeField1TSWithoutSDA::writeLL(med_idt fid, const MEDFileWritable& opts, const MEDFileFieldNameScope& nasc) const
{
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::writeLL : empty field !");
  if(_field_per_mesh.size()>1)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::writeLL : In MED3.0 mode in writting mode only ONE underlying mesh supported !");
  _field_per_mesh[0]->copyOptionsFrom(opts);
  _field_per_mesh[0]->writeLL(fid,nasc);
}

/*!
 * MED file does not support ' ' at the end of the field name. This method corrects the possibly invalid input \a nonCorrectFieldName to a correct one by right stripping input.
 */
std::string MEDFileAnyTypeField1TSWithoutSDA::FieldNameToMEDFileConvention(const std::string& nonCorrectFieldName)
{
  std::string::size_type pos0(nonCorrectFieldName.find_last_not_of(' '));
  if(pos0==std::string::npos)
    return nonCorrectFieldName;
  if(pos0+1==nonCorrectFieldName.length())
    return nonCorrectFieldName;
  return nonCorrectFieldName.substr(0,pos0+1);
}

/*!
 * This methods returns true is the allocation has been needed leading to a modification of state in \a this->_nb_of_tuples_to_be_allocated.
 * If false is returned the memory allocation is not required.
 */
bool MEDFileAnyTypeField1TSWithoutSDA::allocIfNecessaryTheArrayToReceiveDataFromFile()
{
  if(_nb_of_tuples_to_be_allocated>=0)
    {
      getOrCreateAndGetArray()->alloc(_nb_of_tuples_to_be_allocated,getNumberOfComponents());
      _nb_of_tuples_to_be_allocated=-2;
      return true;
    }
  if(_nb_of_tuples_to_be_allocated==-2 || _nb_of_tuples_to_be_allocated==-3)
    return false;
  if(_nb_of_tuples_to_be_allocated==-1)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::allocIfNecessaryTheArrayToReceiveDataFromFile : trying to read from a file an empty instance ! Need to prepare the structure before !");
  if(_nb_of_tuples_to_be_allocated<-3)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::allocIfNecessaryTheArrayToReceiveDataFromFile : internal error !");
  throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::allocIfNecessaryTheArrayToReceiveDataFromFile : internal error !");
}

void MEDFileAnyTypeField1TSWithoutSDA::loadOnlyStructureOfDataRecursively(med_idt fid, const MEDFileFieldNameScope& nasc, const MEDFileMeshes *ms, const MEDFileEntities *entities)
{
  med_int numdt,numit;
  med_float dt;
  med_int meshnumdt,meshnumit;
  MEDFILESAFECALLERRD0(MEDfieldComputingStepInfo,(fid,nasc.getName().c_str(),_csit,&numdt,&numit,&_dt));
  {
    med_bool localMesh;
    med_int nmesh;
    INTERP_KERNEL::AutoPtr<char> meshName(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
    MEDFILESAFECALLERRD0(MEDfield23ComputingStepMeshInfo,(fid,nasc.getName().c_str(),_csit,&numdt,&numit,&dt,&nmesh,meshName,&localMesh,&meshnumdt,&meshnumit)); // to check with Adrien for legacy MED files
  }
  //MEDFILESAFECALLERRD0(MEDfieldComputingStepMeshInfo,(fid,nasc.getName().c_str(),_csit,&numdt,&numit,&_dt,&meshnumdt,&meshnumit));
  if(_iteration!=numdt || _order!=numit)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::loadBigArraysRecursively : unexpected exception internal error !");
  _field_per_mesh.resize(1);
  //
  MEDFileMesh *mm(0);
  if(ms)
    {
      mm=ms->getMeshWithName(getMeshName());
    }
  //
  _field_per_mesh[0]=MEDFileFieldPerMesh::NewOnRead(fid,this,0,meshnumdt,meshnumit,nasc,mm,entities);
  _nb_of_tuples_to_be_allocated=0;
  _field_per_mesh[0]->loadOnlyStructureOfDataRecursively(fid,_nb_of_tuples_to_be_allocated,nasc);
}

void MEDFileAnyTypeField1TSWithoutSDA::loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc)
{
  allocIfNecessaryTheArrayToReceiveDataFromFile();
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    (*it)->loadBigArraysRecursively(fid,nasc);
}

void MEDFileAnyTypeField1TSWithoutSDA::loadBigArraysRecursivelyIfNecessary(med_idt fid, const MEDFileFieldNameScope& nasc)
{
  if(allocIfNecessaryTheArrayToReceiveDataFromFile())
    for(std::vector< MCAuto< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
      (*it)->loadBigArraysRecursively(fid,nasc);
}

void MEDFileAnyTypeField1TSWithoutSDA::loadStructureAndBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc, const MEDFileMeshes *ms, const MEDFileEntities *entities)
{
  loadOnlyStructureOfDataRecursively(fid,nasc,ms,entities);
  loadBigArraysRecursively(fid,nasc);
}

void MEDFileAnyTypeField1TSWithoutSDA::unloadArrays()
{
  DataArray *thisArr(getUndergroundDataArray());
  if(thisArr && thisArr->isAllocated())
    {
      _nb_of_tuples_to_be_allocated=thisArr->getNumberOfTuples();
      thisArr->desallocate();
    }
}

std::size_t MEDFileAnyTypeField1TSWithoutSDA::getHeapMemorySizeWithoutChildren() const
{
  return _mesh_name.capacity()+_dt_unit.capacity()+_field_per_mesh.capacity()*sizeof(MCAuto< MEDFileFieldPerMesh >);
}

std::vector<const BigMemoryObject *> MEDFileAnyTypeField1TSWithoutSDA::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  if(getUndergroundDataArray())
    ret.push_back(getUndergroundDataArray());
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    ret.push_back((const MEDFileFieldPerMesh *)*it);
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
void MEDFileAnyTypeField1TSWithoutSDA::setFieldNoProfileSBT(const TimeHolder *th, const MEDCouplingFieldTemplate *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc)
{
  const MEDCouplingMesh *mesh(field->getMesh());
  //
  TypeOfField type(field->getTypeOfField());
  std::vector<DataArrayInt *> dummy;
  if(mesh)
    setMeshName(mesh->getName());
  int start(copyTinyInfoFrom(th,field,arr));
  int pos(addNewEntryIfNecessary(mesh));
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
void MEDFileAnyTypeField1TSWithoutSDA::setFieldProfile(const TimeHolder *th, const MEDCouplingFieldTemplate *field, const DataArray *arrOfVals, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::setFieldProfile : input field is null !");
  if(!arrOfVals || !arrOfVals->isAllocated())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::setFieldProfile : input array is null or not allocated !");
  TypeOfField type=field->getTypeOfField();
  std::vector<DataArrayInt *> idsInPflPerType;
  std::vector<DataArrayInt *> idsPerType;
  std::vector<int> code,code2;
  MCAuto<MEDCouplingMesh> m(mesh->getMeshAtLevel(meshDimRelToMax));
  if(type!=ON_NODES)
    {
      m->splitProfilePerType(profile,code,idsInPflPerType,idsPerType);
      std::vector< MCAuto<DataArrayInt> > idsInPflPerType2(idsInPflPerType.size()); std::copy(idsInPflPerType.begin(),idsInPflPerType.end(),idsInPflPerType2.begin());
      std::vector< MCAuto<DataArrayInt> > idsPerType2(idsPerType.size()); std::copy(idsPerType.begin(),idsPerType.end(),idsPerType2.begin()); 
      std::vector<const DataArrayInt *> idsPerType3(idsPerType.size()); std::copy(idsPerType.begin(),idsPerType.end(),idsPerType3.begin());
      // start of check
      MCAuto<MEDCouplingFieldTemplate> field2=field->clone(false);
      int nbOfTuplesExp=field2->getNumberOfTuplesExpectedRegardingCode(code,idsPerType3);
      if(nbOfTuplesExp!=arrOfVals->getNumberOfTuples())
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeField1TSWithoutSDA::setFieldProfile : The array is expected to have " << nbOfTuplesExp << " tuples ! It has " << arrOfVals->getNumberOfTuples() << " !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      // end of check
      int start(copyTinyInfoFrom(th,field,arrOfVals));
      code2=m->getDistributionOfTypes();
      //
      int pos=addNewEntryIfNecessary(m);
      _field_per_mesh[pos]->assignFieldProfile(start,profile,code,code2,idsInPflPerType,idsPerType,field,arrOfVals,m,glob,nasc);
    }
  else
    {
      if(!profile || !profile->isAllocated() || profile->getNumberOfComponents()!=1)
        throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::setFieldProfile : input profile is null, not allocated or with number of components != 1 !");
      std::vector<int> v(3); v[0]=-1; v[1]=profile->getNumberOfTuples(); v[2]=0;
      std::vector<const DataArrayInt *> idsPerType3(1); idsPerType3[0]=profile;
      int nbOfTuplesExp=field->getNumberOfTuplesExpectedRegardingCode(v,idsPerType3);
      if(nbOfTuplesExp!=arrOfVals->getNumberOfTuples())
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeField1TSWithoutSDA::setFieldProfile : For node field, the array is expected to have " << nbOfTuplesExp << " tuples ! It has " << arrOfVals->getNumberOfTuples() << " !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      int start(copyTinyInfoFrom(th,field,arrOfVals));
      int pos(addNewEntryIfNecessary(m));
      _field_per_mesh[pos]->assignNodeFieldProfile(start,profile,field,arrOfVals,glob,nasc);
    }
}

/*!
 * \param [in] newNbOfTuples - The new nb of tuples to be allocated.
 */
void MEDFileAnyTypeField1TSWithoutSDA::allocNotFromFile(int newNbOfTuples)
{
  if(_nb_of_tuples_to_be_allocated>=0)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::allocNotFromFile : the object is expected to be appended to a data coming from a file but not loaded ! Load before appending data !");
  DataArray *arr(getOrCreateAndGetArray());
  arr->alloc(newNbOfTuples,arr->getNumberOfComponents());
  _nb_of_tuples_to_be_allocated=-3;
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
int MEDFileAnyTypeField1TSWithoutSDA::copyTinyInfoFrom(const TimeHolder *th, const MEDCouplingFieldTemplate *field, const DataArray *arr)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::copyTinyInfoFrom : input field is NULL !");
  std::string name(field->getName());
  setName(name.c_str());
  if(field->getMesh())
    setMeshName(field->getMesh()->getName());
  setDtUnit(th->getTimeUnit());
  if(name.empty())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::copyTinyInfoFrom : unsupported fields with no name in MED file !");
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::copyTinyInfoFrom : no array set !");
  if(!arr->isAllocated())
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::copyTinyInfoFrom : array is not allocated !");
  _dt=th->getTime(_iteration,_order);
  getOrCreateAndGetArray()->setInfoAndChangeNbOfCompo(arr->getInfoOnComponents());
  if(!getOrCreateAndGetArray()->isAllocated())
    {
      allocNotFromFile(arr->getNumberOfTuples());
      return 0;
    }
  else
    {
      int oldNbOfTuples=getOrCreateAndGetArray()->getNumberOfTuples();
      int newNbOfTuples=oldNbOfTuples+arr->getNumberOfTuples();
      getOrCreateAndGetArray()->reAlloc(newNbOfTuples);
      _nb_of_tuples_to_be_allocated=-3;
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
void MEDFileAnyTypeField1TSWithoutSDA::setInfo(const std::vector<std::string>& infos)
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

bool MEDFileAnyTypeField1TSWithoutSDA::presenceOfMultiDiscPerGeoType() const
{
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    {
      const MEDFileFieldPerMesh *fpm(*it);
      if(!fpm)
        continue;
      if(fpm->presenceOfMultiDiscPerGeoType())
        return true;
    }
  return false;
}

bool MEDFileAnyTypeField1TSWithoutSDA::presenceOfStructureElements() const
{
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    if((*it).isNotNull())
      if((*it)->presenceOfStructureElements())
        return true;
  return false;
}

bool MEDFileAnyTypeField1TSWithoutSDA::onlyStructureElements() const
{
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    if((*it).isNotNull())
      if(!(*it)->onlyStructureElements())
        return false;
  return true;
}

void MEDFileAnyTypeField1TSWithoutSDA::killStructureElements()
{
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    if((*it).isNotNull())
      (*it)->killStructureElements();
}

void MEDFileAnyTypeField1TSWithoutSDA::keepOnlyStructureElements()
{
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    if((*it).isNotNull())
      (*it)->keepOnlyStructureElements();
}

void MEDFileAnyTypeField1TSWithoutSDA::keepOnlyOnSE(const std::string& seName)
{
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    if((*it).isNotNull())
      (*it)->keepOnlyOnSE(seName);
}

void MEDFileAnyTypeField1TSWithoutSDA::getMeshSENames(std::vector< std::pair<std::string,std::string> >& ps) const
{
  for(std::vector< MCAuto< MEDFileFieldPerMesh > >::const_iterator it=_field_per_mesh.begin();it!=_field_per_mesh.end();it++)
    if((*it).isNotNull())
      (*it)->getMeshSENames(ps);
}

MEDCouplingFieldDouble *MEDFileAnyTypeField1TSWithoutSDA::fieldOnMesh(const MEDFileFieldGlobsReal *glob, const MEDFileMesh *mesh, MCAuto<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const
{
  static const char MSG0[]="MEDFileAnyTypeField1TSWithoutSDA::fieldOnMesh : the field is too complex to be able to be extracted with  \"field\" method ! Call getFieldOnMeshAtLevel method instead to deal with complexity !";
  if(_field_per_mesh.empty())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::fieldOnMesh : the field is empty ! Nothing to extract !");
  if(_field_per_mesh.size()>1)
    throw INTERP_KERNEL::Exception(MSG0);
  if(_field_per_mesh[0].isNull())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::fieldOnMesh : the field is inconsistent !");
  const MEDFileFieldPerMesh *pm(_field_per_mesh[0]);
  std::set<TypeOfField> types;
  pm->fillTypesOfFieldAvailable(types);
  if(types.size()!=1)
    throw INTERP_KERNEL::Exception(MSG0);
  TypeOfField type(*types.begin());
  int meshDimRelToMax(0);
  if(type==ON_NODES)
    meshDimRelToMax=0;
  else
    {
      int myDim(std::numeric_limits<int>::max());
      bool isUnique(pm->isUniqueLevel(myDim));
      if(!isUnique)
        throw INTERP_KERNEL::Exception(MSG0);
      meshDimRelToMax=myDim-mesh->getMeshDimension();
      if(meshDimRelToMax>0)
        throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TSWithoutSDA::fieldOnMesh : the mesh attached to field is not compatible with the field !");
    }
  return MEDFileAnyTypeField1TSWithoutSDA::getFieldOnMeshAtLevel(type,meshDimRelToMax,0/*renumPol*/,glob,mesh,arrOut,nasc);
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
MEDCouplingFieldDouble *MEDFileAnyTypeField1TSWithoutSDA::getFieldAtLevel(TypeOfField type, int meshDimRelToMax, const std::string& mName, int renumPol, const MEDFileFieldGlobsReal *glob, MCAuto<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const
{
  MCAuto<MEDFileMesh> mm;
  if(mName.empty())
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
MEDCouplingFieldDouble *MEDFileAnyTypeField1TSWithoutSDA::getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol, const MEDFileFieldGlobsReal *glob, const MEDFileMesh *mesh, MCAuto<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const
{
  MCAuto<MEDCouplingMesh> m(mesh->getMeshAtLevel(meshDimRelToMax,false));
  const DataArrayInt *d(mesh->getNumberFieldAtLevel(meshDimRelToMax)),*e(mesh->getNumberFieldAtLevel(1));
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
MEDCouplingFieldDouble *MEDFileAnyTypeField1TSWithoutSDA::getFieldAtTopLevel(TypeOfField type, const std::string& mName, int renumPol, const MEDFileFieldGlobsReal *glob, MCAuto<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const
{
  MCAuto<MEDFileMesh> mm;
  if(mName.empty())
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
MEDCouplingFieldDouble *MEDFileAnyTypeField1TSWithoutSDA::getFieldOnMeshAtLevel(TypeOfField type, int renumPol, const MEDFileFieldGlobsReal *glob, const MEDCouplingMesh *mesh, const DataArrayInt *cellRenum, const DataArrayInt *nodeRenum, MCAuto<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const
{
  static const char msg1[]="MEDFileField1TSWithoutSDA::getFieldOnMeshAtLevel : request for a renumbered field following mesh numbering whereas it is a profile field !";
  bool isPfl=false;
  MCAuto<MEDCouplingFieldDouble> ret=_field_per_mesh[0]->getFieldOnMeshAtLevel(type,glob,mesh,isPfl,arrOut,nasc);
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
                throw INTERP_KERNEL::Exception(oss.str());
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
                throw INTERP_KERNEL::Exception(oss.str());
              }
            MCAuto<DataArrayInt> nodeRenumSafe=nodeRenum->checkAndPreparePermutation();
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
DataArray *MEDFileAnyTypeField1TSWithoutSDA::getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl, const MEDFileFieldGlobsReal *glob, const MEDFileFieldNameScope& nasc) const
{
  MCAuto<MEDCouplingMesh> m(mesh->getMeshAtLevel(meshDimRelToMax));
  MCAuto<DataArray> ret=_field_per_mesh[0]->getFieldOnMeshAtLevelWithPfl(type,m,pfl,glob,nasc);
  ret->setName(nasc.getName().c_str());
  return ret.retn();
}

//= MEDFileField1TSWithoutSDA

/*!
 * Throws if a given value is not a valid (non-extended) relative dimension.
 *  \param [in] meshDimRelToMax - the relative dimension value.
 *  \throw If \a meshDimRelToMax > 0.
 */
void MEDFileField1TSWithoutSDA::CheckMeshDimRel(int meshDimRelToMax)
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
std::vector<int> MEDFileField1TSWithoutSDA::CheckSBTMesh(const MEDCouplingMesh *mesh)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::CheckSBTMesh : input mesh is NULL !");
  std::set<INTERP_KERNEL::NormalizedCellType> geoTypes=mesh->getAllGeoTypes();
  int nbOfTypes=geoTypes.size();
  std::vector<int> code(3*nbOfTypes);
  MCAuto<DataArrayInt> arr1=DataArrayInt::New();
  arr1->alloc(nbOfTypes,1);
  int *arrPtr=arr1->getPointer();
  std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator it=geoTypes.begin();
  for(int i=0;i<nbOfTypes;i++,it++)
    arrPtr[i]=std::distance(typmai2,std::find(typmai2,typmai2+MED_N_CELL_FIXED_GEO,*it));
  MCAuto<DataArrayInt> arr2=arr1->checkAndPreparePermutation();
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

MEDFileField1TSWithoutSDA *MEDFileField1TSWithoutSDA::New(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order, const std::vector<std::string>& infos)
{
  return new MEDFileField1TSWithoutSDA(fieldName,meshName,csit,iteration,order,infos);
}

/*!
 * Returns all attributes and values of parts of \a this field lying on a given mesh.
 * Each part differs from other ones by a type of supporting mesh entity. The _i_-th
 * item of every of returned sequences refers to the _i_-th part of \a this field.
 * Thus all sequences returned by this method are of the same length equal to number
 * of different types of supporting entities.<br>
 * A field part can include sub-parts with several different spatial discretizations,
 * \ref MEDCoupling::ON_CELLS "ON_CELLS" and \ref MEDCoupling::ON_GAUSS_PT "ON_GAUSS_PT"
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
 *          \ref MEDCoupling::ON_CELLS "ON_CELLS" and 
 *          \ref MEDCoupling::ON_GAUSS_PT "ON_GAUSS_PT" for example.
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
std::vector< std::vector<DataArrayDouble *> > MEDFileField1TSWithoutSDA::getFieldSplitedByType2(const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const
{
  if(mname.empty())
    if(_field_per_mesh.empty())
      throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getFieldSplitedByType : This is empty !");
  std::vector< std::vector< std::pair<int,int> > > ret0=_field_per_mesh[0]->getFieldSplitedByType(types,typesF,pfls,locs);
  int nbOfRet=ret0.size();
  std::vector< std::vector<DataArrayDouble *> > ret(nbOfRet);
  for(int i=0;i<nbOfRet;i++)
    {
      const std::vector< std::pair<int,int> >& p=ret0[i];
      int nbOfRet1=p.size();
      ret[i].resize(nbOfRet1);
      for(int j=0;j<nbOfRet1;j++)
        {
          DataArrayDouble *tmp=_arr->selectByTupleIdSafeSlice(p[j].first,p[j].second,1);
          ret[i][j]=tmp;
        }
    }
  return ret;
}

const char *MEDFileField1TSWithoutSDA::getTypeStr() const
{
  return TYPE_STR;
}

MEDFileIntField1TSWithoutSDA *MEDFileField1TSWithoutSDA::convertToInt() const
{
  MCAuto<MEDFileIntField1TSWithoutSDA> ret(new MEDFileIntField1TSWithoutSDA);
  ret->MEDFileAnyTypeField1TSWithoutSDA::operator =(*this);
  ret->deepCpyLeavesFrom(*this);
  const DataArrayDouble *arr(_arr);
  if(arr)
    {
      MCAuto<DataArrayInt> arr2(arr->convertToIntArr());
      ret->setArray(arr2);
    }
  return ret.retn();
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
DataArray *MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const
{
  return getUndergroundDataArrayTemplateExt(entries);
}

MEDFileField1TSWithoutSDA::MEDFileField1TSWithoutSDA(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order, const std::vector<std::string>& infos):MEDFileField1TSTemplateWithoutSDA<double>(fieldName,meshName,csit,iteration,order)
{
  DataArrayDouble *arr(getOrCreateAndGetArrayTemplate());
  arr->setInfoAndChangeNbOfCompo(infos);
}

MEDFileField1TSWithoutSDA::MEDFileField1TSWithoutSDA():MEDFileField1TSTemplateWithoutSDA<double>()
{
}

MEDFileField1TSWithoutSDA *MEDFileField1TSWithoutSDA::shallowCpy() const
{
  MCAuto<MEDFileField1TSWithoutSDA> ret(new MEDFileField1TSWithoutSDA(*this));
  ret->deepCpyLeavesFrom(*this);
  return ret.retn();
}

MEDFileField1TSWithoutSDA *MEDFileField1TSWithoutSDA::deepCopy() const
{
  MCAuto<MEDFileField1TSWithoutSDA> ret(shallowCpy());
  if(_arr.isNotNull())
    ret->_arr=_arr->deepCopy();
  return ret.retn();
}

//= MEDFileIntField1TSWithoutSDA

MEDFileIntField1TSWithoutSDA *MEDFileIntField1TSWithoutSDA::New(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order, const std::vector<std::string>& infos)
{
  return new MEDFileIntField1TSWithoutSDA(fieldName,meshName,csit,iteration,order,infos);
}

MEDFileIntField1TSWithoutSDA::MEDFileIntField1TSWithoutSDA()
{
}

MEDFileIntField1TSWithoutSDA::MEDFileIntField1TSWithoutSDA(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order,
                                                           const std::vector<std::string>& infos):MEDFileField1TSNDTemplateWithoutSDA<int>(fieldName,meshName,csit,iteration,order,infos)
{
  DataArrayInt *arr(getOrCreateAndGetArrayTemplate());
  arr->setInfoAndChangeNbOfCompo(infos);
}

const char *MEDFileIntField1TSWithoutSDA::getTypeStr() const
{
  return TYPE_STR;
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
DataArray *MEDFileIntField1TSWithoutSDA::getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const
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
DataArrayInt *MEDFileIntField1TSWithoutSDA::getUndergroundDataArrayIntExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const
{
  if(_field_per_mesh.size()!=1)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt : field lies on several meshes, this method has no sense !");
  if(_field_per_mesh[0]==0)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt : no field specified !");
  _field_per_mesh[0]->getUndergroundDataArrayExt(entries);
  return getUndergroundDataArrayTemplate();
}

MEDFileIntField1TSWithoutSDA *MEDFileIntField1TSWithoutSDA::shallowCpy() const
{
  MCAuto<MEDFileIntField1TSWithoutSDA> ret(new MEDFileIntField1TSWithoutSDA(*this));
  ret->deepCpyLeavesFrom(*this);
  return ret.retn();
}

MEDFileIntField1TSWithoutSDA *MEDFileIntField1TSWithoutSDA::deepCopy() const
{
  MCAuto<MEDFileIntField1TSWithoutSDA> ret(shallowCpy());
  if(_arr.isNotNull())
    ret->_arr=_arr->deepCopy();
  return ret.retn();
}

//= MEDFileFloatField1TSWithoutSDA

MEDFileFloatField1TSWithoutSDA *MEDFileFloatField1TSWithoutSDA::New(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order, const std::vector<std::string>& infos)
{
  return new MEDFileFloatField1TSWithoutSDA(fieldName,meshName,csit,iteration,order,infos);
}

MEDFileFloatField1TSWithoutSDA::MEDFileFloatField1TSWithoutSDA()
{
}

MEDFileFloatField1TSWithoutSDA::MEDFileFloatField1TSWithoutSDA(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order,
                                                               const std::vector<std::string>& infos):MEDFileField1TSNDTemplateWithoutSDA<float>(fieldName,meshName,csit,iteration,order,infos)
{
  DataArrayFloat *arr(getOrCreateAndGetArrayTemplate());
  arr->setInfoAndChangeNbOfCompo(infos);
}

const char *MEDFileFloatField1TSWithoutSDA::getTypeStr() const
{
  return TYPE_STR;
}

/*!
 * Returns a pointer to the underground DataArrayFloat instance and a
 * sequence describing parameters of a support of each part of \a this field. The
 * caller should not decrRef() the returned DataArrayFloat. This method allows for a
 * direct access to the field values. This method is intended for the field lying on one
 * mesh only.
 *  \param [in,out] entries - the sequence describing parameters of a support of each
 *         part of \a this field. Each item of this sequence consists of two parts. The
 *         first part describes a type of mesh entity and an id of discretization of a
 *         current field part. The second part describes a range of values [begin,end)
 *         within the returned array relating to the current field part.
 *  \return DataArrayFloat * - the pointer to the field values array.
 *  \throw If the number of underlying meshes is not equal to 1.
 *  \throw If no field values are available.
 *  \sa getUndergroundDataArray()
 */
DataArray *MEDFileFloatField1TSWithoutSDA::getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const
{
  return getUndergroundDataArrayFloatExt(entries);
}

/*!
 * Returns a pointer to the underground DataArrayFloat instance and a
 * sequence describing parameters of a support of each part of \a this field. The
 * caller should not decrRef() the returned DataArrayFloat. This method allows for a
 * direct access to the field values. This method is intended for the field lying on one
 * mesh only.
 *  \param [in,out] entries - the sequence describing parameters of a support of each
 *         part of \a this field. Each item of this sequence consists of two parts. The
 *         first part describes a type of mesh entity and an id of discretization of a
 *         current field part. The second part describes a range of values [begin,end)
 *         within the returned array relating to the current field part.
 *  \return DataArrayFloat * - the pointer to the field values array.
 *  \throw If the number of underlying meshes is not equal to 1.
 *  \throw If no field values are available.
 *  \sa getUndergroundDataArray()
 */
DataArrayFloat *MEDFileFloatField1TSWithoutSDA::getUndergroundDataArrayFloatExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const
{
  if(_field_per_mesh.size()!=1)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt : field lies on several meshes, this method has no sense !");
  if(_field_per_mesh[0]==0)
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt : no field specified !");
  _field_per_mesh[0]->getUndergroundDataArrayExt(entries);
  return getUndergroundDataArrayTemplate();
}

MEDFileFloatField1TSWithoutSDA *MEDFileFloatField1TSWithoutSDA::shallowCpy() const
{
  MCAuto<MEDFileFloatField1TSWithoutSDA> ret(new MEDFileFloatField1TSWithoutSDA(*this));
  ret->deepCpyLeavesFrom(*this);
  return ret.retn();
}

MEDFileFloatField1TSWithoutSDA *MEDFileFloatField1TSWithoutSDA::deepCopy() const
{
  MCAuto<MEDFileFloatField1TSWithoutSDA> ret(shallowCpy());
  if(_arr.isNotNull())
    ret->_arr=_arr->deepCopy();
  return ret.retn();
}

//= MEDFileAnyTypeField1TS

MEDFileAnyTypeField1TS::MEDFileAnyTypeField1TS()
{
}

MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeField1TS::BuildContentFrom(med_idt fid, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
{
  med_field_type typcha;
  //
  std::vector<std::string> infos;
  std::string dtunit,fieldName,meshName;
  LocateField2(fid,0,true,fieldName,typcha,infos,dtunit,meshName);
  MCAuto<MEDFileAnyTypeField1TSWithoutSDA> ret;
  switch(typcha)
  {
    case MED_FLOAT64:
      {
        ret=MEDFileField1TSWithoutSDA::New(fieldName,meshName,-1,-1/*iteration*/,-1/*order*/,std::vector<std::string>());
        break;
      }
    case MED_INT32:
      {
        ret=MEDFileIntField1TSWithoutSDA::New(fieldName,meshName,-1,-1/*iteration*/,-1/*order*/,std::vector<std::string>());
        break;
      }
    case MED_FLOAT32:
      {
        ret=MEDFileFloatField1TSWithoutSDA::New(fieldName,meshName,-1,-1/*iteration*/,-1/*order*/,std::vector<std::string>());
        break;
      }
    case MED_INT:
      {
        if(sizeof(med_int)==sizeof(int))
          {
            ret=MEDFileIntField1TSWithoutSDA::New(fieldName,meshName,-1,-1/*iteration*/,-1/*order*/,std::vector<std::string>());
            break;
          }
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileAnyTypeField1TS::BuildContentFrom(fid) : file \'" << FileNameFromFID(fid) << "\' contains field with name \'" << fieldName << "\' but the type of the first field is not in [MED_FLOAT64, MED_INT32] !";
        throw INTERP_KERNEL::Exception(oss.str());
      }
  }
  ret->setDtUnit(dtunit.c_str());
  ret->getOrCreateAndGetArray()->setInfoAndChangeNbOfCompo(infos);
  //
  med_int numdt,numit;
  med_float dt;
  MEDFILESAFECALLERRD0(MEDfieldComputingStepInfo,(fid,fieldName.c_str(),1,&numdt,&numit,&dt));
  ret->setTime(numdt,numit,dt);
  ret->_csit=1;
  if(loadAll)
    ret->loadStructureAndBigArraysRecursively(fid,*((const MEDFileAnyTypeField1TSWithoutSDA*)ret),ms,entities);
  else
    ret->loadOnlyStructureOfDataRecursively(fid,*((const MEDFileAnyTypeField1TSWithoutSDA*)ret),ms,entities);
  return ret.retn();
}

MEDFileAnyTypeField1TS::MEDFileAnyTypeField1TS(med_idt fid, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
try:MEDFileFieldGlobsReal(fid)
{
  _content=BuildContentFrom(fid,loadAll,ms,entities);
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
{
    throw e;
}

MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeField1TS::BuildContentFrom(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
{
  med_field_type typcha;
  std::vector<std::string> infos;
  std::string dtunit,meshName;
  int nbSteps(0);
  {
    int iii=-1;
    nbSteps=LocateField(fid,fieldName,iii,typcha,infos,dtunit,meshName);
  }
  MCAuto<MEDFileAnyTypeField1TSWithoutSDA> ret;
  switch(typcha)
  {
    case MED_FLOAT64:
      {
        ret=MEDFileField1TSWithoutSDA::New(fieldName,meshName,-1,-1/*iteration*/,-1/*order*/,std::vector<std::string>());
        break;
      }
    case MED_INT32:
      {
        ret=MEDFileIntField1TSWithoutSDA::New(fieldName,meshName,-1,-1/*iteration*/,-1/*order*/,std::vector<std::string>());
        break;
      }
    case MED_FLOAT32:
      {
        ret=MEDFileFloatField1TSWithoutSDA::New(fieldName,meshName,-1,-1/*iteration*/,-1/*order*/,std::vector<std::string>());
        break;
      }
    case MED_INT:
      {
        if(sizeof(med_int)==sizeof(int))
          {
            ret=MEDFileIntField1TSWithoutSDA::New(fieldName,meshName,-1,-1/*iteration*/,-1/*order*/,std::vector<std::string>());
            break;
          }
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileAnyTypeField1TS::BuildContentFrom(fid,fieldName) : file \'" << FileNameFromFID(fid) << "\' contains field with name \'" << fieldName << "\' but the type of field is not in [MED_FLOAT64, MED_INT32, MED_FLOAT32] !";
        throw INTERP_KERNEL::Exception(oss.str());
      }
  }
  ret->setMeshName(meshName);
  ret->setDtUnit(dtunit.c_str());
  ret->getOrCreateAndGetArray()->setInfoAndChangeNbOfCompo(infos);
  //
  if(nbSteps<1)
    {
      std::ostringstream oss; oss << "MEDFileField1TS(fid,fieldName) : file \'" << FileNameFromFID(fid) << "\' contains field with name \'" << fieldName << "\' but there is no time steps on it !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  //
  med_int numdt,numit;
  med_float dt;
  MEDFILESAFECALLERRD0(MEDfieldComputingStepInfo,(fid,fieldName.c_str(),1,&numdt,&numit,&dt));
  ret->setTime(numdt,numit,dt);
  ret->_csit=1;
  if(loadAll)
    ret->loadStructureAndBigArraysRecursively(fid,*((const MEDFileAnyTypeField1TSWithoutSDA*)ret),ms,entities);
  else
    ret->loadOnlyStructureOfDataRecursively(fid,*((const MEDFileAnyTypeField1TSWithoutSDA*)ret),ms,entities);
  return ret.retn();
}

MEDFileAnyTypeField1TS::MEDFileAnyTypeField1TS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
try:MEDFileFieldGlobsReal(fid)
{
  _content=BuildContentFrom(fid,fieldName,loadAll,ms,entities);
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
{
    throw e;
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::BuildNewInstanceFromContent(MEDFileAnyTypeField1TSWithoutSDA *c)
{
  if(!c)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TS::BuildNewInstanceFromContent : empty content in input : unable to build a new instance !");
  if(dynamic_cast<const MEDFileField1TSWithoutSDA *>(c))
    {
      MCAuto<MEDFileField1TS> ret(MEDFileField1TS::New());
      ret->_content=c; c->incrRef();
      return ret.retn();
    }
  if(dynamic_cast<const MEDFileIntField1TSWithoutSDA *>(c))
    {
      MCAuto<MEDFileIntField1TS> ret(MEDFileIntField1TS::New());
      ret->_content=c; c->incrRef();
      return ret.retn();
    }
  if(dynamic_cast<const MEDFileFloatField1TSWithoutSDA *>(c))
    {
      MCAuto<MEDFileFloatField1TS> ret(MEDFileFloatField1TS::New());
      ret->_content=c; c->incrRef();
      return ret.retn();
    }
  throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TS::BuildNewInstanceFromContent : internal error ! a content of type different from FLOAT64 FLOAT32 and INT32 has been built but not intercepted !");
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::BuildNewInstanceFromContent(MEDFileAnyTypeField1TSWithoutSDA *c, med_idt fid)
{
  MEDFileAnyTypeField1TS *ret(BuildNewInstanceFromContent(c));
  ret->setFileName(FileNameFromFID(fid));
  return ret;
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::New(const std::string& fileName, bool loadAll)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return New(fid,loadAll);
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::New(med_idt fid, bool loadAll)
{
  MCAuto<MEDFileAnyTypeField1TSWithoutSDA> c(BuildContentFrom(fid,loadAll,0,0));
  MCAuto<MEDFileAnyTypeField1TS> ret(BuildNewInstanceFromContent(c,fid));
  ret->loadGlobals(fid);
  return ret.retn();
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::New(const std::string& fileName, const std::string& fieldName, bool loadAll)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return New(fid,fieldName,loadAll);
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::New(med_idt fid, const std::string& fieldName, bool loadAll)
{
  MCAuto<MEDFileAnyTypeField1TSWithoutSDA> c(BuildContentFrom(fid,fieldName,loadAll,0,0));
  MCAuto<MEDFileAnyTypeField1TS> ret(BuildNewInstanceFromContent(c,fid));
  ret->loadGlobals(fid);
  return ret.retn();
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::New(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return New(fid,fieldName,iteration,order,loadAll);
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::New(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll)
{
  MCAuto<MEDFileAnyTypeField1TSWithoutSDA> c(BuildContentFrom(fid,fieldName,iteration,order,loadAll,0,0));
  MCAuto<MEDFileAnyTypeField1TS> ret(BuildNewInstanceFromContent(c,fid));
  ret->loadGlobals(fid);
  return ret.retn();
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::NewAdv(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileEntities *entities)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return NewAdv(fid,fieldName,iteration,order,loadAll,entities);
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::NewAdv(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileEntities *entities)
{
  MCAuto<MEDFileAnyTypeField1TSWithoutSDA> c(BuildContentFrom(fid,fieldName,iteration,order,loadAll,0,entities));
  MCAuto<MEDFileAnyTypeField1TS> ret(BuildNewInstanceFromContent(c,fid));
  ret->loadGlobals(fid);
  return ret.retn();
}

MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeField1TS::BuildContentFrom(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
{
  med_field_type typcha;
  std::vector<std::string> infos;
  std::string dtunit,meshName;
  int iii(-1);
  int nbOfStep2(LocateField(fid,fieldName,iii,typcha,infos,dtunit,meshName));
  MCAuto<MEDFileAnyTypeField1TSWithoutSDA> ret;
  switch(typcha)
  {
    case MED_FLOAT64:
      {
        ret=MEDFileField1TSWithoutSDA::New(fieldName,meshName,-1,iteration,order,std::vector<std::string>());
        break;
      }
    case MED_INT32:
      {
        ret=MEDFileIntField1TSWithoutSDA::New(fieldName,meshName,-1,iteration,order,std::vector<std::string>());
        break;
      }
    case MED_FLOAT32:
      {
        ret=MEDFileFloatField1TSWithoutSDA::New(fieldName,meshName,-1,iteration,order,std::vector<std::string>());
        break;
      }
    case MED_INT:
      {
        if(sizeof(med_int)==sizeof(int))
          {
            ret=MEDFileIntField1TSWithoutSDA::New(fieldName,meshName,-1,iteration,order,std::vector<std::string>());
            break;
          }
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileAnyTypeField1TS::BuildContentFrom(fid,fieldName,iteration,order) : file \'" << FileNameFromFID(fid) << "\' contains field with name \'" << fieldName << "\' but the type of field is not in [MED_FLOAT64, MED_INT32, MED_FLOAT32] !";
        throw INTERP_KERNEL::Exception(oss.str());
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
      MEDFILESAFECALLERRD0(MEDfieldComputingStepInfo,(fid,fieldName.c_str(),i+1,&numdt,&numit,&dt));
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
      std::ostringstream oss; oss << "No such iteration (" << iteration << "," << order << ") in existing field '" << fieldName << "' in file '" << FileNameFromFID(fid) << "' ! Available iterations are : ";
      for(std::vector< std::pair<int,int> >::const_iterator iter=dtits.begin();iter!=dtits.end();iter++)
        oss << "(" << (*iter).first << "," << (*iter).second << "), ";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  if(loadAll)
    ret->loadStructureAndBigArraysRecursively(fid,*((const MEDFileAnyTypeField1TSWithoutSDA*)ret),ms,entities);
  else
    ret->loadOnlyStructureOfDataRecursively(fid,*((const MEDFileAnyTypeField1TSWithoutSDA*)ret),ms,entities);
  return ret.retn();
}

MEDFileAnyTypeField1TS::MEDFileAnyTypeField1TS(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
try:MEDFileFieldGlobsReal(fid)
{
  _content=BuildContentFrom(fid,fieldName,iteration,order,loadAll,ms,entities);
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

int MEDFileAnyTypeField1TS::LocateField2(med_idt fid, int fieldIdCFormat, bool checkFieldId, std::string& fieldName, med_field_type& typcha, std::vector<std::string>& infos, std::string& dtunitOut, std::string& meshName)
{
  if(checkFieldId)
    {
      int nbFields=MEDnField(fid);
      if(fieldIdCFormat>=nbFields)
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeField1TS::LocateField2(fileName) : in file \'" << FileNameFromFID(fid) << "\' number of fields is " << nbFields << " ! Trying to request for id " << fieldIdCFormat << " !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  int ncomp(MEDfieldnComponent(fid,fieldIdCFormat+1));
  INTERP_KERNEL::AutoPtr<char> comp(MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> unit(MEDLoaderBase::buildEmptyString(ncomp*MED_SNAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> dtunit(MEDLoaderBase::buildEmptyString(MED_LNAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> nomcha(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> nomMaa(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  med_bool localMesh;
  int nbOfStep;
  MEDFILESAFECALLERRD0(MEDfieldInfo,(fid,fieldIdCFormat+1,nomcha,nomMaa,&localMesh,&typcha,comp,unit,dtunit,&nbOfStep));
  fieldName=MEDLoaderBase::buildStringFromFortran(nomcha,MED_NAME_SIZE);
  dtunitOut=MEDLoaderBase::buildStringFromFortran(dtunit,MED_LNAME_SIZE);
  meshName=MEDLoaderBase::buildStringFromFortran(nomMaa,MED_NAME_SIZE);
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
int MEDFileAnyTypeField1TS::LocateField(med_idt fid, const std::string& fieldName, int& posCFormat, med_field_type& typcha, std::vector<std::string>& infos, std::string& dtunitOut, std::string& meshName)
{
  int nbFields=MEDnField(fid);
  bool found=false;
  std::vector<std::string> fns(nbFields);
  int nbOfStep2(-1);
  for(int i=0;i<nbFields && !found;i++)
    {
      std::string tmp,tmp2;
      nbOfStep2=LocateField2(fid,i,false,tmp,typcha,infos,dtunitOut,tmp2);
      fns[i]=tmp;
      found=(tmp==fieldName);
      if(found)
        {
          posCFormat=i;
          meshName=tmp2;
        }
    }
  if(!found)
    {
      std::ostringstream oss; oss << "No such field '" << fieldName << "' in file '" << FileNameFromFID(fid) << "' ! Available fields are : ";
      for(std::vector<std::string>::const_iterator it=fns.begin();it!=fns.end();it++)
        oss << "\"" << *it << "\" ";
      throw INTERP_KERNEL::Exception(oss.str());
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
void MEDFileAnyTypeField1TS::setProfileNameOnLeaf(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const std::string& newPflName, bool forceRenameOnGlob)
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
      throw INTERP_KERNEL::Exception(oss.str());
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
void MEDFileAnyTypeField1TS::setLocNameOnLeaf(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const std::string& newLocName, bool forceRenameOnGlob)
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
      throw INTERP_KERNEL::Exception(oss.str());
    }
}

MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeField1TS::contentNotNullBase()
{
  MEDFileAnyTypeField1TSWithoutSDA *ret=_content;
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TS : content is expected to be not null !");
  return ret;
}

const MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeField1TS::contentNotNullBase() const
{
  const MEDFileAnyTypeField1TSWithoutSDA *ret=_content;
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TS : const content is expected to be not null !");
  return ret;
}

/*!
 * This method alloc the arrays and load potentially huge arrays contained in this field.
 * This method should be called when a MEDFileAnyTypeField1TS::New constructor has been with false as the last parameter.
 * This method can be also called to refresh or reinit values from a file.
 * 
 * \throw If the fileName is not set or points to a non readable MED file.
 * \sa MEDFileAnyTypeField1TS::loadArraysIfNecessary
 */
void MEDFileAnyTypeField1TS::loadArrays()
{
  if(getFileName().empty())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TS::loadArrays : the structure does not come from a file !");
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(getFileName()));
  contentNotNullBase()->loadBigArraysRecursively(fid,*contentNotNullBase());
}

/*!
 * This method behaves as MEDFileAnyTypeField1TS::loadArrays does, the first call, if \a this was built using a file without loading big arrays.
 * But once data loaded once, this method does nothing. Contrary to MEDFileAnyTypeField1TS::loadArrays and MEDFileAnyTypeField1TS::unloadArrays
 * this method does not throw if \a this does not come from file read.
 * 
 * \sa MEDFileAnyTypeField1TS::loadArrays, MEDFileAnyTypeField1TS::unloadArrays
 */
void MEDFileAnyTypeField1TS::loadArraysIfNecessary()
{
  if(!getFileName().empty())
    {
      MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(getFileName()));
      contentNotNullBase()->loadBigArraysRecursivelyIfNecessary(fid,*contentNotNullBase());
    }
}

/*!
 * This method releases potentially big data arrays and so returns to the same heap memory than status loaded with 'loadAll' parameter set to false.
 * \b WARNING, this method does release arrays even if \a this does not come from a load of a MED file.
 * So this method can lead to a loss of data. If you want to unload arrays safely call MEDFileAnyTypeField1TS::unloadArraysWithoutDataLoss instead.
 * 
 * \sa MEDFileAnyTypeField1TS::loadArrays, MEDFileAnyTypeField1TS::loadArraysIfNecessary, MEDFileAnyTypeField1TS::unloadArraysWithoutDataLoss
 */
void MEDFileAnyTypeField1TS::unloadArrays()
{
  contentNotNullBase()->unloadArrays();
}

/*!
 * This method potentially releases big data arrays if \a this is coming from a file. If \a this has been built from scratch this method will have no effect.
 * This method is the symetrical method of MEDFileAnyTypeField1TS::loadArraysIfNecessary.
 * This method is useful to reduce \b safely amount of heap memory necessary for \a this by using MED file as database.
 * 
 * \sa MEDFileAnyTypeField1TS::loadArraysIfNecessary
 */
void MEDFileAnyTypeField1TS::unloadArraysWithoutDataLoss()
{
  if(!getFileName().empty())
    contentNotNullBase()->unloadArrays();
}

void MEDFileAnyTypeField1TS::writeLL(med_idt fid) const
{
  int nbComp(getNumberOfComponents());
  INTERP_KERNEL::AutoPtr<char> comp(MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> unit(MEDLoaderBase::buildEmptyString(nbComp*MED_SNAME_SIZE));
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
  MEDFILESAFECALLERWR0(MEDfieldCr,(fid,getName().c_str(),getMEDFileFieldType(),nbComp,comp,unit,getDtUnit().c_str(),getMeshName().c_str()));
  writeGlobals(fid,*this);
  contentNotNullBase()->writeLL(fid,*this,*contentNotNullBase());
}

std::size_t MEDFileAnyTypeField1TS::getHeapMemorySizeWithoutChildren() const
{
  return MEDFileFieldGlobsReal::getHeapMemorySizeWithoutChildren();
}

std::vector<const BigMemoryObject *> MEDFileAnyTypeField1TS::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(MEDFileFieldGlobsReal::getDirectChildrenWithNull());
  ret.push_back((const MEDFileAnyTypeField1TSWithoutSDA *)_content);
  return ret;
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

void MEDFileAnyTypeField1TS::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  contentNotNullBase()->changePflsRefsNamesGen2(mapOfModif);
}

void MEDFileAnyTypeField1TS::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
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

void MEDFileAnyTypeField1TS::setName(const std::string& name)
{
  contentNotNullBase()->setName(name);
}

void MEDFileAnyTypeField1TS::simpleRepr(int bkOffset, std::ostream& oss, int f1tsId) const
{
  contentNotNullBase()->simpleRepr(bkOffset,oss,f1tsId);
}

std::string MEDFileAnyTypeField1TS::getDtUnit() const
{
  return contentNotNullBase()->getDtUnit();
}

void MEDFileAnyTypeField1TS::setDtUnit(const std::string& dtUnit)
{
  contentNotNullBase()->setDtUnit(dtUnit);
}

std::string MEDFileAnyTypeField1TS::getMeshName() const
{
  return contentNotNullBase()->getMeshName();
}

void MEDFileAnyTypeField1TS::setMeshName(const std::string& newMeshName)
{
  contentNotNullBase()->setMeshName(newMeshName);
}

bool MEDFileAnyTypeField1TS::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab)
{
  return contentNotNullBase()->changeMeshNames(modifTab);
}

int MEDFileAnyTypeField1TS::getMeshIteration() const
{
  return contentNotNullBase()->getMeshIteration();
}

int MEDFileAnyTypeField1TS::getMeshOrder() const
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

void MEDFileAnyTypeField1TS::fillTypesOfFieldAvailable(std::vector<TypeOfField>& types) const
{
  contentNotNullBase()->fillTypesOfFieldAvailable(types);
}

void MEDFileAnyTypeField1TS::setInfo(const std::vector<std::string>& infos)
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

bool MEDFileAnyTypeField1TS::presenceOfMultiDiscPerGeoType() const
{
  return contentNotNullBase()->presenceOfMultiDiscPerGeoType();
}

MEDFileFieldPerMeshPerTypePerDisc *MEDFileAnyTypeField1TS::getLeafGivenMeshAndTypeAndLocId(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId)
{
  return contentNotNullBase()->getLeafGivenMeshAndTypeAndLocId(mName,typ,locId);
}

const MEDFileFieldPerMeshPerTypePerDisc *MEDFileAnyTypeField1TS::getLeafGivenMeshAndTypeAndLocId(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId) const
{
  return contentNotNullBase()->getLeafGivenMeshAndTypeAndLocId(mName,typ,locId);
}

int MEDFileAnyTypeField1TS::getNonEmptyLevels(const std::string& mname, std::vector<int>& levs) const
{
  return contentNotNullBase()->getNonEmptyLevels(mname,levs);
}

void MEDFileAnyTypeField1TS::convertMedBallIntoClassic()
{
  return contentNotNullBase()->convertMedBallIntoClassic();
}

void MEDFileAnyTypeField1TS::makeReduction(INTERP_KERNEL::NormalizedCellType ct, TypeOfField tof, const DataArrayInt *pfl)
{
  return contentNotNullBase()->makeReduction(ct,tof,pfl);
}

std::vector<TypeOfField> MEDFileAnyTypeField1TS::getTypesOfFieldAvailable() const
{
  return contentNotNullBase()->getTypesOfFieldAvailable();
}

std::vector< std::vector<std::pair<int,int> > > MEDFileAnyTypeField1TS::getFieldSplitedByType(const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF,
                                                                                              std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const
{
  return contentNotNullBase()->getFieldSplitedByType(mname,types,typesF,pfls,locs);
}

/*!
 * This method returns as MEDFileAnyTypeField1TS new instances as number of components in \a this.
 * The returned instances are deep copy of \a this except that for globals that are share with those contained in \a this.
 * ** WARNING ** do no forget to rename the ouput instances to avoid to write n-times in the same MED file field !
 */
std::vector< MCAuto< MEDFileAnyTypeField1TS > > MEDFileAnyTypeField1TS::splitComponents() const
{
  const MEDFileAnyTypeField1TSWithoutSDA *content(_content);
  if(!content)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TS::splitComponents : no content in this ! Unable to split components !");
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > contentsSplit=content->splitComponents();
  std::size_t sz(contentsSplit.size());
  std::vector< MCAuto< MEDFileAnyTypeField1TS > > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=shallowCpy();
      ret[i]->_content=contentsSplit[i];
    }
  return ret;
}

/*!
 * This method returns as MEDFileAnyTypeField1TS new instances as number of spatial discretizations in \a this.
 * The returned instances are shallowed copied of \a this except that for globals that are share with those contained in \a this.
 */
std::vector< MCAuto< MEDFileAnyTypeField1TS > > MEDFileAnyTypeField1TS::splitDiscretizations() const
{
  const MEDFileAnyTypeField1TSWithoutSDA *content(_content);
  if(!content)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TS::splitDiscretizations : no content in this ! Unable to split discretization !");
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > contentsSplit(content->splitDiscretizations());
  std::size_t sz(contentsSplit.size());
  std::vector< MCAuto< MEDFileAnyTypeField1TS > > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=shallowCpy();
      ret[i]->_content=contentsSplit[i];
    }
  return ret;
}

/*!
 * This method returns as MEDFileAnyTypeField1TS new instances as number of maximal number of discretization in \a this.
 * The returned instances are shallowed copied of \a this except that for globals that are share with those contained in \a this.
 */
std::vector< MCAuto< MEDFileAnyTypeField1TS > > MEDFileAnyTypeField1TS::splitMultiDiscrPerGeoTypes() const
{
  const MEDFileAnyTypeField1TSWithoutSDA *content(_content);
  if(!content)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeField1TS::splitMultiDiscrPerGeoTypes : no content in this ! Unable to split discretization !");
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > contentsSplit(content->splitMultiDiscrPerGeoTypes());
  std::size_t sz(contentsSplit.size());
  std::vector< MCAuto< MEDFileAnyTypeField1TS > > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=shallowCpy();
      ret[i]->_content=contentsSplit[i];
    }
  return ret;
}

MEDFileAnyTypeField1TS *MEDFileAnyTypeField1TS::deepCopy() const
{
  MCAuto<MEDFileAnyTypeField1TS> ret=shallowCpy();
  if((const MEDFileAnyTypeField1TSWithoutSDA *)_content)
    ret->_content=_content->deepCopy();
  ret->deepCpyGlobs(*this);
  return ret.retn();
}

int MEDFileAnyTypeField1TS::copyTinyInfoFrom(const MEDCouplingFieldDouble *field, const DataArray *arr)
{
  MCAuto<MEDCouplingFieldTemplate> ft(MEDCouplingFieldTemplate::New(*field));
  return copyTinyInfoFrom(field->timeDiscrSafe(),ft,arr);
}

int MEDFileAnyTypeField1TS::copyTinyInfoFrom(const TimeHolder *th, const MEDCouplingFieldTemplate *field, const DataArray *arr)
{
  return contentNotNullBase()->copyTinyInfoFrom(th,field,arr);
}

//= MEDFileField1TS

/*!
 * This method performs a copy with datatype modification ( float64->int32 ) of \a this. The globals information are copied
 * following the given input policy.
 *
 * \param [in] isDeepCpyGlobs - a boolean that indicates the behaviour concerning globals (profiles and localizations)
 *                            By default (true) the globals are deeply copied.
 * \return MEDFileIntField1TS * - a new object that is the result of the conversion of \a this to int32 field.
 */
MEDFileIntField1TS *MEDFileField1TS::convertToInt(bool isDeepCpyGlobs) const
{
  MCAuto<MEDFileIntField1TS> ret;
  const MEDFileAnyTypeField1TSWithoutSDA *content(_content);
  if(content)
    {
      const MEDFileField1TSWithoutSDA *contc=dynamic_cast<const MEDFileField1TSWithoutSDA *>(content);
      if(!contc)
        throw INTERP_KERNEL::Exception("MEDFileField1TS::convertToInt : the content inside this is not FLOAT64 ! This is incoherent !");
      MCAuto<MEDFileIntField1TSWithoutSDA> newc(contc->convertToInt());
      ret=static_cast<MEDFileIntField1TS *>(MEDFileAnyTypeField1TS::BuildNewInstanceFromContent((MEDFileIntField1TSWithoutSDA *)newc));
    }
  else
    ret=MEDFileIntField1TS::New();
  if(isDeepCpyGlobs)
    ret->deepCpyGlobs(*this);
  else
    ret->shallowCpyGlobs(*this);
  return ret.retn();
}

MEDFileField1TS::MEDFileField1TS(med_idt fid, bool loadAll, const MEDFileMeshes *ms)
try:MEDFileTemplateField1TS<double>(fid,loadAll,ms)
{
}
catch(INTERP_KERNEL::Exception& e)
{ throw e; }

MEDFileField1TS::MEDFileField1TS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms)
try:MEDFileTemplateField1TS<double>(fid,fieldName,loadAll,ms)
{
}
catch(INTERP_KERNEL::Exception& e)
{ throw e; }

MEDFileField1TS::MEDFileField1TS(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileMeshes *ms)
try:MEDFileTemplateField1TS<double>(fid,fieldName,iteration,order,loadAll,ms)
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
try:MEDFileTemplateField1TS<double>(other,shallowCopyOfContent)
{
}
catch(INTERP_KERNEL::Exception& e)
{ throw e; }

MEDFileField1TS *MEDFileField1TS::shallowCpy() const
{
  return new MEDFileField1TS(*this);
}

std::vector< std::vector<DataArrayDouble *> > MEDFileField1TS::getFieldSplitedByType2(const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF,
                                                                                      std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const
{
  return contentNotNull()->getFieldSplitedByType2(mname,types,typesF,pfls,locs);
}

//= MEDFileIntField1TS

MCAuto<MEDCouplingFieldDouble> MEDFileIntField1TS::ConvertFieldIntToFieldDouble(const MEDCouplingFieldInt *f)
{
  if(!f)
    throw INTERP_KERNEL::Exception("MEDFileIntField1TS::ConvertFieldIntToFieldDouble : null input field !");
  int t1,t2;
  double t0(f->getTime(t1,t2));
  std::string tu(f->getTimeUnit());
  MCAuto<MEDCouplingFieldTemplate> ft(MEDCouplingFieldTemplate::New(*f));
  MCAuto<MEDCouplingFieldDouble> ret(MEDCouplingFieldDouble::New(*ft));
  ret->setTime(t0,t1,t2); ret->setTimeUnit(tu);
  return ret;
}

//= MEDFileFloatField1TS

//= MEDFileFloatField1TS

//= MEDFileAnyTypeFieldMultiTSWithoutSDA

MEDFileAnyTypeFieldMultiTSWithoutSDA::MEDFileAnyTypeFieldMultiTSWithoutSDA()
{
}

MEDFileAnyTypeFieldMultiTSWithoutSDA::MEDFileAnyTypeFieldMultiTSWithoutSDA(const std::string& fieldName, const std::string& meshName):MEDFileFieldNameScope(fieldName,meshName)
{
}

/*!
 * \param [in] fieldId field id in C mode
 */
MEDFileAnyTypeFieldMultiTSWithoutSDA::MEDFileAnyTypeFieldMultiTSWithoutSDA(med_idt fid, int fieldId, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
{
  med_field_type typcha;
  std::string dtunitOut,meshName;
  int nbOfStep(MEDFileAnyTypeField1TS::LocateField2(fid,fieldId,false,_name,typcha,_infos,dtunitOut,meshName));
  setMeshName(meshName);
  setDtUnit(dtunitOut.c_str());
  loadStructureOrStructureAndBigArraysRecursively(fid,nbOfStep,typcha,loadAll,ms,entities);
}

MEDFileAnyTypeFieldMultiTSWithoutSDA::MEDFileAnyTypeFieldMultiTSWithoutSDA(med_idt fid, const std::string& fieldName, const std::string& meshName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
try:MEDFileFieldNameScope(fieldName,meshName),_infos(infos)
{
  setDtUnit(dtunit.c_str());
  loadStructureOrStructureAndBigArraysRecursively(fid,nbOfStep,fieldTyp,loadAll,ms,entities);
}
catch(INTERP_KERNEL::Exception& e)
{
    throw e;
}

std::size_t MEDFileAnyTypeFieldMultiTSWithoutSDA::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(_mesh_name.capacity()+_name.capacity()+_infos.capacity()*sizeof(std::string)+_time_steps.capacity()*sizeof(MCAuto<MEDFileField1TSWithoutSDA>));
  for(std::vector<std::string>::const_iterator it=_infos.begin();it!=_infos.end();it++)
    ret+=(*it).capacity();
  return ret;
}

std::vector<const BigMemoryObject *> MEDFileAnyTypeFieldMultiTSWithoutSDA::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    ret.push_back((const MEDFileAnyTypeField1TSWithoutSDA *)*it);
  return ret;
}

/*!
 * If one of the id in [ \a startIds , \a endIds ) points to a null element, there is not throw. Simply, this empty element is added as if it were not
 * NULL.
 */
MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::buildFromTimeStepIds(const int *startIds, const int *endIds) const
{
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> ret=createNew();
  ret->setInfo(_infos);
  int sz=(int)_time_steps.size();
  for(const int *id=startIds;id!=endIds;id++)
    {
      if(*id>=0 && *id<sz)
        {
          const MEDFileAnyTypeField1TSWithoutSDA *tse=_time_steps[*id];
          MCAuto<MEDFileAnyTypeField1TSWithoutSDA> tse2;
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
          throw INTERP_KERNEL::Exception(oss.str());
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
MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::buildFromTimeStepIds2(int bg, int end, int step) const
{
  static const char msg[]="MEDFileAnyTypeFieldMultiTSWithoutSDA::buildFromTimeStepIds2";
  int nbOfEntriesToKeep=DataArrayInt::GetNumberOfItemGivenBESRelative(bg,end,step,msg);
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> ret=createNew();
  ret->setInfo(_infos);
  int sz=(int)_time_steps.size();
  int j=bg;
  for(int i=0;i<nbOfEntriesToKeep;i++,j+=step)
    {
      if(j>=0 && j<sz)
        {
          const MEDFileAnyTypeField1TSWithoutSDA *tse=_time_steps[j];
          MCAuto<MEDFileAnyTypeField1TSWithoutSDA> tse2;
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
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  if(ret->getNumberOfTS()>0)
    ret->synchronizeNameScope();
  ret->copyNameScope(*this);
  return ret.retn();
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::partOfThisLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const
{
  int id=0;
  MCAuto<DataArrayInt> ids=DataArrayInt::New(); ids->alloc(0,1);
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,id++)
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

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::partOfThisNotLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const
{
  int id=0;
  MCAuto<DataArrayInt> ids=DataArrayInt::New(); ids->alloc(0,1);
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,id++)
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

bool MEDFileAnyTypeFieldMultiTSWithoutSDA::presenceOfStructureElements() const
{
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it).isNotNull())
      if((*it)->presenceOfStructureElements())
        return true;
  return false;
}

bool MEDFileAnyTypeFieldMultiTSWithoutSDA::onlyStructureElements() const
{
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it).isNotNull())
      if(!(*it)->onlyStructureElements())
        return false;
  return true;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::killStructureElements()
{
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > ret;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
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
  _time_steps=ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::keepOnlyStructureElements()
{
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > ret;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it).isNotNull())
      {
        if((*it)->presenceOfStructureElements())
          {
            if(!(*it)->onlyStructureElements())
              (*it)->keepOnlyStructureElements();
            ret.push_back(*it);
          }
      }
  _time_steps=ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::keepOnlyOnSE(const std::string& seName)
{
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > ret;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it).isNotNull())
      (*it)->keepOnlyOnSE(seName);
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::getMeshSENames(std::vector< std::pair<std::string,std::string> >& ps) const
{
  std::vector< std::pair<std::string,std::string> > ps2;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it).isNotNull())
      {
        (*it)->getMeshSENames(ps2);
        break;
      }
  if(ps2.empty())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::getMeshSENames : this appears to not contain SE only !");
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it).isNotNull())
      {
        std::vector< std::pair<std::string,std::string> > ps3;
        (*it)->getMeshSENames(ps3);
        if(ps2!=ps3)
          throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::getMeshSENames : For the moment only homogeneous SE def through time managed !");
      }
  for(std::vector< std::pair<std::string,std::string> >::const_iterator it=ps2.begin();it!=ps2.end();it++)
    {
      std::vector< std::pair<std::string,std::string> >::iterator it2(std::find(ps.begin(),ps.end(),*it));
      if(it2==ps.end())
        ps.push_back(*it);
    }
}

bool MEDFileAnyTypeFieldMultiTSWithoutSDA::presenceOfMultiDiscPerGeoType() const
{
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *cur(*it);
      if(!cur)
        continue;
      if(cur->presenceOfMultiDiscPerGeoType())
        return true;
    }
  return false;
}

const std::vector<std::string>& MEDFileAnyTypeFieldMultiTSWithoutSDA::getInfo() const
{
  return _infos;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::setInfo(const std::vector<std::string>& info)
{
  _infos=info;
}

int MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepPos(int iteration, int order) const
{
  int ret=0;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,ret++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *pt(*it);
      if(pt->isDealingTS(iteration,order))
        return ret;
    }
  std::ostringstream oss; oss << "MEDFileFieldMultiTS::getTimeStepPos : Muli timestep field on time (" << iteration << "," << order << ") does not exist ! Available (iteration,order) are :\n";
  std::vector< std::pair<int,int> > vp=getIterations();
  for(std::vector< std::pair<int,int> >::const_iterator it2=vp.begin();it2!=vp.end();it2++)
    oss << "(" << (*it2).first << "," << (*it2).second << ") ";
  throw INTERP_KERNEL::Exception(oss.str());
}

const MEDFileAnyTypeField1TSWithoutSDA& MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepEntry(int iteration, int order) const
{
  return *_time_steps[getTimeStepPos(iteration,order)];
}

MEDFileAnyTypeField1TSWithoutSDA& MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepEntry(int iteration, int order)
{
  return *_time_steps[getTimeStepPos(iteration,order)];
}

bool MEDFileAnyTypeFieldMultiTSWithoutSDA::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab)
{
  bool ret(false);
  for(std::vector< std::pair<std::string,std::string> >::const_iterator it=modifTab.begin();it!=modifTab.end();it++)
    {
      if((*it).first==getMeshName())
        {
          setMeshName((*it).second);
          ret=true;
        }
    }
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
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
DataArray *MEDFileAnyTypeFieldMultiTSWithoutSDA::getUndergroundDataArray(int iteration, int order) const
{
  return getTimeStepEntry(iteration,order).getUndergroundDataArray();
}

/*!
 * See doc at MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt
 */
DataArray *MEDFileAnyTypeFieldMultiTSWithoutSDA::getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const
{
  return getTimeStepEntry(iteration,order).getUndergroundDataArrayExt(entries);
}

bool MEDFileAnyTypeFieldMultiTSWithoutSDA::renumberEntitiesLyingOnMesh(const std::string& meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N,
                                                                       MEDFileFieldGlobsReal& glob)
{
  bool ret=false;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileAnyTypeField1TSWithoutSDA *f1ts(*it);
      if(f1ts)
        ret=f1ts->renumberEntitiesLyingOnMesh(meshName,oldCode,newCode,renumO2N,glob) || ret;
    }
  return ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::accept(MEDFileFieldVisitor& visitor) const
{
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    if((*it).isNotNull())
      {
        visitor.newTimeStepEntry(*it);
        (*it)->accept(visitor);
        visitor.endTimeStepEntry(*it);
      }
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
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,i++)
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

std::vector< std::pair<int,int> > MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeSteps(std::vector<double>& ret1) const
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
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  return ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::pushBackTimeStep(MCAuto<MEDFileAnyTypeField1TSWithoutSDA>& tse)
{
  MEDFileAnyTypeField1TSWithoutSDA *tse2(tse);
  if(!tse2)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::pushBackTimeStep : input content object is null !");
  checkCoherencyOfType(tse2);
  if(_time_steps.empty())
    {
      setName(tse2->getName());
      setMeshName(tse2->getMeshName());
      setInfo(tse2->getInfo());
    }
  checkThatComponentsMatch(tse2->getInfo());
  if(getDtUnit().empty() && !tse->getDtUnit().empty())
    setDtUnit(tse->getDtUnit());
  _time_steps.push_back(tse);
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::synchronizeNameScope()
{
  std::size_t nbOfCompo=_infos.size();
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileAnyTypeField1TSWithoutSDA *cur=(*it);
      if(cur)
        {
          if((cur->getInfo()).size()!=nbOfCompo)
            {
              std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::synchronizeNameScope : Mismatch in the number of components of parts ! Should be " << nbOfCompo;
              oss << " ! but the field at iteration=" << cur->getIteration() << " order=" << cur->getOrder() << " has " << (cur->getInfo()).size() << " components !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
          cur->copyNameScope(*this);
        }
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::loadStructureOrStructureAndBigArraysRecursively(med_idt fid, int nbPdt, med_field_type fieldTyp, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
{
  _time_steps.resize(nbPdt);
  for(int i=0;i<nbPdt;i++)
    {
      std::vector< std::pair<int,int> > ts;
      med_int numdt=0,numo=0;
      med_float dt=0.0;
      MEDFILESAFECALLERRD0(MEDfieldComputingStepInfo,(fid,_name.c_str(),i+1,&numdt,&numo,&dt));
      switch(fieldTyp)
      {
        case MED_FLOAT64:
          {
            _time_steps[i]=MEDFileField1TSWithoutSDA::New(getName(),getMeshName(),i+1,numdt,numo,_infos);
            break;
          }
        case MED_INT32:
          {
            _time_steps[i]=MEDFileIntField1TSWithoutSDA::New(getName(),getMeshName(),i+1,numdt,numo,_infos);
            break;
          }
        case MED_FLOAT32:
          {
            _time_steps[i]=MEDFileFloatField1TSWithoutSDA::New(getName(),getMeshName(),i+1,numdt,numo,_infos);
            break;
          }
        case MED_INT:
          {
            if(sizeof(med_int)==sizeof(int))
              {
                _time_steps[i]=MEDFileIntField1TSWithoutSDA::New(getName(),getMeshName(),i+1,numdt,numo,_infos);
                break;
              }
          }
        default:
          throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::loadStructureOrStructureAndBigArraysRecursively : managed field type are : FLOAT64, INT32, FLOAT32 !");
      }
      if(loadAll)
        _time_steps[i]->loadStructureAndBigArraysRecursively(fid,*this,ms,entities);
      else
        _time_steps[i]->loadOnlyStructureOfDataRecursively(fid,*this,ms,entities);
      synchronizeNameScope();
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::writeLL(med_idt fid, const MEDFileWritable& opts) const
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
  MEDFILESAFECALLERWR0(MEDfieldCr,(fid,_name.c_str(),getMEDFileFieldType(),nbComp,comp,unit,getDtUnit().c_str(),getMeshName().c_str()));
  int nbOfTS=_time_steps.size();
  for(int i=0;i<nbOfTS;i++)
    _time_steps[i]->writeLL(fid,opts,*this);
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc)
{
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileAnyTypeField1TSWithoutSDA *elt(*it);
      if(elt)
        elt->loadBigArraysRecursively(fid,nasc);
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::loadBigArraysRecursivelyIfNecessary(med_idt fid, const MEDFileFieldNameScope& nasc)
{
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileAnyTypeField1TSWithoutSDA *elt(*it);
      if(elt)
        elt->loadBigArraysRecursivelyIfNecessary(fid,nasc);
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::unloadArrays()
{
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      MEDFileAnyTypeField1TSWithoutSDA *elt(*it);
      if(elt)
        elt->unloadArrays();
    }
}

int MEDFileAnyTypeFieldMultiTSWithoutSDA::getNumberOfTS() const
{
  return _time_steps.size();
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::eraseEmptyTS()
{
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA>  > newTS;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *tmp=(*it);
      if(tmp)
        newTS.push_back(*it);
    }
  _time_steps=newTS;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::eraseTimeStepIds(const int *startIds, const int *endIds)
{
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > newTS;
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
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  for(int iii=0;iii<maxId;iii++)
    if(idsToDel.find(iii)==idsToDel.end())
      newTS.push_back(_time_steps[iii]);
  _time_steps=newTS;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::eraseTimeStepIds2(int bg, int end, int step)
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
  std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > newTS;
  for(std::size_t i=0;i<sz;i++)
    if(b[i])
      newTS.push_back(_time_steps[i]);
  _time_steps=newTS;
}

int MEDFileAnyTypeFieldMultiTSWithoutSDA::getPosOfTimeStep(int iteration, int order) const
{
  int ret=0;
  std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::getPosOfTimeStep : No such time step (" << iteration << "," << order << ") !\nPossibilities are : "; 
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,ret++)
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
  throw INTERP_KERNEL::Exception(oss.str());
}

int MEDFileAnyTypeFieldMultiTSWithoutSDA::getPosGivenTime(double time, double eps) const
{
  int ret=0;
  std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::getPosGivenTime : No such time step " << time << "! \nPossibilities are : ";
  oss.precision(15);
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA>  >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,ret++)
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
  throw INTERP_KERNEL::Exception(oss.str());
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
int MEDFileAnyTypeFieldMultiTSWithoutSDA::getNonEmptyLevels(int iteration, int order, const std::string& mname, std::vector<int>& levs) const
{
  return getTimeStepEntry(iteration,order).getNonEmptyLevels(mname,levs);
}

const MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2(int pos) const
{
  if(pos<0 || pos>=(int)_time_steps.size())
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << " whereas should be in [0," << _time_steps.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  const MEDFileAnyTypeField1TSWithoutSDA *item=_time_steps[pos];
  if(item==0)
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << ", this pos id exists but the underlying Field1TS is null !";
      oss << "\nTry to use following method eraseEmptyTS !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  return item;
}

MEDFileAnyTypeField1TSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2(int pos)
{
  if(pos<0 || pos>=(int)_time_steps.size())
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << " whereas should be in [0," << _time_steps.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  MEDFileAnyTypeField1TSWithoutSDA *item=_time_steps[pos];
  if(item==0)
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::getTimeStepAtPos2 : request for pos #" << pos << ", this pos id exists but the underlying Field1TS is null !";
      oss << "\nTry to use following method eraseEmptyTS !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  return item;
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTSWithoutSDA::getPflsReallyUsed2() const
{
  std::vector<std::string> ret;
  std::set<std::string> ret2;
  for(std::vector< MCAuto< MEDFileAnyTypeField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
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
  for(std::vector< MCAuto< MEDFileAnyTypeField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
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
  for(std::vector< MCAuto< MEDFileAnyTypeField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getPflsReallyUsedMulti2();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

std::vector<std::string> MEDFileAnyTypeFieldMultiTSWithoutSDA::getLocsReallyUsedMulti2() const
{
  std::vector<std::string> ret;
  for(std::vector< MCAuto< MEDFileAnyTypeField1TSWithoutSDA > >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    {
      std::vector<std::string> tmp=(*it)->getLocsReallyUsedMulti2();
      ret.insert(ret.end(),tmp.begin(),tmp.end());
    }
  return ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::changePflsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  for(std::vector< MCAuto< MEDFileAnyTypeField1TSWithoutSDA > >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    (*it)->changePflsRefsNamesGen2(mapOfModif);
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::changeLocsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  for(std::vector< MCAuto< MEDFileAnyTypeField1TSWithoutSDA > >::iterator it=_time_steps.begin();it!=_time_steps.end();it++)
    (*it)->changeLocsRefsNamesGen2(mapOfModif);
}

std::vector< std::vector<TypeOfField> > MEDFileAnyTypeFieldMultiTSWithoutSDA::getTypesOfFieldAvailable() const
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
std::vector< std::vector< std::pair<int,int> > > MEDFileAnyTypeFieldMultiTSWithoutSDA::getFieldSplitedByType(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const
{
  return getTimeStepEntry(iteration,order).getFieldSplitedByType(mname,types,typesF,pfls,locs);
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTSWithoutSDA::deepCopy() const
{
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> ret=shallowCpy();
  std::size_t i=0;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,i++)
    {
      if((const MEDFileAnyTypeField1TSWithoutSDA *)*it)
        ret->_time_steps[i]=(*it)->deepCopy();
    }
  return ret.retn();
}

std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > MEDFileAnyTypeFieldMultiTSWithoutSDA::splitComponents() const
{
  std::size_t sz(_infos.size()),sz2(_time_steps.size());
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > ret(sz);
  std::vector< std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > > ts(sz2);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=shallowCpy();
      ret[i]->_infos.resize(1); ret[i]->_infos[0]=_infos[i];
    }
  for(std::size_t i=0;i<sz2;i++)
    {
      std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > ret1=_time_steps[i]->splitComponents();
      if(ret1.size()!=sz)
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::splitComponents : At rank #" << i << " number of components is " << ret1.size() << " whereas it should be for all time steps " << sz << " !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      ts[i]=ret1;
    }
  for(std::size_t i=0;i<sz;i++)
    for(std::size_t j=0;j<sz2;j++)
      ret[i]->_time_steps[j]=ts[j][i];
  return ret;
}

/*!
 * This method splits into discretization each time steps in \a this.
 * ** WARNING ** the returned instances are not compulsary defined on the same time steps series !
 */
std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > MEDFileAnyTypeFieldMultiTSWithoutSDA::splitDiscretizations() const
{
  std::size_t sz(_time_steps.size());
  std::vector< std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > > items(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *timeStep(_time_steps[i]);
      if(!timeStep)
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::splitDiscretizations : time step #" << i << " is null !"; 
          throw INTERP_KERNEL::Exception(oss.str());
        }
      items[i]=timeStep->splitDiscretizations();  
    }
  //
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > ret;
  std::vector< std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > > ret2;
  std::vector< TypeOfField > types;
  for(std::vector< std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > >::const_iterator it0=items.begin();it0!=items.end();it0++)
    for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it1=(*it0).begin();it1!=(*it0).end();it1++)
      {
        std::vector<TypeOfField> ts=(*it1)->getTypesOfFieldAvailable();
        if(ts.size()!=1)
          throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::splitDiscretizations : it appears that the splitting of MEDFileAnyTypeField1TSWithoutSDA::splitDiscretizations has returned invalid result !");
        std::vector< TypeOfField >::iterator it2=std::find(types.begin(),types.end(),ts[0]);
        if(it2==types.end())
          types.push_back(ts[0]);
      }
  ret.resize(types.size()); ret2.resize(types.size());
  for(std::vector< std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > >::const_iterator it0=items.begin();it0!=items.end();it0++)
    for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it1=(*it0).begin();it1!=(*it0).end();it1++)
      {
        TypeOfField typ=(*it1)->getTypesOfFieldAvailable()[0];
        std::size_t pos=std::distance(types.begin(),std::find(types.begin(),types.end(),typ));
        ret2[pos].push_back(*it1);
      }
  for(std::size_t i=0;i<types.size();i++)
    {
      MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> elt(createNew());
      for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::iterator it1=ret2[i].begin();it1!=ret2[i].end();it1++)
        elt->pushBackTimeStep(*it1);//also updates infos in elt
      ret[i]=elt;
      elt->MEDFileFieldNameScope::operator=(*this);
    }
  return ret;
}

/*!
 * Contrary to splitDiscretizations method this method makes the hypothesis that the times series are **NOT** impacted by the splitting of multi discretization.
 */
std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > MEDFileAnyTypeFieldMultiTSWithoutSDA::splitMultiDiscrPerGeoTypes() const
{
  std::size_t sz(_time_steps.size());
  std::vector< std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > > items(sz);
  std::size_t szOut(std::numeric_limits<std::size_t>::max());
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *timeStep(_time_steps[i]);
      if(!timeStep)
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::splitMultiDiscrPerGeoTypes : time step #" << i << " is null !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      items[i]=timeStep->splitMultiDiscrPerGeoTypes();
      if(szOut==std::numeric_limits<std::size_t>::max())
        szOut=items[i].size();
      else
        if(items[i].size()!=szOut)
          throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::splitMultiDiscrPerGeoTypes : The splitting per discretization is expected to be same among time steps !");
    }
  if(szOut==std::numeric_limits<std::size_t>::max())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::splitMultiDiscrPerGeoTypes : empty field !");
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > ret(szOut);
  for(std::size_t i=0;i<szOut;i++)
    {
      MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> elt(createNew());
      for(std::size_t j=0;j<sz;j++)
        elt->pushBackTimeStep(items[j][i]);
      ret[i]=elt;
      elt->MEDFileFieldNameScope::operator=(*this);
    }
  return ret;
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::copyTinyInfoFrom(const MEDCouplingFieldDouble *field, const DataArray *arr)
{
  setName(field->getName());
  if(field->getMesh())
    setMeshName(field->getMesh()->getName());
  if(_name.empty())
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::copyTinyInfoFrom : unsupported fields with no name in MED file !");
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::copyTinyInfoFrom : no array set !");
  _infos=arr->getInfoOnComponents();
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::checkCoherencyOfTinyInfo(const MEDCouplingFieldDouble *field, const DataArray *arr) const
{
  static const char MSG[]="MEDFileFieldMultiTSWithoutSDA::checkCoherencyOfTinyInfo : invalid ";
  if(_name!=field->getName())
    {
      std::ostringstream oss; oss << MSG << "name ! should be \"" << _name;
      oss << "\" and it is set in input field to \"" << field->getName() << "\" !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::checkCoherencyOfTinyInfo : no array set !");
  checkThatComponentsMatch(arr->getInfoOnComponents());
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::checkThatComponentsMatch(const std::vector<std::string>& compos) const
{
  static const char MSG[]="MEDFileFieldMultiTSWithoutSDA::checkThatComponentsMatch : ";
  if(getInfo().size()!=compos.size())
    {
      std::ostringstream oss; oss << MSG << "mismatch of number of components between this (" << getInfo().size() << ") and ";
      oss << " number of components of element to append (" << compos.size() << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  if(_infos!=compos)
    {
      std::ostringstream oss; oss << MSG << "components have same size but are different ! should be \"";
      std::copy(_infos.begin(),_infos.end(),std::ostream_iterator<std::string>(oss,", "));
      oss << " But compo in input fields are : ";
      std::copy(compos.begin(),compos.end(),std::ostream_iterator<std::string>(oss,", "));
      oss << " !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::checkThatNbOfCompoOfTSMatchThis() const
{
  std::size_t sz=_infos.size();
  int j=0;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,j++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *elt(*it);
      if(elt)
        if(elt->getInfo().size()!=sz)
          {
            std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::checkThatNbOfCompoOfTSMatchThis : At pos #" << j << " the number of components is equal to ";
            oss << elt->getInfo().size() << " whereas it is expected to be equal to " << sz << " !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
    }
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::appendFieldNoProfileSBT : input field is NULL !");
  if(!_time_steps.empty())
    checkCoherencyOfTinyInfo(field,arr);
  MEDFileAnyTypeField1TSWithoutSDA *objC(createNew1TSWithoutSDAEmptyInstance());
  MCAuto<MEDFileAnyTypeField1TSWithoutSDA> obj(objC);
  {
    MCAuto<MEDCouplingFieldTemplate> ft(MEDCouplingFieldTemplate::New(*field));
    objC->setFieldNoProfileSBT(field->timeDiscrSafe(),ft,arr,glob,*this);
  }
  copyTinyInfoFrom(field,arr);
  _time_steps.push_back(obj);
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::appendFieldProfile(const MEDCouplingFieldDouble *field, const DataArray *arr, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob)
{
  if(!field)
    throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTSWithoutSDA::appendFieldNoProfileSBT : input field is NULL !");
  if(!_time_steps.empty())
    checkCoherencyOfTinyInfo(field,arr);
  MEDFileAnyTypeField1TSWithoutSDA *objC=createNew1TSWithoutSDAEmptyInstance();
  MCAuto<MEDFileAnyTypeField1TSWithoutSDA> obj(objC);
  {
    MCAuto<MEDCouplingFieldTemplate> ft(MEDCouplingFieldTemplate::NewWithoutCheck(*field));
    objC->setFieldProfile(field->timeDiscrSafe(),ft,arr,mesh,meshDimRelToMax,profile,glob,*this);
  }
  copyTinyInfoFrom(field,arr);
  setMeshName(objC->getMeshName());
  _time_steps.push_back(obj);
}

void MEDFileAnyTypeFieldMultiTSWithoutSDA::setIteration(int i, MCAuto<MEDFileAnyTypeField1TSWithoutSDA> ts)
{
  int sz=(int)_time_steps.size();
  if(i<0 || i>=sz)
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::setIteration : trying to set element at place #" << i << " should be in [0," << sz << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  const MEDFileAnyTypeField1TSWithoutSDA *tsPtr(ts);
  if(tsPtr)
    {
      if(tsPtr->getNumberOfComponents()!=(int)_infos.size())
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTSWithoutSDA::setIteration : trying to set element with " << tsPtr->getNumberOfComponents() << " components ! Should be " << _infos.size() <<  " !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  _time_steps[i]=ts;
}

//= MEDFileFieldMultiTSWithoutSDA

/*!
 * entry point for users that want to iterate into MEDFile DataStructure with a reduced overhead because output arrays are extracted (created) specially
 * for the call of this method. That's why the DataArrayDouble instance in returned vector of vector should be dealed by the caller.
 */
std::vector< std::vector<DataArrayDouble *> > MEDFileFieldMultiTSWithoutSDA::getFieldSplitedByType2(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const
{
  const MEDFileAnyTypeField1TSWithoutSDA& myF1TS=getTimeStepEntry(iteration,order);
  const MEDFileField1TSWithoutSDA *myF1TSC=dynamic_cast<const MEDFileField1TSWithoutSDA *>(&myF1TS);
  if(!myF1TSC)
    throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::getFieldSplitedByType2 : mismatch of type of field expecting FLOAT64 !");
  return myF1TSC->getFieldSplitedByType2(mname,types,typesF,pfls,locs);
}

MEDFileIntFieldMultiTSWithoutSDA *MEDFileFieldMultiTSWithoutSDA::convertToInt() const
{
  MCAuto<MEDFileIntFieldMultiTSWithoutSDA> ret(new MEDFileIntFieldMultiTSWithoutSDA);
  ret->MEDFileAnyTypeFieldMultiTSWithoutSDA::operator =(*this);
  int i=0;
  for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=_time_steps.begin();it!=_time_steps.end();it++,i++)
    {
      const MEDFileAnyTypeField1TSWithoutSDA *eltToConv(*it);
      if(eltToConv)
        {
          const MEDFileField1TSWithoutSDA *eltToConvC=dynamic_cast<const MEDFileField1TSWithoutSDA *>(eltToConv);
          if(!eltToConvC)
            throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::convertToInt : presence of an invalid 1TS type ! Should be of type FLOAT64 !");
          MCAuto<MEDFileAnyTypeField1TSWithoutSDA> elt=eltToConvC->convertToInt();
          ret->setIteration(i,elt);
        }
    }
  return ret.retn();
}

//= MEDFileAnyTypeFieldMultiTS

MEDFileAnyTypeFieldMultiTS::MEDFileAnyTypeFieldMultiTS()
{
}

MEDFileAnyTypeFieldMultiTS::MEDFileAnyTypeFieldMultiTS(med_idt fid, bool loadAll, const MEDFileMeshes *ms)
try:MEDFileFieldGlobsReal(fid)
{
  _content=BuildContentFrom(fid,loadAll,ms);
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
{
    throw e;
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTS::BuildContentFrom(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
{
  med_field_type typcha;
  std::vector<std::string> infos;
  std::string dtunit;
  std::string meshName;
  int i(-1);
  MEDFileAnyTypeField1TS::LocateField(fid,fieldName,i,typcha,infos,dtunit,meshName);
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> ret;
  switch(typcha)
  {
    case MED_FLOAT64:
      {
        ret=new MEDFileFieldMultiTSWithoutSDA(fid,i,loadAll,ms,entities);
        break;
      }
    case MED_INT32:
      {
        ret=new MEDFileIntFieldMultiTSWithoutSDA(fid,i,loadAll,ms,entities);
        break;
      }
    case MED_FLOAT32:
      {
        ret=new MEDFileFloatFieldMultiTSWithoutSDA(fid,i,loadAll,ms,entities);
        break;
      }
    case MED_INT:
      {
        if(sizeof(med_int)==sizeof(int))
          {
            ret=new MEDFileIntFieldMultiTSWithoutSDA(fid,i,loadAll,ms,entities);
            break;
          }
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::BuildContentFrom(fid,fieldName) : file \'" << FileNameFromFID(fid) << "\' contains field with name \'" << fieldName << "\' but the type of field is not in [MED_FLOAT64, MED_INT32, MED_FLOAT32] !";
        throw INTERP_KERNEL::Exception(oss.str());
      }
  }
  ret->setMeshName(meshName);
  ret->setDtUnit(dtunit.c_str());
  return ret.retn();
}

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTS::BuildContentFrom(med_idt fid, bool loadAll, const MEDFileMeshes *ms)
{
  med_field_type typcha;
  //
  std::vector<std::string> infos;
  std::string dtunit,fieldName,meshName;
  MEDFileAnyTypeField1TS::LocateField2(fid,0,true,fieldName,typcha,infos,dtunit,meshName);
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> ret;
  switch(typcha)
  {
    case MED_FLOAT64:
      {
        ret=new MEDFileFieldMultiTSWithoutSDA(fid,0,loadAll,ms,0);
        break;
      }
    case MED_INT32:
      {
        ret=new MEDFileIntFieldMultiTSWithoutSDA(fid,0,loadAll,ms,0);
        break;
      }
    case MED_FLOAT32:
      {
        ret=new MEDFileFloatFieldMultiTSWithoutSDA(fid,0,loadAll,ms,0);
        break;
      }
    case MED_INT:
      {
        if(sizeof(med_int)==sizeof(int))
          {
            ret=new MEDFileIntFieldMultiTSWithoutSDA(fid,0,loadAll,ms,0);
            break;
          }
      }
    default:
      {
        std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::BuildContentFrom(fid) : file \'" << FileNameFromFID(fid) << "\' contains field with name \'" << fieldName << "\' but the type of the first field is not in [MED_FLOAT64, MED_INT32, MED_FLOAT32] !";
        throw INTERP_KERNEL::Exception(oss.str());
      }
  }
  ret->setMeshName(meshName);
  ret->setDtUnit(dtunit.c_str());
  return ret.retn();
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::BuildNewInstanceFromContent(MEDFileAnyTypeFieldMultiTSWithoutSDA *c)
{
  if(!c)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::BuildNewInstanceFromContent : empty content in input : unable to build a new instance !");
  if(dynamic_cast<const MEDFileFieldMultiTSWithoutSDA *>(c))
    {
      MCAuto<MEDFileFieldMultiTS> ret(MEDFileFieldMultiTS::New());
      ret->_content=c;  c->incrRef();
      return ret.retn();
    }
  if(dynamic_cast<const MEDFileIntFieldMultiTSWithoutSDA *>(c))
    {
      MCAuto<MEDFileIntFieldMultiTS> ret(MEDFileIntFieldMultiTS::New());
      ret->_content=c;  c->incrRef();
      return ret.retn();
    }
  if(dynamic_cast<const MEDFileFloatFieldMultiTSWithoutSDA *>(c))
    {
      MCAuto<MEDFileFloatFieldMultiTS> ret(MEDFileFloatFieldMultiTS::New());
      ret->_content=c;  c->incrRef();
      return ret.retn();
    }
  throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::BuildNewInstanceFromContent : internal error ! a content of type different from FLOAT64 FLOAT32 and INT32 has been built but not intercepted !");
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::BuildNewInstanceFromContent(MEDFileAnyTypeFieldMultiTSWithoutSDA *c, med_idt fid)
{
  MEDFileAnyTypeFieldMultiTS *ret(BuildNewInstanceFromContent(c));
  std::string fileName(FileNameFromFID(fid));
  ret->setFileName(fileName);
  return ret;
}

MEDFileAnyTypeFieldMultiTS::MEDFileAnyTypeFieldMultiTS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
try:MEDFileFieldGlobsReal(fid)
{
  _content=BuildContentFrom(fid,fieldName,loadAll,ms,entities);
  loadGlobals(fid);
}
catch(INTERP_KERNEL::Exception& e)
{
    throw e;
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
MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::New(const std::string& fileName, bool loadAll)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return New(fid,loadAll);
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::New(med_idt fid, bool loadAll)
{
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> c(BuildContentFrom(fid,loadAll,0));
  MCAuto<MEDFileAnyTypeFieldMultiTS> ret(BuildNewInstanceFromContent(c,fid));
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
MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::New(const std::string& fileName, const std::string& fieldName, bool loadAll)
{
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
  return New(fid,fieldName,loadAll);
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::New(med_idt fid, const std::string& fieldName, bool loadAll)
{
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> c(BuildContentFrom(fid,fieldName,loadAll,0,0));
  MCAuto<MEDFileAnyTypeFieldMultiTS> ret(BuildNewInstanceFromContent(c,fid));
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

MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTS::contentNotNullBase()
{
  MEDFileAnyTypeFieldMultiTSWithoutSDA *ret=_content;
  if(!ret)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS : content is expected to be not null !");
  return ret;
}

const MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileAnyTypeFieldMultiTS::contentNotNullBase() const
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

void MEDFileAnyTypeFieldMultiTS::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  contentNotNullBase()->changePflsRefsNamesGen2(mapOfModif);
}

void MEDFileAnyTypeFieldMultiTS::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif)
{
  contentNotNullBase()->changeLocsRefsNamesGen2(mapOfModif);
}

int MEDFileAnyTypeFieldMultiTS::getNumberOfTS() const
{
  return contentNotNullBase()->getNumberOfTS();
}

void MEDFileAnyTypeFieldMultiTS::eraseEmptyTS()
{
  contentNotNullBase()->eraseEmptyTS();
}

void MEDFileAnyTypeFieldMultiTS::eraseTimeStepIds(const int *startIds, const int *endIds)
{
  contentNotNullBase()->eraseTimeStepIds(startIds,endIds);
}

void MEDFileAnyTypeFieldMultiTS::eraseTimeStepIds2(int bg, int end, int step)
{
  contentNotNullBase()->eraseTimeStepIds2(bg,end,step);
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::buildSubPart(const int *startIds, const int *endIds) const
{
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> c=contentNotNullBase()->buildFromTimeStepIds(startIds,endIds);
  MCAuto<MEDFileAnyTypeFieldMultiTS> ret=shallowCpy();
  ret->_content=c;
  return ret.retn();
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::buildSubPartSlice(int bg, int end, int step) const
{
  MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> c=contentNotNullBase()->buildFromTimeStepIds2(bg,end,step);
  MCAuto<MEDFileAnyTypeFieldMultiTS> ret=shallowCpy();
  ret->_content=c;
  return ret.retn();
}

std::vector< std::pair<int,int> > MEDFileAnyTypeFieldMultiTS::getIterations() const
{
  return contentNotNullBase()->getIterations();
}

void MEDFileAnyTypeFieldMultiTS::pushBackTimeSteps(const std::vector<MEDFileAnyTypeField1TS *>& f1ts)
{
  for(std::vector<MEDFileAnyTypeField1TS *>::const_iterator it=f1ts.begin();it!=f1ts.end();it++)
    pushBackTimeStep(*it);
}

void MEDFileAnyTypeFieldMultiTS::pushBackTimeSteps(MEDFileAnyTypeFieldMultiTS *fmts)
{
  if(!fmts)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::pushBackTimeSteps : Input fmts is NULL !");
  int nbOfTS(fmts->getNumberOfTS());
  for(int i=0;i<nbOfTS;i++)
    {
      MCAuto<MEDFileAnyTypeField1TS> elt(fmts->getTimeStepAtPos(i));
      pushBackTimeStep(elt);
    }
}

void MEDFileAnyTypeFieldMultiTS::pushBackTimeStep(MEDFileAnyTypeField1TS *f1ts)
{
  if(!f1ts)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::pushBackTimeStep : input pointer is NULL !");
  checkCoherencyOfType(f1ts);
  f1ts->incrRef();
  MCAuto<MEDFileAnyTypeField1TS> f1tsSafe(f1ts);
  MEDFileAnyTypeField1TSWithoutSDA *c=f1ts->contentNotNullBase();
  c->incrRef();
  MCAuto<MEDFileAnyTypeField1TSWithoutSDA> cSafe(c);
  if(!((MEDFileAnyTypeFieldMultiTSWithoutSDA *)_content))
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTSWithoutSDA::pushBackTimeStep : no content in this !");
  _content->pushBackTimeStep(cSafe);
  appendGlobs(*f1ts,1e-12);
}

void MEDFileAnyTypeFieldMultiTS::synchronizeNameScope()
{
  contentNotNullBase()->synchronizeNameScope();
}

int MEDFileAnyTypeFieldMultiTS::getPosOfTimeStep(int iteration, int order) const
{
  return contentNotNullBase()->getPosOfTimeStep(iteration,order);
}

int MEDFileAnyTypeFieldMultiTS::getPosGivenTime(double time, double eps) const
{
  return contentNotNullBase()->getPosGivenTime(time,eps);
}

int MEDFileAnyTypeFieldMultiTS::getNonEmptyLevels(int iteration, int order, const std::string& mname, std::vector<int>& levs) const
{
  return contentNotNullBase()->getNonEmptyLevels(iteration,order,mname,levs);
}

std::vector< std::vector<TypeOfField> > MEDFileAnyTypeFieldMultiTS::getTypesOfFieldAvailable() const
{
  return contentNotNullBase()->getTypesOfFieldAvailable();
}

std::vector< std::vector< std::pair<int,int> > > MEDFileAnyTypeFieldMultiTS::getFieldSplitedByType(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const
{
  return contentNotNullBase()->getFieldSplitedByType(iteration,order,mname,types,typesF,pfls,locs);
}

std::string MEDFileAnyTypeFieldMultiTS::getName() const
{
  return contentNotNullBase()->getName();
}

void MEDFileAnyTypeFieldMultiTS::setName(const std::string& name)
{
  contentNotNullBase()->setName(name);
}

std::string MEDFileAnyTypeFieldMultiTS::getDtUnit() const
{
  return contentNotNullBase()->getDtUnit();
}

void MEDFileAnyTypeFieldMultiTS::setDtUnit(const std::string& dtUnit)
{
  contentNotNullBase()->setDtUnit(dtUnit);
}

void MEDFileAnyTypeFieldMultiTS::simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const
{
  contentNotNullBase()->simpleRepr(bkOffset,oss,fmtsId);
}

std::vector< std::pair<int,int> > MEDFileAnyTypeFieldMultiTS::getTimeSteps(std::vector<double>& ret1) const
{
  return contentNotNullBase()->getTimeSteps(ret1);
}

std::string MEDFileAnyTypeFieldMultiTS::getMeshName() const
{
  return contentNotNullBase()->getMeshName();
}

void MEDFileAnyTypeFieldMultiTS::setMeshName(const std::string& newMeshName)
{
  contentNotNullBase()->setMeshName(newMeshName);
}

bool MEDFileAnyTypeFieldMultiTS::changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab)
{
  return contentNotNullBase()->changeMeshNames(modifTab);
}

const std::vector<std::string>& MEDFileAnyTypeFieldMultiTS::getInfo() const
{
  return contentNotNullBase()->getInfo();
}

bool MEDFileAnyTypeFieldMultiTS::presenceOfMultiDiscPerGeoType() const
{
  return contentNotNullBase()->presenceOfMultiDiscPerGeoType();
}

void MEDFileAnyTypeFieldMultiTS::setInfo(const std::vector<std::string>& info)
{
  return contentNotNullBase()->setInfo(info);
}

int MEDFileAnyTypeFieldMultiTS::getNumberOfComponents() const
{
  const std::vector<std::string> ret=getInfo();
  return (int)ret.size();
}

void MEDFileAnyTypeFieldMultiTS::writeLL(med_idt fid) const
{
  writeGlobals(fid,*this);
  contentNotNullBase()->writeLL(fid,*this);
}

/*!
 * This method alloc the arrays and load potentially huge arrays contained in this field.
 * This method should be called when a MEDFileAnyTypeFieldMultiTS::New constructor has been with false as the last parameter.
 * This method can be also called to refresh or reinit values from a file.
 * 
 * \throw If the fileName is not set or points to a non readable MED file.
 */
void MEDFileAnyTypeFieldMultiTS::loadArrays()
{
  if(getFileName().empty())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::loadArrays : the structure does not come from a file !");
  MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(getFileName()));
  contentNotNullBase()->loadBigArraysRecursively(fid,*contentNotNullBase());
}

/*!
 * This method behaves as MEDFileAnyTypeFieldMultiTS::loadArrays does, the first call, if \a this was built using a file without loading big arrays.
 * But once data loaded once, this method does nothing.
 * 
 * \throw If the fileName is not set or points to a non readable MED file.
 * \sa MEDFileAnyTypeFieldMultiTS::loadArrays, MEDFileAnyTypeFieldMultiTS::unloadArrays
 */
void MEDFileAnyTypeFieldMultiTS::loadArraysIfNecessary()
{
  if(!getFileName().empty())
    {
      MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(getFileName()));
      contentNotNullBase()->loadBigArraysRecursivelyIfNecessary(fid,*contentNotNullBase());
    }
}

/*!
 * This method releases potentially big data arrays and so returns to the same heap memory than status loaded with 'loadAll' parameter set to false.
 * \b WARNING, this method does release arrays even if \a this does not come from a load of a MED file.
 * So this method can lead to a loss of data. If you want to unload arrays safely call MEDFileAnyTypeFieldMultiTS::unloadArraysWithoutDataLoss instead.
 * 
 * \sa MEDFileAnyTypeFieldMultiTS::loadArrays, MEDFileAnyTypeFieldMultiTS::loadArraysIfNecessary, MEDFileAnyTypeFieldMultiTS::unloadArraysWithoutDataLoss
 */
void MEDFileAnyTypeFieldMultiTS::unloadArrays()
{
  contentNotNullBase()->unloadArrays();
}

/*!
 * This method potentially releases big data arrays if \a this is coming from a file. If \a this has been built from scratch this method will have no effect.
 * This method is the symetrical method of MEDFileAnyTypeFieldMultiTS::loadArraysIfNecessary.
 * This method is useful to reduce \b safely amount of heap memory necessary for \a this by using MED file as database.
 * 
 * \sa MEDFileAnyTypeFieldMultiTS::loadArraysIfNecessary
 */
void MEDFileAnyTypeFieldMultiTS::unloadArraysWithoutDataLoss()
{
  if(!getFileName().empty())
    contentNotNullBase()->unloadArrays();
}

std::string MEDFileAnyTypeFieldMultiTS::simpleRepr() const
{
  std::ostringstream oss;
  contentNotNullBase()->simpleRepr(0,oss,-1);
  simpleReprGlobs(oss);
  return oss.str();
}

std::size_t MEDFileAnyTypeFieldMultiTS::getHeapMemorySizeWithoutChildren() const
{
  return MEDFileFieldGlobsReal::getHeapMemorySizeWithoutChildren();
}

std::vector<const BigMemoryObject *> MEDFileAnyTypeFieldMultiTS::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(MEDFileFieldGlobsReal::getDirectChildrenWithNull());
  ret.push_back((const MEDFileAnyTypeFieldMultiTSWithoutSDA *)_content);
  return ret;
}

/*!
 * This method returns as MEDFileAnyTypeFieldMultiTS new instances as number of components in \a this.
 * The returned instances are deep copy of \a this except that for globals that are share with those contained in \a this.
 * ** WARNING ** do no forget to rename the ouput instances to avoid to write n-times in the same MED file field !
 */
std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > MEDFileAnyTypeFieldMultiTS::splitComponents() const
{
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *content(_content);
  if(!content)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::splitComponents : no content in this ! Unable to split components !");
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > contentsSplit=content->splitComponents();
  std::size_t sz(contentsSplit.size());
  std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=shallowCpy();
      ret[i]->_content=contentsSplit[i];
    }
  return ret;
}

/*!
 * This method returns as MEDFileAnyTypeFieldMultiTS new instances as number of discretizations over time steps in \a this.
 * The returned instances are shallow copied of \a this included globals that are share with those contained in \a this.
 */
std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > MEDFileAnyTypeFieldMultiTS::splitDiscretizations() const
{
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *content(_content);
  if(!content)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::splitDiscretizations : no content in this ! Unable to split discretizations !");
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > contentsSplit(content->splitDiscretizations());
  std::size_t sz(contentsSplit.size());
  std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=shallowCpy();
      ret[i]->_content=contentsSplit[i];
    }
  return ret;
}

/*!
 * This method returns as MEDFileAnyTypeFieldMultiTS new instances as number of sub-discretizations over time steps in \a this.
 * The returned instances are shallow copied of \a this included globals that are share with those contained in \a this.
 */
std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > MEDFileAnyTypeFieldMultiTS::splitMultiDiscrPerGeoTypes() const
{
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *content(_content);
  if(!content)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::splitMultiDiscrPerGeoTypes : no content in this ! Unable to split discretizations !");
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > contentsSplit(content->splitMultiDiscrPerGeoTypes());
  std::size_t sz(contentsSplit.size());
  std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=shallowCpy();
      ret[i]->_content=contentsSplit[i];
    }
  return ret;
}

MEDFileAnyTypeFieldMultiTS *MEDFileAnyTypeFieldMultiTS::deepCopy() const
{
  MCAuto<MEDFileAnyTypeFieldMultiTS> ret=shallowCpy();
  if((const MEDFileAnyTypeFieldMultiTSWithoutSDA *)_content)
    ret->_content=_content->deepCopy();
  ret->deepCpyGlobs(*this);
  return ret.retn();
}

MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> MEDFileAnyTypeFieldMultiTS::getContent()
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
MEDFileAnyTypeField1TS *MEDFileAnyTypeFieldMultiTS::getTimeStep(int iteration, int order) const
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
MEDFileAnyTypeField1TS *MEDFileAnyTypeFieldMultiTS::getTimeStepGivenTime(double time, double eps) const
{
  int pos=getPosGivenTime(time,eps);
  return getTimeStepAtPos(pos);
}

/*!
 * This method groups not null items in \a vectFMTS per time step series. Two time series are considered equal if the list of their pair of integers iteration,order are equal.
 * The float64 value of time attached to the pair of integers are not considered here.
 * WARNING the returned pointers are not incremented. The caller is \b not responsible to deallocate them ! This method only reorganizes entries in \a vectFMTS.
 *
 * \param [in] vectFMTS - vector of not null fields defined on a same global data pointer.
 * \throw If there is a null pointer in \a vectFMTS.
 */
std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > MEDFileAnyTypeFieldMultiTS::SplitIntoCommonTimeSeries(const std::vector<MEDFileAnyTypeFieldMultiTS *>& vectFMTS)
{
  static const char msg[]="MEDFileAnyTypeFieldMultiTS::SplitIntoCommonTimeSeries : presence of null instance in input vector !";
  std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > ret;
  std::list<MEDFileAnyTypeFieldMultiTS *> lstFMTS(vectFMTS.begin(),vectFMTS.end());
  while(!lstFMTS.empty())
    {
      std::list<MEDFileAnyTypeFieldMultiTS *>::iterator it(lstFMTS.begin());
      MEDFileAnyTypeFieldMultiTS *curIt(*it);
      if(!curIt)
        throw INTERP_KERNEL::Exception(msg);
      std::vector< std::pair<int,int> > refIts=curIt->getIterations();
      std::vector<MEDFileAnyTypeFieldMultiTS *> elt;
      elt.push_back(curIt); it=lstFMTS.erase(it);
      while(it!=lstFMTS.end())
        {
          curIt=*it;
          if(!curIt)
            throw INTERP_KERNEL::Exception(msg);
          std::vector< std::pair<int,int> > curIts=curIt->getIterations();
          if(refIts==curIts)
            { elt.push_back(curIt); it=lstFMTS.erase(it); }
          else
            it++;
        }
      ret.push_back(elt);
    }
  return ret;
}

/*!
 * This method splits the input list \a vectFMTS considering the aspect of the geometrical support over time.
 * All returned instances in a subvector can be safely loaded, rendered along time
 * All items must be defined on the same time step ids ( see MEDFileAnyTypeFieldMultiTS::SplitIntoCommonTimeSeries method ).
 * Each item in \a vectFMTS is expected to have one and exactly one spatial discretization along time.
 * All items in \a vectFMTS must lie on the mesh (located by meshname and time step) and compatible with the input mesh \a mesh (having the same name than those in items).
 * All items in \a vectFMTS whose spatial discretization is not ON_NODES will appear once.
 * For items in \a vectFMTS that are ON_NODES it is possible to appear several times (more than once or once) in the returned vector.
 *
 * \param [in] vectFMTS - list of multi times step part all defined each on a same spatial discretization along time and pointing to a mesh whose name is equal to \c mesh->getName().
 * \param [in] mesh - the mesh shared by all items in \a vectFMTS across time.
 * \param [out] fsc - A vector having same size than returned vector. It specifies the support comporator of the corresponding vector of MEDFileAnyTypeFieldMultiTS in returned vector of vector.
 * \return - A vector of vector of objects that contains the same pointers (objects) than thoose in \a vectFMTS except that there are organized differently. So pointers included in returned vector of vector should \b not been dealt by the caller.
 *
 * \throw If an element in \a vectFMTS has not only one spatial discretization set.
 * \throw If an element in \a vectFMTS change of spatial discretization along time.
 * \throw If an element in \a vectFMTS lies on a mesh with meshname different from those in \a mesh.
 * \thorw If some elements in \a vectFMTS do not have the same times steps.
 * \throw If mesh is null.
 * \throw If an element in \a vectFMTS is null.
 * \sa MEDFileAnyTypeFieldMultiTS::AreOnSameSupportAcrossTime
 */
std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupport(const std::vector<MEDFileAnyTypeFieldMultiTS *>& vectFMTS, const MEDFileMesh *mesh, std::vector< MCAuto<MEDFileFastCellSupportComparator> >& fsc)
{
  static const char msg[]="MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupport : presence of a null instance in the input vector !";
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupport : input mesh is null !");
  std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > ret;
  if(vectFMTS.empty())
    return ret;
  std::vector<MEDFileAnyTypeFieldMultiTS *>::const_iterator it(vectFMTS.begin());
  MEDFileAnyTypeFieldMultiTS *frstElt(*it);
  if(!frstElt)
    throw INTERP_KERNEL::Exception(msg);
  std::size_t i=0;
  std::vector<MEDFileAnyTypeFieldMultiTS *> vectFMTSNotNodes;
  std::vector<MEDFileAnyTypeFieldMultiTS *> vectFMTSNodes;
  for(;it!=vectFMTS.end();it++,i++)
    {
      if(!(*it))
        throw INTERP_KERNEL::Exception(msg);
      TypeOfField tof0,tof1;
      if(CheckSupportAcrossTime(frstElt,*it,mesh,tof0,tof1)>0)
        {
          if(tof1!=ON_NODES)
            vectFMTSNotNodes.push_back(*it);
          else
            vectFMTSNodes.push_back(*it);
        }
      else
        vectFMTSNotNodes.push_back(*it);
    }
  std::vector< MCAuto<MEDFileFastCellSupportComparator> > cmps;
  std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > retCell=SplitPerCommonSupportNotNodesAlg(vectFMTSNotNodes,mesh,cmps);
  ret=retCell;
  for(std::vector<MEDFileAnyTypeFieldMultiTS *>::const_iterator it2=vectFMTSNodes.begin();it2!=vectFMTSNodes.end();it2++)
    {
      i=0;
      bool isFetched(false);
      for(std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> >::const_iterator it0=retCell.begin();it0!=retCell.end();it0++,i++)
        {
          if((*it0).empty())
            throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupport : internal error !");
          if(cmps[i]->isCompatibleWithNodesDiscr(*it2))
            { ret[i].push_back(*it2); isFetched=true; }
        }
      if(!isFetched)
        {
          std::vector<MEDFileAnyTypeFieldMultiTS *> tmp(1,*it2);
          MCAuto<MEDFileMeshStruct> tmp2(MEDFileMeshStruct::New(mesh));
          ret.push_back(tmp); retCell.push_back(tmp); cmps.push_back(MEDFileFastCellSupportComparator::New(tmp2,*it2));
        }
    }
  fsc=cmps;
  return ret;
}

/*!
 * WARNING no check here. The caller must be sure that all items in vectFMTS are coherent each other in time steps, only one same spatial discretization and not ON_NODES.
 * \param [out] cmps - same size than the returned vector.
 */
std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupportNotNodesAlg(const std::vector<MEDFileAnyTypeFieldMultiTS *>& vectFMTS, const MEDFileMesh *mesh, std::vector< MCAuto<MEDFileFastCellSupportComparator> >& cmps)
{
  std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > ret;
  std::list<MEDFileAnyTypeFieldMultiTS *> lstFMTS(vectFMTS.begin(),vectFMTS.end());
  while(!lstFMTS.empty())
    {
      std::list<MEDFileAnyTypeFieldMultiTS *>::iterator it(lstFMTS.begin());
      MEDFileAnyTypeFieldMultiTS *ref(*it);
      std::vector<MEDFileAnyTypeFieldMultiTS *> elt;
      elt.push_back(ref); it=lstFMTS.erase(it);
      MCAuto<MEDFileMeshStruct> mst(MEDFileMeshStruct::New(mesh));
      MCAuto<MEDFileFastCellSupportComparator> cmp(MEDFileFastCellSupportComparator::New(mst,ref));
      while(it!=lstFMTS.end())
        {
          MEDFileAnyTypeFieldMultiTS *curIt(*it);
          if(cmp->isEqual(curIt))
            { elt.push_back(curIt); it=lstFMTS.erase(it); }
          else
            it++;
        }
      ret.push_back(elt); cmps.push_back(cmp);
    }
  return ret;
}

/*!
 * This method scan the two main structs along time of \a f0 and \a f1 to see if there are all lying on the same mesh along time than those in \a mesh.
 * \a f0 and \a f1 must be defined each only on a same spatial discretization even if this can be different each other.
 *
 * \throw If \a f0 or \a f1 has not only one spatial discretization set.
 * \throw If \a f0 or \a f1 change of spatial discretization along time.
 * \throw If \a f0 or \a f1 on a mesh with meshname different from those in \a mesh.
 * \thorw If \a f0 and \a f1 do not have the same times steps.
 * \throw If mesh is null.
 * \throw If \a f0 or \a f1 is null.
 * \sa MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupport
 */
int MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime(MEDFileAnyTypeFieldMultiTS *f0, MEDFileAnyTypeFieldMultiTS *f1, const MEDFileMesh *mesh, TypeOfField& tof0, TypeOfField& tof1)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : input mesh is null !");
  if(!f0 || !f1)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : presence of null instance in fields over time !");
  if(f0->getMeshName()!=mesh->getName())
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : first field points to mesh \""<< f0->getMeshName() << "\" and input mesh to compare has name \"" << mesh->getName() << "\" !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  if(f1->getMeshName()!=mesh->getName())
    {
      std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : second field points to mesh \""<< f1->getMeshName() << "\" and input mesh to compare has name \"" << mesh->getName() << "\" !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  int nts=f0->getNumberOfTS();
  if(nts!=f1->getNumberOfTS())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : number of time steps are not the same !");
  if(nts==0)
    return nts;
  for(int i=0;i<nts;i++)
    {
      MCAuto<MEDFileAnyTypeField1TS> f0cur=f0->getTimeStepAtPos(i);
      MCAuto<MEDFileAnyTypeField1TS> f1cur=f1->getTimeStepAtPos(i);
      std::vector<TypeOfField> tofs0(f0cur->getTypesOfFieldAvailable()),tofs1(f1cur->getTypesOfFieldAvailable());
      if(tofs0.size()!=1 || tofs1.size()!=1)
        throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : All time steps must be defined on only one spatial discretization !");
      if(i!=0)
        {
          if(tof0!=tofs0[0] || tof1!=tofs1[0])
            throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : Across times steps MEDFileAnyTypeFieldMultiTS instances have to keep the same unique spatial discretization !");
        }
      else
        { tof0=tofs0[0]; tof1=tofs1[0]; }
      if(f0cur->getMeshIteration()!=mesh->getIteration() || f0cur->getMeshOrder()!=mesh->getOrder())
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : first field points to mesh time step (" << f0cur->getMeshIteration() << ","<< f0cur->getMeshOrder() << ") whereas input mesh points to time step (" << mesh->getIteration() << "," << mesh->getOrder() << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      if(f1cur->getMeshIteration()!=mesh->getIteration() || f1cur->getMeshOrder()!=mesh->getOrder())
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : second field points to mesh time step (" << f1cur->getMeshIteration() << ","<< f1cur->getMeshOrder() << ") whereas input mesh points to time step (" << mesh->getIteration() << "," << mesh->getOrder() << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      if(f0cur->getIteration()!=f1cur->getIteration() || f0cur->getOrder()!=f1cur->getOrder())
        {
          std::ostringstream oss; oss << "MEDFileAnyTypeFieldMultiTS::CheckSupportAcrossTime : all the time steps must be the same ! it is not the case (" << f0cur->getIteration() << "," << f0cur->getOrder() << ")!=(" << f1cur->getIteration() << "," << f1cur->getOrder() << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  return nts;
}

template<class T>
MCAuto<MEDFileAnyTypeField1TS> AggregateHelperF1TS(const std::vector< typename MLFieldTraits<T>::F1TSType const *>& f1tss, const std::vector< std::vector< std::pair<int,int> > >& dts)
{
  MCAuto< typename MLFieldTraits<T>::F1TSType > ret(MLFieldTraits<T>::F1TSType::New());
  if(f1tss.empty())
    throw INTERP_KERNEL::Exception("AggregateHelperF1TS : empty vector !");
  std::size_t sz(f1tss.size()),i(0);
  std::vector< typename MLFieldTraits<T>::F1TSWSDAType const *> f1tsw(sz);
  for(typename std::vector< typename MLFieldTraits<T>::F1TSType const *>::const_iterator it=f1tss.begin();it!=f1tss.end();it++,i++)
    {
      typename MLFieldTraits<T>::F1TSType const *elt(*it);
      if(!elt)
        throw INTERP_KERNEL::Exception("AggregateHelperF1TS : presence of a null pointer !");
      f1tsw[i]=dynamic_cast<typename MLFieldTraits<T>::F1TSWSDAType const *>(elt->contentNotNullBase());
    }
  typename MLFieldTraits<T>::F1TSWSDAType *retc(dynamic_cast<typename MLFieldTraits<T>::F1TSWSDAType *>(ret->contentNotNullBase()));
  if(!retc)
    throw INTERP_KERNEL::Exception("AggregateHelperF1TS : internal error 1 !");
  retc->aggregate(f1tsw,dts);
  ret->setDtUnit(f1tss[0]->getDtUnit());
  return DynamicCast<typename MLFieldTraits<T>::F1TSType , MEDFileAnyTypeField1TS>(ret);
}

template<class T>
MCAuto< MEDFileAnyTypeFieldMultiTS > AggregateHelperFMTS(const std::vector< typename MLFieldTraits<T>::FMTSType const *>& fmtss, const std::vector< std::vector< std::pair<int,int> > >& dts)
{
  MCAuto< typename MLFieldTraits<T>::FMTSType > ret(MLFieldTraits<T>::FMTSType::New());
  if(fmtss.empty())
    throw INTERP_KERNEL::Exception("AggregateHelperFMTS : empty vector !");
  std::size_t sz(fmtss.size());
  for(typename std::vector< typename MLFieldTraits<T>::FMTSType const *>::const_iterator it=fmtss.begin();it!=fmtss.end();it++)
    {
      typename MLFieldTraits<T>::FMTSType const *elt(*it);
      if(!elt)
        throw INTERP_KERNEL::Exception("AggregateHelperFMTS : presence of null pointer !");
    }
  int nbTS(fmtss[0]->getNumberOfTS());
  for(typename std::vector< typename MLFieldTraits<T>::FMTSType const *>::const_iterator it=fmtss.begin();it!=fmtss.end();it++)
    if((*it)->getNumberOfTS()!=nbTS)
      throw INTERP_KERNEL::Exception("AggregateHelperFMTS : all fields must have the same number of TS !");
  for(int iterTS=0;iterTS<nbTS;iterTS++)
    {
      std::size_t i(0);
      std::vector< typename MLFieldTraits<T>::F1TSType const *> f1tss(sz);
      std::vector< MCAuto<typename MLFieldTraits<T>::F1TSType> > f1tss2(sz);
      for(typename std::vector< typename MLFieldTraits<T>::FMTSType const *>::const_iterator it=fmtss.begin();it!=fmtss.end();it++,i++)
        { f1tss2[i]=(*it)->getTimeStepAtPos(iterTS); f1tss[i]=f1tss2[i]; }
      MCAuto<MEDFileAnyTypeField1TS> f1ts(AggregateHelperF1TS<T>(f1tss,dts));
      ret->pushBackTimeStep(f1ts);
      ret->setDtUnit(f1ts->getDtUnit());
    }
  return DynamicCast<typename MLFieldTraits<T>::FMTSType , MEDFileAnyTypeFieldMultiTS>(ret);
}

/*!
 * \a dts and \a ftmss are expected to have same size.
 */
MCAuto<MEDFileAnyTypeFieldMultiTS> MEDFileAnyTypeFieldMultiTS::Aggregate(const std::vector<const MEDFileAnyTypeFieldMultiTS *>& fmtss, const std::vector< std::vector< std::pair<int,int> > >& dts)
{
  if(fmtss.empty())
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::Aggregate : input vector is empty !");
  std::size_t sz(fmtss.size());
  std::vector<const MEDFileFieldMultiTS *> fmtss1;
  std::vector<const MEDFileIntFieldMultiTS *> fmtss2;
  for(std::vector<const MEDFileAnyTypeFieldMultiTS *>::const_iterator it=fmtss.begin();it!=fmtss.end();it++)
    {
      if(!(*it))
        throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::Aggregate : presence of null instance in input vector !");
      const MEDFileFieldMultiTS *elt1(dynamic_cast<const MEDFileFieldMultiTS *>(*it));
      if(elt1)
        {
          fmtss1.push_back(elt1);
          continue;
        }
      const MEDFileIntFieldMultiTS *elt2(dynamic_cast<const MEDFileIntFieldMultiTS *>(*it));
      if(elt2)
        {
          fmtss2.push_back(elt2);
          continue;
        }
      throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::Aggregate : not recognized type !");
    }
  if(fmtss1.size()!=sz && fmtss2.size()!=sz)
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::Aggregate : type of data is not homogeneous !");
  if(fmtss1.size()==sz)
    return AggregateHelperFMTS<double>(fmtss1,dts);
  if(fmtss2.size()!=sz)
    return AggregateHelperFMTS<int>(fmtss2,dts);
  throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::Aggregate : not implemented yet !");
}

MEDFileAnyTypeFieldMultiTSIterator *MEDFileAnyTypeFieldMultiTS::iterator()
{
  return new MEDFileAnyTypeFieldMultiTSIterator(this);
}

//= MEDFileFieldMultiTS

MEDFileAnyTypeFieldMultiTS *MEDFileFieldMultiTS::shallowCpy() const
{
  return new MEDFileFieldMultiTS(*this);
}

/*!
 * This method performs a copy with datatype modification ( float64->int32 ) of \a this. The globals information are copied
 * following the given input policy.
 *
 * \param [in] isDeepCpyGlobs - a boolean that indicates the behaviour concerning globals (profiles and localizations)
 *                            By default (true) the globals are deeply copied.
 * \return MEDFileIntFieldMultiTS * - a new object that is the result of the conversion of \a this to int32 field.
 */
MEDFileIntFieldMultiTS *MEDFileFieldMultiTS::convertToInt(bool isDeepCpyGlobs) const
{
  MCAuto<MEDFileIntFieldMultiTS> ret;
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *content(_content);
  if(content)
    {
      const MEDFileFieldMultiTSWithoutSDA *contc=dynamic_cast<const MEDFileFieldMultiTSWithoutSDA *>(content);
      if(!contc)
        throw INTERP_KERNEL::Exception("MEDFileFieldMultiTS::convertToInt : the content inside this is not FLOAT64 ! This is incoherent !");
      MCAuto<MEDFileIntFieldMultiTSWithoutSDA> newc(contc->convertToInt());
      ret=static_cast<MEDFileIntFieldMultiTS *>(MEDFileAnyTypeFieldMultiTS::BuildNewInstanceFromContent((MEDFileIntFieldMultiTSWithoutSDA *)newc));
    }
  else
    ret=MEDFileIntFieldMultiTS::New();
  if(isDeepCpyGlobs)
    ret->deepCpyGlobs(*this);
  else
    ret->shallowCpyGlobs(*this);
  return ret.retn();
}

MEDFileFieldMultiTS::MEDFileFieldMultiTS(med_idt fid, bool loadAll, const MEDFileMeshes *ms)
try:MEDFileTemplateFieldMultiTS<double>(fid,loadAll,ms)
{
}
catch(INTERP_KERNEL::Exception& e)
{ throw e; }

MEDFileFieldMultiTS::MEDFileFieldMultiTS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
try:MEDFileTemplateFieldMultiTS<double>(fid,fieldName,loadAll,ms,entities)
{
}
catch(INTERP_KERNEL::Exception& e)
{ throw e; }

MEDFileFieldMultiTS::MEDFileFieldMultiTS(const MEDFileFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent):MEDFileTemplateFieldMultiTS<double>(other,shallowCopyOfContent)
{
}

std::vector< std::vector<DataArrayDouble *> > MEDFileFieldMultiTS::getFieldSplitedByType2(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const
{
  return contentNotNull()->getFieldSplitedByType2(iteration,order,mname,types,typesF,pfls,locs);
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

MEDFileAnyTypeField1TS *MEDFileAnyTypeFieldMultiTSIterator::nextt()
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
 * This method is the symetrical method of MEDFileFields::loadArraysIfNecessary.
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
 * \return If true a renumbering has been performed. The structure in \a this has been modified. If false, nothing has been done: it is typically the case if \a meshName is not refered by any 
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
