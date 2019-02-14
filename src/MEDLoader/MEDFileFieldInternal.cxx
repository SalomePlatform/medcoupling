// Copyright (C) 2017-2019  CEA/DEN, EDF R&D
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

#include "MEDFileFieldInternal.hxx"
#include "MEDFileField.hxx"
#include "MEDFileFieldVisitor.hxx"
#include "MEDFileStructureElement.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDFileSafeCaller.txx"
#include "MEDFileEntities.hxx"

#include "MEDCouplingGaussLocalization.hxx"
#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include "CellModel.hxx"

extern med_geometry_type typmai[MED_N_CELL_FIXED_GEO];
extern INTERP_KERNEL::NormalizedCellType typmai2[MED_N_CELL_FIXED_GEO];
extern med_geometry_type typmai3[34];

using namespace MEDCoupling;

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
  MEDFILESAFECALLERRD0(MEDlocalizationRd,(fid,locName.c_str(),MED_FULL_INTERLACE,&_ref_coo[0],&_gs_coo[0],&_w[0]));
}

MEDFileFieldLoc::MEDFileFieldLoc(med_idt fid, int id, const MEDFileEntities *entities)
{
  med_geometry_type geotype;
  med_geometry_type sectiongeotype;
  int nsectionmeshcell;
  INTERP_KERNEL::AutoPtr<char> locName(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> geointerpname(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  INTERP_KERNEL::AutoPtr<char> sectionmeshname(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
  MEDFILESAFECALLERRD0(MEDlocalizationInfo,(fid,id+1,locName,&geotype,&_dim,&_nb_gauss_pt,geointerpname,sectionmeshname,&nsectionmeshcell,&sectiongeotype));
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
        MEDFILESAFECALLERRD0(MEDmeshGeotypeParameter,(fid,geotype,&dummy,&_nb_node_per_cell));
      }
    }
  _ref_coo.resize(_dim*_nb_node_per_cell);
  _gs_coo.resize(_dim*_nb_gauss_pt);
  _w.resize(_nb_gauss_pt);
  MEDFILESAFECALLERRD0(MEDlocalizationRd,(fid,locName,MED_FULL_INTERLACE,&_ref_coo[0],&_gs_coo[0],&_w[0]));
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
  MEDFILESAFECALLERWR0(MEDlocalizationWr,(fid,_name.c_str(),typmai3[(int)getGeoType()],_dim,&_ref_coo[0],MED_FULL_INTERLACE,_nb_gauss_pt,&_gs_coo[0],&_w[0],MED_NO_INTERPOLATION,MED_NO_MESH_SUPPORT));
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
 * Leaf method of field with profile assignment. This method is the most general one. No optimization is done here.
 * \param [in] pflName input containing name of profile if any. 0 if no profile (except for GAUSS_PT where a no profile can hide a profile when split by loc_id).
 * \param [in] multiTypePfl is the end user profile specified in high level API
 * \param [in] idsInPfl is the selection into the \a multiTypePfl whole profile that corresponds to the current geometric type.
 * \param [in] locIds is the profile needed to be created for MED file format. It can be null if all cells of current geometric type are fetched in \a multiTypePfl.
 *             \b WARNING if not null the MED file profile can be subdivided again in case of Gauss points.
 * \param [in] mesh is the mesh coming from the MEDFileMesh instance in correspondence with the MEDFileField. The mesh inside the \a field is simply ignored.
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
          MEDFILESAFECALLERRD0(MEDfilterBlockOfEntityCr,(fid,/*nentity*/overallNval,/*nvaluesperentity*/nbi,/*nconstituentpervalue*/nbOfCompo,
                                                         MED_ALL_CONSTITUENT,MED_FULL_INTERLACE,MED_COMPACT_STMODE,MED_NO_PROFILE,
                                                         /*start*/start+1,/*stride*/step,/*count*/1,/*blocksize*/nbOfEltsToLoad,
                                                         /*lastblocksize=useless because count=1*/0,&filter));
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
            MEDFILESAFECALLERRD0(MEDfilterBlockOfEntityCr,(fid,/*nentity*/overallNval,/*nvaluesperentity*/nbi,/*nconstituentpervalue*/nbOfCompo,
                                                           MED_ALL_CONSTITUENT,MED_FULL_INTERLACE,MED_COMPACT_STMODE,MED_NO_PROFILE,
                                                           /*start*/a+1,/*stride*/1,/*count*/1,/*blocksize*/nbOfEltsToLoad,
                                                           /*lastblocksize=useless because count=1*/0,&filter));
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
      if(_localization!="MED_GAUSS_ELNO")//For compatibility with MED2.3
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

void MEDFileFieldPerMeshPerTypePerDisc::incrementNbOfVals(int deltaNbVal)
{
  int nbi((_end-_start)/_nval);
  _nval+=deltaNbVal;
  _end+=nbi*deltaNbVal;
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
 * \param [in] offset the offset id used to take into account that \a result is not compulsory empty in input
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
 * \param [in] mesh is the mesh coming from the MEDFileMesh instance in correspondence with the MEDFileField. The mesh inside the \a field is simply ignored.
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

int MEDFileFieldPerMeshPerTypeCommon::locIdOfLeaf(const MEDFileFieldPerMeshPerTypePerDisc *leaf) const
{
  int ret(0);
  for(std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >::const_iterator it=_field_pm_pt_pd.begin();it!=_field_pm_pt_pd.end();it++,ret++)
    {
      const MEDFileFieldPerMeshPerTypePerDisc *cand(*it);
      if(cand==leaf)
        return ret;
    }
  throw INTERP_KERNEL::Exception("MEDFileFieldPerMeshPerTypeCommon::locIdOfLeaf : not found such a leaf in this !");
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
 * \param [in] mesh is the mesh coming from the MEDFileMesh instance in correspondence with the MEDFileField. The mesh inside the \a field is simply ignored.
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
  oss << "Possibilities are : ";
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
  oss << "Possibilities are : ";
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
 * No check of this is performed. 'da' array contains an array in old2New style to be applied to mesh to obtain the right support.
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
      oss << "So it is impossible to return a well defined MEDCouplingFieldDouble instance on specified mesh on a specified meshDim !" << std::endl;
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
