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

#include "MEDFileField1TS.hxx"
#include "MEDFileFieldVisitor.hxx"
#include "MEDFileSafeCaller.txx"
#include "MEDLoaderBase.hxx"
#include "MEDFileField.txx"

#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingFieldDouble.hxx"

using namespace MEDCoupling;

extern INTERP_KERNEL::NormalizedCellType typmai2[MED_N_CELL_FIXED_GEO];

template class MEDCoupling::MEDFileField1TSTemplateWithoutSDA<int>;
template class MEDCoupling::MEDFileField1TSTemplateWithoutSDA<float>;
template class MEDCoupling::MEDFileField1TSTemplateWithoutSDA<double>;
template class MEDCoupling::MEDFileField1TSNDTemplateWithoutSDA<int>;
template class MEDCoupling::MEDFileField1TSNDTemplateWithoutSDA<float>;
template class MEDCoupling::MEDFileTemplateField1TS<int>;
template class MEDCoupling::MEDFileTemplateField1TS<float>;
template class MEDCoupling::MEDFileTemplateField1TS<double>;
template class MEDCoupling::MEDFileNDTemplateField1TS<int>;
template class MEDCoupling::MEDFileNDTemplateField1TS<float>;

const char MEDFileField1TSWithoutSDA::TYPE_STR[]="FLOAT64";
const char MEDFileIntField1TSWithoutSDA::TYPE_STR[]="INT32";
const char MEDFileFloatField1TSWithoutSDA::TYPE_STR[]="FLOAT32";

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
 * This method performs a sub splitting as splitDiscretizations does but finer. This is the finest splitting level that can be done.
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
    throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::writeLL : In MED3.0 mode in writing mode only ONE underlying mesh supported !");
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
 *  \param [in] nasc - the name scope used to assign names. Depends on the caller on the top call stack
 *  \param [in] smartPflKiller - specifies if this method tries at most to avoid profiles
 *  \throw If either \a field or \a mesh or \a profile has an empty name.
 *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in \a mesh.
 *  \throw If the data array of \a field is not set.
 *  \throw If \a this->_arr is already allocated but has different number of components
 *         than \a field.
 *  \throw If elements in \a mesh are not in the order suitable for writing to the MED file.
 *  \sa setFieldNoProfileSBT()
 */
void MEDFileAnyTypeField1TSWithoutSDA::setFieldProfile(const TimeHolder *th, const MEDCouplingFieldTemplate *field, const DataArray *arrOfVals, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc, bool smartPflKiller)
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
      m->splitProfilePerType(profile,code,idsInPflPerType,idsPerType,smartPflKiller);
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
  arr->setInfoOnComponents(infos);//will throw an exception if number of components mismatches
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
 * It is the responsibility of the caller to invoke MEDFileFieldGlobs::appendProfile or MEDFileFieldGlobs::getProfile
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
 * It is the responsibility of the caller to invoke MEDFileFieldGlobs::appendProfile or MEDFileFieldGlobs::getProfile
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
 * This method is the symmetrical method of MEDFileAnyTypeField1TS::loadArraysIfNecessary.
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
 * The returned instances are deep copy of \a this except that for globals that are shared with those contained in \a this.
 * ** WARNING ** do no forget to rename the output instances to avoid to write n-times in the same MED file field !
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
 * The returned instances are shallowed copied of \a this except that for globals that are shared with those contained in \a this.
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
 * The returned instances are shallowed copied of \a this except that for globals that are shared with those contained in \a this.
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
