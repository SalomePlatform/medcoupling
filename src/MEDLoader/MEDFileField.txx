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

#ifndef __MEDFILEFIELD_TXX__
#define __MEDFILEFIELD_TXX__

#include "MEDFileField.hxx"
#include "MEDCouplingTraits.hxx"

namespace MEDCoupling
{
  template<class T>
  void MEDFileField1TSTemplateWithoutSDA<T>::setArray(DataArray *arr)
  {
    if(!arr)
      {
        _nb_of_tuples_to_be_allocated=-1;
        _arr=0;
        return ;
      }
    typename Traits<T>::ArrayType *arrC=dynamic_cast<typename Traits<T>::ArrayType *>(arr);
    if(!arrC)
      throw INTERP_KERNEL::Exception("MEDFileField1TSTemplateWithoutSDA::setArray : the input not null array is not of type DataArrayDouble !");
    else
      _nb_of_tuples_to_be_allocated=-3;
    arrC->incrRef();
    _arr=arrC;
  }
  
  /*!
   * Returns a pointer to the underground DataArrayTemplate<T>  instance. So the
   * caller should not decrRef() it. This method allows for a direct access to the field
   * values. This method is quite unusable if there is more than a nodal field or a cell
   * field on single geometric cell type. 
   *  \return DataArrayTemplate<T> * - the pointer to the field values array.
   */
  template<class T>
  DataArray *MEDFileField1TSTemplateWithoutSDA<T>::getOrCreateAndGetArray()
  {
    return getOrCreateAndGetArrayTemplate();
  }
  
  template<class T>
  const DataArray *MEDFileField1TSTemplateWithoutSDA<T>::getOrCreateAndGetArray() const
  {
    return getOrCreateAndGetArrayTemplate();
  }
  
  template<class T>
  DataArray *MEDFileField1TSTemplateWithoutSDA<T>::createNewEmptyDataArrayInstance() const
  {
    return Traits<T>::ArrayType::New();
  }
  
  template<class T>
  typename Traits<T>::ArrayType const *MEDFileField1TSTemplateWithoutSDA<T>::getOrCreateAndGetArrayTemplate() const
  {
    typename Traits<T>::ArrayType const *ret(_arr);
    if(ret)
      return ret;
    (const_cast< MEDFileField1TSTemplateWithoutSDA<T> *>(this))->_arr=Traits<T>::ArrayType::New();
    return _arr;
  }
  
  template<class T>
  typename Traits<T>::ArrayType *MEDFileField1TSTemplateWithoutSDA<T>::getOrCreateAndGetArrayTemplate()
  {
    typename Traits<T>::ArrayType *ret(_arr);
    if(ret)
      return ret;
    _arr=Traits<T>::ArrayType::New();
    return _arr;
  }

  /*!
   * Returns a pointer to the underground DataArrayTemplate<T> instance. So the
   * caller should not decrRef() it. This method allows for a direct access to the field
   * values. This method is quite unusable if there is more than a nodal field or a cell
   * field on single geometric cell type. 
   *  \return DataArrayTemplate<T> * - the pointer to the field values array.
   */
  template<class T>
  typename Traits<T>::ArrayType *MEDFileField1TSTemplateWithoutSDA<T>::getUndergroundDataArrayTemplate() const
  {
    typename Traits<T>::ArrayType const *ret(_arr);
    if(ret)
      return const_cast<typename Traits<T>::ArrayType *>(ret);
    else
      return 0;
  }
  
  /*!
   * Returns a pointer to the underground DataArrayDouble instance. So the
   * caller should not decrRef() it. This method allows for a direct access to the field
   * values. This method is quite unusable if there is more than a nodal field or a cell
   * field on single geometric cell type. 
   *  \return DataArrayDouble * - the pointer to the field values array.
   */
  template<class T>
  DataArray *MEDFileField1TSTemplateWithoutSDA<T>::getUndergroundDataArray() const
  {
    return getUndergroundDataArrayTemplate();
  }
  
  template<class T>
  void MEDFileField1TSTemplateWithoutSDA<T>::aggregate(const std::vector<typename MLFieldTraits<T>::F1TSWSDAType const *>& f1tss, const std::vector< std::vector< std::pair<int,int> > >& dts)
  {
    if(f1tss.empty())
      throw INTERP_KERNEL::Exception("MEDFileField1TSTemplateWithoutSDA::aggregate : empty vector !");
    std::size_t sz(f1tss.size()),ii(0);
    std::vector<const MEDFileFieldPerMesh *> pms;
    std::vector<const DataArray *> das(sz);
    for(typename std::vector<typename MLFieldTraits<T>::F1TSWSDAType const *>::const_iterator it=f1tss.begin();it!=f1tss.end();it++,ii++)
      {
        if(!*it)
          throw INTERP_KERNEL::Exception("MEDFileField1TSTemplateWithoutSDA::aggregate : presence of null pointer in input vector !");
        if((*it)->_field_per_mesh.empty())
          throw INTERP_KERNEL::Exception("MEDFileField1TSTemplateWithoutSDA::aggregate : no info !");
        const typename Traits<T>::ArrayType *arr((*it)->getUndergroundDataArrayTemplate());
        if(!arr)
          throw INTERP_KERNEL::Exception("MEDFileField1TSTemplateWithoutSDA::aggregate : presence of null array !");
        das[ii]=arr;
        pms.push_back((*it)->_field_per_mesh[0]);
      }
    typename MLFieldTraits<T>::F1TSWSDAType const *refPt(f1tss[0]);
    setName(refPt->getName());
    
    const DataArray *arr(refPt->getUndergroundDataArray());
    int nbCompo(arr->getNumberOfComponents());
    for(typename std::vector<typename MLFieldTraits<T>::F1TSWSDAType const *>::const_iterator it=f1tss.begin();it!=f1tss.end();it++)
      {
        const typename Traits<T>::ArrayType *myArr((*it)->getUndergroundDataArrayTemplate());
        if(myArr->getNumberOfComponents()!=nbCompo)
          throw INTERP_KERNEL::Exception("MEDFileField1TSTemplateWithoutSDA::aggregate : arrays must have same number of components !");
      }
    std::vector<std::pair< int, std::pair<int,int> > > extractInfo;
    int start(0);
    MCAuto<MEDFileFieldPerMesh> fpm(MEDFileFieldPerMesh::Aggregate(start,pms,dts,this,extractInfo));
    _field_per_mesh.push_back(fpm);
    int iteration,order;
    double tv(f1tss[0]->getTime(iteration,order));
    _iteration=iteration; _order=order; _dt=tv;
    _arr=Traits<T>::ArrayType::New();
    _arr->alloc(start,nbCompo); _arr->copyStringInfoFrom(*arr);
    start=0;
    for(std::vector<std::pair< int, std::pair<int,int> > >::const_iterator it=extractInfo.begin();it!=extractInfo.end();it++)
      {
        const DataArray *zeArr(das[(*it).first]);
        _arr->setContigPartOfSelectedValuesSlice(start,zeArr,(*it).second.first,(*it).second.second,1);
        start+=(*it).second.second-(*it).second.first;
      }
  }
}

#endif
