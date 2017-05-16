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
   *  \sa getUndergroundDataArrayTemplate()
   */
  template<class T>
  typename Traits<T>::ArrayType *MEDFileField1TSTemplateWithoutSDA<T>::getUndergroundDataArrayTemplateExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const
  {
    if(this->_field_per_mesh.size()!=1)
      throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt : field lies on several meshes, this method has no sense !");
    if(this->_field_per_mesh[0]==0)
      throw INTERP_KERNEL::Exception("MEDFileField1TSWithoutSDA::getUndergroundDataArrayExt : no field specified !");
    this->_field_per_mesh[0]->getUndergroundDataArrayExt(entries);
    return getUndergroundDataArrayTemplate();
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

  template<class T>
  MEDFileTemplateField1TS<T>::MEDFileTemplateField1TS()
  {
    _content=new typename MLFieldTraits<T>::F1TSWSDAType;
  }

  /*!
   * Returns a new empty instance of MEDFileField1TS.
   *  \return MEDFileField1TS * - a new instance of MEDFileField1TS. The caller
   *          is to delete this field using decrRef() as it is no more needed.
   */
  template<class T>
  typename MLFieldTraits<T>::F1TSType *MEDFileTemplateField1TS<T>::New()
  {
    MCAuto<typename MLFieldTraits<T>::F1TSType> ret(new typename MLFieldTraits<T>::F1TSType);
    ret->contentNotNull();
    return ret.retn();
  }

  /*!
   * Returns a new instance of MEDFileField1TS holding data of the first time step of 
   * the first field that has been read from a specified MED file.
   *  \param [in] fileName - the name of the MED file to read.
   *  \return MEDFileField1TS * - a new instance of MEDFileFieldMultiTS. The caller
   *          is to delete this field using decrRef() as it is no more needed.
   *  \throw If reading the file fails.
   */
  template<class T>
  typename MLFieldTraits<T>::F1TSType *MEDFileTemplateField1TS<T>::New(const std::string& fileName, bool loadAll)
  {
    MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
    return New(fid,loadAll);
  }
  
  template<class T>
  typename MLFieldTraits<T>::F1TSType *MEDFileTemplateField1TS<T>::New(med_idt fid, bool loadAll)
  {
    MCAuto<typename MLFieldTraits<T>::F1TSType> ret(new typename MLFieldTraits<T>::F1TSType(fid,loadAll,0));
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
  template<class T>
  typename MLFieldTraits<T>::F1TSType *MEDFileTemplateField1TS<T>::New(const std::string& fileName, const std::string& fieldName, bool loadAll)
  {
    MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
    return New(fid,fieldName,loadAll);
  }

  template<class T>
  typename MLFieldTraits<T>::F1TSType *MEDFileTemplateField1TS<T>::New(med_idt fid, const std::string& fieldName, bool loadAll)
  {
    MCAuto<typename MLFieldTraits<T>::F1TSType> ret(new typename MLFieldTraits<T>::F1TSType(fid,fieldName,loadAll,0));
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
  template<class T>
  typename MLFieldTraits<T>::F1TSType *MEDFileTemplateField1TS<T>::New(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll)
  {
    MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
    return New(fid,fieldName,iteration,order,loadAll);
  }
  
  template<class T>
  typename MLFieldTraits<T>::F1TSType *MEDFileTemplateField1TS<T>::New(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll)
  {
    MCAuto<typename MLFieldTraits<T>::F1TSType> ret(new typename MLFieldTraits<T>::F1TSType(fid,fieldName,iteration,order,loadAll,0));
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
  template<class T>
  typename MLFieldTraits<T>::F1TSType *MEDFileTemplateField1TS<T>::New(const typename MLFieldTraits<T>::F1TSWSDAType& other, bool shallowCopyOfContent)
  {
    MCAuto<typename MLFieldTraits<T>::F1TSType> ret(new typename MLFieldTraits<T>::F1TSType(other,shallowCopyOfContent));
    ret->contentNotNull();
    return ret.retn();
  }
  
  template<class T>
  const typename MLFieldTraits<T>::F1TSWSDAType *MEDFileTemplateField1TS<T>::contentNotNull() const
  {
    const MEDFileAnyTypeField1TSWithoutSDA *pt(_content);
    if(!pt)
      throw INTERP_KERNEL::Exception("MEDFileTemplateField1TS<T>::contentNotNull : the content pointer is null !");
    const typename MLFieldTraits<T>::F1TSWSDAType *ret(dynamic_cast<const typename MLFieldTraits<T>::F1TSWSDAType *>(pt));
    if(!ret)
      {
        std::ostringstream oss; oss << "MEDFileTemplateField1TS<T>::contentNotNull : the content pointer is not null but it is not of type double ! Reason is maybe that the read field has not the type " << MLFieldTraits<T>::F1TSWSDAType::TYPE_STR;
        throw INTERP_KERNEL::Exception(oss.str());
      }
    return ret;
  }
  
  template<class T>
  typename MLFieldTraits<T>::F1TSWSDAType *MEDFileTemplateField1TS<T>::contentNotNull()
  {
    MEDFileAnyTypeField1TSWithoutSDA *pt(_content);
    if(!pt)
      throw INTERP_KERNEL::Exception("MEDFileTemplateField1TS<T>::contentNotNull : the non const content pointer is null !");
    typename MLFieldTraits<T>::F1TSWSDAType *ret(dynamic_cast<typename MLFieldTraits<T>::F1TSWSDAType *>(pt));
    if(!ret)
      {
        std::ostringstream oss; oss << "MEDFileTemplateField1TS<T>::contentNotNull : the non const content pointer is not null but it is not of type double ! Reason is maybe that the read field has not the type " << MLFieldTraits<T>::F1TSWSDAType::TYPE_STR;
        throw INTERP_KERNEL::Exception(oss.str());
      }
    return ret;
  }
  
  template<class T>
  typename Traits<T>::ArrayType *MEDFileTemplateField1TS<T>::ReturnSafelyTypedDataArray(MCAuto<DataArray>& arr)
  {
    if(arr.isNull())
      throw INTERP_KERNEL::Exception("MEDFileField1TS::ReturnSafelyTypedDataArray : no array !");
    typename Traits<T>::ArrayType *arrOutC(dynamic_cast<typename Traits<T>::ArrayType *>((DataArray*)arr));
    if(!arrOutC)
      throw INTERP_KERNEL::Exception("MEDFileField1TS::ReturnSafelyTypedDataArray : mismatch between dataArrays type and MEDFileField1TS ! Expected double !");
    arrOutC->incrRef();
    return arrOutC;
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
  template<class T>
  typename Traits<T>::ArrayType *MEDFileTemplateField1TS<T>::getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const
  {
    MCAuto<DataArray> arr(contentNotNull()->getFieldWithProfile(type,meshDimRelToMax,mesh,pfl,this,*contentNotNull()));
    return ReturnSafelyTypedDataArray(arr);
  }

  template<class T>
  typename Traits<T>::ArrayType *MEDFileTemplateField1TS<T>::getUndergroundDataArray() const
  {
    return contentNotNull()->getUndergroundDataArrayTemplate();
  }
  
  template<class T>
  typename Traits<T>::ArrayType *MEDFileTemplateField1TS<T>::getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const
  {
    return contentNotNull()->getUndergroundDataArrayTemplateExt(entries);
  }
}

#endif
