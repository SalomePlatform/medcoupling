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
#include "MEDCouplingFieldInt.hxx"
#include "MEDCouplingFieldFloat.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldTemplate.hxx"

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

  ///////////////////////////////////////////////////////

  template<class T>
  MEDFileField1TSWithoutSDA *MEDFileField1TSNDTemplateWithoutSDA<T>::convertToDouble() const
  {
    MCAuto<MEDFileField1TSWithoutSDA> ret(new MEDFileField1TSWithoutSDA);
    ret->MEDFileAnyTypeField1TSWithoutSDA::operator =(*this);
    ret->deepCpyLeavesFrom(*this);
    if(this->_arr.isNotNull())
      {
        MCAuto<DataArrayDouble> arr2(this->_arr->convertToDblArr());
        ret->setArray(arr2);
      }
    return ret.retn();
  }
  
  ///////////////////////////////////////////////////////

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

  template<class T>
  MCAuto<typename Traits<T>::FieldType> MEDFileTemplateField1TS<T>::SetDataArrayInField(MEDCouplingFieldDouble *f, MCAuto<DataArray>& arr)
  {
    if(!f)
      throw INTERP_KERNEL::Exception("MEDFileTemplateField1TS<T>::SetDataArrayInField : input field is NULL !");
    if(arr.isNull())
      throw INTERP_KERNEL::Exception("MEDFileTemplateField1TS<T>::SetDataArrayInField : no array !");
    int t1,t2;
    double t0(f->getTime(t1,t2));
    std::string tu(f->getTimeUnit());
    MCAuto<typename Traits<T>::ArrayType> arr2(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(arr));
    MCAuto<MEDCouplingFieldTemplate> ft(MEDCouplingFieldTemplate::New(*f));
    MCAuto<typename Traits<T>::FieldType> ret(Traits<T>::FieldType::New(*ft));
    ret->setTime(t0,t1,t2); ret->setArray(arr2); ret->setTimeUnit(tu);
    return ret.retn();
  }

  template<class T>
  MCAuto<MEDCouplingFieldDouble> MEDFileTemplateField1TS<T>::ToFieldTemplateWithTime(const typename Traits<T>::FieldType *f)
  {
    int t1,t2;
    double t0(f->getTime(t1,t2));
    std::string tu(f->getTimeUnit());
    MCAuto<MEDCouplingFieldTemplate> ft(MEDCouplingFieldTemplate::NewWithoutCheck(*f));
    MCAuto<MEDCouplingFieldDouble> ret(MEDCouplingFieldDouble::New(*ft));
    ret->setTime(t0,t1,t2); ret->setTimeUnit(tu);
    return ret.retn();
  }

  /*!
   * This is the simplest version to fetch a field for MED structure. One drawback : if \a this is a complex field (multi spatial discretization inside a same field) this method will throw exception and more advance
   * method should be called (getFieldOnMeshAtLevel for example).
   * But for normal usage of field in MED file world this method is the most efficient to fetch data.
   *
   * \param [in] mesh - the mesh the field is lying on
   * \return typename Traits<T>::FieldType * - a new instance of typename Traits<T>::FieldType. The
   *          caller is to delete this field using decrRef() as it is no more needed. 
   */
  template<class T>
  typename Traits<T>::FieldType *MEDFileTemplateField1TS<T>::field(const MEDFileMesh *mesh) const
  {
    MCAuto<DataArray> arrOut;
    MCAuto<MEDCouplingFieldDouble> ret(contentNotNull()->fieldOnMesh(this,mesh,arrOut,*contentNotNull()));
    MCAuto<typename Traits<T>::FieldType> ret2(MEDFileTemplateField1TS<T>::SetDataArrayInField(ret,arrOut));
    return ret2.retn();
  }

  /*!
   * Returns a new typename Traits<T>::FieldType of a given type lying on
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
   *  \return typename Traits<T>::FieldType * - a new instance of typename Traits<T>::FieldType. The
   *          caller is to delete this field using decrRef() as it is no more needed. 
   *  \throw If \a this field has not been constructed via file reading.
   *  \throw If the MED file is not readable.
   *  \throw If there is no mesh in the MED file.
   *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in the mesh.
   *  \throw If no field values of the given \a type or given \a meshDimRelToMax are available.
   *  \sa getFieldOnMeshAtLevel()
   */
  template<class T>
  typename Traits<T>::FieldType *MEDFileTemplateField1TS<T>::getFieldAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol) const
  {
    if(getFileName().empty())
      throw INTERP_KERNEL::Exception("MEDFileTemplateField1TS<T>::getFieldAtLevel : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtLevel method instead !");
    MCAuto<DataArray> arrOut;
    MCAuto<MEDCouplingFieldDouble> ret(contentNotNull()->getFieldAtLevel(type,meshDimRelToMax,std::string(),renumPol,this,arrOut,*contentNotNull()));
    MCAuto<typename Traits<T>::FieldType> ret2(MEDFileTemplateField1TS<T>::SetDataArrayInField(ret,arrOut));
    return ret2.retn();
  }

  /*!
   * Returns a new typename Traits<T>::FieldType of a given type lying on
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
   *  \return typename Traits<T>::FieldType * - a new instance of typename Traits<T>::FieldType. The
   *          caller is to delete this field using decrRef() as it is no more needed. 
   *  \throw If \a this field has not been constructed via file reading.
   *  \throw If the MED file is not readable.
   *  \throw If there is no mesh in the MED file.
   *  \throw If no field values of the given \a type.
   *  \throw If no field values lying on the top level support.
   *  \sa getFieldAtLevel()
   */
  template<class T>
  typename Traits<T>::FieldType *MEDFileTemplateField1TS<T>::getFieldAtTopLevel(TypeOfField type, int renumPol) const
  {
    if(getFileName().empty())
      throw INTERP_KERNEL::Exception("MEDFileTemplateField1TS<T>::getFieldAtTopLevel : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtTopLevel method instead !");
    MCAuto<DataArray> arrOut;
    MCAuto<MEDCouplingFieldDouble> ret(contentNotNull()->getFieldAtTopLevel(type,std::string(),renumPol,this,arrOut,*contentNotNull()));
    MCAuto<typename Traits<T>::FieldType> ret2(MEDFileTemplateField1TS<T>::SetDataArrayInField(ret,arrOut));
    return ret2.retn();
  }

  /*!
   * Returns a new typename Traits<T>::FieldType of given type lying on a given mesh.
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
   *  \return typename Traits<T>::FieldType * - a new instance of typename Traits<T>::FieldType. The
   *          caller is to delete this field using decrRef() as it is no more needed. 
   *  \throw If no field of \a this is lying on \a mesh.
   *  \throw If the mesh is empty.
   *  \throw If no field values of the given \a type are available.
   *  \sa getFieldAtLevel()
   *  \sa getFieldOnMeshAtLevel() 
   */
  template<class T>
  typename Traits<T>::FieldType *MEDFileTemplateField1TS<T>::getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, int renumPol) const
  {
    MCAuto<DataArray> arrOut;
    MCAuto<MEDCouplingFieldDouble> ret(contentNotNull()->getFieldOnMeshAtLevel(type,renumPol,this,mesh,0,0,arrOut,*contentNotNull()));
    MCAuto<typename Traits<T>::FieldType> ret2(MEDFileTemplateField1TS<T>::SetDataArrayInField(ret,arrOut));
    return ret2.retn();
  }

  /*!
   * Returns a new typename Traits<T>::FieldType of a given type lying on a given support.
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
   *  \return typename Traits<T>::FieldType * - a new instance of typename Traits<T>::FieldType. The
   *          caller is to delete this field using decrRef() as it is no more needed. 
   *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in the mesh.
   *  \throw If no field of \a this is lying on \a mesh.
   *  \throw If no field values of the given \a type or given \a meshDimRelToMax are available.
   *  \sa getFieldAtLevel()
   *  \sa getFieldOnMeshAtLevel() 
   */
  template<class T>
  typename Traits<T>::FieldType *MEDFileTemplateField1TS<T>::getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol) const
  {
    MCAuto<DataArray> arrOut;
    MCAuto<MEDCouplingFieldDouble> ret(contentNotNull()->getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,this,mesh,arrOut,*contentNotNull()));
    MCAuto<typename Traits<T>::FieldType> ret2(MEDFileTemplateField1TS<T>::SetDataArrayInField(ret,arrOut));
    return ret2.retn();
  }

  /*!
   * Returns a new typename Traits<T>::FieldType of a given type lying on a given support.
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
   *  \return typename Traits<T>::FieldType * - a new instance of typename Traits<T>::FieldType. The
   *          caller is to delete this field using decrRef() as it is no more needed. 
   *  \throw If the MED file is not readable.
   *  \throw If there is no mesh named \a mName in the MED file.
   *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in the mesh.
   *  \throw If \a this field has not been constructed via file reading.
   *  \throw If no field of \a this is lying on the mesh named \a mName.
   *  \throw If no field values of the given \a type or given \a meshDimRelToMax are available.
   *  \sa getFieldAtLevel()
   */
  template<class T>
  typename Traits<T>::FieldType *MEDFileTemplateField1TS<T>::getFieldAtLevelOld(TypeOfField type, const std::string& mname, int meshDimRelToMax, int renumPol) const
  {
    if(getFileName().empty())
      throw INTERP_KERNEL::Exception("MEDFileTemplateField1TS<T>::getFieldAtLevelOld : Request for a method that can be used for instances coming from file loading ! Use getFieldOnMeshAtLevel method instead !");
    MCAuto<DataArray> arrOut;
    MCAuto<MEDCouplingFieldDouble> ret(contentNotNull()->getFieldAtLevel(type,meshDimRelToMax,mname,renumPol,this,arrOut,*contentNotNull()));
    MCAuto<typename Traits<T>::FieldType> ret2(MEDFileTemplateField1TS<T>::SetDataArrayInField(ret,arrOut));
    return ret2.retn();
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
  template<class T>
  void MEDFileTemplateField1TS<T>::setFieldNoProfileSBT(const typename Traits<T>::FieldType *field)
  {
    setFileName("");
    MCAuto<MEDCouplingFieldTemplate> ft(MEDCouplingFieldTemplate::New(*field));
    contentNotNull()->setFieldNoProfileSBT(field->timeDiscrSafe(),ft,field->getArray(),*this,*contentNotNull());
  }

  /*!
   * Adds a MEDCouplingFieldDouble to \a this. As described in \ref MEDLoaderMainC a field in MED file sense
   * can be an aggregation of several MEDCouplingFieldDouble instances.
   * The mesh support of input parameter \a field is ignored here, it can be NULL.
   * The support of field \a field is expected to be those computed with the input parameter \a mesh, \a meshDimRelToMax,
   * and \a profile.
   *
   * This method will check that the field based on the computed support is coherent. If not an exception will be thrown.
   * A new profile is added only if no equal profile is missing.
   * For more info, see \ref AdvMEDLoaderAPIFieldRW
   *  \param [in] field - the field to add to \a this. The mesh support of field is ignored.
   *  \param [in] mesh - the supporting mesh of \a field.
   *  \param [in] meshDimRelToMax - a relative dimension of mesh entities \a field lies on (useless if field spatial discretization is ON_NODES).
   *  \param [in] profile - ids of mesh entities on which corresponding field values lie.
   *  \throw If either \a field or \a mesh or \a profile has an empty name.
   *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in \a mesh.
   *  \throw If the data array of \a field is not set.
   *  \throw If the data array of \a this is already allocated but has different number of
   *         components than \a field.
   *  \throw If elements in \a mesh are not in the order suitable for writing to the MED file.
   *  \sa setFieldNoProfileSBT()
   */
  template<class T>
  void MEDFileTemplateField1TS<T>::setFieldProfile(const typename Traits<T>::FieldType *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile)
  {
    setFileName("");
    MCAuto<MEDCouplingFieldTemplate> ft(MEDCouplingFieldTemplate::NewWithoutCheck(*field));
    contentNotNull()->setFieldProfile(field->timeDiscrSafe(),ft,field->getArray(),mesh,meshDimRelToMax,profile,*this,*contentNotNull());
  }

  /*!
   * Return an extraction of \a this using \a extractDef map to specify the extraction.
   * The keys of \a extractDef is level relative to max ext of \a mm mesh.
   *
   * \return A new object that the caller is responsible to deallocate.
   * \sa MEDFileUMesh::deduceNodeSubPartFromCellSubPart , MEDFileUMesh::extractPart
   */
  template<class T>
  typename MLFieldTraits<T>::F1TSType *MEDFileTemplateField1TS<T>::extractPartImpl(const std::map<int, MCAuto<DataArrayInt> >& extractDef, MEDFileMesh *mm) const
  {
    if(!mm)
      throw INTERP_KERNEL::Exception("MEDFileField1TS::extractPart : input mesh is NULL !");
    MCAuto<typename MLFieldTraits<T>::F1TSType> ret(MLFieldTraits<T>::F1TSType::New());
    std::vector<TypeOfField> tof(getTypesOfFieldAvailable());
    for(std::vector<TypeOfField>::const_iterator it0=tof.begin();it0!=tof.end();it0++)
      {
        if((*it0)!=ON_NODES)
          {
            std::vector<int> levs;
            getNonEmptyLevels(mm->getName(),levs);
            for(std::vector<int>::const_iterator lev=levs.begin();lev!=levs.end();lev++)
              {
                std::map<int, MCAuto<DataArrayInt> >::const_iterator it2(extractDef.find(*lev));
                if(it2!=extractDef.end())
                  {
                    MCAuto<DataArrayInt> t((*it2).second);
                    if(t.isNull())
                      throw INTERP_KERNEL::Exception("MEDFileField1TS::extractPart : presence of a value with null pointer 1 !");
                    MCAuto<typename Traits<T>::FieldType> f(getFieldOnMeshAtLevel(ON_CELLS,(*lev),mm));
                    MCAuto<typename Traits<T>::FieldType> fOut(f->buildSubPart(t));
                    ret->setFieldNoProfileSBT(fOut);
                  }
              }
          }
        else
          {
            std::map<int, MCAuto<DataArrayInt> >::const_iterator it2(extractDef.find(1));
            if(it2==extractDef.end())
              throw INTERP_KERNEL::Exception("MEDFileField1TS::extractPart : presence of a NODE field and no extract array available for NODE !");
            MCAuto<DataArrayInt> t((*it2).second);
            if(t.isNull())
              throw INTERP_KERNEL::Exception("MEDFileField1TS::extractPart : presence of a value with null pointer 1 !");
            MCAuto<typename Traits<T>::FieldType> f(getFieldOnMeshAtLevel(ON_NODES,0,mm));
            MCAuto<typename Traits<T>::FieldType> fOut(f->deepCopy());
            typename Traits<T>::ArrayType *arr(f->getArray());
            MCAuto<typename Traits<T>::ArrayType> newArr(arr->selectByTupleIdSafe(t->begin(),t->end()));
            fOut->setArray(newArr);
            ret->setFieldNoProfileSBT(fOut);
          }
      }
    return ret.retn();
  }

  //////////////////////////

  /*!
   * This method performs a copy with datatype modification ( int32->float64 ) of \a this. The globals information are copied
   * following the given input policy.
   *
   * \param [in] isDeepCpyGlobs - a boolean that indicates the behaviour concerning globals (profiles and localizations)
   *                            By default (true) the globals are deeply copied.
   * \return MEDFileField1TS * - a new object that is the result of the conversion of \a this to float64 field.
   */
  template<class T>
  MEDFileField1TS *MEDFileNDTemplateField1TS<T>::convertToDouble(bool isDeepCpyGlobs) const
  {
    MCAuto<MEDFileField1TS> ret;
    const MEDFileAnyTypeField1TSWithoutSDA *content(this->_content);
    if(content)
      {
        const typename MLFieldTraits<T>::F1TSWSDAType *contc(dynamic_cast<const typename MLFieldTraits<T>::F1TSWSDAType *>(content));
        if(!contc)
          {
            std::ostringstream oss; oss << "MEDFileNDTemplateField1TS<T>::convertToDouble : the content inside this is not " << MLFieldTraits<T>::F1TSWSDAType::TYPE_STR << " ! This is incoherent !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
        MCAuto<MEDFileField1TSWithoutSDA> newc(contc->convertToDouble());
        ret=static_cast<MEDFileField1TS *>(MEDFileAnyTypeField1TS::BuildNewInstanceFromContent((MEDFileField1TSWithoutSDA *)newc));
      }
    else
      ret=MEDFileField1TS::New();
    if(isDeepCpyGlobs)
      ret->deepCpyGlobs(*this);
    else
      ret->shallowCpyGlobs(*this);
    return ret.retn();
  }

  //////////////////////////

  template<class T>
  typename MLFieldTraits<T>::FMTSWSDAType *MEDFileTemplateFieldMultiTSWithoutSDA<T>::New(med_idt fid, const std::string& fieldName, const std::string& meshName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities)
  {
    return new typename MLFieldTraits<T>::FMTSWSDAType(fid,fieldName,meshName,fieldTyp,infos,nbOfStep,dtunit,loadAll,ms,entities);
  }
  
  template<class T>
  void MEDFileTemplateFieldMultiTSWithoutSDA<T>::checkCoherencyOfType(const MEDFileAnyTypeField1TSWithoutSDA *f1ts) const
  {
    if(!f1ts)
      throw INTERP_KERNEL::Exception("MEDFileFieldMultiTSWithoutSDA::checkCoherencyOfType : input field1TS is NULL ! Impossible to check !");
    const typename MLFieldTraits<T>::F1TSWSDAType *f1tsC(dynamic_cast<const typename MLFieldTraits<T>::F1TSWSDAType *>(f1ts));
    if(!f1tsC)
      {
        std::ostringstream oss; oss << "MEDFileFieldMultiTSWithoutSDA::checkCoherencyOfType : the input field1TS is not a " << MLFieldTraits<T>::F1TSWSDAType::TYPE_STR << " type !";
        throw INTERP_KERNEL::Exception(oss.str());
      }
  }
  
  template<class T>
  const char *MEDFileTemplateFieldMultiTSWithoutSDA<T>::getTypeStr() const
  {
    return MLFieldTraits<T>::F1TSWSDAType::TYPE_STR;
  }
  
  template<class T>
  MEDFileAnyTypeFieldMultiTSWithoutSDA *MEDFileTemplateFieldMultiTSWithoutSDA<T>::createNew() const
  {
    return new typename MLFieldTraits<T>::FMTSWSDAType;
  }
  
  template<class T>
  MEDFileAnyTypeField1TSWithoutSDA *MEDFileTemplateFieldMultiTSWithoutSDA<T>::createNew1TSWithoutSDAEmptyInstance() const
  {
    return new typename MLFieldTraits<T>::F1TSWSDAType;
  }
  
  //////////////////////////

  template<class T>
  MEDFileFieldMultiTSWithoutSDA *MEDFileNDTemplateFieldMultiTSWithoutSDA<T>::convertToDouble() const
  {
    MCAuto<MEDFileFieldMultiTSWithoutSDA> ret(new MEDFileFieldMultiTSWithoutSDA);
    ret->MEDFileAnyTypeFieldMultiTSWithoutSDA::operator =(*this);
    int i=0;
    for(std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> >::const_iterator it=this->_time_steps.begin();it!=this->_time_steps.end();it++,i++)
      {
        const MEDFileAnyTypeField1TSWithoutSDA *eltToConv(*it);
        if(eltToConv)
          {
            const typename MLFieldTraits<T>::F1TSWSDAType *eltToConvC(dynamic_cast<const typename MLFieldTraits<T>::F1TSWSDAType *>(eltToConv));
            if(!eltToConvC)
              throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTSWithoutSDA::convertToInt : presence of an invalid 1TS type ! Should be of type INT32 !");
          MCAuto<MEDFileAnyTypeField1TSWithoutSDA> elt(eltToConvC->convertToDouble());
          ret->setIteration(i,elt);
          }
      }
    return ret.retn();
  }
  
  //////////////////////////
  
  template<class T>
  MEDFileTemplateFieldMultiTS<T>::MEDFileTemplateFieldMultiTS()
  {
    _content=new typename MLFieldTraits<T>::FMTSWSDAType;
  }

  template<class T>
  MEDFileTemplateFieldMultiTS<T>::MEDFileTemplateFieldMultiTS(med_idt fid, bool loadAll, const MEDFileMeshes *ms):MEDFileAnyTypeFieldMultiTS(fid,loadAll,ms)
  {
  }
  
  template<class T>
  MEDFileTemplateFieldMultiTS<T>::MEDFileTemplateFieldMultiTS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities):MEDFileAnyTypeFieldMultiTS(fid,fieldName,loadAll,ms,entities)
  {
  }
  
  template<class T>
  MEDFileTemplateFieldMultiTS<T>::MEDFileTemplateFieldMultiTS(const typename MLFieldTraits<T>::FMTSWSDAType& other, bool shallowCopyOfContent):MEDFileAnyTypeFieldMultiTS(other,shallowCopyOfContent)
  {
  }

  /*!
   * Return an extraction of \a this using \a extractDef map to specify the extraction.
   * The keys of \a extractDef is level relative to max ext of \a mm mesh.
   *
   * \return A new object that the caller is responsible to deallocate.
   */
  template<class T>
  typename MLFieldTraits<T>::FMTSType *MEDFileTemplateFieldMultiTS<T>::extractPartImpl(const std::map<int, MCAuto<DataArrayInt> >& extractDef, MEDFileMesh *mm) const
  {
    if(!mm)
      throw INTERP_KERNEL::Exception("MEDFileTemplateFieldMultiTS<T>::extractPart : mesh is null !");
    MCAuto<typename MLFieldTraits<T>::FMTSType> fmtsOut(MLFieldTraits<T>::FMTSType::New());
    int nbTS(getNumberOfTS());
    for(int i=0;i<nbTS;i++)
      {
        MCAuto<MEDFileAnyTypeField1TS> f1ts(getTimeStepAtPos(i));
        MCAuto<typename MLFieldTraits<T>::F1TSType> f1ts2(DynamicCastSafe<MEDFileAnyTypeField1TS,typename MLFieldTraits<T>::F1TSType>(f1ts));
        MCAuto<typename MLFieldTraits<T>::F1TSType> f1tsOut(f1ts2->extractPartImpl(extractDef,mm));
        fmtsOut->pushBackTimeStep(f1tsOut);
      }
    return fmtsOut.retn();
  }

  /*!
   * Returns a new empty instance of MEDFileFieldMultiTS.
   *  \return MEDFileFieldMultiTS * - a new instance of MEDFileFieldMultiTS. The caller
   *          is to delete this field using decrRef() as it is no more needed.
   */
  template<class T>
  typename MLFieldTraits<T>::FMTSType *MEDFileTemplateFieldMultiTS<T>::New()
  {
    return new typename MLFieldTraits<T>::FMTSType;
  }

  /*!
   * Returns a new instance of MEDFileTemplateFieldMultiTS<T> holding data of the first field
   * that has been read from a specified MED file.
   *  \param [in] fileName - the name of the MED file to read.
   *  \return MEDFileTemplateFieldMultiTS<T> * - a new instance of MEDFileTemplateFieldMultiTS<T>. The caller
   *          is to delete this field using decrRef() as it is no more needed.
   *  \throw If reading the file fails.
   */
  template<class T>
  typename MLFieldTraits<T>::FMTSType *MEDFileTemplateFieldMultiTS<T>::New(const std::string& fileName, bool loadAll)
  {
    MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
    return New(fid,loadAll);
  }

  template<class T>
  typename MLFieldTraits<T>::FMTSType *MEDFileTemplateFieldMultiTS<T>::New(med_idt fid, bool loadAll)
  {
    MCAuto<typename MLFieldTraits<T>::FMTSType> ret(new typename MLFieldTraits<T>::FMTSType(fid,loadAll,0));
    ret->contentNotNull();//to check that content type matches with \a this type.
    return ret.retn();
  }

  /*!
   * Returns a new instance of MEDFileFieldMultiTS holding data of a given field
   * that has been read from a specified MED file.
   *  \param [in] fileName - the name of the MED file to read.
   *  \param [in] fieldName - the name of the field to read.
   *  \return MEDFileTemplateFieldMultiTS<T> * - a new instance of MEDFileTemplateFieldMultiTS<T>. The caller
   *          is to delete this field using decrRef() as it is no more needed.
   *  \throw If reading the file fails.
   *  \throw If there is no field named \a fieldName in the file.
   */
  template<class T>
  typename MLFieldTraits<T>::FMTSType *MEDFileTemplateFieldMultiTS<T>::New(const std::string& fileName, const std::string& fieldName, bool loadAll)
  {
    MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
    return New(fid,fieldName,loadAll);
  }

  template<class T>
  typename MLFieldTraits<T>::FMTSType *MEDFileTemplateFieldMultiTS<T>::New(med_idt fid, const std::string& fieldName, bool loadAll)
  {
    MCAuto<typename MLFieldTraits<T>::FMTSType> ret(new typename MLFieldTraits<T>::FMTSType(fid,fieldName,loadAll,0));
    ret->contentNotNull();//to check that content type matches with \a this type.
    return ret.retn();
  }

  /*!
   * Returns a new instance of MEDFileFieldMultiTS. If \a shallowCopyOfContent is true the content of \a other is shallow copied.
   * If \a shallowCopyOfContent is false, \a other is taken to be the content of \a this.
   *
   * Returns a new instance of MEDFileTemplateFieldMultiTS<T> holding either a shallow copy
   * of a given MEDFileTemplateFieldMultiTS<T>WithoutSDA ( \a other ) or \a other itself.
   * \warning this is a shallow copy constructor
   *  \param [in] other - a MEDFileField1TSWithoutSDA to copy.
   *  \param [in] shallowCopyOfContent - if \c true, a shallow copy of \a other is created.
   *  \return MEDFileTemplateFieldMultiTS<T> * - a new instance of MEDFileTemplateFieldMultiTS<T>. The caller
   *          is to delete this field using decrRef() as it is no more needed.
   */
  template<class T>
  typename MLFieldTraits<T>::FMTSType *MEDFileTemplateFieldMultiTS<T>::New(const typename MLFieldTraits<T>::FMTSWSDAType& other, bool shallowCopyOfContent)
  {
    return new typename MLFieldTraits<T>::FMTSType(other,shallowCopyOfContent);
  }

  template<class T>
  typename MLFieldTraits<T>::FMTSType *MEDFileTemplateFieldMultiTS<T>::LoadSpecificEntities(const std::string& fileName, const std::string& fieldName, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> >& entities, bool loadAll)
  {
    MEDFileUtilities::AutoFid fid(OpenMEDFileForRead(fileName));
    INTERP_KERNEL::AutoCppPtr<MEDFileEntities> ent(new MEDFileStaticEntities(entities));
    MCAuto<typename MLFieldTraits<T>::FMTSType> ret(new typename MLFieldTraits<T>::FMTSType(fid,fieldName,loadAll,0,ent));
    ret->contentNotNull();//to check that content type matches with \a this type.
    return ret.retn();
  }

  /*!
   * This is the simplest version to fetch a field for MED structure. One drawback : if \a this is a complex field (multi spatial discretization inside a same field) this method will throw exception and more advance
   * method should be called (getFieldOnMeshAtLevel for example).
   * But for normal usage of field in MED file world this method is the most efficient to fetch data.
   *
   * \param [in] iteration - the iteration number of a required time step.
   * \param [in] order - the iteration order number of required time step.
   * \param [in] mesh - the mesh the field is lying on
   * \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
   *          caller is to delete this field using decrRef() as it is no more needed. 
   */
  template<class T>
  typename Traits<T>::FieldType *MEDFileTemplateFieldMultiTS<T>::field(int iteration, int order, const MEDFileMesh *mesh) const
  {
    const MEDFileAnyTypeField1TSWithoutSDA& myF1TS(contentNotNullBase()->getTimeStepEntry(iteration,order));
    MCAuto<DataArray> arrOut;
    MCAuto<MEDCouplingFieldDouble> ret(myF1TS.fieldOnMesh(this,mesh,arrOut,*contentNotNullBase()));
    MCAuto<typename Traits<T>::FieldType> ret2(MEDFileTemplateField1TS<T>::SetDataArrayInField(ret,arrOut));
    return ret2.retn();
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
  template<class T>
  typename Traits<T>::FieldType *MEDFileTemplateFieldMultiTS<T>::getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol) const
  {
    const MEDFileAnyTypeField1TSWithoutSDA& myF1TS(contentNotNullBase()->getTimeStepEntry(iteration,order));
    const typename MLFieldTraits<T>::F1TSWSDAType *myF1TSC(dynamic_cast<const typename MLFieldTraits<T>::F1TSWSDAType *>(&myF1TS));
    if(!myF1TSC)
      throw INTERP_KERNEL::Exception("MEDFileTemplateFieldMultiTS<T>::getFieldAtLevel : mismatch of type of field expecting FLOAT64 !");
    MCAuto<DataArray> arrOut;
    MCAuto<MEDCouplingFieldDouble> ret(myF1TSC->getFieldAtLevel(type,meshDimRelToMax,std::string(),renumPol,this,arrOut,*contentNotNullBase()));
    MCAuto<typename Traits<T>::FieldType> ret2(MEDFileTemplateField1TS<T>::SetDataArrayInField(ret,arrOut));
    return ret2.retn();
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
  template<class T>
  typename Traits<T>::FieldType *MEDFileTemplateFieldMultiTS<T>::getFieldAtTopLevel(TypeOfField type, int iteration, int order, int renumPol) const
  {
    const MEDFileAnyTypeField1TSWithoutSDA& myF1TS(contentNotNullBase()->getTimeStepEntry(iteration,order));
    const typename MLFieldTraits<T>::F1TSWSDAType *myF1TSC(dynamic_cast<const typename MLFieldTraits<T>::F1TSWSDAType *>(&myF1TS));
    if(!myF1TSC)
      throw INTERP_KERNEL::Exception("MEDFileTemplateFieldMultiTS<T>::getFieldAtTopLevel : mismatch of type of field !");
    MCAuto<DataArray> arrOut;
    MCAuto<MEDCouplingFieldDouble> ret(myF1TSC->getFieldAtTopLevel(type,std::string(),renumPol,this,arrOut,*contentNotNullBase()));
    MCAuto<typename Traits<T>::FieldType> ret2(MEDFileTemplateField1TS<T>::SetDataArrayInField(ret,arrOut));
    return ret2.retn();
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
  template<class T>
  typename Traits<T>::FieldType *MEDFileTemplateFieldMultiTS<T>::getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol) const
  {
    const MEDFileAnyTypeField1TSWithoutSDA& myF1TS(contentNotNullBase()->getTimeStepEntry(iteration,order));
    const typename MLFieldTraits<T>::F1TSWSDAType *myF1TSC(dynamic_cast<const typename MLFieldTraits<T>::F1TSWSDAType *>(&myF1TS));
    if(!myF1TSC)
      throw INTERP_KERNEL::Exception("MEDFileTemplateFieldMultiTS<T>::getFieldOnMeshAtLevel : mismatch of type of field !");
    MCAuto<DataArray> arrOut;
    MCAuto<MEDCouplingFieldDouble> ret(myF1TSC->getFieldOnMeshAtLevel(type,meshDimRelToMax,renumPol,this,mesh,arrOut,*contentNotNullBase()));
    MCAuto<typename Traits<T>::FieldType> ret2(MEDFileTemplateField1TS<T>::SetDataArrayInField(ret,arrOut));
    return ret2.retn();
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
  template<class T>
  typename Traits<T>::FieldType *MEDFileTemplateFieldMultiTS<T>::getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, int renumPol) const
  {
    const MEDFileAnyTypeField1TSWithoutSDA& myF1TS(contentNotNullBase()->getTimeStepEntry(iteration,order));
    const typename MLFieldTraits<T>::F1TSWSDAType *myF1TSC(dynamic_cast<const typename MLFieldTraits<T>::F1TSWSDAType *>(&myF1TS));
    if(!myF1TSC)
      throw INTERP_KERNEL::Exception("MEDFileTemplateFieldMultiTS<T>::getFieldOnMeshAtLevel : mismatch of type of field !");
    MCAuto<DataArray> arrOut;
    MCAuto<MEDCouplingFieldDouble> ret(myF1TSC->getFieldOnMeshAtLevel(type,renumPol,this,mesh,0,0,arrOut,*contentNotNullBase()));
    MCAuto<typename Traits<T>::FieldType> ret2(MEDFileTemplateField1TS<T>::SetDataArrayInField(ret,arrOut));
    return ret2.retn();
  }

  /*!
   * This method has a close behaviour than MEDFileFieldMultiTS::getFieldAtLevel.
   * This method is called 'old' because the user should give the mesh name he wants to use for it's field.
   * This method is useful for MED2 file format when field on different mesh was autorized.
   */
  template<class T>
  typename Traits<T>::FieldType *MEDFileTemplateFieldMultiTS<T>::getFieldAtLevelOld(TypeOfField type, int iteration, int order, const std::string& mname, int meshDimRelToMax, int renumPol) const
  {
    const MEDFileAnyTypeField1TSWithoutSDA& myF1TS(contentNotNullBase()->getTimeStepEntry(iteration,order));
    const typename MLFieldTraits<T>::F1TSWSDAType *myF1TSC(dynamic_cast<const typename MLFieldTraits<T>::F1TSWSDAType *>(&myF1TS));
    if(!myF1TSC)
      throw INTERP_KERNEL::Exception("MEDFileTemplateFieldMultiTS<T>::getFieldAtLevelOld : mismatch of type of field !");
    MCAuto<DataArray> arrOut;
    MCAuto<MEDCouplingFieldDouble> ret(myF1TSC->getFieldAtLevel(type,meshDimRelToMax,mname,renumPol,this,arrOut,*contentNotNullBase()));
    MCAuto<typename Traits<T>::FieldType> ret2(MEDFileTemplateField1TS<T>::SetDataArrayInField(ret,arrOut));
    return ret2.retn();
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
  template<class T>
  typename Traits<T>::ArrayType *MEDFileTemplateFieldMultiTS<T>::getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const
  {
    const MEDFileAnyTypeField1TSWithoutSDA& myF1TS(contentNotNullBase()->getTimeStepEntry(iteration,order));
    const typename MLFieldTraits<T>::F1TSWSDAType *myF1TSC(dynamic_cast<const typename MLFieldTraits<T>::F1TSWSDAType *>(&myF1TS));
    if(!myF1TSC)
      throw INTERP_KERNEL::Exception("MEDFileTemplateFieldMultiTS<T>::getFieldWithProfile : mismatch of type of field !");
    MCAuto<DataArray> ret(myF1TSC->getFieldWithProfile(type,meshDimRelToMax,mesh,pfl,this,*contentNotNullBase()));
    return MEDFileTemplateField1TS<T>::ReturnSafelyTypedDataArray(ret);
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
  template<class T>
  void MEDFileTemplateFieldMultiTS<T>::appendFieldNoProfileSBT(const typename Traits<T>::FieldType *field)
  {
    const typename Traits<T>::ArrayType *arr(NULL);
    if(field)
      arr=field->getArray();
    MCAuto<MEDCouplingFieldDouble> field2(MEDFileTemplateField1TS<T>::ToFieldTemplateWithTime(field));
    contentNotNull()->appendFieldNoProfileSBT(field2,arr,*this);
  }

  /*!
   * Adds a MEDCouplingFieldDouble to \a this as another time step.
   * The mesh support of input parameter \a field is ignored here, it can be NULL.
   * The support of field \a field is expected to be those computed with the input parameter \a mesh, \a meshDimRelToMax,
   * and \a profile.
   *
   * This method will check that the field based on the computed support is coherent. If not an exception will be thrown.
   * A new profile is added only if no equal profile is missing.
   * For more info, see \ref AdvMEDLoaderAPIFieldRW
   *  \param [in] field - the field to add to \a this. The mesh support of field is ignored.
   *  \param [in] mesh - the supporting mesh of \a field.
   *  \param [in] meshDimRelToMax - a relative dimension of mesh entities \a field lies on (useless if field spatial discretization is ON_NODES).
   *  \param [in] profile - ids of mesh entities on which corresponding field values lie.
   *  \throw If either \a field or \a mesh or \a profile has an empty name.
   *  \throw If there are no mesh entities of \a meshDimRelToMax dimension in \a mesh.
   *  \throw If the data array of \a field is not set.
   *  \throw If the data array of \a this is already allocated but has different number of
   *         components than \a field.
   *  \throw If elements in \a mesh are not in the order suitable for writing to the MED file.
   *  \sa setFieldNoProfileSBT()
   */
  template<class T>
  void MEDFileTemplateFieldMultiTS<T>::appendFieldProfile(const typename Traits<T>::FieldType *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile)
  {
    const typename Traits<T>::ArrayType *arr(NULL);
    if(field)
      arr=field->getArray();
    MCAuto<MEDCouplingFieldDouble> field2(MEDFileTemplateField1TS<T>::ToFieldTemplateWithTime(field));
    contentNotNull()->appendFieldProfile(field2,arr,mesh,meshDimRelToMax,profile,*this);
  }

  template<class T>
  const typename MLFieldTraits<T>::FMTSWSDAType *MEDFileTemplateFieldMultiTS<T>::contentNotNull() const
  {
    const MEDFileAnyTypeFieldMultiTSWithoutSDA *pt(_content);
    if(!pt)
      throw INTERP_KERNEL::Exception("MEDFileTemplateFieldMultiTS<T>::contentNotNull : the content pointer is null !");
    const typename MLFieldTraits<T>::FMTSWSDAType *ret=dynamic_cast<const typename MLFieldTraits<T>::FMTSWSDAType *>(pt);
    if(!ret)
      throw INTERP_KERNEL::Exception("MEDFileTemplateFieldMultiTS<T>::contentNotNull : the content pointer is not null but it is not of type double ! Reason is maybe that the read field has not the type FLOAT64 !");
    return ret;
  }

  template<class T>
  typename MLFieldTraits<T>::FMTSWSDAType *MEDFileTemplateFieldMultiTS<T>::contentNotNull()
  {
    MEDFileAnyTypeFieldMultiTSWithoutSDA *pt(_content);
    if(!pt)
      throw INTERP_KERNEL::Exception("MEDFileTemplateFieldMultiTS<T>::contentNotNull : the non const content pointer is null !");
    typename MLFieldTraits<T>::FMTSWSDAType *ret(dynamic_cast<typename MLFieldTraits<T>::FMTSWSDAType *>(pt));
    if(!ret)
      throw INTERP_KERNEL::Exception("MEDFileTemplateFieldMultiTS<T>::contentNotNull : the non const content pointer is not null but it is not of type double ! Reason is maybe that the read field has not the type FLOAT64 !");
    return ret;
  }

  /*!
   * Returns a new MEDFileField1TS holding data of a given time step of \a this field.
   *  \param [in] pos - a time step id.
   *  \return MEDFileField1TS * - a new instance of MEDFileField1TS. The caller is to
   *          delete this field using decrRef() as it is no more needed.
   *  \throw If \a pos is not a valid time step id.
   */
  template<class T>
  typename MLFieldTraits<T>::F1TSType *MEDFileTemplateFieldMultiTS<T>::getTimeStepAtPos(int pos) const
  {
    const MEDFileAnyTypeField1TSWithoutSDA *item(contentNotNullBase()->getTimeStepAtPos2(pos));
    if(!item)
      {
        std::ostringstream oss; oss << "MEDFileFieldMultiTS::getTimeStepAtPos : field at pos #" << pos << " is null !";
        throw INTERP_KERNEL::Exception(oss.str());
      }
    const typename MLFieldTraits<T>::F1TSWSDAType *itemC=dynamic_cast<const typename MLFieldTraits<T>::F1TSWSDAType *>(item);
    if(itemC)
      {
        MCAuto<typename MLFieldTraits<T>::F1TSType> ret(MLFieldTraits<T>::F1TSType::New(*itemC,false));
        ret->shallowCpyGlobs(*this);
        return ret.retn();
      }
    std::ostringstream oss; oss << "MEDFileFieldMultiTS::getTimeStepAtPos : type of field at pos #" << pos << " is not " << MLFieldTraits<T>::F1TSWSDAType::TYPE_STR << " !";
    throw INTERP_KERNEL::Exception(oss.str());
  }

  template<class T>
  typename Traits<T>::ArrayType *MEDFileTemplateFieldMultiTS<T>::getUndergroundDataArray(int iteration, int order) const
  {
    DataArray *ret(contentNotNull()->getUndergroundDataArray(iteration,order));
    if(!ret)
      return NULL;
    typename Traits<T>::ArrayType *ret2(dynamic_cast<typename Traits<T>::ArrayType *>(ret));
    if(!ret2)
      {
        std::ostringstream oss; oss << "MEDFileTemplateFieldMultiTS<T>::getUndergroundDataArray : invalid type of data dectected ! Expecting " << MLFieldTraits<T>::F1TSWSDAType::TYPE_STR;
        throw INTERP_KERNEL::Exception(oss.str());
      }
    return ret2;
  }

  template<class T>
  typename Traits<T>::ArrayType *MEDFileTemplateFieldMultiTS<T>::getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const
  {
    DataArray *ret(contentNotNull()->getUndergroundDataArrayExt(iteration,order,entries));
    if(!ret)
      return NULL;
    typename Traits<T>::ArrayType *ret2(dynamic_cast<typename Traits<T>::ArrayType *>(ret));
    if(!ret2)
      {
        std::ostringstream oss; oss << "MEDFileTemplateFieldMultiTS<T>::getUndergroundDataArrayExt : invalid type of data dectected ! Expecting " << MLFieldTraits<T>::F1TSWSDAType::TYPE_STR;
        throw INTERP_KERNEL::Exception(oss.str());
      }
    return ret2;
  }
  
  template<class T>
  typename MLFieldTraits<T>::FMTSType *MEDFileTemplateFieldMultiTS<T>::buildNewEmptyImpl() const
  {
    return MLFieldTraits<T>::FMTSType::New();
  }
  
  template<class T>
  void MEDFileTemplateFieldMultiTS<T>::checkCoherencyOfType(const MEDFileAnyTypeField1TS *f1ts) const
  {
    if(!f1ts)
      throw INTERP_KERNEL::Exception("MEDFileTemplateFieldMultiTS<T>::checkCoherencyOfType : input field1TS is NULL ! Impossible to check !");
    const typename MLFieldTraits<T>::F1TSType *f1tsC=dynamic_cast<const typename MLFieldTraits<T>::F1TSType *>(f1ts);
    if(!f1tsC)
      {
        std::ostringstream oss; oss << "MEDFileTemplateFieldMultiTS<T>::checkCoherencyOfType : the input field1TS is not a " << MLFieldTraits<T>::F1TSWSDAType::TYPE_STR << " type !";
        throw INTERP_KERNEL::Exception(oss.str());
      }
  }

  //////////////////////////

  /*!
   * This method performs a copy with datatype modification ( int32->float64 ) of \a this. The globals information are copied
   * following the given input policy.
   *
   * \param [in] isDeepCpyGlobs - a boolean that indicates the behaviour concerning globals (profiles and localizations)
   *                            By default (true) the globals are deeply copied.
   * \return MEDFileFieldMultiTS * - a new object that is the result of the conversion of \a this to float64 field.
   */
  template<class T>
  MEDFileFieldMultiTS *MEDFileNDTemplateFieldMultiTS<T>::convertToDouble(bool isDeepCpyGlobs) const
  {
    MCAuto<MEDFileFieldMultiTS> ret;
    const MEDFileAnyTypeFieldMultiTSWithoutSDA *content(this->_content);
    if(content)
      {
        const typename MLFieldTraits<T>::FMTSWSDAType *contc=dynamic_cast<const typename MLFieldTraits<T>::FMTSWSDAType *>(content);
        if(!contc)
          throw INTERP_KERNEL::Exception("MEDFileIntFieldMultiTS::convertToInt : the content inside this is not INT32 ! This is incoherent !");
        MCAuto<MEDFileFieldMultiTSWithoutSDA> newc(contc->convertToDouble());
        ret=static_cast<MEDFileFieldMultiTS *>(MEDFileAnyTypeFieldMultiTS::BuildNewInstanceFromContent((MEDFileFieldMultiTSWithoutSDA *)newc));
      }
    else
      ret=MEDFileFieldMultiTS::New();
    if(isDeepCpyGlobs)
      ret->deepCpyGlobs(*this);
    else
      ret->shallowCpyGlobs(*this);
    return ret.retn();
  }
}

#endif
