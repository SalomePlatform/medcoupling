// Copyright (C) 2016  CEA/DEN, EDF R&D
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

#ifndef __MEDCOUPLINGFIELDT_TXX__
#define __MEDCOUPLINGFIELDT_TXX__

#include "MEDCouplingTimeDiscretization.hxx"
#include "MEDCouplingMesh.hxx"

namespace MEDCoupling
{
  template<class T>
  MEDCouplingFieldT<T>::MEDCouplingFieldT(const MEDCouplingFieldT<T>& other, bool deepCopy):MEDCouplingField(other,deepCopy),_time_discr(other._time_discr->performCopyOrIncrRef(deepCopy))
  {
  }
  
  /*!
   * Checks if \a this field is correctly defined, else an exception is thrown.
   *  \throw If the mesh is not set.
   *  \throw If the data array is not set.
   *  \throw If the spatial discretization of \a this field is NULL.
   *  \throw If \a this->getTimeTolerance() < 0.
   *  \throw If the temporal discretization data is incorrect.
   *  \throw If mesh data does not correspond to field data.
   */
  template<class T>
  void MEDCouplingFieldT<T>::checkConsistencyLight() const
  {
    MEDCouplingField::checkConsistencyLight();
    _time_discr->checkConsistencyLight();
    _type->checkCoherencyBetween(_mesh,getArray());
  }

  template<class T>
  MEDCouplingFieldT<T>::MEDCouplingFieldT(const MEDCouplingField& other, MEDCouplingTimeDiscretizationTemplate<T> *timeDiscr, bool deepCopy):MEDCouplingField(other,deepCopy),_time_discr(timeDiscr)
  {
  }
  
  template<class T>
  MEDCouplingFieldT<T>::MEDCouplingFieldT(TypeOfField type, MEDCouplingTimeDiscretizationTemplate<T> *timeDiscr):MEDCouplingField(type),_time_discr(timeDiscr)
  {
  }

  template<class T>
  MEDCouplingFieldT<T>::MEDCouplingFieldT(MEDCouplingFieldDiscretization *type, NatureOfField n, MEDCouplingTimeDiscretizationTemplate<T> *timeDiscr):MEDCouplingField(type,n),_time_discr(timeDiscr)
  {
  }

  template<class T>
  MEDCouplingFieldT<T>::~MEDCouplingFieldT()
  {
    delete _time_discr;
  }

  /*!
   * This method synchronizes time information (time, iteration, order, time unit) regarding the information in \c this->_mesh.
   * \throw If no mesh is set in this. Or if \a this is not compatible with time setting (typically NO_TIME)
   */
  template<class T>
  void MEDCouplingFieldT<T>::synchronizeTimeWithMesh()
  {
    if(!_mesh)
      throw INTERP_KERNEL::Exception("MEDCouplingFieldT::synchronizeTimeWithMesh : no mesh set in this !");
    int it(-1),ordr(-1);
    double val(_mesh->getTime(it,ordr));
    std::string timeUnit(_mesh->getTimeUnit());
    setTime(val,it,ordr);
    setTimeUnit(timeUnit);
  }

  /*!
   * Returns a new MEDCouplingFieldDouble which is a copy of \a this one. The data
   * of \a this field is copied either deep or shallow depending on \a recDeepCpy
   * parameter. But the underlying mesh is always deep copied.
   * Data that can be copied either deeply or shallow are:
   * - \ref MEDCouplingTemporalDisc "temporal discretization" data that holds array(s)
   * of field values,
   * - \ref MEDCouplingSpatialDisc "a spatial discretization".
   * 
   * This method behaves exactly like clone() except that here the underlying **mesh is
   * always deeply duplicated**, whatever the value \a recDeepCpy parameter.
   * The result of \c cloneWithMesh(true) is exactly the same as that of deepCopy().
   * So the resulting field can not be used together with \a this one in the methods
   * like operator+(), operator*() etc. To avoid deep copying the underlying mesh,
   * the user can call clone().
   *  \param [in] recDeepCpy - if \c true, the copy of the underlying data arrays is
   *         deep, else all data arrays of \a this field are shared by the new field.
   *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
   *         caller is to delete this field using decrRef() as it is no more needed.
   * \sa clone()
   */
  template<class T>
  typename Traits<T>::FieldType *MEDCouplingFieldT<T>::cloneWithMesh(bool recDeepCpy) const
  {
    MCAuto<  typename Traits<T>::FieldType > ret(clone(recDeepCpy));
    if(_mesh)
      {
        MCAuto<MEDCouplingMesh> mCpy(_mesh->deepCopy());
        ret->setMesh(mCpy);
      }
    return ret.retn();
  }

  template<class T>
  bool MEDCouplingFieldT<T>::isEqual(const MEDCouplingFieldT<T> *other, double meshPrec, T valsPrec) const
  {
    std::string tmp;
    return isEqualIfNotWhy(other,meshPrec,valsPrec,tmp);
  }
  
  template<class T>
  bool MEDCouplingFieldT<T>::isEqualIfNotWhy(const MEDCouplingFieldT<T> *other, double meshPrec, T valsPrec, std::string& reason) const
  {
    if(!other)
      throw INTERP_KERNEL::Exception("MEDCouplingFieldT::isEqualIfNotWhy : other instance is NULL !");
    if(!isEqualIfNotWhyProtected(other,meshPrec,reason))
      return false;
    if(!_time_discr->isEqualIfNotWhy(other->_time_discr,T(valsPrec),reason))
      {
        reason.insert(0,"In FieldT time discretizations differ :");
        return false;
      }
    return true;
  }
  
  /*!
   * Checks equality of \a this and \a other field. Only numeric data is considered,
   * i.e. names, description etc are not compared.
   *  \param [in] other - the field to compare with.
   *  \param [in] meshPrec - a precision used to compare node coordinates of meshes.
   *  \param [in] valsPrec - a precision used to compare data arrays of the two fields.
   *  \return bool - \c true if the two fields are equal, \c false else.
   *  \throw If \a other == NULL.
   *  \throw If the spatial discretization of \a this field is NULL.
   */
  template<class T>
  bool MEDCouplingFieldT<T>::isEqualWithoutConsideringStr(const MEDCouplingFieldT<T> *other, double meshPrec, T valsPrec) const
  {
    if(!other)
      return false;
    if(!isEqualWithoutConsideringStrProtected(other,meshPrec))
      return false;
    if(!_time_discr->isEqualWithoutConsideringStr(other->_time_discr,valsPrec))
      return false;
    return true;
  }
  
  /*!
   * Copies tiny info (component names, name and description) from an \a other field to
   * \a this one.
   * \warning The underlying mesh is not renamed (for safety reason).
   *  \param [in] other - the field to copy the tiny info from.
   *  \throw If \a this->getNumberOfComponents() != \a other->getNumberOfComponents()
   */
  template<class T>
  void MEDCouplingFieldT<T>::copyTinyStringsFrom(const MEDCouplingField *other)
  {
    MEDCouplingField::copyTinyStringsFrom(other);
    const MEDCouplingFieldT<T> *otherC(dynamic_cast<const MEDCouplingFieldT<T> *>(other));
    if(otherC)
      {
        _time_discr->copyTinyStringsFrom(*otherC->_time_discr);
      }
  }

  /*!
   * Copies only times, order and iteration from an \a other field to
   * \a this one. The underlying mesh is not impacted by this method.
   * Arrays are not impacted neither.
   *  \param [in] other - the field to tiny attributes from.
   *  \throw If \a this->getNumberOfComponents() != \a other->getNumberOfComponents()
   */
  template<class T>
  void MEDCouplingFieldT<T>::copyTinyAttrFrom(const MEDCouplingFieldT<T> *other)
  {
    if(other)
      {
        _time_discr->copyTinyAttrFrom(*other->_time_discr);
      }
  }
  
  template<class T>
  void MEDCouplingFieldT<T>::copyAllTinyAttrFrom(const MEDCouplingFieldT<T> *other)
  {
    copyTinyStringsFrom(other);
    copyTinyAttrFrom(other);
  }
  
  /*!
   * Permutes values of \a this field according to a given permutation array for cells
   * renumbering. The underlying mesh is deeply copied and its cells are also permuted. 
   * The number of cells remains the same; for that the permutation array \a old2NewBg
   * should not contain equal ids.
   * ** Warning, this method modifies the mesh aggreagated by \a this (by performing a deep copy ) **.
   *
   *  \param [in] old2NewBg - the permutation array in "Old to New" mode. Its length is
   *         to be equal to \a this->getMesh()->getNumberOfCells().
   *  \param [in] check - if \c true, \a old2NewBg is transformed to a new permutation
   *         array, so that its maximal cell id to correspond to (be less than) the number
   *         of cells in mesh. This new array is then used for the renumbering. If \a 
   *         check == \c false, \a old2NewBg is used as is, that is less secure as validity 
   *         of ids in \a old2NewBg is not checked.
   *  \throw If the mesh is not set.
   *  \throw If the spatial discretization of \a this field is NULL.
   *  \throw If \a check == \c true and \a old2NewBg contains equal ids.
   *  \throw If mesh nature does not allow renumbering (e.g. structured mesh).
   * 
   *  \if ENABLE_EXAMPLES
   *  \ref cpp_mcfielddouble_renumberCells "Here is a C++ example".<br>
   *  \ref  py_mcfielddouble_renumberCells "Here is a Python example".
   *  \endif
   */
  template<class T>
  void MEDCouplingFieldT<T>::renumberCells(const int *old2NewBg, bool check)
  {
    renumberCellsWithoutMesh(old2NewBg,check);
    MCAuto<MEDCouplingMesh> m(_mesh->deepCopy());
    m->renumberCells(old2NewBg,check);
    setMesh(m);
    updateTime();
  }

  /*!
   * Permutes values of \a this field according to a given permutation array for cells
   * renumbering. The underlying mesh is \b not permuted. 
   * The number of cells remains the same; for that the permutation array \a old2NewBg
   * should not contain equal ids.
   * This method performs a part of job of renumberCells(). The reasonable use of this
   * method is only for multi-field instances lying on the same mesh to avoid a
   * systematic duplication and renumbering of _mesh attribute. 
   * \warning Use this method with a lot of care!
   *  \param [in] old2NewBg - the permutation array in "Old to New" mode. Its length is
   *         to be equal to \a this->getMesh()->getNumberOfCells().
   *  \param [in] check - if \c true, \a old2NewBg is transformed to a new permutation
   *         array, so that its maximal cell id to correspond to (be less than) the number
   *         of cells in mesh. This new array is then used for the renumbering. If \a 
   *         check == \c false, \a old2NewBg is used as is, that is less secure as validity 
   *         of ids in \a old2NewBg is not checked.
   *  \throw If the mesh is not set.
   *  \throw If the spatial discretization of \a this field is NULL.
   *  \throw If \a check == \c true and \a old2NewBg contains equal ids.
   *  \throw If mesh nature does not allow renumbering (e.g. structured mesh).
   */
  template<class T>
  void MEDCouplingFieldT<T>::renumberCellsWithoutMesh(const int *old2NewBg, bool check)
  {
    if(!_mesh)
      throw INTERP_KERNEL::Exception("Expecting a defined mesh to be able to operate a renumbering !");
    if(_type.isNull())
      throw INTERP_KERNEL::Exception("Expecting a spatial discretization to be able to operate a renumbering !");
    //
    _type->renumberCells(old2NewBg,check);
    std::vector< typename MEDCoupling::Traits<T>::ArrayType *> arrays;
    timeDiscrSafe()->getArrays(arrays);
    std::vector<DataArray *> arrays2(arrays.size()); std::copy(arrays.begin(),arrays.end(),arrays2.begin());
    _type->renumberArraysForCell(_mesh,arrays2,old2NewBg,check);
    //
    updateTime();
  }
  
  /*!
   * This method is more strict than MEDCouplingField::areCompatibleForMerge method.
   * This method is used for operation on fields to operate a first check before attempting operation.
   */
  template<class T>
  bool MEDCouplingFieldT<T>::areStrictlyCompatible(const MEDCouplingField *other) const
  {
    std::string tmp;
    if(!MEDCouplingField::areStrictlyCompatible(other))
      return false;
    const MEDCouplingFieldT<T> *otherC(dynamic_cast<const MEDCouplingFieldT<T> *>(other));
    if(!otherC)
      return false;
    if(!_time_discr->areStrictlyCompatible(otherC->_time_discr,tmp))
      return false;
    return true;
  }
  
  template<class T>
  bool MEDCouplingFieldT<T>::areStrictlyCompatibleForMulDiv(const MEDCouplingField *other) const
  {
    if(!MEDCouplingField::areStrictlyCompatibleForMulDiv(other))
      return false;
    const MEDCouplingFieldT<T> *otherC(dynamic_cast<const MEDCouplingFieldT<T> *>(other));
    if(!otherC)
      return false;
    if(!_time_discr->areStrictlyCompatibleForDiv(otherC->_time_discr))
      return false;
    return true;
  }

  /*!
   * Method with same principle than MEDCouplingFieldDouble::areStrictlyCompatibleForMulDiv method except that
   * number of components between \a this and 'other' can be different here (for operator/).
   */
  template<class T>
  bool MEDCouplingFieldT<T>::areCompatibleForDiv(const MEDCouplingField *other) const
  {
    if(!MEDCouplingField::areStrictlyCompatibleForMulDiv(other))
      return false;
    const MEDCouplingFieldT<T> *otherC(dynamic_cast<const MEDCouplingFieldT<T> *>(other));
    if(!otherC)
      return false;
    if(!_time_discr->areStrictlyCompatibleForDiv(otherC->_time_discr))
      return false;
    return true;
  }
  
  template<class T>
  bool MEDCouplingFieldT<T>::areCompatibleForMul(const MEDCouplingField *other) const
  {
    if(!MEDCouplingField::areStrictlyCompatibleForMulDiv(other))
      return false;
    const MEDCouplingFieldT<T> *otherC(dynamic_cast<const MEDCouplingFieldT<T> *>(other));
    if(!otherC)
      return false;
    if(!_time_discr->areStrictlyCompatibleForMul(otherC->_time_discr))
      return false;
    return true;
  }

  /*!
   * Returns a string describing \a this field. This string is outputted by \c print
   * Python command. The string includes info on
   * - name,
   * - description,
   * - \ref MEDCouplingSpatialDisc "spatial discretization",
   * - \ref MEDCouplingTemporalDisc "time discretization",
   * - \ref NatureOfField,
   * - components,
   * - mesh.
   *
   *  \return std::string - the string describing \a this field.
   */
  template<class T>
  std::string MEDCouplingFieldT<T>::simpleRepr() const
  {
    std::ostringstream ret;
    ret << Traits<T>::FieldTypeName << " with name : \"" << getName() << "\"\n";
    ret << "Description of field is : \"" << getDescription() << "\"\n";
    if(_type)
      { ret << Traits<T>::FieldTypeName << " space discretization is : " << _type->getStringRepr() << "\n"; }
    else
      { ret << Traits<T>::FieldTypeName << " has no spatial discretization !\n"; }
    if(_time_discr)
      { ret << Traits<T>::FieldTypeName << " time discretization is : " << _time_discr->getStringRepr() << "\n"; }
    else
      { ret << Traits<T>::FieldTypeName << " has no time discretization !\n"; }
    ret << Traits<T>::FieldTypeName << " nature of field is : \"" << MEDCouplingNatureOfField::GetReprNoThrow(_nature) << "\"\n";
    if(getArray())
      {
        if(getArray()->isAllocated())
          {
            int nbOfCompo=getArray()->getNumberOfComponents();
            ret << Traits<T>::FieldTypeName << " default array has " << nbOfCompo << " components and " << getArray()->getNumberOfTuples() << " tuples.\n";
            ret << Traits<T>::FieldTypeName << " default array has following info on components : ";
            for(int i=0;i<nbOfCompo;i++)
              ret << "\"" << getArray()->getInfoOnComponent(i) << "\" ";
            ret << "\n";
          }
        else
          {
            ret << "Array set but not allocated !\n";
          }
      }
    if(_mesh)
      ret << "Mesh support information :\n__________________________\n" << _mesh->simpleRepr();
    else
      ret << "Mesh support information : No mesh set !\n";
    return ret.str();
  }
  
  template<class T>
  void MEDCouplingFieldT<T>::reprQuickOverview(std::ostream& stream) const
  {
    stream << Traits<T>::FieldTypeName << " C++ instance at " << this << ". Name : \"" << _name << "\"." << std::endl;
    const char *nat(0);
    try
      {
        nat=MEDCouplingNatureOfField::GetRepr(_nature);
        stream << "Nature of field : " << nat << ".\n";
      }
    catch(INTERP_KERNEL::Exception& e)
      {  }
    const MEDCouplingFieldDiscretization *fd(_type);
    if(!fd)
      stream << "No spatial discretization set !";
    else
      fd->reprQuickOverview(stream);
    stream << std::endl;
    if(!_mesh)
      stream << "\nNo mesh support defined !";
    else
      {
        std::ostringstream oss;
        _mesh->reprQuickOverview(oss);
        std::string tmp(oss.str());
        stream << "\nMesh info : " << tmp.substr(0,tmp.find('\n'));
      }
    if(_time_discr)
      {
        const typename Traits<T>::ArrayType *arr(_time_discr->getArray());
        if(arr)
          {
            stream << "\n\nArray info : ";
            arr->reprQuickOverview(stream);
          }
        else
          {
            stream << "\n\nNo data array set !";
          }
      }
  }

  /*!
   * Returns a type of \ref MEDCouplingTemporalDisc "time discretization" of \a this field.
   *  \return MEDCoupling::TypeOfTimeDiscretization - an enum item describing the time
   *          discretization type.
   */
  template<class T>
  TypeOfTimeDiscretization MEDCouplingFieldT<T>::getTimeDiscretization() const
  {
    return _time_discr->getEnum();
  }

  /*!
   * Builds a newly created field, that the caller will have the responsability to deal with.
   * \n This method makes the assumption that \a this field is correctly defined when this method is called (\a this->checkConsistencyLight() returns without any exception thrown), **no check of this will be done**.
   * \n This method returns a restriction of \a this so that only tuple ids specified in [ \a partBg , \a partEnd ) will be contained in the returned field.
   * \n Parameter [\a partBg, \a partEnd ) specifies **cell ids whatever the spatial discretization** of \a this (
   * \ref MEDCoupling::ON_CELLS "ON_CELLS", 
   * \ref MEDCoupling::ON_NODES "ON_NODES",
   * \ref MEDCoupling::ON_GAUSS_PT "ON_GAUSS_PT", 
   * \ref MEDCoupling::ON_GAUSS_NE "ON_GAUSS_NE",
   * \ref MEDCoupling::ON_NODES_KR "ON_NODES_KR").
   *
   * For example, \a this is a field on cells lying on a mesh that have 10 cells, \a partBg contains the following cell ids [3,7,6].
   * Then the returned field will lie on mesh having 3 cells and will contain 3 tuples.
   *- Tuple #0 of the result field will refer to the cell #0 of returned mesh. The cell #0 of returned mesh will be equal to the cell #3 of \a this->getMesh().
   *- Tuple #1 of the result field will refer to the cell #1 of returned mesh. The cell #1 of returned mesh will be equal to the cell #7 of \a this->getMesh().
   *- Tuple #2 of the result field will refer to the cell #2 of returned mesh. The cell #2 of returned mesh will be equal to the cell #6 of \a this->getMesh().
   *
   * Let, for example, \a this be a field on nodes lying on a mesh that have 10 cells and 11 nodes, and \a partBg contains following cellIds [3,7,6].
   * Thus \a this currently contains 11 tuples. If the restriction of mesh to 3 cells leads to a mesh with 6 nodes, then the returned field
   * will contain 6 tuples and \a this field will lie on this restricted mesh. 
   *
   * \param [in] partBg - start (included) of input range of cell ids to select [ \a partBg, \a partEnd )
   * \param [in] partEnd - end (not included) of input range of cell ids to select [ \a partBg, \a partEnd )
   * \return a newly allocated field the caller should deal with.
   * 
   * \throw if there is presence of an invalid cell id in [ \a partBg, \a partEnd ) regarding the number of cells of \a this->getMesh().
   *
   * \if ENABLE_EXAMPLES
   * \ref cpp_mcfielddouble_subpart1 "Here a C++ example."<br>
   * \ref py_mcfielddouble_subpart1 "Here a Python example."
   * \endif
   * \sa MEDCoupling::MEDCouplingFieldDouble::buildSubPart(const DataArrayInt *) const, MEDCouplingFieldDouble::buildSubPartRange
   */
  template<class T>
  typename Traits<T>::FieldType *MEDCouplingFieldT<T>::buildSubPart(const int *partBg, const int *partEnd) const
  {
    if(_type.isNull())
      throw INTERP_KERNEL::Exception("MEDCouplingFieldT::buildSubPart : Expecting a not NULL spatial discretization !");
    DataArrayInt *arrSelect;
    MCAuto<MEDCouplingMesh> m=_type->buildSubMeshData(_mesh,partBg,partEnd,arrSelect);
    MCAuto<DataArrayInt> arrSelect2(arrSelect);
    MCAuto< typename Traits<T>::FieldType > ret(clone(false));//quick shallow copy.
    const MEDCouplingFieldDiscretization *disc=getDiscretization();
    if(disc)
      ret->setDiscretization(MCAuto<MEDCouplingFieldDiscretization>(disc->clonePart(partBg,partEnd)));
    ret->setMesh(m);
    std::vector<typename Traits<T>::ArrayType *> arrays;
    timeDiscrSafe()->getArrays(arrays);
    std::vector<typename Traits<T>::ArrayType *> arrs;
    std::vector< MCAuto< typename Traits<T>::ArrayType > > arrsSafe;
    const int *arrSelBg=arrSelect->begin();
    const int *arrSelEnd=arrSelect->end();
    for(typename std::vector<typename Traits<T>::ArrayType *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
      {
        typename Traits<T>::ArrayType *arr(0);
        if(*iter)
          arr=(*iter)->selectByTupleIdSafe(arrSelBg,arrSelEnd);
        arrs.push_back(arr); arrsSafe.push_back(arr);
      }
    ret->timeDiscrSafe()->setArrays(arrs,0);
    return ret.retn();
  }

  /*!
   * Builds a newly created field, that the caller will have the responsability to deal with (decrRef()).
   * This method makes the assumption that the field is correctly defined when this method is called, no check of this will be done.
   * This method returns a restriction of \a this so that only tuples with ids specified in \a part will be contained in the returned field.
   * Parameter \a part specifies **cell ids whatever the spatial discretization of this** (
   * \ref MEDCoupling::ON_CELLS "ON_CELLS", 
   * \ref MEDCoupling::ON_NODES "ON_NODES",
   * \ref MEDCoupling::ON_GAUSS_PT "ON_GAUSS_PT", 
   * \ref MEDCoupling::ON_GAUSS_NE "ON_GAUSS_NE",
   * \ref MEDCoupling::ON_NODES_KR "ON_NODES_KR").
   *
   * For example, \a this is a field on cells lying on a mesh that have 10 cells, \a part contains following cell ids [3,7,6].
   * Then the returned field will lie on mesh having 3 cells and the returned field will contain 3 tuples.<br>
   * Tuple #0 of the result field will refer to the cell #0 of returned mesh. The cell #0 of returned mesh will be equal to the cell #3 of \a this->getMesh().<br>
   * Tuple #1 of the result field will refer to the cell #1 of returned mesh. The cell #1 of returned mesh will be equal to the cell #7 of \a this->getMesh().<br>
   * Tuple #2 of the result field will refer to the cell #2 of returned mesh. The cell #2 of returned mesh will be equal to the cell #6 of \a this->getMesh().
   *
   * Let, for example, \a this be a field on nodes lying on a mesh that have 10 cells and 11 nodes, and \a part contains following cellIds [3,7,6].
   * Thus \a this currently contains 11 tuples. If the restriction of mesh to 3 cells leads to a mesh with 6 nodes, then the returned field
   * will contain 6 tuples and \a this field will lie on this restricted mesh. 
   *
   *  \param [in] part - an array of cell ids to include to the result field.
   *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The caller is to delete this field using decrRef() as it is no more needed.
   *
   *  \if ENABLE_EXAMPLES
   *  \ref cpp_mcfielddouble_subpart1 "Here is a C++ example".<br>
   *  \ref  py_mcfielddouble_subpart1 "Here is a Python example".
   *  \endif
   *  \sa MEDCouplingFieldDouble::buildSubPartRange
   */
  template<class T>
  typename Traits<T>::FieldType *MEDCouplingFieldT<T>::buildSubPart(const DataArrayInt *part) const
  {
    if(part==0)
      throw INTERP_KERNEL::Exception("MEDCouplingFieldT::buildSubPart : not empty array must be passed to this method !");
    return buildSubPart(part->begin(),part->end());
  }
  
  /*!
   * This method is equivalent to MEDCouplingFieldDouble::buildSubPart, the only difference is that the input range of cell ids is
   * given using a range given \a begin, \a end and \a step to optimize the part computation.
   * 
   * \sa MEDCouplingFieldDouble::buildSubPart
   */
  template<class T>
  typename Traits<T>::FieldType *MEDCouplingFieldT<T>::buildSubPartRange(int begin, int end, int step) const
  {
    if(_type.isNull())
      throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::buildSubPart : Expecting a not NULL spatial discretization !");
    DataArrayInt *arrSelect;
    int beginOut,endOut,stepOut;
    MCAuto<MEDCouplingMesh> m(_type->buildSubMeshDataRange(_mesh,begin,end,step,beginOut,endOut,stepOut,arrSelect));
    MCAuto<DataArrayInt> arrSelect2(arrSelect);
    MCAuto< typename Traits<T>::FieldType > ret(clone(false));//quick shallow copy.
    const MEDCouplingFieldDiscretization *disc=getDiscretization();
    if(disc)
      ret->setDiscretization(MCAuto<MEDCouplingFieldDiscretization>(disc->clonePartRange(begin,end,step)));
    ret->setMesh(m);
    std::vector<typename Traits<T>::ArrayType *> arrays;
    timeDiscrSafe()->getArrays(arrays);
    std::vector<typename Traits<T>::ArrayType *> arrs;
    std::vector< MCAuto< typename Traits<T>::ArrayType > > arrsSafe;
    for(typename std::vector<typename Traits<T>::ArrayType *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
      {
        typename Traits<T>::ArrayType *arr(0);
        if(*iter)
          {
            if(arrSelect)
              {
                const int *arrSelBg=arrSelect->begin();
                const int *arrSelEnd=arrSelect->end();
                arr=(*iter)->selectByTupleIdSafe(arrSelBg,arrSelEnd);
              }
            else
              arr=(*iter)->selectByTupleIdSafeSlice(beginOut,endOut,stepOut);
          }
        arrs.push_back(arr); arrsSafe.push_back(arr);
      }
    ret->timeDiscrSafe()->setArrays(arrs,0);
    return ret.retn();
  }
  
  template<class T>
  const MEDCouplingTimeDiscretizationTemplate<T> *MEDCouplingFieldT<T>::timeDiscrSafe() const
  {
    const MEDCouplingTimeDiscretizationTemplate<T> *ret(_time_discr);
    if(!ret)
      throw INTERP_KERNEL::Exception("const FieldT : Null type of time discr !");
    return ret;
  }
  
  template<class T>
  MEDCouplingTimeDiscretizationTemplate<T> *MEDCouplingFieldT<T>::timeDiscrSafe()
  {
    MEDCouplingTimeDiscretizationTemplate<T> *ret(_time_discr);
    if(!ret)
      throw INTERP_KERNEL::Exception("const FieldT : Null type of time discr !");
    return ret;
  }

  template<class T>
  void MEDCouplingFieldT<T>::getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const
  {
    tinyInfo.clear();
    timeDiscrSafe()->getTinySerializationStrInformation(tinyInfo);
    tinyInfo.push_back(_name);
    tinyInfo.push_back(_desc);
    tinyInfo.push_back(getTimeUnit());
  }

  /*!
   * This method retrieves some critical values to resize and prepare remote instance.
   * The first two elements returned in tinyInfo correspond to the parameters to give in constructor.
   * @param tinyInfo out parameter resized correctly after the call. The length of this vector is tiny.
   */
  template<class T>
  void MEDCouplingFieldT<T>::getTinySerializationIntInformation(std::vector<int>& tinyInfo) const
  {
    if(_type.isNull())
      throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform getTinySerializationIntInformation !");
    tinyInfo.clear();
    tinyInfo.push_back((int)_type->getEnum());
    tinyInfo.push_back((int)timeDiscrSafe()->getEnum());
    tinyInfo.push_back((int)_nature);
    timeDiscrSafe()->getTinySerializationIntInformation(tinyInfo);
    std::vector<int> tinyInfo2;
    _type->getTinySerializationIntInformation(tinyInfo2);
    tinyInfo.insert(tinyInfo.end(),tinyInfo2.begin(),tinyInfo2.end());
    tinyInfo.push_back((int)tinyInfo2.size());
  }

  /*!
   * This method retrieves some critical values to resize and prepare remote instance.
   * @param tinyInfo out parameter resized correctly after the call. The length of this vector is tiny.
   */
  template<class T>
  void MEDCouplingFieldT<T>::getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const
  {
    if(_type.isNull())
      throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform getTinySerializationDbleInformation !");
    tinyInfo.clear();
    timeDiscrSafe()->getTinySerializationDbleInformation(tinyInfo);
    std::vector<double> tinyInfo2;
    _type->getTinySerializationDbleInformation(tinyInfo2);
    tinyInfo.insert(tinyInfo.end(),tinyInfo2.begin(),tinyInfo2.end());
    tinyInfo.push_back((int)tinyInfo2.size());//very bad, lack of time to improve it
  }

  /*!
   * This method has to be called to the new instance filled by CORBA, MPI, File...
   * @param tinyInfoI is the value retrieves from distant result of getTinySerializationIntInformation on source instance to be copied.
   * @param dataInt out parameter. If not null the pointer is already owned by \a this after the call of this method. In this case no decrRef must be applied.
   * @param arrays out parameter is a vector resized to the right size. The pointers in the vector is already owned by \a this after the call of this method.
   *               No decrRef must be applied to every instances in returned vector.
   * \sa checkForUnserialization
   */
  template<class T>
  void MEDCouplingFieldT<T>::resizeForUnserialization(const std::vector<int>& tinyInfoI, DataArrayInt *&dataInt, std::vector<typename Traits<T>::ArrayType *>& arrays)
  {
    if(_type.isNull())
      throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform resizeForUnserialization !");
    dataInt=0;
    std::vector<int> tinyInfoITmp(tinyInfoI);
    int sz=tinyInfoITmp.back();
    tinyInfoITmp.pop_back();
    std::vector<int> tinyInfoITmp2(tinyInfoITmp.begin(),tinyInfoITmp.end()-sz);
    std::vector<int> tinyInfoI2(tinyInfoITmp2.begin()+3,tinyInfoITmp2.end());
    timeDiscrSafe()->resizeForUnserialization(tinyInfoI2,arrays);
    std::vector<int> tinyInfoITmp3(tinyInfoITmp.end()-sz,tinyInfoITmp.end());
    _type->resizeForUnserialization(tinyInfoITmp3,dataInt);
  }

  /*!
   * This method is extremely close to resizeForUnserialization except that here the arrays in \a dataInt and in \a arrays are attached in \a this
   * after having checked that size is correct. This method is used in python pickeling context to avoid copy of data.
   * \sa resizeForUnserialization
   */
  template<class T>
  void MEDCouplingFieldT<T>::checkForUnserialization(const std::vector<int>& tinyInfoI, const DataArrayInt *dataInt, const std::vector<typename Traits<T>::ArrayType *>& arrays)
  {
    if(_type.isNull())
      throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform resizeForUnserialization !");
    std::vector<int> tinyInfoITmp(tinyInfoI);
    int sz=tinyInfoITmp.back();
    tinyInfoITmp.pop_back();
    std::vector<int> tinyInfoITmp2(tinyInfoITmp.begin(),tinyInfoITmp.end()-sz);
    std::vector<int> tinyInfoI2(tinyInfoITmp2.begin()+3,tinyInfoITmp2.end());
    timeDiscrSafe()->checkForUnserialization(tinyInfoI2,arrays);
    std::vector<int> tinyInfoITmp3(tinyInfoITmp.end()-sz,tinyInfoITmp.end());
    _type->checkForUnserialization(tinyInfoITmp3,dataInt);
  }

  template<class T>
  void MEDCouplingFieldT<T>::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS)
  {
    if(_type.isNull())
      throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform finishUnserialization !");
    std::vector<int> tinyInfoI2(tinyInfoI.begin()+3,tinyInfoI.end());
    //
    std::vector<double> tmp(tinyInfoD);
    int sz=(int)tinyInfoD.back();//very bad, lack of time to improve it
    tmp.pop_back();
    std::vector<double> tmp1(tmp.begin(),tmp.end()-sz);
    std::vector<double> tmp2(tmp.end()-sz,tmp.end());
    //
    timeDiscrSafe()->finishUnserialization(tinyInfoI2,tmp1,tinyInfoS);
    _nature=(NatureOfField)tinyInfoI[2];
    _type->finishUnserialization(tmp2);
    int nbOfElemS=(int)tinyInfoS.size();
    _name=tinyInfoS[nbOfElemS-3];
    _desc=tinyInfoS[nbOfElemS-2];
    setTimeUnit(tinyInfoS[nbOfElemS-1]);
  }

  /*!
   * Contrary to MEDCouplingPointSet class the returned arrays are \b not the responsabilities of the caller.
   * The values returned must be consulted only in readonly mode.
   */
  template<class T>
  void MEDCouplingFieldT<T>::serialize(DataArrayInt *&dataInt, std::vector<typename Traits<T>::ArrayType *>& arrays) const
  {
    if(_type.isNull())
      throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform serialize !");
    timeDiscrSafe()->getArrays(arrays);
    _type->getSerializationIntArray(dataInt);
  }
}

#endif
