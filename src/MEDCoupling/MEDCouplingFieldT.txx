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
  bool MEDCouplingFieldT<T>::isEqualIfNotWhy(const MEDCouplingField *other, double meshPrec, double valsPrec, std::string& reason) const
  {
    if(!other)
      throw INTERP_KERNEL::Exception("MEDCouplingFieldT::isEqualIfNotWhy : other instance is NULL !");
    const MEDCouplingFieldT<T> *otherC(dynamic_cast<const MEDCouplingFieldT<T> *>(other));
    if(!otherC)
      {
        reason="field given in input is not castable in MEDCouplingFieldT !";
        return false;
      }
    if(!MEDCouplingField::isEqualIfNotWhy(other,meshPrec,valsPrec,reason))
      return false;
    if(!_time_discr->isEqualIfNotWhy(otherC->_time_discr,T(valsPrec),reason))
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
  bool MEDCouplingFieldT<T>::isEqualWithoutConsideringStr(const MEDCouplingField *other, double meshPrec, double valsPrec) const
  {
    const MEDCouplingFieldT<T> *otherC(dynamic_cast<const MEDCouplingFieldT<T> *>(other));
    if(!otherC)
      return false;
    if(!MEDCouplingField::isEqualWithoutConsideringStr(other,meshPrec,valsPrec))
      return false;
    if(!_time_discr->isEqualWithoutConsideringStr(otherC->_time_discr,T(valsPrec)))
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
}

#endif
