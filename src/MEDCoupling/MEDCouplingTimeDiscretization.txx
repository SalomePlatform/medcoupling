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

#ifndef __MEDCOUPLINGTIMEDISCRETIZATION_TXX__
#define __MEDCOUPLINGTIMEDISCRETIZATION_TXX__

#include "MEDCouplingTimeDiscretization.hxx"

#include <cmath>
#include <sstream>

namespace MEDCoupling
{
  template<class T>
  const double MEDCouplingTimeDiscretizationTemplate<T>::TIME_TOLERANCE_DFT=1.e-12;
  
  template<class T>
  void MEDCouplingTimeDiscretizationTemplate<T>::updateTime() const
  {
    if(_array)
      updateTimeWith(*_array);
  }
  
  template<class T>
  void MEDCouplingTimeDiscretizationTemplate<T>::setArray(typename Traits<T>::ArrayType *array, TimeLabel *owner)
  {
    if(array!=_array)
      {
        if(_array)
        _array->decrRef();
        _array=array;
        if(_array)
          _array->incrRef();
        if(owner)
          owner->declareAsNew();
      }
  }
  
  template<class T>
  void MEDCouplingTimeDiscretizationTemplate<T>::copyTinyAttrFrom(const MEDCouplingTimeDiscretizationTemplate<T>& other)
  {
    _time_tolerance=other._time_tolerance;
    _time_unit=other._time_unit;
  }
  
  template<class T>
  void MEDCouplingTimeDiscretizationTemplate<T>::copyTinyStringsFrom(const MEDCouplingTimeDiscretizationTemplate<T>& other)
  {
    _time_unit=other._time_unit;
    if(_array && other._array)
      _array->copyStringInfoFrom(*other._array);
  }

  template<class T>
  std::size_t MEDCouplingTimeDiscretizationTemplate<T>::getHeapMemorySizeWithoutChildren() const
  {
    return _time_unit.capacity();
  }
  
  template<class T>
  void MEDCouplingTimeDiscretizationTemplate<T>::checkConsistencyLight() const
  {
    if(!_array)
      throw INTERP_KERNEL::Exception("Field invalid because no values set !");
    if(_time_tolerance<0.)
      throw INTERP_KERNEL::Exception("time tolerance is expected to be greater than 0. !");
  }
  
  template<class T>
  std::vector<const BigMemoryObject *> MEDCouplingTimeDiscretizationTemplate<T>::getDirectChildrenWithNull() const
  {
    std::vector<const BigMemoryObject *> ret;
    ret.push_back(_array);
    return ret;
  }
  
  template<class T>
  bool MEDCouplingTimeDiscretizationTemplate<T>::areStrictlyCompatible(const MEDCouplingTimeDiscretizationTemplate<T> *other, std::string& reason) const
  {
    std::ostringstream oss; oss.precision(15);
    if(_time_unit!=other->_time_unit)
      {
        oss << "Field discretizations differ : this time unit = \"" << _time_unit << "\" and other time unit = \"" << other->_time_unit << "\" !";
        reason=oss.str();
        return false;
      }
    if(std::fabs(_time_tolerance-other->_time_tolerance)>1.e-16)
      {
        oss << "Field discretizations differ : this time tolerance = \"" << _time_tolerance << "\" and other time tolerance = \"" << other->_time_tolerance << "\" !";
        reason=oss.str();
        return false;
      }
    if(_array==0 && other->_array==0)
      return true;
    if(_array==0 || other->_array==0)
      {
        reason="Field discretizations differ : Only one timediscretization between the two this and other has a DataArrayDouble for values defined";
        return false;
      }
    if(_array->getNumberOfComponents()!=other->_array->getNumberOfComponents())
      return false;
    if(_array->getNumberOfTuples()!=other->_array->getNumberOfTuples())
      return false;
    return true;
  }
  
  template<class T>
  bool MEDCouplingTimeDiscretizationTemplate<T>::areCompatible(const MEDCouplingTimeDiscretizationTemplate<T> *other) const
  {
    if(std::fabs(_time_tolerance-other->_time_tolerance)>1.e-16)
      return false;
    if(_array==0 && other->_array==0)
      return true;
    if(_array==0 || other->_array==0)
      return false;
    if(_array->getNumberOfComponents()!=other->_array->getNumberOfComponents())
      return false;
    return true;
  }
  
  template<class T>
  bool MEDCouplingTimeDiscretizationTemplate<T>::areStrictlyCompatibleForMul(const MEDCouplingTimeDiscretizationTemplate<T> *other) const
  {
    if(std::fabs(_time_tolerance-other->_time_tolerance)>1.e-16)
      return false;
    if(_array==0 && other->_array==0)
      return true;
    if(_array==0 || other->_array==0)
      return false;
    int nbC1(_array->getNumberOfComponents()),nbC2(other->_array->getNumberOfComponents());
    int nbMin(std::min(nbC1,nbC2));
    if(nbC1!=nbC2 && nbMin!=1)
      return false;
    return true;
  }
  
  template<class T>
  bool MEDCouplingTimeDiscretizationTemplate<T>::areStrictlyCompatibleForDiv(const MEDCouplingTimeDiscretizationTemplate<T> *other) const
  {
    if(std::fabs(_time_tolerance-other->_time_tolerance)>1.e-16)
      return false;
    if(_array==0 && other->_array==0)
      return true;
    if(_array==0 || other->_array==0)
      return false;
    int nbC1(_array->getNumberOfComponents()),nbC2(other->_array->getNumberOfComponents());
    if(nbC1!=nbC2 && nbC2!=1)
      return false;
    return true;
  }

  template<class T>
  MEDCouplingTimeDiscretizationTemplate<T>::MEDCouplingTimeDiscretizationTemplate():_time_tolerance(TIME_TOLERANCE_DFT),_array(0)
  {
  }

  template<class T>
  MEDCouplingTimeDiscretizationTemplate<T>::MEDCouplingTimeDiscretizationTemplate(const MEDCouplingTimeDiscretizationTemplate<T>& other, bool deepCopy):_time_unit(other._time_unit),_time_tolerance(other._time_tolerance)
  {
    if(other._array)
      _array=other._array->performCopyOrIncrRef(deepCopy);
    else
      _array=0;
  }
  
  template<class T>
  MEDCouplingTimeDiscretizationTemplate<T>::~MEDCouplingTimeDiscretizationTemplate()
  {
    if(_array)
      _array->decrRef();
  }
  
  template<class T>
  void MEDCouplingTimeDiscretizationTemplate<T>::setEndArray(typename Traits<T>::ArrayType *array, TimeLabel *owner)
  {
    throw INTERP_KERNEL::Exception("setEndArray not available for this type of time discretization !");
  }
  
  template<class T>
  void MEDCouplingTimeDiscretizationTemplate<T>::setArrays(const std::vector<typename Traits<T>::ArrayType *>& arrays, TimeLabel *owner)
  {
    if(arrays.size()!=1)
      throw INTERP_KERNEL::Exception("MEDCouplingTimeDiscretization::setArrays : number of arrays must be one.");
    setArray(arrays.back(),owner);
  }
  
  template<class T>
  const typename Traits<T>::ArrayType *MEDCouplingTimeDiscretizationTemplate<T>::getEndArray() const
  {
    throw INTERP_KERNEL::Exception("getEndArray not available for this type of time discretization !");
  }
  
  template<class T>
  typename Traits<T>::ArrayType *MEDCouplingTimeDiscretizationTemplate<T>::getEndArray()
  {
    throw INTERP_KERNEL::Exception("getEndArray not available for this type of time discretization !");
  }
  
  template<class T>
  void MEDCouplingTimeDiscretizationTemplate<T>::getArrays(std::vector<typename Traits<T>::ArrayType *>& arrays) const
  {
    arrays.resize(1);
    arrays[0]=_array;
  }
}

#endif
