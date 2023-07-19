// Copyright (C) 2020-2023  CEA, EDF
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

#pragma once

#include "MEDCouplingMemArray.hxx"

namespace MEDCoupling
{
  /*!
   * Parallel representation of a DataArray
   *
   * This class is very specific to the requirement of parallel code computations.
   */
  class ParaDataArray : public RefCountObject
  {
  };

  template<class T>
  class ParaDataArrayTemplate : public ParaDataArray
  {
  protected:
    ParaDataArrayTemplate(typename Traits<T>::ArrayType *seqDa);
  protected:
    std::size_t getHeapMemorySizeWithoutChildren() const override;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const override;
    void checkOKOneComponent(const std::string& msg) const;
  protected:
    MCAuto<typename Traits<T>::ArrayType> _seq_da;
  };

  template<class T>
  class ParaDataArrayDiscrete : public ParaDataArrayTemplate<T>
  {
  public:
    DataArrayIdType *buildComplement(T nbOfElems) const;
  protected:
    ParaDataArrayDiscrete(typename Traits<T>::ArrayType *seqDa):ParaDataArrayTemplate<T>(seqDa) { }
  };

  class ParaDataArrayInt32 : public ParaDataArrayDiscrete<Int32>
  {
  public:
    static ParaDataArrayInt32 *New(DataArrayInt32 *seqDa);
  private:
    ParaDataArrayInt32(DataArrayInt32 *seqDa):ParaDataArrayDiscrete<Int32>(seqDa) { }
    std::string getClassName() const override { return "ParaDataArrayInt32"; }
  };

  class ParaDataArrayInt64 : public ParaDataArrayDiscrete<Int64>
  {
  public:
    static ParaDataArrayInt64 *New(DataArrayInt64 *seqDa);
  private:
    ParaDataArrayInt64(DataArrayInt64 *seqDa):ParaDataArrayDiscrete<Int64>(seqDa) { }
    std::string getClassName() const override { return "ParaDataArrayInt64"; }
  };
  
  #ifndef MEDCOUPLING_USE_64BIT_IDS
  using ParaDataArrayIdType = ParaDataArrayInt32;
  #else
  using ParaDataArrayIdType = ParaDataArrayInt64;
  #endif
}
