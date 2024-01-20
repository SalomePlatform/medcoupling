// Copyright (C) 2007-2024  CEA, EDF
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

#include "MEDCoupling.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MCAuto.hxx"

namespace MEDCoupling
{
  class PartDefinition : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT static PartDefinition *New(mcIdType start, mcIdType stop, mcIdType step);
    MEDCOUPLING_EXPORT static PartDefinition *New(DataArrayIdType *listOfIds);
    MEDCOUPLING_EXPORT static PartDefinition *Unserialize(std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >& bigArraysI);
    MEDCOUPLING_EXPORT virtual bool isEqual(const PartDefinition *other, std::string& what) const = 0;
    MEDCOUPLING_EXPORT virtual PartDefinition *deepCopy() const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayIdType *toDAI() const = 0;
    MEDCOUPLING_EXPORT virtual mcIdType getNumberOfElems() const = 0;
    MEDCOUPLING_EXPORT virtual PartDefinition *operator+(const PartDefinition& other) const = 0;
    MEDCOUPLING_EXPORT virtual std::string getRepr() const = 0;
    MEDCOUPLING_EXPORT virtual PartDefinition *composeWith(const PartDefinition *other) const = 0;
    MEDCOUPLING_EXPORT virtual void checkConsistencyLight() const = 0;
    MEDCOUPLING_EXPORT virtual PartDefinition *tryToSimplify() const = 0;
    MEDCOUPLING_EXPORT virtual void serialize(std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >& bigArraysI) const = 0;
  protected:
    virtual ~PartDefinition();
  };

  class SlicePartDefinition;

  class DataArrayPartDefinition : public PartDefinition
  {
  public:
    MEDCOUPLING_EXPORT static DataArrayPartDefinition *New(DataArrayIdType *listOfIds);
    MEDCOUPLING_EXPORT std::string getClassName() const override { return std::string("DataArrayPartDefinition"); }
    MEDCOUPLING_EXPORT bool isEqual(const PartDefinition *other, std::string& what) const;
    MEDCOUPLING_EXPORT DataArrayPartDefinition *deepCopy() const;
    MEDCOUPLING_EXPORT DataArrayIdType *toDAI() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfElems() const;
    MEDCOUPLING_EXPORT PartDefinition *operator+(const PartDefinition& other) const;
    MEDCOUPLING_EXPORT std::string getRepr() const;
    MEDCOUPLING_EXPORT PartDefinition *composeWith(const PartDefinition *other) const;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT PartDefinition *tryToSimplify() const;
    MEDCOUPLING_EXPORT void serialize(std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >& bigArraysI) const;
  private:
    DataArrayPartDefinition(DataArrayIdType *listOfIds);
    void checkInternalArrayOK() const;
    static void CheckInternalArrayOK(const DataArrayIdType *listOfIds);
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    DataArrayPartDefinition *add1(const DataArrayPartDefinition *other) const;
    DataArrayPartDefinition *add2(const SlicePartDefinition *other) const;
    virtual ~DataArrayPartDefinition();
  private:
    MCAuto<DataArrayIdType> _arr;
  };

  class SlicePartDefinition : public PartDefinition
  {
  public:
    MEDCOUPLING_EXPORT static SlicePartDefinition *New(mcIdType start, mcIdType stop, mcIdType step);
    MEDCOUPLING_EXPORT std::string getClassName() const override { return std::string("SlicePartDefinition"); }
    MEDCOUPLING_EXPORT bool isEqual(const PartDefinition *other, std::string& what) const;
    MEDCOUPLING_EXPORT SlicePartDefinition *deepCopy() const;
    MEDCOUPLING_EXPORT DataArrayIdType *toDAI() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfElems() const;
    MEDCOUPLING_EXPORT PartDefinition *operator+(const PartDefinition& other) const;
    MEDCOUPLING_EXPORT std::string getRepr() const;
    MEDCOUPLING_EXPORT PartDefinition *composeWith(const PartDefinition *other) const;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT PartDefinition *tryToSimplify() const;
    MEDCOUPLING_EXPORT void serialize(std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >& bigArraysI) const;
    //specific method
    MEDCOUPLING_EXPORT mcIdType getEffectiveStop() const;
    MEDCOUPLING_EXPORT void getSlice(mcIdType& start, mcIdType& stop, mcIdType& step) const;
  private:
    SlicePartDefinition(mcIdType start, mcIdType stop, mcIdType step);
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    DataArrayPartDefinition *add1(const DataArrayPartDefinition *other) const;
    PartDefinition *add2(const SlicePartDefinition *other) const;
    virtual ~SlicePartDefinition();
  private:
    mcIdType _start;
    mcIdType _stop;
    mcIdType _step;
  };
}
