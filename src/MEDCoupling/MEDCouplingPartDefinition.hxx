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

#ifndef __PARAMEDMEM_MEDCOUPLINGPARTDEFINITION_HXX__
#define __PARAMEDMEM_MEDCOUPLINGPARTDEFINITION_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MCAuto.hxx"

namespace MEDCoupling
{
  class PartDefinition : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT static PartDefinition *New(int start, int stop, int step);
    MEDCOUPLING_EXPORT static PartDefinition *New(DataArrayInt *listOfIds);
    MEDCOUPLING_EXPORT static PartDefinition *Unserialize(std::vector<int>& tinyInt, std::vector< MCAuto<DataArrayInt> >& bigArraysI);
    MEDCOUPLING_EXPORT virtual bool isEqual(const PartDefinition *other, std::string& what) const = 0;
    MEDCOUPLING_EXPORT virtual PartDefinition *deepCopy() const = 0;
    MEDCOUPLING_EXPORT virtual DataArrayInt *toDAI() const = 0;
    MEDCOUPLING_EXPORT virtual int getNumberOfElems() const = 0;
    MEDCOUPLING_EXPORT virtual PartDefinition *operator+(const PartDefinition& other) const = 0;
    MEDCOUPLING_EXPORT virtual std::string getRepr() const = 0;
    MEDCOUPLING_EXPORT virtual PartDefinition *composeWith(const PartDefinition *other) const = 0;
    MEDCOUPLING_EXPORT virtual void checkConsistencyLight() const = 0;
    MEDCOUPLING_EXPORT virtual PartDefinition *tryToSimplify() const = 0;
    MEDCOUPLING_EXPORT virtual void serialize(std::vector<int>& tinyInt, std::vector< MCAuto<DataArrayInt> >& bigArraysI) const = 0;
  protected:
    virtual ~PartDefinition();
  };

  class SlicePartDefinition;

  class DataArrayPartDefinition : public PartDefinition
  {
  public:
    MEDCOUPLING_EXPORT static DataArrayPartDefinition *New(DataArrayInt *listOfIds);
    MEDCOUPLING_EXPORT bool isEqual(const PartDefinition *other, std::string& what) const;
    MEDCOUPLING_EXPORT DataArrayPartDefinition *deepCopy() const;
    MEDCOUPLING_EXPORT DataArrayInt *toDAI() const;
    MEDCOUPLING_EXPORT int getNumberOfElems() const;
    MEDCOUPLING_EXPORT PartDefinition *operator+(const PartDefinition& other) const;
    MEDCOUPLING_EXPORT std::string getRepr() const;
    MEDCOUPLING_EXPORT PartDefinition *composeWith(const PartDefinition *other) const;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT PartDefinition *tryToSimplify() const;
    MEDCOUPLING_EXPORT void serialize(std::vector<int>& tinyInt, std::vector< MCAuto<DataArrayInt> >& bigArraysI) const;
  private:
    DataArrayPartDefinition(DataArrayInt *listOfIds);
    void checkInternalArrayOK() const;
    static void CheckInternalArrayOK(const DataArrayInt *listOfIds);
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    DataArrayPartDefinition *add1(const DataArrayPartDefinition *other) const;
    DataArrayPartDefinition *add2(const SlicePartDefinition *other) const;
    virtual ~DataArrayPartDefinition();
  private:
    MCAuto<DataArrayInt> _arr;
  };

  class SlicePartDefinition : public PartDefinition
  {
  public:
    MEDCOUPLING_EXPORT static SlicePartDefinition *New(int start, int stop, int step);
    MEDCOUPLING_EXPORT bool isEqual(const PartDefinition *other, std::string& what) const;
    MEDCOUPLING_EXPORT SlicePartDefinition *deepCopy() const;
    MEDCOUPLING_EXPORT DataArrayInt *toDAI() const;
    MEDCOUPLING_EXPORT int getNumberOfElems() const;
    MEDCOUPLING_EXPORT PartDefinition *operator+(const PartDefinition& other) const;
    MEDCOUPLING_EXPORT std::string getRepr() const;
    MEDCOUPLING_EXPORT PartDefinition *composeWith(const PartDefinition *other) const;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT PartDefinition *tryToSimplify() const;
    MEDCOUPLING_EXPORT void serialize(std::vector<int>& tinyInt, std::vector< MCAuto<DataArrayInt> >& bigArraysI) const;
    //specific method
    MEDCOUPLING_EXPORT int getEffectiveStop() const;
    MEDCOUPLING_EXPORT void getSlice(int& start, int& stop, int& step) const;
  private:
    SlicePartDefinition(int start, int stop, int step);
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    DataArrayPartDefinition *add1(const DataArrayPartDefinition *other) const;
    PartDefinition *add2(const SlicePartDefinition *other) const;
    virtual ~SlicePartDefinition();
  private:
    int _start;
    int _stop;
    int _step;
  };
}

#endif
