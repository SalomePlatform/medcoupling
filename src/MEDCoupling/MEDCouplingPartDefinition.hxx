// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

namespace ParaMEDMEM
{
  class PartDefinition : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT static PartDefinition *New(int start, int stop, int step);
    MEDCOUPLING_EXPORT static PartDefinition *New(DataArrayInt *listOfIds);
    MEDCOUPLING_EXPORT virtual DataArrayInt *toDAI() const = 0;
    MEDCOUPLING_EXPORT virtual int getNumberOfElems() const = 0;
    MEDCOUPLING_EXPORT virtual PartDefinition *operator+(const PartDefinition& other) const = 0;
    MEDCOUPLING_EXPORT virtual std::string getRepr() const = 0;
  protected:
    virtual ~PartDefinition();
  };

  class SlicePartDefinition;

  class DataArrayPartDefinition : public PartDefinition
  {
  public:
    MEDCOUPLING_EXPORT static DataArrayPartDefinition *New(DataArrayInt *listOfIds);
    MEDCOUPLING_EXPORT DataArrayInt *toDAI() const;
    MEDCOUPLING_EXPORT int getNumberOfElems() const;
    MEDCOUPLING_EXPORT PartDefinition *operator+(const PartDefinition& other) const;
    MEDCOUPLING_EXPORT std::string getRepr() const;
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
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _arr;
  };

  class SlicePartDefinition : public PartDefinition
  {
  public:
    MEDCOUPLING_EXPORT static SlicePartDefinition *New(int start, int stop, int step);
    MEDCOUPLING_EXPORT DataArrayInt *toDAI() const;
    MEDCOUPLING_EXPORT int getNumberOfElems() const;
    MEDCOUPLING_EXPORT PartDefinition *operator+(const PartDefinition& other) const;
    MEDCOUPLING_EXPORT std::string getRepr() const;
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
