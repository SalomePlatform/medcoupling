// Copyright (C) 2020-2023  CEA/DEN, EDF R&D
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

#include "MEDCouplingSkyLineArray.hxx"
#include "ProcessorGroup.hxx"
#include "MEDCouplingMemArray.hxx"

#include <string>
#include <vector>

namespace MEDCoupling
{
  /*!
   * Parallel representation of a SkyLineArray
   *
   * This class is very specific to the requirement of parallel code computations.
   */
  class ParaSkyLineArray : public RefCountObject
  {
  public:
    static ParaSkyLineArray *New(MEDCouplingSkyLineArray *ska, DataArrayIdType *globalIds);
    MCAuto<ParaSkyLineArray> equiRedistribute(mcIdType nbOfEntities) const;
    MEDCouplingSkyLineArray *getSkyLineArray() const;
    DataArrayIdType *getGlobalIdsArray() const;
    virtual ~ParaSkyLineArray() { }
  private:
    ParaSkyLineArray(MEDCouplingSkyLineArray *ska, DataArrayIdType *globalIds);
  protected:
    std::string getClassName() const override { return "ParaSkyLineArray"; }
    std::size_t getHeapMemorySizeWithoutChildren() const override;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const override;
  private:
    MCAuto<MEDCouplingSkyLineArray> _ska;
    MCAuto<DataArrayIdType> _global_ids;
  };
}
