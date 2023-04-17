// Copyright (C) 2007-2022  CEA/DEN, EDF R&D
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
// Author : Anida Khizar (CEA/DES)

#ifndef __MEDFILTERENTITY_HXX__
#define __MEDFILTERENTITY_HXX__

#include "MEDCouplingPartDefinition.hxx"
#include "med.h"
#include <memory>

namespace MEDCoupling
{

  /*!
   *
   * This class encapsulates the med_filter object to create the appropriate filter based on the partition given as input:
   * if the partition represents a slice of values, then it's more efficient to have a block filter (to treat a block of data)
   * otherwise, a generic filter is necessary
   *
   */
  class MEDFilterEntity
  {
  public:
    inline MEDFilterEntity();
    ~MEDFilterEntity() { if (_filter != nullptr) MEDfilterClose(_filter.get()); }

    inline void fill(med_idt fid, mcIdType nbOfEntity, mcIdType nbOfValuesPerEntity, mcIdType nbOfConstituentPerValue,
                     const med_int constituentSelect, const med_switch_mode switchMode, const med_storage_mode storageMode, const char * const profileName,
                     const PartDefinition* pd);
    const med_filter *getPtr() const { return _filter.get(); }

  private:
    std::shared_ptr<med_filter> _filter;
  };

  MEDFilterEntity::MEDFilterEntity() : _filter(std::make_shared<med_filter>())
  {
    med_filter& ref = *_filter.get();

    // gcc < 9.x compilers are not able to assign to the shared_ptr th med_filter structure
#if defined(WIN32)
    ref = MED_FILTER_INIT;
#else
    ref = (med_filter)MED_FILTER_INIT;
#endif // WIN32
  }

  void MEDFilterEntity::fill(med_idt fid, mcIdType nbOfEntity, mcIdType nbOfValuesPerEntity, mcIdType nbOfConstituentPerValue,
                             const med_int constituentSelect, const med_switch_mode switchMode, const med_storage_mode storageMode, const char * const profileName,
                             const PartDefinition* pd)
  {
    const SlicePartDefinition *spd(dynamic_cast<const SlicePartDefinition *>(pd));
    if(spd)
      {
        //Here, pd contains a slice, so it's more efficient to define a filter of block
        //(which will load contiguous values)
        mcIdType nbOfEltsToLoad = spd->getNumberOfElems();
        mcIdType strt,end,step;
        spd->getSlice(strt,end,step);
        if(strt<0)
          throw INTERP_KERNEL::Exception("MEDFilterEntity::fill : start pos is negative !");
        if(end>nbOfEntity)
          throw INTERP_KERNEL::Exception("MEDFilterEntity::fill : end is after the authorized range !");
        MEDfilterBlockOfEntityCr(fid,ToMedInt(nbOfEntity),ToMedInt(nbOfValuesPerEntity),ToMedInt(nbOfConstituentPerValue),
                                 constituentSelect,switchMode,storageMode,profileName,
                                 /*start*/ToMedInt(strt+1),/*stride*/ToMedInt(step),/*count*/1,/*blocksize*/ToMedInt(nbOfEltsToLoad),
                                 /*lastblocksize=useless because count=1*/0,_filter.get());
        return;
      }
    const DataArrayPartDefinition *dpd(dynamic_cast<const DataArrayPartDefinition *>(pd));
    if(dpd)
      {
        mcIdType nbOfEltsToLoad = dpd->getNumberOfElems();

      //convert to fortran indexing
      std::vector<mcIdType> dpdPlus1;
      MCAuto<DataArrayIdType> partition(pd->toDAI());
      std::copy(partition->begin(), partition->end(), std::back_inserter(dpdPlus1));
      std::for_each(dpdPlus1.begin(), dpdPlus1.end(), [](mcIdType &node){ node+=1; });

        //Here, pd contains a random selection of non-contiguous values:
        //we need to use a more generic filter (less efficient)
        MEDfilterEntityCr(fid,ToMedInt(nbOfEntity),ToMedInt(nbOfValuesPerEntity),ToMedInt(nbOfConstituentPerValue),
                          constituentSelect,switchMode,storageMode,profileName,
                          ToMedInt(nbOfEltsToLoad), dpdPlus1.data(),
                          _filter.get());
        return;
      }
    throw INTERP_KERNEL::Exception("MEDFilterEntity::fill : empty part definition !");
  }

} // namespace

#endif
