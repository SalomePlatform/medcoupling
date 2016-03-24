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
// Author : Anthony Geay (CEA/DEN)

#ifndef __PARAMEDMEM_MEDCOUPLINGMULTIFIELDS_HXX__
#define __PARAMEDMEM_MEDCOUPLINGMULTIFIELDS_HXX__

#include "MEDCouplingRefCountObject.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MCAuto.hxx"

#include "InterpKernelException.hxx"

#include <vector>

namespace MEDCoupling
{
  class MEDCouplingMesh;
  class DataArrayDouble;
  class MEDCouplingFieldDouble;
  class MEDCouplingFieldTemplate;

  class MEDCouplingMultiFields : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingMultiFields *New(const std::vector<MEDCouplingFieldDouble *>& fs);
    MEDCOUPLING_EXPORT static MEDCouplingMultiFields *New();
    MEDCOUPLING_EXPORT MEDCouplingMultiFields *deepCopy() const;
    MEDCOUPLING_EXPORT std::string getName() const;
    MEDCOUPLING_EXPORT std::string getDescription() const;
    MEDCOUPLING_EXPORT std::string getTimeUnit() const;
    MEDCOUPLING_EXPORT double getTimeResolution() const;
    MEDCOUPLING_EXPORT virtual std::string simpleRepr() const;
    MEDCOUPLING_EXPORT virtual std::string advancedRepr() const;
    MEDCOUPLING_EXPORT virtual bool isEqual(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    MEDCOUPLING_EXPORT virtual bool isEqualWithoutConsideringStr(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    MEDCOUPLING_EXPORT const MEDCouplingFieldDouble *getFieldWithId(int id) const;
    MEDCOUPLING_EXPORT std::vector<const MEDCouplingFieldDouble *> getFields() const;
    MEDCOUPLING_EXPORT int getNumberOfFields() const;
    MEDCOUPLING_EXPORT const MEDCouplingFieldDouble *getFieldAtPos(int id) const;
    MEDCOUPLING_EXPORT virtual std::vector<MEDCouplingMesh *> getMeshes() const;
    MEDCOUPLING_EXPORT virtual std::vector<MEDCouplingMesh *> getDifferentMeshes(std::vector<int>& refs) const;
    MEDCOUPLING_EXPORT virtual std::vector<DataArrayDouble *> getArrays() const;
    MEDCOUPLING_EXPORT virtual std::vector<DataArrayDouble *> getDifferentArrays(std::vector< std::vector<int> >& refs) const;
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT void getTinySerializationInformation(std::vector<int>& tinyInfo, std::vector<double>& tinyInfo2, int& nbOfDiffMeshes, int& nbOfDiffArr) const;
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD,
                                                  const std::vector<MEDCouplingFieldTemplate *>& ft, const std::vector<MEDCouplingMesh *>& ms,
                                                  const std::vector<DataArrayDouble *>& das);
    MEDCOUPLING_EXPORT virtual void checkConsistencyLight() const;
  protected:
    MEDCOUPLING_EXPORT MEDCouplingMultiFields(const std::vector<MEDCouplingFieldDouble *>& fs);
    MEDCOUPLING_EXPORT MEDCouplingMultiFields(const MEDCouplingMultiFields& other);
    MEDCOUPLING_EXPORT MEDCouplingMultiFields();
  protected:
    std::vector< MCAuto<MEDCouplingFieldDouble> > _fs;
  };
}

#endif

