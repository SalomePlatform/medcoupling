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
// Author : Anthony Geay (CEA/DEN)

#ifndef __PARAMEDMEM_MEDCOUPLINGDEFINITIONTIME_HXX__
#define __PARAMEDMEM_MEDCOUPLINGDEFINITIONTIME_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "MCAuto.hxx"


#include <ostream>
#include <cstddef>
#include <string>
#include <vector>

namespace MEDCoupling
{
  class MEDCouplingFieldDouble;

  class MEDCouplingDefinitionTimeSlice : public RefCountObject
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingDefinitionTimeSlice *New(const MEDCouplingFieldDouble *f, int meshId, const std::vector<int>& arrId, int fieldId);
    MEDCOUPLING_EXPORT static MEDCouplingDefinitionTimeSlice *New(TypeOfTimeDiscretization type, const std::vector<int>& tiI, const std::vector<double>& tiD);
    MEDCOUPLING_EXPORT int getArrayId() const { return _array_id; }
    MEDCOUPLING_EXPORT virtual MEDCouplingDefinitionTimeSlice *copy() const = 0;
    MEDCOUPLING_EXPORT virtual bool isEqual(const MEDCouplingDefinitionTimeSlice& other, double eps) const;
    MEDCOUPLING_EXPORT virtual void getHotSpotsTime(std::vector<double>& ret) const = 0;
    MEDCOUPLING_EXPORT virtual void getIdsOnTime(double tm, double eps, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const = 0;
    MEDCOUPLING_EXPORT virtual bool isContaining(double tmp, double eps) const = 0;
    MEDCOUPLING_EXPORT virtual int getStartId() const;
    MEDCOUPLING_EXPORT virtual int getEndId() const;
    MEDCOUPLING_EXPORT virtual void appendRepr(std::ostream& stream) const;
    MEDCOUPLING_EXPORT virtual double getStartTime() const = 0;
    MEDCOUPLING_EXPORT virtual double getEndTime() const = 0;
    MEDCOUPLING_EXPORT virtual void getTinySerializationInformation(std::vector<int>& tiI, std::vector<double>& tiD) const = 0;
    MEDCOUPLING_EXPORT virtual TypeOfTimeDiscretization getTimeType() const = 0;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const override;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const override;
    MEDCOUPLING_EXPORT bool isFullyIncludedInMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const;
    MEDCOUPLING_EXPORT bool isOverllapingWithMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const;
    MEDCOUPLING_EXPORT bool isAfterMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const;
    MEDCOUPLING_EXPORT bool isBeforeMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const;
  protected:
    MEDCOUPLING_EXPORT MEDCouplingDefinitionTimeSlice() { }
    MEDCOUPLING_EXPORT MEDCouplingDefinitionTimeSlice(const MEDCouplingFieldDouble *f, int meshId, int arrId, int fieldId);
  protected:
    int _mesh_id;
    int _array_id;
    int _field_id;
  };

  class MEDCouplingDefinitionTimeSliceInst : public MEDCouplingDefinitionTimeSlice
  {
  public:
    static MEDCouplingDefinitionTimeSliceInst *New(const std::vector<int>& tiI, const std::vector<double>& tiD);
    MEDCOUPLING_EXPORT std::string getClassName() const override { return std::string("MEDCouplingDefinitionTimeSliceInst"); }
    MEDCouplingDefinitionTimeSlice *copy() const override;
    bool isEqual(const MEDCouplingDefinitionTimeSlice& other, double eps) const override;
    void getHotSpotsTime(std::vector<double>& ret) const override;
    void getIdsOnTime(double tm, double eps, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const override;
    bool isContaining(double tmp, double eps) const override;
    void appendRepr(std::ostream& stream) const override;
    double getStartTime() const override;
    double getEndTime() const override;
    void getTinySerializationInformation(std::vector<int>& tiI, std::vector<double>& tiD) const override;
    void unserialize(const std::vector<int>& tiI, const std::vector<double>& tiD);
    TypeOfTimeDiscretization getTimeType() const override;
  public:
    MEDCouplingDefinitionTimeSliceInst(const MEDCouplingFieldDouble *f, int meshId, int arrId, int fieldId);
  protected:
    MEDCouplingDefinitionTimeSliceInst() { }
  protected:
    double _instant;
  };

  class MEDCouplingDefinitionTimeSliceCstOnTI : public  MEDCouplingDefinitionTimeSlice
  {
  public:
    static MEDCouplingDefinitionTimeSliceCstOnTI *New(const std::vector<int>& tiI, const std::vector<double>& tiD);
    MEDCOUPLING_EXPORT std::string getClassName() const override { return std::string("MEDCouplingDefinitionTimeSliceCstOnTI"); }
    MEDCouplingDefinitionTimeSlice *copy() const override;
    bool isEqual(const MEDCouplingDefinitionTimeSlice& other, double eps) const override;
    void getHotSpotsTime(std::vector<double>& ret) const override;
    void getIdsOnTime(double tm, double eps, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const override;
    bool isContaining(double tmp, double eps) const override;
    void appendRepr(std::ostream& stream) const override;
    double getStartTime() const override;
    double getEndTime() const override;
    void getTinySerializationInformation(std::vector<int>& tiI, std::vector<double>& tiD) const override;
    void unserialize(const std::vector<int>& tiI, const std::vector<double>& tiD);
    TypeOfTimeDiscretization getTimeType() const override;
  public:
    MEDCouplingDefinitionTimeSliceCstOnTI(const MEDCouplingFieldDouble *f, int meshId, int arrId, int fieldId);
  protected:
    MEDCouplingDefinitionTimeSliceCstOnTI() { }
  protected:
    double _start;
    double _end;
  };

  class MEDCouplingDefinitionTimeSliceLT : public MEDCouplingDefinitionTimeSlice
  {
  public:
    static MEDCouplingDefinitionTimeSliceLT *New(const std::vector<int>& tiI, const std::vector<double>& tiD);
    std::string getClassName() const override { return std::string("MEDCouplingDefinitionTimeSliceLT"); }
    MEDCouplingDefinitionTimeSlice *copy() const override;
    bool isEqual(const MEDCouplingDefinitionTimeSlice& other, double eps) const override;
    void getHotSpotsTime(std::vector<double>& ret) const override;
    void getIdsOnTime(double tm, double eps, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const override;
    bool isContaining(double tmp, double eps) const override;
    void appendRepr(std::ostream& stream) const override;
    double getStartTime() const override;
    double getEndTime() const override;
    int getEndId() const override;
    void getTinySerializationInformation(std::vector<int>& tiI, std::vector<double>& tiD) const override;
    void unserialize(const std::vector<int>& tiI, const std::vector<double>& tiD);
    TypeOfTimeDiscretization getTimeType() const override;
  public:
    MEDCouplingDefinitionTimeSliceLT(const MEDCouplingFieldDouble *f, int meshId, int arrId, int arr2Id, int fieldId);
  protected:
    MEDCouplingDefinitionTimeSliceLT() { }
  protected:
    int _array_id_end;
    double _start;
    double _end;
  };

  class MEDCouplingDefinitionTime
  {
  public:
    MEDCOUPLING_EXPORT MEDCouplingDefinitionTime();
    MEDCOUPLING_EXPORT MEDCouplingDefinitionTime(const std::vector<const MEDCouplingFieldDouble *>& fs, const std::vector<int>& meshRefs, const std::vector<std::vector<int> >& arrRefs);
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT void assign(const MEDCouplingDefinitionTime& other);
    MEDCOUPLING_EXPORT bool isEqual(const MEDCouplingDefinitionTime& other) const;
    MEDCOUPLING_EXPORT double getTimeResolution() const { return _eps; }
    MEDCOUPLING_EXPORT void getIdsOnTimeRight(double tm, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const;
    MEDCOUPLING_EXPORT void getIdsOnTimeLeft(double tm, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const;
    MEDCOUPLING_EXPORT void getIdsOnTime(double tm, std::vector<int>& meshIds, std::vector<int>& arrIds, std::vector<int>& arrIdsInField, std::vector<int>& fieldIds) const;
    MEDCOUPLING_EXPORT std::vector<double> getHotSpotsTime() const;
    MEDCOUPLING_EXPORT void appendRepr(std::ostream& stream) const;
  public:
    MEDCOUPLING_EXPORT void getTinySerializationInformation(std::vector<int>& tinyInfoI, std::vector<double>& tinyInfoD) const;
    MEDCOUPLING_EXPORT void unserialize(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD);
  private:
    double _eps;
    std::vector< MCAuto<MEDCouplingDefinitionTimeSlice> > _slices;
    static const double EPS_DFT;
  };
}

#endif
