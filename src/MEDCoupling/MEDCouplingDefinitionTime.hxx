//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#ifndef __PARAMEDMEM_MEDCOUPLINGDEFINITIONTIME_HXX__
#define __PARAMEDMEM_MEDCOUPLINGDEFINITIONTIME_HXX__

#include "MEDCouplingRefCountObject.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include "InterpKernelException.hxx"

#include <vector>
#include <sstream>

namespace ParaMEDMEM
{
  class MEDCouplingFieldDouble;

  class MEDCouplingDefinitionTimeSlice : public RefCountObject
  {
  public:
    static MEDCouplingDefinitionTimeSlice *New(const MEDCouplingFieldDouble *f, int meshId, const std::vector<int>& arrId, int fieldId) throw(INTERP_KERNEL::Exception);
    static MEDCouplingDefinitionTimeSlice *New(TypeOfTimeDiscretization type, const std::vector<int>& tiI, const std::vector<double>& tiD) throw(INTERP_KERNEL::Exception);
    int getArrayId() const { return _array_id; }
    virtual MEDCouplingDefinitionTimeSlice *copy() const = 0;
    virtual bool isEqual(const MEDCouplingDefinitionTimeSlice& other, double eps) const;
    virtual void getHotSpotsTime(std::vector<double>& ret) const = 0;
    virtual void getIdsOnTime(double tm, double eps, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const throw(INTERP_KERNEL::Exception) = 0;
    virtual bool isContaining(double tmp, double eps) const = 0;
    virtual int getStartId() const;
    virtual int getEndId() const;
    virtual void appendRepr(std::ostream& stream) const;
    virtual double getStartTime() const = 0;
    virtual double getEndTime() const = 0;
    virtual void getTinySerializationInformation(std::vector<int>& tiI, std::vector<double>& tiD) const = 0;
    virtual TypeOfTimeDiscretization getTimeType() const = 0;
    bool isFullyIncludedInMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const;
    bool isOverllapingWithMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const;
    bool isAfterMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const;
    bool isBeforeMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const;
  protected:
    MEDCouplingDefinitionTimeSlice() { }
    MEDCouplingDefinitionTimeSlice(const MEDCouplingFieldDouble *f, int meshId, int arrId, int fieldId) throw(INTERP_KERNEL::Exception);
  protected:
    int _mesh_id;
    int _array_id;
    int _field_id;
  };

  class MEDCouplingDefinitionTimeSliceInst : public MEDCouplingDefinitionTimeSlice
  {
  public:
    static MEDCouplingDefinitionTimeSliceInst *New(const std::vector<int>& tiI, const std::vector<double>& tiD);
    MEDCouplingDefinitionTimeSlice *copy() const;
    bool isEqual(const MEDCouplingDefinitionTimeSlice& other, double eps) const;
    void getHotSpotsTime(std::vector<double>& ret) const;
    void getIdsOnTime(double tm, double eps, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const throw(INTERP_KERNEL::Exception);
    bool isContaining(double tmp, double eps) const;
    void appendRepr(std::ostream& stream) const;
    double getStartTime() const;
    double getEndTime() const;
    void getTinySerializationInformation(std::vector<int>& tiI, std::vector<double>& tiD) const;
    void unserialize(const std::vector<int>& tiI, const std::vector<double>& tiD);
    TypeOfTimeDiscretization getTimeType() const;
  public:
    MEDCouplingDefinitionTimeSliceInst(const MEDCouplingFieldDouble *f, int meshId, int arrId, int fieldId) throw(INTERP_KERNEL::Exception);
  protected:
    MEDCouplingDefinitionTimeSliceInst() { }
  protected:
    double _instant;
  };

  class MEDCouplingDefinitionTimeSliceCstOnTI : public  MEDCouplingDefinitionTimeSlice
  {
  public:
    static MEDCouplingDefinitionTimeSliceCstOnTI *New(const std::vector<int>& tiI, const std::vector<double>& tiD);
    MEDCouplingDefinitionTimeSlice *copy() const;
    bool isEqual(const MEDCouplingDefinitionTimeSlice& other, double eps) const;
    void getHotSpotsTime(std::vector<double>& ret) const;
    void getIdsOnTime(double tm, double eps, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const throw(INTERP_KERNEL::Exception);
    bool isContaining(double tmp, double eps) const;
    void appendRepr(std::ostream& stream) const;
    double getStartTime() const;
    double getEndTime() const;
    void getTinySerializationInformation(std::vector<int>& tiI, std::vector<double>& tiD) const;
    void unserialize(const std::vector<int>& tiI, const std::vector<double>& tiD);
    TypeOfTimeDiscretization getTimeType() const;
  public:
    MEDCouplingDefinitionTimeSliceCstOnTI(const MEDCouplingFieldDouble *f, int meshId, int arrId, int fieldId) throw(INTERP_KERNEL::Exception);
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
    MEDCouplingDefinitionTimeSlice *copy() const;
    bool isEqual(const MEDCouplingDefinitionTimeSlice& other, double eps) const;
    void getHotSpotsTime(std::vector<double>& ret) const;
    void getIdsOnTime(double tm, double eps, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const throw(INTERP_KERNEL::Exception);
    bool isContaining(double tmp, double eps) const;
    void appendRepr(std::ostream& stream) const;
    double getStartTime() const;
    double getEndTime() const;
    int getEndId() const;
    void getTinySerializationInformation(std::vector<int>& tiI, std::vector<double>& tiD) const;
    void unserialize(const std::vector<int>& tiI, const std::vector<double>& tiD);
    TypeOfTimeDiscretization getTimeType() const;
  public:
    MEDCouplingDefinitionTimeSliceLT(const MEDCouplingFieldDouble *f, int meshId, int arrId, int arr2Id, int fieldId) throw(INTERP_KERNEL::Exception);
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
    MEDCouplingDefinitionTime();
    MEDCouplingDefinitionTime(const std::vector<const MEDCouplingFieldDouble *>& fs, const std::vector<int>& meshRefs, const std::vector<std::vector<int> >& arrRefs) throw(INTERP_KERNEL::Exception);
    void assign(const MEDCouplingDefinitionTime& other);
    bool isEqual(const MEDCouplingDefinitionTime& other) const;
    double getTimeResolution() const { return _eps; }
    void getIdsOnTimeRight(double tm, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const throw(INTERP_KERNEL::Exception);
    void getIdsOnTimeLeft(double tm, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const throw(INTERP_KERNEL::Exception);
    void getIdsOnTime(double tm, std::vector<int>& meshIds, std::vector<int>& arrIds, std::vector<int>& arrIdsInField, std::vector<int>& fieldIds) const throw(INTERP_KERNEL::Exception);
    std::vector<double> getHotSpotsTime() const;
    void appendRepr(std::ostream& stream) const;
  public:
    void getTinySerializationInformation(std::vector<int>& tinyInfoI, std::vector<double>& tinyInfoD) const;
    void unserialize(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD);
  private:
    double _eps;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingDefinitionTimeSlice> > _slices;
    static const double EPS_DFT;
  };
}

#endif
