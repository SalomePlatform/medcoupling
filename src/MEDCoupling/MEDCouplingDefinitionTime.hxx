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
    int getArrayId() const { return _array_id; }
    virtual bool isContaining(double tmp, double eps) const = 0;
    virtual int getStartId() const;
    virtual int getEndId() const;
    virtual void appendRepr(std::ostream& stream) const;
    virtual double getStartTime() const = 0;
    virtual double getEndTime() const = 0;
    bool isFullyIncludedInMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const;
    bool isOverllapingWithMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const;
    bool isAfterMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const;
    bool isBeforeMe(const MEDCouplingDefinitionTimeSlice *other, double eps) const;
  protected:
    MEDCouplingDefinitionTimeSlice(const MEDCouplingFieldDouble *f, int meshId, int arrId, int fieldId) throw(INTERP_KERNEL::Exception);
  protected:
    int _mesh_id;
    int _array_id;
    int _field_id;
  };

  class MEDCouplingDefinitionTimeSliceInst : public MEDCouplingDefinitionTimeSlice
  {
  public:
    bool isContaining(double tmp, double eps) const;
    void appendRepr(std::ostream& stream) const;
    double getStartTime() const;
    double getEndTime() const;
  public:
    MEDCouplingDefinitionTimeSliceInst(const MEDCouplingFieldDouble *f, int meshId, int arrId, int fieldId) throw(INTERP_KERNEL::Exception);
  protected:
    double _instant;
  };

  class MEDCouplingDefinitionTimeSliceCstOnTI : public  MEDCouplingDefinitionTimeSlice
  {
  public:
    bool isContaining(double tmp, double eps) const;
    void appendRepr(std::ostream& stream) const;
    double getStartTime() const;
    double getEndTime() const;
  public:
    MEDCouplingDefinitionTimeSliceCstOnTI(const MEDCouplingFieldDouble *f, int meshId, int arrId, int fieldId) throw(INTERP_KERNEL::Exception);
  protected:
    double _start;
    double _end;
  };
  

  class MEDCouplingDefinitionTimeSliceLT : public MEDCouplingDefinitionTimeSlice
  {
  public:
    bool isContaining(double tmp, double eps) const;
    void appendRepr(std::ostream& stream) const;
    double getStartTime() const;
    double getEndTime() const;
    int getEndId() const;
  public:
    MEDCouplingDefinitionTimeSliceLT(const MEDCouplingFieldDouble *f, int meshId, int arrId, int arr2Id, int fieldId) throw(INTERP_KERNEL::Exception);
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
    double getTimeResolution() const { return _eps; }
    void getIdsOnTimeRight(double tm, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const throw(INTERP_KERNEL::Exception);
    void getIdsOnTimeLeft(double tm, int& meshId, int& arrId, int& arrIdInField, int& fieldId) const throw(INTERP_KERNEL::Exception);
    void appendRepr(std::ostream& stream) const;
  private:
    void getIdsOnTime(double tm, std::vector<int>& meshIds, std::vector<int>& arrIds, std::vector<int>& arrIdsInField, std::vector<int>& fieldIds) const throw(INTERP_KERNEL::Exception);
  private:
    double _eps;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingDefinitionTimeSlice> > _slices;
  };
}

#endif
