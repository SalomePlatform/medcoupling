//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#ifndef __MEDCOUPLINGTIMEDISCRETIZATION_HXX__
#define __MEDCOUPLINGTIMEDISCRETIZATION_HXX__

#include "MEDCoupling.hxx"
#include "RefCountObject.hxx"
#include "InterpKernelException.hxx"

#include <vector>

namespace ParaMEDMEM
{
  class DataArrayDouble;
  class TimeLabel;

  class MEDCOUPLING_EXPORT MEDCouplingTimeDiscretization
  {
  protected:
    MEDCouplingTimeDiscretization();
    MEDCouplingTimeDiscretization(const MEDCouplingTimeDiscretization& other, bool deepCpy);
  public:
    static MEDCouplingTimeDiscretization *New(TypeOfTimeDiscretization type);
    virtual MEDCouplingTimeDiscretization *performCpy(bool deepCpy) const = 0;
    void setTimeTolerance(double val);
    double getTimeTolerance() const { return _time_tolerance; }
    virtual void checkNoTimePresence() const throw(INTERP_KERNEL::Exception) = 0;
    virtual void checkTimePresence(double time) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void setArray(DataArrayDouble *array, TimeLabel *owner);
    virtual void setArrays(const std::vector<DataArrayDouble *>& arrays, TimeLabel *owner) throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getArray() const { return _array; }
    virtual DataArrayDouble *getEndArray() const { return _array; }
    //! Warning contrary to getArray method this method returns an object to deal with.
    virtual DataArrayDouble *getArrayOnTime(double time) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void getArrays(std::vector<DataArrayDouble *>& arrays) const;
    virtual bool isBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    virtual bool isStrictlyBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    double getTime(int& dt, int& it) const throw(INTERP_KERNEL::Exception) { return getStartTime(dt,it); }
    virtual double getStartTime(int& dt, int& it) const throw(INTERP_KERNEL::Exception) = 0;
    virtual double getEndTime(int& dt, int& it) const throw(INTERP_KERNEL::Exception) = 0;
    void setTime(double time, int dt, int it) throw(INTERP_KERNEL::Exception) { setStartTime(time,dt,it); }
    virtual void setStartTime(double time, int dt, int it) throw(INTERP_KERNEL::Exception) = 0;
    virtual void setEndTime(double time, int dt, int it) throw(INTERP_KERNEL::Exception) = 0;
    virtual void getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void getValueOnDiscTime(int eltId, int dt, int it, double *value) const throw(INTERP_KERNEL::Exception) = 0;
    virtual ~MEDCouplingTimeDiscretization();
  protected:
    double _time_tolerance;
    DataArrayDouble *_array;
  protected:
    static const double TIME_TOLERANCE_DFT;
  };

  class MEDCOUPLING_EXPORT MEDCouplingNoTimeLabel : public MEDCouplingTimeDiscretization
  {
  public:
    MEDCouplingNoTimeLabel();
    MEDCouplingNoTimeLabel(const MEDCouplingTimeDiscretization& other, bool deepCpy);
    MEDCouplingTimeDiscretization *performCpy(bool deepCpy) const;
    void checkNoTimePresence() const throw(INTERP_KERNEL::Exception) { }
    void checkTimePresence(double time) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getArrayOnTime(double time) const throw(INTERP_KERNEL::Exception);
    bool isBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    bool isStrictlyBefore(const MEDCouplingTimeDiscretization *other) const throw(INTERP_KERNEL::Exception);
    double getStartTime(int& dt, int& it) const throw(INTERP_KERNEL::Exception);
    double getEndTime(int& dt, int& it) const throw(INTERP_KERNEL::Exception);
    void setStartTime(double time, int dt, int it) throw(INTERP_KERNEL::Exception);
    void setEndTime(double time, int dt, int it) throw(INTERP_KERNEL::Exception);
    void getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception);
    void getValueOnDiscTime(int eltId, int dt, int it, double *value) const throw(INTERP_KERNEL::Exception);
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=NO_TIME;
  private:
    static const char EXCEPTION_MSG[];
  };

  class MEDCOUPLING_EXPORT MEDCouplingWithTimeStep : public MEDCouplingTimeDiscretization
  {
  protected:
    MEDCouplingWithTimeStep(const MEDCouplingWithTimeStep& other, bool deepCpy);
  public:
    MEDCouplingWithTimeStep();
    MEDCouplingTimeDiscretization *performCpy(bool deepCpy) const;
    void checkNoTimePresence() const throw(INTERP_KERNEL::Exception);
    void checkTimePresence(double time) const throw(INTERP_KERNEL::Exception);
    void setStartTime(double time, int dt, int it) throw(INTERP_KERNEL::Exception) { _time=time; _dt=dt; _it=it; }
    void setEndTime(double time, int dt, int it) throw(INTERP_KERNEL::Exception) { _time=time; _dt=dt; _it=it; }
    double getStartTime(int& dt, int& it) const throw(INTERP_KERNEL::Exception) { dt=_dt; it=_it; return _time; }
    double getEndTime(int& dt, int& it) const throw(INTERP_KERNEL::Exception) { dt=_dt; it=_it; return _time; }
    DataArrayDouble *getArrayOnTime(double time) const throw(INTERP_KERNEL::Exception);
    void getValueOnTime(int eltId, double time, double *value) const throw(INTERP_KERNEL::Exception);
    void getValueOnDiscTime(int eltId, int dt, int it, double *value) const throw(INTERP_KERNEL::Exception);
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=ONE_TIME;
  private:
    static const char EXCEPTION_MSG[];
  protected:
    double _time;
    int _dt;
    int _it;
  };

  class MEDCOUPLING_EXPORT MEDCouplingTwoTimeSteps : public MEDCouplingTimeDiscretization
  {
  protected:
    MEDCouplingTwoTimeSteps();
    ~MEDCouplingTwoTimeSteps();
  public:
    void checkNoTimePresence() const throw(INTERP_KERNEL::Exception);
    void checkTimePresence(double time) const throw(INTERP_KERNEL::Exception);
    void getArrays(std::vector<DataArrayDouble *>& arrays) const;
    DataArrayDouble *getEndArray() const { return _end_array; }
    void setStartTime(double time, int dt, int it) throw(INTERP_KERNEL::Exception) { _start_time=time; _start_dt=dt; _start_it=it; }
    void setEndTime(double time, int dt, int it) throw(INTERP_KERNEL::Exception) { _end_time=time; _end_dt=dt; _end_it=it; }
    double getStartTime(int& dt, int& it) const throw(INTERP_KERNEL::Exception) { dt=_start_dt; it=_start_it; return _start_time; }
    double getEndTime(int& dt, int& it) const throw(INTERP_KERNEL::Exception) { dt=_end_dt; it=_end_it; return _end_time; }
  protected:
    double _start_time;
    double _end_time;
    int _start_dt;
    int _end_dt;
    int _start_it;
    int _end_it;
    DataArrayDouble *_end_array;
  };

  class MEDCOUPLING_EXPORT MEDCouplingLinearTime : public MEDCouplingTwoTimeSteps
  {
  public:
    static const TypeOfTimeDiscretization DISCRETIZATION=LINEAR_TIME;
  };
}

#endif
