// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#ifndef __MEDFILEPARAMETER_HXX__
#define __MEDFILEPARAMETER_HXX__

#include "MEDLoaderDefines.hxx"
#include "MEDFileUtilities.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

namespace ParaMEDMEM
{
  class MEDLOADER_EXPORT MEDFileParameter1TS : public RefCountObject
  {
  public:
    virtual MEDFileParameter1TS *deepCpy() const throw(INTERP_KERNEL::Exception) = 0;
    virtual bool isEqual(const MEDFileParameter1TS *other, double eps, std::string& what) const;
    virtual void simpleRepr2(int bkOffset, std::ostream& oss) const = 0;
    virtual void readValue(med_idt fid, const std::string& name) throw(INTERP_KERNEL::Exception) = 0;
    virtual void writeLL(med_idt fid, const std::string& name, const MEDFileWritable& mw) const throw(INTERP_KERNEL::Exception) = 0;
  public:
    void setIteration(int it) { _iteration=it; }
    int getIteration() const { return _iteration; }
    void setOrder(int order) { _order=order; }
    int getOrder() const { return _order; }
    void setTimeValue(double time) { _time=time; }
    void setTime(int dt, int it, double time) { _time=time; _iteration=dt; _order=it; }
    double getTime(int& dt, int& it) { dt=_iteration; it=_order; return _time; }
    double getTimeValue() const { return _time; }
  protected:
    MEDFileParameter1TS(int iteration, int order, double time);
    MEDFileParameter1TS();
  protected:
    int _iteration;
    int _order;
    double _time;    
  };
  
  class MEDLOADER_EXPORT MEDFileParameterDouble1TSWTI : public MEDFileParameter1TS
  {
  public:
    static MEDFileParameterDouble1TSWTI *New(int iteration, int order, double time);
    MEDFileParameter1TS *deepCpy() const throw(INTERP_KERNEL::Exception);
    void setValue(double val) throw(INTERP_KERNEL::Exception) { _arr=val; }
    double getValue() const throw(INTERP_KERNEL::Exception) { return _arr; }
    bool isEqual(const MEDFileParameter1TS *other, double eps, std::string& what) const;
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildren() const;
    void readValue(med_idt fid, const std::string& name) throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
  protected:
    MEDFileParameterDouble1TSWTI();
    MEDFileParameterDouble1TSWTI(int iteration, int order, double time);
    void simpleRepr2(int bkOffset, std::ostream& oss) const;
    void finishLoading(med_idt fid, const std::string& name, int dt, int it, int nbOfSteps) throw(INTERP_KERNEL::Exception);
    void finishLoading(med_idt fid, const std::string& name, int timeStepId) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid, const std::string& name, const MEDFileWritable& mw) const throw(INTERP_KERNEL::Exception);
  protected:
    double _arr;
  };

  class MEDLOADER_EXPORT MEDFileParameterTinyInfo : public MEDFileWritable
  {
  public:
    void setDescription(const char *name) { _desc_name=name; }
    std::string getDescription() const { return _desc_name; }
    void setTimeUnit(const char *unit) { _dt_unit=unit; }
    std::string getTimeUnit() const { return _dt_unit; }
    std::size_t getHeapMemSizeOfStrings() const;
    bool isEqualStrings(const MEDFileParameterTinyInfo& other, std::string& what) const;
  protected:
    void writeLLHeader(med_idt fid, med_parameter_type typ) const throw(INTERP_KERNEL::Exception);
    void mainRepr(int bkOffset, std::ostream& oss) const;
  protected:
    std::string _dt_unit;
    std::string _name;
    std::string _desc_name;
  };

  class MEDLOADER_EXPORT MEDFileParameterDouble1TS : public MEDFileParameterDouble1TSWTI, public MEDFileParameterTinyInfo
  {
  public:
    static MEDFileParameterDouble1TS *New();
    static MEDFileParameterDouble1TS *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileParameterDouble1TS *New(const char *fileName, const char *paramName) throw(INTERP_KERNEL::Exception);
    static MEDFileParameterDouble1TS *New(const char *fileName, const char *paramName, int dt, int it) throw(INTERP_KERNEL::Exception);
    virtual MEDFileParameter1TS *deepCpy() const throw(INTERP_KERNEL::Exception);
    virtual bool isEqual(const MEDFileParameter1TS *other, double eps, std::string& what) const;
    virtual std::string simpleRepr() const;
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildren() const;
    void setName(const char *name) throw(INTERP_KERNEL::Exception) { _name=name; }
    std::string getName() const throw(INTERP_KERNEL::Exception) { return _name; }
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
  private:
    MEDFileParameterDouble1TS();
    MEDFileParameterDouble1TS(const char *fileName) throw(INTERP_KERNEL::Exception);
    MEDFileParameterDouble1TS(const char *fileName, const char *paramName) throw(INTERP_KERNEL::Exception);
    MEDFileParameterDouble1TS(const char *fileName, const char *paramName, int dt, int it) throw(INTERP_KERNEL::Exception);
  };

  class MEDLOADER_EXPORT MEDFileParameterMultiTS : public RefCountObject, public MEDFileParameterTinyInfo
  {
  public:
    static MEDFileParameterMultiTS *New();
    static MEDFileParameterMultiTS *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileParameterMultiTS *New(const char *fileName, const char *paramName) throw(INTERP_KERNEL::Exception);
    std::string getName() const { return _name; }
    void setName(const char *name) { _name=name; }
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDFileParameterMultiTS *deepCpy() const throw(INTERP_KERNEL::Exception);
    bool isEqual(const MEDFileParameterMultiTS *other, double eps, std::string& what) const;
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid, const MEDFileWritable& mw) const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    void appendValue(int dt, int it, double time, double val) throw(INTERP_KERNEL::Exception);
    double getDoubleValue(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    int getPosOfTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    int getPosGivenTime(double time, double eps=1e-8) const throw(INTERP_KERNEL::Exception);
    MEDFileParameter1TS *getTimeStepAtPos(int posId) const throw(INTERP_KERNEL::Exception);
    void eraseTimeStepIds(const int *startIds, const int *endIds) throw(INTERP_KERNEL::Exception);
    std::vector< std::pair<int,int> > getIterations() const throw(INTERP_KERNEL::Exception);
    std::vector< std::pair<int,int> > getTimeSteps(std::vector<double>& ret1) const throw(INTERP_KERNEL::Exception);
    void simpleRepr2(int bkOffset, std::ostream& oss) const;
  protected:
    MEDFileParameterMultiTS();
    MEDFileParameterMultiTS(const MEDFileParameterMultiTS& other, bool deepCopy);
    MEDFileParameterMultiTS(const char *fileName) throw(INTERP_KERNEL::Exception);
    MEDFileParameterMultiTS(const char *fileName, const char *paramName) throw(INTERP_KERNEL::Exception);
    void finishLoading(med_idt fid, med_parameter_type typ, int nbOfSteps) throw(INTERP_KERNEL::Exception);
  protected:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileParameter1TS> > _param_per_ts;
  };

  class MEDLOADER_EXPORT MEDFileParameters : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileParameters *New();
    static MEDFileParameters *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDFileParameters *deepCpy() const throw(INTERP_KERNEL::Exception);
    bool isEqual(const MEDFileParameters *other, double eps, std::string& what) const;
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getParamsNames() const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    void simpleReprWithoutHeader(std::ostream& oss) const;
    void resize(int newSize) throw(INTERP_KERNEL::Exception);
    void pushParam(MEDFileParameterMultiTS *param) throw(INTERP_KERNEL::Exception);
    void setParamAtPos(int i, MEDFileParameterMultiTS *param) throw(INTERP_KERNEL::Exception);
    MEDFileParameterMultiTS *getParamAtPos(int i) const throw(INTERP_KERNEL::Exception);
    MEDFileParameterMultiTS *getParamWithName(const char *paramName) const throw(INTERP_KERNEL::Exception);
    void destroyParamAtPos(int i) throw(INTERP_KERNEL::Exception);
    int getPosFromParamName(const char *paramName) const throw(INTERP_KERNEL::Exception);
    int getNumberOfParams() const throw(INTERP_KERNEL::Exception);
  protected:
    void simpleRepr2(int bkOffset, std::ostream& oss) const;
    MEDFileParameters(const char *fileName) throw(INTERP_KERNEL::Exception);
    MEDFileParameters(const MEDFileParameters& other, bool deepCopy);
    MEDFileParameters();
  protected:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileParameterMultiTS> > _params;
  };
}

#endif
