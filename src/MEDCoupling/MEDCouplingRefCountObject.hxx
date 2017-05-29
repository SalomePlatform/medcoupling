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

#ifndef __PARAMEDMEM_MEDCOUPLINGREFCOUNTOBJECT_HXX__
#define __PARAMEDMEM_MEDCOUPLINGREFCOUNTOBJECT_HXX__

#include "MEDCoupling.hxx"

#include <set>
#include <map>
#include <vector>
#include <string>
#include <cstddef>

namespace MEDCoupling
{
  typedef enum
  {
    C_DEALLOC = 2,
    CPP_DEALLOC = 3
  } DeallocType;

  //! The various spatial discretization of a field
  typedef enum
  {
    ON_CELLS = 0,
    ON_NODES = 1,
    ON_GAUSS_PT = 2,
    ON_GAUSS_NE = 3,
    ON_NODES_KR = 4
  } TypeOfField;

  //! The various temporal discretization of a field
  typedef enum
  {
    NO_TIME = 4,
    ONE_TIME = 5,
    LINEAR_TIME = 6,
    CONST_ON_TIME_INTERVAL = 7
  } TypeOfTimeDiscretization;

  typedef bool (*FunctionToEvaluate)(const double *pos, double *res);

  MEDCOUPLING_EXPORT const char *MEDCouplingVersionStr();
  MEDCOUPLING_EXPORT int MEDCouplingVersion();
  MEDCOUPLING_EXPORT void MEDCouplingVersionMajMinRel(int& maj, int& minor, int& releas);
  MEDCOUPLING_EXPORT int MEDCouplingSizeOfVoidStar();
  MEDCOUPLING_EXPORT bool MEDCouplingByteOrder();
  MEDCOUPLING_EXPORT const char *MEDCouplingByteOrderStr();

  class BigMemoryObject
  {
  public:
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySize() const;
    MEDCOUPLING_EXPORT std::string getHeapMemorySizeStr() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getAllTheProgeny() const;
    MEDCOUPLING_EXPORT bool isObjectInTheProgeny(const BigMemoryObject *obj) const;
    MEDCOUPLING_EXPORT static std::size_t GetHeapMemorySizeOfObjs(const std::vector<const BigMemoryObject *>& objs);
    MEDCOUPLING_EXPORT virtual std::size_t getHeapMemorySizeWithoutChildren() const = 0;
    MEDCOUPLING_EXPORT virtual std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const = 0;
    MEDCOUPLING_EXPORT virtual ~BigMemoryObject();
  private:
    static std::size_t GetHeapMemoryOfSet(std::set<const BigMemoryObject *>& s1, std::set<const BigMemoryObject *>& s2);
  };

  class RefCountObjectOnly
  {
  protected:
    MEDCOUPLING_EXPORT RefCountObjectOnly();
    MEDCOUPLING_EXPORT RefCountObjectOnly(const RefCountObjectOnly& other);
  public:
    MEDCOUPLING_EXPORT bool decrRef() const;
    MEDCOUPLING_EXPORT void incrRef() const;
    MEDCOUPLING_EXPORT int getRCValue() const;
    MEDCOUPLING_EXPORT RefCountObjectOnly& operator=(const RefCountObjectOnly& other);
  protected:
    virtual ~RefCountObjectOnly();
  private:
    mutable int _cnt;
  };

  class RefCountObject : public RefCountObjectOnly, public BigMemoryObject
  {
  protected:
    MEDCOUPLING_EXPORT RefCountObject();
    MEDCOUPLING_EXPORT RefCountObject(const RefCountObject& other);
    MEDCOUPLING_EXPORT virtual ~RefCountObject();
  };

  class GlobalDict
  {
  public:
    MEDCOUPLING_EXPORT static GlobalDict *GetInstance();
    MEDCOUPLING_EXPORT bool hasKey(const std::string& key) const;
    MEDCOUPLING_EXPORT std::string value(const std::string& key) const;
    MEDCOUPLING_EXPORT std::vector<std::string> keys() const;
    MEDCOUPLING_EXPORT void erase(const std::string& key);
    MEDCOUPLING_EXPORT void clear();
    MEDCOUPLING_EXPORT void setKeyValue(const std::string& key, const std::string& val);
    MEDCOUPLING_EXPORT void setKeyValueForce(const std::string& key, const std::string& val);
    MEDCOUPLING_EXPORT std::string printSelf() const;
  private:
    GlobalDict() { }
  private:
    static GlobalDict *UNIQUE_INSTANCE;
  private:
    std::map<std::string, std::string> _my_map;
  };
}

#endif
