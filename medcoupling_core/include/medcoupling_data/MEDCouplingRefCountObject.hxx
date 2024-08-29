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

#pragma once

#include "MEDCoupling.hxx"

#include <set>
#include <map>
#include <atomic>
#include <vector>
#include <string>
#include <cstddef>

namespace MEDCoupling
{
  enum class DeallocType
  {
    C_DEALLOC = 2,
    CPP_DEALLOC = 3,
    C_DEALLOC_WITH_OFFSET = 4
  };

  //! The various spatial discretization of a field
  typedef enum
  {
    ON_CELLS = 0,
    ON_NODES = 1,
    ON_GAUSS_PT = 2,
    ON_GAUSS_NE = 3,
    ON_NODES_KR = 4,
    ON_NODES_FE = 5
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
  MEDCOUPLING_EXPORT std::size_t MEDCouplingSizeOfIDs();
  MEDCOUPLING_EXPORT bool MEDCouplingByteOrder();
  MEDCOUPLING_EXPORT const char *MEDCouplingByteOrderStr();
  MEDCOUPLING_EXPORT bool IsCXX11Compiled();
  
  class MEDCOUPLING_EXPORT BigMemoryObject
  {
  public:
    std::size_t getHeapMemorySize() const;
    std::string getHeapMemorySizeStr() const;
    std::vector<const BigMemoryObject *> getDirectChildren() const;
    std::vector<const BigMemoryObject *> getAllTheProgeny() const;
    bool isObjectInTheProgeny(const BigMemoryObject *obj) const;
    static std::size_t GetHeapMemorySizeOfObjs(const std::vector<const BigMemoryObject *>& objs);
    virtual std::string getClassName() const { return "BigMemoryObject"; }
    virtual std::size_t getHeapMemorySizeWithoutChildren() const = 0;
    virtual std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const = 0;
    std::string debugHeapMemorySize() const;
    virtual ~BigMemoryObject();
  private:
    static std::size_t GetHeapMemoryOfSet(std::set<const BigMemoryObject *>& s1, std::set<const BigMemoryObject *>& s2);
  };

  class MEDCOUPLING_EXPORT RefCountObjectOnly
  {
  protected:
    RefCountObjectOnly();
    RefCountObjectOnly(const RefCountObjectOnly& other);
  public:
    bool decrRef() const;
    void incrRef() const;
    int getRCValue() const;
    RefCountObjectOnly& operator=(const RefCountObjectOnly& other);
  protected:
    virtual ~RefCountObjectOnly();
  private:
    mutable std::atomic<int> _cnt;
  };

  class MEDCOUPLING_EXPORT RefCountObject : public RefCountObjectOnly, public BigMemoryObject
  {
  protected:
    RefCountObject();
    RefCountObject(const RefCountObject& other);
    virtual ~RefCountObject();
  };

  class MEDCOUPLING_EXPORT GlobalDict
  {
  public:
    static GlobalDict *GetInstance();
    bool hasKey(const std::string& key) const;
    std::string value(const std::string& key) const;
    std::vector<std::string> keys() const;
    void erase(const std::string& key);
    void clear();
    void setKeyValue(const std::string& key, const std::string& val);
    void setKeyValueForce(const std::string& key, const std::string& val);
    std::string printSelf() const;
  private:
    GlobalDict() { }
  private:
    static GlobalDict *UNIQUE_INSTANCE;
  private:
    std::map<std::string, std::string> _my_map;
  };
}

