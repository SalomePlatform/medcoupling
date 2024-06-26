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

namespace INTERP_KERNEL
{
  class Exception
  {
  public:
    Exception(const char* what);
    ~Exception() throw ();
    const char *what() const throw ();
    %extend
    {
      std::string __str__() const
      {
        return std::string(self->what());
      }
    }
  };
}

/*
 * ABN: Install default exception handler: this catches all INTERP_KERNEL::Exception (even when no
 * except declaration was added to the function declaration) and propagates it to the Python level.
 */
%exception {
  try {
    $action
  }
  catch (INTERP_KERNEL::Exception& _e) {
    // Reraise with SWIG_Python_Raise
    SWIG_Python_Raise(SWIG_NewPointerObj((new INTERP_KERNEL::Exception(static_cast< const INTERP_KERNEL::Exception& >(_e))),SWIGTYPE_p_INTERP_KERNEL__Exception,SWIG_POINTER_OWN), "INTERP_KERNEL::Exception", SWIGTYPE_p_INTERP_KERNEL__Exception);
#ifdef MEDCOUPLING_SWIG4_COMPAT
    return nullptr;
#else
    SWIG_fail;
#endif
  }
}

namespace MEDCoupling
{
  class TimeLabel
  {
  public:
    void declareAsNew() const;
    virtual void updateTime() const;
    unsigned int getTimeOfThis() const;
  protected:
    ~TimeLabel();
  };
}

namespace MEDCoupling
{
  enum class DeallocType
    {
      C_DEALLOC = 2,
      CPP_DEALLOC = 3
    };

  const char *MEDCouplingVersionStr();
  int MEDCouplingVersion();
  int MEDCouplingSizeOfVoidStar();
  int MEDCouplingSizeOfIDs();
  bool MEDCouplingByteOrder();
  const char *MEDCouplingByteOrderStr();
  bool IsCXX11Compiled();
  
  class BigMemoryObject
  {
  public:
    std::size_t getHeapMemorySize() const;
    std::string getHeapMemorySizeStr() const;
    bool isObjectInTheProgeny(const BigMemoryObject *obj) const;
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::string debugHeapMemorySize() const;
    std::string getClassName() const;
    virtual ~BigMemoryObject();
    %extend
    {
      virtual PyObject *getDirectChildren() const
      {
        std::vector<const BigMemoryObject *> c(self->getDirectChildren());
        PyObject *ret(PyList_New(c.size()));
        for(std::size_t i=0;i<c.size();i++)
          PyList_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(c[i]),SWIGTYPE_p_MEDCoupling__BigMemoryObject, 0 | 0 ));
        return ret;
      }

      PyObject *getAllTheProgeny() const
      {
        std::vector<const BigMemoryObject *> c(self->getAllTheProgeny());
        PyObject *ret(PyList_New(c.size()));
        for(std::size_t i=0;i<c.size();i++)
          PyList_SetItem(ret,i,SWIG_NewPointerObj(SWIG_as_voidptr(c[i]),SWIGTYPE_p_MEDCoupling__BigMemoryObject, 0 | 0 ));
        return ret;
      }

      static std::size_t GetHeapMemorySizeOfObjs(PyObject *objs)
      {
        std::vector<const BigMemoryObject *> cppObjs;
        convertFromPyObjVectorOfObj<const MEDCoupling::BigMemoryObject *>(objs,SWIGTYPE_p_MEDCoupling__BigMemoryObject,"BigMemoryObject",cppObjs);
        return BigMemoryObject::GetHeapMemorySizeOfObjs(cppObjs);
      }
    }
  };
  
  class RefCountObjectOnly
  {
  public:
    bool decrRef() const;
    void incrRef() const;
    int getRCValue() const;
  protected:
    ~RefCountObjectOnly();
  };

  class RefCountObject : public RefCountObjectOnly, public BigMemoryObject
  {
  protected:
    ~RefCountObject();
  };

  class GlobalDict
  {
  public:
    static GlobalDict *GetInstance();
    bool hasKey(const std::string& key) const;
    std::string value(const std::string& key) const;
    std::vector<std::string> keys() const;
    void erase(const std::string& key);
    void clear();
    void setKeyValue(const std::string& key, const std::string& value);
    void setKeyValueForce(const std::string& key, const std::string& value);
  private:
    GlobalDict();
  public:
    %extend
    {
      std::string __str__() const
      {
        return self->printSelf();
      }
    }
  };
}

%inline
{
  PyObject *MEDCouplingVersionMajMinRel()
  {
    int tmp0=0,tmp1=0,tmp2=0;
    MEDCouplingVersionMajMinRel(tmp0,tmp1,tmp2);
    PyObject *res = PyList_New(3);
    PyList_SetItem(res,0,SWIG_From_int(tmp0));
    PyList_SetItem(res,1,SWIG_From_int(tmp1));
    PyList_SetItem(res,2,SWIG_From_int(tmp2));
    return res;
  }

  bool MEDCouplingHasNumPyBindings()
  {
#ifdef WITH_NUMPY
    return true;
#else
    return false;
#endif
  }

  bool MEDCouplingHasSciPyBindings()
  {
#ifdef WITH_SCIPY
    return true;
#else
    return false;
#endif
  }
  
  bool MEDCouplingUse64BitIDs()
  {
#ifndef MEDCOUPLING_USE_64BIT_IDS
    return false;
#else
    return true;
#endif
  }

  std::string MEDCouplingCompletionScript()
  {
    static const char script[]="import rlcompleter,readline\nreadline.parse_and_bind('tab:complete')";
    std::ostringstream oss; oss << "MEDCouplingCompletionScript : error when trying to activate completion ! readline not present ?\nScript is :\n" << script;
    if(PyRun_SimpleString(script)!=0)
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    return std::string(script);
  }
}

%pythoncode %{
def INTERPKERNELExceptionReduceFunct(a,b):
    ret=InterpKernelException.__new__(a)
    ret.__init__(*b)
    return ret
def INTERPKERNELExceptionReduce(self):
    return INTERPKERNELExceptionReduceFunct,(InterpKernelException,(self.what(),))
%}
