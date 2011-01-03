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

#ifdef WITH_NUMPY2
#include <numpy/arrayobject.h>
#endif

static PyObject* convertMesh(ParaMEDMEM::MEDCouplingMesh* mesh, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(dynamic_cast<ParaMEDMEM::MEDCouplingUMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,owner);
  if(dynamic_cast<ParaMEDMEM::MEDCouplingExtrudedMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDCouplingExtrudedMesh,owner);
  if(dynamic_cast<ParaMEDMEM::MEDCouplingCMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDCouplingCMesh,owner);
  if(!ret)
    {
      const char msg[]="Not recognized type of mesh on downcast !";
      PyErr_SetString(PyExc_TypeError,msg);
      throw INTERP_KERNEL::Exception(msg);
    }
  return ret;
}

static PyObject *convertIntArrToPyList(const int *ptr, int size) throw(INTERP_KERNEL::Exception)
{
#ifndef WITH_NUMPY2
  PyObject *ret=PyList_New(size);
  for(int i=0;i<size;i++)
    PyList_SetItem(ret,i,PyInt_FromLong(ptr[i]));
  return ret;
#else
  npy_intp dim = (npy_intp) size;
  int *tmp=new int[size];
  std::copy(ptr,ptr+size,tmp);
  return PyArray_SimpleNewFromData(1,&dim,NPY_INT,const_cast<int *>(tmp));
#endif
}

static PyObject *convertIntArrToPyList2(const std::vector<int>& v) throw(INTERP_KERNEL::Exception)
{
#ifndef WITH_NUMPY2
  int size=v.size();
  PyObject *ret=PyList_New(size);
  for(int i=0;i<size;i++)
    PyList_SetItem(ret,i,PyInt_FromLong(v[i]));
  return ret;
#else
  npy_intp dim = (npy_intp) v.size();
  int *tmp=new int[v.size()];
  std::copy(v.begin(),v.end(),tmp);
  return PyArray_SimpleNewFromData(1,&dim,NPY_INT,tmp);
#endif
}

static PyObject *convertIntArrToPyList3(const std::set<int>& v) throw(INTERP_KERNEL::Exception)
{
  int size=v.size();
  PyObject *ret=PyList_New(size);
  std::set<int>::const_iterator it=v.begin();
  for(int i=0;i<size;i++,it++)
    PyList_SetItem(ret,i,PyInt_FromLong(*it));
  return ret;
}

static PyObject *convertIntArrToPyListOfTuple(const int *vals, int nbOfComp, int nbOfTuples) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=PyList_New(nbOfTuples);
  for(int i=0;i<nbOfTuples;i++)
    {
      PyObject *t=PyTuple_New(nbOfComp);
      for(int j=0;j<nbOfComp;j++)
        PyTuple_SetItem(t,j,PyInt_FromLong(vals[i*nbOfComp+j]));
      PyList_SetItem(ret,i,t);
    }
  return ret;
}

static int *convertPyToNewIntArr2(PyObject *pyLi, int *size) throw(INTERP_KERNEL::Exception)
{
  if(PyList_Check(pyLi))
    {
      *size=PyList_Size(pyLi);
      int *tmp=new int[*size];
      for(int i=0;i<*size;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(PyInt_Check(o))
            {
              int val=(int)PyInt_AS_LONG(o);
              tmp[i]=val;
            }
          else
            {
              delete [] tmp;
              const char msg[]="list must contain integers only";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
        }
      return tmp;
    }
  else if(PyTuple_Check(pyLi))
    {
      *size=PyTuple_Size(pyLi);
      int *tmp=new int[*size];
      for(int i=0;i<*size;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          if(PyInt_Check(o))
            {
              int val=(int)PyInt_AS_LONG(o);
              tmp[i]=val;
            }
          else
            {
              delete [] tmp;
              const char msg[]="tuple must contain integers only";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
        }
      return tmp;
    }
  else
    {
#ifndef WITH_NUMPY2
      const char msg[]="convertPyToNewIntArr2 : not a list";
      PyErr_SetString(PyExc_TypeError,msg);
      throw INTERP_KERNEL::Exception(msg);
#else
      if(PyArray_Check(pyLi))
        {
          npy_intp mySize = PyArray_SIZE(pyLi);
          int *ret=(int *)PyArray_BYTES(pyLi);
          *size=mySize;
          return ret;
        }
      else
        {
          const char msg[]="convertPyToNewIntArr2 : not a list nor PyArray";
          PyErr_SetString(PyExc_TypeError,msg);
          throw INTERP_KERNEL::Exception(msg);
        }
#endif
    }
}

static void convertPyToNewIntArr3(PyObject *pyLi, std::vector<int>& arr) throw(INTERP_KERNEL::Exception)
{
  if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      arr.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(PyInt_Check(o))
            {
              int val=(int)PyInt_AS_LONG(o);
              arr[i]=val;
            }
          else
            {
              const char msg[]="list must contain integers only";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
        }
    }
  else if(PyTuple_Check(pyLi))
    {
      int size=PyTuple_Size(pyLi);
      arr.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          if(PyInt_Check(o))
            {
              int val=(int)PyInt_AS_LONG(o);
              arr[i]=val;
            }
          else
            {
              const char msg[]="tuple must contain integers only";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
        }
    }
  else
    {
#ifndef WITH_NUMPY2
      const char msg[]="convertPyToNewIntArr3 : not a list";
      PyErr_SetString(PyExc_TypeError,msg);
      throw INTERP_KERNEL::Exception(msg);
#else
      if(PyArray_Check(pyLi))
        {
          npy_intp mySize = PyArray_SIZE(pyLi);
          int *ret=(int *)PyArray_BYTES(pyLi);
          arr.resize(mySize);
          std::copy(ret,ret+mySize,arr.begin());
          return ;
        }
      else
        {
          const char msg[]="convertPyToNewIntArr3 : not a list nor PyArray";
          PyErr_SetString(PyExc_TypeError,msg);
          throw INTERP_KERNEL::Exception(msg);
        }
#endif
    }
}

static void fillArrayWithPyListInt(PyObject *pyLi, int *arrToFill, int sizeOfArray, int dftVal) throw(INTERP_KERNEL::Exception)
{
  if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(PyInt_Check(o))
            {
              int val=(int)PyInt_AS_LONG(o);
              if(i<sizeOfArray)
                arrToFill[i]=val;
            }
          else
            {
              const char msg[]="list must contain integers only";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
        }
      for(int i=size;i<sizeOfArray;i++)
        arrToFill[i]=dftVal;
      return;
      
    }
  else if(PyTuple_Check(pyLi))
    {
      int size=PyTuple_Size(pyLi);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          if(PyInt_Check(o))
            {
              int val=(int)PyInt_AS_LONG(o);
              if(i<sizeOfArray)
                arrToFill[i]=val;
            }
          else
            {
              const char msg[]="tuple must contain integers only";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
        }
      for(int i=size;i<sizeOfArray;i++)
        arrToFill[i]=dftVal;
      return;
    }
  else
    {
      const char msg[]="fillArrayWithPyListInt : not a list";
      PyErr_SetString(PyExc_TypeError,msg);
      throw INTERP_KERNEL::Exception(msg);
    }
}

static PyObject *convertDblArrToPyList(const double *ptr, int size) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=PyList_New(size);
  for(int i=0;i<size;i++)
    PyList_SetItem(ret,i,PyFloat_FromDouble(ptr[i]));
  return ret;
}

static PyObject *convertDblArrToPyList2(const std::vector<double>& v) throw(INTERP_KERNEL::Exception)
{
  int size=v.size();
  PyObject *ret=PyList_New(size);
  for(int i=0;i<size;i++)
    PyList_SetItem(ret,i,PyFloat_FromDouble(v[i]));
  return ret;
}

static PyObject *convertDblArrToPyListOfTuple(const double *vals, int nbOfComp, int nbOfTuples) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=PyList_New(nbOfTuples);
  for(int i=0;i<nbOfTuples;i++)
    {
      PyObject *t=PyTuple_New(nbOfComp);
      for(int j=0;j<nbOfComp;j++)
        PyTuple_SetItem(t,j,PyFloat_FromDouble(vals[i*nbOfComp+j]));
      PyList_SetItem(ret,i,t);
    }
  return ret;
}

static double *convertPyToNewDblArr2(PyObject *pyLi, int *size) throw(INTERP_KERNEL::Exception)
{
  if(PyList_Check(pyLi))
    {
      *size=PyList_Size(pyLi);
      double *tmp=new double[*size];
      for(int i=0;i<*size;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(PyFloat_Check(o))
            {
              double val=PyFloat_AS_DOUBLE(o);
              tmp[i]=val;
            }
          else if(PyInt_Check(o))
            {
              long val0=PyInt_AS_LONG(o);
              double val=val0;
              tmp[i]=val;
            }
          else
            {
              delete [] tmp;
              const char msg[]="convertPyToNewDblArr2 : list must contain floats/integers only";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
        }
      return tmp;
    }
  else if(PyTuple_Check(pyLi))
    {
      *size=PyTuple_Size(pyLi);
      double *tmp=new double[*size];
      for(int i=0;i<*size;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          if(PyFloat_Check(o))
            {
              double val=PyFloat_AS_DOUBLE(o);
              tmp[i]=val;
            }
          else if(PyInt_Check(o))
            {
              long val0=PyInt_AS_LONG(o);
              double val=val0;
              tmp[i]=val;
            }
          else
            {
              delete [] tmp;
              const char msg[]="convertPyToNewDblArr2 : tuple must contain floats/integers only";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
        }
      return tmp;
    }
  else
    {
      const char msg[]="convertPyToNewDblArr2 : not a list";
      PyErr_SetString(PyExc_TypeError,msg);
      throw INTERP_KERNEL::Exception(msg);
    }
}

static void fillArrayWithPyListDbl(PyObject *pyLi, double *arrToFill, int sizeOfArray, double dftVal) throw(INTERP_KERNEL::Exception)
{
  if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(PyFloat_Check(o))
            {
              double val=PyFloat_AS_DOUBLE(o);
              if(i<sizeOfArray)
                arrToFill[i]=val;
            }
          else if(PyInt_Check(o))
            {
              long val0=PyInt_AS_LONG(o);
              double val=val0;
              if(i<sizeOfArray)
                arrToFill[i]=val;
            }
          else
            {
              const char msg[]="fillArrayWithPyListDbl : list must contain floats/integers only";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
        }
      for(int i=size;i<sizeOfArray;i++)
        arrToFill[i]=dftVal;
      return;
    }
  else if(PyTuple_Check(pyLi))
    {
      int size=PyTuple_Size(pyLi);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          if(PyFloat_Check(o))
            {
              double val=PyFloat_AS_DOUBLE(o);
              arrToFill[i]=val;
            }
          else if(PyInt_Check(o))
            {
              long val0=PyInt_AS_LONG(o);
              double val=val0;
              arrToFill[i]=val;
            }
          else
            {
              const char msg[]="fillArrayWithPyListDbl : tuple must contain floats/integers only";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
        }
      for(int i=size;i<sizeOfArray;i++)
        arrToFill[i]=dftVal;
      return ;
    }
  else
    {
      const char msg[]="convertPyToNewIntArr : not a list";
      PyErr_SetString(PyExc_TypeError,msg);
      throw INTERP_KERNEL::Exception(msg);
    }
}

void convertPyObjToVecUMeshesCst(PyObject *ms, std::vector<const ParaMEDMEM::MEDCouplingUMesh *>& v) throw(INTERP_KERNEL::Exception)
{
  if(PyList_Check(ms))
    {
      int size=PyList_Size(ms);
      v.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *obj=PyList_GetItem(ms,i);
          void *argp;
          int status=SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,0|0);
          if(!SWIG_IsOK(status))
            {
              const char msg[]="list must contain only MEDCouplingUMesh";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
          const ParaMEDMEM::MEDCouplingUMesh *arg=reinterpret_cast< const ParaMEDMEM::MEDCouplingUMesh * >(argp);
          v[i]=arg;
        }
    }
  else
    {
      const char msg[]="convertPyObjToVecUMeshesCst : not a list";
      PyErr_SetString(PyExc_TypeError,msg);
      throw INTERP_KERNEL::Exception(msg);
    }
}

void convertPyObjToVecDataArrayDblCst(PyObject *ms, std::vector<const ParaMEDMEM::DataArrayDouble *>& v) throw(INTERP_KERNEL::Exception)
{
  if(PyList_Check(ms))
    {
      int size=PyList_Size(ms);
      v.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *obj=PyList_GetItem(ms,i);
          void *argp;
          int status=SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble,0|0);
          if(!SWIG_IsOK(status))
            {
              const char msg[]="list must contain only DataArrayDouble";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
          const ParaMEDMEM::DataArrayDouble *arg=reinterpret_cast< const ParaMEDMEM::DataArrayDouble * >(argp);
          v[i]=arg;
        }
    }
  else
    {
      const char msg[]="convertPyObjToVecDataArrayDblCst : not a list";
      PyErr_SetString(PyExc_TypeError,msg);
      throw INTERP_KERNEL::Exception(msg);
    }
}

void convertPyObjToVecFieldDblCst(PyObject *ms, std::vector<const ParaMEDMEM::MEDCouplingFieldDouble *>& v) throw(INTERP_KERNEL::Exception)
{
  if(PyList_Check(ms))
    {
      int size=PyList_Size(ms);
      v.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *obj=PyList_GetItem(ms,i);
          void *argp;
          int status=SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,0|0);
          if(!SWIG_IsOK(status))
            {
              const char msg[]="list must contain only MEDCouplingFieldDouble";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
          const ParaMEDMEM::MEDCouplingFieldDouble *arg=reinterpret_cast< const ParaMEDMEM::MEDCouplingFieldDouble * >(argp);
          v[i]=arg;
        }
    }
  else
    {
      const char msg[]="convertPyObjToVecFieldDblCst : not a list";
      PyErr_SetString(PyExc_TypeError,msg);
      throw INTERP_KERNEL::Exception(msg);
    }
}

void convertPyObjToVecDataArrayIntCst(PyObject *ms, std::vector<const ParaMEDMEM::DataArrayInt *>& v) throw(INTERP_KERNEL::Exception)
{
  if(PyList_Check(ms))
    {
      int size=PyList_Size(ms);
      v.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *obj=PyList_GetItem(ms,i);
          void *argp;
          int status=SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,0|0);
          if(!SWIG_IsOK(status))
            {
              const char msg[]="list must contain only DataArrayInt";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
          ParaMEDMEM::DataArrayInt *arg=reinterpret_cast< ParaMEDMEM::DataArrayInt * >(argp);
          v[i]=arg;
        }
    }
  else
    {
      const char msg[]="convertPyObjToVecDataArrayInt : not a list";
      PyErr_SetString(PyExc_TypeError,msg);
      throw INTERP_KERNEL::Exception(msg);
    }
}
