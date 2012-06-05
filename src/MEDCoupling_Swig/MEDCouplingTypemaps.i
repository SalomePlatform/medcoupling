// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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

static PyObject* convertMultiFields(ParaMEDMEM::MEDCouplingMultiFields *mfs, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(dynamic_cast<ParaMEDMEM::MEDCouplingFieldOverTime *>(mfs))
    ret=SWIG_NewPointerObj((void*)mfs,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldOverTime,owner);
  else
    ret=SWIG_NewPointerObj((void*)mfs,SWIGTYPE_p_ParaMEDMEM__MEDCouplingMultiFields,owner);
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

static void convertPyToVectorPairInt(PyObject *pyLi, std::vector< std::pair<int,int> >& arr) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="list must contain tuples of 2 integers only or tuple must contain tuples of 2 integers only !";
  if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      arr.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(PyTuple_Check(o))
            {
              int sz2=PyTuple_Size(o);
              if(sz2!=2)
                throw INTERP_KERNEL::Exception(msg);
              PyObject *o_0=PyTuple_GetItem(o,0);
              if(!PyInt_Check(o_0))
                throw INTERP_KERNEL::Exception(msg);
              PyObject *o_1=PyTuple_GetItem(o,1);
              if(!PyInt_Check(o_1))
                throw INTERP_KERNEL::Exception(msg);
              arr[i].first=(int)PyInt_AS_LONG(o_0);
              arr[i].second=(int)PyInt_AS_LONG(o_1);
            }
          else
            throw INTERP_KERNEL::Exception(msg);
        }
    }
  else if(PyTuple_Check(pyLi))
    {
      int size=PyTuple_Size(pyLi);
      arr.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          if(PyTuple_Check(o))
            {
              int sz2=PyTuple_Size(o);
              if(sz2!=2)
                throw INTERP_KERNEL::Exception(msg);
              PyObject *o_0=PyTuple_GetItem(o,0);
              if(!PyInt_Check(o_0))
                throw INTERP_KERNEL::Exception(msg);
              PyObject *o_1=PyTuple_GetItem(o,1);
              if(!PyInt_Check(o_1))
                throw INTERP_KERNEL::Exception(msg);
              arr[i].first=(int)PyInt_AS_LONG(o_0);
              arr[i].second=(int)PyInt_AS_LONG(o_1);
            }
          else
            throw INTERP_KERNEL::Exception(msg);
        }
    }
  else
    throw INTERP_KERNEL::Exception(msg);
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

static void fillArrayWithPyListInt(PyObject *pyLi, int *arrToFill, int sizeOfArray, int dftVal, bool chckSize) throw(INTERP_KERNEL::Exception)
{
  if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      if(chckSize)
        if(size!=sizeOfArray)
          {
            std::ostringstream oss; oss << "fillArrayWithPyListInt : List expected to be of size " << sizeOfArray << " but the size is " << size << " !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
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
            throw INTERP_KERNEL::Exception("fillArrayWithPyListInt : List must contain integers only !");
        }
      for(int i=size;i<sizeOfArray;i++)
        arrToFill[i]=dftVal;
      return;
      
    }
  else if(PyTuple_Check(pyLi))
    {
      int size=PyTuple_Size(pyLi);
      if(chckSize)
        if(size!=sizeOfArray)
          {
            std::ostringstream oss; oss << "fillArrayWithPyListInt : Tuple expected to be of size " << sizeOfArray << " but the size is " << size << " !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
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

static void fillArrayWithPyListDbl(PyObject *pyLi, double *arrToFill, int sizeOfArray, double dftVal, bool chckSize) throw(INTERP_KERNEL::Exception)
{
  if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      if(chckSize)
        if(size!=sizeOfArray)
          {
            std::ostringstream oss; oss << "fillArrayWithPyListDbl : List expected to be of size " << sizeOfArray << " but the size is " << size << " !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
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
      if(chckSize)
        if(size!=sizeOfArray)
          {
            std::ostringstream oss; oss << "fillArrayWithPyListDbl : Tuple expected to be of size " << sizeOfArray << " but the size is " << size << " !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
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
              const char msg[]="list must contain only instance of MEDCouplingUMesh";
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

void convertPyObjToVecUMeshes(PyObject *ms, std::vector<ParaMEDMEM::MEDCouplingUMesh *>& v) throw(INTERP_KERNEL::Exception)
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
              const char msg[]="list must contain only instance of MEDCouplingUMesh";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
          ParaMEDMEM::MEDCouplingUMesh *arg=reinterpret_cast< ParaMEDMEM::MEDCouplingUMesh * >(argp);
          v[i]=arg;
        }
    }
  else
    {
      const char msg[]="convertPyObjToVecUMeshes : not a list";
      PyErr_SetString(PyExc_TypeError,msg);
      throw INTERP_KERNEL::Exception(msg);
    }
}

void convertPyObjToVecMeshesCst(PyObject *ms, std::vector<const ParaMEDMEM::MEDCouplingMesh *>& v) throw(INTERP_KERNEL::Exception)
{
  if(PyList_Check(ms))
    {
      int size=PyList_Size(ms);
      v.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *obj=PyList_GetItem(ms,i);
          void *argp;
          int status=SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingMesh,0|0);
          if(!SWIG_IsOK(status))
            {
              const char msg[]="list must contain only instance of MEDCouplingMesh";
              PyErr_SetString(PyExc_TypeError,msg);
              throw INTERP_KERNEL::Exception(msg);
            }
          const ParaMEDMEM::MEDCouplingMesh *arg=reinterpret_cast< const ParaMEDMEM::MEDCouplingMesh * >(argp);
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
              const char msg[]="list must contain only instance of MEDCouplingFieldDouble";
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
              const char msg[]="list must contain only instance of DataArrayInt";
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

/*!
 * if python int -> cpp int sw=1
 * if python list[int] -> cpp vector<int> sw=2
 * if python tuple[int] -> cpp vector<int> sw=2
 * if python DataArrayInt -> cpp DataArrayInt sw=3
 * if python DataArrayIntTuple -> cpp DataArrayIntTuple sw=4
 *
 * switch between (int,vector<int>,DataArrayInt)
 */
static void convertObjToPossibleCpp1(PyObject *value, int& sw, int& iTyypp, std::vector<int>& stdvecTyypp, ParaMEDMEM::DataArrayInt *& daIntTyypp, ParaMEDMEM::DataArrayIntTuple *&daIntTuple) throw(INTERP_KERNEL::Exception)
{
  sw=-1;
  if(PyInt_Check(value))
    {
      iTyypp=(int)PyInt_AS_LONG(value);
      sw=1;
      return;
    }
  if(PyTuple_Check(value))
    {
      int size=PyTuple_Size(value);
      stdvecTyypp.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(value,i);
          if(PyInt_Check(o))
            stdvecTyypp[i]=(int)PyInt_AS_LONG(o);
          else
            {
              std::ostringstream oss; oss << "Tuple as been detected but element #" << i << " is not integer ! only tuples of integers accepted !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      sw=2;
      return;
    }
  if(PyList_Check(value))
    {
      int size=PyList_Size(value);
      stdvecTyypp.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(value,i);
          if(PyInt_Check(o))
            stdvecTyypp[i]=(int)PyInt_AS_LONG(o);
          else
            {
              std::ostringstream oss; oss << "List as been detected but element #" << i << " is not integer ! only lists of integers accepted !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      sw=2;
      return;
    }
  void *argp;
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,0|0);
  if(SWIG_IsOK(status))
    {
      daIntTyypp=reinterpret_cast< ParaMEDMEM::DataArrayInt * >(argp);
      sw=3;
      return;
    }
  status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayIntTuple,0|0);
  if(SWIG_IsOK(status))
    {  
      daIntTuple=reinterpret_cast< ParaMEDMEM::DataArrayIntTuple * >(argp);
      sw=4;
      return ;
    }
  throw INTERP_KERNEL::Exception("5 types accepted : integer, tuple of integer, list of integer, DataArrayInt, DataArrayIntTuple");
}

/*!
 * if python double -> cpp double sw=1
 * if python int -> cpp double sw=1
 * if python list[double] -> cpp vector<double> sw=2
 * if python list[int] -> cpp vector<double> sw=2
 * if python tuple[double] -> cpp vector<double> sw=2
 * if python tuple[int] -> cpp vector<double> sw=2
 * if python DataArrayDouble -> cpp DataArrayDouble sw=3
 *
 * switch between (int,vector<int>,DataArrayInt)
 */
static void convertObjToPossibleCpp4(PyObject *value, int& sw, double& iTyypp, std::vector<double>& stdvecTyypp, ParaMEDMEM::DataArrayDouble *& daIntTyypp) throw(INTERP_KERNEL::Exception)
{
  sw=-1;
  if(PyFloat_Check(value))
    {
      iTyypp=PyFloat_AS_DOUBLE(value);
      sw=1;
      return;
    }
  if(PyInt_Check(value))
    {
      iTyypp=(double)PyInt_AS_LONG(value);
      sw=1;
      return;
    }
  if(PyTuple_Check(value))
    {
      int size=PyTuple_Size(value);
      stdvecTyypp.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(value,i);
          if(PyFloat_Check(o))
            stdvecTyypp[i]=PyFloat_AS_DOUBLE(o);
          else if(PyInt_Check(o))
            stdvecTyypp[i]=(double)PyInt_AS_LONG(o);
          else
            {
              std::ostringstream oss; oss << "Tuple as been detected but element #" << i << " is not double ! only tuples of doubles accepted or integer !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      sw=2;
      return;
    }
  if(PyList_Check(value))
    {
      int size=PyList_Size(value);
      stdvecTyypp.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(value,i);
          if(PyFloat_Check(o))
            stdvecTyypp[i]=PyFloat_AS_DOUBLE(o);
          else if(PyInt_Check(o))
            stdvecTyypp[i]=(double)PyInt_AS_LONG(o);
          else
            {
              std::ostringstream oss; oss << "List as been detected but element #" << i << " is not double ! only lists of doubles accepted or integer !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      sw=2;
      return;
    }
  void *argp;
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble,0|0);
  if(!SWIG_IsOK(status))
    throw INTERP_KERNEL::Exception("5 types accepted : double float, integer, tuple of double float or int, list of double float or int, DataArrayDouble");
  daIntTyypp=reinterpret_cast< ParaMEDMEM::DataArrayDouble * >(argp);
  sw=3;
}

/*!
 * if python double -> cpp double sw=1
 * if python int -> cpp double sw=1
 * if python list[double] -> cpp vector<double> sw=2
 * if python list[int] -> cpp vector<double> sw=2
 * if python tuple[double] -> cpp vector<double> sw=2
 * if python tuple[int] -> cpp vector<double> sw=2
 * if python DataArrayDoubleTuple -> cpp DataArrayDoubleTuple sw=3
 *
 * switch between (int,vector<int>,DataArrayInt)
 */
static void convertObjToPossibleCpp44(PyObject *value, int& sw, double& iTyypp, std::vector<double>& stdvecTyypp, ParaMEDMEM::DataArrayDoubleTuple *& daIntTyypp) throw(INTERP_KERNEL::Exception)
{
  sw=-1;
  if(PyFloat_Check(value))
    {
      iTyypp=PyFloat_AS_DOUBLE(value);
      sw=1;
      return;
    }
  if(PyInt_Check(value))
    {
      iTyypp=(double)PyInt_AS_LONG(value);
      sw=1;
      return;
    }
  if(PyTuple_Check(value))
    {
      int size=PyTuple_Size(value);
      stdvecTyypp.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(value,i);
          if(PyFloat_Check(o))
            stdvecTyypp[i]=PyFloat_AS_DOUBLE(o);
          else if(PyInt_Check(o))
            stdvecTyypp[i]=(double)PyInt_AS_LONG(o);
          else
            {
              std::ostringstream oss; oss << "Tuple as been detected but element #" << i << " is not double ! only tuples of doubles accepted or integer !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      sw=2;
      return;
    }
  if(PyList_Check(value))
    {
      int size=PyList_Size(value);
      stdvecTyypp.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(value,i);
          if(PyFloat_Check(o))
            stdvecTyypp[i]=PyFloat_AS_DOUBLE(o);
          else if(PyInt_Check(o))
            stdvecTyypp[i]=(double)PyInt_AS_LONG(o);
          else
            {
              std::ostringstream oss; oss << "List as been detected but element #" << i << " is not double ! only lists of doubles accepted or integer !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      sw=2;
      return;
    }
  void *argp;
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayDoubleTuple,0|0);
  if(!SWIG_IsOK(status))
    throw INTERP_KERNEL::Exception("5 types accepted : double float, integer, tuple of double float or int, list of double float or int, DataArrayDoubleTuple");
  daIntTyypp=reinterpret_cast< ParaMEDMEM::DataArrayDoubleTuple * >(argp);
  sw=3;
}

/*!
 * if python int -> cpp int sw=1
 * if python list[int] -> cpp vector<int> sw=2
 * if python tuple[int] -> cpp vector<int> sw=2
 * if python slicp -> cpp pair sw=3
 * if python DataArrayInt -> cpp DataArrayInt sw=4
 *
 * switch between (int,vector<int>,DataArrayInt)
 */
static void convertObjToPossibleCpp2(PyObject *value, int nbelem, int& sw, int& iTyypp, std::vector<int>& stdvecTyypp, std::pair<int, std::pair<int,int> >& p, ParaMEDMEM::DataArrayInt *& daIntTyypp) throw(INTERP_KERNEL::Exception)
{
  const char *msg="5 types accepted : integer, tuple of integer, list of integer, slice, DataArrayInt, DataArrayIntTuple";
  sw=-1;
  if(PyInt_Check(value))
    {
      iTyypp=(int)PyInt_AS_LONG(value);
      sw=1;
      return;
    }
  if(PyTuple_Check(value))
    {
      int size=PyTuple_Size(value);
      stdvecTyypp.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(value,i);
          if(PyInt_Check(o))
            stdvecTyypp[i]=(int)PyInt_AS_LONG(o);
          else
            {
              std::ostringstream oss; oss << "Tuple as been detected but element #" << i << " is not integer ! only tuples of integers accepted !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      sw=2;
      return;
    }
  if(PyList_Check(value))
    {
      int size=PyList_Size(value);
      stdvecTyypp.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(value,i);
          if(PyInt_Check(o))
            stdvecTyypp[i]=(int)PyInt_AS_LONG(o);
          else
            {
              std::ostringstream oss; oss << "List as been detected but element #" << i << " is not integer ! only lists of integers accepted !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      sw=2;
      return;
    }
  if(PySlice_Check(value))
    {
      Py_ssize_t strt,stp,step;
      PySliceObject *oC=reinterpret_cast<PySliceObject *>(value);
      if(PySlice_GetIndices(oC,nbelem,&strt,&stp,&step)!=0)
        {
          std::ostringstream oss; oss << "Slice in subscriptable object DataArray invalid : number of elemnts is : " << nbelem;
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      p.first=strt;
      p.second.first=stp;
      p.second.second=step;
      sw=3;
      return ;
    }
  void *argp;
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,0|0);
  if(SWIG_IsOK(status))
    {
      daIntTyypp=reinterpret_cast< ParaMEDMEM::DataArrayInt * >(argp);
      if(!daIntTyypp)
        {
          std::ostringstream oss; oss << msg << " Instance in null !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      sw=4;
      return ;
    }
  status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayIntTuple,0|0);;
  if(SWIG_IsOK(status))
    {
      ParaMEDMEM::DataArrayIntTuple *tmp=reinterpret_cast< ParaMEDMEM::DataArrayIntTuple * >(argp);
      if(!tmp)
        {
          std::ostringstream oss; oss << msg << " Instance in null !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      stdvecTyypp.resize(tmp->getNumberOfCompo());
      std::copy(tmp->getConstPointer(),tmp->getConstPointer()+tmp->getNumberOfCompo(),stdvecTyypp.begin());
      sw=2;
      return ;
    }
  throw INTERP_KERNEL::Exception(msg);
}

static void convertObjToPossibleCpp22(PyObject *value, int nbelem, int& sw, int& iTyypp, std::vector<int>& stdvecTyypp, std::pair<int, std::pair<int,int> >& p, ParaMEDMEM::DataArrayIntTuple *& daIntTyypp) throw(INTERP_KERNEL::Exception)
{
  sw=-1;
  if(PyInt_Check(value))
    {
      iTyypp=(int)PyInt_AS_LONG(value);
      sw=1;
      return;
    }
  if(PyTuple_Check(value))
    {
      int size=PyTuple_Size(value);
      stdvecTyypp.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(value,i);
          if(PyInt_Check(o))
            stdvecTyypp[i]=(int)PyInt_AS_LONG(o);
          else
            {
              std::ostringstream oss; oss << "Tuple as been detected but element #" << i << " is not integer ! only tuples of integers accepted !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      sw=2;
      return;
    }
  if(PyList_Check(value))
    {
      int size=PyList_Size(value);
      stdvecTyypp.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(value,i);
          if(PyInt_Check(o))
            stdvecTyypp[i]=(int)PyInt_AS_LONG(o);
          else
            {
              std::ostringstream oss; oss << "List as been detected but element #" << i << " is not integer ! only lists of integers accepted !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      sw=2;
      return;
    }
  if(PySlice_Check(value))
    {
      Py_ssize_t strt,stp,step;
      PySliceObject *oC=reinterpret_cast<PySliceObject *>(value);
      if(PySlice_GetIndices(oC,nbelem,&strt,&stp,&step)!=0)
        {
          std::ostringstream oss; oss << "Slice in subscriptable object DataArray invalid : number of elemnts is : " << nbelem;
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      p.first=strt;
      p.second.first=stp;
      p.second.second=step;
      sw=3;
      return ;
    }
  void *argp;
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayIntTuple,0|0);
  if(!SWIG_IsOK(status))
    throw INTERP_KERNEL::Exception("4 types accepted : integer, tuple of integer, list of integer, slice, DataArrayIntTuple");
  daIntTyypp=reinterpret_cast< ParaMEDMEM::DataArrayIntTuple * >(argp);
  sw=4;
}

/*!
 * if value int -> cpp it sw=1
 * if value list[int] -> vt sw=2
 * if value tuple[int] -> vt sw=2
 * if value slice -> pt sw=3
 * if value DataArrayInt -> dt sw=4
 * if value tuple [int,int] -> cpp it,ip sw=5
 * if value tuple [list[int],int] -> cpp vt,ip sw=6
 * if value tuple [tuple[int],int] -> cpp vt,ip sw=6
 * if value tuple [slice,int] -> cpp pt,ip sw=7
 * if value tuple [DaI,int] -> cpp dt,ip sw=8
 * if value tuple [int,list[int]] -> cpp it,vc sw=9
 * if value tuple [list[int],list[int]] -> cpp vt,vc sw=10
 * if value tuple [tuple[int],list[int]] -> cpp vt,vc sw=10
 * if value tuple [slice,list[int]] -> cpp pt,vc sw=11
 * if value tuple [DaI,list[int]] -> cpp dt,vc sw=12
 * if value tuple [int,tuple[int]] -> cpp it,vc sw=9
 * if value tuple [list[int],tuple[int]] -> cpp vt,vc sw=10
 * if value tuple [tuple[int],tuple[int]] -> cpp vt,vc sw=10
 * if value tuple [slice,tuple[int]] -> cpp pt,vc sw=11
 * if value tuple [DaI,tuple[int]] -> cpp dt,vc sw=12
 * if value tuple [int,slice] -> cpp it,pc sw=13
 * if value tuple [list[int],slice] -> cpp vt,pc sw=14
 * if value tuple [tuple[int],slice] -> cpp vt,pc sw=14
 * if value tuple [slice,slice] -> cpp pt,pc sw=15
 * if value tuple [DaI,slice] -> cpp dt,pc sw=16
 *
 * switch between (int,vector<int>,DataArrayInt)
 */
static void convertObjToPossibleCpp3(PyObject *value, int nbTuple, int nbCompo, int& sw, int& it, int& ic, std::vector<int>& vt, std::vector<int>& vc,
                                     std::pair<int, std::pair<int,int> >& pt, std::pair<int, std::pair<int,int> >& pc,
                                     ParaMEDMEM::DataArrayInt *&dt, ParaMEDMEM::DataArrayInt *&dc) throw(INTERP_KERNEL::Exception)
{
  if(!PyTuple_Check(value))
    {
      convertObjToPossibleCpp2(value,nbTuple,sw,it,vt,pt,dt);
      return ;
    }
  else
    {
      int sz=PyTuple_Size(value);
      if(sz!=2)
        throw INTERP_KERNEL::Exception("Unexpected nb of slice element : 1 or 2 expected !\n1st is for tuple selection, 2nd for component selection !");
      PyObject *ob0=PyTuple_GetItem(value,0);
      int sw1,sw2;
      convertObjToPossibleCpp2(ob0,nbTuple,sw1,it,vt,pt,dt);
      PyObject *ob1=PyTuple_GetItem(value,1);
      convertObjToPossibleCpp2(ob1,nbCompo,sw2,ic,vc,pc,dc);
      sw=4*sw2+sw1;
    }
}

/*!
 * if value int -> cpp val sw=1
 * if value double -> cpp val sw=1
 * if value DataArrayDouble -> cpp DataArrayDouble sw=2
 * if value DataArrayDoubleTuple -> cpp DataArrayDoubleTuple sw=3
 * if value list[int,double] -> cpp std::vector<double> sw=4
 * if value tuple[int,double] -> cpp std::vector<double> sw=4
 */
static void convertObjToPossibleCpp5(PyObject *value, int& sw, double& val, ParaMEDMEM::DataArrayDouble *&d, ParaMEDMEM::DataArrayDoubleTuple *&e, std::vector<double>& f)
{
  sw=-1;
  if(PyFloat_Check(value))
    {
      val=PyFloat_AS_DOUBLE(value);
      sw=1;
      return;
    }
  if(PyInt_Check(value))
    {
      val=(double)PyInt_AS_LONG(value);
      sw=1;
      return;
    }
  if(PyTuple_Check(value))
    {
      int size=PyTuple_Size(value);
      f.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(value,i);
          if(PyFloat_Check(o))
            f[i]=PyFloat_AS_DOUBLE(o);
          else if(PyInt_Check(o))
            f[i]=(double)PyInt_AS_LONG(o);
          else
            {
              std::ostringstream oss; oss << "Tuple as been detected but element #" << i << " is not double ! only tuples of doubles accepted or integer !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      sw=4;
      return;
    }
  if(PyList_Check(value))
    {
      int size=PyList_Size(value);
      f.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(value,i);
          if(PyFloat_Check(o))
            f[i]=PyFloat_AS_DOUBLE(o);
          else if(PyInt_Check(o))
            f[i]=(double)PyInt_AS_LONG(o);
          else
            {
              std::ostringstream oss; oss << "List as been detected but element #" << i << " is not double ! only lists of doubles accepted or integer !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      sw=4;
      return;
    }
  void *argp;
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble,0|0);
  if(SWIG_IsOK(status))
    {  
      d=reinterpret_cast< ParaMEDMEM::DataArrayDouble * >(argp);
      sw=2;
      return ;
    }
  status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayDoubleTuple,0|0);
  if(SWIG_IsOK(status))
    {  
      e=reinterpret_cast< ParaMEDMEM::DataArrayDoubleTuple * >(argp);
      sw=3;
      return ;
    }
  throw INTERP_KERNEL::Exception("4 types accepted : integer, double, DataArrayDouble, DataArrayDoubleTuple");
}
