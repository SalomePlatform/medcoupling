// Copyright (C) 2007-2017  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (EDF R&D)

#ifndef __MEDCOUPLINGDATAARRAYTRAITS_HXX__
#define __MEDCOUPLINGDATAARRAYTRAITS_HXX__

#include "MEDCouplingMemArray.hxx"

#include <Python.h>

#ifdef WITH_NUMPY
#include <numpy/arrayobject.h>
#if NPY_API_VERSION <= 0x00000006
#  define MED_NUMPY_OWNDATA NPY_OWNDATA
#else
#  define MED_NUMPY_OWNDATA NPY_ARRAY_OWNDATA
#endif
#endif

#ifdef WITH_NUMPY
// specific DataArray deallocator callback. This deallocator is used both in the constructor of DataArray and in the toNumPyArr
// method. This dellocator uses weakref to determine if the linked numArr is still alive or not. If alive the ownership is given to it.
// if no more alive the "standart" DataArray deallocator is called.
void numarrdeal(void *pt, void *wron)
{
  void **wronc=(void **)wron;
  PyObject *weakRefOnOwner=reinterpret_cast<PyObject *>(wronc[0]);
  PyObject *obj=PyWeakref_GetObject(weakRefOnOwner);
  if(obj!=Py_None)
    {
      Py_XINCREF(obj);
      PyArrayObject *objC=reinterpret_cast<PyArrayObject *>(obj);
      objC->flags|=MED_NUMPY_OWNDATA;
      Py_XDECREF(weakRefOnOwner);
      Py_XDECREF(obj);
    }
  else
    {
      typedef void (*MyDeallocator)(void *,void *);
      MyDeallocator deall=(MyDeallocator)wronc[1];
      deall(pt,NULL);
      Py_XDECREF(weakRefOnOwner);
    }
  delete [] wronc;
}
#endif

template<class MCData>
struct PyCallBackDataArraySt {
    PyObject_HEAD
    MCData *_pt_mc;
};

typedef struct PyCallBackDataArraySt<MEDCoupling::DataArrayByte> PyCallBackDataArrayChar;
typedef struct PyCallBackDataArraySt<MEDCoupling::DataArrayInt> PyCallBackDataArrayInt;
typedef struct PyCallBackDataArraySt<MEDCoupling::DataArrayFloat> PyCallBackDataArrayFloat;
typedef struct PyCallBackDataArraySt<MEDCoupling::DataArrayDouble> PyCallBackDataArrayDouble;

extern "C"
{
  static int callbackmcdataarray___init__(PyObject *self, PyObject *args, PyObject *kwargs) { return 0; }
  
  static PyObject *callbackmcdataarraychar___new__(PyTypeObject *type, PyObject *args, PyObject *kwargs)
  {
    PyCallBackDataArrayChar *self = (PyCallBackDataArrayChar *) ( type->tp_alloc(type, 0) );
    return (PyObject *)self;
  }

  static PyObject *callbackmcdataarrayint___new__(PyTypeObject *type, PyObject *args, PyObject *kwargs)
  {
    PyCallBackDataArrayInt *self = (PyCallBackDataArrayInt *) ( type->tp_alloc(type, 0) );
    return (PyObject *)self;
  }
  
  static PyObject *callbackmcdataarrayfloat___new__(PyTypeObject *type, PyObject *args, PyObject *kwargs)
  {
    PyCallBackDataArrayFloat *self = (PyCallBackDataArrayFloat *) ( type->tp_alloc(type, 0) );
    return (PyObject *)self;
  }
  
  static PyObject *callbackmcdataarraydouble___new__(PyTypeObject *type, PyObject *args, PyObject *kwargs)
  {
    PyCallBackDataArrayDouble *self = (PyCallBackDataArrayDouble *) ( type->tp_alloc(type, 0) );
    return (PyObject *)self;
  }
  
  static void callbackmcdataarray_dealloc(PyObject *self)
  {
    Py_TYPE(self)->tp_free(self);
  }

  
  // real callback called when a numpy arr having more than one DataArray instance client on it is destroyed.
  // In this case, all the "weak" clients, except the first one, invoke this call back that desable the content of these "weak" clients.
  static PyObject *callbackmcdataarraychar_call(PyCallBackDataArrayChar *self, PyObject *args, PyObject *kw)
  {
    if(self->_pt_mc)
      {
        MEDCoupling::MemArray<char>& mma=self->_pt_mc->accessToMemArray();
        mma.destroy();
      }
    Py_XINCREF(Py_None);
    return Py_None;
  }

  // real callback called when a numpy arr having more than one DataArray instance client on it is destroyed.
  // In this case, all the "weak" clients, except the first one, invoke this call back that desable the content of these "weak" clients.
  static PyObject *callbackmcdataarrayint_call(PyCallBackDataArrayInt *self, PyObject *args, PyObject *kw)
  {
    if(self->_pt_mc)
      {
        MEDCoupling::MemArray<int>& mma=self->_pt_mc->accessToMemArray();
        mma.destroy();
      }
    Py_XINCREF(Py_None);
    return Py_None;
  }

  // real callback called when a numpy arr having more than one DataArray instance client on it is destroyed.
  // In this case, all the "weak" clients, except the first one, invoke this call back that desable the content of these "weak" clients.
  static PyObject *callbackmcdataarrayfloat_call(PyCallBackDataArrayFloat *self, PyObject *args, PyObject *kw)
  {
    if(self->_pt_mc)
      {
        MEDCoupling::MemArray<float>& mma=self->_pt_mc->accessToMemArray();
        mma.destroy();
      }
    Py_XINCREF(Py_None);
    return Py_None;
  }
  
  // real callback called when a numpy arr having more than one DataArray instance client on it is destroyed.
  // In this case, all the "weak" clients, except the first one, invoke this call back that desable the content of these "weak" clients.
  static PyObject *callbackmcdataarraydouble_call(PyCallBackDataArrayDouble *self, PyObject *args, PyObject *kw)
  {
    if(self->_pt_mc)
      {
        MEDCoupling::MemArray<double>& mma=self->_pt_mc->accessToMemArray();
        mma.destroy();
      }
    Py_XINCREF(Py_None);
    return Py_None;
  }
}

PyTypeObject PyCallBackDataArrayChar_RefType = {
  PyVarObject_HEAD_INIT(&PyType_Type, 0)
  "callbackmcdataarraychar",
  sizeof(PyCallBackDataArrayChar),
  0,
  callbackmcdataarray_dealloc,            /*tp_dealloc*/
  0,                          /*tp_print*/
  0,                          /*tp_getattr*/
  0,                          /*tp_setattr*/
  0,                          /*tp_compare*/
  0,                          /*tp_repr*/
  0,                          /*tp_as_number*/
  0,                          /*tp_as_sequence*/
  0,                          /*tp_as_mapping*/
  0,                          /*tp_hash*/
  (ternaryfunc)callbackmcdataarraychar_call,  /*tp_call*/
  0,                          /*tp_str*/
  0,                          /*tp_getattro*/
  0,                          /*tp_setattro*/
  0,                          /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_BASETYPE,  /*tp_flags*/
  0,                          /*tp_doc*/
  0,                          /*tp_traverse*/
  0,                          /*tp_clear*/
  0,                          /*tp_richcompare*/
  0,                          /*tp_weaklistoffset*/
  0,                          /*tp_iter*/
  0,                          /*tp_iternext*/
  0,                          /*tp_methods*/
  0,                          /*tp_members*/
  0,                          /*tp_getset*/
  0,                          /*tp_base*/
  0,                          /*tp_dict*/
  0,                          /*tp_descr_get*/
  0,                          /*tp_descr_set*/
  0,                          /*tp_dictoffset*/
  callbackmcdataarray___init__,           /*tp_init*/
  PyType_GenericAlloc,        /*tp_alloc*/
  callbackmcdataarraychar___new__,            /*tp_new*/
  PyObject_GC_Del,            /*tp_free*/
};


PyTypeObject PyCallBackDataArrayInt_RefType = {
  PyVarObject_HEAD_INIT(&PyType_Type, 0)
  "callbackmcdataarrayint",
  sizeof(PyCallBackDataArrayInt),
  0,
  callbackmcdataarray_dealloc,            /*tp_dealloc*/
  0,                          /*tp_print*/
  0,                          /*tp_getattr*/
  0,                          /*tp_setattr*/
  0,                          /*tp_compare*/
  0,                          /*tp_repr*/
  0,                          /*tp_as_number*/
  0,                          /*tp_as_sequence*/
  0,                          /*tp_as_mapping*/
  0,                          /*tp_hash*/
  (ternaryfunc)callbackmcdataarrayint_call,  /*tp_call*/
  0,                          /*tp_str*/
  0,                          /*tp_getattro*/
  0,                          /*tp_setattro*/
  0,                          /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_BASETYPE,  /*tp_flags*/
  0,                          /*tp_doc*/
  0,                          /*tp_traverse*/
  0,                          /*tp_clear*/
  0,                          /*tp_richcompare*/
  0,                          /*tp_weaklistoffset*/
  0,                          /*tp_iter*/
  0,                          /*tp_iternext*/
  0,                          /*tp_methods*/
  0,                          /*tp_members*/
  0,                          /*tp_getset*/
  0,                          /*tp_base*/
  0,                          /*tp_dict*/
  0,                          /*tp_descr_get*/
  0,                          /*tp_descr_set*/
  0,                          /*tp_dictoffset*/
  callbackmcdataarray___init__,           /*tp_init*/
  PyType_GenericAlloc,        /*tp_alloc*/
  callbackmcdataarrayint___new__,            /*tp_new*/
  PyObject_GC_Del,            /*tp_free*/
};

PyTypeObject PyCallBackDataArrayFloat_RefType = {
  PyVarObject_HEAD_INIT(&PyType_Type, 0)
  "callbackmcdataarraydouble",
  sizeof(PyCallBackDataArrayFloat),
  0,
  callbackmcdataarray_dealloc,            /*tp_dealloc*/
  0,                          /*tp_print*/
  0,                          /*tp_getattr*/
  0,                          /*tp_setattr*/
  0,                          /*tp_compare*/
  0,                          /*tp_repr*/
  0,                          /*tp_as_number*/
  0,                          /*tp_as_sequence*/
  0,                          /*tp_as_mapping*/
  0,                          /*tp_hash*/
  (ternaryfunc)callbackmcdataarraydouble_call,  /*tp_call*/
  0,                          /*tp_str*/
  0,                          /*tp_getattro*/
  0,                          /*tp_setattro*/
  0,                          /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_BASETYPE,  /*tp_flags*/
  0,                          /*tp_doc*/
  0,                          /*tp_traverse*/
  0,                          /*tp_clear*/
  0,                          /*tp_richcompare*/
  0,                          /*tp_weaklistoffset*/
  0,                          /*tp_iter*/
  0,                          /*tp_iternext*/
  0,                          /*tp_methods*/
  0,                          /*tp_members*/
  0,                          /*tp_getset*/
  0,                          /*tp_base*/
  0,                          /*tp_dict*/
  0,                          /*tp_descr_get*/
  0,                          /*tp_descr_set*/
  0,                          /*tp_dictoffset*/
  callbackmcdataarray___init__,           /*tp_init*/
  PyType_GenericAlloc,        /*tp_alloc*/
  callbackmcdataarrayfloat___new__,            /*tp_new*/
  PyObject_GC_Del,            /*tp_free*/
};

PyTypeObject PyCallBackDataArrayDouble_RefType = {
  PyVarObject_HEAD_INIT(&PyType_Type, 0)
  "callbackmcdataarraydouble",
  sizeof(PyCallBackDataArrayDouble),
  0,
  callbackmcdataarray_dealloc,            /*tp_dealloc*/
  0,                          /*tp_print*/
  0,                          /*tp_getattr*/
  0,                          /*tp_setattr*/
  0,                          /*tp_compare*/
  0,                          /*tp_repr*/
  0,                          /*tp_as_number*/
  0,                          /*tp_as_sequence*/
  0,                          /*tp_as_mapping*/
  0,                          /*tp_hash*/
  (ternaryfunc)callbackmcdataarraydouble_call,  /*tp_call*/
  0,                          /*tp_str*/
  0,                          /*tp_getattro*/
  0,                          /*tp_setattro*/
  0,                          /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_BASETYPE,  /*tp_flags*/
  0,                          /*tp_doc*/
  0,                          /*tp_traverse*/
  0,                          /*tp_clear*/
  0,                          /*tp_richcompare*/
  0,                          /*tp_weaklistoffset*/
  0,                          /*tp_iter*/
  0,                          /*tp_iternext*/
  0,                          /*tp_methods*/
  0,                          /*tp_members*/
  0,                          /*tp_getset*/
  0,                          /*tp_base*/
  0,                          /*tp_dict*/
  0,                          /*tp_descr_get*/
  0,                          /*tp_descr_set*/
  0,                          /*tp_dictoffset*/
  callbackmcdataarray___init__,           /*tp_init*/
  PyType_GenericAlloc,        /*tp_alloc*/
  callbackmcdataarraydouble___new__,            /*tp_new*/
  PyObject_GC_Del,            /*tp_free*/
};

#ifdef WITH_NUMPY
template<class T>
struct NPYTraits
{
};

template<>
struct NPYTraits<double>
{
  static const int NPYObjectType=NPY_DOUBLE;
  static PyTypeObject *NPYFunc;
  static PyObject *Array_SWIGTYPE;
};

template<>
struct NPYTraits<float>
{
  static const int NPYObjectType=NPY_FLOAT;
  static PyTypeObject *NPYFunc;
  static PyObject *Array_SWIGTYPE;
};
#endif

#endif
