// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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

#include "InterpKernelAutoPtr.hxx"

/*!
 * This method is an extention of PySlice_GetIndices but less
 * open than PySlice_GetIndicesEx that accepts too many situations.
 */
void GetIndicesOfSlice(PySliceObject *slice, Py_ssize_t length, Py_ssize_t *start, Py_ssize_t *stop, Py_ssize_t *step, const char *msgInCaseOfFailure)
{
  int ret(PySlice_GetIndices(slice,length,start,stop,step));
  if(ret==0)
    return ;
  if(*step>0 && *start==*stop && length==*start)
    return ;
  throw INTERP_KERNEL::Exception(msgInCaseOfFailure);
}

/*!
 * This method allows to retrieve slice info from \a slice.
 */
void GetIndicesOfSliceExplicitely(PySliceObject *slice, Py_ssize_t *start, Py_ssize_t *stop, Py_ssize_t *step, const char *msgInCaseOfFailure)
{
  int ret(PySlice_GetIndices(slice,std::numeric_limits<int>::max(),start,stop,step));
  if(ret==0)
    {
      if(*start!=std::numeric_limits<int>::max() && *stop!=std::numeric_limits<int>::max())
        return ;
      std::ostringstream oss;
      oss << msgInCaseOfFailure << " The input slice contains some unknowns that can't be determined in static method ! The input slice must be explicit here !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  throw INTERP_KERNEL::Exception(msgInCaseOfFailure);
}

int InterpreteNegativeInt(int val, int nbelem)
{
  if(val<0)
    {
      int newVal(nbelem+val);
      if(newVal<0)
        {
          std::ostringstream oss; oss << "interpreteNegativeInt : request for negative int=" << val << " but number of elems is equal to " << nbelem << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      return newVal;
    }
  else
    return val;
}

#ifdef WITH_NUMPY
#include <numpy/arrayobject.h>

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
      objC->flags|=NPY_OWNDATA;
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

template<class MCData>
struct PyCallBackDataArraySt {
    PyObject_HEAD
    MCData *_pt_mc;
};

typedef struct PyCallBackDataArraySt<ParaMEDMEM::DataArrayInt> PyCallBackDataArrayInt;
typedef struct PyCallBackDataArraySt<ParaMEDMEM::DataArrayDouble> PyCallBackDataArrayDouble;

extern "C"
{
  static int callbackmcdataarray___init__(PyObject *self, PyObject *args, PyObject *kwargs) { return 0; }
  
  static PyObject *callbackmcdataarrayint___new__(PyTypeObject *type, PyObject *args, PyObject *kwargs)
  {
    PyCallBackDataArrayInt *self = (PyCallBackDataArrayInt *) ( type->tp_alloc(type, 0) );
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
  static PyObject *callbackmcdataarrayint_call(PyCallBackDataArrayInt *self, PyObject *args, PyObject *kw)
  {
    if(self->_pt_mc)
      {
        ParaMEDMEM::MemArray<int>& mma=self->_pt_mc->accessToMemArray();
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
        ParaMEDMEM::MemArray<double>& mma=self->_pt_mc->accessToMemArray();
        mma.destroy();
      }
    Py_XINCREF(Py_None);
    return Py_None;
  }
}

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

// this is the second type of specific deallocator, only valid for the constructor of DataArrays taking numpy array
// in input when an another DataArray is already client of this.
template<class MCData>
void numarrdeal2(void *pt, void *obj)
{
  typedef struct PyCallBackDataArraySt<MCData> PyCallBackDataArray;
  void **obj1=(void **)obj;
  PyCallBackDataArray *cbdaic=reinterpret_cast<PyCallBackDataArray *>(obj1[0]);
  PyObject *weakRefOnOwner=reinterpret_cast<PyObject *>(obj1[1]);
  cbdaic->_pt_mc=0;
  Py_XDECREF(weakRefOnOwner);
  Py_XDECREF(cbdaic);
  delete [] obj1;
}

template<class MCData, class T>
MCData *BuildNewInstance(PyObject *elt0, int npyObjectType, PyTypeObject *pytype, const char *msg)
{
  int ndim=PyArray_NDIM(elt0);
  if(ndim!=1 && ndim!=2)
    throw INTERP_KERNEL::Exception("Input numpy array should have dimension equal to 1 or 2 !");
  if(PyArray_DESCR(elt0)->type_num != npyObjectType)
    {
      std::ostringstream oss; oss << "Input numpy array has not the type " << msg << "!";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  npy_intp sz0=PyArray_DIM(elt0,0);
  npy_intp sz1=ndim==2?PyArray_DIM(elt0,1):1;
  //
  int itemSize=PyArray_ITEMSIZE(elt0);
  if(itemSize!=sizeof(T))
    {
      std::ostringstream oss; oss << "Input numpy array has not itemSize set to " << sizeof(T) << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(itemSize*sz1!=PyArray_STRIDE(elt0,0))
    throw INTERP_KERNEL::Exception("Input numpy array has stride that mismatches the item size ! Data are not packed in the right way for DataArrays !");
  if(ndim==2)
    if(itemSize!=PyArray_STRIDE(elt0,1))
      throw INTERP_KERNEL::Exception("Input numpy array has stride that mismatches the item size ! Data are not packed in the right way for DataArrays for component #1 !");
  const char *data=PyArray_BYTES(elt0);
  typename ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<MCData> ret=MCData::New();
  if(PyArray_ISBEHAVED(elt0))//aligned and writeable and in machine byte-order
    {
      PyArrayObject *elt0C=reinterpret_cast<PyArrayObject *>(elt0);
      PyArrayObject *eltOwning=(PyArray_FLAGS(elt0C) & NPY_OWNDATA)?elt0C:NULL;
      int mask=NPY_OWNDATA; mask=~mask;
      elt0C->flags&=mask;
      PyObject *deepestObj=elt0;
      PyObject *base=elt0C->base;
      if(base) deepestObj=base;
      while(base)
        {
          if(PyArray_Check(base))
            {
              PyArrayObject *baseC=reinterpret_cast<PyArrayObject *>(base);
              eltOwning=(PyArray_FLAGS(baseC) & NPY_OWNDATA)?baseC:eltOwning;
              baseC->flags&=mask;
              base=baseC->base;
              if(base) deepestObj=base;
            }
          else
            break;
        }
      typename ParaMEDMEM::MemArray<T>& mma=ret->accessToMemArray();
      if(eltOwning==NULL)
        {
          PyCallBackDataArraySt<MCData> *cb=PyObject_GC_New(PyCallBackDataArraySt<MCData>,pytype);
          cb->_pt_mc=ret;
          ret->useArray(reinterpret_cast<const T *>(data),true,ParaMEDMEM::C_DEALLOC,sz0,sz1);
          PyObject *ref=PyWeakref_NewRef(deepestObj,(PyObject *)cb);
          void **objs=new void *[2]; objs[0]=cb; objs[1]=ref;
          mma.setParameterForDeallocator(objs);
          mma.setSpecificDeallocator(numarrdeal2<MCData>);
          //"Impossible to share this numpy array chunk of data, because already shared by an another non numpy array object (maybe an another DataArrayInt instance) ! Release it, or perform a copy on the input array !");
        }
      else
        {
          ret->useArray(reinterpret_cast<const T *>(data),true,ParaMEDMEM::C_DEALLOC,sz0,sz1);
          PyObject *ref=PyWeakref_NewRef(reinterpret_cast<PyObject *>(eltOwning),NULL);
          typename ParaMEDMEM::MemArray<T>::Deallocator tmp(ParaMEDMEM::MemArray<T>::CDeallocator);
          void **tmp2 = reinterpret_cast<void**>(&tmp); // MSVC2010 does not support constructor()
          void **objs=new void *[2]; objs[0]=ref; objs[1]=*tmp2;
          mma.setParameterForDeallocator(objs);
          mma.setSpecificDeallocator(numarrdeal);
        }
    }
  else if(PyArray_ISBEHAVED_RO(elt0))
    ret->useArray(reinterpret_cast<const T *>(data),false,ParaMEDMEM::CPP_DEALLOC,sz0,sz1);
  return ret.retn();
}


int NumpyArrSetBaseObjectExt(PyArrayObject *arr, PyObject *obj)
{
    if (obj == NULL) {
        PyErr_SetString(PyExc_ValueError,
                "Cannot set the NumPy array 'base' "
                "dependency to NULL after initialization");
        return -1;
    }
    /*
     * Allow the base to be set only once. Once the object which
     * owns the data is set, it doesn't make sense to change it.
     */
    if (PyArray_BASE(arr) != NULL) {
        Py_DECREF(obj);
        PyErr_SetString(PyExc_ValueError,
                "Cannot set the NumPy array 'base' "
                "dependency more than once");
        return -1;
    }

    /*
     * Don't allow infinite chains of views, always set the base
     * to the first owner of the data.  
     * That is, either the first object which isn't an array, 
     * or the first object which owns its own data.
     */

    while (PyArray_Check(obj) && (PyObject *)arr != obj) {
        PyArrayObject *obj_arr = (PyArrayObject *)obj;
        PyObject *tmp;
 

        /* If this array owns its own data, stop collapsing */
        if (PyArray_CHKFLAGS(obj_arr, NPY_OWNDATA)) {
            break;
        }   

        tmp = PyArray_BASE(obj_arr);
        /* If there's no base, stop collapsing */
        if (tmp == NULL) {
            break;
        }
        /* Stop the collapse new base when the would not be of the same 
         * type (i.e. different subclass).
         */
        if (Py_TYPE(tmp) != Py_TYPE(arr)) {
            break;
        }


        Py_INCREF(tmp);
        Py_DECREF(obj);
        obj = tmp;
    }

    /* Disallow circular references */
    if ((PyObject *)arr == obj) {
        Py_DECREF(obj);
        PyErr_SetString(PyExc_ValueError,
                "Cannot create a circular NumPy array 'base' dependency");
        return -1;
    }

    arr->base = obj;

    return 0;
}

template<class MCData, class T>
PyObject *ToNumPyArray(MCData *self, int npyObjectType, const char *MCDataStr)
{
  if(!self->isAllocated())
    {
      std::ostringstream oss; oss << MCDataStr << "::toNumPyArray : this is not allocated !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  ParaMEDMEM::MemArray<T>& mem=self->accessToMemArray();
  int nbComp=self->getNumberOfComponents();
  if(nbComp==0)
    {
      std::ostringstream oss; oss << MCDataStr << "::toNumPyArray : number of components of this is 0 ! Should be > 0 !"; 
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  int nbDims=nbComp==1?1:2;
  npy_intp dim[2];
  dim[0]=(npy_intp)self->getNumberOfTuples(); dim[1]=nbComp;
  const T *bg=self->getConstPointer();
  PyObject *ret(PyArray_SimpleNewFromData(nbDims,dim,npyObjectType,const_cast<T *>(bg)));
  if(mem.isDeallocatorCalled())
    {
      if(mem.getDeallocator()!=numarrdeal)
        {// case for the first call of toNumPyArray
          PyObject *ref(PyWeakref_NewRef(ret,NULL));
          typename ParaMEDMEM::MemArray<T>::Deallocator tmp(mem.getDeallocator());
          void **tmp2 = reinterpret_cast<void**>(&tmp); // MSVC2010 does not support constructor()
          void **objs=new void *[2]; objs[0]=reinterpret_cast<void*>(ref); objs[1]=*tmp2;
          mem.setParameterForDeallocator(objs);
          mem.setSpecificDeallocator(numarrdeal);
          return ret;
        }
      else
        {// case for the second and other call of toNumPyArray
          void **objs=(void **)mem.getParameterForDeallocator();
          PyObject *weakRefOnOwner=(PyObject *)objs[0];
          PyObject *obj=PyWeakref_GetObject(weakRefOnOwner);
          if(obj!=Py_None)
            {//the previous numArray exists let numpy deals the numpy array each other by declaring the still alive instance as base
              Py_XINCREF(obj);
              NumpyArrSetBaseObjectExt((PyArrayObject*)ret,obj);
            }
          else
            {//the previous numArray no more exists -> declare the newly created numpy array as the first one.
              Py_XDECREF(weakRefOnOwner);
              PyObject *ref=PyWeakref_NewRef(ret,NULL);
              objs[0]=ref;
            }
        }
    }
  return ret;
}

SWIGINTERN PyObject *ParaMEDMEM_DataArrayInt_toNumPyArray(ParaMEDMEM::DataArrayInt *self);
SWIGINTERN PyObject *ParaMEDMEM_DataArrayDouble_toNumPyArray(ParaMEDMEM::DataArrayDouble *self);

PyObject *ToCSRMatrix(const std::vector<std::map<int,double> >& m, int nbCols) throw(INTERP_KERNEL::Exception)
{
  int nbRows((int)m.size());
  ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayInt> indPtr(ParaMEDMEM::DataArrayInt::New()),indices(ParaMEDMEM::DataArrayInt::New());
  ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> data(ParaMEDMEM::DataArrayDouble::New());
  indPtr->alloc(nbRows+1,1);
  int *intPtr_ptr(indPtr->getPointer()); intPtr_ptr[0]=0; intPtr_ptr++;
  int sz2(0);
  for(std::vector<std::map<int,double> >::const_iterator it0=m.begin();it0!=m.end();it0++,intPtr_ptr++)
    {
      sz2+=(int)(*it0).size();
      *intPtr_ptr=sz2;
    }
  indices->alloc(sz2,1); data->alloc(sz2,1);
  int *indices_ptr(indices->getPointer());
  double *data_ptr(data->getPointer());
  for(std::vector<std::map<int,double> >::const_iterator it0=m.begin();it0!=m.end();it0++)
    for(std::map<int,double>::const_iterator it1=(*it0).begin();it1!=(*it0).end();it1++,indices_ptr++,data_ptr++)
      {
        *indices_ptr=(*it1).first;
        *data_ptr=(*it1).second;
      }
  PyObject *a(ParaMEDMEM_DataArrayDouble_toNumPyArray(data)),*b(ParaMEDMEM_DataArrayInt_toNumPyArray(indices)),*c(ParaMEDMEM_DataArrayInt_toNumPyArray(indPtr));
  //
  PyObject *args(PyTuple_New(1)),*args0(PyTuple_New(3)),*kw(PyDict_New()),*kw1(PyTuple_New(2));
  PyTuple_SetItem(args0,0,a); PyTuple_SetItem(args0,1,b); PyTuple_SetItem(args0,2,c); PyTuple_SetItem(args,0,args0);
  PyTuple_SetItem(kw1,0,PyInt_FromLong(nbRows)); PyTuple_SetItem(kw1,1,PyInt_FromLong(nbCols));
  PyObject *tmp1(PyString_FromString("shape"));
  PyDict_SetItem(kw,tmp1,kw1); Py_DECREF(tmp1); Py_DECREF(kw1);
  PyObject* pdict=PyDict_New();
  PyDict_SetItemString(pdict, "__builtins__", PyEval_GetBuiltins());
  PyObject *tmp(PyRun_String("from scipy.sparse import csr_matrix", Py_single_input, pdict, pdict));
  if(!tmp)
    throw INTERP_KERNEL::Exception("Problem during loading csr_matrix in scipy.sparse ! Is Scipy module available in present ?");
  PyObject *csrMatrixCls=PyDict_GetItemString(pdict,"csr_matrix");
  if(!csrMatrixCls)
    throw INTERP_KERNEL::Exception("csr_matrix not found in scipy.sparse ! Is Scipy module available in present ?");
  PyObject *ret(PyObject_Call(csrMatrixCls,args,kw));
  Py_DECREF(pdict); Py_XDECREF(tmp); Py_DECREF(args); Py_DECREF(kw);
  return ret;
}

#endif

static PyObject *convertDataArrayChar(ParaMEDMEM::DataArrayChar *dac, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!dac)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<ParaMEDMEM::DataArrayByte *>(dac))
    ret=SWIG_NewPointerObj((void*)dac,SWIGTYPE_p_ParaMEDMEM__DataArrayByte,owner);
  if(dynamic_cast<ParaMEDMEM::DataArrayAsciiChar *>(dac))
    ret=SWIG_NewPointerObj((void*)dac,SWIGTYPE_p_ParaMEDMEM__DataArrayAsciiChar,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of DataArrayChar on downcast !");
  return ret;
}

static PyObject *convertDataArray(ParaMEDMEM::DataArray *dac, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!dac)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<ParaMEDMEM::DataArrayDouble *>(dac))
    ret=SWIG_NewPointerObj((void*)dac,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble,owner);
  if(dynamic_cast<ParaMEDMEM::DataArrayInt *>(dac))
    ret=SWIG_NewPointerObj((void*)dac,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,owner);
  if(dynamic_cast<ParaMEDMEM::DataArrayByte *>(dac))
    ret=SWIG_NewPointerObj((void*)dac,SWIGTYPE_p_ParaMEDMEM__DataArrayByte,owner);
  if(dynamic_cast<ParaMEDMEM::DataArrayAsciiChar *>(dac))
    ret=SWIG_NewPointerObj((void*)dac,SWIGTYPE_p_ParaMEDMEM__DataArrayAsciiChar,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of DataArray on downcast !");
  return ret;
}

static PyObject *convertIntArrToPyList(const int *ptr, int size) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=PyList_New(size);
  for(int i=0;i<size;i++)
    PyList_SetItem(ret,i,PyInt_FromLong(ptr[i]));
  return ret;
}

static PyObject *convertIntArrToPyList2(const std::vector<int>& v) throw(INTERP_KERNEL::Exception)
{
  int size=v.size();
  PyObject *ret=PyList_New(size);
  for(int i=0;i<size;i++)
    PyList_SetItem(ret,i,PyInt_FromLong(v[i]));
  return ret;
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
              throw INTERP_KERNEL::Exception("list must contain integers only");
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
              throw INTERP_KERNEL::Exception("tuple must contain integers only");
            }
        }
      return tmp;
    }
  else
    {
      throw INTERP_KERNEL::Exception("convertPyToNewIntArr2 : not a list");
    }
}

static PyObject *convertFromVectorPairInt(const std::vector< std::pair<int,int> >& arr) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=PyList_New(arr.size());
  for(std::size_t i=0;i<arr.size();i++)
    {
      PyObject *t=PyTuple_New(2);
      PyTuple_SetItem(t,0,PyInt_FromLong(arr[i].first));
      PyTuple_SetItem(t,1,PyInt_FromLong(arr[i].second));
      PyList_SetItem(ret,i,t);
    }
  return ret;
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
            throw INTERP_KERNEL::Exception("list must contain integers only");
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
            throw INTERP_KERNEL::Exception("tuple must contain integers only");
        }
    }
  else
    {
      throw INTERP_KERNEL::Exception("convertPyToNewIntArr3 : not a list nor a tuple");
    }
}

static void convertPyToNewIntArr4(PyObject *pyLi, int recurseLev, int nbOfSubPart, std::vector<int>& arr) throw(INTERP_KERNEL::Exception)
{
  if(recurseLev<0)
    throw INTERP_KERNEL::Exception("convertPyToNewIntArr4 : invalid list of integers level of recursion !");
  arr.clear();
  if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(PyInt_Check(o))
            {
              int val=(int)PyInt_AS_LONG(o);
              arr.push_back(val);
            }
          else
            {
              std::vector<int> arr2;
              convertPyToNewIntArr4(o,recurseLev-1,nbOfSubPart,arr2);
              if(nbOfSubPart>=1 && nbOfSubPart!=(int)arr2.size())
                  {
                    std::ostringstream oss; oss << "convertPyToNewIntArr4 : input list at lev " <<  recurseLev << " invalid nb of subpart elts expected " << nbOfSubPart << " having " << arr2.size() << " !";
                    throw INTERP_KERNEL::Exception(oss.str().c_str());
                  }
              arr.insert(arr.end(),arr2.begin(),arr2.end());
            }
        }
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
              arr.push_back(val);
            }
          else
            {
              std::vector<int> arr2;
              convertPyToNewIntArr4(o,recurseLev-1,nbOfSubPart,arr2);
              if(nbOfSubPart>=1 && nbOfSubPart!=(int)arr2.size())
                  {
                    std::ostringstream oss; oss << "convertPyToNewIntArr4 : input list at lev " <<  recurseLev << " invalid nb of subpart elts expected " << nbOfSubPart << " having " << arr2.size() << " !";
                    throw INTERP_KERNEL::Exception(oss.str().c_str());
                  }
              arr.insert(arr.end(),arr2.begin(),arr2.end());
            }
        }
    }
  else
    throw INTERP_KERNEL::Exception("convertPyToNewIntArr4 : not a list nor a tuple recursively !");
}

static void checkFillArrayWithPyList(int size1, int size2, int& nbOfTuples, int& nbOfComp) throw(INTERP_KERNEL::Exception)
{
  if(nbOfTuples==-1)
    {
      if(nbOfComp==-1) { nbOfTuples=size1; nbOfComp=size2; }
      else { if(nbOfComp==size2) { nbOfTuples=size1; } else
          {
            std::ostringstream oss; oss << "fillArrayWithPyListDbl2 : mismatch between nb of elemts : Input has " << size1 << " tuples and " << size2 << " components";
            oss << " whereas nb of components expected is " << nbOfComp << " !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          } }
    }
  else
    {
      if(nbOfComp!=-1)
        {
          if((nbOfTuples!=size1 || nbOfComp!=size2))
            {
              if(size2!=1 || size1!=nbOfComp*nbOfTuples)
                {
                  std::ostringstream oss; oss << "fillArrayWithPyListDbl2 : mismatch between nb of elemts : Input has " << size1 << " tuples and " << size2 << " components";
                  oss << " whereas nb of tuples expected is " << nbOfTuples << " and number of components expected is " << nbOfComp << " !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
            }
        }
      else
        {
          if(nbOfTuples==size1)
            nbOfComp=size2;
          else
            {
              std::ostringstream oss; oss << "fillArrayWithPyListDbl2 : mismatch between nb of elemts : Input has " << size1 << " tuples and " << size2 << " components";
              oss << " whereas nb of tuples expected is " << nbOfTuples << " !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
    }
}

static void fillArrayWithPyListInt3(PyObject *pyLi, int& nbOfElt, std::vector<int>& ret)
{
  static const char MSG[]="fillArrayWithPyListInt3 : It appears that the input list or tuple is composed by elts having different sizes !";
  if(PyInt_Check(pyLi))
    {
      long val=PyInt_AS_LONG(pyLi);
      if(nbOfElt==-1)
        nbOfElt=1;
      else
        if(nbOfElt!=1)
          throw INTERP_KERNEL::Exception(MSG);
      ret.push_back(val);
    }
  else if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      int tmp=0;
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          int tmp1=-1;
          fillArrayWithPyListInt3(o,tmp1,ret);
          tmp+=tmp1;
        }
      if(nbOfElt==-1)
        nbOfElt=tmp;
      else
        {
          if(nbOfElt!=tmp)
            throw INTERP_KERNEL::Exception(MSG);
        }
    }
  else if(PyTuple_Check(pyLi))
    {
      int size=PyTuple_Size(pyLi);
      int tmp=0;
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          int tmp1=-1;
          fillArrayWithPyListInt3(o,tmp1,ret);
          tmp+=tmp1;
        }
      if(nbOfElt==-1)
        nbOfElt=tmp;
      else
        {
          if(nbOfElt!=tmp)
            throw INTERP_KERNEL::Exception(MSG);
        }
    }
  else
    throw INTERP_KERNEL::Exception("fillArrayWithPyListInt3 : Unrecognized type ! Should be a composition of tuple,list,int !");
}

static std::vector<int> fillArrayWithPyListInt2(PyObject *pyLi, int& nbOfTuples, int& nbOfComp) throw(INTERP_KERNEL::Exception)
{
  std::vector<int> ret;
  int size1=-1,size2=-1;
  if(PyList_Check(pyLi))
    {
      size1=PyList_Size(pyLi);
      for(int i=0;i<size1;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          fillArrayWithPyListInt3(o,size2,ret);
        }
      if(size1==0)
        size2=1;
    }
  else if(PyTuple_Check(pyLi))
    {
      size1=PyTuple_Size(pyLi);
      for(int i=0;i<size1;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          fillArrayWithPyListInt3(o,size2,ret);
        }
      if(size1==0)
        size2=1;
    }
  else
    throw INTERP_KERNEL::Exception("fillArrayWithPyListInt2 : Unrecognized type ! Should be a tuple or a list !");
  //
  checkFillArrayWithPyList(size1,size2,nbOfTuples,nbOfComp);
  return ret;
}

static bool fillStringVector(PyObject *pyLi, std::vector<std::string>& vec) throw(INTERP_KERNEL::Exception)
{
  if(PyList_Check(pyLi))
    {
      Py_ssize_t sz=PyList_Size(pyLi);
      vec.resize(sz);
      for(int i=0;i<sz;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(PyString_Check(o))
            vec[i]=PyString_AsString(o);
          else
            return false;
        }
      return true;
    }
  else if(PyTuple_Check(pyLi))
    {
      Py_ssize_t sz=PyTuple_Size(pyLi);
      vec.resize(sz);
      for(int i=0;i<sz;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          if(PyString_Check(o))
            vec[i]=PyString_AsString(o);
          else
            return false;
        }
      return true;
    }
  else
    return false;
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

static PyObject *convertCharArrToPyListOfTuple(const char *vals, int nbOfComp, int nbOfTuples) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=PyList_New(nbOfTuples);
  INTERP_KERNEL::AutoPtr<char> tmp=new char[nbOfComp+1]; tmp[nbOfComp]='\0';
  for(int i=0;i<nbOfTuples;i++)
    {
      std::copy(vals+i*nbOfComp,vals+(i+1)*nbOfComp,(char *)tmp);
      PyList_SetItem(ret,i,PyString_FromString(tmp));
    }
  return ret;
}

static double *convertPyToNewDblArr2(PyObject *pyLi, int *size) throw(INTERP_KERNEL::Exception)
{
  if(PyList_Check(pyLi))
    {
      *size=PyList_Size(pyLi);
      double *tmp=(double *)malloc((*size)*sizeof(double));
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
              free(tmp);
              throw INTERP_KERNEL::Exception("convertPyToNewDblArr2 : list must contain floats/integers only");
            }
        }
      return tmp;
    }
  else if(PyTuple_Check(pyLi))
    {
      *size=PyTuple_Size(pyLi);
      double *tmp=(double *)malloc((*size)*sizeof(double));
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
              free(tmp);
              throw INTERP_KERNEL::Exception("convertPyToNewDblArr2 : tuple must contain floats/integers only");
            }
        }
      return tmp;
    }
  else
    throw INTERP_KERNEL::Exception("convertPyToNewDblArr2 : not a list");
}

static void fillArrayWithPyListDbl3(PyObject *pyLi, int& nbOfElt, std::vector<double>& ret)
{
  static const char MSG[]="fillArrayWithPyListDbl3 : It appears that the input list or tuple is composed by elts having different sizes !";
  if(PyFloat_Check(pyLi))
    {
      if(nbOfElt==-1)
        nbOfElt=1;
      else
        if(nbOfElt!=1)
          throw INTERP_KERNEL::Exception(MSG);
      double val=PyFloat_AS_DOUBLE(pyLi);
      ret.push_back(val);
    }
  else if(PyInt_Check(pyLi))
    {
      long val0=PyInt_AS_LONG(pyLi);
      double val=val0;
      if(nbOfElt==-1)
        nbOfElt=1;
      else
        if(nbOfElt!=1)
          throw INTERP_KERNEL::Exception(MSG);
      ret.push_back(val);
    }
  else if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      int tmp=0;
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          int tmp1=-1;
          fillArrayWithPyListDbl3(o,tmp1,ret);
          tmp+=tmp1;
        }
      if(nbOfElt==-1)
        nbOfElt=tmp;
      else
        {
          if(nbOfElt!=tmp)
            throw INTERP_KERNEL::Exception(MSG);
        }
    }
  else if(PyTuple_Check(pyLi))
    {
      int size=PyTuple_Size(pyLi);
      int tmp=0;
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          int tmp1=-1;
          fillArrayWithPyListDbl3(o,tmp1,ret);
          tmp+=tmp1;
        }
      if(nbOfElt==-1)
        nbOfElt=tmp;
      else
        {
          if(nbOfElt!=tmp)
            throw INTERP_KERNEL::Exception(MSG);
        }
    }
  else
    throw INTERP_KERNEL::Exception("fillArrayWithPyListDbl3 : Unrecognized type ! Should be a composition of tuple,list,int and float !");
}

static std::vector<double> fillArrayWithPyListDbl2(PyObject *pyLi, int& nbOfTuples, int& nbOfComp) throw(INTERP_KERNEL::Exception)
{
  std::vector<double> ret;
  int size1=-1,size2=-1;
  if(PyList_Check(pyLi))
    {
      size1=PyList_Size(pyLi);
      for(int i=0;i<size1;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          fillArrayWithPyListDbl3(o,size2,ret);
        }
      if(size1==0)
        size2=1;
    }
  else if(PyTuple_Check(pyLi))
    {
      size1=PyTuple_Size(pyLi);
      for(int i=0;i<size1;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          fillArrayWithPyListDbl3(o,size2,ret);
        }
      if(size1==0)
        size2=1;
    }
  else
    throw INTERP_KERNEL::Exception("fillArrayWithPyListDbl2 : Unrecognized type ! Should be a tuple or a list !");
  //
  checkFillArrayWithPyList(size1,size2,nbOfTuples,nbOfComp);
  return ret;
}

//convertFromPyObjVectorOfObj<const ParaMEDMEM::MEDCouplingUMesh *>(pyLi,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,"MEDCouplingUMesh")
template<class T>
static void convertFromPyObjVectorOfObj(PyObject *pyLi, swig_type_info *ty, const char *typeStr, typename std::vector<T>& ret)
{
  void *argp=0;
  if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      ret.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *obj=PyList_GetItem(pyLi,i);
          int status=SWIG_ConvertPtr(obj,&argp,ty,0|0);
          if(!SWIG_IsOK(status))
            {
              std::ostringstream oss; oss << "convertFromPyObjVectorOfObj : list is excepted to contain only " << typeStr << " instances !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          T arg=reinterpret_cast< T >(argp);
          ret[i]=arg;
        }
    }
  else if(PyTuple_Check(pyLi))
    {
      int size=PyTuple_Size(pyLi);
      ret.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *obj=PyTuple_GetItem(pyLi,i);
          int status=SWIG_ConvertPtr(obj,&argp,ty,0|0);
          if(!SWIG_IsOK(status))
            {
              std::ostringstream oss; oss << "convertFromPyObjVectorOfObj : tuple is excepted to contain only " << typeStr << " instances !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          T arg=reinterpret_cast< T >(argp);
          ret[i]=arg;
        }
    }
  else if(SWIG_IsOK(SWIG_ConvertPtr(pyLi,&argp,ty,0|0)))
    {
      ret.resize(1);
      T arg=reinterpret_cast< T >(argp);
      ret[0]=arg;
    }
  else
    throw INTERP_KERNEL::Exception("convertFromPyObjVectorOfObj : not a list nor a tuple");
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
 * if python slicp -> cpp pair sw=3 (begin,end,step)
 * if python DataArrayInt -> cpp DataArrayInt sw=4 . The returned pointer cannot be the null pointer ! If null an exception is thrown.
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
      Py_ssize_t strt=2,stp=2,step=2;
      PySliceObject *oC=reinterpret_cast<PySliceObject *>(value);
      GetIndicesOfSlice(oC,nbelem,&strt,&stp,&step,"Slice in subscriptable object DataArray invalid !");
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
  status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayIntTuple,0|0);
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

/*!
 * Idem than convertObjToPossibleCpp2
 */
static void convertObjToPossibleCpp2WithNegIntInterp(PyObject *value, int nbelem, int& sw, int& iTyypp, std::vector<int>& stdvecTyypp, std::pair<int, std::pair<int,int> >& p, ParaMEDMEM::DataArrayInt *& daIntTyypp) throw(INTERP_KERNEL::Exception)
{
  convertObjToPossibleCpp2(value,nbelem,sw,iTyypp,stdvecTyypp,p,daIntTyypp);
  if(sw==1)
    {
      iTyypp=InterpreteNegativeInt(iTyypp,nbelem);
    }
}

/*!
 * if python int -> cpp int sw=1
 * if python tuple[int] -> cpp vector<int> sw=2
 * if python list[int] -> cpp vector<int> sw=2
 * if python slice -> cpp pair sw=3
 * if python DataArrayIntTuple -> cpp DataArrayIntTuple sw=4 . WARNING The returned pointer can be the null pointer !
 */
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
      Py_ssize_t strt=2,stp=2,step=2;
      PySliceObject *oC=reinterpret_cast<PySliceObject *>(value);
      GetIndicesOfSlice(oC,nbelem,&strt,&stp,&step,"Slice in subscriptable object DataArray invalid !");
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
 * if python string with size one -> cpp char sw=1
 * if python string with size different from one -> cpp string sw=2
 * if python tuple[string] or list[string] -> vector<string> sw=3
 * if python not null pointer of DataArrayChar -> cpp DataArrayChar sw=4
 * switch between (int,string,vector<string>,DataArrayChar)
 */
static void convertObjToPossibleCpp6(PyObject *value, int& sw, char& cTyp, std::string& sType, std::vector<std::string>& vsType, ParaMEDMEM::DataArrayChar *& dacType) throw(INTERP_KERNEL::Exception)
{
  const char *msg="4 types accepted : string, list or tuple of strings having same size, not null DataArrayChar instance.";
  sw=-1;
  if(PyString_Check(value))
    {
      const char *pt=PyString_AsString(value);
      Py_ssize_t sz=PyString_Size(value);
      if(sz==1)
        {
          cTyp=pt[0];
          sw=1;
          return;
        }
      else
        {
          sType=pt;
          sw=2;
          return;
        }
    }
  if(PyTuple_Check(value))
    {
      int size=PyTuple_Size(value);
      vsType.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(value,i);
          if(PyString_Check(o))
            vsType[i]=PyString_AsString(o);
          else
            {
              std::ostringstream oss; oss << "Tuple as been detected but element #" << i << " is not a string ! only tuples of strings accepted !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      sw=3;
      return;
    }
  if(PyList_Check(value))
    {
      int size=PyList_Size(value);
      vsType.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(value,i);
          if(PyString_Check(o))
            vsType[i]=PyString_AsString(o);
          else
            {
              std::ostringstream oss; oss << "List as been detected but element #" << i << " is not string ! only lists of strings accepted !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      sw=3;
      return;
    }
  void *argp;
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayChar,0|0);
  if(SWIG_IsOK(status))
    {
      dacType=reinterpret_cast< ParaMEDMEM::DataArrayChar * >(argp);
      if(!dacType)
        {
          std::ostringstream oss; oss << msg << " Instance in null !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      sw=4;
      return ;
    }
  throw INTERP_KERNEL::Exception(msg);
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
      convertObjToPossibleCpp2WithNegIntInterp(value,nbTuple,sw,it,vt,pt,dt);
      return ;
    }
  else
    {
      int sz=PyTuple_Size(value);
      if(sz!=2)
        throw INTERP_KERNEL::Exception("Unexpected nb of slice element : 1 or 2 expected !\n1st is for tuple selection, 2nd for component selection !");
      PyObject *ob0=PyTuple_GetItem(value,0);
      int sw1,sw2;
      convertObjToPossibleCpp2WithNegIntInterp(ob0,nbTuple,sw1,it,vt,pt,dt);
      PyObject *ob1=PyTuple_GetItem(value,1);
      convertObjToPossibleCpp2WithNegIntInterp(ob1,nbCompo,sw2,ic,vc,pc,dc);
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

/*!
 * if value int -> cpp val sw=1
 * if value double -> cpp val sw=1
 * if value DataArrayDouble -> cpp DataArrayDouble sw=2
 * if value DataArrayDoubleTuple -> cpp DataArrayDoubleTuple sw=3
 * if value list[int,double] -> cpp std::vector<double> sw=4
 * if value tuple[int,double] -> cpp std::vector<double> sw=4
 */
static const double *convertObjToPossibleCpp5_Safe(PyObject *value, int& sw, double& val, ParaMEDMEM::DataArrayDouble *&d, ParaMEDMEM::DataArrayDoubleTuple *&e, std::vector<double>& f,
                                                   const char *msg, int nbTuplesExpected, int nbCompExpected, bool throwIfNullPt) throw(INTERP_KERNEL::Exception)
{
  sw=-1;
  if(PyFloat_Check(value))
    {
      val=PyFloat_AS_DOUBLE(value);
      sw=1;
      if(nbTuplesExpected*nbCompExpected!=1)
        {
          std::ostringstream oss; oss << msg << "dimension expected to be " << nbTuplesExpected*nbCompExpected << " , and your data in input has dimension one (single PyFloat) !"; 
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      return &val;
    }
  if(PyInt_Check(value))
    {
      val=(double)PyInt_AS_LONG(value);
      sw=1;
      if(nbTuplesExpected*nbCompExpected!=1)
        {
          std::ostringstream oss; oss << msg << "dimension expected to be " << nbTuplesExpected*nbCompExpected << " , and your data in input has dimension one (single PyInt) !"; 
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      return &val;
    }
  if(PyTuple_Check(value) || PyList_Check(value))
    {
      try
        {
          int tmp1=nbTuplesExpected,tmp2=nbCompExpected;
          std::vector<double> ret=fillArrayWithPyListDbl2(value,tmp1,tmp2);
          sw=4;
          f=ret;
          return &f[0];
        }
      catch(INTERP_KERNEL::Exception& exc) { throw exc; }
    }
  void *argp;
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble,0|0);
  if(SWIG_IsOK(status))
    {  
      d=reinterpret_cast< ParaMEDMEM::DataArrayDouble * >(argp);
      sw=2;
      if(d)
        {
          if(d->getNumberOfTuples()==nbTuplesExpected)
            {
              if(d->getNumberOfComponents()==nbCompExpected)
                {
                  return d->getConstPointer();
                }
              else
                {
                  std::ostringstream oss; oss << msg << "nb of components expected to be " << nbCompExpected << " , and input has " << d->getNumberOfComponents() << " components !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
            }
          else
            {
              std::ostringstream oss; oss << msg << " input DataArrayDouble should have a number of tuples equal to " << nbTuplesExpected << " and there are " << d->getNumberOfTuples() << " tuples !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          if(throwIfNullPt)
            {
              std::ostringstream oss; oss << msg << " null pointer not accepted!";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          else
            return 0;
        }
    }
  status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayDoubleTuple,0|0);
  if(SWIG_IsOK(status))
    {  
      e=reinterpret_cast< ParaMEDMEM::DataArrayDoubleTuple * >(argp);
      sw=3;
      if(e->getNumberOfCompo()==nbCompExpected)
        {
          if(nbTuplesExpected==1)
            return e->getConstPointer();
          else
            {
              std::ostringstream oss; oss << msg << "nb of tuples expected to be " << nbTuplesExpected << " , and input DataArrayDoubleTuple has always one tuple by contruction !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << msg << "nb of components expected to be " <<  nbCompExpected << " , and input DataArrayDoubleTuple has " << e->getNumberOfCompo() << " components !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  throw INTERP_KERNEL::Exception("4 types accepted : integer, double, DataArrayDouble, DataArrayDoubleTuple");
}

/*!
 * if value int -> cpp val sw=1
 * if value double -> cpp val sw=1
 * if value DataArrayDouble -> cpp DataArrayDouble sw=2
 * if value DataArrayDoubleTuple -> cpp DataArrayDoubleTuple sw=3
 * if value list[int,double] -> cpp std::vector<double> sw=4
 * if value tuple[int,double] -> cpp std::vector<double> sw=4
 */
static const double *convertObjToPossibleCpp5_Safe2(PyObject *value, int& sw, double& val, ParaMEDMEM::DataArrayDouble *&d, ParaMEDMEM::DataArrayDoubleTuple *&e, std::vector<double>& f,
                                                    const char *msg, int nbCompExpected, bool throwIfNullPt, int& nbTuples) throw(INTERP_KERNEL::Exception)
{
  sw=-1;
  if(PyFloat_Check(value))
    {
      val=PyFloat_AS_DOUBLE(value);
      sw=1;
      if(nbCompExpected!=1)
        {
          std::ostringstream oss; oss << msg << "dimension expected to be " << nbCompExpected << " , and your data in input has dimension one (single PyFloat) !"; 
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      nbTuples=1;
      return &val;
    }
  if(PyInt_Check(value))
    {
      val=(double)PyInt_AS_LONG(value);
      sw=1;
      if(nbCompExpected!=1)
        {
          std::ostringstream oss; oss << msg << "dimension expected to be " << nbCompExpected << " , and your data in input has dimension one (single PyInt) !"; 
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      nbTuples=1;
      return &val;
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
      if(size%nbCompExpected!=0)
        {
          std::ostringstream oss; oss << msg << "dimension expected to be a multiple of " << nbCompExpected << " , and your data in input has dimension " << f.size() << " !"; 
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      nbTuples=size/nbCompExpected;
      return &f[0];
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
      if(size%nbCompExpected!=0)
        {
          std::ostringstream oss; oss << msg << "dimension expected to be a multiple of " << nbCompExpected << " , and your data in input has dimension " << f.size() << " !"; 
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      nbTuples=size/nbCompExpected;
      return &f[0];
    }
  void *argp;
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble,0|0);
  if(SWIG_IsOK(status))
    {  
      d=reinterpret_cast< ParaMEDMEM::DataArrayDouble * >(argp);
      sw=2;
      if(d)
        {
          if(d->getNumberOfComponents()==nbCompExpected)
            {
              nbTuples=d->getNumberOfTuples();
              return d->getConstPointer();
            }
          else
            {
              std::ostringstream oss; oss << msg << "nb of components expected to be a multiple of " << nbCompExpected << " , and input has " << d->getNumberOfComponents() << " components !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          if(throwIfNullPt)
            {
              std::ostringstream oss; oss << msg << " null pointer not accepted!";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          else
            { nbTuples=0; return 0; }
        }
    }
  status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayDoubleTuple,0|0);
  if(SWIG_IsOK(status))
    {  
      e=reinterpret_cast< ParaMEDMEM::DataArrayDoubleTuple * >(argp);
      sw=3;
      if(e)
        {
          if(e->getNumberOfCompo()==nbCompExpected)
            {
              nbTuples=1;
              return e->getConstPointer();
            }
          else
            {
              std::ostringstream oss; oss << msg << "nb of components expected to be " <<  nbCompExpected << " , and input DataArrayDoubleTuple has " << e->getNumberOfCompo() << " components !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          if(throwIfNullPt)
            {
              std::ostringstream oss; oss << msg << " null pointer not accepted!";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          else
            { nbTuples=0; return 0; }
        }
    }
  throw INTERP_KERNEL::Exception("4 types accepted : integer, double, DataArrayDouble, DataArrayDoubleTuple");
}

/*!
 * if value int -> cpp val sw=1
 * if value double -> cpp val sw=1
 * if value DataArrayDouble -> cpp DataArrayDouble sw=2
 * if value DataArrayDoubleTuple -> cpp DataArrayDoubleTuple sw=3
 * if value list[int,double] -> cpp std::vector<double> sw=4
 * if value tuple[int,double] -> cpp std::vector<double> sw=4
 */
static const double *convertObjToPossibleCpp5_SingleCompo(PyObject *value, int& sw, double& val, std::vector<double>& f,
                                                          const char *msg, bool throwIfNullPt, int& nbTuples) throw(INTERP_KERNEL::Exception)
{
  ParaMEDMEM::DataArrayDouble *d=0;
  ParaMEDMEM::DataArrayDoubleTuple *e=0;
  sw=-1;
  if(PyFloat_Check(value))
    {
      val=PyFloat_AS_DOUBLE(value);
      sw=1;
      nbTuples=1;
      return &val;
    }
  if(PyInt_Check(value))
    {
      val=(double)PyInt_AS_LONG(value);
      sw=1;
      nbTuples=1;
      return &val;
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
      nbTuples=size;
      return &f[0];
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
      nbTuples=size;
      return &f[0];
    }
  void *argp;
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble,0|0);
  if(SWIG_IsOK(status))
    {  
      d=reinterpret_cast< ParaMEDMEM::DataArrayDouble * >(argp);
      sw=2;
      if(d)
        {
          if(d->getNumberOfComponents()==1)
            {
              nbTuples=d->getNumberOfTuples();
              return d->getConstPointer();
            }
          else
            {
              std::ostringstream oss; oss << msg << "nb of components expected to be one, and input has " << d->getNumberOfComponents() << " components !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          if(throwIfNullPt)
            {
              std::ostringstream oss; oss << msg << " null pointer not accepted!";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          else
            { nbTuples=0; return 0; }
        }
    }
  status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayDoubleTuple,0|0);
  if(SWIG_IsOK(status))
    {  
      e=reinterpret_cast< ParaMEDMEM::DataArrayDoubleTuple * >(argp);
      sw=3;
      if(e)
        {
          nbTuples=e->getNumberOfCompo();
          return e->getConstPointer();
        }
      else
        {
          if(throwIfNullPt)
            {
              std::ostringstream oss; oss << msg << " null pointer not accepted!";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          else
            { nbTuples=0; return 0; }
        }
    }
  throw INTERP_KERNEL::Exception("4 types accepted : integer, double, DataArrayDouble, DataArrayDoubleTuple");
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
static const int *convertObjToPossibleCpp1_Safe(PyObject *value, int& sw, int& sz, int& iTyypp, std::vector<int>& stdvecTyypp) throw(INTERP_KERNEL::Exception)
{
  sw=-1;
  if(PyInt_Check(value))
    {
      iTyypp=(int)PyInt_AS_LONG(value);
      sw=1; sz=1;
      return &iTyypp;
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
      sw=2; sz=size;
      return &stdvecTyypp[0];
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
      sw=2; sz=size;
      return &stdvecTyypp[0];
    }
  void *argp;
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,0|0);
  if(SWIG_IsOK(status))
    {
      ParaMEDMEM::DataArrayInt *daIntTyypp=reinterpret_cast< ParaMEDMEM::DataArrayInt * >(argp);
      if(daIntTyypp)
        {
          sw=3; sz=daIntTyypp->getNbOfElems();
          return daIntTyypp->begin();
        }
      else
        {
          sz=0;
          return 0;
        }
    }
  status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_ParaMEDMEM__DataArrayIntTuple,0|0);
  if(SWIG_IsOK(status))
    {  
      ParaMEDMEM::DataArrayIntTuple *daIntTuple=reinterpret_cast< ParaMEDMEM::DataArrayIntTuple * >(argp);
      sw=4; sz=daIntTuple->getNumberOfCompo();
      return daIntTuple->getConstPointer();
    }
  throw INTERP_KERNEL::Exception("5 types accepted : integer, tuple of integer, list of integer, DataArrayInt, DataArrayIntTuple");
}

static ParaMEDMEM::DataArray *CheckAndRetrieveDataArrayInstance(PyObject *obj, const char *msg)
{
  void *aBasePtrVS=0;
  int status=SWIG_ConvertPtr(obj,&aBasePtrVS,SWIGTYPE_p_ParaMEDMEM__DataArray,0|0);
  if(!SWIG_IsOK(status))
    {
      status=SWIG_ConvertPtr(obj,&aBasePtrVS,SWIGTYPE_p_ParaMEDMEM__DataArrayDouble,0|0);
      if(!SWIG_IsOK(status))
        {
          status=SWIG_ConvertPtr(obj,&aBasePtrVS,SWIGTYPE_p_ParaMEDMEM__DataArrayInt,0|0);
          if(!SWIG_IsOK(status))
            {
              status=SWIG_ConvertPtr(obj,&aBasePtrVS,SWIGTYPE_p_ParaMEDMEM__DataArrayAsciiChar,0|0);
              if(!SWIG_IsOK(status))
                {
                  status=SWIG_ConvertPtr(obj,&aBasePtrVS,SWIGTYPE_p_ParaMEDMEM__DataArrayByte,0|0);
                  std::ostringstream oss; oss << msg << " ! Accepted instances are DataArrayDouble, DataArrayInt, DataArrayAsciiChar, DataArrayByte !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
            }
        }
    }
  return reinterpret_cast< ParaMEDMEM::DataArray * >(aBasePtrVS);
}
