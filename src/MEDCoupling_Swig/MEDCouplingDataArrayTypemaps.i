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

#ifndef __MEDCOUPLINGDATAARRAYTYPEMAPS_I__
#define __MEDCOUPLINGDATAARRAYTYPEMAPS_I__

#if PY_VERSION_HEX >= 0x03000000
#define PyInt_AS_LONG PyLong_AS_LONG
#endif

#include "InterpKernelAutoPtr.hxx"
#include "MEDCouplingDataArrayTraits.hxx"

#include <sstream>

static PyObject *convertArray(MEDCoupling::DataArray *array, int owner)
{
  PyObject *ret(NULL);
  if(!array)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDCoupling::DataArrayDouble *>(array))
    ret=SWIG_NewPointerObj((void*)array,SWIGTYPE_p_MEDCoupling__DataArrayDouble,owner);
  if(dynamic_cast<MEDCoupling::DataArrayInt *>(array))
    ret=SWIG_NewPointerObj((void*)array,SWIGTYPE_p_MEDCoupling__DataArrayInt,owner);
  if(dynamic_cast<MEDCoupling::DataArrayFloat *>(array))
    ret=SWIG_NewPointerObj((void*)array,SWIGTYPE_p_MEDCoupling__DataArrayFloat,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of array on downcast !");
  return ret;
}

/*!
 * This method is an extention of PySlice_GetIndices but less
 * open than PySlice_GetIndicesEx that accepts too many situations.
 */
void GetIndicesOfSlice(PyObject *slice, Py_ssize_t length, Py_ssize_t *start, Py_ssize_t *stop, Py_ssize_t *step, const char *msgInCaseOfFailure)
{
  int ret(PySlice_GetIndices(
#if PY_VERSION_HEX >= 0x03000000
        slice,
#else
        reinterpret_cast<PySliceObject *>(slice),
#endif
        length,start,stop,step));
  if(ret==0)
    return ;
  if(*step>0 && *start==*stop && length==*start)
    return ;
  throw INTERP_KERNEL::Exception(msgInCaseOfFailure);
}

/*!
 * This method allows to retrieve slice info from \a slice.
 */
void GetIndicesOfSliceExplicitely(PyObject *slice, Py_ssize_t *start, Py_ssize_t *stop, Py_ssize_t *step, const char *msgInCaseOfFailure)
{
  int ret(PySlice_GetIndices(
#if PY_VERSION_HEX >= 0x03000000
        slice,
#else
        reinterpret_cast<PySliceObject *>(slice),
#endif
        std::numeric_limits<int>::max(),start,stop,step));
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
  typename MEDCoupling::MCAuto<MCData> ret=MCData::New();
  if(PyArray_ISBEHAVED(elt0))//aligned and writeable and in machine byte-order
    {
      PyArrayObject *elt0C=reinterpret_cast<PyArrayObject *>(elt0);
      PyArrayObject *eltOwning=(PyArray_FLAGS(elt0C) & MED_NUMPY_OWNDATA)?elt0C:NULL;
      int mask=MED_NUMPY_OWNDATA; mask=~mask;
      elt0C->flags&=mask;
      PyObject *deepestObj=elt0;
      PyObject *base=elt0C->base;
      if(base) deepestObj=base;
      bool isSpetialCase(false);
      while(base)
        {
          if(PyArray_Check(base))
            {
              PyArrayObject *baseC=reinterpret_cast<PyArrayObject *>(base);
              eltOwning=(PyArray_FLAGS(baseC) & MED_NUMPY_OWNDATA)?baseC:eltOwning;
              baseC->flags&=mask;
              base=baseC->base;
              if(base) deepestObj=base;
            }
          else
            {
              isSpetialCase=true;
              break;
            }
        }
      if(isSpetialCase)
        {// this case is present for numpy arrayint coming from load of pickelized string. The owner of elt0 is not an array -> A copy is requested.
          std::size_t nbOfElems(sz0*sz1);
          T *dataCpy=(T*)malloc(sizeof(T)*nbOfElems);
          std::copy(reinterpret_cast<const T*>(data),reinterpret_cast<const T*>(data)+nbOfElems,dataCpy);
          ret->useArray(dataCpy,true,MEDCoupling::C_DEALLOC,sz0,sz1);
          return ret.retn();
        }
      typename MEDCoupling::MemArray<T>& mma=ret->accessToMemArray();
      if(eltOwning==NULL)
        {
          PyCallBackDataArraySt<MCData> *cb=PyObject_GC_New(PyCallBackDataArraySt<MCData>,pytype);
          cb->_pt_mc=ret;
          ret->useArray(reinterpret_cast<const T *>(data),true,MEDCoupling::C_DEALLOC,sz0,sz1);
          PyObject *ref=PyWeakref_NewRef(deepestObj,(PyObject *)cb);
          void **objs=new void *[2]; objs[0]=cb; objs[1]=ref;
          mma.setParameterForDeallocator(objs);
          mma.setSpecificDeallocator(numarrdeal2<MCData>);
          //"Impossible to share this numpy array chunk of data, because already shared by an another non numpy array object (maybe an another DataArrayInt instance) ! Release it, or perform a copy on the input array !");
        }
      else
        {
          ret->useArray(reinterpret_cast<const T *>(data),true,MEDCoupling::C_DEALLOC,sz0,sz1);
          PyObject *ref=PyWeakref_NewRef(reinterpret_cast<PyObject *>(eltOwning),NULL);
          typename MEDCoupling::MemArray<T>::Deallocator tmp(MEDCoupling::MemArray<T>::CDeallocator);
          void **tmp2 = reinterpret_cast<void**>(&tmp); // MSVC2010 does not support constructor()
          void **objs=new void *[2]; objs[0]=ref; objs[1]=*tmp2;
          mma.setParameterForDeallocator(objs);
          mma.setSpecificDeallocator(numarrdeal);
        }
    }
  else if(PyArray_ISBEHAVED_RO(elt0))
    ret->useArray(reinterpret_cast<const T *>(data),false,MEDCoupling::CPP_DEALLOC,sz0,sz1);
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
        if (PyArray_CHKFLAGS(obj_arr, MED_NUMPY_OWNDATA )) {
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
PyObject *ToNumPyArrayUnderground(MCData *self, int npyObjectType, const char *MCDataStr, int nbTuples, int nbComp)
{
  if(!self->isAllocated())
    {
      std::ostringstream oss; oss << MCDataStr << "::toNumPyArray : this is not allocated !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDCoupling::MemArray<T>& mem=self->accessToMemArray();
  if(nbComp==0)
    {
      std::ostringstream oss; oss << MCDataStr << "::toNumPyArray : number of components of this is 0 ! Should be > 0 !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  int nbDims=nbComp==1?1:2;
  npy_intp dim[2];
  dim[0]=(npy_intp)nbTuples; dim[1]=nbComp;
  const T *bg=self->getConstPointer();
  PyObject *ret(PyArray_SimpleNewFromData(nbDims,dim,npyObjectType,const_cast<T *>(bg)));
  if(mem.isDeallocatorCalled())
    {
      if(mem.getDeallocator()!=numarrdeal)
        {// case for the first call of toNumPyArray
          PyObject *ref(PyWeakref_NewRef(ret,NULL));
          typename MEDCoupling::MemArray<T>::Deallocator tmp(mem.getDeallocator());
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

template<class MCData, class T>
PyObject *ToNumPyArray(MCData *self, int npyObjectType, const char *MCDataStr)
{
  return ToNumPyArrayUnderground<MCData,T>(self,npyObjectType,MCDataStr,self->getNumberOfTuples(),self->getNumberOfComponents());
}

SWIGINTERN PyObject *MEDCoupling_DataArrayInt_toNumPyArray(MEDCoupling::DataArrayInt *self);
SWIGINTERN PyObject *MEDCoupling_DataArrayDouble_toNumPyArray(MEDCoupling::DataArrayDouble *self);

#endif

#ifdef WITH_SCIPY
PyObject *ToCSRMatrix(const std::vector<std::map<int,double> >& m, int nbCols)
{
  int nbRows((int)m.size());
  MEDCoupling::MCAuto<MEDCoupling::DataArrayInt> indPtr(MEDCoupling::DataArrayInt::New()),indices(MEDCoupling::DataArrayInt::New());
  MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> data(MEDCoupling::DataArrayDouble::New());
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
  PyObject *a(MEDCoupling_DataArrayDouble_toNumPyArray(data)),*b(MEDCoupling_DataArrayInt_toNumPyArray(indices)),*c(MEDCoupling_DataArrayInt_toNumPyArray(indPtr));
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

static PyObject *convertDataArrayChar(MEDCoupling::DataArrayChar *dac, int owner)
{
  PyObject *ret=0;
  if(!dac)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDCoupling::DataArrayByte *>(dac))
    ret=SWIG_NewPointerObj((void*)dac,SWIGTYPE_p_MEDCoupling__DataArrayByte,owner);
  if(dynamic_cast<MEDCoupling::DataArrayAsciiChar *>(dac))
    ret=SWIG_NewPointerObj((void*)dac,SWIGTYPE_p_MEDCoupling__DataArrayAsciiChar,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of DataArrayChar on downcast !");
  return ret;
}

static PyObject *convertDataArray(MEDCoupling::DataArray *dac, int owner)
{
  PyObject *ret=0;
  if(!dac)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDCoupling::DataArrayDouble *>(dac))
    ret=SWIG_NewPointerObj((void*)dac,SWIGTYPE_p_MEDCoupling__DataArrayDouble,owner);
  if(dynamic_cast<MEDCoupling::DataArrayInt *>(dac))
    ret=SWIG_NewPointerObj((void*)dac,SWIGTYPE_p_MEDCoupling__DataArrayInt,owner);
  if(dynamic_cast<MEDCoupling::DataArrayFloat *>(dac))
    ret=SWIG_NewPointerObj((void*)dac,SWIGTYPE_p_MEDCoupling__DataArrayFloat,owner);
  if(dynamic_cast<MEDCoupling::DataArrayByte *>(dac))
    ret=SWIG_NewPointerObj((void*)dac,SWIGTYPE_p_MEDCoupling__DataArrayByte,owner);
  if(dynamic_cast<MEDCoupling::DataArrayAsciiChar *>(dac))
    ret=SWIG_NewPointerObj((void*)dac,SWIGTYPE_p_MEDCoupling__DataArrayAsciiChar,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of DataArray on downcast !");
  return ret;
}

static PyObject *convertIntArrToPyList(const int *ptr, int size)
{
  PyObject *ret=PyList_New(size);
  for(int i=0;i<size;i++)
    PyList_SetItem(ret,i,PyInt_FromLong(ptr[i]));
  return ret;
}

static PyObject *convertIntArrToPyList2(const std::vector<int>& v)
{
  int size=v.size();
  PyObject *ret=PyList_New(size);
  for(int i=0;i<size;i++)
    PyList_SetItem(ret,i,PyInt_FromLong(v[i]));
  return ret;
}

static PyObject *convertIntArrToPyList3(const std::set<int>& v)
{
  int size=v.size();
  PyObject *ret=PyList_New(size);
  std::set<int>::const_iterator it=v.begin();
  for(int i=0;i<size;i++,it++)
    PyList_SetItem(ret,i,PyInt_FromLong(*it));
  return ret;
}

static bool convertPyObjectToStrNT(PyObject *obj, std::string& ret)
{
  if(PyString_Check(obj))
    {
      ret=PyString_AsString(obj);
      return true;
    }
#if PY_VERSION_HEX >= 0x03000000
  else if(PyUnicode_Check(obj))
    {
      ret=PyUnicode_AsUTF8(obj);
      return true;
    }
#endif
  return false;
}

static std::string convertPyObjectToStr(PyObject *obj, const char *msg=NULL)
{
  std::string ret;
  if(PyString_Check(obj))
    ret=PyString_AsString(obj);
#if PY_VERSION_HEX >= 0x03000000
  else if(PyUnicode_Check(obj))
    ret=PyUnicode_AsUTF8(obj);
#endif
  else
    {
      std::ostringstream oss;
      if(msg)
        oss << msg;
      else
        oss << "PyWrap convertPyObjectToStr : expect a string like py object !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  return ret;
}

static PyObject *convertIntArrToPyListOfTuple(const int *vals, int nbOfComp, int nbOfTuples)
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

static int *convertPyToNewIntArr2(PyObject *pyLi, int *size)
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

static PyObject *convertFromVectorPairInt(const std::vector< std::pair<int,int> >& arr)
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

static void convertPyToVectorPairInt(PyObject *pyLi, std::vector< std::pair<int,int> >& arr)
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

static void convertPyToVectorPairStringInt(PyObject *pyLi, std::vector< std::pair<std::string,int> >& arr)
{
  const char msg[]="convertPyToVectorPairStringInt : list must contain tuples of 2 integers only or tuple must contain tuples of 1 string and 1 integer only !";
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
              PyObject *o_1=PyTuple_GetItem(o,1);
              arr[i].first=convertPyObjectToStr(o_0,msg);
              if(!PyInt_Check(o_1))
                throw INTERP_KERNEL::Exception(msg);
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
              PyObject *o_1=PyTuple_GetItem(o,1);
              arr[i].first=convertPyObjectToStr(o_0,msg);
              if(!PyInt_Check(o_1))
                throw INTERP_KERNEL::Exception(msg);
              arr[i].second=(int)PyInt_AS_LONG(o_1);
            }
          else
            throw INTERP_KERNEL::Exception(msg);
        }
    }
  else
    throw INTERP_KERNEL::Exception(msg);
}

static void convertPyToNewIntArr3(PyObject *pyLi, std::vector<int>& arr)
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

static void convertPyToNewIntArr4(PyObject *pyLi, int recurseLev, int nbOfSubPart, std::vector<int>& arr)
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

static void checkFillArrayWithPyList(int size1, int size2, int& nbOfTuples, int& nbOfComp)
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

static std::vector<int> fillArrayWithPyListInt2(PyObject *pyLi, int& nbOfTuples, int& nbOfComp)
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

static bool fillStringVector(PyObject *pyLi, std::vector<std::string>& vec)
{
  if(PyList_Check(pyLi))
    {
      Py_ssize_t sz=PyList_Size(pyLi);
      vec.resize(sz);
      for(int i=0;i<sz;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(!convertPyObjectToStrNT(o,vec[i]))
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
          if(!convertPyObjectToStrNT(o,vec[i]))
            return false;
        }
      return true;
    }
  else
    return false;
}
static void convertPyToVectorOfVectorOfString(PyObject *pyLi, std::vector< std::vector<std::string> >& arr)
{
  const char msg[]="convertPyToVectorOfVectorOfString : expecting list of list of strings !";
  if(PyList_Check(pyLi))
    {
      Py_ssize_t sz=PyList_Size(pyLi);
      arr.resize(sz);
      for(int i=0;i<sz;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(!fillStringVector(o,arr[i]))
            throw INTERP_KERNEL::Exception(msg);
        }
    }
  else if(PyTuple_Check(pyLi))
    {
      Py_ssize_t sz=PyTuple_Size(pyLi);
      arr.resize(sz);
      for(int i=0;i<sz;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          if(!fillStringVector(o,arr[i]))
            throw INTERP_KERNEL::Exception(msg);
        }
    }
  else
    throw INTERP_KERNEL::Exception(msg);
}

static bool fillIntVector(PyObject *pyLi, std::vector<int>& vec)
{
  if(PyList_Check(pyLi))
    {
      Py_ssize_t sz=PyList_Size(pyLi);
      vec.resize(sz);
      for(int i=0;i<sz;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(PyInt_Check(o))
            vec[i]=PyInt_AS_LONG(o);
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
          if(PyInt_Check(o))
            vec[i]=PyInt_AS_LONG(o);
          else
            return false;
        }
      return true;
    }
  else
    return false;
}

static void convertPyToVectorOfVectorOfInt(PyObject *pyLi, std::vector< std::vector<int> >& arr)
{
  const char msg[]="convertPyToVectorOfVectorOfInt : expecting list of list of strings !";
  if(PyList_Check(pyLi))
    {
      Py_ssize_t sz=PyList_Size(pyLi);
      arr.resize(sz);
      for(int i=0;i<sz;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(!fillIntVector(o,arr[i]))
            throw INTERP_KERNEL::Exception(msg);
        }
    }
  else if(PyTuple_Check(pyLi))
    {
      Py_ssize_t sz=PyTuple_Size(pyLi);
      arr.resize(sz);
      for(int i=0;i<sz;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          if(!fillIntVector(o,arr[i]))
            throw INTERP_KERNEL::Exception(msg);
        }
    }
  else
    throw INTERP_KERNEL::Exception(msg);
}

static void convertPyToVectorPairStringVecString(PyObject *pyLi, std::vector< std::pair<std::string, std::vector<std::string> > >& arr)
{
  const char msg[]="convertPyToVectorPairStringVecString : expecting list of tuples containing each exactly 2 items : one string and one vector of string !";
  if(PyList_Check(pyLi))
    {
      Py_ssize_t sz=PyList_Size(pyLi);
      arr.resize(sz);
      for(int i=0;i<sz;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(PyTuple_Check(o))
            {
              int sz2=PyTuple_Size(o);
              if(sz2!=2)
                throw INTERP_KERNEL::Exception(msg);
              std::pair<std::string, std::vector<std::string> > item;
              PyObject *o_0=PyTuple_GetItem(o,0);
              item.first=convertPyObjectToStr(o_0,msg);
              PyObject *o_1=PyTuple_GetItem(o,1);
              if(!fillStringVector(o_1,item.second))
                throw INTERP_KERNEL::Exception(msg);
              arr[i]=item;
            }
          else
            throw INTERP_KERNEL::Exception(msg);
        }
    }
  else if(PyTuple_Check(pyLi))
    {
      Py_ssize_t sz=PyTuple_Size(pyLi);
      arr.resize(sz);
      for(int i=0;i<sz;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          if(PyTuple_Check(o))
            {
              int sz2=PyTuple_Size(o);
              if(sz2!=2)
                throw INTERP_KERNEL::Exception(msg);
              std::pair<std::string, std::vector<std::string> > item;
              PyObject *o_0=PyTuple_GetItem(o,0);
              item.first=convertPyObjectToStr(o_0,msg);
              PyObject *o_1=PyTuple_GetItem(o,1);
              if(!fillStringVector(o_1,item.second))
                throw INTERP_KERNEL::Exception(msg);
              arr[i]=item;
            }
          else
            throw INTERP_KERNEL::Exception(msg);
        }
    }
  else
    throw INTERP_KERNEL::Exception(msg);
}

template<class T>
PyObject *convertDblArrToPyList(const T *ptr, int size)
{
  PyObject *ret(PyList_New(size));
  for(int i=0;i<size;i++)
    PyList_SetItem(ret,i,PyFloat_FromDouble(ptr[i]));
  return ret;
}

static PyObject *convertDblArrToPyList2(const std::vector<double>& v)
{
  int size(v.size());
  PyObject *ret(PyList_New(size));
  for(int i=0;i<size;i++)
    PyList_SetItem(ret,i,PyFloat_FromDouble(v[i]));
  return ret;
}

template<class T>
PyObject *convertDblArrToPyListOfTuple(const T *vals, int nbOfComp, int nbOfTuples)
{
  PyObject *ret(PyList_New(nbOfTuples));
  for(int i=0;i<nbOfTuples;i++)
    {
      PyObject *t=PyTuple_New(nbOfComp);
      for(int j=0;j<nbOfComp;j++)
        PyTuple_SetItem(t,j,PyFloat_FromDouble(vals[i*nbOfComp+j]));
      PyList_SetItem(ret,i,t);
    }
  return ret;
}

static PyObject *convertCharArrToPyListOfTuple(const char *vals, int nbOfComp, int nbOfTuples)
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

static double *convertPyToNewDblArr2(PyObject *pyLi, int *size)
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

static std::vector<double> fillArrayWithPyListDbl2(PyObject *pyLi, int& nbOfTuples, int& nbOfComp)
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

//convertFromPyObjVectorOfObj<const MEDCoupling::MEDCouplingUMesh *>(pyLi,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,"MEDCouplingUMesh")
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
static void convertIntStarLikePyObjToCpp(PyObject *value, int& sw, int& iTyypp, std::vector<int>& stdvecTyypp, MEDCoupling::DataArrayInt *& daIntTyypp, MEDCoupling::DataArrayIntTuple *&daIntTuple)
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
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayInt,0|0);
  if(SWIG_IsOK(status))
    {
      daIntTyypp=reinterpret_cast< MEDCoupling::DataArrayInt * >(argp);
      sw=3;
      return;
    }
  status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayIntTuple,0|0);
  if(SWIG_IsOK(status))
    {
      daIntTuple=reinterpret_cast< MEDCoupling::DataArrayIntTuple * >(argp);
      sw=4;
      return ;
    }
  throw INTERP_KERNEL::Exception("5 types accepted : integer, tuple of integer, list of integer, DataArrayInt, DataArrayIntTuple");
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
static const int *convertIntStarLikePyObjToCppIntStar(PyObject *value, int& sw, int& sz, int& iTyypp, std::vector<int>& stdvecTyypp)
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
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayInt,0|0);
  if(SWIG_IsOK(status))
    {
      MEDCoupling::DataArrayInt *daIntTyypp=reinterpret_cast< MEDCoupling::DataArrayInt * >(argp);
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
  status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayIntTuple,0|0);
  if(SWIG_IsOK(status))
    {
      MEDCoupling::DataArrayIntTuple *daIntTuple=reinterpret_cast< MEDCoupling::DataArrayIntTuple * >(argp);
      sw=4; sz=daIntTuple->getNumberOfCompo();
      return daIntTuple->getConstPointer();
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
template<class T>
void considerPyObjAsATStarLikeObject(PyObject *value, int& sw, T& iTyypp, std::vector<T>& stdvecTyypp, typename MEDCoupling::Traits<T>::ArrayType *& daIntTyypp, swig_type_info *ti)
{
  sw=-1;
  if(PyFloat_Check(value))
    {
      iTyypp=(T)PyFloat_AS_DOUBLE(value);
      sw=1;
      return;
    }
  if(PyInt_Check(value))
    {
      iTyypp=(T)PyInt_AS_LONG(value);
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
            stdvecTyypp[i]=(T)PyFloat_AS_DOUBLE(o);
          else if(PyInt_Check(o))
            stdvecTyypp[i]=(T)PyInt_AS_LONG(o);
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
            stdvecTyypp[i]=(T)PyFloat_AS_DOUBLE(o);
          else if(PyInt_Check(o))
            stdvecTyypp[i]=(T)PyInt_AS_LONG(o);
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
  int status=SWIG_ConvertPtr(value,&argp,ti,0|0);
  if(!SWIG_IsOK(status))
    throw INTERP_KERNEL::Exception("5 types accepted : double float, integer, tuple of double float or int, list of double float or int, DataArrayDouble");
  daIntTyypp=reinterpret_cast< typename MEDCoupling::Traits<T>::ArrayType * >(argp);
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
static void convertDoubleStarLikePyObjToCpp(PyObject *value, int& sw, double& iTyypp, std::vector<double>& stdvecTyypp, MEDCoupling::DataArrayDoubleTuple *& daIntTyypp)
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
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayDoubleTuple,0|0);
  if(!SWIG_IsOK(status))
    throw INTERP_KERNEL::Exception("5 types accepted : double float, integer, tuple of double float or int, list of double float or int, DataArrayDoubleTuple");
  daIntTyypp=reinterpret_cast< MEDCoupling::DataArrayDoubleTuple * >(argp);
  sw=3;
}

template<class T>
void convertFPStarLikePyObjToCpp_2(PyObject *value, int& sw, T& val, typename MEDCoupling::Traits<T>::ArrayType *&d, typename MEDCoupling::Traits<T>::ArrayTuple *&e, std::vector<T>& f, swig_type_info *ti_da, swig_type_info *ti_tuple)
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
      val=(T)PyInt_AS_LONG(value);
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
            f[i]=(T)PyInt_AS_LONG(o);
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
            f[i]=(T)PyInt_AS_LONG(o);
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
  int status=SWIG_ConvertPtr(value,&argp,ti_da,0|0);
  if(SWIG_IsOK(status))
    {
      d=reinterpret_cast< typename MEDCoupling::Traits<T>::ArrayType * >(argp);
      sw=2;
      return ;
    }
  status=SWIG_ConvertPtr(value,&argp,ti_tuple,0|0);
  if(SWIG_IsOK(status))
    {
      e=reinterpret_cast< typename MEDCoupling::Traits<T>::ArrayTuple * >(argp);
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
static void convertDoubleStarLikePyObjToCpp_2(PyObject *value, int& sw, double& val, MEDCoupling::DataArrayDouble *&d, MEDCoupling::DataArrayDoubleTuple *&e, std::vector<double>& f)
{
  convertFPStarLikePyObjToCpp_2<double>(value,sw,val,d,e,f,SWIGTYPE_p_MEDCoupling__DataArrayDouble,SWIGTYPE_p_MEDCoupling__DataArrayDoubleTuple);
}

/*!
 * if value int -> cpp val sw=1
 * if value double -> cpp val sw=1
 * if value DataArrayDouble -> cpp DataArrayDouble sw=2
 * if value DataArrayDoubleTuple -> cpp DataArrayDoubleTuple sw=3
 * if value list[int,double] -> cpp std::vector<double> sw=4
 * if value tuple[int,double] -> cpp std::vector<double> sw=4
 */
static void convertFloatStarLikePyObjToCpp_2(PyObject *value, int& sw, float& val, MEDCoupling::DataArrayFloat *&d, MEDCoupling::DataArrayFloatTuple *&e, std::vector<float>& f)
{
  convertFPStarLikePyObjToCpp_2<float>(value,sw,val,d,e,f,SWIGTYPE_p_MEDCoupling__DataArrayFloat,SWIGTYPE_p_MEDCoupling__DataArrayFloatTuple);
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
static void convertIntStarOrSliceLikePyObjToCpp(PyObject *value, int nbelem, int& sw, int& iTyypp, std::vector<int>& stdvecTyypp, std::pair<int, std::pair<int,int> >& p, MEDCoupling::DataArrayInt *& daIntTyypp)
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
      GetIndicesOfSlice(value,nbelem,&strt,&stp,&step,"Slice in subscriptable object DataArray invalid !");
      p.first=strt;
      p.second.first=stp;
      p.second.second=step;
      sw=3;
      return ;
    }
  void *argp;
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayInt,0|0);
  if(SWIG_IsOK(status))
    {
      daIntTyypp=reinterpret_cast< MEDCoupling::DataArrayInt * >(argp);
      if(!daIntTyypp)
        {
          std::ostringstream oss; oss << msg << " Instance in null !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      sw=4;
      return ;
    }
  status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayIntTuple,0|0);
  if(SWIG_IsOK(status))
    {
      MEDCoupling::DataArrayIntTuple *tmp=reinterpret_cast< MEDCoupling::DataArrayIntTuple * >(argp);
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
 * Idem than convertIntStarOrSliceLikePyObjToCpp
 */
static void convertIntStarOrSliceLikePyObjToCppWithNegIntInterp(PyObject *value, int nbelem, int& sw, int& iTyypp, std::vector<int>& stdvecTyypp, std::pair<int, std::pair<int,int> >& p, MEDCoupling::DataArrayInt *& daIntTyypp)
{
  convertIntStarOrSliceLikePyObjToCpp(value,nbelem,sw,iTyypp,stdvecTyypp,p,daIntTyypp);
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
static void convertObjToPossibleCpp22(PyObject *value, int nbelem, int& sw, int& iTyypp, std::vector<int>& stdvecTyypp, std::pair<int, std::pair<int,int> >& p, MEDCoupling::DataArrayIntTuple *& daIntTyypp)
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
      GetIndicesOfSlice(value,nbelem,&strt,&stp,&step,"Slice in subscriptable object DataArray invalid !");
      p.first=strt;
      p.second.first=stp;
      p.second.second=step;
      sw=3;
      return ;
    }
  void *argp;
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayIntTuple,0|0);
  if(!SWIG_IsOK(status))
    throw INTERP_KERNEL::Exception("4 types accepted : integer, tuple of integer, list of integer, slice, DataArrayIntTuple");
  daIntTyypp=reinterpret_cast< MEDCoupling::DataArrayIntTuple * >(argp);
  sw=4;
}

/*!
 * if python string with size one -> cpp char sw=1
 * if python string with size different from one -> cpp string sw=2
 * if python tuple[string] or list[string] -> vector<string> sw=3
 * if python not null pointer of DataArrayChar -> cpp DataArrayChar sw=4
 * switch between (int,string,vector<string>,DataArrayChar)
 */
static void convertObjToPossibleCpp6(PyObject *value, int& sw, char& cTyp, std::string& sType, std::vector<std::string>& vsType, MEDCoupling::DataArrayChar *& dacType)
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
#if PY_VERSION_HEX >= 0x03000000
  if(PyUnicode_Check(value))
    {
      Py_ssize_t sz;
      const char *pt = PyUnicode_AsUTF8AndSize(value, &sz);
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
#endif
  if(PyTuple_Check(value))
    {
      int size=PyTuple_Size(value);
      vsType.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(value,i);
          try
            {
              vsType[i]=convertPyObjectToStr(o);
            }
          catch(INTERP_KERNEL::Exception& e)
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
          try
            {
              vsType[i]=convertPyObjectToStr(o);
            }
          catch(INTERP_KERNEL::Exception& e)
            {
              std::ostringstream oss; oss << "List as been detected but element #" << i << " is not a string ! only tuples of strings accepted !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      sw=3;
      return;
    }
  void *argp;
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayChar,0|0);
  if(SWIG_IsOK(status))
    {
      dacType=reinterpret_cast< MEDCoupling::DataArrayChar * >(argp);
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
                                     MEDCoupling::DataArrayInt *&dt, MEDCoupling::DataArrayInt *&dc)
{
  if(!PyTuple_Check(value))
    {
      convertIntStarOrSliceLikePyObjToCppWithNegIntInterp(value,nbTuple,sw,it,vt,pt,dt);
      return ;
    }
  else
    {
      int sz=PyTuple_Size(value);
      if(sz!=2)
        throw INTERP_KERNEL::Exception("Unexpected nb of slice element : 1 or 2 expected !\n1st is for tuple selection, 2nd for component selection !");
      PyObject *ob0=PyTuple_GetItem(value,0);
      int sw1,sw2;
      convertIntStarOrSliceLikePyObjToCppWithNegIntInterp(ob0,nbTuple,sw1,it,vt,pt,dt);
      PyObject *ob1=PyTuple_GetItem(value,1);
      convertIntStarOrSliceLikePyObjToCppWithNegIntInterp(ob1,nbCompo,sw2,ic,vc,pc,dc);
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
static const double *convertObjToPossibleCpp5_Safe(PyObject *value, int& sw, double& val, MEDCoupling::DataArrayDouble *&d, MEDCoupling::DataArrayDoubleTuple *&e, std::vector<double>& f,
                                                   const char *msg, int nbTuplesExpected, int nbCompExpected, bool throwIfNullPt)
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
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayDouble,0|0);
  if(SWIG_IsOK(status))
    {
      d=reinterpret_cast< MEDCoupling::DataArrayDouble * >(argp);
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
  status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayDoubleTuple,0|0);
  if(SWIG_IsOK(status))
    {
      e=reinterpret_cast< MEDCoupling::DataArrayDoubleTuple * >(argp);
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
static const double *convertObjToPossibleCpp5_Safe2(PyObject *value, int& sw, double& val, MEDCoupling::DataArrayDouble *&d, MEDCoupling::DataArrayDoubleTuple *&e, std::vector<double>& f,
                                                    const char *msg, int nbCompExpected, bool throwIfNullPt, int& nbTuples)
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
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayDouble,0|0);
  if(SWIG_IsOK(status))
    {
      d=reinterpret_cast< MEDCoupling::DataArrayDouble * >(argp);
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
  status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayDoubleTuple,0|0);
  if(SWIG_IsOK(status))
    {
      e=reinterpret_cast< MEDCoupling::DataArrayDoubleTuple * >(argp);
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
                                                          const char *msg, bool throwIfNullPt, int& nbTuples)
{
  MEDCoupling::DataArrayDouble *d=0;
  MEDCoupling::DataArrayDoubleTuple *e=0;
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
  int status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayDouble,0|0);
  if(SWIG_IsOK(status))
    {
      d=reinterpret_cast< MEDCoupling::DataArrayDouble * >(argp);
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
  status=SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayDoubleTuple,0|0);
  if(SWIG_IsOK(status))
    {
      e=reinterpret_cast< MEDCoupling::DataArrayDoubleTuple * >(argp);
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

static MEDCoupling::DataArray *CheckAndRetrieveDataArrayInstance(PyObject *obj, const char *msg)
{
  void *aBasePtrVS=0;
  int status=SWIG_ConvertPtr(obj,&aBasePtrVS,SWIGTYPE_p_MEDCoupling__DataArray,0|0);
  if(!SWIG_IsOK(status))
    {
      status=SWIG_ConvertPtr(obj,&aBasePtrVS,SWIGTYPE_p_MEDCoupling__DataArrayDouble,0|0);
      if(!SWIG_IsOK(status))
        {
          status=SWIG_ConvertPtr(obj,&aBasePtrVS,SWIGTYPE_p_MEDCoupling__DataArrayInt,0|0);
          if(!SWIG_IsOK(status))
            {
              status=SWIG_ConvertPtr(obj,&aBasePtrVS,SWIGTYPE_p_MEDCoupling__DataArrayAsciiChar,0|0);
              if(!SWIG_IsOK(status))
                {
                  status=SWIG_ConvertPtr(obj,&aBasePtrVS,SWIGTYPE_p_MEDCoupling__DataArrayByte,0|0);
                  std::ostringstream oss; oss << msg << " ! Accepted instances are DataArrayDouble, DataArrayInt, DataArrayAsciiChar, DataArrayByte !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
            }
        }
    }
  return reinterpret_cast< MEDCoupling::DataArray * >(aBasePtrVS);
}

static PyObject *NewMethWrapCallInitOnlyIfEmptyDictInInput(PyObject *cls, PyObject *args, const char *clsName)
{
  if(!PyTuple_Check(args))
    {
      std::ostringstream oss; oss << clsName << ".__new__ : the args in input is expected to be a tuple !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  PyObject *builtinsd(PyEval_GetBuiltins());//borrowed
  PyObject *obj(PyDict_GetItemString(builtinsd,"object"));//borrowed
  PyObject *selfMeth(PyObject_GetAttrString(obj,"__new__"));
  //
  PyObject *tmp0(PyTuple_New(1));
  PyTuple_SetItem(tmp0,0,cls); Py_XINCREF(cls);
  PyObject *instance(PyObject_CallObject(selfMeth,tmp0));
  Py_DECREF(tmp0);
  Py_DECREF(selfMeth);
  if(PyTuple_Size(args)==2 && PyDict_Check(PyTuple_GetItem(args,1)) && PyDict_Size(PyTuple_GetItem(args,1))==0 )
    {// NOT general case. only true if in unpickeling context ! call __init__. Because for all other cases, __init__ is called right after __new__ !
      PyObject *initMeth(PyObject_GetAttrString(instance,"__init__"));
      PyObject *tmp3(PyTuple_New(0));
      PyObject *tmp2(PyObject_CallObject(initMeth,tmp3));
      Py_XDECREF(tmp2);
      Py_DECREF(tmp3);
      Py_DECREF(initMeth);
    }
  return instance;
}

template<class T>
static PyObject *NewMethWrapCallInitOnlyIfDictWithSingleEltInInputGeneral(PyObject *cls, PyObject *args, const char *clsName)
{
  if(!PyTuple_Check(args))
    {
      std::ostringstream oss; oss << clsName << ".__new__ : the args in input is expected to be a tuple !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  PyObject *builtinsd(PyEval_GetBuiltins());//borrowed
  PyObject *obj(PyDict_GetItemString(builtinsd,"object"));//borrowed
  PyObject *selfMeth(PyObject_GetAttrString(obj,"__new__"));
  //
  PyObject *tmp0(PyTuple_New(1));
  PyTuple_SetItem(tmp0,0,cls); Py_XINCREF(cls);
  PyObject *instance(PyObject_CallObject(selfMeth,tmp0));
  Py_DECREF(tmp0);
  Py_DECREF(selfMeth);
  if(PyTuple_Size(args)==2 && PyDict_Check(PyTuple_GetItem(args,1)) && PyDict_Size(PyTuple_GetItem(args,1))==1 )
    {// NOT general case. only true if in unpickeling context ! call __init__. Because for all other cases, __init__ is called right after __new__ !
      PyObject *initMeth(PyObject_GetAttrString(instance,"__init__"));
      PyObject *zeNumpyRepr(0);
      {
        PyObject *tmp1(PyInt_FromLong(0));
       zeNumpyRepr=PyDict_GetItem(PyTuple_GetItem(args,1),tmp1);//borrowed
        Py_DECREF(tmp1);
      }
      if(!zeNumpyRepr)
        {
          std::ostringstream oss; oss << clsName << ".__new__ : the args in input is expected to be a tuple !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      T tt;
      {
        PyObject *tmp3(0);
        try
          {
            tmp3=tt(zeNumpyRepr);
          }
        catch(INTERP_KERNEL::Exception& e)
          {
            std::ostringstream oss; oss << clsName << ".__new__ : Invalid type in input " << " : " << e.what();
            throw INTERP_KERNEL::Exception(oss.str());
          }
        {
          PyObject *tmp2(PyObject_CallObject(initMeth,tmp3));
          Py_XDECREF(tmp2);
        }
        Py_DECREF(tmp3);
      }
      Py_DECREF(initMeth);
    }
  return instance;
}

struct SinglePyObjToBePutInATuple
{
  PyObject *operator()(PyObject *zeNumpyRepr)
  {
    PyObject *tmp3(PyTuple_New(1));
    PyTuple_SetItem(tmp3,0,zeNumpyRepr); Py_XINCREF(zeNumpyRepr);
    return tmp3;
  }
};

struct SinglePyObjExpectToBeAListOfSz2
{
  PyObject *operator()(PyObject *uniqueElt)
  {
    if(!PyTuple_Check(uniqueElt) || PyTuple_Size(uniqueElt)!=2)
      throw INTERP_KERNEL::Exception("Not a tuple of size 2 !");
    Py_XINCREF(uniqueElt);
    return uniqueElt;
  }
};

static PyObject *NewMethWrapCallInitOnlyIfDictWithSingleEltInInput(PyObject *cls, PyObject *args, const char *clsName)
{
  return NewMethWrapCallInitOnlyIfDictWithSingleEltInInputGeneral<SinglePyObjToBePutInATuple>(cls,args,clsName);
}

static PyObject *convertPartDefinition(MEDCoupling::PartDefinition *pd, int owner)
{
  PyObject *ret=0;
  if(!pd)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDCoupling::DataArrayPartDefinition *>(pd))
    ret=SWIG_NewPointerObj((void*)pd,SWIGTYPE_p_MEDCoupling__DataArrayPartDefinition,owner);
  else
    ret=SWIG_NewPointerObj((void*)pd,SWIGTYPE_p_MEDCoupling__SlicePartDefinition,owner);
  return ret;
}

template<class T>
static typename MEDCoupling::Traits<T>::ArrayType *DataArrayT_New(PyObject *elt0, PyObject *nbOfTuples, PyObject *elt2)
{
  const char *msgBase="MEDCoupling::DataArrayDouble::New : Available API are : \n-DataArrayDouble.New()\n-DataArrayDouble.New([1.,3.,4.])\n-DataArrayDouble.New([1.,3.,4.],3)\n-DataArrayDouble.New([1.,3.,4.,5.],2,2)\n-DataArrayDouble.New([1.,3.,4.,5.,7,8.],3,2)\n-DataArrayDouble.New([(1.,3.),(4.,5.),(7,8.)])\n-DataArrayDouble.New(5)\n-DataArrayDouble.New(5,2)";
  std::string msg(msgBase);
#ifdef WITH_NUMPY
  msg+="\n-DataArrayDouble.New(numpy array with dtype=float64)";
#endif
  msg+=" !";
  if(PyList_Check(elt0) || PyTuple_Check(elt0))
    {
      if(nbOfTuples)
        {
          if(PyInt_Check(nbOfTuples))
            {
              int nbOfTuples1=PyInt_AS_LONG(nbOfTuples);
              if(nbOfTuples1<0)
                throw INTERP_KERNEL::Exception("DataArrayDouble::New : should be a positive set of allocated memory !");
              if(elt2)
                {
                  if(PyInt_Check(elt2))
                    {//DataArrayDouble.New([1.,3.,4.,5.],2,2)
                      int nbOfCompo=PyInt_AS_LONG(elt2);
                      if(nbOfCompo<0)
                        throw INTERP_KERNEL::Exception("DataArrayDouble::New : should be a positive number of components !");
                      MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > ret(MEDCoupling::Traits<T>::ArrayType::New());
                      std::vector<double> tmp(fillArrayWithPyListDbl2(elt0,nbOfTuples1,nbOfCompo));
                      ret->alloc(nbOfTuples1,nbOfCompo); std::copy(tmp.begin(),tmp.end(),ret->getPointer());
                      return ret.retn();
                    }
                  else
                    throw INTERP_KERNEL::Exception(msg.c_str());
                }
              else
                {//DataArrayDouble.New([1.,3.,4.],3)
                  MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > ret(MEDCoupling::Traits<T>::ArrayType::New());
                  int tmpp1(-1);
                  std::vector<double> tmp(fillArrayWithPyListDbl2(elt0,nbOfTuples1,tmpp1));
                  ret->alloc(nbOfTuples1,tmpp1); std::copy(tmp.begin(),tmp.end(),ret->getPointer());
                  return ret.retn();
                }
            }
          else
            throw INTERP_KERNEL::Exception(msg.c_str());
        }
      else
        {// DataArrayDouble.New([1.,3.,4.])
          MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > ret(MEDCoupling::Traits<T>::ArrayType::New());
          int tmpp1(-1),tmpp2(-1);
          std::vector<double> tmp=fillArrayWithPyListDbl2(elt0,tmpp1,tmpp2);
          ret->alloc(tmpp1,tmpp2); std::copy(tmp.begin(),tmp.end(),ret->getPointer());
          return ret.retn();
        }
    }
  else if(PyInt_Check(elt0))
    {
      int nbOfTuples1(PyInt_AS_LONG(elt0));
      if(nbOfTuples1<0)
        throw INTERP_KERNEL::Exception("DataArrayDouble::New : should be a positive set of allocated memory !");
      if(nbOfTuples)
        {
          if(!elt2)
            {
              if(PyInt_Check(nbOfTuples))
                {//DataArrayDouble.New(5,2)
                  int nbOfCompo=PyInt_AS_LONG(nbOfTuples);
                  if(nbOfCompo<0)
                    throw INTERP_KERNEL::Exception("DataArrayDouble::New : should be a positive number of components !");
                  MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > ret(MEDCoupling::Traits<T>::ArrayType::New());
                  ret->alloc(nbOfTuples1,nbOfCompo);
                  return ret.retn();
                }
              else
                throw INTERP_KERNEL::Exception(msg.c_str());
            }
          else
            throw INTERP_KERNEL::Exception(msg.c_str());
        }
      else
        {//DataArrayDouble.New(5)
          MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > ret(MEDCoupling::Traits<T>::ArrayType::New());
          ret->alloc(nbOfTuples1,1);
          return ret.retn();
        }
    }
#ifdef WITH_NUMPY
  else if(PyArray_Check(elt0) && nbOfTuples==NULL && elt2==NULL)
    {//DataArrayDouble.New(numpyArray)
      return BuildNewInstance< typename MEDCoupling::Traits<T>::ArrayType , T >(elt0,NPYTraits<T>::NPYObjectType,NPYTraits<T>::NPYFunc,MEDCoupling::Traits<T>::NPYStr);
    }
#endif
  else
    throw INTERP_KERNEL::Exception(msg.c_str());
  throw INTERP_KERNEL::Exception(msg.c_str());//to make g++ happy
}

template<class T>
typename MEDCoupling::Traits<T>::ArrayType *DataArrayT__setitem__internal(typename MEDCoupling::Traits<T>::ArrayType *self, PyObject *obj, PyObject *value, swig_type_info *ti)
{
  self->checkAllocated();
  const char msg[]="Unexpected situation in DataArrayDouble::__setitem__ !";
  int nbOfTuples(self->getNumberOfTuples()),nbOfComponents(self->getNumberOfComponents());
  int sw1,sw2;
  T i1;
  std::vector<T> v1;
  typename MEDCoupling::Traits<T>::ArrayType *d1=0;
  considerPyObjAsATStarLikeObject<T>(value,sw1,i1,v1,d1,ti);
  int it1,ic1;
  std::vector<int> vt1,vc1;
  std::pair<int, std::pair<int,int> > pt1,pc1;
  MEDCoupling::DataArrayInt *dt1=0,*dc1=0;
  convertObjToPossibleCpp3(obj,nbOfTuples,nbOfComponents,sw2,it1,ic1,vt1,vc1,pt1,pc1,dt1,dc1);
  MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > tmp;
  switch(sw2)
    {
    case 1:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple1(i1,it1,it1+1,1,0,nbOfComponents,1);
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues1(tmp,it1,it1+1,1,0,nbOfComponents,1,false);
            return self;
          case 3:
            self->setPartOfValues1(d1,it1,it1+1,1,0,nbOfComponents,1);
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    case 2:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple3(i1,&vt1[0],&vt1[0]+vt1.size(),0,nbOfComponents,1);
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues3(tmp,&vt1[0],&vt1[0]+vt1.size(),0,nbOfComponents,1,false);
            return self;
          case 3:
            self->setPartOfValues3(d1,&vt1[0],&vt1[0]+vt1.size(),0,nbOfComponents,1);
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    case 3:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple1(i1,pt1.first,pt1.second.first,pt1.second.second,0,nbOfComponents,1);
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues1(tmp,pt1.first,pt1.second.first,pt1.second.second,0,nbOfComponents,1,false);
            return self;
          case 3:
            self->setPartOfValues1(d1,pt1.first,pt1.second.first,pt1.second.second,0,nbOfComponents,1);
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    case 4:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple3(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),0,nbOfComponents,1);
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues3(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),0,nbOfComponents,1,false);
            return self;
          case 3:
            self->setPartOfValues3(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),0,nbOfComponents,1);
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    case 5:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple1(i1,it1,it1+1,1,ic1,ic1+1,1);
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues1(tmp,it1,it1+1,1,ic1,ic1+1,1,false);
            return self;
          case 3:
            self->setPartOfValues1(d1,it1,it1+1,1,ic1,ic1+1,1);
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    case 6:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple3(i1,&vt1[0],&vt1[0]+vt1.size(),ic1,ic1+1,1);
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues3(tmp,&vt1[0],&vt1[0]+vt1.size(),ic1,ic1+1,1,false);
            return self;
          case 3:
            self->setPartOfValues3(d1,&vt1[0],&vt1[0]+vt1.size(),ic1,ic1+1,1);
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    case 7:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple1(i1,pt1.first,pt1.second.first,pt1.second.second,ic1,ic1+1,1);
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues1(tmp,pt1.first,pt1.second.first,pt1.second.second,ic1,ic1+1,1,false);
            return self;
          case 3:
            self->setPartOfValues1(d1,pt1.first,pt1.second.first,pt1.second.second,ic1,ic1+1,1);
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    case 8:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple3(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),ic1,ic1+1,1);
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues3(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),ic1,ic1+1,1,false);
            return self;
          case 3:
            self->setPartOfValues3(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),ic1,ic1+1,1);
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    case 9:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple2(i1,&it1,&it1+1,&vc1[0],&vc1[0]+vc1.size());
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues2(tmp,&it1,&it1+1,&vc1[0],&vc1[0]+vc1.size(),false);
            return self;
          case 3:
            self->setPartOfValues2(d1,&it1,&it1+1,&vc1[0],&vc1[0]+vc1.size());
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    case 10:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple2(i1,&vt1[0],&vt1[0]+vt1.size(),&vc1[0],&vc1[0]+vc1.size());
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues2(tmp,&vt1[0],&vt1[0]+vt1.size(),&vc1[0],&vc1[0]+vc1.size(),false);
            return self;
          case 3:
            self->setPartOfValues2(d1,&vt1[0],&vt1[0]+vt1.size(),&vc1[0],&vc1[0]+vc1.size());
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    case 11:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple4(i1,pt1.first,pt1.second.first,pt1.second.second,&vc1[0],&vc1[0]+vc1.size());
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues4(tmp,pt1.first,pt1.second.first,pt1.second.second,&vc1[0],&vc1[0]+vc1.size(),false);
            return self;
          case 3:
            self->setPartOfValues4(d1,pt1.first,pt1.second.first,pt1.second.second,&vc1[0],&vc1[0]+vc1.size());
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    case 12:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple2(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),&vc1[0],&vc1[0]+vc1.size());
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues2(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),&vc1[0],&vc1[0]+vc1.size(),false);
            return self;
          case 3:
            self->setPartOfValues2(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),&vc1[0],&vc1[0]+vc1.size());
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    case 13:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple1(i1,it1,it1+1,1,pc1.first,pc1.second.first,pc1.second.second);
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues1(tmp,it1,it1+1,1,pc1.first,pc1.second.first,pc1.second.second,false);
            return self;
          case 3:
            self->setPartOfValues1(d1,it1,it1+1,1,pc1.first,pc1.second.first,pc1.second.second);
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    case 14:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple3(i1,&vt1[0],&vt1[0]+vt1.size(),pc1.first,pc1.second.first,pc1.second.second);
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues3(tmp,&vt1[0],&vt1[0]+vt1.size(),pc1.first,pc1.second.first,pc1.second.second,false);
            return self;
          case 3:
            self->setPartOfValues3(d1,&vt1[0],&vt1[0]+vt1.size(),pc1.first,pc1.second.first,pc1.second.second);
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    case 15:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple1(i1,pt1.first,pt1.second.first,pt1.second.second,pc1.first,pc1.second.first,pc1.second.second);
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues1(tmp,pt1.first,pt1.second.first,pt1.second.second,pc1.first,pc1.second.first,pc1.second.second,false);
            return self;
          case 3:
            self->setPartOfValues1(d1,pt1.first,pt1.second.first,pt1.second.second,pc1.first,pc1.second.first,pc1.second.second);
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    case 16:
      {
        switch(sw1)
          {
          case 1:
            self->setPartOfValuesSimple3(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),pc1.first,pc1.second.first,pc1.second.second);
            return self;
          case 2:
            tmp=MEDCoupling::Traits<T>::ArrayType::New();
            tmp->useArray(&v1[0],false,MEDCoupling::CPP_DEALLOC,1,v1.size());
            self->setPartOfValues3(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),pc1.first,pc1.second.first,pc1.second.second,false);
            return self;
          case 3:
            self->setPartOfValues3(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),pc1.first,pc1.second.first,pc1.second.second);
            return self;
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        break;
      }
    default:
      throw INTERP_KERNEL::Exception(msg);
    }
  return self;
}

template<class T>
PyObject *DataArrayT__getitem__internal(const typename MEDCoupling::Traits<T>::ArrayType *self, PyObject *obj, swig_type_info *ti)
{
  const char msg[]="Unexpected situation in DataArrayDouble::__getitem__ !";
  const char msg2[]="DataArrayDouble::__getitem__ : Mismatch of slice values in 2nd parameter (components) !";
  self->checkAllocated();
  int nbOfTuples(self->getNumberOfTuples()),nbOfComponents(self->getNumberOfComponents());
  int it1,ic1;
  std::vector<int> vt1,vc1;
  std::pair<int, std::pair<int,int> > pt1,pc1;
  MEDCoupling::DataArrayInt *dt1=0,*dc1=0;
  int sw;
  convertObjToPossibleCpp3(obj,nbOfTuples,nbOfComponents,sw,it1,ic1,vt1,vc1,pt1,pc1,dt1,dc1);
  MEDCoupling::MCAuto<typename MEDCoupling::Traits<T>::ArrayType > ret;
  switch(sw)
    {
    case 1:
      if(nbOfComponents==1)
        return PyFloat_FromDouble((T)self->getIJSafe(it1,0));
      return SWIG_NewPointerObj(SWIG_as_voidptr(self->selectByTupleIdSafe(&it1,&it1+1)),ti, SWIG_POINTER_OWN | 0 );
    case 2:
      return SWIG_NewPointerObj(SWIG_as_voidptr(self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size())),ti, SWIG_POINTER_OWN | 0 );
    case 3:
      return SWIG_NewPointerObj(SWIG_as_voidptr(self->selectByTupleIdSafeSlice(pt1.first,pt1.second.first,pt1.second.second)),ti, SWIG_POINTER_OWN | 0 );
    case 4:
      return SWIG_NewPointerObj(SWIG_as_voidptr(self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems())),ti, SWIG_POINTER_OWN | 0 );
    case 5:
      return PyFloat_FromDouble((T)self->getIJSafe(it1,ic1));
    case 6:
      {
        ret=self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size());
        std::vector<int> v2(1,ic1);
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(v2)),ti, SWIG_POINTER_OWN | 0 );
      }
    case 7:
      {
        ret=self->selectByTupleIdSafeSlice(pt1.first,pt1.second.first,pt1.second.second);
        std::vector<int> v2(1,ic1);
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(v2)),ti, SWIG_POINTER_OWN | 0 );
      }
    case 8:
      {
        ret=self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems());
        std::vector<int> v2(1,ic1);
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(v2)),ti, SWIG_POINTER_OWN | 0 );
      }
    case 9:
      {
        ret=self->selectByTupleIdSafe(&it1,&it1+1);
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(vc1)),ti, SWIG_POINTER_OWN | 0 );
      }
    case 10:
      {
        ret=self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size());
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(vc1)),ti, SWIG_POINTER_OWN | 0 );
      }
    case 11:
      {
        ret=self->selectByTupleIdSafeSlice(pt1.first,pt1.second.first,pt1.second.second);
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(vc1)),ti, SWIG_POINTER_OWN | 0 );
      }
    case 12:
      {
        ret=self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems());
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(vc1)),ti, SWIG_POINTER_OWN | 0 );
      }
    case 13:
      {
        ret=self->selectByTupleIdSafe(&it1,&it1+1);
        int nbOfComp(MEDCoupling::DataArray::GetNumberOfItemGivenBESRelative(pc1.first,pc1.second.first,pc1.second.second,msg2));
        std::vector<int> v2(nbOfComp);
        for(int i=0;i<nbOfComp;i++)
          v2[i]=pc1.first+i*pc1.second.second;
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(v2)),ti, SWIG_POINTER_OWN | 0 );
      }
    case 14:
      {
        ret=self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size());
        int nbOfComp(MEDCoupling::DataArray::GetNumberOfItemGivenBESRelative(pc1.first,pc1.second.first,pc1.second.second,msg2));
        std::vector<int> v2(nbOfComp);
        for(int i=0;i<nbOfComp;i++)
          v2[i]=pc1.first+i*pc1.second.second;
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(v2)),ti, SWIG_POINTER_OWN | 0 );
      }
    case 15:
      {
        ret=self->selectByTupleIdSafeSlice(pt1.first,pt1.second.first,pt1.second.second);
        int nbOfComp(MEDCoupling::DataArray::GetNumberOfItemGivenBESRelative(pc1.first,pc1.second.first,pc1.second.second,msg2));
        std::vector<int> v2(nbOfComp);
        for(int i=0;i<nbOfComp;i++)
          v2[i]=pc1.first+i*pc1.second.second;
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(v2)),ti, SWIG_POINTER_OWN | 0 );
      }
    case 16:
      {
        ret=self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems());
        int nbOfComp(MEDCoupling::DataArray::GetNumberOfItemGivenBESRelative(pc1.first,pc1.second.first,pc1.second.second,msg2));
        std::vector<int> v2(nbOfComp);
        for(int i=0;i<nbOfComp;i++)
          v2[i]=pc1.first+i*pc1.second.second;
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(v2)),ti, SWIG_POINTER_OWN | 0 );
      }
    default:
      throw INTERP_KERNEL::Exception(msg);
    }
}

template<class T>
PyObject *DataArrayT_imul__internal(PyObject *trueSelf, PyObject *obj, typename MEDCoupling::Traits<T>::ArrayType *self, swig_type_info *ti_da, swig_type_info *ti_tuple)
{
  const char msg[]="Unexpected situation in __imul__ !";
  T val;
  typename MEDCoupling::Traits<T>::ArrayType *a;
  typename MEDCoupling::Traits<T>::ArrayTuple *aa;
  std::vector<T> bb;
  int sw;
  convertFPStarLikePyObjToCpp_2<T>(obj,sw,val,a,aa,bb,ti_da,ti_tuple);
  switch(sw)
    {
    case 1:
      {
        self->applyLin(val,0.);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    case 2:
      {
        self->multiplyEqual(a);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    case 3:
      {
        MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > aaa(aa->buildDA(1,self->getNumberOfComponents()));
        self->multiplyEqual(aaa);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    case 4:
      {
        MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > aaa(MEDCoupling::Traits<T>::ArrayType::New()); aaa->useArray(&bb[0],false,MEDCoupling::CPP_DEALLOC,1,(int)bb.size());
        self->multiplyEqual(aaa);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    default:
      throw INTERP_KERNEL::Exception(msg);
    }
}

template<class T>
PyObject *DataArrayT_idiv__internal(PyObject *trueSelf, PyObject *obj, typename MEDCoupling::Traits<T>::ArrayType *self, swig_type_info *ti_da, swig_type_info *ti_tuple)
{
  const char msg[]="Unexpected situation in __idiv__ !";
  T val;
  typename MEDCoupling::Traits<T>::ArrayType *a;
  typename MEDCoupling::Traits<T>::ArrayTuple *aa;
  std::vector<T> bb;
  int sw;
  convertFPStarLikePyObjToCpp_2<T>(obj,sw,val,a,aa,bb,ti_da,ti_tuple);
  switch(sw)
    {
    case 1:
      {
        if(val==0.)
          throw INTERP_KERNEL::Exception("DataArrayDouble::__div__ : trying to divide by zero !");
        self->applyLin(1./val,0.);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    case 2:
      {
        self->divideEqual(a);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    case 3:
      {
        MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > aaa(aa->buildDA(1,self->getNumberOfComponents()));
        self->divideEqual(aaa);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    case 4:
      {
        MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > aaa(MEDCoupling::Traits<T>::ArrayType::New()); aaa->useArray(&bb[0],false,MEDCoupling::CPP_DEALLOC,1,(int)bb.size());
        self->divideEqual(aaa);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    default:
      throw INTERP_KERNEL::Exception(msg);
    }
}

template<class T>
PyObject *DataArrayT_iadd__internal(PyObject *trueSelf, PyObject *obj, typename MEDCoupling::Traits<T>::ArrayType *self, swig_type_info *ti_da, swig_type_info *ti_tuple)
{
  const char msg[]="Unexpected situation in __iadd__ !";
  T val;
  typename MEDCoupling::Traits<T>::ArrayType *a;
  typename MEDCoupling::Traits<T>::ArrayTuple *aa;
  std::vector<T> bb;
  int sw;
  convertFPStarLikePyObjToCpp_2<T>(obj,sw,val,a,aa,bb,ti_da,ti_tuple);
  switch(sw)
    {
    case 1:
      {
        self->applyLin(1.,val);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    case 2:
      {
        self->addEqual(a);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    case 3:
      {
        MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > aaa(aa->buildDA(1,self->getNumberOfComponents()));
        self->addEqual(aaa);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    case 4:
      {
        MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > aaa(MEDCoupling::Traits<T>::ArrayType::New()); aaa->useArray(&bb[0],false,MEDCoupling::CPP_DEALLOC,1,(int)bb.size());
        self->addEqual(aaa);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    default:
      throw INTERP_KERNEL::Exception(msg);
    }
}

template<class T>
PyObject *DataArrayT_isub__internal(PyObject *trueSelf, PyObject *obj, typename MEDCoupling::Traits<T>::ArrayType *self, swig_type_info *ti_da, swig_type_info *ti_tuple)
{
  const char msg[]="Unexpected situation in __isub__ !";
  T val;
  typename MEDCoupling::Traits<T>::ArrayType *a;
  typename MEDCoupling::Traits<T>::ArrayTuple *aa;
  std::vector<T> bb;
  int sw;
  convertFPStarLikePyObjToCpp_2<T>(obj,sw,val,a,aa,bb,ti_da,ti_tuple);
  switch(sw)
    {
    case 1:
      {
        self->applyLin(1.,-val);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    case 2:
      {
        self->substractEqual(a);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    case 3:
      {
        MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > aaa(aa->buildDA(1,self->getNumberOfComponents()));
        self->substractEqual(aaa);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    case 4:
      {
        MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > aaa(MEDCoupling::Traits<T>::ArrayType::New()); aaa->useArray(&bb[0],false,MEDCoupling::CPP_DEALLOC,1,(int)bb.size());
        self->substractEqual(aaa);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
    default:
      throw INTERP_KERNEL::Exception(msg);
    }
}

template<class T>
struct SWIGTITraits
{ };

template<>
struct SWIGTITraits<double>
{ static swig_type_info *TI; static swig_type_info *TI_TUPLE; };

template<>
struct SWIGTITraits<float>
{ static swig_type_info *TI; static swig_type_info *TI_TUPLE; };

template<>
struct SWIGTITraits<int>
{ static swig_type_info *TI; static swig_type_info *TI_TUPLE; };

swig_type_info *SWIGTITraits<double>::TI=NULL;//unfortunately SWIGTYPE_p_MEDCoupling__DataArrayDouble is null when called here ! Postpone initialization at inlined initializeMe()
swig_type_info *SWIGTITraits<float>::TI=NULL;//unfortunately SWIGTYPE_p_MEDCoupling__DataArrayFloat is null when called here ! Postpone initialization at inlined initializeMe()
swig_type_info *SWIGTITraits<int>::TI=NULL;//unfortunately SWIGTYPE_p_MEDCoupling__DataArrayFloat is null when called here ! Postpone initialization at inlined initializeMe()
swig_type_info *SWIGTITraits<double>::TI_TUPLE=NULL;//unfortunately SWIGTYPE_p_MEDCoupling__DataArrayDouble is null when called here ! Postpone initialization at inlined initializeMe()
swig_type_info *SWIGTITraits<float>::TI_TUPLE=NULL;//unfortunately SWIGTYPE_p_MEDCoupling__DataArrayFloat is null when called here ! Postpone initialization at inlined initializeMe()
swig_type_info *SWIGTITraits<int>::TI_TUPLE=NULL;//unfortunately SWIGTYPE_p_MEDCoupling__DataArrayFloat is null when called here ! Postpone initialization at inlined initializeMe()

#ifdef WITH_NUMPY
PyTypeObject *NPYTraits<double>::NPYFunc=&PyCallBackDataArrayDouble_RefType;

PyTypeObject *NPYTraits<float>::NPYFunc=&PyCallBackDataArrayFloat_RefType;
#endif

template<class T>
typename MEDCoupling::Traits<T>::ArrayType *DataArrayT__setitem__(typename MEDCoupling::Traits<T>::ArrayType *self, PyObject *obj, PyObject *value)
{
  return DataArrayT__setitem__internal<T>(self,obj,value,SWIGTITraits<T>::TI);
}

template<class T>
PyObject *DataArrayT__getitem(const typename MEDCoupling::Traits<T>::ArrayType *self, PyObject *obj)
{
  return DataArrayT__getitem__internal<T>(self,obj,SWIGTITraits<T>::TI);
}

template<class T>
PyObject *DataArrayT_imul(PyObject *trueSelf, PyObject *obj, typename MEDCoupling::Traits<T>::ArrayType *self)
{
  return DataArrayT_imul__internal<T>(trueSelf,obj,self,SWIGTITraits<T>::TI,SWIGTITraits<T>::TI_TUPLE);
}

template<class T>
PyObject *DataArrayT_idiv(PyObject *trueSelf, PyObject *obj, typename MEDCoupling::Traits<T>::ArrayType *self)
{
  return DataArrayT_idiv__internal<T>(trueSelf,obj,self,SWIGTITraits<T>::TI,SWIGTITraits<T>::TI_TUPLE);
}

template<class T>
PyObject *DataArrayT_iadd(PyObject *trueSelf, PyObject *obj, typename MEDCoupling::Traits<T>::ArrayType *self)
{
  return DataArrayT_iadd__internal<T>(trueSelf,obj,self,SWIGTITraits<T>::TI,SWIGTITraits<T>::TI_TUPLE);
}

template<class T>
PyObject *DataArrayT_isub(PyObject *trueSelf, PyObject *obj, typename MEDCoupling::Traits<T>::ArrayType *self)
{
  return DataArrayT_isub__internal<T>(trueSelf,obj,self,SWIGTITraits<T>::TI,SWIGTITraits<T>::TI_TUPLE);
}

template<class T>
typename MEDCoupling::Traits<T>::ArrayType *DataArrayFPT_rmul(typename MEDCoupling::Traits<T>::ArrayType *self, PyObject *obj)
{
  const char msg[]="Unexpected situation in __rmul__ !";
  T val;
  typename MEDCoupling::Traits<T>::ArrayType *a;
  typename MEDCoupling::Traits<T>::ArrayTuple *aa;
  std::vector<T> bb;
  int sw;
  convertFPStarLikePyObjToCpp_2<T>(obj,sw,val,a,aa,bb,SWIGTITraits<T>::TI,SWIGTITraits<T>::TI_TUPLE);
  switch(sw)
    {
    case 1:
      {
        typename MEDCoupling::MCAuto<typename MEDCoupling::Traits<T>::ArrayType> ret(self->deepCopy());
        ret->applyLin(val,0.);
        return ret.retn();
      }
    case 3:
      {
        typename MEDCoupling::MCAuto<typename MEDCoupling::Traits<T>::ArrayType> aaa(aa->buildDA(1,self->getNumberOfComponents()));
        return MEDCoupling::Traits<T>::ArrayType::Multiply(self,aaa);
      }
    case 4:
      {
        typename MEDCoupling::MCAuto<typename MEDCoupling::Traits<T>::ArrayType> aaa(MEDCoupling::Traits<T>::ArrayType::New()); aaa->useArray(&bb[0],false,MEDCoupling::CPP_DEALLOC,1,(int)bb.size());
        return MEDCoupling::Traits<T>::ArrayType::Multiply(self,aaa);
      }
    default:
      throw INTERP_KERNEL::Exception(msg);
    }
}

#endif
