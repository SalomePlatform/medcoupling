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
// Author : Anthony Geay (EDF R&D)

#include <vector>

static PyObject *convertMEDFileMesh(MEDCoupling::MEDFileMesh* mesh, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!mesh)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDCoupling::MEDFileUMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_MEDCoupling__MEDFileUMesh,owner);
  if(dynamic_cast<MEDCoupling::MEDFileCMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_MEDCoupling__MEDFileCMesh,owner);
  if(dynamic_cast<MEDCoupling::MEDFileCurveLinearMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_MEDCoupling__MEDFileCurveLinearMesh,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of MEDFileMesh on downcast !");
  return ret;
}

static PyObject *convertMEDFileParameter1TS(MEDCoupling::MEDFileParameter1TS* p1ts, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!p1ts)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDFileParameterDouble1TS *>(p1ts))
    ret=SWIG_NewPointerObj((void*)p1ts,SWIGTYPE_p_MEDCoupling__MEDFileParameterDouble1TS,owner);
  if(dynamic_cast<MEDFileParameterDouble1TSWTI *>(p1ts))
    ret=SWIG_NewPointerObj((void*)p1ts,SWIGTYPE_p_MEDCoupling__MEDFileParameterDouble1TSWTI,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of MEDFileParameter1TS on downcast !");
  return ret;
}

static PyObject *convertMEDFileField1TS(MEDCoupling::MEDFileAnyTypeField1TS *p, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!p)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDFileField1TS *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_MEDCoupling__MEDFileField1TS,owner);
  if(dynamic_cast<MEDFileIntField1TS *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_MEDCoupling__MEDFileIntField1TS,owner);
  if(dynamic_cast<MEDFileFloatField1TS *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_MEDCoupling__MEDFileFloatField1TS,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of MEDFileAnyTypeField1TS on downcast !");
  return ret;
}

static PyObject *convertMEDFileFieldMultiTS(MEDCoupling::MEDFileAnyTypeFieldMultiTS *p, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!p)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDFileFieldMultiTS *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_MEDCoupling__MEDFileFieldMultiTS,owner);
  if(dynamic_cast<MEDFileIntFieldMultiTS *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_MEDCoupling__MEDFileIntFieldMultiTS,owner);
  if(dynamic_cast<MEDFileFloatFieldMultiTS *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_MEDCoupling__MEDFileFloatFieldMultiTS,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of MEDFileAnyTypeFieldMultiTS on downcast !");
  return ret;
}

static PyObject *convertMEDMeshMultiLev(MEDCoupling::MEDMeshMultiLev *p, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!p)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDUMeshMultiLev *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_MEDCoupling__MEDUMeshMultiLev,owner);
  if(dynamic_cast<MEDCMeshMultiLev *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_MEDCoupling__MEDCMeshMultiLev,owner);
  if(dynamic_cast<MEDCurveLinearMeshMultiLev *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_MEDCoupling__MEDCurveLinearMeshMultiLev,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of MEDMeshMultiLev on downcast !");
  return ret;
}

static std::vector<std::pair<int,int> > convertTimePairIdsFromPy(PyObject *pyLi) throw(INTERP_KERNEL::Exception)
{
  std::vector<std::pair<int,int> > ret;
  if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      ret.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(PyTuple_Check(o))
            {
              std::pair<int,int> p;
              int size2=PyTuple_Size(o);
              if(size2!=2)
                throw INTERP_KERNEL::Exception("tuples in list must be of size 2 (dt,it) !");
              PyObject *o0=PyTuple_GetItem(o,0);
              if(PyInt_Check(o0))
                p.first=(int)PyInt_AS_LONG(o0);
              else
                throw INTERP_KERNEL::Exception("First elem of tuples in list must be integer : dt !");
              PyObject *o1=PyTuple_GetItem(o,1);
              if(PyInt_Check(o1))
                p.second=(int)PyInt_AS_LONG(o1);
              else
                throw INTERP_KERNEL::Exception("Second elem of tuples in list must be integer : dt !");
              ret[i]=p;
            }
          else
            throw INTERP_KERNEL::Exception("list must contain tuples only");
        }
    }
  else
    throw INTERP_KERNEL::Exception("convertTimePairIdsFromPy : not a list");
  return ret;
}

static std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > convertVecPairIntToVecPairTOFCT(const std::vector<std::pair<int,int> >& tmp)
{
  std::size_t sz(tmp.size());
  std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > entitiesCpp(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      entitiesCpp[i].first=(TypeOfField)tmp[i].first;
      entitiesCpp[i].second=(INTERP_KERNEL::NormalizedCellType)tmp[i].second;
    }
  return entitiesCpp;
}

static void converPyListToVecString(PyObject *pyLi, std::vector<std::string>& v)
{
  static const char msg0[]="In list passed in argument some elements are NOT strings ! Expected a list containing only strings !";
  static const char msg1[]="In tuple passed in argument some elements are NOT strings ! Expected a list containing only strings !";
  static const char msg2[]="Unrecognized python argument : expected a list of string or tuple of string or string !";
  if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      v.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          v[i]=convertPyObjectToStr(o,msg0);
        }
      return ;
    }
  else if(PyTuple_Check(pyLi))
    {
      int size=PyTuple_Size(pyLi);
      v.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          v[i]=convertPyObjectToStr(o,msg1);
        }
      return ;
    }
  v.resize(1);
  v[0]=convertPyObjectToStr(pyLi,msg2);
}

static PyObject *convertFieldDoubleVecToPy(const std::vector<MEDCoupling::MEDCouplingFieldDouble *>& li)
{
  int sz=li.size();
  PyObject *ret=PyList_New(sz);
  for(int i=0;i<sz;i++)
    {
      PyObject *o=SWIG_NewPointerObj((void*)li[i],SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,SWIG_POINTER_OWN | 0);
      PyList_SetItem(ret,i,o);
    }
  return ret;
}

PyObject *convertVecPairIntToPy(const std::vector< std::pair<int,int> >& vec)
{
  PyObject *ret(PyList_New(vec.size()));
  int rk=0;
  for(std::vector< std::pair<int,int> >::const_iterator iter=vec.begin();iter!=vec.end();iter++,rk++)
    {
      PyObject *elt=PyTuple_New(2);
      PyTuple_SetItem(elt,0,SWIG_From_int((*iter).first));
      PyTuple_SetItem(elt,1,SWIG_From_int((*iter).second));
      PyList_SetItem(ret,rk,elt);
    }
  return ret;
}

PyObject *convertVecPairVecStToPy(const std::vector< std::pair<std::vector<std::string>, std::string > >& vec)
{
  int sz=(int)vec.size();
  PyObject *ret=PyList_New(sz);
  for(int i=0;i<sz;i++)
    {
      PyObject *t=PyTuple_New(2);
      int sz2=(int)vec[i].first.size();
      PyObject *ll=PyList_New(sz2);
      for(int j=0;j<sz2;j++)
        PyList_SetItem(ll,j,PyString_FromString(vec[i].first[j].c_str()));
      PyTuple_SetItem(t,0,ll);
      PyTuple_SetItem(t,1,PyString_FromString(vec[i].second.c_str()));
      PyList_SetItem(ret,i,t);
    }
  return ret;
}

PyObject *convertVectPairStToPy(const std::vector< std::pair<std::string, std::string > >& vec)
{
  int sz=(int)vec.size();
  PyObject *ret=PyList_New(sz);
  for(int i=0;i<sz;i++)
    {
      PyObject *t=PyTuple_New(2);
      PyTuple_SetItem(t,0,PyString_FromString(vec[i].first.c_str()));
      PyTuple_SetItem(t,1,PyString_FromString(vec[i].second.c_str()));
      PyList_SetItem(ret,i,t);
    }
  return ret;
}

std::vector< std::pair<std::string, std::string > > convertVecPairStStFromPy(PyObject *pyLi)
{
  std::vector< std::pair<std::string, std::string > > ret;
  const char *msg="convertVecPairStStFromPy : Expecting PyList of Tuples of size 2 ! The first elt in tuple is one string and the 2nd one a string !";
  if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      ret.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(PyTuple_Check(o))
            {
              std::pair<std::string, std::string> p;
              int size2=PyTuple_Size(o);
              if(size2!=2)
                throw INTERP_KERNEL::Exception(msg);
              PyObject *o0=PyTuple_GetItem(o,0);
              p.first=convertPyObjectToStr(o0,msg);
              PyObject *o1=PyTuple_GetItem(o,1);
              p.second=convertPyObjectToStr(o1,msg);
              ret[i]=p;
            }
          else
            throw INTERP_KERNEL::Exception(msg);
        }
      return ret;
    }
  throw INTERP_KERNEL::Exception(msg);
}

std::vector< std::pair<std::vector<std::string>, std::string > > convertVecPairVecStFromPy(PyObject *pyLi)
{
  std::vector< std::pair<std::vector<std::string>, std::string > > ret;
  const char *msg="convertVecPairVecStFromPy : Expecting PyList of Tuples of size 2 ! The first elt in tuple is a list of strings and the 2nd one a string !";
  if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      ret.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(PyTuple_Check(o))
            {
              std::pair<std::vector<std::string>, std::string> p;
              int size2=PyTuple_Size(o);
              if(size2!=2)
                throw INTERP_KERNEL::Exception(msg);
              PyObject *o0=PyTuple_GetItem(o,0);
              if(PyList_Check(o0))
                {
                  int size3=PyList_Size(o0);
                  p.first.resize(size3);
                  for(int j=0;j<size3;j++)
                    {
                      PyObject *o0j=PyList_GetItem(o0,j);
                      p.first[j]=convertPyObjectToStr(o0j,msg);
                    }
                }
              else
                throw INTERP_KERNEL::Exception(msg);
              PyObject *o1=PyTuple_GetItem(o,1);
              p.second=convertPyObjectToStr(o1,msg);
              ret[i]=p;
            }
          else
            throw INTERP_KERNEL::Exception(msg);
        }
      return ret;
    }
  throw INTERP_KERNEL::Exception(msg);
}

/*!
 * Called by MEDFileAnyTypeFieldMultiTS::__getitem__ when \a elt0 is neither a list nor a slice.
 * In this case a MEDFileAnyTypeField1TS object is returned.
 */
int MEDFileAnyTypeFieldMultiTSgetitemSingleTS__(const MEDFileAnyTypeFieldMultiTS *self, PyObject *elt0) throw(INTERP_KERNEL::Exception)
{
  if(elt0 && PyInt_Check(elt0))
    {//fmts[3]
      return InterpreteNegativeInt(PyInt_AS_LONG(elt0),self->getNumberOfTS());
    }
  else if(elt0 && PyTuple_Check(elt0))
    {
      if(PyTuple_Size(elt0)==2)
        {
          PyObject *o0=PyTuple_GetItem(elt0,0);
          PyObject *o1=PyTuple_GetItem(elt0,1);
          if(PyInt_Check(o0) && PyInt_Check(o1))
            {//fmts(1,-1)
              int iter=PyInt_AS_LONG(o0);
              int order=PyInt_AS_LONG(o1);
              return self->getPosOfTimeStep(iter,order);
            }
          else
            throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::__getitem__ : invalid input param ! input is a tuple of size 2 but two integers are expected in this tuple to request a time steps !");
        }
      else
        throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::__getitem__ : invalid input param ! input is a tuple of size != 2 ! two integers are expected in this tuple to request a time steps !");
    }
  else if(elt0 && PyFloat_Check(elt0))
    {
      double val=PyFloat_AS_DOUBLE(elt0);
      return self->getPosGivenTime(val);
    }
  else
    throw INTERP_KERNEL::Exception("MEDFileAnyTypeFieldMultiTS::__getitem__ : invalid input params ! expected fmts[int], fmts[int,int], or fmts[double] to request one time step ! To request a series of time steps invoke fmts[slice], fmts[list of int], fmts[list of double], or fmts[list of int,int] !");
}

/*!
 * Called by MEDFileAnyTypeFieldMultiTS::__getitem__ when \a obj is neither a list nor a slice.
 * In this case a MEDFileAnyTypeField1TS object is returned.
 */
int MEDFileFieldsgetitemSingleTS__(const MEDFileFields *self, PyObject *obj) throw(INTERP_KERNEL::Exception)
{
  static const char msg[]="MEDFileFields::__getitem__ : only integer or string with fieldname supported !";
  if(PyInt_Check(obj))
    {
      return InterpreteNegativeInt((int)PyInt_AS_LONG(obj),self->getNumberOfFields());
    }
  return self->getPosFromFieldName(convertPyObjectToStr(obj,msg));
}

void convertToMapIntDataArrayInt(PyObject *pyMap, std::map<int, MCAuto<DataArrayInt> >& cppMap)
{
  if(!PyDict_Check(pyMap))
    throw INTERP_KERNEL::Exception("convertToMapIntDataArrayInt : input is not a python map !");
  PyObject *key, *value;
  Py_ssize_t pos(0);
  cppMap.clear();
  while (PyDict_Next(pyMap,&pos,&key,&value))
    {
      if(!PyInt_Check(key))
        throw INTERP_KERNEL::Exception("convertToMapIntDataArrayInt : keys in map must be PyInt !");
      long k(PyInt_AS_LONG(key));
      void *argp(0);
      int status(SWIG_ConvertPtr(value,&argp,SWIGTYPE_p_MEDCoupling__DataArrayInt,0|0));
      if(!SWIG_IsOK(status))
        {
          std::ostringstream oss; oss << "convertToMapIntDataArrayInt : values in map must be DataArrayInt !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      DataArrayInt *arg(reinterpret_cast<DataArrayInt*>(argp));
      MCAuto<DataArrayInt> arg2(arg);
      if(arg)
        arg->incrRef();
      cppMap[k]=arg2;
    }
}

template<class T>
PyObject *MEDFileField1TS_getFieldWithProfile(const typename MLFieldTraits<T>::F1TSType *self, TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh) 
{
  DataArrayInt *ret1(NULL);
  typename MEDCoupling::Traits<T>::ArrayType *ret0(self->getFieldWithProfile(type,meshDimRelToMax,mesh,ret1));
  PyObject *ret(PyTuple_New(2));
  PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTITraits<T>::TI, SWIG_POINTER_OWN | 0 ));
  PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
  return ret;
}

