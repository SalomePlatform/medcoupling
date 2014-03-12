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

#include <vector>

static PyObject *convertMEDFileMesh(ParaMEDMEM::MEDFileMesh* mesh, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!mesh)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<ParaMEDMEM::MEDFileUMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDFileUMesh,owner);
  if(dynamic_cast<ParaMEDMEM::MEDFileCMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDFileCMesh,owner);
  if(dynamic_cast<ParaMEDMEM::MEDFileCurveLinearMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDFileCurveLinearMesh,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of MEDFileMesh on downcast !");
  return ret;
}

static PyObject *convertMEDFileParameter1TS(ParaMEDMEM::MEDFileParameter1TS* p1ts, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!p1ts)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDFileParameterDouble1TS *>(p1ts))
    ret=SWIG_NewPointerObj((void*)p1ts,SWIGTYPE_p_ParaMEDMEM__MEDFileParameterDouble1TS,owner);
  if(dynamic_cast<MEDFileParameterDouble1TSWTI *>(p1ts))
    ret=SWIG_NewPointerObj((void*)p1ts,SWIGTYPE_p_ParaMEDMEM__MEDFileParameterDouble1TSWTI,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of MEDFileParameter1TS on downcast !");
  return ret;
}

static PyObject *convertMEDFileField1TS(ParaMEDMEM::MEDFileAnyTypeField1TS *p, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!p)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDFileField1TS *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_ParaMEDMEM__MEDFileField1TS,owner);
  if(dynamic_cast<MEDFileIntField1TS *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_ParaMEDMEM__MEDFileIntField1TS,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of MEDFileAnyTypeField1TS on downcast !");
  return ret;
}

static PyObject *convertMEDFileFieldMultiTS(ParaMEDMEM::MEDFileAnyTypeFieldMultiTS *p, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!p)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDFileFieldMultiTS *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_ParaMEDMEM__MEDFileFieldMultiTS,owner);
  if(dynamic_cast<MEDFileIntFieldMultiTS *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_ParaMEDMEM__MEDFileIntFieldMultiTS,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of MEDFileAnyTypeFieldMultiTS on downcast !");
  return ret;
}

static PyObject *convertMEDMeshMultiLev(ParaMEDMEM::MEDMeshMultiLev *p, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!p)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDUMeshMultiLev *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_ParaMEDMEM__MEDUMeshMultiLev,owner);
  if(dynamic_cast<MEDCMeshMultiLev *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_ParaMEDMEM__MEDCMeshMultiLev,owner);
  if(dynamic_cast<MEDCurveLinearMeshMultiLev *>(p))
    ret=SWIG_NewPointerObj((void*)p,SWIGTYPE_p_ParaMEDMEM__MEDCurveLinearMeshMultiLev,owner);
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

static void converPyListToVecString(PyObject *pyLi, std::vector<std::string>& v)
{
  if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      v.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(!PyString_Check(o))
            throw INTERP_KERNEL::Exception("In list passed in argument some elements are NOT strings ! Expected a list containing only strings !");
          const char *st=PyString_AsString(o);
          v[i]=std::string(st);
        }
    }
  else if(PyTuple_Check(pyLi))
    {
      int size=PyTuple_Size(pyLi);
      v.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyTuple_GetItem(pyLi,i);
          if(!PyString_Check(o))
            throw INTERP_KERNEL::Exception("In tuple passed in argument some elements are NOT strings ! Expected a tuple containing only strings !");
          const char *st=PyString_AsString(o);
          v[i]=std::string(st);
        }
    }
  else if(PyString_Check(pyLi))
    {
      v.resize(1);
      v[0]=std::string((const char *)PyString_AsString(pyLi));
    }
  else
    {
      throw INTERP_KERNEL::Exception("Unrecognized python argument : expected a list of string or tuple of string or string !");
    }
}

static PyObject *convertFieldDoubleVecToPy(const std::vector<ParaMEDMEM::MEDCouplingFieldDouble *>& li)
{
  int sz=li.size();
  PyObject *ret=PyList_New(sz);
  for(int i=0;i<sz;i++)
    {
      PyObject *o=SWIG_NewPointerObj((void*)li[i],SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,SWIG_POINTER_OWN | 0);
      PyList_SetItem(ret,i,o);
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
              if(PyString_Check(o0))
                p.first=std::string(PyString_AsString(o0));
              else
                throw INTERP_KERNEL::Exception(msg);
              PyObject *o1=PyTuple_GetItem(o,1);
              if(PyString_Check(o1))
                p.second=std::string(PyString_AsString(o1));
              else
                throw INTERP_KERNEL::Exception(msg);
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
                      if(PyString_Check(o0j))
                        {
                          p.first[j]=std::string(PyString_AsString(o0j));
                        }
                      else
                        throw INTERP_KERNEL::Exception(msg);
                    }
                }
              else
                throw INTERP_KERNEL::Exception(msg);
              PyObject *o1=PyTuple_GetItem(o,1);
              if(PyString_Check(o1))
                p.second=std::string(PyString_AsString(o1));
              else
                throw INTERP_KERNEL::Exception(msg);
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
  if(PyInt_Check(obj))
    {
      return InterpreteNegativeInt((int)PyInt_AS_LONG(obj),self->getNumberOfFields());
    }
  else if(PyString_Check(obj))
    {
      return self->getPosFromFieldName(PyString_AsString(obj));
    }
  else
    throw INTERP_KERNEL::Exception("MEDFileFields::__getitem__ : only integer or string with fieldname supported !");
}
