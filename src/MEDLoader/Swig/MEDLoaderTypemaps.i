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

#include <vector>

static PyObject* convertMEDFileMesh(ParaMEDMEM::MEDFileMesh* mesh, int owner)
{
  PyObject *ret=0;
  if(dynamic_cast<ParaMEDMEM::MEDFileUMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDFileUMesh,owner);
  if(dynamic_cast<ParaMEDMEM::MEDFileCMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDFileCMesh,owner);
  if(!ret)
    {
      PyErr_SetString(PyExc_TypeError,"Not recognized type of MEDFileMesh on downcast !");
      PyErr_Print();
    }
  return ret;
}

static std::vector<std::pair<int,int> > convertTimePairIdsFromPy(PyObject *pyLi)
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
                {
                  const char msg[]="tuples in list must be of size 2 (dt,it) !";
                  ret.clear();
                  PyErr_SetString(PyExc_TypeError,msg);
                  PyErr_Print();
                  throw INTERP_KERNEL::Exception(msg);
                }
              PyObject *o0=PyTuple_GetItem(o,0);
              if(PyInt_Check(o0))
                p.first=(int)PyInt_AS_LONG(o0);
              else
                {
                  const char msg[]="First elem of tuples in list must be integer : dt !";
                  ret.clear();
                  PyErr_SetString(PyExc_TypeError,msg);
                  PyErr_Print();
                  throw INTERP_KERNEL::Exception(msg);
                }
              PyObject *o1=PyTuple_GetItem(o,1);
              if(PyInt_Check(o1))
                p.second=(int)PyInt_AS_LONG(o1);
              else
                {
                  const char msg[]="Second elem of tuples in list must be integer : dt !";
                  ret.clear();
                  PyErr_SetString(PyExc_TypeError,msg);
                  PyErr_Print();
                  throw INTERP_KERNEL::Exception(msg);
                }
              ret[i]=p;
            }
          else
            {
              const char msg[]="list must contain tuples only";
              ret.clear();
              PyErr_SetString(PyExc_TypeError,msg);
              PyErr_Print();
              throw INTERP_KERNEL::Exception(msg);
              return ret;
            }
        }
    }
  else
    {
      ret.clear();
      const char msg[]="convertTimePairIdsFromPy : not a list";
      PyErr_SetString(PyExc_TypeError,msg);
      PyErr_Print();
      throw INTERP_KERNEL::Exception(msg);
    }
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

static std::vector<const ParaMEDMEM::MEDCouplingUMesh *> convertUMeshVecFromPy(PyObject *pyLi)
{
  std::vector<const ParaMEDMEM::MEDCouplingUMesh *> ret;
  if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      ret.resize(size);
      for(int i=0;i<size;i++)
        {
          PyObject *obj=PyList_GetItem(pyLi,i);
          void *argp;
          int status=SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,0|0);
          if(!SWIG_IsOK(status))
            {
              const char msg[]="list must contain only MEDCouplingUMesh";
              PyErr_SetString(PyExc_TypeError,msg);
              PyErr_Print();
              throw INTERP_KERNEL::Exception(msg);
            }
          ParaMEDMEM::MEDCouplingUMesh *arg=reinterpret_cast< ParaMEDMEM::MEDCouplingUMesh * >(argp);
          ret[i]=arg;
        }
    }
  else
    {
      ret.clear();
      const char msg[]="convertFieldDoubleVecFromPy : not a list";
      PyErr_SetString(PyExc_TypeError,msg);
      PyErr_Print();
      throw INTERP_KERNEL::Exception(msg);
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

std::vector< std::pair<std::vector<std::string>, std::string > > convertVecPairVecStFromPy(PyObject *pyLi)
{
  std::vector< std::pair<std::vector<std::string>, std::string > > ret;
  const char *msg="convertVecPairVecStFromPy : Expecting PyList of Tuples of size 2 ! The first elt in tupe is a list of strings and the 2nd one a string !";
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
                  int size2=PyList_Size(o0);
                  p.first.resize(size2);
                  for(int j=0;j<size2;j++)
                    {
                      PyObject *o0j=PyTuple_GetItem(o0,j);
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
