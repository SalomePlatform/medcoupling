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

#include <vector>

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

static std::vector<ParaMEDMEM::MEDCouplingUMesh *> convertFieldDoubleVecFromPy(PyObject *pyLi)
{
  std::vector<ParaMEDMEM::MEDCouplingUMesh *> ret;
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
