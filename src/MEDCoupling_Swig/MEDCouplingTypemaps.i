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

static PyObject* convertMesh(ParaMEDMEM::MEDCouplingMesh* mesh, int owner)
{
  PyObject *ret;
  if(dynamic_cast<ParaMEDMEM::MEDCouplingUMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,owner);
  if(dynamic_cast<ParaMEDMEM::MEDCouplingExtrudedMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDCouplingExtrudedMesh,owner);
  return ret;
}

static PyObject *convertIntArrToPyList(const int *ptr, int size)
{
  PyObject *ret=PyList_New(size);
  for(int i=0;i<size;i++)
    PyList_SetItem(ret,i,PyInt_FromLong(ptr[i]));
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
              PyErr_SetString(PyExc_TypeError,"list must contain integers only");
              return NULL;
            }
        }
      return tmp;
    }
  else
    {
      PyErr_SetString(PyExc_TypeError,"convertPyToNewIntArr : not a list");
      return 0;
    }
}

static PyObject *convertDblArrToPyList(const double *ptr, int size)
{
  PyObject *ret=PyList_New(size);
  for(int i=0;i<size;i++)
    PyList_SetItem(ret,i,PyFloat_FromDouble(ptr[i]));
  return ret;
}

static double *convertPyToNewDblArr2(PyObject *pyLi)
{
  if(PyList_Check(pyLi))
    {
      int size=PyList_Size(pyLi);
      double *tmp=new double[size];
      for(int i=0;i<size;i++)
        {
          PyObject *o=PyList_GetItem(pyLi,i);
          if(PyFloat_Check(o))
            {
              double val=PyFloat_AS_DOUBLE(o);
              tmp[i]=val;
            }
          else
            {
              PyErr_SetString(PyExc_TypeError,"list must contain floats only");
              PyErr_Print();
              return NULL;
            }
        }
      return tmp;
    }
  else
    {
      PyErr_SetString(PyExc_TypeError,"convertPyToNewIntArr : not a list");
      PyErr_Print();
      return 0;
    }
}
