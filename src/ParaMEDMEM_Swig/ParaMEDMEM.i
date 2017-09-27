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

%module ParaMEDMEM

%include "ParaMEDMEM.typemap"
%include "MEDCouplingCommon.i"

%{
#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "Topology.hxx"
#include "MPIProcessorGroup.hxx"
#include "DEC.hxx"
#include "InterpKernelDEC.hxx"
#include "NonCoincidentDEC.hxx"
#include "StructuredCoincidentDEC.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "ICoCoMEDField.hxx"
#include "ComponentTopology.hxx"

#include <mpi.h>

using namespace MEDCoupling;
using namespace ICoCo;
      
enum mpi_constants { mpi_comm_world, mpi_comm_self, mpi_double, mpi_int };
%}

%include "CommInterface.hxx"
%include "ProcessorGroup.hxx"
%include "DECOptions.hxx"
%include "ParaMESH.hxx"
%include "ParaFIELD.hxx"
%include "MPIProcessorGroup.hxx"
%include "ComponentTopology.hxx"
%include "DEC.hxx"
%include "InterpKernelDEC.hxx"
%include "StructuredCoincidentDEC.hxx"

%rename(ICoCoMEDField) ICoCo::MEDField;
%include "ICoCoMEDField.hxx"

%nodefaultctor;

/* This object can be used only if MED_ENABLE_FVM is defined*/
#ifdef MED_ENABLE_FVM
class NonCoincidentDEC : public DEC
{
public:
  NonCoincidentDEC(ProcessorGroup& source, ProcessorGroup& target);
};
#endif

%extend MEDCoupling::ParaMESH
{
  PyObject *getGlobalNumberingCell2() const
  {
    const int *tmp=self->getGlobalNumberingCell();
    int size=self->getCellMesh()->getNumberOfCells();
    PyObject *ret=PyList_New(size);
    for(int i=0;i<size;i++)
      PyList_SetItem(ret,i,PyInt_FromLong(tmp[i])); 
    return ret;
  }

  PyObject *getGlobalNumberingFace2() const
  {
    const int *tmp=self->getGlobalNumberingFace();
    int size=self->getFaceMesh()->getNumberOfCells();
    PyObject *ret=PyList_New(size);
    for(int i=0;i<size;i++)
      PyList_SetItem(ret,i,PyInt_FromLong(tmp[i])); 
    return ret;
  }

  PyObject *getGlobalNumberingNode2() const
  {
    const int *tmp=self->getGlobalNumberingNode();
    int size=self->getCellMesh()->getNumberOfNodes();
    PyObject *ret=PyList_New(size);
    for(int i=0;i<size;i++)
      PyList_SetItem(ret,i,PyInt_FromLong(tmp[i])); 
    return ret;
  }
}

//=============================================================================================
// Interface for MPI-realization-specific constants like MPI_COMM_WORLD.
//
// Type and values of constants like MPI_COMM_WORLD depends on MPI realization
// and usually such constants actually are macros. To have such symbols in python
// and translate them into correct values we use the following technique.
// We define some constants (enum mpi_constants) and map them into real MPI values
// using typemaps, and we create needed python symbols equal to 'mpi_constants'
// via %pythoncode directive.

// Constants corresponding to similar MPI definitions
enum mpi_constants { mpi_comm_world, mpi_comm_self, mpi_double, mpi_int };

// Map mpi_comm_world and mpi_comm_self -> MPI_COMM_WORLD and MPI_COMM_SELF
%typemap(in) MPI_Comm
{ 
  switch (PyInt_AsLong($input))
    {
    case mpi_comm_world: $1 = MPI_COMM_WORLD; break;
    case mpi_comm_self:  $1 = MPI_COMM_SELF;  break;
    default:
      PyErr_SetString(PyExc_TypeError,"unexpected value of MPI_Comm");
      return NULL;
    }
}
// Map mpi_double and mpi_int -> MPI_DOUBLE and MPI_INT
%typemap(in) MPI_Datatype
{ 
  switch (PyInt_AsLong($input))
    {
    case mpi_double:     $1 = MPI_DOUBLE;     break;
    case mpi_int:        $1 = MPI_INT;        break;
    default:
      PyErr_SetString(PyExc_TypeError,"unexpected value of MPI_Datatype");
      return NULL;
    }
}
// The following code gets inserted into the result python file:
// create needed python symbols
%pythoncode %{
MPI_COMM_WORLD = mpi_comm_world
MPI_COMM_SELF  = mpi_comm_self
MPI_DOUBLE     = mpi_double
MPI_INT        = mpi_int
%}
//=============================================================================================

// ==============
// MPI_Comm_size
// ==============
%inline %{ PyObject* MPI_Comm_size(MPI_Comm comm)
  {
    int res = 0;
    int err = MPI_Comm_size(comm, &res);
    if ( err != MPI_SUCCESS )
      {
        PyErr_SetString(PyExc_RuntimeError,"Erorr in MPI_Comm_size()");
        return NULL;
      }
    return PyInt_FromLong( res );
  } %}

// ==============
// MPI_Comm_rank
// ==============
%inline %{ PyObject* MPI_Comm_rank(MPI_Comm comm)
  {
    int res = 0;
    int err = MPI_Comm_rank(comm, &res);
    if ( err != MPI_SUCCESS )
      {
        PyErr_SetString(PyExc_RuntimeError,"Erorr in MPI_Comm_rank()");
        return NULL;
      }
    return PyInt_FromLong( res );
  } 
  %}

int MPI_Init(int *argc, char ***argv );
int MPI_Barrier(MPI_Comm comm);
int MPI_Finalize();

// ==========
// MPI_Bcast
// ==========

%inline %{ PyObject* MPI_Bcast(PyObject* buffer, int nb, MPI_Datatype type, int root, MPI_Comm c)
  {
    // buffer must be a list
    if (!PyList_Check(buffer))
      {
        PyErr_SetString(PyExc_TypeError, "buffer is expected to be a list");
        return NULL;
      }
    // check list size
    int aSize = PyList_Size(buffer);
    if ( aSize != nb )
      {
        std::ostringstream stream; stream << "buffer is expected to be of size " << nb;
        PyErr_SetString(PyExc_ValueError, stream.str().c_str());
        return NULL;
      }
    // allocate and fill a buffer
    void* aBuf = 0;
    int* intBuf = 0;
    double* dblBuf = 0;
    if ( type == MPI_DOUBLE )
      {
        aBuf = (void*) ( dblBuf = new double[ nb ] );
        for ( int i = 0; i < aSize; ++i )
          dblBuf[i] = PyFloat_AS_DOUBLE( PyList_GetItem( buffer, i ));
      }
    else if ( type == MPI_INT )
      {
        aBuf = (void*) ( intBuf = new int[ nb ] );
        for ( int i = 0; i < aSize; ++i )
          intBuf[i] = int( PyInt_AS_LONG( PyList_GetItem( buffer, i )));
      }
    else
      {
        PyErr_SetString(PyExc_TypeError, "Only MPI_DOUBLE and MPI_INT supported");
        return NULL;
      }
    // call MPI_Bcast
    int err = MPI_Bcast(aBuf, nb, type, root, c);
    // treat error
    if ( err != MPI_SUCCESS )
      {
        PyErr_SetString(PyExc_RuntimeError,"Erorr in MPI_Bcast()");
        delete [] intBuf; delete [] dblBuf;
        return NULL;
      }
    // put recieved data into the list
    int pyerr = 0;
    if ( type == MPI_DOUBLE )
      {
        for ( int i = 0; i < aSize && !pyerr; ++i )
          pyerr = PyList_SetItem(buffer, i, PyFloat_FromDouble( dblBuf[i] ));
        delete [] dblBuf;
      }
    else
      {
        for ( int i = 0; i < aSize && !pyerr; ++i )
          pyerr = PyList_SetItem(buffer, i, PyInt_FromLong( intBuf[i] ));
        delete [] intBuf;
      }
    if ( pyerr )
      {
        PyErr_SetString(PyExc_RuntimeError, "Error of PyList_SetItem()");
        return NULL;
      }
    return PyInt_FromLong( err );

  }
  %}

%pythoncode %{
def MEDCouplingDataArrayDoublenew(cls,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayDouble____new___(cls,args)
def MEDCouplingDataArrayDoubleIadd(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayDouble____iadd___(self, self, *args)
def MEDCouplingDataArrayDoubleIsub(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayDouble____isub___(self, self, *args)
def MEDCouplingDataArrayDoubleImul(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayDouble____imul___(self, self, *args)
def MEDCouplingDataArrayDoubleIdiv(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayDouble____idiv___(self, self, *args)
def MEDCouplingDataArrayDoubleIpow(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayDouble____ipow___(self, self, *args)
def MEDCouplingFieldDoublenew(cls,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.MEDCouplingFieldDouble____new___(cls,args)
def MEDCouplingFieldDoubleIadd(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.MEDCouplingFieldDouble____iadd___(self, self, *args)
def MEDCouplingFieldDoubleIsub(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.MEDCouplingFieldDouble____isub___(self, self, *args)
def MEDCouplingFieldDoubleImul(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.MEDCouplingFieldDouble____imul___(self, self, *args)
def MEDCouplingFieldDoubleIdiv(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.MEDCouplingFieldDouble____idiv___(self, self, *args)
def MEDCouplingFieldDoubleIpow(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.MEDCouplingFieldDouble____ipow___(self, self, *args)
def MEDCouplingFieldIntnew(cls,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.MEDCouplingFieldInt____new___(cls,args)
def MEDCouplingFieldFloatnew(cls,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.MEDCouplingFieldFloat____new___(cls,args)
def MEDCouplingDataArrayIntnew(cls,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayInt____new___(cls,args)
def MEDCouplingDataArrayIntIadd(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayInt____iadd___(self, self, *args)
def MEDCouplingDataArrayIntIsub(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayInt____isub___(self, self, *args)
def MEDCouplingDataArrayIntImul(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayInt____imul___(self, self, *args)
def MEDCouplingDataArrayIntIdiv(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayInt____idiv___(self, self, *args)
def MEDCouplingDataArrayIntImod(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayInt____imod___(self, self, *args)
def MEDCouplingDataArrayIntIpow(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayInt____ipow___(self, self, *args)
def MEDCouplingDataArrayBytenew(cls,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayByte____new___(cls,args)
def MEDCouplingDataArrayFloatnew(cls,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayFloat____new___(cls,args)
def MEDCouplingDataArrayFloatIadd(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayFloat____iadd___(self, self, *args)
def MEDCouplingDataArrayFloatIsub(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayFloat____isub___(self, self, *args)
def MEDCouplingDataArrayFloatImul(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayFloat____imul___(self, self, *args)
def MEDCouplingDataArrayFloatIdiv(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayFloat____idiv___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIadd(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayDoubleTuple____iadd___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIsub(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayDoubleTuple____isub___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleImul(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayDoubleTuple____imul___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIdiv(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayDoubleTuple____idiv___(self, self, *args)
def MEDCouplingDataArrayIntTupleIadd(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayIntTuple____iadd___(self, self, *args)
def MEDCouplingDataArrayIntTupleIsub(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayIntTuple____isub___(self, self, *args)
def MEDCouplingDataArrayIntTupleImul(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayIntTuple____imul___(self, self, *args)
def MEDCouplingDataArrayIntTupleIdiv(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayIntTuple____idiv___(self, self, *args)
def MEDCouplingDataArrayIntTupleImod(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DataArrayIntTuple____imod___(self, self, *args)
def ParaMEDMEMDenseMatrixIadd(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DenseMatrix____iadd___(self, self, *args)
def ParaMEDMEMDenseMatrixIsub(self,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.DenseMatrix____isub___(self, self, *args)
def MEDCouplingUMeshnew(cls,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.MEDCouplingUMesh____new___(cls,args)
def MEDCoupling1DGTUMeshnew(cls,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.MEDCoupling1DGTUMesh____new___(cls,args)
def MEDCoupling1SGTUMeshnew(cls,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.MEDCoupling1SGTUMesh____new___(cls,args)
def MEDCouplingCurveLinearMeshnew(cls,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.MEDCouplingCurveLinearMesh____new___(cls,args)
def MEDCouplingCMeshnew(cls,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.MEDCouplingCMesh____new___(cls,args)
def MEDCouplingIMeshnew(cls,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.MEDCouplingIMesh____new___(cls,args)
def MEDCouplingExtrudedMeshnew(cls,*args):
    import _ParaMEDMEM
    return _ParaMEDMEM.MEDCouplingMappedExtrudedMesh____new___(cls,args)
%}

%include "MEDCouplingFinalize.i"
