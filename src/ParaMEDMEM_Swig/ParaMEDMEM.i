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

#define MEDCOUPLING_EXPORT
#define INTERPKERNEL_EXPORT

%include "MEDCouplingCommon.i"

%include std_set.i
%include std_string.i

%template() std::set<int>;

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

using namespace INTERP_KERNEL;
using namespace MEDCoupling;
using namespace ICoCo;
%}

%include "InterpolationOptions.hxx"
%include "CommInterface.hxx"
%include "ProcessorGroup.hxx"
%include "DECOptions.hxx"
%include "ParaMESH.hxx"
%include "ParaFIELD.hxx"
%include "MPIProcessorGroup.hxx"
%include "ComponentTopology.hxx"
%include "DEC.hxx"
%include "DisjointDEC.hxx"
%include "InterpKernelDEC.hxx"
%include "StructuredCoincidentDEC.hxx"

%include "ICoCoField.hxx"
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

%pythoncode %{
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
%}

%include "MEDCouplingFinalize.i"
