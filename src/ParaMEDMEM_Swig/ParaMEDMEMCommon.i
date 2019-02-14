// Copyright (C) 2017-2019  CEA/DEN, EDF R&D
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

%include std_set.i

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
