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
// Author : Anthony Geay (EDF R&D)

%{
#include "RenumberingFactory.hxx"
#include "RENUMBER_Renumbering.hxx"
  
using namespace MED_RENUMBER;
%}

%newobject MED_RENUMBER::RenumberingFactory;

class Renumbering
{
public:
  %extend
  {
    virtual PyObject *renumber(const MEDCoupling::DataArrayIdType *graph, const MEDCoupling::DataArrayIdType *index_graph)
    {
      if(!graph || !index_graph)
        throw INTERP_KERNEL::Exception("wrap of Renumbering::renumber : One of the input arrays is NULL !");
      if(!graph->isAllocated() || !index_graph->isAllocated())
        throw INTERP_KERNEL::Exception("wrap of Renumbering::renumber : One of the input arrays is not allocated !");
      MEDCoupling::DataArrayIdType *out0(0),*out1(0);
      self->renumber(graph->begin(),index_graph->begin(),index_graph->getNumberOfTuples()-1,out0,out1);
      PyObject *ret=PyTuple_New(2);
      PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(out0),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
      PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(out1),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
      return ret;
    }
  }
  virtual ~Renumbering();
};

namespace MED_RENUMBER
{
  Renumbering *RenumberingFactory(const std::string& s);
  std::vector<std::string> RenumberAvailableMethods();
  std::vector<std::string> AllRenumberMethods();
}
