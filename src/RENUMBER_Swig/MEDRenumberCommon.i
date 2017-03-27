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
// Author : Anthony Geay (CEA/DEN)

%module MEDRenumber

%include std_vector.i
%include std_string.i

%{
#include "MEDCouplingMemArray.txx"
#include "MCAuto.hxx"
#include "MEDCouplingDataArrayTypemaps.i"

#include "RenumberingFactory.hxx"
#include "RENUMBER_Renumbering.hxx"

using namespace MEDCoupling;
using namespace INTERP_KERNEL;
 using namespace MED_RENUMBER;
%}

%template(ivec) std::vector<int>;
%template(dvec) std::vector<double>;
%template(svec) std::vector<std::string>;

#ifdef WITH_NUMPY
%init %{ import_array(); %}
#endif

%init %{ initializeMe(); %}

%feature("autodoc", "1");
%feature("docstring");

%newobject MED_RENUMBER::RenumberingFactory;

%nodefaultctor;

%rename (InterpKernelException) INTERP_KERNEL::Exception;

%include "MEDCouplingRefCountObject.i"
%include "MEDCouplingMemArray.i"

%{
  void initializeMe()
  {// AGY : here initialization of C++ traits in MEDCouplingDataArrayTypemaps.i for code factorization. Awful, I know, but no other solutions.
    SWIGTITraits<double>::TI=SWIGTYPE_p_MEDCoupling__DataArrayDouble;
    SWIGTITraits<float>::TI=SWIGTYPE_p_MEDCoupling__DataArrayFloat;
  }
%}

class Renumbering
{
public:
  %extend
  {
    virtual PyObject *renumber(const MEDCoupling::DataArrayInt *graph, const MEDCoupling::DataArrayInt *index_graph) throw(INTERP_KERNEL::Exception)
    {
      if(!graph || !index_graph)
        throw INTERP_KERNEL::Exception("wrap of Renumbering::renumber : One of the input arrays is NULL !");
      if(!graph->isAllocated() || !index_graph->isAllocated())
        throw INTERP_KERNEL::Exception("wrap of Renumbering::renumber : One of the input arrays is not allocated !");
      MEDCoupling::DataArrayInt *out0(0),*out1(0);
      self->renumber(graph->begin(),index_graph->begin(),index_graph->getNumberOfTuples()-1,out0,out1);
      PyObject *ret=PyTuple_New(2);
      PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(out0),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
      PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(out1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
      return ret;
    }
  }
  virtual ~Renumbering();
};

namespace MED_RENUMBER
{
  Renumbering *RenumberingFactory(const std::string& s) throw(INTERP_KERNEL::Exception);
}

%inline
{
  std::vector<std::string> RenumberAvailableMethods()throw(INTERP_KERNEL::Exception)
  {
    std::vector<std::string> ret;
#ifdef HAS_BOOST
    ret.push_back(std::string("BOOST"));
#endif
#ifdef HAS_METIS
    ret.push_back(std::string("METIS"));
#endif
    return ret;
  }
}

%pythoncode %{
import os
__filename=os.environ.get('PYTHONSTARTUP')
if __filename and os.path.isfile(__filename):
  exec(open(__filename).read())
  pass
%}
