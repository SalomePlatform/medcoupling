// Copyright (C) 2007-2019  CEA/DEN, EDF R&D
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

%include std_vector.i
%include std_string.i

%include "MEDCouplingCommon.i"

%{
#ifndef WIN32
#include "MEDCouplingMemArray.txx"
#endif
#include "MCAuto.hxx"
#include "MEDCouplingDataArrayTypemaps.i"

using namespace MEDCoupling;
using namespace INTERP_KERNEL;
%}

%template(ivec) std::vector<int>;
%template(dvec) std::vector<double>;
%template(svec) std::vector<std::string>;

#ifdef WITH_NUMPY
%init %{ import_array(); %}
#endif

%init %{ initializeMe_renumber(); %}

%feature("autodoc", "1");
%feature("docstring");

%nodefaultctor;

%rename (InterpKernelException) INTERP_KERNEL::Exception;

%include "MEDCouplingRefCountObject.i"
%include "MEDCouplingMemArray.i"

%{
  void initializeMe_renumber()
  {// AGY : here initialization of C++ traits in MEDCouplingDataArrayTypemaps.i for code factorization. Awful, I know, but no other solutions.
    SWIGTITraits<double>::TI=SWIGTYPE_p_MEDCoupling__DataArrayDouble;
    SWIGTITraits<float>::TI=SWIGTYPE_p_MEDCoupling__DataArrayFloat;
  }
%}

%include "MEDRenumberImpl.i"

%pythoncode %{
import os
__filename=os.environ.get('PYTHONSTARTUP')
if __filename and os.path.isfile(__filename):
  exec(open(__filename).read())
  pass
%}
