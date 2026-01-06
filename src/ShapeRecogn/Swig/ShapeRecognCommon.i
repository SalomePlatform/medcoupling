// Copyright (C) 2024-2026  CEA, EDF
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

%include "std_string.i"
%include "MEDCouplingCommon.i"

%{
#include "MEDCouplingMemArray.txx"
#include "MCAuto.hxx"
#include "MEDCouplingDataArrayTypemaps.i"

using namespace MEDCoupling;
using namespace INTERP_KERNEL;
%}

#ifdef WITH_NUMPY
%init %{ import_array(); %}
#endif

%init %{ initializeMe_shape_recogn(); %}

%feature("autodoc", "1");
%feature("docstring");

%nodefaultctor;

%rename (InterpKernelException) INTERP_KERNEL::Exception;

%{
  void initializeMe_shape_recogn()
  {// AGY : here initialization of C++ traits in MEDCouplingDataArrayTypemaps.i for code factorization. Awful, I know, but no other solutions.
    SWIGTITraits<double>::TI=SWIGTYPE_p_MEDCoupling__DataArrayDouble;
    SWIGTITraits<float>::TI=SWIGTYPE_p_MEDCoupling__DataArrayFloat;
  }
%}

%include "ShapeRecognImpl.i"

%pythoncode %{
def AreasListOfPropertiesForPrimitive(cls, primitiveType):
  # Cone
  if primitiveType == 3:
    return ["Radius","Angle","Axis","AxisPoint","Apex"]
  # Plan
  elif primitiveType == 0:
    return ["Normal","AffinePoint"]
  # Cylinder
  elif primitiveType == 2:
    return ["Radius","AxisPoint","Axis"]
  # Tore
  elif primitiveType == 4:
    return ["Center","Radius","MinorRadius"]
  # Sphere
  elif primitiveType == 1:
    return ["Center","Radius"]
  raise RuntimeError(f"Unrecognized primitive {primitiveType}")

def AreasReprOfArea(self, areaId ):
  def evaluateProps(self, areaId, listOfProps ):
    def evaluateProp(self, areaId, prop):
      return eval(f"self.get{prop}({areaId})")
    return [ f"{prop} = {evaluateProp(self,areaId,prop)}" for prop in listOfProps]
  primitiveType = self.getPrimitiveType( areaId )
  ret = [ ConvertPrimitiveToString( primitiveType ) ]
  ret += evaluateProps(self,areaId,Areas.ListOfPropertiesForPrimitive(primitiveType))
  return "\n".join( ret )

Areas.ListOfPropertiesForPrimitive = classmethod( AreasListOfPropertiesForPrimitive )
Areas.reprOfArea = AreasReprOfArea

import os
__filename=os.environ.get('PYTHONSTARTUP')
if __filename and os.path.isfile(__filename):
  with open(__filename) as __fp:
        exec(__fp.read())
%}
