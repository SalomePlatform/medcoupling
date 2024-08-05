%module ShapeRecogn

%include "std_string.i"
%include "MEDCouplingCommon.i"

%{
#include "ShapeRecognMesh.hxx"
#include "ShapeRecognMeshBuilder.hxx"
using namespace MEDCoupling;
%}

%ignore getAreas() const;
%ignore getNodes() const;
%include "ShapeRecognMesh.hxx"
%include "ShapeRecognMeshBuilder.hxx"
