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

#ifndef __MEDCOUPLINGTYPEMAPS_I__
#define __MEDCOUPLINGTYPEMAPS_I__

#include "MEDCouplingDataArrayTypemaps.i"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingIMesh.hxx"
#include "MEDCouplingCurveLinearMesh.hxx"
#include "MEDCouplingMappedExtrudedMesh.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
#include "MEDCouplingMultiFields.hxx"
#include "MEDCouplingPartDefinition.hxx"
#include "MEDCouplingCartesianAMRMesh.hxx"

static PyObject *convertMesh(MEDCoupling::MEDCouplingMesh *mesh, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!mesh)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDCoupling::MEDCouplingUMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_MEDCoupling__MEDCouplingUMesh,owner);
  if(dynamic_cast<MEDCoupling::MEDCoupling1SGTUMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_MEDCoupling__MEDCoupling1SGTUMesh,owner);
  if(dynamic_cast<MEDCoupling::MEDCoupling1DGTUMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_MEDCoupling__MEDCoupling1DGTUMesh,owner);
  if(dynamic_cast<MEDCoupling::MEDCouplingMappedExtrudedMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_MEDCoupling__MEDCouplingMappedExtrudedMesh,owner);
  if(dynamic_cast<MEDCoupling::MEDCouplingCMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_MEDCoupling__MEDCouplingCMesh,owner);
  if(dynamic_cast<MEDCoupling::MEDCouplingCurveLinearMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_MEDCoupling__MEDCouplingCurveLinearMesh,owner);
  if(dynamic_cast<MEDCoupling::MEDCouplingIMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_MEDCoupling__MEDCouplingIMesh,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of mesh on downcast !");
  return ret;
}

static PyObject *convertFieldDiscretization(MEDCoupling::MEDCouplingFieldDiscretization *fd, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!fd)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDCoupling::MEDCouplingFieldDiscretizationP0 *>(fd))
    ret=SWIG_NewPointerObj(reinterpret_cast<void*>(fd),SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDiscretizationP0,owner);
  if(dynamic_cast<MEDCoupling::MEDCouplingFieldDiscretizationP1 *>(fd))
    ret=SWIG_NewPointerObj(reinterpret_cast<void*>(fd),SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDiscretizationP1,owner);
  if(dynamic_cast<MEDCoupling::MEDCouplingFieldDiscretizationGauss *>(fd))
    ret=SWIG_NewPointerObj(reinterpret_cast<void*>(fd),SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDiscretizationGauss,owner);
  if(dynamic_cast<MEDCoupling::MEDCouplingFieldDiscretizationGaussNE *>(fd))
    ret=SWIG_NewPointerObj(reinterpret_cast<void*>(fd),SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDiscretizationGaussNE,owner);
  if(dynamic_cast<MEDCoupling::MEDCouplingFieldDiscretizationKriging *>(fd))
    ret=SWIG_NewPointerObj(reinterpret_cast<void*>(fd),SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDiscretizationKriging,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of field discretization on downcast !");
  return ret;
}

static PyObject* convertMultiFields(MEDCoupling::MEDCouplingMultiFields *mfs, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!mfs)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDCoupling::MEDCouplingFieldOverTime *>(mfs))
    ret=SWIG_NewPointerObj((void*)mfs,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldOverTime,owner);
  else
    ret=SWIG_NewPointerObj((void*)mfs,SWIGTYPE_p_MEDCoupling__MEDCouplingMultiFields,owner);
  return ret;
}

static PyObject *convertPartDefinition(MEDCoupling::PartDefinition *pd, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!pd)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDCoupling::DataArrayPartDefinition *>(pd))
    ret=SWIG_NewPointerObj((void*)pd,SWIGTYPE_p_MEDCoupling__DataArrayPartDefinition,owner);
  else
    ret=SWIG_NewPointerObj((void*)pd,SWIGTYPE_p_MEDCoupling__SlicePartDefinition,owner);
  return ret;
}

static PyObject *convertCartesianAMRMesh(MEDCoupling::MEDCouplingCartesianAMRMeshGen *mesh, int owner) throw(INTERP_KERNEL::Exception)
{
  if(!mesh)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDCoupling::MEDCouplingCartesianAMRMeshSub *>(mesh))
    {
      return SWIG_NewPointerObj(reinterpret_cast<void*>(mesh),SWIGTYPE_p_MEDCoupling__MEDCouplingCartesianAMRMeshSub,owner);
    }
  if(dynamic_cast<MEDCoupling::MEDCouplingCartesianAMRMesh *>(mesh))
    {
      return SWIG_NewPointerObj(reinterpret_cast<void*>(mesh),SWIGTYPE_p_MEDCoupling__MEDCouplingCartesianAMRMesh,owner);
    }
  throw INTERP_KERNEL::Exception("convertCartesianAMRMesh wrap : unrecognized type of cartesian AMR mesh !");
}

static PyObject *convertDataForGodFather(MEDCoupling::MEDCouplingDataForGodFather *data, int owner) throw(INTERP_KERNEL::Exception)
{
  if(!data)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDCoupling::MEDCouplingAMRAttribute *>(data))
    {
      return SWIG_NewPointerObj(reinterpret_cast<void*>(data),SWIGTYPE_p_MEDCoupling__MEDCouplingAMRAttribute,owner);
    }
  throw INTERP_KERNEL::Exception("convertDataForGodFather wrap : unrecognized data type for AMR !");
}

static PyObject *convertCartesianAMRPatch(MEDCoupling::MEDCouplingCartesianAMRPatchGen *patch, int owner) throw(INTERP_KERNEL::Exception)
{
  if(!patch)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<MEDCoupling::MEDCouplingCartesianAMRPatchGF *>(patch))
    {
      return SWIG_NewPointerObj(reinterpret_cast<void*>(patch),SWIGTYPE_p_MEDCoupling__MEDCouplingCartesianAMRPatchGF,owner);
    }
  if(dynamic_cast<MEDCoupling::MEDCouplingCartesianAMRPatch *>(patch))
    {
      return SWIG_NewPointerObj(reinterpret_cast<void*>(patch),SWIGTYPE_p_MEDCoupling__MEDCouplingCartesianAMRPatch,owner);
    }
  throw INTERP_KERNEL::Exception("convertCartesianAMRPatch wrap : unrecognized type of cartesian AMR patch !");
}

static MEDCoupling::MEDCouplingFieldDouble *MEDCoupling_MEDCouplingFieldDouble___add__Impl(MEDCoupling::MEDCouplingFieldDouble *self, PyObject *obj) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__add__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
  const char msg2[]="in MEDCouplingFieldDouble.__add__ : self field has no Array of values set !";
  void *argp;
  //
  if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,0|0)))
    {
      MEDCoupling::MEDCouplingFieldDouble *other=reinterpret_cast< MEDCoupling::MEDCouplingFieldDouble * >(argp);
      if(other)
        return (*self)+(*other);
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  //
  double val;
  MEDCoupling::DataArrayDouble *a;
  MEDCoupling::DataArrayDoubleTuple *aa;
  std::vector<double> bb;
  int sw;
  convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
  switch(sw)
    {
    case 1:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=self->getArray()->deepCopy();
        ret->applyLin(1.,val);
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 2:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=MEDCoupling::DataArrayDouble::Add(self->getArray(),a);
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 3:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=MEDCoupling::DataArrayDouble::Add(self->getArray(),aaa);
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 4:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> aaa=MEDCoupling::DataArrayDouble::New(); aaa->useArray(&bb[0],false,MEDCoupling::CPP_DEALLOC,1,(int)bb.size());
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=MEDCoupling::DataArrayDouble::Add(self->getArray(),aaa);
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    default:
      { throw INTERP_KERNEL::Exception(msg); }
    }
}

static MEDCoupling::MEDCouplingFieldDouble *MEDCoupling_MEDCouplingFieldDouble___radd__Impl(MEDCoupling::MEDCouplingFieldDouble *self, PyObject *obj) throw(INTERP_KERNEL::Exception)
{
  return MEDCoupling_MEDCouplingFieldDouble___add__Impl(self,obj);
}

static MEDCoupling::MEDCouplingFieldDouble *MEDCoupling_MEDCouplingFieldDouble___rsub__Impl(MEDCoupling::MEDCouplingFieldDouble *self, PyObject *obj) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__rsub__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
  const char msg2[]="in MEDCouplingFieldDouble.__rsub__ : self field has no Array of values set !";
  void *argp;
  //
  if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,0|0)))
    {
      MEDCoupling::MEDCouplingFieldDouble *other=reinterpret_cast< MEDCoupling::MEDCouplingFieldDouble * >(argp);
      if(other)
        return (*other)-(*self);
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  //
  double val;
  MEDCoupling::DataArrayDouble *a;
  MEDCoupling::DataArrayDoubleTuple *aa;
  std::vector<double> bb;
  int sw;
  convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
  switch(sw)
    {
    case 1:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=self->getArray()->deepCopy();
        ret->applyLin(-1.,val);
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 2:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=MEDCoupling::DataArrayDouble::Substract(a,self->getArray());
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 3:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=MEDCoupling::DataArrayDouble::Substract(aaa,self->getArray());
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 4:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> aaa=MEDCoupling::DataArrayDouble::New(); aaa->useArray(&bb[0],false,MEDCoupling::CPP_DEALLOC,1,(int)bb.size());
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=MEDCoupling::DataArrayDouble::Substract(aaa,self->getArray());
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    default:
      { throw INTERP_KERNEL::Exception(msg); }
    }
}

static MEDCoupling::MEDCouplingFieldDouble *MEDCoupling_MEDCouplingFieldDouble___mul__Impl(MEDCoupling::MEDCouplingFieldDouble *self, PyObject *obj) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__mul__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
  const char msg2[]="in MEDCouplingFieldDouble.__mul__ : self field has no Array of values set !";
  void *argp;
  //
  if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,0|0)))
    {
      MEDCoupling::MEDCouplingFieldDouble *other=reinterpret_cast< MEDCoupling::MEDCouplingFieldDouble * >(argp);
      if(other)
        return (*self)*(*other);
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  //
  double val;
  MEDCoupling::DataArrayDouble *a;
  MEDCoupling::DataArrayDoubleTuple *aa;
  std::vector<double> bb;
  int sw;
  convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
  switch(sw)
    {
    case 1:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=self->getArray()->deepCopy();
        ret->applyLin(val,0.);
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 2:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=MEDCoupling::DataArrayDouble::Multiply(self->getArray(),a);
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 3:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=MEDCoupling::DataArrayDouble::Multiply(self->getArray(),aaa);
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 4:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> aaa=MEDCoupling::DataArrayDouble::New(); aaa->useArray(&bb[0],false,MEDCoupling::CPP_DEALLOC,1,(int)bb.size());
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=MEDCoupling::DataArrayDouble::Multiply(self->getArray(),aaa);
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    default:
      { throw INTERP_KERNEL::Exception(msg); }
    }
}

MEDCoupling::MEDCouplingFieldDouble *MEDCoupling_MEDCouplingFieldDouble___rmul__Impl(MEDCoupling::MEDCouplingFieldDouble *self, PyObject *obj) throw(INTERP_KERNEL::Exception)
{
  return MEDCoupling_MEDCouplingFieldDouble___mul__Impl(self,obj);
}

MEDCoupling::MEDCouplingFieldDouble *MEDCoupling_MEDCouplingFieldDouble___rdiv__Impl(MEDCoupling::MEDCouplingFieldDouble *self, PyObject *obj) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__rdiv__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
  const char msg2[]="in MEDCouplingFieldDouble.__div__ : self field has no Array of values set !";
  void *argp;
  //
  if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_MEDCoupling__MEDCouplingFieldDouble,0|0)))
    {
      MEDCoupling::MEDCouplingFieldDouble *other=reinterpret_cast< MEDCoupling::MEDCouplingFieldDouble * >(argp);
      if(other)
        return (*other)/(*self);
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  //
  double val;
  MEDCoupling::DataArrayDouble *a;
  MEDCoupling::DataArrayDoubleTuple *aa;
  std::vector<double> bb;
  int sw;
  convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
  switch(sw)
    {
    case 1:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=self->getArray()->deepCopy();
        ret->applyInv(val);
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 2:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=MEDCoupling::DataArrayDouble::Divide(a,self->getArray());
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 3:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=MEDCoupling::DataArrayDouble::Divide(aaa,self->getArray());
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 4:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> aaa=MEDCoupling::DataArrayDouble::New(); aaa->useArray(&bb[0],false,MEDCoupling::CPP_DEALLOC,1,(int)bb.size());
        MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> ret=MEDCoupling::DataArrayDouble::Divide(aaa,self->getArray());
        MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    default:
      { throw INTERP_KERNEL::Exception(msg); }
    }
}

#endif
