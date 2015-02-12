// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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

#include "MEDCouplingDataArrayTypemaps.i"

static PyObject *convertMesh(ParaMEDMEM::MEDCouplingMesh *mesh, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!mesh)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<ParaMEDMEM::MEDCouplingUMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDCouplingUMesh,owner);
  if(dynamic_cast<ParaMEDMEM::MEDCoupling1SGTUMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDCoupling1SGTUMesh,owner);
  if(dynamic_cast<ParaMEDMEM::MEDCoupling1DGTUMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDCoupling1DGTUMesh,owner);
  if(dynamic_cast<ParaMEDMEM::MEDCouplingExtrudedMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDCouplingExtrudedMesh,owner);
  if(dynamic_cast<ParaMEDMEM::MEDCouplingCMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDCouplingCMesh,owner);
  if(dynamic_cast<ParaMEDMEM::MEDCouplingCurveLinearMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDCouplingCurveLinearMesh,owner);
  if(dynamic_cast<ParaMEDMEM::MEDCouplingIMesh *>(mesh))
    ret=SWIG_NewPointerObj((void*)mesh,SWIGTYPE_p_ParaMEDMEM__MEDCouplingIMesh,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of mesh on downcast !");
  return ret;
}

static PyObject *convertFieldDiscretization(ParaMEDMEM::MEDCouplingFieldDiscretization *fd, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!fd)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<ParaMEDMEM::MEDCouplingFieldDiscretizationP0 *>(fd))
    ret=SWIG_NewPointerObj(reinterpret_cast<void*>(fd),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDiscretizationP0,owner);
  if(dynamic_cast<ParaMEDMEM::MEDCouplingFieldDiscretizationP1 *>(fd))
    ret=SWIG_NewPointerObj(reinterpret_cast<void*>(fd),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDiscretizationP1,owner);
  if(dynamic_cast<ParaMEDMEM::MEDCouplingFieldDiscretizationGauss *>(fd))
    ret=SWIG_NewPointerObj(reinterpret_cast<void*>(fd),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDiscretizationGauss,owner);
  if(dynamic_cast<ParaMEDMEM::MEDCouplingFieldDiscretizationGaussNE *>(fd))
    ret=SWIG_NewPointerObj(reinterpret_cast<void*>(fd),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDiscretizationGaussNE,owner);
  if(dynamic_cast<ParaMEDMEM::MEDCouplingFieldDiscretizationKriging *>(fd))
    ret=SWIG_NewPointerObj(reinterpret_cast<void*>(fd),SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDiscretizationKriging,owner);
  if(!ret)
    throw INTERP_KERNEL::Exception("Not recognized type of field discretization on downcast !");
  return ret;
}

static PyObject* convertMultiFields(ParaMEDMEM::MEDCouplingMultiFields *mfs, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!mfs)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<ParaMEDMEM::MEDCouplingFieldOverTime *>(mfs))
    ret=SWIG_NewPointerObj((void*)mfs,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldOverTime,owner);
  else
    ret=SWIG_NewPointerObj((void*)mfs,SWIGTYPE_p_ParaMEDMEM__MEDCouplingMultiFields,owner);
  return ret;
}

static PyObject *convertPartDefinition(ParaMEDMEM::PartDefinition *pd, int owner) throw(INTERP_KERNEL::Exception)
{
  PyObject *ret=0;
  if(!pd)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<ParaMEDMEM::DataArrayPartDefinition *>(pd))
    ret=SWIG_NewPointerObj((void*)pd,SWIGTYPE_p_ParaMEDMEM__DataArrayPartDefinition,owner);
  else
    ret=SWIG_NewPointerObj((void*)pd,SWIGTYPE_p_ParaMEDMEM__SlicePartDefinition,owner);
  return ret;
}

static PyObject *convertCartesianAMRMesh(ParaMEDMEM::MEDCouplingCartesianAMRMeshGen *mesh, int owner) throw(INTERP_KERNEL::Exception)
{
  if(!mesh)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<ParaMEDMEM::MEDCouplingCartesianAMRMeshSub *>(mesh))
    {
      return SWIG_NewPointerObj(reinterpret_cast<void*>(mesh),SWIGTYPE_p_ParaMEDMEM__MEDCouplingCartesianAMRMeshSub,owner);
    }
  if(dynamic_cast<ParaMEDMEM::MEDCouplingCartesianAMRMesh *>(mesh))
    {
      return SWIG_NewPointerObj(reinterpret_cast<void*>(mesh),SWIGTYPE_p_ParaMEDMEM__MEDCouplingCartesianAMRMesh,owner);
    }
  throw INTERP_KERNEL::Exception("convertCartesianAMRMesh wrap : unrecognized type of cartesian AMR mesh !");
}

static PyObject *convertDataForGodFather(ParaMEDMEM::MEDCouplingDataForGodFather *data, int owner) throw(INTERP_KERNEL::Exception)
{
  if(!data)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<ParaMEDMEM::MEDCouplingAMRAttribute *>(data))
    {
      return SWIG_NewPointerObj(reinterpret_cast<void*>(data),SWIGTYPE_p_ParaMEDMEM__MEDCouplingAMRAttribute,owner);
    }
  throw INTERP_KERNEL::Exception("convertDataForGodFather wrap : unrecognized data type for AMR !");
}

static PyObject *convertCartesianAMRPatch(ParaMEDMEM::MEDCouplingCartesianAMRPatchGen *patch, int owner) throw(INTERP_KERNEL::Exception)
{
  if(!patch)
    {
      Py_XINCREF(Py_None);
      return Py_None;
    }
  if(dynamic_cast<ParaMEDMEM::MEDCouplingCartesianAMRPatchGF *>(patch))
    {
      return SWIG_NewPointerObj(reinterpret_cast<void*>(patch),SWIGTYPE_p_ParaMEDMEM__MEDCouplingCartesianAMRPatchGF,owner);
    }
  if(dynamic_cast<ParaMEDMEM::MEDCouplingCartesianAMRPatch *>(patch))
    {
      return SWIG_NewPointerObj(reinterpret_cast<void*>(patch),SWIGTYPE_p_ParaMEDMEM__MEDCouplingCartesianAMRPatch,owner);
    }
  throw INTERP_KERNEL::Exception("convertCartesianAMRPatch wrap : unrecognized type of cartesian AMR patch !");
}

static ParaMEDMEM::MEDCouplingFieldDouble *ParaMEDMEM_MEDCouplingFieldDouble___add__Impl(ParaMEDMEM::MEDCouplingFieldDouble *self, PyObject *obj) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__add__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
  const char msg2[]="in MEDCouplingFieldDouble.__add__ : self field has no Array of values set !";
  void *argp;
  //
  if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,0|0)))
    {
      ParaMEDMEM::MEDCouplingFieldDouble *other=reinterpret_cast< ParaMEDMEM::MEDCouplingFieldDouble * >(argp);
      if(other)
        return (*self)+(*other);
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  //
  double val;
  ParaMEDMEM::DataArrayDouble *a;
  ParaMEDMEM::DataArrayDoubleTuple *aa;
  std::vector<double> bb;
  int sw;
  convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
  switch(sw)
    {
    case 1:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=self->getArray()->deepCpy();
        ret->applyLin(1.,val);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 2:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=ParaMEDMEM::DataArrayDouble::Add(self->getArray(),a);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 3:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=ParaMEDMEM::DataArrayDouble::Add(self->getArray(),aaa);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 4:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> aaa=ParaMEDMEM::DataArrayDouble::New(); aaa->useArray(&bb[0],false,ParaMEDMEM::CPP_DEALLOC,1,(int)bb.size());
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=ParaMEDMEM::DataArrayDouble::Add(self->getArray(),aaa);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    default:
      { throw INTERP_KERNEL::Exception(msg); }
    }
}

static ParaMEDMEM::MEDCouplingFieldDouble *ParaMEDMEM_MEDCouplingFieldDouble___radd__Impl(ParaMEDMEM::MEDCouplingFieldDouble *self, PyObject *obj) throw(INTERP_KERNEL::Exception)
{
  return ParaMEDMEM_MEDCouplingFieldDouble___add__Impl(self,obj);
}

static ParaMEDMEM::MEDCouplingFieldDouble *ParaMEDMEM_MEDCouplingFieldDouble___rsub__Impl(ParaMEDMEM::MEDCouplingFieldDouble *self, PyObject *obj) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__rsub__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
  const char msg2[]="in MEDCouplingFieldDouble.__rsub__ : self field has no Array of values set !";
  void *argp;
  //
  if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,0|0)))
    {
      ParaMEDMEM::MEDCouplingFieldDouble *other=reinterpret_cast< ParaMEDMEM::MEDCouplingFieldDouble * >(argp);
      if(other)
        return (*other)-(*self);
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  //
  double val;
  ParaMEDMEM::DataArrayDouble *a;
  ParaMEDMEM::DataArrayDoubleTuple *aa;
  std::vector<double> bb;
  int sw;
  convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
  switch(sw)
    {
    case 1:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=self->getArray()->deepCpy();
        ret->applyLin(-1.,val);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 2:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=ParaMEDMEM::DataArrayDouble::Substract(a,self->getArray());
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 3:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=ParaMEDMEM::DataArrayDouble::Substract(aaa,self->getArray());
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 4:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> aaa=ParaMEDMEM::DataArrayDouble::New(); aaa->useArray(&bb[0],false,ParaMEDMEM::CPP_DEALLOC,1,(int)bb.size());
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=ParaMEDMEM::DataArrayDouble::Substract(aaa,self->getArray());
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    default:
      { throw INTERP_KERNEL::Exception(msg); }
    }
}

static ParaMEDMEM::MEDCouplingFieldDouble *ParaMEDMEM_MEDCouplingFieldDouble___mul__Impl(ParaMEDMEM::MEDCouplingFieldDouble *self, PyObject *obj) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__mul__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
  const char msg2[]="in MEDCouplingFieldDouble.__mul__ : self field has no Array of values set !";
  void *argp;
  //
  if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,0|0)))
    {
      ParaMEDMEM::MEDCouplingFieldDouble *other=reinterpret_cast< ParaMEDMEM::MEDCouplingFieldDouble * >(argp);
      if(other)
        return (*self)*(*other);
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  //
  double val;
  ParaMEDMEM::DataArrayDouble *a;
  ParaMEDMEM::DataArrayDoubleTuple *aa;
  std::vector<double> bb;
  int sw;
  convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
  switch(sw)
    {
    case 1:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=self->getArray()->deepCpy();
        ret->applyLin(val,0.);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 2:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=ParaMEDMEM::DataArrayDouble::Multiply(self->getArray(),a);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 3:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=ParaMEDMEM::DataArrayDouble::Multiply(self->getArray(),aaa);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 4:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> aaa=ParaMEDMEM::DataArrayDouble::New(); aaa->useArray(&bb[0],false,ParaMEDMEM::CPP_DEALLOC,1,(int)bb.size());
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=ParaMEDMEM::DataArrayDouble::Multiply(self->getArray(),aaa);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    default:
      { throw INTERP_KERNEL::Exception(msg); }
    }
}

ParaMEDMEM::MEDCouplingFieldDouble *ParaMEDMEM_MEDCouplingFieldDouble___rmul__Impl(ParaMEDMEM::MEDCouplingFieldDouble *self, PyObject *obj) throw(INTERP_KERNEL::Exception)
{
  return ParaMEDMEM_MEDCouplingFieldDouble___mul__Impl(self,obj);
}

ParaMEDMEM::MEDCouplingFieldDouble *ParaMEDMEM_MEDCouplingFieldDouble___rdiv__Impl(ParaMEDMEM::MEDCouplingFieldDouble *self, PyObject *obj) throw(INTERP_KERNEL::Exception)
{
  const char msg[]="Unexpected situation in MEDCouplingFieldDouble.__rdiv__ ! Expecting a not null MEDCouplingFieldDouble or DataArrayDouble or DataArrayDoubleTuple instance, or a list of double, or a double.";
  const char msg2[]="in MEDCouplingFieldDouble.__div__ : self field has no Array of values set !";
  void *argp;
  //
  if(SWIG_IsOK(SWIG_ConvertPtr(obj,&argp,SWIGTYPE_p_ParaMEDMEM__MEDCouplingFieldDouble,0|0)))
    {
      ParaMEDMEM::MEDCouplingFieldDouble *other=reinterpret_cast< ParaMEDMEM::MEDCouplingFieldDouble * >(argp);
      if(other)
        return (*other)/(*self);
      else
        throw INTERP_KERNEL::Exception(msg);
    }
  //
  double val;
  ParaMEDMEM::DataArrayDouble *a;
  ParaMEDMEM::DataArrayDoubleTuple *aa;
  std::vector<double> bb;
  int sw;
  convertObjToPossibleCpp5(obj,sw,val,a,aa,bb);
  switch(sw)
    {
    case 1:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=self->getArray()->deepCpy();
        ret->applyInv(val);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 2:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=ParaMEDMEM::DataArrayDouble::Divide(a,self->getArray());
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 3:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> aaa=aa->buildDADouble(1,self->getNumberOfComponents());
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=ParaMEDMEM::DataArrayDouble::Divide(aaa,self->getArray());
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    case 4:
      {
        if(!self->getArray())
          throw INTERP_KERNEL::Exception(msg2);
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> aaa=ParaMEDMEM::DataArrayDouble::New(); aaa->useArray(&bb[0],false,ParaMEDMEM::CPP_DEALLOC,1,(int)bb.size());
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::DataArrayDouble> ret=ParaMEDMEM::DataArrayDouble::Divide(aaa,self->getArray());
        ParaMEDMEM::MEDCouplingAutoRefCountObjectPtr<ParaMEDMEM::MEDCouplingFieldDouble> ret2=self->clone(false);
        ret2->setArray(ret);
        return ret2.retn();
      }
    default:
      { throw INTERP_KERNEL::Exception(msg); }
    }
}
