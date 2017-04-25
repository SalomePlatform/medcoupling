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

template<class T>
typename MEDCoupling::Traits<T>::FieldType *fieldT_buildSubPart(const MEDCoupling::MEDCouplingFieldT<T> *self, PyObject *li)
{
  int sw;
  int singleVal;
  std::vector<int> multiVal;
  std::pair<int, std::pair<int,int> > slic;
  MEDCoupling::DataArrayInt *daIntTyypp=0;
  const MEDCoupling::MEDCouplingMesh *mesh=self->getMesh();
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::buildSubPart : field lies on a null mesh !");
  int nbc=mesh->getNumberOfCells();
  convertObjToPossibleCpp2(li,nbc,sw,singleVal,multiVal,slic,daIntTyypp);
  switch(sw)
    {
    case 1:
      {
        if(singleVal>=nbc)
          {
            std::ostringstream oss;
            oss << "Requesting for cell id " << singleVal << " having only " << nbc << " cells !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
        if(singleVal>=0)
          return self->buildSubPart(&singleVal,&singleVal+1);
        else
          {
            if(nbc+singleVal>0)
              {
                int tmp=nbc+singleVal;
                return self->buildSubPart(&tmp,&tmp+1);
              }
            else
              {
                std::ostringstream oss;
                oss << "Requesting for cell id " << singleVal << " having only " << nbc << " cells !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
          }
      }
    case 2:
      {
        return self->buildSubPart(&multiVal[0],&multiVal[0]+multiVal.size());
      }
    case 3:
      {
        return self->buildSubPartRange(slic.first,slic.second.first,slic.second.second);
      }
    case 4:
      {
        if(!daIntTyypp)
          throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::buildSubPart : null instance has been given in input !");
        daIntTyypp->checkAllocated();
        return self->buildSubPart(daIntTyypp->begin(),daIntTyypp->end());
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::buildSubPart : unrecognized type in input ! Possibilities are : int, list or tuple of int DataArrayInt instance !");
    }
}

template<class T>
typename MEDCoupling::Traits<T>::FieldType *fieldT__getitem__(const MEDCoupling::MEDCouplingFieldT<T> *self, PyObject *li)
{
  const char msg[]="MEDCouplingFieldDouble::__getitem__ : invalid call  Available API are : \n-myField[dataArrayInt]\n-myField[slice]\n-myField[pythonListOfCellIds]\n-myField[integer]\n-myField[dataArrayInt,1]\n-myField[slice,1]\n-myField[pythonListOfCellIds,1]\n-myField[integer,1]\n";
  if(PyTuple_Check(li))
    {
      Py_ssize_t sz=PyTuple_Size(li);
      if(sz!=2)
        throw INTERP_KERNEL::Exception(msg);
      PyObject *elt0=PyTuple_GetItem(li,0),*elt1=PyTuple_GetItem(li,1);
      int sw;
      int singleVal;
      std::vector<int> multiVal;
      std::pair<int, std::pair<int,int> > slic;
      MEDCoupling::DataArrayInt *daIntTyypp=0;
      if(!self->getArray())
        throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::__getitem__ : no array set on field to deduce number of components !");
      try
        { convertObjToPossibleCpp2(elt1,self->getArray()->getNumberOfComponents(),sw,singleVal,multiVal,slic,daIntTyypp); }
      catch(INTERP_KERNEL::Exception& e)
        { std::ostringstream oss; oss << "MEDCouplingFieldDouble::__getitem__ : invalid type in 2nd parameter (compo) !" << e.what(); throw INTERP_KERNEL::Exception(oss.str().c_str()); }
      typename MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::FieldType > ret0(fieldT_buildSubPart<T>(self,elt0));
      typename MEDCoupling::Traits<T>::ArrayType *ret0Arr=ret0->getArray();
      if(!ret0Arr)
        throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::__getitem__ : no array exists to apply restriction on component on it !");
      switch(sw)
        {
        case 1:
          {
            std::vector<int> v2(1,singleVal);
            MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > aarr(ret0Arr->keepSelectedComponents(v2));
            ret0->setArray(aarr);
            return ret0.retn();
          }
        case 2:
          {
            MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > aarr(ret0Arr->keepSelectedComponents(multiVal));
            ret0->setArray(aarr);
            return ret0.retn();
          }
        case 3:
          {
            int nbOfComp(MEDCoupling::DataArray::GetNumberOfItemGivenBESRelative(slic.first,slic.second.first,slic.second.second,"MEDCouplingFieldDouble::__getitem__ : invalid range in 2nd parameter (components) !"));
            std::vector<int> v2(nbOfComp);
            for(int i=0;i<nbOfComp;i++)
              v2[i]=slic.first+i*slic.second.second;
            MEDCoupling::MCAuto< typename MEDCoupling::Traits<T>::ArrayType > aarr(ret0Arr->keepSelectedComponents(v2));
            ret0->setArray(aarr);
            return ret0.retn();
          }
        default:
          throw INTERP_KERNEL::Exception(msg);
        }
    }
  else
    return fieldT_buildSubPart<T>(self,li);
}


#endif
