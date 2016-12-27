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

#include "MEDCouplingRemapper.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
#include "MEDCouplingMappedExtrudedMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingNormalizedUnstructuredMesh.txx"
#include "MEDCouplingNormalizedCartesianMesh.txx"

#include "Interpolation1D.txx"
#include "Interpolation2DCurve.hxx"
#include "Interpolation2D.txx"
#include "Interpolation3D.txx"
#include "Interpolation3DSurf.hxx"
#include "Interpolation2D1D.txx"
#include "Interpolation2D3D.txx"
#include "Interpolation3D1D.txx"
#include "InterpolationCU.txx"
#include "InterpolationCC.txx"

using namespace MEDCoupling;

MEDCouplingRemapper::MEDCouplingRemapper():_src_ft(0),_target_ft(0),_interp_matrix_pol(IK_ONLY_PREFERED),_nature_of_deno(NoNature),_time_deno_update(0)
{
}

MEDCouplingRemapper::~MEDCouplingRemapper()
{
  releaseData(false);
}

/*!
 * This method is the second step of the remapping process. The remapping process works in three phases :
 *
 * - Set remapping options appropriately
 * - The computation of remapping matrix
 * - Apply the matrix vector multiply to obtain the result of the remapping
 * 
 * This method performs the second step (computation of remapping matrix) which may be CPU-time consuming. This phase is also the most critical (where the most tricky algorithm) in the remapping process.
 * Strictly speaking to perform the computation of the remapping matrix the field templates source-side and target-side is required (which is the case of MEDCouplingRemapper::prepareEx).
 * So this method is less precise but a more user friendly way to compute a remapping matrix.
 *
 * \param [in] srcMesh the source mesh
 * \param [in] targetMesh the target mesh
 * \param [in] method A string obtained by aggregation of the spatial discretisation string representation of source field and target field. The string representation is those returned by MEDCouplingFieldDiscretization::getStringRepr.
 *             Example : "P0" is for cell discretization. "P1" is for node discretization. So "P0P1" for \a method parameter means from a source cell field (lying on \a srcMesh) to a target node field (lying on \a targetMesh).
 *
 * \sa MEDCouplingRemapper::prepareEx
 */
int MEDCouplingRemapper::prepare(const MEDCouplingMesh *srcMesh, const MEDCouplingMesh *targetMesh, const std::string& method)
{
  if(!srcMesh || !targetMesh)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepare : presence of NULL input pointer !");
  std::string srcMethod,targetMethod;
  INTERP_KERNEL::Interpolation<INTERP_KERNEL::Interpolation3D>::CheckAndSplitInterpolationMethod(method,srcMethod,targetMethod);
  MCAuto<MEDCouplingFieldTemplate> src=MEDCouplingFieldTemplate::New(MEDCouplingFieldDiscretization::GetTypeOfFieldFromStringRepr(srcMethod));
  src->setMesh(srcMesh);
  MCAuto<MEDCouplingFieldTemplate> target=MEDCouplingFieldTemplate::New(MEDCouplingFieldDiscretization::GetTypeOfFieldFromStringRepr(targetMethod));
  target->setMesh(targetMesh);
  return prepareEx(src,target);
}

/*!
 * This method is the generalization of MEDCouplingRemapper::prepare. Indeed, MEDCouplingFieldTemplate instances gives all required information to compute the remapping matrix.
 * This method must be used instead of MEDCouplingRemapper::prepare if Gauss point to Gauss point must be applied.
 *
 * \param [in] src is the field template source side.
 * \param [in] target is the field template target side.
 *
 * \sa MEDCouplingRemapper::prepare
 */
int MEDCouplingRemapper::prepareEx(const MEDCouplingFieldTemplate *src, const MEDCouplingFieldTemplate *target)
{
  if(!src || !target)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareEx : presence of NULL input pointer !");
  if(!src->getMesh() || !target->getMesh())
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareEx : presence of NULL mesh pointer in given field template !");
  releaseData(true);
  _src_ft=const_cast<MEDCouplingFieldTemplate *>(src); _src_ft->incrRef();
  _target_ft=const_cast<MEDCouplingFieldTemplate *>(target); _target_ft->incrRef();
  if(isInterpKernelOnlyOrNotOnly())
    return prepareInterpKernelOnly();
  else
    return prepareNotInterpKernelOnly();
}

int MEDCouplingRemapper::prepareInterpKernelOnly()
{
  int meshInterpType=((int)_src_ft->getMesh()->getType()*16)+(int)_target_ft->getMesh()->getType();
  // *** Remember:
//  typedef enum
//    {
//      UNSTRUCTURED = 5,
//      CARTESIAN = 7,
//      EXTRUDED = 8,
//      CURVE_LINEAR = 9,
//      SINGLE_STATIC_GEO_TYPE_UNSTRUCTURED = 10,
//      SINGLE_DYNAMIC_GEO_TYPE_UNSTRUCTURED = 11,
//      IMAGE_GRID = 12
//    } MEDCouplingMeshType;

  switch(meshInterpType)
  {
    case 90:   // UNSTRUCTURED - SINGLE_STATIC_GEO_TYPE_UNSTRUCTURED
    case 91:   // UNSTRUCTURED - SINGLE_DYNAMIC_GEO_TYPE_UNSTRUCTURED
    case 165:  // SINGLE_STATIC_GEO_TYPE_UNSTRUCTURED - UNSTRUCTURED
    case 181:  // SINGLE_DYNAMIC_GEO_TYPE_UNSTRUCTURED - UNSTRUCTURED
    case 170:  // SINGLE_STATIC_GEO_TYPE_UNSTRUCTURED - SINGLE_STATIC_GEO_TYPE_UNSTRUCTURED
    case 171:  // SINGLE_STATIC_GEO_TYPE_UNSTRUCTURED - SINGLE_DYNAMIC_GEO_TYPE_UNSTRUCTURED
    case 186:  // SINGLE_DYNAMIC_GEO_TYPE_UNSTRUCTURED - SINGLE_STATIC_GEO_TYPE_UNSTRUCTURED
    case 187:  // SINGLE_DYNAMIC_GEO_TYPE_UNSTRUCTURED - SINGLE_DYNAMIC_GEO_TYPE_UNSTRUCTURED
    case 85:   // UNSTRUCTURED - UNSTRUCTURED
      return prepareInterpKernelOnlyUU();
    case 167:  // SINGLE_STATIC_GEO_TYPE_UNSTRUCTURED - CARTESIAN
    case 183:  // SINGLE_DYNAMIC_GEO_TYPE_UNSTRUCTURED - CARTESIAN
    case 87:   // UNSTRUCTURED - CARTESIAN
      return prepareInterpKernelOnlyUC();
    case 122:  // CARTESIAN - SINGLE_STATIC_GEO_TYPE_UNSTRUCTURED
    case 123:  // CARTESIAN - SINGLE_DYNAMIC_GEO_TYPE_UNSTRUCTURED
    case 117:  // CARTESIAN - UNSTRUCTURED
      return prepareInterpKernelOnlyCU();
    case 119:  // CARTESIAN - CARTESIAN
      return prepareInterpKernelOnlyCC();
    case 136:  // EXTRUDED - EXTRUDED
      return prepareInterpKernelOnlyEE();
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareInterpKernelOnly : Not managed type of meshes ! Dealt meshes type are : Unstructured<->Unstructured, Unstructured<->Cartesian, Cartesian<->Cartesian, Extruded<->Extruded !");
  }
}

int MEDCouplingRemapper::prepareNotInterpKernelOnly()
{
  std::string srcm,trgm,method;
  method=checkAndGiveInterpolationMethodStr(srcm,trgm);
  switch(CheckInterpolationMethodManageableByNotOnlyInterpKernel(method))
  {
    case 0:
      return prepareNotInterpKernelOnlyGaussGauss();
    default:
      {
        std::ostringstream oss; oss << "MEDCouplingRemapper::prepareNotInterpKernelOnly : INTERNAL ERROR ! the method \"" << method << "\" declared as managed bu not implemented !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  }
}

/*!
 * This method performs the operation source to target using matrix computed in MEDCoupling::MEDCouplingRemapper::prepare method.
 * If meshes of \b srcField and \b targetField do not match exactly those given into \ref MEDCoupling::MEDCouplingRemapper::prepare "prepare method" an exception will be thrown.
 * 
 * \param [in] srcField is the source field from which the interpolation will be done. The mesh into \b srcField should be the same than those specified on MEDCoupling::MEDCouplingRemapper::prepare.
 * \param [in/out] targetField the destination field with the allocated array in which all tuples will be overwritten.
 * \param [in] dftValue is the value that will be assigned in the targetField to each entity of target mesh (entity depending on the method selected on prepare invocation) that is not intercepted by any entity of source mesh.
 *             For example in "P0P0" case (cell-cell) if a cell in target mesh is not overlapped by any source cell the \a dftValue value will be attached on that cell in the returned \a targetField. In some cases a target
 *             cell not intercepted by any source cell is a bug so in this case it is advised to set a huge value (1e300 for example) to \a dftValue to quickly point to the problem. But for users doing parallelism a target cell can
 *             be intercepted by a source cell on a different process. In this case 0. assigned to \a dftValue is more appropriate.
 *
 * \sa transferField
 */
void MEDCouplingRemapper::transfer(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField, double dftValue)
{
  if(!srcField || !targetField)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::transfer : input field must be both not NULL !");
  transferUnderground(srcField,targetField,true,dftValue);
}

/*!
 * This method is equivalent to MEDCoupling::MEDCouplingRemapper::transfer except that here \b targetField is a in/out parameter.
 * If an entity (cell for example) in targetField is not fetched by any entity (cell for example) of \b srcField, the value in targetField is
 * let unchanged.
 * This method requires that \b targetField was fully defined and allocated. If the array is not allocated an exception will be thrown.
 * 
 * \param [in] srcField is the source field from which the interpolation will be done. The mesh into \b srcField should be the same than those specified on MEDCoupling::MEDCouplingRemapper::prepare.
 * \param [in,out] targetField the destination field with the allocated array in which only tuples whose entities are fetched by interpolation will be overwritten only.
 */
void MEDCouplingRemapper::partialTransfer(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField)
{
  if(!srcField || !targetField)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::partialTransfer : input field must be both not NULL !");
  transferUnderground(srcField,targetField,false,std::numeric_limits<double>::max());
}

void MEDCouplingRemapper::reverseTransfer(MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *targetField, double dftValue)
{
  if(!srcField || !targetField)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::reverseTransfer : input fields must be both not NULL !");
  checkPrepare();
  targetField->checkConsistencyLight();
  if(_src_ft->getDiscretization()->getStringRepr()!=srcField->getDiscretization()->getStringRepr())
    throw INTERP_KERNEL::Exception("Incoherency with prepare call for source field");
  if(_target_ft->getDiscretization()->getStringRepr()!=targetField->getDiscretization()->getStringRepr())
    throw INTERP_KERNEL::Exception("Incoherency with prepare call for target field");
  if(srcField->getNature()!=targetField->getNature())
    throw INTERP_KERNEL::Exception("Natures of fields mismatch !");
  if(targetField->getNumberOfTuplesExpected()!=_target_ft->getNumberOfTuplesExpected())
    {
      std::ostringstream oss;
      oss << "MEDCouplingRemapper::reverseTransfer : in given source field the number of tuples required is " << _target_ft->getNumberOfTuplesExpected() << " (on prepare) and number of tuples in given target field is " << targetField->getNumberOfTuplesExpected();
      oss << " ! It appears that the target support is not the same between the prepare and the transfer !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  DataArrayDouble *array(srcField->getArray());
  int trgNbOfCompo=targetField->getNumberOfComponents();
  if(array)
    {
      srcField->checkConsistencyLight();
      if(trgNbOfCompo!=srcField->getNumberOfTuplesExpected())
        throw INTERP_KERNEL::Exception("Number of components mismatch !");
    }
  else
    {
      MCAuto<DataArrayDouble > tmp(DataArrayDouble::New());
      tmp->alloc(srcField->getNumberOfTuplesExpected(),trgNbOfCompo);
      srcField->setArray(tmp);
    }
  computeDeno(srcField->getNature(),srcField,targetField);
  double *resPointer(srcField->getArray()->getPointer());
  const double *inputPointer=targetField->getArray()->getConstPointer();
  computeReverseProduct(inputPointer,trgNbOfCompo,dftValue,resPointer);
}

/*!
 * This method performs the operation source to target using matrix computed in MEDCoupling::MEDCouplingRemapper::prepare method.
 * If mesh of \b srcField does not match exactly those given into \ref MEDCoupling::MEDCouplingRemapper::prepare "prepare method" an exception will be thrown.
 * 
 * \param [in] srcField is the source field from which the interpolation will be done. The mesh into \b srcField should be the same than those specified on MEDCoupling::MEDCouplingRemapper::prepare.
 * \param [in] dftValue is the value that will be assigned in the targetField to each entity of target mesh (entity depending on the method selected on prepare invocation) that is not intercepted by any entity of source mesh.
 *             For example in "P0P0" case (cell-cell) if a cell in target mesh is not overlapped by any source cell the \a dftValue value will be attached on that cell in the returned \a targetField. In some cases a target
 *             cell not intercepted by any source cell is a bug so in this case it is advised to set a huge value (1e300 for example) to \a dftValue to quickly point to the problem. But for users doing parallelism a target cell can
 *             be intercepted by a source cell on a different process. In this case 0. assigned to \a dftValue is more appropriate.
 * \return destination field to be deallocated by the caller.
 *
 * \sa transfer
 */
MEDCouplingFieldDouble *MEDCouplingRemapper::transferField(const MEDCouplingFieldDouble *srcField, double dftValue)
{
  checkPrepare();
  if(!srcField)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::transferField : input srcField is NULL !");
  srcField->checkConsistencyLight();
  if(_src_ft->getDiscretization()->getStringRepr()!=srcField->getDiscretization()->getStringRepr())
    throw INTERP_KERNEL::Exception("Incoherency with prepare call for source field");
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(*_target_ft,srcField->getTimeDiscretization());
  ret->setNature(srcField->getNature());
  transfer(srcField,ret,dftValue);
  ret->copyAllTinyAttrFrom(srcField);//perform copy of tiny strings after and not before transfer because the array will be created on transfer
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingRemapper::reverseTransferField(const MEDCouplingFieldDouble *targetField, double dftValue)
{
  if(!targetField)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::transferField : input targetField is NULL !");
  targetField->checkConsistencyLight();
  checkPrepare();
  if(_target_ft->getDiscretization()->getStringRepr()!=targetField->getDiscretization()->getStringRepr())
    throw INTERP_KERNEL::Exception("Incoherency with prepare call for target field");
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(*_src_ft,targetField->getTimeDiscretization());
  ret->setNature(targetField->getNature());
  reverseTransfer(ret,targetField,dftValue);
  ret->copyAllTinyAttrFrom(targetField);//perform copy of tiny strings after and not before reverseTransfer because the array will be created on reverseTransfer
  return ret;
}

/*!
 * This method does nothing more than inherited INTERP_KERNEL::InterpolationOptions::setOptionInt method. This method
 * is here only for automatic CORBA generators.
 */
bool MEDCouplingRemapper::setOptionInt(const std::string& key, int value)
{
  return INTERP_KERNEL::InterpolationOptions::setOptionInt(key,value);
}

/*!
 * This method does nothing more than inherited INTERP_KERNEL::InterpolationOptions::setOptionInt method. This method
 * is here only for automatic CORBA generators.
 */
bool MEDCouplingRemapper::setOptionDouble(const std::string& key, double value)
{
  return INTERP_KERNEL::InterpolationOptions::setOptionDouble(key,value);
}

/*!
 * This method does nothing more than inherited INTERP_KERNEL::InterpolationOptions::setOptionInt method. This method
 * is here only for automatic CORBA generators.
 */
bool MEDCouplingRemapper::setOptionString(const std::string& key, const std::string& value)
{
  return INTERP_KERNEL::InterpolationOptions::setOptionString(key,value);
}

/*!
 * This method returns the interpolation matrix policy. This policy specifies which interpolation matrix method to keep or prefered.
 * If interpolation matrix policy is :
 *
 * - set to IK_ONLY_PREFERED (0) (the default) : the INTERP_KERNEL only method is prefered. That is to say, if it is possible to treat the case
 *   regarding spatial discretization of source and target with INTERP_KERNEL only method, INTERP_KERNEL only method will be performed.
 *   If not, the \b not only INTERP_KERNEL method will be attempt.
 * 
 * - set to NOT_IK_ONLY_PREFERED (1) : the \b NOT only INTERP_KERNEL method is prefered. That is to say, if it is possible to treat the case
 *   regarding spatial discretization of source and target with \b NOT only INTERP_KERNEL method, \b NOT only INTERP_KERNEL method, will be performed.
 *   If not, the INTERP_KERNEL only method will be attempt.
 * 
 * - IK_ONLY_FORCED (2) : Only INTERP_KERNEL only method will be launched.
 *
 * - NOT_IK_ONLY_FORCED (3) : Only \b NOT INTERP_KERNEL only method will be launched.
 * 
 * \sa MEDCouplingRemapper::setInterpolationMatrixPolicy
 */
int MEDCouplingRemapper::getInterpolationMatrixPolicy() const
{
  return _interp_matrix_pol;
}

/*!
 * This method sets a new interpolation matrix policy. The default one is IK_PREFERED (0). The input is of type \c int to be dealt by standard Salome
 * CORBA component generators. This method throws an INTERP_KERNEL::Exception if a the input integer is not in the available possibilities, that is to say not in
 * [0 (IK_PREFERED) , 1 (NOT_IK_PREFERED), 2 (IK_ONLY_FORCED), 3 (NOT_IK_ONLY_FORCED)].
 *
 * If interpolation matrix policy is :
 *
 * - set to IK_ONLY_PREFERED (0) (the default) : the INTERP_KERNEL only method is prefered. That is to say, if it is possible to treat the case
 *   regarding spatial discretization of source and target with INTERP_KERNEL only method, INTERP_KERNEL only method will be performed.
 *   If not, the \b not only INTERP_KERNEL method will be attempt.
 * 
 * - set to NOT_IK_ONLY_PREFERED (1) : the \b NOT only INTERP_KERNEL method is prefered. That is to say, if it is possible to treat the case
 *   regarding spatial discretization of source and target with \b NOT only INTERP_KERNEL method, \b NOT only INTERP_KERNEL method, will be performed.
 *   If not, the INTERP_KERNEL only method will be attempt.
 * 
 * - IK_ONLY_FORCED (2) : Only INTERP_KERNEL only method will be launched.
 *
 * - NOT_IK_ONLY_FORCED (3) : Only \b NOT INTERP_KERNEL only method will be launched.
 * 
 * \input newInterpMatPol the new interpolation matrix method policy. This parameter is of type \c int and not of type \c MEDCoupling::InterpolationMatrixPolicy
 *                        for automatic generation of CORBA component.
 * 
 * \sa MEDCouplingRemapper::getInterpolationMatrixPolicy
 */
void MEDCouplingRemapper::setInterpolationMatrixPolicy(int newInterpMatPol)
{
  switch(newInterpMatPol)
  {
    case 0:
      _interp_matrix_pol=IK_ONLY_PREFERED;
      break;
    case 1:
      _interp_matrix_pol=NOT_IK_ONLY_PREFERED;
      break;
    case 2:
      _interp_matrix_pol=IK_ONLY_FORCED;
      break;
    case 3:
      _interp_matrix_pol=NOT_IK_ONLY_FORCED;
      break;
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingRemapper::setInterpolationMatrixPolicy : invalid input integer value ! Should be in [0 (IK_PREFERED) , 1 (NOT_IK_PREFERED), 2 (IK_ONLY_FORCED), 3 (NOT_IK_ONLY_FORCED)] ! For information, the default is IK_PREFERED=0 !");
  }
}

int MEDCouplingRemapper::prepareInterpKernelOnlyUU()
{
  const MEDCouplingPointSet *src_mesh=static_cast<const MEDCouplingPointSet *>(_src_ft->getMesh());
  const MEDCouplingPointSet *target_mesh=static_cast<const MEDCouplingPointSet *>(_target_ft->getMesh());
  std::string srcMeth,trgMeth;
  std::string method(checkAndGiveInterpolationMethodStr(srcMeth,trgMeth));
  const int srcMeshDim=src_mesh->getMeshDimension();
  int srcSpaceDim=-1;
  if(srcMeshDim!=-1)
    srcSpaceDim=src_mesh->getSpaceDimension();
  const int trgMeshDim=target_mesh->getMeshDimension();
  int trgSpaceDim=-1;
  if(trgMeshDim!=-1)
    trgSpaceDim=target_mesh->getSpaceDimension();
  if(trgSpaceDim!=srcSpaceDim)
    if(trgSpaceDim!=-1 && srcSpaceDim!=-1)
      throw INTERP_KERNEL::Exception("Incoherent space dimension detected between target and source.");
  int nbCols;
  if(srcMeshDim==1 && trgMeshDim==1 && srcSpaceDim==1)
    {
      MEDCouplingNormalizedUnstructuredMesh<1,1> source_mesh_wrapper(src_mesh);
      MEDCouplingNormalizedUnstructuredMesh<1,1> target_mesh_wrapper(target_mesh);
      INTERP_KERNEL::Interpolation1D interpolation(*this);
      nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
    }
  else if(srcMeshDim==1 && trgMeshDim==1 && srcSpaceDim==2)
    {
      MEDCouplingNormalizedUnstructuredMesh<2,1> source_mesh_wrapper(src_mesh);
      MEDCouplingNormalizedUnstructuredMesh<2,1> target_mesh_wrapper(target_mesh);
      INTERP_KERNEL::Interpolation2DCurve interpolation(*this);
      nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
    }
  else if(srcMeshDim==2 && trgMeshDim==2 && srcSpaceDim==2)
    {
      MEDCouplingNormalizedUnstructuredMesh<2,2> source_mesh_wrapper(src_mesh);
      MEDCouplingNormalizedUnstructuredMesh<2,2> target_mesh_wrapper(target_mesh);
      INTERP_KERNEL::Interpolation2D interpolation(*this);
      nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
    }
  else if(srcMeshDim==3 && trgMeshDim==3 && srcSpaceDim==3)
    {
      MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(src_mesh);
      MEDCouplingNormalizedUnstructuredMesh<3,3> target_mesh_wrapper(target_mesh);
      INTERP_KERNEL::Interpolation3D interpolation(*this);
      nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
    }
  else if(srcMeshDim==2 && trgMeshDim==2 && srcSpaceDim==3)
    {
      MEDCouplingNormalizedUnstructuredMesh<3,2> source_mesh_wrapper(src_mesh);
      MEDCouplingNormalizedUnstructuredMesh<3,2> target_mesh_wrapper(target_mesh);
      INTERP_KERNEL::Interpolation3DSurf interpolation(*this);
      nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
    }
  else if(srcMeshDim==3 && trgMeshDim==1 && srcSpaceDim==3)
    {
      if(getIntersectionType()!=INTERP_KERNEL::PointLocator)
        throw INTERP_KERNEL::Exception("Invalid interpolation requested between 3D and 1D ! Select PointLocator as intersection type !");
      MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(src_mesh);
      MEDCouplingNormalizedUnstructuredMesh<3,3> target_mesh_wrapper(target_mesh);
      INTERP_KERNEL::Interpolation3D1D interpolation(*this);
      nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
    }
  else if(srcMeshDim==1 && trgMeshDim==3 && srcSpaceDim==3)
    {
      if(getIntersectionType()!=INTERP_KERNEL::PointLocator)
        throw INTERP_KERNEL::Exception("Invalid interpolation requested between 3D and 1D ! Select PointLocator as intersection type !");
      MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(src_mesh);
      MEDCouplingNormalizedUnstructuredMesh<3,3> target_mesh_wrapper(target_mesh);
      INTERP_KERNEL::Interpolation3D1D interpolation(*this);
      std::vector<std::map<int,double> > matrixTmp;
      std::string revMethod(BuildMethodFrom(trgMeth,srcMeth));
      nbCols=interpolation.interpolateMeshes(target_mesh_wrapper,source_mesh_wrapper,matrixTmp,revMethod);
      ReverseMatrix(matrixTmp,nbCols,_matrix);
      nbCols=matrixTmp.size();
    }
  else if(srcMeshDim==2 && trgMeshDim==1 && srcSpaceDim==2)
    {
      if(getIntersectionType()==INTERP_KERNEL::PointLocator)
        {
          MEDCouplingNormalizedUnstructuredMesh<2,2> source_mesh_wrapper(src_mesh);
          MEDCouplingNormalizedUnstructuredMesh<2,2> target_mesh_wrapper(target_mesh);
          INTERP_KERNEL::Interpolation2D interpolation(*this);
          nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
        }
      else
        {
          MEDCouplingNormalizedUnstructuredMesh<2,2> source_mesh_wrapper(src_mesh);
          MEDCouplingNormalizedUnstructuredMesh<2,2> target_mesh_wrapper(target_mesh);
          INTERP_KERNEL::Interpolation2D1D interpolation(*this);
          std::vector<std::map<int,double> > matrixTmp;
          std::string revMethod(BuildMethodFrom(trgMeth,srcMeth));
          nbCols=interpolation.interpolateMeshes(target_mesh_wrapper,source_mesh_wrapper,matrixTmp,revMethod);
          ReverseMatrix(matrixTmp,nbCols,_matrix);
          nbCols=matrixTmp.size();
          INTERP_KERNEL::Interpolation2D1D::DuplicateFacesType duplicateFaces=interpolation.retrieveDuplicateFaces();
          if(!duplicateFaces.empty())
            {
              std::ostringstream oss; oss << "An unexpected situation happened ! For the following 1D Cells are part of edges shared by 2D cells :\n";
              for(std::map<int,std::set<int> >::const_iterator it=duplicateFaces.begin();it!=duplicateFaces.end();it++)
                {
                  oss << "1D Cell #" << (*it).first << " is part of common edge of following 2D cells ids : ";
                  std::copy((*it).second.begin(),(*it).second.end(),std::ostream_iterator<int>(oss," "));
                  oss << std::endl;
                }
            }
        }
    }
  else if(srcMeshDim==1 && trgMeshDim==2 && srcSpaceDim==2)
    {
      if(getIntersectionType()==INTERP_KERNEL::PointLocator)
        {
          MEDCouplingNormalizedUnstructuredMesh<2,2> source_mesh_wrapper(src_mesh);
          MEDCouplingNormalizedUnstructuredMesh<2,2> target_mesh_wrapper(target_mesh);
          INTERP_KERNEL::Interpolation2D interpolation(*this);
          std::vector<std::map<int,double> > matrixTmp;
          std::string revMethod(BuildMethodFrom(trgMeth,srcMeth));
          nbCols=interpolation.interpolateMeshes(target_mesh_wrapper,source_mesh_wrapper,matrixTmp,revMethod);
          ReverseMatrix(matrixTmp,nbCols,_matrix);
          nbCols=matrixTmp.size();
        }
      else
        {
          MEDCouplingNormalizedUnstructuredMesh<2,2> source_mesh_wrapper(src_mesh);
          MEDCouplingNormalizedUnstructuredMesh<2,2> target_mesh_wrapper(target_mesh);
          INTERP_KERNEL::Interpolation2D1D interpolation(*this);
          nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
          INTERP_KERNEL::Interpolation2D1D::DuplicateFacesType duplicateFaces=interpolation.retrieveDuplicateFaces();
          if(!duplicateFaces.empty())
            {
              std::ostringstream oss; oss << "An unexpected situation happend ! For the following 1D Cells are part of edges shared by 2D cells :\n";
              for(std::map<int,std::set<int> >::const_iterator it=duplicateFaces.begin();it!=duplicateFaces.end();it++)
                {
                  oss << "1D Cell #" << (*it).first << " is part of common edge of following 2D cells ids : ";
                  std::copy((*it).second.begin(),(*it).second.end(),std::ostream_iterator<int>(oss," "));
                  oss << std::endl;
                }
            }
        }
    }
  else if(srcMeshDim==2 && trgMeshDim==3 && srcSpaceDim==3)
    {
      MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(src_mesh);
      MEDCouplingNormalizedUnstructuredMesh<3,3> target_mesh_wrapper(target_mesh);
      INTERP_KERNEL::Interpolation2D3D interpolation(*this);
      nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
      INTERP_KERNEL::Interpolation2D3D::DuplicateFacesType duplicateFaces=interpolation.retrieveDuplicateFaces();
      if(!duplicateFaces.empty())
        {
          std::ostringstream oss; oss << "An unexpected situation happend ! For the following 2D Cells are part of edges shared by 3D cells :\n";
          for(std::map<int,std::set<int> >::const_iterator it=duplicateFaces.begin();it!=duplicateFaces.end();it++)
            {
              oss << "2D Cell #" << (*it).first << " is part of common face of following 3D cells ids : ";
              std::copy((*it).second.begin(),(*it).second.end(),std::ostream_iterator<int>(oss," "));
              oss << std::endl;
            }
        }
    }
  else if(srcMeshDim==3 && trgMeshDim==2 && srcSpaceDim==3)
    {
      if(getIntersectionType()==INTERP_KERNEL::PointLocator)
        {
          MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(src_mesh);
          MEDCouplingNormalizedUnstructuredMesh<3,3> target_mesh_wrapper(target_mesh);
          INTERP_KERNEL::Interpolation3D interpolation(*this);
          nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
        }
      else
        {
          MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(src_mesh);
          MEDCouplingNormalizedUnstructuredMesh<3,3> target_mesh_wrapper(target_mesh);
          INTERP_KERNEL::Interpolation2D3D interpolation(*this);
          std::vector<std::map<int,double> > matrixTmp;
          std::string revMethod(BuildMethodFrom(trgMeth,srcMeth));
          nbCols=interpolation.interpolateMeshes(target_mesh_wrapper,source_mesh_wrapper,matrixTmp,revMethod);
          ReverseMatrix(matrixTmp,nbCols,_matrix);
          nbCols=matrixTmp.size();
          INTERP_KERNEL::Interpolation2D3D::DuplicateFacesType duplicateFaces=interpolation.retrieveDuplicateFaces();
          if(!duplicateFaces.empty())
            {
              std::ostringstream oss; oss << "An unexpected situation happend ! For the following 2D Cells are part of edges shared by 3D cells :\n";
              for(std::map<int,std::set<int> >::const_iterator it=duplicateFaces.begin();it!=duplicateFaces.end();it++)
                {
                  oss << "2D Cell #" << (*it).first << " is part of common face of following 3D cells ids : ";
                  std::copy((*it).second.begin(),(*it).second.end(),std::ostream_iterator<int>(oss," "));
                  oss << std::endl;
                }
            }
        }
    }
  else if(trgMeshDim==-1)
    {
      if(srcMeshDim==2 && srcSpaceDim==2)
        {
          MEDCouplingNormalizedUnstructuredMesh<2,2> source_mesh_wrapper(src_mesh);
          INTERP_KERNEL::Interpolation2D interpolation(*this);
          nbCols=interpolation.toIntegralUniform(source_mesh_wrapper,_matrix,srcMeth);
        }
      else if(srcMeshDim==3 && srcSpaceDim==3)
        {
          MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(src_mesh);
          INTERP_KERNEL::Interpolation3D interpolation(*this);
          nbCols=interpolation.toIntegralUniform(source_mesh_wrapper,_matrix,srcMeth);
        }
      else if(srcMeshDim==2 && srcSpaceDim==3)
        {
          MEDCouplingNormalizedUnstructuredMesh<3,2> source_mesh_wrapper(src_mesh);
          INTERP_KERNEL::Interpolation3DSurf interpolation(*this);
          nbCols=interpolation.toIntegralUniform(source_mesh_wrapper,_matrix,srcMeth);
        }
      else
        throw INTERP_KERNEL::Exception("No interpolation available for the given mesh and space dimension of source mesh to -1D targetMesh");
    }
  else if(srcMeshDim==-1)
    {
      if(trgMeshDim==2 && trgSpaceDim==2)
        {
          MEDCouplingNormalizedUnstructuredMesh<2,2> source_mesh_wrapper(target_mesh);
          INTERP_KERNEL::Interpolation2D interpolation(*this);
          nbCols=interpolation.fromIntegralUniform(source_mesh_wrapper,_matrix,trgMeth);
        }
      else if(trgMeshDim==3 && trgSpaceDim==3)
        {
          MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(target_mesh);
          INTERP_KERNEL::Interpolation3D interpolation(*this);
          nbCols=interpolation.fromIntegralUniform(source_mesh_wrapper,_matrix,trgMeth);
        }
      else if(trgMeshDim==2 && trgSpaceDim==3)
        {
          MEDCouplingNormalizedUnstructuredMesh<3,2> source_mesh_wrapper(target_mesh);
          INTERP_KERNEL::Interpolation3DSurf interpolation(*this);
          nbCols=interpolation.fromIntegralUniform(source_mesh_wrapper,_matrix,trgMeth);
        }
      else
        throw INTERP_KERNEL::Exception("No interpolation available for the given mesh and space dimension of source mesh from -1D sourceMesh");
    }
  else
    throw INTERP_KERNEL::Exception("No interpolation available for the given mesh and space dimension");
  _deno_multiply.clear();
  _deno_multiply.resize(_matrix.size());
  _deno_reverse_multiply.clear();
  _deno_reverse_multiply.resize(nbCols);
  declareAsNew();
  return 1;
}

int MEDCouplingRemapper::prepareInterpKernelOnlyEE()
{
  std::string srcMeth,trgMeth;
  std::string methC=checkAndGiveInterpolationMethodStr(srcMeth,trgMeth);
  const MEDCouplingMappedExtrudedMesh *src_mesh=static_cast<const MEDCouplingMappedExtrudedMesh *>(_src_ft->getMesh());
  const MEDCouplingMappedExtrudedMesh *target_mesh=static_cast<const MEDCouplingMappedExtrudedMesh *>(_target_ft->getMesh());
  if(methC!="P0P0")
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareInterpKernelOnlyEE : Only P0P0 method implemented for Extruded/Extruded meshes !");
  MCAuto<MEDCouplingUMesh> src2D(src_mesh->getMesh2D()->clone(false)); src2D->changeSpaceDimension(2,0.);
  MCAuto<MEDCouplingUMesh> trg2D(target_mesh->getMesh2D()->clone(false)); trg2D->changeSpaceDimension(2,0.);
  MEDCouplingNormalizedUnstructuredMesh<2,2> source_mesh_wrapper(src2D);
  MEDCouplingNormalizedUnstructuredMesh<2,2> target_mesh_wrapper(trg2D);
  INTERP_KERNEL::Interpolation2D interpolation2D(*this);
  std::vector<std::map<int,double> > matrix2D;
  int nbCols2D=interpolation2D.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,matrix2D,methC);
  MEDCouplingUMesh *s1D,*t1D;
  double v[3];
  MEDCouplingMappedExtrudedMesh::Project1DMeshes(src_mesh->getMesh1D(),target_mesh->getMesh1D(),getPrecision(),s1D,t1D,v);
  MEDCouplingNormalizedUnstructuredMesh<1,1> s1DWrapper(s1D);
  MEDCouplingNormalizedUnstructuredMesh<1,1> t1DWrapper(t1D);
  std::vector<std::map<int,double> > matrix1D;
  INTERP_KERNEL::Interpolation1D interpolation1D(*this);
  if(interpolation1D.getIntersectionType()==INTERP_KERNEL::Geometric2D)// For intersection type of 1D, Geometric2D do not deal with it ! -> make interpolation1D not inherite from this
    interpolation1D.setIntersectionType(INTERP_KERNEL::Triangulation);//
  int nbCols1D=interpolation1D.interpolateMeshes(s1DWrapper,t1DWrapper,matrix1D,methC);
  s1D->decrRef();
  t1D->decrRef();
  buildFinalInterpolationMatrixByConvolution(matrix1D,matrix2D,src_mesh->getMesh3DIds()->getConstPointer(),nbCols2D,nbCols1D,
                                             target_mesh->getMesh3DIds()->getConstPointer());
  //
  _deno_multiply.clear();
  _deno_multiply.resize(_matrix.size());
  _deno_reverse_multiply.clear();
  _deno_reverse_multiply.resize(nbCols2D*nbCols1D);
  declareAsNew();
  return 1;
}

int MEDCouplingRemapper::prepareInterpKernelOnlyUC()
{
  std::string srcMeth,trgMeth;
  std::string methodCpp=checkAndGiveInterpolationMethodStr(srcMeth,trgMeth);
  if(methodCpp!="P0P0")
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareInterpKernelOnlyUC : only P0P0 interpolation supported for the moment !");
  const MEDCouplingUMesh *src_mesh=static_cast<const MEDCouplingUMesh *>(_src_ft->getMesh());
  const MEDCouplingCMesh *target_mesh=static_cast<const MEDCouplingCMesh *>(_target_ft->getMesh());
  const int srcMeshDim=src_mesh->getMeshDimension();
  const int srcSpceDim=src_mesh->getSpaceDimension();
  const int trgMeshDim=target_mesh->getMeshDimension();
  if(srcMeshDim!=srcSpceDim || srcMeshDim!=trgMeshDim)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareInterpKernelOnlyUC : space dim of src unstructured should be equal to mesh dim of src unstructured and should be equal also equal to trg cartesian dimension !");
  std::vector<std::map<int,double> > res;
  switch(srcMeshDim)
  {
    case 1:
      {
        MEDCouplingNormalizedCartesianMesh<1> targetWrapper(target_mesh);
        MEDCouplingNormalizedUnstructuredMesh<1,1> sourceWrapper(src_mesh);
        INTERP_KERNEL::InterpolationCU myInterpolator(*this);
        myInterpolator.interpolateMeshes(targetWrapper,sourceWrapper,res,"P0P0");
        break;
      }
    case 2:
      {
        MEDCouplingNormalizedCartesianMesh<2> targetWrapper(target_mesh);
        MEDCouplingNormalizedUnstructuredMesh<2,2> sourceWrapper(src_mesh);
        INTERP_KERNEL::InterpolationCU myInterpolator(*this);
        myInterpolator.interpolateMeshes(targetWrapper,sourceWrapper,res,"P0P0");
        break;
      }
    case 3:
      {
        MEDCouplingNormalizedCartesianMesh<3> targetWrapper(target_mesh);
        MEDCouplingNormalizedUnstructuredMesh<3,3> sourceWrapper(src_mesh);
        INTERP_KERNEL::InterpolationCU myInterpolator(*this);
        myInterpolator.interpolateMeshes(targetWrapper,sourceWrapper,res,"P0P0");
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareInterpKernelOnlyUC : only dimension 1 2 or 3 supported !");
  }
  ReverseMatrix(res,target_mesh->getNumberOfCells(),_matrix);
  nullifiedTinyCoeffInCrudeMatrixAbs(0.);
  //
  _deno_multiply.clear();
  _deno_multiply.resize(_matrix.size());
  _deno_reverse_multiply.clear();
  _deno_reverse_multiply.resize(src_mesh->getNumberOfCells());
  declareAsNew();
  return 1;
}

int MEDCouplingRemapper::prepareInterpKernelOnlyCU()
{
  std::string srcMeth,trgMeth;
  std::string methodCpp=checkAndGiveInterpolationMethodStr(srcMeth,trgMeth);
  if(methodCpp!="P0P0")
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareInterpKernelOnlyCU : only P0P0 interpolation supported for the moment !");
  const MEDCouplingCMesh *src_mesh=static_cast<const MEDCouplingCMesh *>(_src_ft->getMesh());
  const MEDCouplingUMesh *target_mesh=static_cast<const MEDCouplingUMesh *>(_target_ft->getMesh());
  const int srcMeshDim=src_mesh->getMeshDimension();
  const int trgMeshDim=target_mesh->getMeshDimension();
  const int trgSpceDim=target_mesh->getSpaceDimension();
  if(trgMeshDim!=trgSpceDim || trgMeshDim!=srcMeshDim)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareInterpKernelOnlyCU : space dim of target unstructured should be equal to mesh dim of target unstructured and should be equal also equal to source cartesian dimension !");
  switch(srcMeshDim)
  {
    case 1:
      {
        MEDCouplingNormalizedCartesianMesh<1> sourceWrapper(src_mesh);
        MEDCouplingNormalizedUnstructuredMesh<1,1> targetWrapper(target_mesh);
        INTERP_KERNEL::InterpolationCU myInterpolator(*this);
        myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,_matrix,"P0P0");
        break;
      }
    case 2:
      {
        MEDCouplingNormalizedCartesianMesh<2> sourceWrapper(src_mesh);
        MEDCouplingNormalizedUnstructuredMesh<2,2> targetWrapper(target_mesh);
        INTERP_KERNEL::InterpolationCU myInterpolator(*this);
        myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,_matrix,"P0P0");
        break;
      }
    case 3:
      {
        MEDCouplingNormalizedCartesianMesh<3> sourceWrapper(src_mesh);
        MEDCouplingNormalizedUnstructuredMesh<3,3> targetWrapper(target_mesh);
        INTERP_KERNEL::InterpolationCU myInterpolator(*this);
        myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,_matrix,"P0P0");
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareInterpKernelOnlyCU : only dimension 1 2 or 3 supported !");
  }
  nullifiedTinyCoeffInCrudeMatrixAbs(0.);
  //
  _deno_multiply.clear();
  _deno_multiply.resize(_matrix.size());
  _deno_reverse_multiply.clear();
  _deno_reverse_multiply.resize(src_mesh->getNumberOfCells());
  declareAsNew();
  return 1;
}

int MEDCouplingRemapper::prepareInterpKernelOnlyCC()
{
  std::string srcMeth,trgMeth;
  std::string methodCpp=checkAndGiveInterpolationMethodStr(srcMeth,trgMeth);
  if(methodCpp!="P0P0")
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareInterpKernelOnlyCC : only P0P0 interpolation supported for the moment !");
  const MEDCouplingCMesh *src_mesh=static_cast<const MEDCouplingCMesh *>(_src_ft->getMesh());
  const MEDCouplingCMesh *target_mesh=static_cast<const MEDCouplingCMesh *>(_target_ft->getMesh());
  const int srcMeshDim=src_mesh->getMeshDimension();
  const int trgMeshDim=target_mesh->getMeshDimension();
  if(trgMeshDim!=srcMeshDim)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareInterpKernelOnlyCC : dim of target cartesian should be equal to dim of source cartesian dimension !");
  switch(srcMeshDim)
  {
    case 1:
      {
        MEDCouplingNormalizedCartesianMesh<1> sourceWrapper(src_mesh);
        MEDCouplingNormalizedCartesianMesh<1> targetWrapper(target_mesh);
        INTERP_KERNEL::InterpolationCC myInterpolator(*this);
        myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,_matrix,"P0P0");
        break;
      }
    case 2:
      {
        MEDCouplingNormalizedCartesianMesh<2> sourceWrapper(src_mesh);
        MEDCouplingNormalizedCartesianMesh<2> targetWrapper(target_mesh);
        INTERP_KERNEL::InterpolationCC myInterpolator(*this);
        myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,_matrix,"P0P0");
        break;
      }
    case 3:
      {
        MEDCouplingNormalizedCartesianMesh<3> sourceWrapper(src_mesh);
        MEDCouplingNormalizedCartesianMesh<3> targetWrapper(target_mesh);
        INTERP_KERNEL::InterpolationCC myInterpolator(*this);
        myInterpolator.interpolateMeshes(sourceWrapper,targetWrapper,_matrix,"P0P0");
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareInterpKernelOnlyCC : only dimension 1 2 or 3 supported !");
  }
  nullifiedTinyCoeffInCrudeMatrixAbs(0.);
  //
  _deno_multiply.clear();
  _deno_multiply.resize(_matrix.size());
  _deno_reverse_multiply.clear();
  _deno_reverse_multiply.resize(src_mesh->getNumberOfCells());
  declareAsNew();
  return 1;
}

int MEDCouplingRemapper::prepareNotInterpKernelOnlyGaussGauss()
{
  if(getIntersectionType()!=INTERP_KERNEL::PointLocator)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareNotInterpKernelOnlyGaussGauss : The intersection type is not supported ! Only PointLocator is supported for Gauss->Gauss interpolation ! Please invoke setIntersectionType(PointLocator) on the MEDCouplingRemapper instance !");
  MCAuto<DataArrayDouble> trgLoc=_target_ft->getLocalizationOfDiscr();
  const double *trgLocPtr=trgLoc->begin();
  int trgSpaceDim=trgLoc->getNumberOfComponents();
  MCAuto<DataArrayInt> srcOffsetArr=_src_ft->getDiscretization()->getOffsetArr(_src_ft->getMesh());
  if(trgSpaceDim!=_src_ft->getMesh()->getSpaceDimension())
    {
      std::ostringstream oss; oss << "MEDCouplingRemapper::prepareNotInterpKernelOnlyGaussGauss : space dimensions mismatch between source and target !";
      oss << " Target discretization localization has dimension " << trgSpaceDim << ", whereas the space dimension of source is equal to ";
      oss << _src_ft->getMesh()->getSpaceDimension() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const int *srcOffsetArrPtr=srcOffsetArr->begin();
  MCAuto<DataArrayDouble> srcLoc=_src_ft->getLocalizationOfDiscr();
  const double *srcLocPtr=srcLoc->begin();
  MCAuto<DataArrayInt> eltsArr,eltsIndexArr;
  int trgNbOfGaussPts=trgLoc->getNumberOfTuples();
  _matrix.resize(trgNbOfGaussPts);
  _src_ft->getMesh()->getCellsContainingPoints(trgLoc->begin(),trgNbOfGaussPts,getPrecision(),eltsArr,eltsIndexArr);
  const int *elts(eltsArr->begin()),*eltsIndex(eltsIndexArr->begin());
  MCAuto<DataArrayInt> nbOfSrcCellsShTrgPts(eltsIndexArr->deltaShiftIndex());
  MCAuto<DataArrayInt> ids0=nbOfSrcCellsShTrgPts->findIdsNotEqual(0);
  for(const int *trgId=ids0->begin();trgId!=ids0->end();trgId++)
    {
      const double *ptTrg=trgLocPtr+trgSpaceDim*(*trgId);
      int srcCellId=elts[eltsIndex[*trgId]];
      double dist=std::numeric_limits<double>::max();
      int srcEntry=-1;
      for(int srcId=srcOffsetArrPtr[srcCellId];srcId<srcOffsetArrPtr[srcCellId+1];srcId++)
        {
          const double *ptSrc=srcLocPtr+trgSpaceDim*srcId;
          double tmp=0.;
          for(int i=0;i<trgSpaceDim;i++)
            tmp+=(ptTrg[i]-ptSrc[i])*(ptTrg[i]-ptSrc[i]);
          if(tmp<dist)
            { dist=tmp; srcEntry=srcId; }
        }
      _matrix[*trgId][srcEntry]=1.;
    }
  if(ids0->getNumberOfTuples()!=trgNbOfGaussPts)
    {
      MCAuto<DataArrayInt> orphanTrgIds=nbOfSrcCellsShTrgPts->findIdsEqual(0);
      MCAuto<DataArrayDouble> orphanTrg=trgLoc->selectByTupleId(orphanTrgIds->begin(),orphanTrgIds->end());
      MCAuto<DataArrayInt> srcIdPerTrg=srcLoc->findClosestTupleId(orphanTrg);
      const int *srcIdPerTrgPtr=srcIdPerTrg->begin();
      for(const int *orphanTrgId=orphanTrgIds->begin();orphanTrgId!=orphanTrgIds->end();orphanTrgId++,srcIdPerTrgPtr++)
        _matrix[*orphanTrgId][*srcIdPerTrgPtr]=2.;
    }
  _deno_multiply.clear();
  _deno_multiply.resize(_matrix.size());
  _deno_reverse_multiply.clear();
  _deno_reverse_multiply.resize(srcLoc->getNumberOfTuples());
  declareAsNew();
  return 1;
}

/*!
 * This method checks that the input interpolation \a method is managed by not INTERP_KERNEL only methods.
 * If no an INTERP_KERNEL::Exception will be thrown. If yes, a magic number will be returned to switch in the MEDCouplingRemapper::prepareNotInterpKernelOnly method.
 */
int MEDCouplingRemapper::CheckInterpolationMethodManageableByNotOnlyInterpKernel(const std::string& method)
{
  if(method=="GAUSSGAUSS")
    return 0;
  std::ostringstream oss; oss << "MEDCouplingRemapper::CheckInterpolationMethodManageableByNotOnlyInterpKernel : ";
  oss << "The method \"" << method << "\" is not manageable by not INTERP_KERNEL only method.";
  oss << " Not only INTERP_KERNEL methods dealed are : GAUSSGAUSS !";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

/*!
 * This method determines regarding \c _interp_matrix_pol attribute ( set by MEDCouplingRemapper::setInterpolationMatrixPolicy and by default equal
 * to IK_ONLY_PREFERED, which method will be applied. If \c true is returned the INTERP_KERNEL only method should be applied to \c false the \b not
 * only INTERP_KERNEL method should be applied.
 */
bool MEDCouplingRemapper::isInterpKernelOnlyOrNotOnly() const
{
  std::string srcm,trgm,method;
  method=checkAndGiveInterpolationMethodStr(srcm,trgm);
  switch(_interp_matrix_pol)
  {
    case IK_ONLY_PREFERED:
      {
        try
        {
            std::string tmp1,tmp2;
            INTERP_KERNEL::Interpolation<INTERP_KERNEL::Interpolation3D>::CheckAndSplitInterpolationMethod(method,tmp1,tmp2);
            return true;
        }
        catch(INTERP_KERNEL::Exception& /*e*/)
        {
            return false;
        }
      }
    case NOT_IK_ONLY_PREFERED:
      {
        try
        {
            CheckInterpolationMethodManageableByNotOnlyInterpKernel(method);
            return false;
        }
        catch(INTERP_KERNEL::Exception& /*e*/)
        {
            return true;
        }
      }
    case IK_ONLY_FORCED:
      return true;
    case NOT_IK_ONLY_FORCED:
      return false;
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingRemapper::isInterpKernelOnlyOrNotOnly : internal error ! The interpolation matrix policy is not managed ! Try to change it using MEDCouplingRemapper::setInterpolationMatrixPolicy !");
  }
}

void MEDCouplingRemapper::updateTime() const
{
}

void MEDCouplingRemapper::checkPrepare() const
{
  const MEDCouplingFieldTemplate *s(_src_ft),*t(_target_ft);
  if(!s || !t)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::checkPrepare : it appears that MEDCouplingRemapper::prepare(Ex) has not been called !");
  if(!s->getMesh() || !t->getMesh())
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::checkPrepare : it appears that no all field templates have their mesh set !");
}

/*!
 * This method builds a code considering already set field discretization int \a this : \a _src_ft and \a _target_ft.
 * This method returns 3 informations (2 in ouput parameters and 1 in return).
 * 
 * \param [out] srcMeth the string code of the discretization of source field template
 * \param [out] trgMeth the string code of the discretization of target field template
 * \return the standardized string code (compatible with INTERP_KERNEL) for matrix of numerators (in \a _matrix)
 */
std::string MEDCouplingRemapper::checkAndGiveInterpolationMethodStr(std::string& srcMeth, std::string& trgMeth) const
{
  const MEDCouplingFieldTemplate *s(_src_ft),*t(_target_ft);
  if(!s || !t)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::checkAndGiveInterpolationMethodStr : it appears that no all field templates have been set !");
  if(!s->getMesh() || !t->getMesh())
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::checkAndGiveInterpolationMethodStr : it appears that no all field templates have their mesh set !");
  srcMeth=_src_ft->getDiscretization()->getRepr();
  trgMeth=_target_ft->getDiscretization()->getRepr();
  return BuildMethodFrom(srcMeth,trgMeth);
}

std::string MEDCouplingRemapper::BuildMethodFrom(const std::string& meth1, const std::string& meth2)
{
  std::string method(meth1); method+=meth2;
  return method;
}

void MEDCouplingRemapper::releaseData(bool matrixSuppression)
{
  _src_ft=0;
  _target_ft=0;
  if(matrixSuppression)
    {
      _matrix.clear();
      _deno_multiply.clear();
      _deno_reverse_multiply.clear();
    }
}

void MEDCouplingRemapper::transferUnderground(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField, bool isDftVal, double dftValue)
{
  if(!srcField || !targetField)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::transferUnderground : srcField or targetField is NULL !");
  srcField->checkConsistencyLight();
  checkPrepare();
  if(_src_ft->getDiscretization()->getStringRepr()!=srcField->getDiscretization()->getStringRepr())
    throw INTERP_KERNEL::Exception("Incoherency with prepare call for source field");
  if(_target_ft->getDiscretization()->getStringRepr()!=targetField->getDiscretization()->getStringRepr())
    throw INTERP_KERNEL::Exception("Incoherency with prepare call for target field");
  if(srcField->getNature()!=targetField->getNature())
    throw INTERP_KERNEL::Exception("Natures of fields mismatch !");
  if(srcField->getNumberOfTuplesExpected()!=_src_ft->getNumberOfTuplesExpected())
    {
      std::ostringstream oss;
      oss << "MEDCouplingRemapper::transferUnderground : in given source field the number of tuples required is " << _src_ft->getNumberOfTuplesExpected() << " (on prepare) and number of tuples in given source field is " << srcField->getNumberOfTuplesExpected();
      oss << " ! It appears that the source support is not the same between the prepare and the transfer !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  DataArrayDouble *array(targetField->getArray());
  int srcNbOfCompo(srcField->getNumberOfComponents());
  if(array)
    {
      targetField->checkConsistencyLight();
      if(srcNbOfCompo!=targetField->getNumberOfComponents())
        throw INTERP_KERNEL::Exception("Number of components mismatch !");
    }
  else
    {
      if(!isDftVal)
        throw INTERP_KERNEL::Exception("MEDCouplingRemapper::partialTransfer : This method requires that the array of target field exists ! Allocate it or call MEDCouplingRemapper::transfer instead !");
      MCAuto<DataArrayDouble> tmp(DataArrayDouble::New());
      tmp->alloc(targetField->getNumberOfTuples(),srcNbOfCompo);
      targetField->setArray(tmp);
    }
  computeDeno(srcField->getNature(),srcField,targetField);
  double *resPointer(targetField->getArray()->getPointer());
  const double *inputPointer(srcField->getArray()->getConstPointer());
  computeProduct(inputPointer,srcNbOfCompo,isDftVal,dftValue,resPointer);
}

void MEDCouplingRemapper::computeDeno(NatureOfField nat, const MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *trgField)
{
  if(nat==NoNature)
    return computeDenoFromScratch(nat,srcField,trgField);
  else if(nat!=_nature_of_deno)
    return computeDenoFromScratch(nat,srcField,trgField);
  else if(nat==_nature_of_deno && _time_deno_update!=getTimeOfThis())
    return computeDenoFromScratch(nat,srcField,trgField);
}

void MEDCouplingRemapper::computeDenoFromScratch(NatureOfField nat, const MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *trgField)
{
  _nature_of_deno=nat;
  _time_deno_update=getTimeOfThis();
  switch(_nature_of_deno)
  {
    case IntensiveMaximum:
      {
        ComputeRowSumAndColSum(_matrix,_deno_multiply,_deno_reverse_multiply);
        break;
      }
    case ExtensiveMaximum:
      {
        MEDCouplingFieldDouble *deno=srcField->getDiscretization()->getMeasureField(srcField->getMesh(),getMeasureAbsStatus());
        MEDCouplingFieldDouble *denoR=trgField->getDiscretization()->getMeasureField(trgField->getMesh(),getMeasureAbsStatus());
        const double *denoPtr=deno->getArray()->getConstPointer();
        const double *denoRPtr=denoR->getArray()->getConstPointer();
        if(trgField->getMesh()->getMeshDimension()==-1)
          {
            double *denoRPtr2=denoR->getArray()->getPointer();
            denoRPtr2[0]=std::accumulate(denoPtr,denoPtr+deno->getNumberOfTuples(),0.);
          }
        if(srcField->getMesh()->getMeshDimension()==-1)
          {
            double *denoPtr2=deno->getArray()->getPointer();
            denoPtr2[0]=std::accumulate(denoRPtr,denoRPtr+denoR->getNumberOfTuples(),0.);
          }
        int idx=0;
        for(std::vector<std::map<int,double> >::const_iterator iter1=_matrix.begin();iter1!=_matrix.end();iter1++,idx++)
          for(std::map<int,double>::const_iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
            {
              _deno_multiply[idx][(*iter2).first]=denoPtr[(*iter2).first];
              _deno_reverse_multiply[(*iter2).first][idx]=denoRPtr[idx];
            }
        deno->decrRef();
        denoR->decrRef();
        break;
      }
    case ExtensiveConservation:
      {
        ComputeColSumAndRowSum(_matrix,_deno_multiply,_deno_reverse_multiply);
        break;
      }
    case IntensiveConservation:
      {
        MEDCouplingFieldDouble *deno=trgField->getDiscretization()->getMeasureField(trgField->getMesh(),getMeasureAbsStatus());
        MEDCouplingFieldDouble *denoR=srcField->getDiscretization()->getMeasureField(srcField->getMesh(),getMeasureAbsStatus());
        const double *denoPtr=deno->getArray()->getConstPointer();
        const double *denoRPtr=denoR->getArray()->getConstPointer();
        if(trgField->getMesh()->getMeshDimension()==-1)
          {
            double *denoRPtr2=denoR->getArray()->getPointer();
            denoRPtr2[0]=std::accumulate(denoPtr,denoPtr+deno->getNumberOfTuples(),0.);
          }
        if(srcField->getMesh()->getMeshDimension()==-1)
          {
            double *denoPtr2=deno->getArray()->getPointer();
            denoPtr2[0]=std::accumulate(denoRPtr,denoRPtr+denoR->getNumberOfTuples(),0.);
          }
        int idx=0;
        for(std::vector<std::map<int,double> >::const_iterator iter1=_matrix.begin();iter1!=_matrix.end();iter1++,idx++)
          for(std::map<int,double>::const_iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
            {
              _deno_multiply[idx][(*iter2).first]=denoPtr[idx];
              _deno_reverse_multiply[(*iter2).first][idx]=denoRPtr[(*iter2).first];
            }
        deno->decrRef();
        denoR->decrRef();
        break;
      }
    case NoNature:
      throw INTERP_KERNEL::Exception("No nature specified ! Select one !");
  }
}

void MEDCouplingRemapper::computeProduct(const double *inputPointer, int inputNbOfCompo, bool isDftVal, double dftValue, double *resPointer)
{
  int idx=0;
  double *tmp=new double[inputNbOfCompo];
  for(std::vector<std::map<int,double> >::const_iterator iter1=_matrix.begin();iter1!=_matrix.end();iter1++,idx++)
    {
      if((*iter1).empty())
        {
          if(isDftVal)
            std::fill(resPointer+idx*inputNbOfCompo,resPointer+(idx+1)*inputNbOfCompo,dftValue);
          continue;
        }
      else
        std::fill(resPointer+idx*inputNbOfCompo,resPointer+(idx+1)*inputNbOfCompo,0.);
      std::map<int,double>::const_iterator iter3=_deno_multiply[idx].begin();
      for(std::map<int,double>::const_iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++,iter3++)
        {
          std::transform(inputPointer+(*iter2).first*inputNbOfCompo,inputPointer+((*iter2).first+1)*inputNbOfCompo,tmp,std::bind2nd(std::multiplies<double>(),(*iter2).second/(*iter3).second));
          std::transform(tmp,tmp+inputNbOfCompo,resPointer+idx*inputNbOfCompo,resPointer+idx*inputNbOfCompo,std::plus<double>());
        }
    }
  delete [] tmp;
}

void MEDCouplingRemapper::computeReverseProduct(const double *inputPointer, int inputNbOfCompo, double dftValue, double *resPointer)
{
  std::vector<bool> isReached(_deno_reverse_multiply.size(),false);
  int idx=0;
  double *tmp=new double[inputNbOfCompo];
  std::fill(resPointer,resPointer+inputNbOfCompo*_deno_reverse_multiply.size(),0.);
  for(std::vector<std::map<int,double> >::const_iterator iter1=_matrix.begin();iter1!=_matrix.end();iter1++,idx++)
    {
      for(std::map<int,double>::const_iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
        {
          isReached[(*iter2).first]=true;
          std::transform(inputPointer+idx*inputNbOfCompo,inputPointer+(idx+1)*inputNbOfCompo,tmp,std::bind2nd(std::multiplies<double>(),(*iter2).second/_deno_reverse_multiply[(*iter2).first][idx]));
          std::transform(tmp,tmp+inputNbOfCompo,resPointer+((*iter2).first)*inputNbOfCompo,resPointer+((*iter2).first)*inputNbOfCompo,std::plus<double>());
        }
    }
  delete [] tmp;
  idx=0;
  for(std::vector<bool>::const_iterator iter3=isReached.begin();iter3!=isReached.end();iter3++,idx++)
    if(!*iter3)
      std::fill(resPointer+idx*inputNbOfCompo,resPointer+(idx+1)*inputNbOfCompo,dftValue);
}

void MEDCouplingRemapper::ReverseMatrix(const std::vector<std::map<int,double> >& matIn, int nbColsMatIn, std::vector<std::map<int,double> >& matOut)
{
  matOut.resize(nbColsMatIn);
  int id=0;
  for(std::vector<std::map<int,double> >::const_iterator iter1=matIn.begin();iter1!=matIn.end();iter1++,id++)
    for(std::map<int,double>::const_iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
      matOut[(*iter2).first][id]=(*iter2).second;
}

void MEDCouplingRemapper::ComputeRowSumAndColSum(const std::vector<std::map<int,double> >& matrixDeno,
                                                 std::vector<std::map<int,double> >& deno, std::vector<std::map<int,double> >& denoReverse)
{
  std::map<int,double> values;
  int idx=0;
  for(std::vector<std::map<int,double> >::const_iterator iter1=matrixDeno.begin();iter1!=matrixDeno.end();iter1++,idx++)
    {
      double sum=0.;
      for(std::map<int,double>::const_iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
        {
          sum+=(*iter2).second;
          values[(*iter2).first]+=(*iter2).second;
        }
      for(std::map<int,double>::const_iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
        deno[idx][(*iter2).first]=sum;
    }
  idx=0;
  for(std::vector<std::map<int,double> >::const_iterator iter1=matrixDeno.begin();iter1!=matrixDeno.end();iter1++,idx++)
    {
      for(std::map<int,double>::const_iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
        denoReverse[(*iter2).first][idx]=values[(*iter2).first];
    }
}

void MEDCouplingRemapper::ComputeColSumAndRowSum(const std::vector<std::map<int,double> >& matrixDeno,
                                                 std::vector<std::map<int,double> >& deno, std::vector<std::map<int,double> >& denoReverse)
{
  std::map<int,double> values;
  int idx=0;
  for(std::vector<std::map<int,double> >::const_iterator iter1=matrixDeno.begin();iter1!=matrixDeno.end();iter1++,idx++)
    {
      double sum=0.;
      for(std::map<int,double>::const_iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
        {
          sum+=(*iter2).second;
          values[(*iter2).first]+=(*iter2).second;
        }
      for(std::map<int,double>::const_iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
        denoReverse[(*iter2).first][idx]=sum;
    }
  idx=0;
  for(std::vector<std::map<int,double> >::const_iterator iter1=matrixDeno.begin();iter1!=matrixDeno.end();iter1++,idx++)
    {
      for(std::map<int,double>::const_iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
        deno[idx][(*iter2).first]=values[(*iter2).first];
    }
}

void MEDCouplingRemapper::buildFinalInterpolationMatrixByConvolution(const std::vector< std::map<int,double> >& m1D,
                                                                     const std::vector< std::map<int,double> >& m2D,
                                                                     const int *corrCellIdSrc, int nbOf2DCellsSrc, int nbOf1DCellsSrc,
                                                                     const int *corrCellIdTrg)
{
  int nbOf2DCellsTrg=m2D.size();
  int nbOf1DCellsTrg=m1D.size();
  int nbOf3DCellsTrg=nbOf2DCellsTrg*nbOf1DCellsTrg;
  _matrix.resize(nbOf3DCellsTrg);
  int id2R=0;
  for(std::vector< std::map<int,double> >::const_iterator iter2R=m2D.begin();iter2R!=m2D.end();iter2R++,id2R++)
    {
      for(std::map<int,double>::const_iterator iter2C=(*iter2R).begin();iter2C!=(*iter2R).end();iter2C++)
        {
          int id1R=0;
          for(std::vector< std::map<int,double> >::const_iterator iter1R=m1D.begin();iter1R!=m1D.end();iter1R++,id1R++)
            {
              for(std::map<int,double>::const_iterator iter1C=(*iter1R).begin();iter1C!=(*iter1R).end();iter1C++)
                {
                  _matrix[corrCellIdTrg[id1R*nbOf2DCellsTrg+id2R]][corrCellIdSrc[(*iter1C).first*nbOf2DCellsSrc+(*iter2C).first]]=(*iter1C).second*((*iter2C).second);
                }
            }
        }
    }
}

void MEDCouplingRemapper::PrintMatrix(const std::vector<std::map<int,double> >& m)
{
  int id=0;
  for(std::vector<std::map<int,double> >::const_iterator iter1=m.begin();iter1!=m.end();iter1++,id++)
    {
      std::cout << "Target Cell # " << id << " : ";
      for(std::map<int,double>::const_iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
        std::cout << "(" << (*iter2).first << "," << (*iter2).second << "), ";
      std::cout << std::endl;
    }
}

const std::vector<std::map<int,double> >& MEDCouplingRemapper::getCrudeMatrix() const
{
  return _matrix;
}

/*!
 * Returns the number of columns of matrix returned by MEDCouplingRemapper::getCrudeMatrix method.
 */
int MEDCouplingRemapper::getNumberOfColsOfMatrix() const
{
  return (int)_deno_reverse_multiply.size();
}

/*!
 * This method is supposed to be called , if needed, right after MEDCouplingRemapper::prepare or MEDCouplingRemapper::prepareEx.
 * If not the behaviour is unpredictable.
 * This method works on precomputed \a this->_matrix. All coefficients in the matrix is lower than \a maxValAbs this coefficient is
 * set to 0. That is to say that its entry disappear from the map storing the corresponding row in the data storage of sparse crude matrix.
 * This method is useful to correct at a high level some problems linked to precision. Indeed, with some \ref NatureOfField "natures of field" some threshold effect
 * can occur.
 *
 * \param [in] maxValAbs is a limit behind which a coefficient is set to 0. \a maxValAbs is expected to be positive, if not this method do nothing.
 * \return a positive value that tells the number of coefficients put to 0. The 0 returned value means that the matrix has remained unchanged.
 * \sa MEDCouplingRemapper::nullifiedTinyCoeffInCrudeMatrix
 */
int MEDCouplingRemapper::nullifiedTinyCoeffInCrudeMatrixAbs(double maxValAbs)
{
  int ret=0;
  std::vector<std::map<int,double> > matrixNew(_matrix.size());
  int i=0;
  for(std::vector<std::map<int,double> >::const_iterator it1=_matrix.begin();it1!=_matrix.end();it1++,i++)
    {
      std::map<int,double>& rowNew=matrixNew[i];
      for(std::map<int,double>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
        {
          if(fabs((*it2).second)>maxValAbs)
            rowNew[(*it2).first]=(*it2).second;
          else
            ret++;
        }
    }
  if(ret>0)
    _matrix=matrixNew;
  return ret;
}

/*!
 * This method is supposed to be called , if needed, right after MEDCouplingRemapper::prepare or MEDCouplingRemapper::prepareEx.
 * If not the behaviour is unpredictable.
 * This method works on precomputed \a this->_matrix. All coefficients in the matrix is lower than delta multiplied by \a scaleFactor this coefficient is
 * set to 0. That is to say that its entry disappear from the map storing the corresponding row in the data storage of sparse crude matrix.
 * delta is the value returned by MEDCouplingRemapper::getMaxValueInCrudeMatrix method.
 * This method is useful to correct at a high level some problems linked to precision. Indeed, with some \ref NatureOfField "natures of field" some threshold effect
 * can occur.
 *
 * \param [in] scaleFactor is the scale factor from which coefficients lower than \a scaleFactor times range width of coefficients are set to zero.
 * \return a positive value that tells the number of coefficients put to 0. The 0 returned value means that the matrix has remained unchanged. If -1 is returned it means
 *         that all coefficients are null.
 * \sa MEDCouplingRemapper::nullifiedTinyCoeffInCrudeMatrixAbs
 */
int MEDCouplingRemapper::nullifiedTinyCoeffInCrudeMatrix(double scaleFactor)
{
  double maxVal=getMaxValueInCrudeMatrix();
  if(maxVal==0.)
    return -1;
  return nullifiedTinyCoeffInCrudeMatrixAbs(scaleFactor*maxVal);
}

/*!
 * This method is supposed to be called , if needed, right after MEDCouplingRemapper::prepare or MEDCouplingRemapper::prepareEx.
 * If not the behaviour is unpredictable.
 * This method returns the maximum of the absolute values of coefficients into the sparse crude matrix.
 * The returned value is positive.
 */
double MEDCouplingRemapper::getMaxValueInCrudeMatrix() const
{
  double ret=0.;
  for(std::vector<std::map<int,double> >::const_iterator it1=_matrix.begin();it1!=_matrix.end();it1++)
    for(std::map<int,double>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
      if(fabs((*it2).second)>ret)
        ret=fabs((*it2).second);
  return ret;
}
