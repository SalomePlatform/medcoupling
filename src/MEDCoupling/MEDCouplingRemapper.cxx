//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#include "MEDCouplingRemapper.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
#include "MEDCouplingNormalizedUnstructuredMesh.txx"

#include "Interpolation2D.txx"
#include "Interpolation3D.txx"
#include "Interpolation3DSurf.txx"

using namespace ParaMEDMEM;

MEDCouplingRemapper::MEDCouplingRemapper():_src_mesh(0),_target_mesh(0),_nature_of_deno(NoNature),_time_deno_update(0)
{
}

MEDCouplingRemapper::~MEDCouplingRemapper()
{
  releaseData(false);
}

int MEDCouplingRemapper::prepare(const MEDCouplingUMesh *srcMesh, const MEDCouplingUMesh *targetMesh, const char *method)
{
  releaseData(true);
  _src_mesh=(MEDCouplingUMesh *)srcMesh; _target_mesh=(MEDCouplingUMesh *)targetMesh;
  _src_mesh->incrRef(); _target_mesh->incrRef();
  INTERP_KERNEL::Interpolation<INTERP_KERNEL::Interpolation3D>::checkAndSplitInterpolationMethod(method,_src_method,_target_method);
  const int srcMeshDim=_src_mesh->getMeshDimension();
  const int srcSpaceDim=_src_mesh->getSpaceDimension();
  const int trgMeshdim=_target_mesh->getMeshDimension();
  const int trgSpaceDim=_target_mesh->getSpaceDimension();
  if(trgSpaceDim!=srcSpaceDim)
    throw INTERP_KERNEL::Exception("Incoherent space dimension detected between target and source.");
  int nbCols;
  if(srcMeshDim==2 && trgMeshdim==2 && srcSpaceDim==2)
    {
      MEDCouplingNormalizedUnstructuredMesh<2,2> source_mesh_wrapper(_src_mesh);
      MEDCouplingNormalizedUnstructuredMesh<2,2> target_mesh_wrapper(_target_mesh);
      INTERP_KERNEL::Interpolation2D interpolation(*this);
      nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
    }
  else if(srcMeshDim==3 && trgMeshdim==3 && srcSpaceDim==3)
    {
      MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(_src_mesh);
      MEDCouplingNormalizedUnstructuredMesh<3,3> target_mesh_wrapper(_target_mesh);
      INTERP_KERNEL::Interpolation3D interpolation(*this);
      nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
    }
  else if(srcMeshDim==2 && trgMeshdim==2 && srcSpaceDim==3)
    {
      MEDCouplingNormalizedUnstructuredMesh<3,2> source_mesh_wrapper(_src_mesh);
      MEDCouplingNormalizedUnstructuredMesh<3,2> target_mesh_wrapper(_target_mesh);
      INTERP_KERNEL::Interpolation3DSurf interpolation(*this);
      nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
    }
  else if(srcMeshDim==3 && trgMeshdim==1 && srcSpaceDim==3)
    {
      if(getIntersectionType()!=INTERP_KERNEL::PointLocator)
        throw INTERP_KERNEL::Exception("Invalid interpolation requested between 3D and 1D ! Select PointLocator as intersection type !");
      MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(_src_mesh);
      MEDCouplingNormalizedUnstructuredMesh<3,3> target_mesh_wrapper(_target_mesh);
      INTERP_KERNEL::Interpolation3D interpolation(*this);
      nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
    }
  else if(srcMeshDim==1 && trgMeshdim==3 && srcSpaceDim==3)
    {
      if(getIntersectionType()!=INTERP_KERNEL::PointLocator)
        throw INTERP_KERNEL::Exception("Invalid interpolation requested between 3D and 1D ! Select PointLocator as intersection type !");
      MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(_src_mesh);
      MEDCouplingNormalizedUnstructuredMesh<3,3> target_mesh_wrapper(_target_mesh);
      INTERP_KERNEL::Interpolation3D interpolation(*this);
      std::vector<std::map<int,double> > matrixTmp;
      nbCols=interpolation.interpolateMeshes(target_mesh_wrapper,source_mesh_wrapper,matrixTmp,method);
      reverseMatrix(matrixTmp,nbCols,_matrix);
      nbCols=matrixTmp.size();
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

void MEDCouplingRemapper::transfer(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField, double dftValue)
{
  if(_src_method!=srcField->getDiscretization()->getStringRepr())
    throw INTERP_KERNEL::Exception("Incoherency with prepare call for source field");
  if(_target_method!=targetField->getDiscretization()->getStringRepr())
    throw INTERP_KERNEL::Exception("Incoherency with prepare call for target field");
  if(srcField->getNature()!=targetField->getNature())
    throw INTERP_KERNEL::Exception("Natures of fields mismatch !");
  DataArrayDouble *array=targetField->getArray();
  int srcNbOfCompo=srcField->getNumberOfComponents();
  if(array)
    {
      if(srcNbOfCompo!=targetField->getNumberOfComponents())
        throw INTERP_KERNEL::Exception("Number of components mismatch !");
    }
  else
    {
      array=DataArrayDouble::New();
      array->alloc(targetField->getNumberOfTuples(),srcNbOfCompo);
      targetField->setArray(array);
      array->decrRef();
    }
  computeDeno(srcField->getNature(),srcField,targetField);
  double *resPointer=array->getPointer();
  const double *inputPointer=srcField->getArray()->getConstPointer();
  computeProduct(inputPointer,srcNbOfCompo,dftValue,resPointer);
}

void MEDCouplingRemapper::reverseTransfer(MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *targetField, double dftValue)
{
  if(_src_method!=srcField->getDiscretization()->getStringRepr())
    throw INTERP_KERNEL::Exception("Incoherency with prepare call for source field");
  if(_target_method!=targetField->getDiscretization()->getStringRepr())
    throw INTERP_KERNEL::Exception("Incoherency with prepare call for target field");
  if(srcField->getNature()!=targetField->getNature())
    throw INTERP_KERNEL::Exception("Natures of fields mismatch !");
  DataArrayDouble *array=srcField->getArray();
  int trgNbOfCompo=targetField->getNumberOfComponents();
  if(array)
    {
      if(trgNbOfCompo!=srcField->getNumberOfComponents())
        throw INTERP_KERNEL::Exception("Number of components mismatch !");
    }
  else
    {
      array=DataArrayDouble::New();
      array->alloc(srcField->getNumberOfTuples(),trgNbOfCompo);
      srcField->setArray(array);
      array->decrRef();
    }
  computeDeno(srcField->getNature(),srcField,targetField);
  double *resPointer=array->getPointer();
  const double *inputPointer=targetField->getArray()->getConstPointer();
  computeReverseProduct(inputPointer,trgNbOfCompo,dftValue,resPointer);
}

MEDCouplingFieldDouble *MEDCouplingRemapper::transferField(const MEDCouplingFieldDouble *srcField, double dftValue)
{
  if(_src_method!=srcField->getDiscretization()->getStringRepr())
    throw INTERP_KERNEL::Exception("Incoherency with prepare call for source field");
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(MEDCouplingFieldDiscretization::getTypeOfFieldFromStringRepr(_target_method.c_str()));
  ret->setNature(srcField->getNature());
  ret->setMesh(_target_mesh);
  transfer(srcField,ret,dftValue);
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingRemapper::reverseTransferField(const MEDCouplingFieldDouble *targetField, double dftValue)
{
  if(_target_method!=targetField->getDiscretization()->getStringRepr())
    throw INTERP_KERNEL::Exception("Incoherency with prepare call for target field");
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(MEDCouplingFieldDiscretization::getTypeOfFieldFromStringRepr(_target_method.c_str()));
  ret->setNature(targetField->getNature());
  ret->setMesh(_src_mesh);
  reverseTransfer(ret,targetField,dftValue);
  return ret;
}

bool MEDCouplingRemapper::setOptionInt(const std::string& key, int value)
{
  return INTERP_KERNEL::InterpolationOptions::setOptionInt(key,value);
}

bool MEDCouplingRemapper::setOptionDouble(const std::string& key, double value)
{
  return INTERP_KERNEL::InterpolationOptions::setOptionDouble(key,value);
}

bool MEDCouplingRemapper::setOptionString(const std::string& key, std::string& value)
{
  return INTERP_KERNEL::InterpolationOptions::setOptionString(key,value);
}

void MEDCouplingRemapper::updateTime()
{
}

void MEDCouplingRemapper::releaseData(bool matrixSuppression)
{
  if(_src_mesh)
    _src_mesh->decrRef();
  if(_target_mesh)
    _target_mesh->decrRef();
  _src_mesh=0;
  _target_mesh=0;
  if(matrixSuppression)
    {
      _matrix.clear();
      _deno_multiply.clear();
      _deno_reverse_multiply.clear();
    }
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
    case ConservativeVolumic:
      {
        computeRowSumAndColSum(_matrix,_deno_multiply,_deno_reverse_multiply);
        break;
      }
    case Integral:
      {
        MEDCouplingFieldDouble *deno=srcField->getDiscretization()->getWeightingField(srcField->getMesh(),true);
        MEDCouplingFieldDouble *denoR=trgField->getDiscretization()->getWeightingField(trgField->getMesh(),true);
        const double *denoPtr=deno->getArray()->getConstPointer();
        const double *denoRPtr=denoR->getArray()->getConstPointer();
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
    case IntegralGlobConstraint:
      {
        computeColSumAndRowSum(_matrix,_deno_multiply,_deno_reverse_multiply);
        break;
      }
    case NoNature:
      throw INTERP_KERNEL::Exception("No nature specified ! Select one !");
    }
}

void MEDCouplingRemapper::computeProduct(const double *inputPointer, int inputNbOfCompo, double dftValue, double *resPointer)
{
  int idx=0;
  double *tmp=new double[inputNbOfCompo];
  for(std::vector<std::map<int,double> >::const_iterator iter1=_matrix.begin();iter1!=_matrix.end();iter1++,idx++)
    {
      if((*iter1).empty())
        {
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

void MEDCouplingRemapper::reverseMatrix(const std::vector<std::map<int,double> >& matIn, int nbColsMatIn, std::vector<std::map<int,double> >& matOut)
{
  matOut.resize(nbColsMatIn);
  int id=0;
  for(std::vector<std::map<int,double> >::const_iterator iter1=matIn.begin();iter1!=matIn.end();iter1++,id++)
    for(std::map<int,double>::const_iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
      matOut[(*iter2).first][id]=(*iter2).second;
}

void MEDCouplingRemapper::computeRowSumAndColSum(const std::vector<std::map<int,double> >& matrixDeno,
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

void MEDCouplingRemapper::computeColSumAndRowSum(const std::vector<std::map<int,double> >& matrixDeno,
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
