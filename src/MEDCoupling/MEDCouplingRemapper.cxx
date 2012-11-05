// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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
#include "MEDCouplingExtrudedMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingNormalizedUnstructuredMesh.txx"
#include "MEDCouplingNormalizedCartesianMesh.txx"

#include "Interpolation1D.txx"
#include "Interpolation2DCurve.hxx"
#include "Interpolation2D.txx"
#include "Interpolation3D.txx"
#include "Interpolation3DSurf.hxx"
#include "Interpolation2D1D.txx"
#include "Interpolation3D2D.txx"
#include "InterpolationCU.txx"
#include "InterpolationCC.txx"

using namespace ParaMEDMEM;

MEDCouplingRemapper::MEDCouplingRemapper():_src_mesh(0),_target_mesh(0),_nature_of_deno(NoNature),_time_deno_update(0)
{
}

MEDCouplingRemapper::~MEDCouplingRemapper()
{
  releaseData(false);
}

int MEDCouplingRemapper::prepare(const MEDCouplingMesh *srcMesh, const MEDCouplingMesh *targetMesh, const char *method) throw(INTERP_KERNEL::Exception)
{
  releaseData(true);
  _src_mesh=const_cast<MEDCouplingMesh *>(srcMesh); _target_mesh=const_cast<MEDCouplingMesh *>(targetMesh);
  _src_mesh->incrRef(); _target_mesh->incrRef();
  int meshInterpType=((int)_src_mesh->getType()*16)+(int)_target_mesh->getType();
  switch(meshInterpType)
    {
    case 85://Unstructured-Unstructured
      return prepareUU(method);
    case 87://Unstructured-Cartesian
      return prepareUC(method);
    case 117://Cartesian-Unstructured
      return prepareCU(method);
    case 119://Cartesian-Cartesian
      return prepareCC(method);
    case 136://Extruded-Extruded
      return prepareEE(method);
    default:
      throw INTERP_KERNEL::Exception("Not managed type of meshes !");
    }
}

int MEDCouplingRemapper::prepareEx(const MEDCouplingFieldTemplate *src, const MEDCouplingFieldTemplate *target) throw(INTERP_KERNEL::Exception)
{
  std::string meth(src->getDiscretization()->getStringRepr());
  meth+=target->getDiscretization()->getStringRepr();
  return prepare(src->getMesh(),target->getMesh(),meth.c_str());
}

/*!
 * This method performs the operation source to target using matrix computed in ParaMEDMEM::MEDCouplingRemapper::prepare method.
 * If meshes of \b srcField and \b targetField do not match exactly those given into \ref ParaMEDMEM::MEDCouplingRemapper::prepare "prepare method" an exception will be thrown.
 * 
 * \param [in] srcField is the source field from which the interpolation will be done. The mesh into \b srcField should be the same than those specified on ParaMEDMEM::MEDCouplingRemapper::prepare.
 * \param [out] targetField the destination field with the allocated array in which all tuples will be overwritten.
 */
void MEDCouplingRemapper::transfer(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField, double dftValue) throw(INTERP_KERNEL::Exception)
{
  transferUnderground(srcField,targetField,true,dftValue);
}

/*!
 * This method is equivalent to ParaMEDMEM::MEDCouplingRemapper::transfer except that here \b targetField is a in/out parameter.
 * If an entity (cell for example) in targetField is not fetched by any entity (cell for example) of \b srcField, the value in targetField is
 * let unchanged.
 * This method requires that \b targetField was fully defined and allocated. If the array is not allocated an exception will be thrown.
 * 
 * \param [in] srcField is the source field from which the interpolation will be done. The mesh into \b srcField should be the same than those specified on ParaMEDMEM::MEDCouplingRemapper::prepare.
 * \param [in,out] targetField the destination field with the allocated array in which only tuples whose entities are fetched by interpolation will be overwritten only.
 */
void MEDCouplingRemapper::partialTransfer(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField) throw(INTERP_KERNEL::Exception)
{
  transferUnderground(srcField,targetField,false,std::numeric_limits<double>::max());
}

void MEDCouplingRemapper::reverseTransfer(MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *targetField, double dftValue) throw(INTERP_KERNEL::Exception)
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

MEDCouplingFieldDouble *MEDCouplingRemapper::transferField(const MEDCouplingFieldDouble *srcField, double dftValue) throw(INTERP_KERNEL::Exception)
{
  if(_src_method!=srcField->getDiscretization()->getStringRepr())
    throw INTERP_KERNEL::Exception("Incoherency with prepare call for source field");
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(MEDCouplingFieldDiscretization::getTypeOfFieldFromStringRepr(_target_method.c_str()),srcField->getTimeDiscretization());
  ret->copyTinyAttrFrom(srcField);
  ret->setNature(srcField->getNature());
  ret->setMesh(_target_mesh);
  transfer(srcField,ret,dftValue);
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingRemapper::reverseTransferField(const MEDCouplingFieldDouble *targetField, double dftValue) throw(INTERP_KERNEL::Exception)
{
  if(_target_method!=targetField->getDiscretization()->getStringRepr())
    throw INTERP_KERNEL::Exception("Incoherency with prepare call for target field");
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(MEDCouplingFieldDiscretization::getTypeOfFieldFromStringRepr(_target_method.c_str()),targetField->getTimeDiscretization());
  ret->copyTinyAttrFrom(targetField);
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

bool MEDCouplingRemapper::setOptionString(const std::string& key, const std::string& value)
{
  return INTERP_KERNEL::InterpolationOptions::setOptionString(key,value);
}

int MEDCouplingRemapper::prepareUU(const char *method) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingUMesh *src_mesh=(MEDCouplingUMesh *)_src_mesh;
  MEDCouplingUMesh *target_mesh=(MEDCouplingUMesh *)_target_mesh;
  INTERP_KERNEL::Interpolation<INTERP_KERNEL::Interpolation3D>::checkAndSplitInterpolationMethod(method,_src_method,_target_method);
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
      INTERP_KERNEL::Interpolation3D interpolation(*this);
      nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
    }
  else if(srcMeshDim==1 && trgMeshDim==3 && srcSpaceDim==3)
    {
      if(getIntersectionType()!=INTERP_KERNEL::PointLocator)
        throw INTERP_KERNEL::Exception("Invalid interpolation requested between 3D and 1D ! Select PointLocator as intersection type !");
      MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(src_mesh);
      MEDCouplingNormalizedUnstructuredMesh<3,3> target_mesh_wrapper(target_mesh);
      INTERP_KERNEL::Interpolation3D interpolation(*this);
      std::vector<std::map<int,double> > matrixTmp;
      nbCols=interpolation.interpolateMeshes(target_mesh_wrapper,source_mesh_wrapper,matrixTmp,method);
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
          nbCols=interpolation.interpolateMeshes(target_mesh_wrapper,source_mesh_wrapper,matrixTmp,method);
          ReverseMatrix(matrixTmp,nbCols,_matrix);
          nbCols=matrixTmp.size();
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
  else if(srcMeshDim==1 && trgMeshDim==2 && srcSpaceDim==2)
    {
      if(getIntersectionType()==INTERP_KERNEL::PointLocator)
        {
          MEDCouplingNormalizedUnstructuredMesh<2,2> source_mesh_wrapper(src_mesh);
          MEDCouplingNormalizedUnstructuredMesh<2,2> target_mesh_wrapper(target_mesh);
          INTERP_KERNEL::Interpolation2D interpolation(*this);
          std::vector<std::map<int,double> > matrixTmp;
          nbCols=interpolation.interpolateMeshes(target_mesh_wrapper,source_mesh_wrapper,matrixTmp,method);
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
      INTERP_KERNEL::Interpolation3D2D interpolation(*this);
      nbCols=interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,_matrix,method);
      INTERP_KERNEL::Interpolation3D2D::DuplicateFacesType duplicateFaces=interpolation.retrieveDuplicateFaces();
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
      MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(src_mesh);
      MEDCouplingNormalizedUnstructuredMesh<3,3> target_mesh_wrapper(target_mesh);
      INTERP_KERNEL::Interpolation3D2D interpolation(*this);
      std::vector<std::map<int,double> > matrixTmp;
      nbCols=interpolation.interpolateMeshes(target_mesh_wrapper,source_mesh_wrapper,matrixTmp,method);
      ReverseMatrix(matrixTmp,nbCols,_matrix);
      nbCols=matrixTmp.size();
      INTERP_KERNEL::Interpolation3D2D::DuplicateFacesType duplicateFaces=interpolation.retrieveDuplicateFaces();
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
  else if(trgMeshDim==-1)
    {
      if(srcMeshDim==2 && srcSpaceDim==2)
        {
          MEDCouplingNormalizedUnstructuredMesh<2,2> source_mesh_wrapper(src_mesh);
          INTERP_KERNEL::Interpolation2D interpolation(*this);
          nbCols=interpolation.toIntegralUniform(source_mesh_wrapper,_matrix,_src_method.c_str());
        }
      else if(srcMeshDim==3 && srcSpaceDim==3)
        {
          MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(src_mesh);
          INTERP_KERNEL::Interpolation3D interpolation(*this);
          nbCols=interpolation.toIntegralUniform(source_mesh_wrapper,_matrix,_src_method.c_str());
        }
      else if(srcMeshDim==2 && srcSpaceDim==3)
        {
          MEDCouplingNormalizedUnstructuredMesh<3,2> source_mesh_wrapper(src_mesh);
          INTERP_KERNEL::Interpolation3DSurf interpolation(*this);
          nbCols=interpolation.toIntegralUniform(source_mesh_wrapper,_matrix,_src_method.c_str());
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
          nbCols=interpolation.fromIntegralUniform(source_mesh_wrapper,_matrix,_target_method.c_str());
        }
      else if(trgMeshDim==3 && trgSpaceDim==3)
        {
          MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(target_mesh);
          INTERP_KERNEL::Interpolation3D interpolation(*this);
          nbCols=interpolation.fromIntegralUniform(source_mesh_wrapper,_matrix,_target_method.c_str());
        }
      else if(trgMeshDim==2 && trgSpaceDim==3)
        {
          MEDCouplingNormalizedUnstructuredMesh<3,2> source_mesh_wrapper(target_mesh);
          INTERP_KERNEL::Interpolation3DSurf interpolation(*this);
          nbCols=interpolation.fromIntegralUniform(source_mesh_wrapper,_matrix,_target_method.c_str());
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

int MEDCouplingRemapper::prepareEE(const char *method) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingExtrudedMesh *src_mesh=(MEDCouplingExtrudedMesh *)_src_mesh;
  MEDCouplingExtrudedMesh *target_mesh=(MEDCouplingExtrudedMesh *)_target_mesh;
  std::string methC(method);
  if(methC!="P0P0")
    throw INTERP_KERNEL::Exception("Only P0P0 method implemented for Extruded/Extruded meshes !");
  INTERP_KERNEL::Interpolation<INTERP_KERNEL::Interpolation3D>::checkAndSplitInterpolationMethod(method,_src_method,_target_method);
  MEDCouplingNormalizedUnstructuredMesh<3,2> source_mesh_wrapper(src_mesh->getMesh2D());
  MEDCouplingNormalizedUnstructuredMesh<3,2> target_mesh_wrapper(target_mesh->getMesh2D());
  INTERP_KERNEL::Interpolation3DSurf interpolation2D(*this);
  std::vector<std::map<int,double> > matrix2D;
  int nbCols2D=interpolation2D.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,matrix2D,method);
  MEDCouplingUMesh *s1D,*t1D;
  double v[3];
  MEDCouplingExtrudedMesh::Project1DMeshes(src_mesh->getMesh1D(),target_mesh->getMesh1D(),getPrecision(),s1D,t1D,v);
  MEDCouplingNormalizedUnstructuredMesh<1,1> s1DWrapper(s1D);
  MEDCouplingNormalizedUnstructuredMesh<1,1> t1DWrapper(t1D);
  std::vector<std::map<int,double> > matrix1D;
  INTERP_KERNEL::Interpolation1D interpolation1D(*this);
  int nbCols1D=interpolation1D.interpolateMeshes(s1DWrapper,t1DWrapper,matrix1D,method);
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

int MEDCouplingRemapper::prepareUC(const char *method) throw(INTERP_KERNEL::Exception)
{
  std::string methodCpp(method);
  if(methodCpp!="P0P0")
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareUC : only P0P0 interpolation supported for the moment !");
  INTERP_KERNEL::Interpolation<INTERP_KERNEL::Interpolation3D>::checkAndSplitInterpolationMethod(method,_src_method,_target_method);
  MEDCouplingUMesh *src_mesh=static_cast<MEDCouplingUMesh *>(_src_mesh);
  MEDCouplingCMesh *target_mesh=static_cast<MEDCouplingCMesh *>(_target_mesh);
  const int srcMeshDim=src_mesh->getMeshDimension();
  const int srcSpceDim=src_mesh->getSpaceDimension();
  const int trgMeshDim=target_mesh->getMeshDimension();
  if(srcMeshDim!=srcSpceDim || srcMeshDim!=trgMeshDim)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareUC : space dim of src unstructured should be equal to mesh dim of src unstructured and should be equal also equal to trg cartesian dimension !");
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
      throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareUC : only dimension 1 2 or 3 supported !");
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

int MEDCouplingRemapper::prepareCU(const char *method) throw(INTERP_KERNEL::Exception)
{
  std::string methodCpp(method);
  if(methodCpp!="P0P0")
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareCU : only P0P0 interpolation supported for the moment !");
  INTERP_KERNEL::Interpolation<INTERP_KERNEL::Interpolation3D>::checkAndSplitInterpolationMethod(method,_src_method,_target_method);
  MEDCouplingCMesh *src_mesh=static_cast<MEDCouplingCMesh *>(_src_mesh);
  MEDCouplingUMesh *target_mesh=static_cast<MEDCouplingUMesh *>(_target_mesh);
  const int srcMeshDim=src_mesh->getMeshDimension();
  const int trgMeshDim=target_mesh->getMeshDimension();
  const int trgSpceDim=target_mesh->getSpaceDimension();
  if(trgMeshDim!=trgSpceDim || trgMeshDim!=srcMeshDim)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareCU : space dim of target unstructured should be equal to mesh dim of target unstructured and should be equal also equal to source cartesian dimension !");
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
      throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareCU : only dimension 1 2 or 3 supported !");
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

int MEDCouplingRemapper::prepareCC(const char *method) throw(INTERP_KERNEL::Exception)
{
  std::string methodCpp(method);
  if(methodCpp!="P0P0")
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareCC : only P0P0 interpolation supported for the moment !");
  INTERP_KERNEL::Interpolation<INTERP_KERNEL::Interpolation3D>::checkAndSplitInterpolationMethod(method,_src_method,_target_method);
  MEDCouplingCMesh *src_mesh=static_cast<MEDCouplingCMesh *>(_src_mesh);
  MEDCouplingCMesh *target_mesh=static_cast<MEDCouplingCMesh *>(_target_mesh);
  const int srcMeshDim=src_mesh->getMeshDimension();
  const int trgMeshDim=target_mesh->getMeshDimension();
  if(trgMeshDim!=srcMeshDim)
    throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareCC : dim of target cartesian should be equal to dim of source cartesian dimension !");
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
      throw INTERP_KERNEL::Exception("MEDCouplingRemapper::prepareCC : only dimension 1 2 or 3 supported !");
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

void MEDCouplingRemapper::updateTime() const
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

void MEDCouplingRemapper::transferUnderground(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField, bool isDftVal, double dftValue) throw(INTERP_KERNEL::Exception)
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
      if(!isDftVal)
        throw INTERP_KERNEL::Exception("MEDCouplingRemapper::partialTransfer : This method requires that the array of target field exists ! Allocate it or call MEDCouplingRemapper::transfer instead !");
      array=DataArrayDouble::New();
      array->alloc(targetField->getNumberOfTuples(),srcNbOfCompo);
      targetField->setArray(array);
      array->decrRef();
    }
  computeDeno(srcField->getNature(),srcField,targetField);
  double *resPointer=array->getPointer();
  const double *inputPointer=srcField->getArray()->getConstPointer();
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

void MEDCouplingRemapper::computeDenoFromScratch(NatureOfField nat, const MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *trgField) throw(INTERP_KERNEL::Exception)
{
  _nature_of_deno=nat;
  _time_deno_update=getTimeOfThis();
  switch(_nature_of_deno)
    {
    case ConservativeVolumic:
      {
        ComputeRowSumAndColSum(_matrix,_deno_multiply,_deno_reverse_multiply);
        break;
      }
    case Integral:
      {
        MEDCouplingFieldDouble *deno=srcField->getDiscretization()->getMeasureField(srcField->getMesh(),true);
        MEDCouplingFieldDouble *denoR=trgField->getDiscretization()->getMeasureField(trgField->getMesh(),true);
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
    case IntegralGlobConstraint:
      {
        ComputeColSumAndRowSum(_matrix,_deno_multiply,_deno_reverse_multiply);
        break;
      }
    case RevIntegral:
      {
        MEDCouplingFieldDouble *deno=trgField->getDiscretization()->getMeasureField(trgField->getMesh(),true);
        MEDCouplingFieldDouble *denoR=srcField->getDiscretization()->getMeasureField(srcField->getMesh(),true);
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
int MEDCouplingRemapper::nullifiedTinyCoeffInCrudeMatrixAbs(double maxValAbs) throw(INTERP_KERNEL::Exception)
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
int MEDCouplingRemapper::nullifiedTinyCoeffInCrudeMatrix(double scaleFactor) throw(INTERP_KERNEL::Exception)
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
double MEDCouplingRemapper::getMaxValueInCrudeMatrix() const throw(INTERP_KERNEL::Exception)
{
  double ret=0.;
  for(std::vector<std::map<int,double> >::const_iterator it1=_matrix.begin();it1!=_matrix.end();it1++)
    for(std::map<int,double>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
      if(fabs((*it2).second)>ret)
        ret=fabs((*it2).second);
  return ret;
}
