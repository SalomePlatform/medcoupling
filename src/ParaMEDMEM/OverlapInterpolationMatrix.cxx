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

#include "OverlapInterpolationMatrix.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "ProcessorGroup.hxx"
#include "TranslationRotationMatrix.hxx"
#include "Interpolation.hxx"
#include "Interpolation1D.txx"
#include "Interpolation2DCurve.hxx"
#include "Interpolation2D.txx"
#include "Interpolation3DSurf.hxx"
#include "Interpolation3D.txx"
#include "Interpolation2D3D.txx"
#include "Interpolation2D1D.txx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingNormalizedUnstructuredMesh.txx"
#include "InterpolationOptions.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "ElementLocator.hxx"
#include "InterpKernelAutoPtr.hxx"

#include <algorithm>

using namespace std;

namespace MEDCoupling
{
  OverlapInterpolationMatrix::OverlapInterpolationMatrix(ParaFIELD *source_field,
                                                         ParaFIELD *target_field,
                                                         const ProcessorGroup& group,
                                                         const DECOptions& dec_options,
                                                         const INTERP_KERNEL::InterpolationOptions& i_opt,
                                                         const OverlapElementLocator & locator):
    INTERP_KERNEL::InterpolationOptions(i_opt),
    DECOptions(dec_options),
    _source_field(source_field),
    _target_field(target_field),
    _source_support(source_field->getSupport()->getCellMesh()),
    _target_support(target_field->getSupport()->getCellMesh()),
    _mapping(group, locator),
    _group(group)
  {
  }

  void OverlapInterpolationMatrix::keepTracksOfSourceIds(int procId, DataArrayInt *ids)
  {
    _mapping.keepTracksOfSourceIds(procId,ids);
  }

  void OverlapInterpolationMatrix::keepTracksOfTargetIds(int procId, DataArrayInt *ids)
  {
    _mapping.keepTracksOfTargetIds(procId,ids);
  }

  OverlapInterpolationMatrix::~OverlapInterpolationMatrix()
  {
  }

  // TODO? Merge with MEDCouplingRemapper::prepareInterpKernelOnlyUU() ?
  /**!
   * Local run (on this proc) of the sequential interpolation algorithm.
   *
   * @param srcIds is null if the source mesh is on the local proc
   * @param trgIds is null if the source mesh is on the local proc
   *
   * One of the 2 is necessarily null (the two can be null together)
   */
  void OverlapInterpolationMatrix::computeLocalIntersection(const MEDCouplingPointSet *src, const DataArrayInt *srcIds, const std::string& srcMeth, int srcProcId,
                                                   const MEDCouplingPointSet *trg, const DataArrayInt *trgIds, const std::string& trgMeth, int trgProcId)
  {
    std::string interpMethod(srcMeth);
    interpMethod+=trgMeth;
    //creating the interpolator structure
    vector<SparseDoubleVec > sparse_matrix_part;
    int colSize=0;
    //computation of the intersection volumes between source and target elements
    const MEDCouplingUMesh *trgC=dynamic_cast<const MEDCouplingUMesh *>(trg);
    const MEDCouplingUMesh *srcC=dynamic_cast<const MEDCouplingUMesh *>(src);
    if ( src->getMeshDimension() == -1 )
      {
        if(trgC->getMeshDimension()==2 && trgC->getSpaceDimension()==2)
          {
            MEDCouplingNormalizedUnstructuredMesh<2,2> target_mesh_wrapper(trgC);
            INTERP_KERNEL::Interpolation2D interpolation(*this);
            colSize=interpolation.fromIntegralUniform(target_mesh_wrapper,sparse_matrix_part,trgMeth);
          }
        else if(trgC->getMeshDimension()==3 && trgC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,3> target_mesh_wrapper(trgC);
            INTERP_KERNEL::Interpolation3D interpolation(*this);
            colSize=interpolation.fromIntegralUniform(target_mesh_wrapper,sparse_matrix_part,trgMeth);
          }
        else if(trgC->getMeshDimension()==2 && trgC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,2> target_mesh_wrapper(trgC);
            INTERP_KERNEL::Interpolation3DSurf interpolation(*this);
            colSize=interpolation.fromIntegralUniform(target_mesh_wrapper,sparse_matrix_part,trgMeth);
          }
        else
          throw INTERP_KERNEL::Exception("No para interpolation available for the given mesh and space dimension of source mesh to -1D targetMesh");
      }
    else if ( trg->getMeshDimension() == -1 )
      {
        if(srcC->getMeshDimension()==2 && srcC->getSpaceDimension()==2)
          {
            MEDCouplingNormalizedUnstructuredMesh<2,2> local_mesh_wrapper(srcC);
            INTERP_KERNEL::Interpolation2D interpolation(*this);
            colSize=interpolation.toIntegralUniform(local_mesh_wrapper,sparse_matrix_part,srcMeth);
          }
        else if(srcC->getMeshDimension()==3 && srcC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,3> local_mesh_wrapper(srcC);
            INTERP_KERNEL::Interpolation3D interpolation(*this);
            colSize=interpolation.toIntegralUniform(local_mesh_wrapper,sparse_matrix_part,srcMeth);
          }
        else if(srcC->getMeshDimension()==2 && srcC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,2> local_mesh_wrapper(srcC);
            INTERP_KERNEL::Interpolation3DSurf interpolation(*this);
            colSize=interpolation.toIntegralUniform(local_mesh_wrapper,sparse_matrix_part,srcMeth);
          }
        else
          throw INTERP_KERNEL::Exception("No para interpolation available for the given mesh and space dimension of distant mesh to -1D sourceMesh");
      }
    else if ( src->getMeshDimension() == 2 && trg->getMeshDimension() == 3
              && trg->getSpaceDimension() == 3 && src->getSpaceDimension() == 3 )
      {
        MEDCouplingNormalizedUnstructuredMesh<3,3> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<3,3> source_wrapper(srcC);
        
        INTERP_KERNEL::Interpolation2D3D interpolator (*this);
        colSize=interpolator.interpolateMeshes(source_wrapper,target_wrapper,sparse_matrix_part,interpMethod);
      }
    else if ( src->getMeshDimension() == 3 && trg->getMeshDimension() == 2
              && trg->getSpaceDimension() == 3 && src->getSpaceDimension() == 3 )
      {
        MEDCouplingNormalizedUnstructuredMesh<3,3> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<3,3> source_wrapper(srcC);
        
        INTERP_KERNEL::Interpolation2D3D interpolator (*this);
        vector<SparseDoubleVec > matrixTranspose;
        colSize=interpolator.interpolateMeshes(target_wrapper,source_wrapper,sparse_matrix_part,interpMethod);//not a bug target in source.
        TransposeMatrix(matrixTranspose,colSize,sparse_matrix_part);
        colSize=matrixTranspose.size();
      }
    else if ( src->getMeshDimension() == 1 && trg->getMeshDimension() == 2
              && trg->getSpaceDimension() == 2 && src->getSpaceDimension() == 2 )
      {
        MEDCouplingNormalizedUnstructuredMesh<2,2> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<2,2> source_wrapper(srcC);
        
        INTERP_KERNEL::Interpolation2D1D interpolator (*this);
        colSize=interpolator.interpolateMeshes(source_wrapper,target_wrapper,sparse_matrix_part,interpMethod);
      }
    else if ( src->getMeshDimension() == 2 && trg->getMeshDimension() == 1
              && trg->getSpaceDimension() == 2 && src->getSpaceDimension() == 2 )
      {
        MEDCouplingNormalizedUnstructuredMesh<2,2> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<2,2> source_wrapper(srcC);
        
        INTERP_KERNEL::Interpolation2D1D interpolator (*this);
        vector<SparseDoubleVec > matrixTranspose;
        colSize=interpolator.interpolateMeshes(target_wrapper,source_wrapper,matrixTranspose,interpMethod);//not a bug target in source.
        TransposeMatrix(matrixTranspose,colSize,sparse_matrix_part);
        colSize=matrixTranspose.size();
      }
    else if (trg->getMeshDimension() != _source_support->getMeshDimension())
      {
        throw INTERP_KERNEL::Exception("local and distant meshes do not have the same space and mesh dimensions");
      }
    else if( src->getMeshDimension() == 1
             && src->getSpaceDimension() == 1 )
      {
        MEDCouplingNormalizedUnstructuredMesh<1,1> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<1,1> source_wrapper(srcC);

        INTERP_KERNEL::Interpolation1D interpolation(*this);
        colSize=interpolation.interpolateMeshes(source_wrapper,target_wrapper,sparse_matrix_part,interpMethod);
      }
    else if( trg->getMeshDimension() == 1
             && trg->getSpaceDimension() == 2 )
      {
        MEDCouplingNormalizedUnstructuredMesh<2,1> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<2,1> source_wrapper(srcC);

        INTERP_KERNEL::Interpolation2DCurve interpolation(*this);
        colSize=interpolation.interpolateMeshes(source_wrapper,target_wrapper,sparse_matrix_part,interpMethod);
      }
    else if ( trg->getMeshDimension() == 2
              && trg->getSpaceDimension() == 3 )
      {
        MEDCouplingNormalizedUnstructuredMesh<3,2> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<3,2> source_wrapper(srcC);

        INTERP_KERNEL::Interpolation3DSurf interpolator (*this);
        colSize=interpolator.interpolateMeshes(source_wrapper,target_wrapper,sparse_matrix_part,interpMethod);
      }
    else if ( trg->getMeshDimension() == 2
              && trg->getSpaceDimension() == 2)
      {
        MEDCouplingNormalizedUnstructuredMesh<2,2> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<2,2> source_wrapper(srcC);

        INTERP_KERNEL::Interpolation2D interpolator (*this);
        colSize=interpolator.interpolateMeshes(source_wrapper,target_wrapper,sparse_matrix_part,interpMethod);
      }
    else if ( trg->getMeshDimension() == 3
              && trg->getSpaceDimension() == 3 )
      {
        MEDCouplingNormalizedUnstructuredMesh<3,3> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<3,3> source_wrapper(srcC);

        INTERP_KERNEL::Interpolation3D interpolator (*this);
        colSize=interpolator.interpolateMeshes(source_wrapper,target_wrapper,sparse_matrix_part,interpMethod);
      }
    else
      {
        throw INTERP_KERNEL::Exception("No interpolator exists for these mesh and space dimensions!");
      }
    /* Fill distributed matrix:
       In sparse_matrix_part rows refer to target, and columns (=first param of map in SparseDoubleVec)
       refer to source.
     */
    _mapping.addContributionST(sparse_matrix_part,srcIds,srcProcId,trgIds,trgProcId);
  }

  /*!
   * 'procsToSendField' gives the list of procs field data has to be sent to.
   * See OverlapElementLocator::computeBoundingBoxesAndTodoList()
   */
  void OverlapInterpolationMatrix::prepare(const std::vector< int >& procsToSendField)
  {
    if(_source_support)
      _mapping.prepare(procsToSendField,_target_field->getField()->getNumberOfTuplesExpected());
    else
      _mapping.prepare(procsToSendField,0);
  }

  void OverlapInterpolationMatrix::computeSurfacesAndDeno()
  {
    if(_target_field->getField()->getNature()==IntensiveMaximum)
      _mapping.computeDenoConservativeVolumic(_target_field->getField()->getNumberOfTuplesExpected());
    else
      throw INTERP_KERNEL::Exception("OverlapDEC: Policy not implemented yet: only IntensiveMaximum!");
//      {
//      if(_target_field->getField()->getNature()==IntensiveConservation)
//        {
//          MCAuto<MEDCouplingFieldDouble> f;
//          int orient = getOrientation(); // From InterpolationOptions inheritance
//          if(orient == 2)  // absolute areas
//             f = _target_support->getMeasureField(true);
//          else
//            if(orient == 0)  // relative areas
//              f = _target_support->getMeasureField(false);
//            else
//              throw INTERP_KERNEL::Exception("OverlapDEC: orientation policy not impl. yet!");
//          _mapping.computeDenoRevIntegral(*f->getArray());
//        }
//      else
//        throw INTERP_KERNEL::Exception("OverlapDEC: Policy not implemented yet: only IntensiveMaximum and IntensiveConservation defined!");
//      }
  }

  void OverlapInterpolationMatrix::multiply(double default_val)
  {
    _mapping.multiply(_source_field->getField(),_target_field->getField(), default_val);
  }

  void OverlapInterpolationMatrix::transposeMultiply()
  {
    _mapping.transposeMultiply(_target_field->getField(),_source_field->getField());
  }
  
  void OverlapInterpolationMatrix::TransposeMatrix(const std::vector<SparseDoubleVec >& matIn,
                                                   int nbColsMatIn, std::vector<SparseDoubleVec >& matOut)
  {
    matOut.resize(nbColsMatIn);
    int id=0;
    for(std::vector<SparseDoubleVec >::const_iterator iter1=matIn.begin();iter1!=matIn.end();iter1++,id++)
      for(SparseDoubleVec::const_iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
        matOut[(*iter2).first][id]=(*iter2).second;
  }
}
