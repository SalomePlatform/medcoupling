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
#include "Interpolation3D2D.txx"
#include "Interpolation2D1D.txx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingNormalizedUnstructuredMesh.txx"
#include "InterpolationOptions.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "ElementLocator.hxx"
#include "InterpKernelAutoPtr.hxx"

#include <algorithm>

using namespace std;

namespace ParaMEDMEM
{
  OverlapInterpolationMatrix::OverlapInterpolationMatrix(ParaFIELD *source_field,
                                                         ParaFIELD *target_field,
                                                         const ProcessorGroup& group,
                                                         const DECOptions& dec_options,
                                                         const INTERP_KERNEL::InterpolationOptions& i_opt):
    INTERP_KERNEL::InterpolationOptions(i_opt),
    DECOptions(dec_options),
    _source_field(source_field),
    _target_field(target_field),
    _source_support(source_field->getSupport()->getCellMesh()),
    _target_support(target_field->getSupport()->getCellMesh()),
    _mapping(group),
    _group(group)
  {
    int nbelems = source_field->getField()->getNumberOfTuples();
    _row_offsets.resize(nbelems+1);
    _coeffs.resize(nbelems);
    _target_volume.resize(nbelems);
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

  void OverlapInterpolationMatrix::addContribution(const MEDCouplingPointSet *src, const DataArrayInt *srcIds, const std::string& srcMeth, int srcProcId,
                                                   const MEDCouplingPointSet *trg, const DataArrayInt *trgIds, const std::string& trgMeth, int trgProcId)
  {
    std::string interpMethod(srcMeth);
    interpMethod+=trgMeth;
    //creating the interpolator structure
    vector<map<int,double> > surfaces;
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
            colSize=interpolation.fromIntegralUniform(target_mesh_wrapper,surfaces,trgMeth);
          }
        else if(trgC->getMeshDimension()==3 && trgC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,3> target_mesh_wrapper(trgC);
            INTERP_KERNEL::Interpolation3D interpolation(*this);
            colSize=interpolation.fromIntegralUniform(target_mesh_wrapper,surfaces,trgMeth);
          }
        else if(trgC->getMeshDimension()==2 && trgC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,2> target_mesh_wrapper(trgC);
            INTERP_KERNEL::Interpolation3DSurf interpolation(*this);
            colSize=interpolation.fromIntegralUniform(target_mesh_wrapper,surfaces,trgMeth);
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
            colSize=interpolation.toIntegralUniform(local_mesh_wrapper,surfaces,srcMeth);
          }
        else if(srcC->getMeshDimension()==3 && srcC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,3> local_mesh_wrapper(srcC);
            INTERP_KERNEL::Interpolation3D interpolation(*this);
            colSize=interpolation.toIntegralUniform(local_mesh_wrapper,surfaces,srcMeth);
          }
        else if(srcC->getMeshDimension()==2 && srcC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,2> local_mesh_wrapper(srcC);
            INTERP_KERNEL::Interpolation3DSurf interpolation(*this);
            colSize=interpolation.toIntegralUniform(local_mesh_wrapper,surfaces,srcMeth);
          }
        else
          throw INTERP_KERNEL::Exception("No para interpolation available for the given mesh and space dimension of distant mesh to -1D sourceMesh");
      }
    else if ( src->getMeshDimension() == 2 && trg->getMeshDimension() == 3
              && trg->getSpaceDimension() == 3 && src->getSpaceDimension() == 3 )
      {
        MEDCouplingNormalizedUnstructuredMesh<3,3> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<3,3> source_wrapper(srcC);
        
        INTERP_KERNEL::Interpolation3D2D interpolator (*this);
        colSize=interpolator.interpolateMeshes(source_wrapper,target_wrapper,surfaces,interpMethod);
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if ( src->getMeshDimension() == 3 && trg->getMeshDimension() == 2
              && trg->getSpaceDimension() == 3 && src->getSpaceDimension() == 3 )
      {
        MEDCouplingNormalizedUnstructuredMesh<3,3> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<3,3> source_wrapper(srcC);
        
        INTERP_KERNEL::Interpolation3D2D interpolator (*this);
        vector<map<int,double> > surfacesTranspose;
        colSize=interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod);//not a bug target in source.
        TransposeMatrix(surfacesTranspose,colSize,surfaces);
        colSize=surfacesTranspose.size();
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if ( src->getMeshDimension() == 1 && trg->getMeshDimension() == 2
              && trg->getSpaceDimension() == 2 && src->getSpaceDimension() == 2 )
      {
        MEDCouplingNormalizedUnstructuredMesh<2,2> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<2,2> source_wrapper(srcC);
        
        INTERP_KERNEL::Interpolation2D1D interpolator (*this);
        colSize=interpolator.interpolateMeshes(source_wrapper,target_wrapper,surfaces,interpMethod);
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if ( src->getMeshDimension() == 2 && trg->getMeshDimension() == 1
              && trg->getSpaceDimension() == 2 && src->getSpaceDimension() == 2 )
      {
        MEDCouplingNormalizedUnstructuredMesh<2,2> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<2,2> source_wrapper(srcC);
        
        INTERP_KERNEL::Interpolation2D1D interpolator (*this);
        vector<map<int,double> > surfacesTranspose;
        colSize=interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfacesTranspose,interpMethod);//not a bug target in source.
        TransposeMatrix(surfacesTranspose,colSize,surfaces);
        colSize=surfacesTranspose.size();
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
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
        colSize=interpolation.interpolateMeshes(source_wrapper,target_wrapper,surfaces,interpMethod);
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if( trg->getMeshDimension() == 1
             && trg->getSpaceDimension() == 2 )
      {
        MEDCouplingNormalizedUnstructuredMesh<2,1> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<2,1> source_wrapper(srcC);

        INTERP_KERNEL::Interpolation2DCurve interpolation(*this);
        colSize=interpolation.interpolateMeshes(source_wrapper,target_wrapper,surfaces,interpMethod);
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if ( trg->getMeshDimension() == 2
              && trg->getSpaceDimension() == 3 )
      {
        MEDCouplingNormalizedUnstructuredMesh<3,2> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<3,2> source_wrapper(srcC);

        INTERP_KERNEL::Interpolation3DSurf interpolator (*this);
        colSize=interpolator.interpolateMeshes(source_wrapper,target_wrapper,surfaces,interpMethod);
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if ( trg->getMeshDimension() == 2
              && trg->getSpaceDimension() == 2)
      {
        MEDCouplingNormalizedUnstructuredMesh<2,2> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<2,2> source_wrapper(srcC);

        INTERP_KERNEL::Interpolation2D interpolator (*this);
        colSize=interpolator.interpolateMeshes(source_wrapper,target_wrapper,surfaces,interpMethod);
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if ( trg->getMeshDimension() == 3
              && trg->getSpaceDimension() == 3 )
      {
        MEDCouplingNormalizedUnstructuredMesh<3,3> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<3,3> source_wrapper(srcC);

        INTERP_KERNEL::Interpolation3D interpolator (*this);
        colSize=interpolator.interpolateMeshes(source_wrapper,target_wrapper,surfaces,interpMethod);
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else
      {
        throw INTERP_KERNEL::Exception("no interpolator exists for these mesh and space dimensions ");
      }
    bool needSourceSurf=isSurfaceComputationNeeded(srcMeth);
    MEDCouplingFieldDouble *source_triangle_surf=0;
    if(needSourceSurf)
      source_triangle_surf=src->getMeasureField(getMeasureAbsStatus());
    //
    fillDistributedMatrix(surfaces,srcIds,srcProcId,trgIds,trgProcId);
    //
    if(needSourceSurf)
      source_triangle_surf->decrRef();
  }

  /*!
   * \b res rows refers to target and column (first param of map) to source.
   */
  void OverlapInterpolationMatrix::fillDistributedMatrix(const std::vector< std::map<int,double> >& res,
                                                         const DataArrayInt *srcIds, int srcProc,
                                                         const DataArrayInt *trgIds, int trgProc)
  {
    _mapping.addContributionST(res,srcIds,srcProc,trgIds,trgProc);
  }

  /*!
   * 'procsInInteraction' gives the global view of interaction between procs.
   * In 'procsInInteraction' for a proc with id i, is in interaction with procs listed in procsInInteraction[i]
   */
  void OverlapInterpolationMatrix::prepare(const std::vector< std::vector<int> >& procsInInteraction)
  {
    if(_source_support)
      _mapping.prepare(procsInInteraction,_target_field->getField()->getNumberOfTuplesExpected());
    else
      _mapping.prepare(procsInInteraction,0);
  }

  void OverlapInterpolationMatrix::computeDeno()
  {
    if(_target_field->getField()->getNature()==ConservativeVolumic)
      _mapping.computeDenoConservativeVolumic(_target_field->getField()->getNumberOfTuplesExpected());
    else
      throw INTERP_KERNEL::Exception("Policy Not implemented yet : only ConservativeVolumic defined !");
  }

  void OverlapInterpolationMatrix::multiply()
  {
    _mapping.multiply(_source_field->getField(),_target_field->getField());
  }

  void OverlapInterpolationMatrix::transposeMultiply()
  {
    _mapping.transposeMultiply(_target_field->getField(),_source_field->getField());
  }
  
  bool OverlapInterpolationMatrix::isSurfaceComputationNeeded(const std::string& method) const
  {
    return method=="P0";
  }

  void OverlapInterpolationMatrix::TransposeMatrix(const std::vector<std::map<int,double> >& matIn, int nbColsMatIn, std::vector<std::map<int,double> >& matOut)
  {
    matOut.resize(nbColsMatIn);
    int id=0;
    for(std::vector<std::map<int,double> >::const_iterator iter1=matIn.begin();iter1!=matIn.end();iter1++,id++)
      for(std::map<int,double>::const_iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
        matOut[(*iter2).first][id]=(*iter2).second;
  }
}
