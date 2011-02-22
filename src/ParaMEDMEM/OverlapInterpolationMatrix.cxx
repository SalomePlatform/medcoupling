//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

  OverlapInterpolationMatrix::~OverlapInterpolationMatrix()
  {
  }

  void OverlapInterpolationMatrix::addContribution(const MEDCouplingPointSet *src, const DataArrayInt *srcIds, const std::string& srcMeth, int srcProcId,
                                                   const MEDCouplingPointSet *trg, const DataArrayInt *trgIds, const std::string& trgMeth, int trgProcId)
  {
    std::string interpMethod(trgMeth);
    interpMethod+=srcMeth;
    //creating the interpolator structure
    vector<map<int,double> > surfaces;
    int colSize=0;
    //computation of the intersection volumes between source and target elements
    const MEDCouplingUMesh *trgC=dynamic_cast<const MEDCouplingUMesh *>(trg);
    const MEDCouplingUMesh *srcC=dynamic_cast<const MEDCouplingUMesh *>(src);
    if ( trg->getMeshDimension() == -1 )
      {
        if(srcC->getMeshDimension()==2 && srcC->getSpaceDimension()==2)
          {
            MEDCouplingNormalizedUnstructuredMesh<2,2> source_mesh_wrapper(srcC);
            INTERP_KERNEL::Interpolation2D interpolation(*this);
            colSize=interpolation.fromIntegralUniform(source_mesh_wrapper,surfaces,srcMeth.c_str());
          }
        else if(srcC->getMeshDimension()==3 && srcC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(srcC);
            INTERP_KERNEL::Interpolation3D interpolation(*this);
            colSize=interpolation.fromIntegralUniform(source_mesh_wrapper,surfaces,srcMeth.c_str());
          }
        else if(srcC->getMeshDimension()==2 && srcC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,2> source_mesh_wrapper(srcC);
            INTERP_KERNEL::Interpolation3DSurf interpolation(*this);
            colSize=interpolation.fromIntegralUniform(source_mesh_wrapper,surfaces,srcMeth.c_str());
          }
        else
          throw INTERP_KERNEL::Exception("No para interpolation available for the given mesh and space dimension of source mesh to -1D targetMesh");
      }
    else if ( srcC->getMeshDimension() == -1 )
      {
        if(trgC->getMeshDimension()==2 && trgC->getSpaceDimension()==2)
          {
            MEDCouplingNormalizedUnstructuredMesh<2,2> distant_mesh_wrapper(trgC);
            INTERP_KERNEL::Interpolation2D interpolation(*this);
            colSize=interpolation.toIntegralUniform(distant_mesh_wrapper,surfaces,srcMeth.c_str());
          }
        else if(trgC->getMeshDimension()==3 && trgC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,3> distant_mesh_wrapper(trgC);
            INTERP_KERNEL::Interpolation3D interpolation(*this);
            colSize=interpolation.toIntegralUniform(distant_mesh_wrapper,surfaces,srcMeth.c_str());
          }
        else if(trgC->getMeshDimension()==2 && trgC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,2> distant_mesh_wrapper(trgC);
            INTERP_KERNEL::Interpolation3DSurf interpolation(*this);
            colSize=interpolation.toIntegralUniform(distant_mesh_wrapper,surfaces,srcMeth.c_str());
          }
        else
          throw INTERP_KERNEL::Exception("No para interpolation available for the given mesh and space dimension of distant mesh to -1D sourceMesh");
      }
    else if (trg->getMeshDimension() != _source_support->getMeshDimension())
      {
        throw INTERP_KERNEL::Exception("local and distant meshes do not have the same space and mesh dimensions");
      }
    else if( trg->getMeshDimension() == 1
             && trg->getSpaceDimension() == 1 )
      {
        MEDCouplingNormalizedUnstructuredMesh<1,1> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<1,1> source_wrapper(srcC);

        INTERP_KERNEL::Interpolation1D interpolation(*this);
        colSize=interpolation.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod.c_str());
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if( trg->getMeshDimension() == 1
             && trg->getSpaceDimension() == 2 )
      {
        MEDCouplingNormalizedUnstructuredMesh<2,1> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<2,1> source_wrapper(srcC);

        INTERP_KERNEL::Interpolation2DCurve interpolation(*this);
        colSize=interpolation.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod.c_str());
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if ( trg->getMeshDimension() == 2
              && trg->getSpaceDimension() == 3 )
      {
        MEDCouplingNormalizedUnstructuredMesh<3,2> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<3,2> source_wrapper(srcC);

        INTERP_KERNEL::Interpolation3DSurf interpolator (*this);
        colSize=interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod.c_str());
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if ( trg->getMeshDimension() == 2
              && trg->getSpaceDimension() == 2)
      {
        MEDCouplingNormalizedUnstructuredMesh<2,2> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<2,2> source_wrapper(srcC);

        INTERP_KERNEL::Interpolation2D interpolator (*this);
        colSize=interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod.c_str());
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if ( trg->getMeshDimension() == 3
              && trg->getSpaceDimension() == 3 )
      {
        MEDCouplingNormalizedUnstructuredMesh<3,3> target_wrapper(trgC);
        MEDCouplingNormalizedUnstructuredMesh<3,3> source_wrapper(srcC);

        INTERP_KERNEL::Interpolation3D interpolator (*this);
        colSize=interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod.c_str());
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else
      {
        throw INTERP_KERNEL::Exception("no interpolator exists for these mesh and space dimensions ");
      }
    bool needTargetSurf=isSurfaceComputationNeeded(trgMeth);
    MEDCouplingFieldDouble *target_triangle_surf=0;
    if(needTargetSurf)
      target_triangle_surf=trg->getMeasureField(getMeasureAbsStatus());
    //
    fillDistributedMatrix(surfaces,srcIds,srcProcId,trgIds,trgProcId);
    //
    if(needTargetSurf)
      target_triangle_surf->decrRef();
  }

  /*!
   * \b WARNING : res rows refers to source and column (first param of map) to target.
   */
  void OverlapInterpolationMatrix::fillDistributedMatrix(const std::vector< std::map<int,double> >& res,
                                                         const DataArrayInt *srcIds, int srcProc,
                                                         const DataArrayInt *trgIds, int trgProc)
  {
    //computing matrix with real ids for target
    int sz=res.size();
    std::vector< std::map<int,double> > res1(sz);
    const int *trgIds2=0;
    int nbTrgIds=_target_field->getField()->getNumberOfTuplesExpected();
    INTERP_KERNEL::AutoPtr<int> tmp2=new int[nbTrgIds];
    if(trgIds)
      trgIds2=trgIds->getConstPointer();
    else
      {
        trgIds2=tmp2;
        for(int i=0;i<nbTrgIds;i++)
          tmp2[i]=i;
      }
    for(int i=0;i<sz;i++)
      {
        std::map<int,double>& m=res1[i];
        const std::map<int,double>& ref=res[i];
        for(std::map<int,double>::const_iterator it=ref.begin();it!=ref.end();it++)
          {
            m[(*it).first]=(*it).second; //if(trgIds2) m[trgIds2[(*it).first]]=(*it).second;
          } 
      }
    //dealing source ids
    if(srcIds)
      _mapping.addContributionST(res1,srcIds->getConstPointer(),trgIds2,nbTrgIds,srcProc,trgProc);
    else
      {
        INTERP_KERNEL::AutoPtr<int> tmp=new int[sz];
        for(int i=0;i<sz;i++)
          tmp[i]=i;
        _mapping.addContributionST(res1,tmp,trgIds2,nbTrgIds,srcProc,trgProc);
      }
  }

  /*!
   * 'procsInInteraction' gives the global view of interaction between procs.
   * In 'procsInInteraction' for a proc with id i, is in interaction with procs listed in procsInInteraction[i]
   */
  void OverlapInterpolationMatrix::prepare(const std::vector< std::vector<int> >& procsInInteraction)
  {
    if(_source_support)
      _mapping.prepare(procsInInteraction,_source_field->getField()->getNumberOfTuplesExpected());
    else
      _mapping.prepare(procsInInteraction,0);
  }

  void OverlapInterpolationMatrix::computeDeno()
  {
    if(_source_field->getField()->getNature()==IntegralGlobConstraint)
      _mapping.computeDenoGlobConstraint();
  }

  void OverlapInterpolationMatrix::transposeMultiply()
  {
    _mapping.transposeMultiply(_target_field->getField(),_source_field->getField());
  }
  
  bool OverlapInterpolationMatrix::isSurfaceComputationNeeded(const std::string& method) const
  {
    return method=="P0";
  }
}
