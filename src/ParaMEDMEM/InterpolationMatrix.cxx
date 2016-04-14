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

#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "ProcessorGroup.hxx"
#include "MxN_Mapping.hxx"
#include "InterpolationMatrix.hxx"
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

#include <algorithm>

using namespace std;

namespace MEDCoupling
{

  /**!
     Creates an empty matrix structure linking two distributed supports.
     The method must be called by all processors belonging to source
     and target groups.
     \param source_support local support
     \param source_group processor group containing the local processors
     \param target_group processor group containing the distant processors
     \param method interpolation method
  */
  InterpolationMatrix::InterpolationMatrix(const MEDCoupling::ParaFIELD *source_field, 
                                           const ProcessorGroup& source_group,
                                           const ProcessorGroup& target_group,
                                           const DECOptions& dec_options,
                                           const INTERP_KERNEL::InterpolationOptions& interp_options):
    INTERP_KERNEL::InterpolationOptions(interp_options),
    DECOptions(dec_options),
    _source_field(source_field),
    _source_support(source_field->getSupport()->getCellMesh()),
    _mapping(source_group, target_group, dec_options),
    _source_group(source_group),
    _target_group(target_group)
  {
    int nbelems = source_field->getField()->getNumberOfTuples();
    _row_offsets.resize(nbelems+1);
    _coeffs.resize(nbelems);
    _target_volume.resize(nbelems);
  }

  InterpolationMatrix::~InterpolationMatrix()
  {
  }


  /*!
     \brief Adds the contribution of a distant subdomain to the*
     interpolation matrix.
     The method adds contribution to the interpolation matrix.
     For each row of the matrix, elements are addded as
     a (column, coeff) pair in the _coeffs array. This column number refers
     to an element on the target side via the _col_offsets array.
     It is made of a series of (iproc, ielem) pairs.
     The number of elements per row is stored in the row_offsets array.

     param distant_support local representation of the distant subdomain
     param iproc_distant id of the distant subdomain (in the distant group)
     param distant_elems mapping between the local representation of
     the subdomain and the actual elem ids on the distant subdomain
   */
  void InterpolationMatrix::addContribution ( MEDCouplingPointSet& distant_support,
                                              int iproc_distant,
                                              const int* distant_elems,
                                              const std::string& srcMeth,
                                              const std::string& targetMeth)
  {
    std::string interpMethod(targetMeth);
    interpMethod+=srcMeth;
    //creating the interpolator structure
    vector<map<int,double> > surfaces;
    //computation of the intersection volumes between source and target elements
    MEDCouplingUMesh *distant_supportC=dynamic_cast<MEDCouplingUMesh *>(&distant_support);
    MEDCouplingUMesh *source_supportC=dynamic_cast<MEDCouplingUMesh *>(_source_support);
    if ( distant_support.getMeshDimension() == -1 )
      {
        if(source_supportC->getMeshDimension()==2 && source_supportC->getSpaceDimension()==2)
          {
            MEDCouplingNormalizedUnstructuredMesh<2,2> source_mesh_wrapper(source_supportC);
            INTERP_KERNEL::Interpolation2D interpolation(*this);
            interpolation.fromIntegralUniform(source_mesh_wrapper,surfaces,srcMeth);
          }
        else if(source_supportC->getMeshDimension()==3 && source_supportC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(source_supportC);
            INTERP_KERNEL::Interpolation3D interpolation(*this);
            interpolation.fromIntegralUniform(source_mesh_wrapper,surfaces,srcMeth);
          }
        else if(source_supportC->getMeshDimension()==2 && source_supportC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,2> source_mesh_wrapper(source_supportC);
            INTERP_KERNEL::Interpolation3DSurf interpolation(*this);
            interpolation.fromIntegralUniform(source_mesh_wrapper,surfaces,srcMeth);
          }
        else
          throw INTERP_KERNEL::Exception("No para interpolation available for the given mesh and space dimension of source mesh to -1D targetMesh");
      }
    else if ( source_supportC->getMeshDimension() == -1 )
      {
        if(distant_supportC->getMeshDimension()==2 && distant_supportC->getSpaceDimension()==2)
          {
            MEDCouplingNormalizedUnstructuredMesh<2,2> distant_mesh_wrapper(distant_supportC);
            INTERP_KERNEL::Interpolation2D interpolation(*this);
            interpolation.toIntegralUniform(distant_mesh_wrapper,surfaces,srcMeth);
          }
        else if(distant_supportC->getMeshDimension()==3 && distant_supportC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,3> distant_mesh_wrapper(distant_supportC);
            INTERP_KERNEL::Interpolation3D interpolation(*this);
            interpolation.toIntegralUniform(distant_mesh_wrapper,surfaces,srcMeth);
          }
        else if(distant_supportC->getMeshDimension()==2 && distant_supportC->getSpaceDimension()==3)
          {
            MEDCouplingNormalizedUnstructuredMesh<3,2> distant_mesh_wrapper(distant_supportC);
            INTERP_KERNEL::Interpolation3DSurf interpolation(*this);
            interpolation.toIntegralUniform(distant_mesh_wrapper,surfaces,srcMeth);
          }
        else
          throw INTERP_KERNEL::Exception("No para interpolation available for the given mesh and space dimension of distant mesh to -1D sourceMesh");
      }
    else if ( distant_support.getMeshDimension() == 2
              && _source_support->getMeshDimension() == 3
              && distant_support.getSpaceDimension() == 3 && _source_support->getSpaceDimension() == 3)
      {
        MEDCouplingNormalizedUnstructuredMesh<3,3> target_wrapper(distant_supportC);
        MEDCouplingNormalizedUnstructuredMesh<3,3> source_wrapper(source_supportC);
        INTERP_KERNEL::Interpolation2D3D interpolator (*this);
        interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod);
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if ( distant_support.getMeshDimension() == 1
              && _source_support->getMeshDimension() == 2
              && distant_support.getSpaceDimension() == 2 && _source_support->getSpaceDimension() == 2)
      {
        MEDCouplingNormalizedUnstructuredMesh<2,2> target_wrapper(distant_supportC);
        MEDCouplingNormalizedUnstructuredMesh<2,2> source_wrapper(source_supportC);
        INTERP_KERNEL::Interpolation2D1D interpolator (*this);
        interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod);
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if ( distant_support.getMeshDimension() == 3
              && _source_support->getMeshDimension() == 1
              && distant_support.getSpaceDimension() == 3 && _source_support->getSpaceDimension() == 3)
      {
        MEDCouplingNormalizedUnstructuredMesh<3,3> target_wrapper(distant_supportC);
        MEDCouplingNormalizedUnstructuredMesh<3,3> source_wrapper(source_supportC);
        INTERP_KERNEL::Interpolation3D interpolator (*this);
        interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod);
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if (distant_support.getMeshDimension() != _source_support->getMeshDimension())
      {
        throw INTERP_KERNEL::Exception("local and distant meshes do not have the same space and mesh dimensions");
      }
    else if( distant_support.getMeshDimension() == 1
             && distant_support.getSpaceDimension() == 1 )
      {
        MEDCouplingNormalizedUnstructuredMesh<1,1> target_wrapper(distant_supportC);
        MEDCouplingNormalizedUnstructuredMesh<1,1> source_wrapper(source_supportC);

        INTERP_KERNEL::Interpolation1D interpolation(*this);
        interpolation.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod);
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if( distant_support.getMeshDimension() == 1
             && distant_support.getSpaceDimension() == 2 )
      {
        MEDCouplingNormalizedUnstructuredMesh<2,1> target_wrapper(distant_supportC);
        MEDCouplingNormalizedUnstructuredMesh<2,1> source_wrapper(source_supportC);

        INTERP_KERNEL::Interpolation2DCurve interpolation(*this);
        interpolation.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod);
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if ( distant_support.getMeshDimension() == 2
              && distant_support.getSpaceDimension() == 3 )
      {
        MEDCouplingNormalizedUnstructuredMesh<3,2> target_wrapper(distant_supportC);
        MEDCouplingNormalizedUnstructuredMesh<3,2> source_wrapper(source_supportC);

        INTERP_KERNEL::Interpolation3DSurf interpolator (*this);
        interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod);
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if ( distant_support.getMeshDimension() == 2
              && distant_support.getSpaceDimension() == 2)
      {
        MEDCouplingNormalizedUnstructuredMesh<2,2> target_wrapper(distant_supportC);
        MEDCouplingNormalizedUnstructuredMesh<2,2> source_wrapper(source_supportC);

        INTERP_KERNEL::Interpolation2D interpolator (*this);
        interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod);
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else if ( distant_support.getMeshDimension() == 3
              && distant_support.getSpaceDimension() == 3 )
      {
        MEDCouplingNormalizedUnstructuredMesh<3,3> target_wrapper(distant_supportC);
        MEDCouplingNormalizedUnstructuredMesh<3,3> source_wrapper(source_supportC);

        INTERP_KERNEL::Interpolation3D interpolator (*this);
        interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod);
        target_wrapper.releaseTempArrays();
        source_wrapper.releaseTempArrays();
      }
    else
      {
        throw INTERP_KERNEL::Exception("no interpolator exists for these mesh and space dimensions ");
      }
    bool needTargetSurf=isSurfaceComputationNeeded(targetMeth);

    MEDCouplingFieldDouble *target_triangle_surf=0;
    if(needTargetSurf)
      target_triangle_surf = distant_support.getMeasureField(getMeasureAbsStatus());
    fillDSFromVM(iproc_distant,distant_elems,surfaces,target_triangle_surf);

    if(needTargetSurf)
      target_triangle_surf->decrRef();
  }

  void InterpolationMatrix::fillDSFromVM(int iproc_distant, const int* distant_elems, const std::vector< std::map<int,double> >& values, MEDCouplingFieldDouble *surf)
  {
    //loop over the elements to build the interpolation
    //matrix structures
    int source_size=values.size();
    for (int ielem=0; ielem < source_size; ielem++) 
      {
        _row_offsets[ielem+1] += values[ielem].size();
        for(map<int,double>::const_iterator iter=values[ielem].begin();iter!=values[ielem].end();iter++)
          {
            int localId;
            if(distant_elems)
              localId=distant_elems[iter->first];
            else
              localId=iter->first;
            //locating the (iproc, itriangle) pair in the list of columns
            map<pair<int,int>,int >::iterator iter2 = _col_offsets.find(make_pair(iproc_distant,localId));
            int col_id;

            if (iter2 == _col_offsets.end())
              {
                //(iproc, itriangle) is not registered in the list
                //of distant elements
                col_id =_col_offsets.size();
                _col_offsets.insert(make_pair(make_pair(iproc_distant,localId),col_id));
                _mapping.addElementFromSource(iproc_distant,localId);
              }
            else 
              {
                col_id = iter2->second;
              }
            //the non zero coefficient is stored 
            //ielem is the row,
            //col_id is the number of the column
            //iter->second is the value of the coefficient
            if(surf)
              {
                double surface = surf->getIJ(iter->first,0);
                _target_volume[ielem].push_back(surface);
              }
            _coeffs[ielem].push_back(make_pair(col_id,iter->second));
          }
      }
  }

  void InterpolationMatrix::serializeMe(std::vector< std::vector< std::map<int,double> > >& data1, std::vector<int>& data2) const
  {
    data1.clear();
    data2.clear();
    const std::vector<std::pair<int,int> >& sendingIds=_mapping.getSendingIds();
    std::set<int> procsS;
    for(std::vector<std::pair<int,int> >::const_iterator iter1=sendingIds.begin();iter1!=sendingIds.end();iter1++)
      procsS.insert((*iter1).first);
    data1.resize(procsS.size());
    data2.resize(procsS.size());
    std::copy(procsS.begin(),procsS.end(),data2.begin());
    std::map<int,int> fastProcAcc;
    int id=0;
    for(std::set<int>::const_iterator iter2=procsS.begin();iter2!=procsS.end();iter2++,id++)
      fastProcAcc[*iter2]=id;
    int nbOfSrcElt=_coeffs.size();
    for(std::vector< std::vector< std::map<int,double> > >::iterator iter3=data1.begin();iter3!=data1.end();iter3++)
      (*iter3).resize(nbOfSrcElt);
    id=0;
    for(std::vector< std::vector< std::pair<int,double> > >::const_iterator iter4=_coeffs.begin();iter4!=_coeffs.end();iter4++,id++)
      {
        for(std::vector< std::pair<int,double> >::const_iterator iter5=(*iter4).begin();iter5!=(*iter4).end();iter5++)
          {
            const std::pair<int,int>& elt=sendingIds[(*iter5).first];
            data1[fastProcAcc[elt.first]][id][elt.second]=(*iter5).second;
          }
      }
  }

  void InterpolationMatrix::initialize()
  {
    int lgth=_coeffs.size();
    _row_offsets.clear(); _row_offsets.resize(lgth+1);
    _coeffs.clear(); _coeffs.resize(lgth);
    _target_volume.clear(); _target_volume.resize(lgth);
    _col_offsets.clear();
    _mapping.initialize();
  }

  void InterpolationMatrix::finishContributionW(ElementLocator& elementLocator)
  {
    NatureOfField nature=elementLocator.getLocalNature();
    switch(nature)
      {
      case IntensiveMaximum:
        computeConservVolDenoW(elementLocator);
        break;
      case ExtensiveMaximum:
        {
          if(!elementLocator.isM1DCorr())
            computeIntegralDenoW(elementLocator);
          else
            computeGlobConstraintDenoW(elementLocator);
          break;
        }
      case ExtensiveConservation:
        computeGlobConstraintDenoW(elementLocator);
        break;
      case IntensiveConservation:
        {
          if(!elementLocator.isM1DCorr())
            computeRevIntegralDenoW(elementLocator);
          else
            computeConservVolDenoW(elementLocator);
          break;
        }
      default:
        throw INTERP_KERNEL::Exception("Not recognized nature of field. Change nature of Field.");
        break;
      }
  }

  void InterpolationMatrix::finishContributionL(ElementLocator& elementLocator)
  {
    NatureOfField nature=elementLocator.getLocalNature();
    switch(nature)
      {
      case IntensiveMaximum:
        computeConservVolDenoL(elementLocator);
        break;
      case ExtensiveMaximum:
        {
          if(!elementLocator.isM1DCorr())
            computeIntegralDenoL(elementLocator);
          else
            computeConservVolDenoL(elementLocator);
          break;
        }
      case ExtensiveConservation:
        //this is not a bug doing like IntensiveMaximum
        computeConservVolDenoL(elementLocator);
        break;
      case IntensiveConservation:
        {
          if(!elementLocator.isM1DCorr())
            computeRevIntegralDenoL(elementLocator);
          else
            computeConservVolDenoL(elementLocator);
          break;
        }
      default:
        throw INTERP_KERNEL::Exception("Not recognized nature of field. Change nature of Field.");
        break;
      }
  }
  
  void InterpolationMatrix::computeConservVolDenoW(ElementLocator& elementLocator)
  {
    computeGlobalColSum(_deno_reverse_multiply);
    computeGlobalRowSum(elementLocator,_deno_multiply,_deno_reverse_multiply);
  }
  
  void InterpolationMatrix::computeConservVolDenoL(ElementLocator& elementLocator)
  {
    int pol1=elementLocator.sendPolicyToWorkingSideL();
    if(pol1==ElementLocator::NO_POST_TREATMENT_POLICY)
      {
        elementLocator.recvFromWorkingSideL();
        elementLocator.sendToWorkingSideL();
      }
    else if(ElementLocator::CUMULATIVE_POLICY)
      {
        //ask for lazy side to deduce ids eventually missing on working side and to send it back.
        elementLocator.recvLocalIdsFromWorkingSideL();
        elementLocator.sendCandidatesGlobalIdsToWorkingSideL();
        elementLocator.recvCandidatesForAddElementsL();
        elementLocator.sendAddElementsToWorkingSideL();
        //Working side has updated its eventually missing ids updates its global ids with lazy side procs contribution
        elementLocator.recvLocalIdsFromWorkingSideL();
        elementLocator.sendGlobalIdsToWorkingSideL();
        //like no post treatment
        elementLocator.recvFromWorkingSideL();
        elementLocator.sendToWorkingSideL();
      }
    else
      throw INTERP_KERNEL::Exception("Not managed policy detected on lazy side : not implemented !");
  }

  void InterpolationMatrix::computeIntegralDenoW(ElementLocator& elementLocator)
  {
    MEDCouplingFieldDouble *source_triangle_surf = _source_support->getMeasureField(getMeasureAbsStatus());
    _deno_multiply.resize(_coeffs.size());
    vector<vector<double> >::iterator iter6=_deno_multiply.begin();
    const double *values=source_triangle_surf->getArray()->getConstPointer();
    for(vector<vector<pair<int,double> > >::const_iterator iter4=_coeffs.begin();iter4!=_coeffs.end();iter4++,iter6++,values++)
      {
        (*iter6).resize((*iter4).size());
        std::fill((*iter6).begin(),(*iter6).end(),*values);
      }
    source_triangle_surf->decrRef();
    _deno_reverse_multiply=_target_volume;
  }

  void InterpolationMatrix::computeRevIntegralDenoW(ElementLocator& elementLocator)
  {
    _deno_multiply=_target_volume;
    MEDCouplingFieldDouble *source_triangle_surf = _source_support->getMeasureField(getMeasureAbsStatus());
    _deno_reverse_multiply.resize(_coeffs.size());
    vector<vector<double> >::iterator iter6=_deno_reverse_multiply.begin();
    const double *values=source_triangle_surf->getArray()->getConstPointer();
    for(vector<vector<pair<int,double> > >::const_iterator iter4=_coeffs.begin();iter4!=_coeffs.end();iter4++,iter6++,values++)
      {
        (*iter6).resize((*iter4).size());
        std::fill((*iter6).begin(),(*iter6).end(),*values);
      }
    source_triangle_surf->decrRef();
  }
  
  /*!
   * Nothing to do because surface computation is on working side.
   */
  void InterpolationMatrix::computeIntegralDenoL(ElementLocator& elementLocator)
  {
  }

  /*!
   * Nothing to do because surface computation is on working side.
   */
  void InterpolationMatrix::computeRevIntegralDenoL(ElementLocator& elementLocator)
  {
  }


  void InterpolationMatrix::computeGlobConstraintDenoW(ElementLocator& elementLocator)
  {
    computeGlobalColSum(_deno_multiply);
    computeGlobalRowSum(elementLocator,_deno_reverse_multiply,_deno_multiply);
  }

  void InterpolationMatrix::computeGlobalRowSum(ElementLocator& elementLocator, std::vector<std::vector<double> >& denoStrorage, std::vector<std::vector<double> >& denoStrorageInv)
  {
    //stores id in distant procs sorted by lazy procs connected with
    vector< vector<int> > rowsPartialSumI;
    //stores for each lazy procs connected with, if global info is available and if it's the case the policy
    vector<int> policyPartial;
    //stores the corresponding values.
    vector< vector<double> > rowsPartialSumD;
    elementLocator.recvPolicyFromLazySideW(policyPartial);
    int pol1=mergePolicies(policyPartial);
    if(pol1==ElementLocator::NO_POST_TREATMENT_POLICY)
      {
        computeLocalRowSum(elementLocator.getDistantProcIds(),rowsPartialSumI,rowsPartialSumD);
        elementLocator.sendSumToLazySideW(rowsPartialSumI,rowsPartialSumD);
        elementLocator.recvSumFromLazySideW(rowsPartialSumD);
      }
    else if(pol1==ElementLocator::CUMULATIVE_POLICY)
      {
        //updateWithNewAdditionnalElements(addingElements);
        //stores for each lazy procs connected with, the ids in global mode if it exists (regarding policyPartial). This array has exactly the size of  rowsPartialSumI,
        //if policyPartial has CUMALATIVE_POLICY in each.
        vector< vector<int> > globalIdsPartial;
        computeLocalRowSum(elementLocator.getDistantProcIds(),rowsPartialSumI,rowsPartialSumD);
        elementLocator.sendLocalIdsToLazyProcsW(rowsPartialSumI);
        elementLocator.recvCandidatesGlobalIdsFromLazyProcsW(globalIdsPartial);
        std::vector< std::vector<int> > addingElements;
        findAdditionnalElements(elementLocator,addingElements,rowsPartialSumI,globalIdsPartial);
        addGhostElements(elementLocator.getDistantProcIds(),addingElements);
        rowsPartialSumI.clear();
        globalIdsPartial.clear();
        computeLocalRowSum(elementLocator.getDistantProcIds(),rowsPartialSumI,rowsPartialSumD);
        elementLocator.sendLocalIdsToLazyProcsW(rowsPartialSumI);
        elementLocator.recvGlobalIdsFromLazyProcsW(rowsPartialSumI,globalIdsPartial);
        //
        elementLocator.sendSumToLazySideW(rowsPartialSumI,rowsPartialSumD);
        elementLocator.recvSumFromLazySideW(rowsPartialSumD);
        mergeRowSum3(globalIdsPartial,rowsPartialSumD);
        mergeCoeffs(elementLocator.getDistantProcIds(),rowsPartialSumI,globalIdsPartial,denoStrorageInv);
      }
    else
      throw INTERP_KERNEL::Exception("Not managed policy detected : not implemented !");
    divideByGlobalRowSum(elementLocator.getDistantProcIds(),rowsPartialSumI,rowsPartialSumD,denoStrorage);
  }

  /*!
   * @param distantProcs in parameter that indicates which lazy procs are concerned.
   * @param resPerProcI out parameter that must be cleared before calling this method. The size of 1st dimension is equal to the size of 'distantProcs'.
   *                    It contains the element ids (2nd dimension) of the corresponding lazy proc.
   * @param  resPerProcD out parameter with the same format than 'resPerProcI'. It contains corresponding sum values.
   */
  void InterpolationMatrix::computeLocalRowSum(const std::vector<int>& distantProcs, std::vector<std::vector<int> >& resPerProcI,
                                               std::vector<std::vector<double> >& resPerProcD) const
  {
    resPerProcI.resize(distantProcs.size());
    resPerProcD.resize(distantProcs.size());
    std::vector<double> res(_col_offsets.size());
    for(vector<vector<pair<int,double> > >::const_iterator iter=_coeffs.begin();iter!=_coeffs.end();iter++)
      for(vector<pair<int,double> >::const_iterator iter3=(*iter).begin();iter3!=(*iter).end();iter3++)
        res[(*iter3).first]+=(*iter3).second;
    set<int> procsSet;
    int id=-1;
    const vector<std::pair<int,int> >& mapping=_mapping.getSendingIds();
    for(vector<std::pair<int,int> >::const_iterator iter2=mapping.begin();iter2!=mapping.end();iter2++)
      {
        std::pair<set<int>::iterator,bool> isIns=procsSet.insert((*iter2).first);
        if(isIns.second)
          id=std::find(distantProcs.begin(),distantProcs.end(),(*iter2).first)-distantProcs.begin();
        resPerProcI[id].push_back((*iter2).second);
        resPerProcD[id].push_back(res[iter2-mapping.begin()]);
      }
  }

  /*!
   * This method is only usable when CUMULATIVE_POLICY detected. This method finds elements ids (typically nodes) lazy side that
   * are not present in columns of 'this' and that should regarding cumulative merge of elements regarding their global ids.
   */
  void InterpolationMatrix::findAdditionnalElements(ElementLocator& elementLocator, std::vector<std::vector<int> >& elementsToAdd,
                                                    const std::vector<std::vector<int> >& resPerProcI, const std::vector<std::vector<int> >& globalIdsPartial)
  {
    std::set<int> globalIds;
    int nbLazyProcs=globalIdsPartial.size();
    for(int i=0;i<nbLazyProcs;i++)
      globalIds.insert(globalIdsPartial[i].begin(),globalIdsPartial[i].end());
    std::vector<int> tmp(globalIds.size());
    std::copy(globalIds.begin(),globalIds.end(),tmp.begin());
    globalIds.clear();
    elementLocator.sendCandidatesForAddElementsW(tmp);
    elementLocator.recvAddElementsFromLazyProcsW(elementsToAdd);
  }

  void InterpolationMatrix::addGhostElements(const std::vector<int>& distantProcs, const std::vector<std::vector<int> >& elementsToAdd)
  {
    std::vector< std::vector< std::map<int,double> > > data1;
    std::vector<int> data2;
    serializeMe(data1,data2);
    initialize();
    int nbOfDistProcs=distantProcs.size();
    for(int i=0;i<nbOfDistProcs;i++)
      {
        int procId=distantProcs[i];
        const std::vector<int>& eltsForThisProc=elementsToAdd[i];
        if(!eltsForThisProc.empty())
          {
            std::vector<int>::iterator iter1=std::find(data2.begin(),data2.end(),procId);
            std::map<int,double> *toFeed=0;
            if(iter1!=data2.end())
              {//to test
                int rank=iter1-data2.begin();
                toFeed=&(data1[rank].back());
              }
            else
              {
                iter1=std::lower_bound(data2.begin(),data2.end(),procId);
                int rank=iter1-data2.begin();
                data2.insert(iter1,procId);
                std::vector< std::map<int,double> > tmp(data1.front().size());
                data1.insert(data1.begin()+rank,tmp);
                toFeed=&(data1[rank].back());
              }
            for(std::vector<int>::const_iterator iter2=eltsForThisProc.begin();iter2!=eltsForThisProc.end();iter2++)
              (*toFeed)[*iter2]=0.;
          }
      }
    //
    nbOfDistProcs=data2.size();
    for(int j=0;j<nbOfDistProcs;j++)
      fillDSFromVM(data2[j],0,data1[j],0);
  }

  int InterpolationMatrix::mergePolicies(const std::vector<int>& policyPartial)
  {
    if(policyPartial.empty())
      return ElementLocator::NO_POST_TREATMENT_POLICY;
    int ref=policyPartial[0];
     std::vector<int>::const_iterator iter1=std::find_if(policyPartial.begin(),policyPartial.end(),std::bind2nd(std::not_equal_to<int>(),ref));
    if(iter1!=policyPartial.end())
      {
        std::ostringstream msg; msg << "Incompatible policies between lazy procs each other : proc # " << iter1-policyPartial.begin();
        throw INTERP_KERNEL::Exception(msg.str().c_str());
      }
    return ref;
  }

  /*!
   * This method introduce global ids aspects in computed 'rowsPartialSumD'.
   * As precondition rowsPartialSumD.size()==policyPartial.size()==globalIdsPartial.size(). Foreach i in [0;rowsPartialSumD.size() ) rowsPartialSumD[i].size()==globalIdsPartial[i].size()
   * @param rowsPartialSumD : in parameter, Partial row sum computed for each lazy procs connected with.
   * @param rowsPartialSumI : in parameter, Corresponding local ids for each lazy procs connected with.
   * @param globalIdsPartial : in parameter, the global numbering, of elements connected with.
   * @param globalIdsLazySideInteraction : out parameter, constituted from all global ids of lazy procs connected with.
   * @para sumCorresponding : out parameter, relative to 'globalIdsLazySideInteraction'
   */
  void InterpolationMatrix::mergeRowSum(const std::vector< std::vector<double> >& rowsPartialSumD, const std::vector< std::vector<int> >& globalIdsPartial,
                                        std::vector<int>& globalIdsLazySideInteraction, std::vector<double>& sumCorresponding)
  {
    std::map<int,double> sumToReturn;
    int nbLazyProcs=rowsPartialSumD.size();
    for(int i=0;i<nbLazyProcs;i++)
      {
        const std::vector<double>& rowSumOfP=rowsPartialSumD[i];
        const std::vector<int>& globalIdsOfP=globalIdsPartial[i];
        std::vector<double>::const_iterator iter1=rowSumOfP.begin();
        std::vector<int>::const_iterator iter2=globalIdsOfP.begin();
        for(;iter1!=rowSumOfP.end();iter1++,iter2++)
          sumToReturn[*iter2]+=*iter1;
      }
    //
    int lgth=sumToReturn.size();
    globalIdsLazySideInteraction.resize(lgth);
    sumCorresponding.resize(lgth);
    std::vector<int>::iterator iter3=globalIdsLazySideInteraction.begin();
    std::vector<double>::iterator iter4=sumCorresponding.begin();
    for(std::map<int,double>::const_iterator iter5=sumToReturn.begin();iter5!=sumToReturn.end();iter5++,iter3++,iter4++)
      {
        *iter3=(*iter5).first;
        *iter4=(*iter5).second;
      }
  }

  /*!
   * This method simply reorganize the result contained in 'sumCorresponding' computed by lazy side into 'rowsPartialSumD' with help of 'globalIdsPartial' and 'globalIdsLazySideInteraction'
   *
   * @param globalIdsPartial : in parameter, global ids sorted by lazy procs
   * @param rowsPartialSumD : in/out parameter, with exactly the same size as 'globalIdsPartial'
   * @param globalIdsLazySideInteraction : in parameter that represents ALL the global ids of every lazy procs in interaction
   * @param sumCorresponding : in parameter with same size as 'globalIdsLazySideInteraction' that stores the corresponding sum of 'globalIdsLazySideInteraction'
   */
  void InterpolationMatrix::mergeRowSum2(const std::vector< std::vector<int> >& globalIdsPartial, std::vector< std::vector<double> >& rowsPartialSumD,
                                         const std::vector<int>& globalIdsLazySideInteraction, const std::vector<double>& sumCorresponding)
  {
    std::map<int,double> acc;
    std::vector<int>::const_iterator iter1=globalIdsLazySideInteraction.begin();
    std::vector<double>::const_iterator iter2=sumCorresponding.begin();
    for(;iter1!=globalIdsLazySideInteraction.end();iter1++,iter2++)
      acc[*iter1]=*iter2;
    //
    int nbLazyProcs=globalIdsPartial.size();
    for(int i=0;i<nbLazyProcs;i++)
      {
        const std::vector<int>& tmp1=globalIdsPartial[i];
        std::vector<double>& tmp2=rowsPartialSumD[i];
        std::vector<int>::const_iterator iter3=tmp1.begin();
        std::vector<double>::iterator iter4=tmp2.begin();
        for(;iter3!=tmp1.end();iter3++,iter4++)
          *iter4=acc[*iter3];
      }
  }
  
  void InterpolationMatrix::mergeRowSum3(const std::vector< std::vector<int> >& globalIdsPartial, std::vector< std::vector<double> >& rowsPartialSumD)
  {
    std::map<int,double> sum;
    std::vector< std::vector<int> >::const_iterator iter1=globalIdsPartial.begin();
    std::vector< std::vector<double> >::iterator iter2=rowsPartialSumD.begin();
    for(;iter1!=globalIdsPartial.end();iter1++,iter2++)
      {
        std::vector<int>::const_iterator iter3=(*iter1).begin();
        std::vector<double>::const_iterator iter4=(*iter2).begin();
        for(;iter3!=(*iter1).end();iter3++,iter4++)
          sum[*iter3]+=*iter4;
      }
    iter2=rowsPartialSumD.begin();
    for(iter1=globalIdsPartial.begin();iter1!=globalIdsPartial.end();iter1++,iter2++)
      {
        std::vector<int>::const_iterator iter3=(*iter1).begin();
        std::vector<double>::iterator iter4=(*iter2).begin();
        for(;iter3!=(*iter1).end();iter3++,iter4++)
          *iter4=sum[*iter3];
      }
  }

  /*!
   * This method updates this->_coeffs attribute in order to take into account hidden (because having the same global number) similar nodes in _coeffs array.
   * If in this->_coeffs two distant element id have the same global id their values will be replaced for each by the sum of the two.
   * @param procsInInteraction input parameter : specifies the procId in absolute of distant lazy procs in interaction with
   * @param rowsPartialSumI input parameter : local ids of distant lazy procs elements in interaction with
   * @param globalIdsPartial input parameter : global ids of distant lazy procs elements in interaction with
   */
  void InterpolationMatrix::mergeCoeffs(const std::vector<int>& procsInInteraction, const std::vector< std::vector<int> >& rowsPartialSumI,
                                        const std::vector<std::vector<int> >& globalIdsPartial, std::vector<std::vector<double> >& denoStrorageInv)
  {
    //preparing fast access structures
    std::map<int,int> procT;
    int localProcId=0;
    for(std::vector<int>::const_iterator iter1=procsInInteraction.begin();iter1!=procsInInteraction.end();iter1++,localProcId++)
      procT[*iter1]=localProcId;
    int size=procsInInteraction.size();
    std::vector<std::map<int,int> > localToGlobal(size);
    for(int i=0;i<size;i++)
      {
        std::map<int,int>& myLocalToGlobal=localToGlobal[i];
        const std::vector<int>& locals=rowsPartialSumI[i];
        const std::vector<int>& globals=globalIdsPartial[i];
        std::vector<int>::const_iterator iter3=locals.begin();
        std::vector<int>::const_iterator iter4=globals.begin();
        for(;iter3!=locals.end();iter3++,iter4++)
          myLocalToGlobal[*iter3]=*iter4;
      }
    //
    const vector<std::pair<int,int> >& mapping=_mapping.getSendingIds();
    std::map<int,double> globalIdVal;
    //accumulate for same global id on lazy part.
    for(vector<vector<pair<int,double> > >::iterator iter1=_coeffs.begin();iter1!=_coeffs.end();iter1++)
      for(vector<pair<int,double> >::iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
        {
          const std::pair<int,int>& distantLocalLazyId=mapping[(*iter2).first];
          int localLazyProcId=procT[distantLocalLazyId.first];
          int globalDistantLazyId=localToGlobal[localLazyProcId][distantLocalLazyId.second];
          globalIdVal[globalDistantLazyId]+=(*iter2).second;
        }
    //perform merge
    std::vector<std::vector<double> >::iterator iter3=denoStrorageInv.begin();
    for(vector<vector<pair<int,double> > >::iterator iter1=_coeffs.begin();iter1!=_coeffs.end();iter1++,iter3++)
      {
        double val=(*iter3).back();
        (*iter3).resize((*iter1).size());
        std::vector<double>::iterator iter4=(*iter3).begin();
        for(vector<pair<int,double> >::iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++,iter4++)
          {
            const std::pair<int,int>& distantLocalLazyId=mapping[(*iter2).first];
            int localLazyProcId=procT[distantLocalLazyId.first];
            int globalDistantLazyId=localToGlobal[localLazyProcId][distantLocalLazyId.second];
            double newVal=globalIdVal[globalDistantLazyId];
            if((*iter2).second!=0.)
              (*iter4)=val*newVal/(*iter2).second;
            else
              (*iter4)=std::numeric_limits<double>::max();
            (*iter2).second=newVal;
          }
      }
  }

  void InterpolationMatrix::divideByGlobalRowSum(const std::vector<int>& distantProcs, const std::vector<std::vector<int> >& resPerProcI,
                                                 const std::vector<std::vector<double> >& resPerProcD, std::vector<std::vector<double> >& deno)
  {
    map<int,double> fastSums;
    int procId=0;
    for(vector<int>::const_iterator iter1=distantProcs.begin();iter1!=distantProcs.end();iter1++,procId++)
      {
        const std::vector<int>& currentProcI=resPerProcI[procId];
        const std::vector<double>& currentProcD=resPerProcD[procId];
        vector<double>::const_iterator iter3=currentProcD.begin();
        for(vector<int>::const_iterator iter2=currentProcI.begin();iter2!=currentProcI.end();iter2++,iter3++)
          fastSums[_col_offsets[std::make_pair(*iter1,*iter2)]]=*iter3;
      }
    deno.resize(_coeffs.size());
    vector<vector<double> >::iterator iter6=deno.begin();
    for(vector<vector<pair<int,double> > >::const_iterator iter4=_coeffs.begin();iter4!=_coeffs.end();iter4++,iter6++)
      {
        (*iter6).resize((*iter4).size());
        vector<double>::iterator iter7=(*iter6).begin();
        for(vector<pair<int,double> >::const_iterator iter5=(*iter4).begin();iter5!=(*iter4).end();iter5++,iter7++)
          *iter7=fastSums[(*iter5).first];
      }
  }

  void InterpolationMatrix::computeGlobalColSum(std::vector<std::vector<double> >& denoStrorage)
  {
    denoStrorage.resize(_coeffs.size());
    vector<vector<double> >::iterator iter2=denoStrorage.begin();
    for(vector<vector<pair<int,double> > >::const_iterator iter1=_coeffs.begin();iter1!=_coeffs.end();iter1++,iter2++)
      {
        (*iter2).resize((*iter1).size());
        double sumOfCurrentRow=0.;
        for(vector<pair<int,double> >::const_iterator iter3=(*iter1).begin();iter3!=(*iter1).end();iter3++)
          sumOfCurrentRow+=(*iter3).second;
        std::fill((*iter2).begin(),(*iter2).end(),sumOfCurrentRow);
      }
  }

  void InterpolationMatrix::resizeGlobalColSum(std::vector<std::vector<double> >& denoStrorage)
  {
    vector<vector<double> >::iterator iter2=denoStrorage.begin();
    for(vector<vector<pair<int,double> > >::const_iterator iter1=_coeffs.begin();iter1!=_coeffs.end();iter1++,iter2++)
      {
        double val=(*iter2).back();
        (*iter2).resize((*iter1).size());
        std::fill((*iter2).begin(),(*iter2).end(),val);
      }
  }


 /**!  The call to this method updates the arrays on the target side
     so that they know which amount of data from which processor they
     should expect.
     That call makes actual interpolations via multiply method
     available.
     */
  void InterpolationMatrix::prepare()
  {
    int nbelems = _source_field->getField()->getNumberOfTuples();
    for (int ielem=0; ielem < nbelems; ielem++)
      {
        _row_offsets[ielem+1]+=_row_offsets[ielem];
      }
    _mapping.prepareSendRecv();
  }



  /*!
     \brief performs t=Ws, where t is the target field, s is the source field

     The call to this method must be called both on the working side
     and on the idle side. On the working side, the vector  T=VT^(-1).(W.S)
     is computed and sent. On the idle side, no computation is done, but the
     result from the working side is received and the field is updated.

     \param field source field on processors involved on the source side,
     target field on processors on the target side
   */
  void InterpolationMatrix::multiply(MEDCouplingFieldDouble& field) const
  {
    int nbcomp = field.getArray()->getNumberOfComponents();
    vector<double> target_value(_col_offsets.size()* nbcomp,0.0);

    //computing the matrix multiply on source side
    if (_source_group.containsMyRank())
      {
        int nbrows = _coeffs.size();

        // performing W.S
        // W is the intersection matrix
        // S is the source vector

        for (int irow=0; irow<nbrows; irow++)
          {
            for (int icomp=0; icomp< nbcomp; icomp++)
              {
                double coeff_row = field.getIJ(irow,icomp);
                for (int icol=_row_offsets[irow]; icol< _row_offsets[irow+1];icol++)
                  {
                    int colid= _coeffs[irow][icol-_row_offsets[irow]].first;
                    double value = _coeffs[irow][icol-_row_offsets[irow]].second;
                    double deno = _deno_multiply[irow][icol-_row_offsets[irow]];
                    target_value[colid*nbcomp+icomp]+=value*coeff_row/deno;
                  }
              }
          }
      }

    if (_target_group.containsMyRank())
      {
        int nbelems = field.getArray()->getNumberOfTuples() ;
        double* value = const_cast<double*> (field.getArray()->getPointer());
        for (int i=0; i<nbelems*nbcomp; i++)
          {
            value[i]=0.0;
          }
      }

    //on source side : sending  T=VT^(-1).(W.S)
    //on target side :: receiving T and storing it in field
    _mapping.sendRecv(&target_value[0],field);
  }
  

  /**!
   \brief performs s=WTt, where t is the target field, s is the source field,
   WT is the transpose matrix from W

     The call to this method must be called both on the working side
     and on the idle side. On the working side, the target vector T is
     received and the vector  S=VS^(-1).(WT.T) is computed to update
     the field.
     On the idle side, no computation is done, but the field is sent.

     param field source field on processors involved on the source side,
     target field on processors on the target side
     */
  void InterpolationMatrix::transposeMultiply(MEDCouplingFieldDouble& field) const
  {
    int nbcomp = field.getArray()->getNumberOfComponents();
    vector<double> source_value(_col_offsets.size()* nbcomp,0.0);
    _mapping.reverseSendRecv(&source_value[0],field);

    //treatment of the transpose matrix multiply on the source side
    if (_source_group.containsMyRank())
      {
        int nbrows    = _coeffs.size();
        double *array = field.getArray()->getPointer() ;

        // Initialization
        std::fill(array, array+nbrows*nbcomp, 0.0) ;

        //performing WT.T
        //WT is W transpose
        //T is the target vector
        for (int irow = 0; irow < nbrows; irow++)
          {
            for (int icol = _row_offsets[irow]; icol < _row_offsets[irow+1]; icol++)
              {
                int colid    = _coeffs[irow][icol-_row_offsets[irow]].first;
                double value = _coeffs[irow][icol-_row_offsets[irow]].second;
                double deno = _deno_reverse_multiply[irow][icol-_row_offsets[irow]];
                for (int icomp=0; icomp<nbcomp; icomp++)
                  {
                    double coeff_row = source_value[colid*nbcomp+icomp];
                    array[irow*nbcomp+icomp] += value*coeff_row/deno;
                  }
              }
          }
      }
  }

  bool InterpolationMatrix::isSurfaceComputationNeeded(const std::string& method) const
  {
    return method=="P0";
  }
}
