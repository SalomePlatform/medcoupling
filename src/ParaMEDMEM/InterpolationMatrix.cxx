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
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "ProcessorGroup.hxx"
#include "MxN_Mapping.hxx"
#include "InterpolationMatrix.hxx"
#include "TranslationRotationMatrix.hxx"
#include "Interpolation.hxx"
#include "Interpolation2D.txx"
#include "Interpolation3DSurf.txx"
#include "Interpolation3D.txx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingNormalizedUnstructuredMesh.txx"
#include "InterpolationOptions.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "GlobalizerMesh.hxx"

// class InterpolationMatrix
// This class enables the storage of an interpolation matrix Wij mapping 
// source field Sj to target field Ti via Ti=Vi^(-1).Wij.Sj.
// The matrix is built and stored on the processors belonging to the source
// group. 

using namespace std;

namespace ParaMEDMEM
{

  //   ====================================================================
  //   Creates an empty matrix structure linking two distributed supports.
  //   The method must be called by all processors belonging to source
  //   and target groups.
  //   param source_support local support
  //   param source_group processor group containing the local processors
  //   param target_group processor group containing the distant processors
  //   param method interpolation method
  //   ====================================================================

  InterpolationMatrix::InterpolationMatrix(const ParaMEDMEM::ParaFIELD *source_field, 
                                           const ProcessorGroup& source_group,
                                           const ProcessorGroup& target_group,
                                           const DECOptions& dec_options,
                                           const INTERP_KERNEL::InterpolationOptions& interp_options):
    _source_field(source_field),
    _source_support(source_field->getSupport()->getCellMesh()),
    _mapping(source_group, target_group, dec_options),
    _source_group(source_group),
    _target_group(target_group),
    DECOptions(dec_options),
    INTERP_KERNEL::InterpolationOptions(interp_options)
  {
    int nbelems = source_field->getField()->getNumberOfTuples();
    _row_offsets.resize(nbelems+1);
    _coeffs.resize(nbelems);
    _target_volume.resize(nbelems);
  }

  InterpolationMatrix::~InterpolationMatrix()
  {
  }


  //   ======================================================================
  //   \brief Adds the contribution of a distant subdomain to the*
  //   interpolation matrix.
  //   The method adds contribution to the interpolation matrix.
  //   For each row of the matrix, elements are addded as
  //   a (column, coeff) pair in the _coeffs array. This column number refers
  //   to an element on the target side via the _col_offsets array.
  //   It is made of a series of (iproc, ielem) pairs. 
  //   The number of elements per row is stored in the row_offsets array.

  //   param distant_support local representation of the distant subdomain
  //   param iproc_distant id of the distant subdomain (in the distant group)
  //   param distant_elems mapping between the local representation of
  //   the subdomain and the actual elem ids on the distant subdomain
  //   ======================================================================

  void InterpolationMatrix::addContribution ( MEDCouplingPointSet& distant_support,
                                              int iproc_distant,
                                              int* distant_elems,
                                              const std::string& srcMeth,
                                              const std::string& targetMeth)
  {
    if (distant_support.getMeshDimension() != _source_support->getMeshDimension())
      {
        throw INTERP_KERNEL::Exception("local and distant meshes do not have the same space and mesh dimensions");
      }
    std::string interpMethod(targetMeth);
    interpMethod+=srcMeth;
    //creating the interpolator structure
    vector<map<int,double> > surfaces;
    int colSize=0;
    //computation of the intersection volumes between source and target elements
    MEDCouplingUMesh *distant_supportC=dynamic_cast<MEDCouplingUMesh *>(&distant_support);
    MEDCouplingUMesh *source_supportC=dynamic_cast<MEDCouplingUMesh *>(_source_support);
    if ( distant_support.getMeshDimension() == 2
         && distant_support.getSpaceDimension() == 3 )
      {
        MEDCouplingNormalizedUnstructuredMesh<3,2> target_wrapper(distant_supportC);
        MEDCouplingNormalizedUnstructuredMesh<3,2> source_wrapper(source_supportC);

        INTERP_KERNEL::Interpolation3DSurf interpolator (*this);
        colSize=interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod.c_str());
        target_wrapper.ReleaseTempArrays();
        source_wrapper.ReleaseTempArrays();
      }
    else if ( distant_support.getMeshDimension() == 2
              && distant_support.getSpaceDimension() == 2)
      {
        MEDCouplingNormalizedUnstructuredMesh<2,2> target_wrapper(distant_supportC);
        MEDCouplingNormalizedUnstructuredMesh<2,2> source_wrapper(source_supportC);

        INTERP_KERNEL::Interpolation2D interpolator (*this);
        colSize=interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod.c_str());
        target_wrapper.ReleaseTempArrays();
        source_wrapper.ReleaseTempArrays();
      }
    else if ( distant_support.getMeshDimension() == 3
              && distant_support.getSpaceDimension() == 3 )
      {
        MEDCouplingNormalizedUnstructuredMesh<3,3> target_wrapper(distant_supportC);
        MEDCouplingNormalizedUnstructuredMesh<3,3> source_wrapper(source_supportC);

        INTERP_KERNEL::Interpolation3D interpolator (*this);
        colSize=interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces,interpMethod.c_str());
        target_wrapper.ReleaseTempArrays();
        source_wrapper.ReleaseTempArrays();
      }
    else
      {
        throw INTERP_KERNEL::Exception("no interpolator exists for these mesh and space dimensions ");
      }
  
    int source_size=surfaces.size();
    bool needTargetSurf=isSurfaceComputationNeeded(targetMeth);

    MEDCouplingFieldDouble *target_triangle_surf;
    if(needTargetSurf)
      target_triangle_surf = distant_support.getMeasureField();

    //loop over the elements to build the interpolation
    //matrix structures
    for (int ielem=0; ielem < source_size; ielem++) 
      {
        _row_offsets[ielem+1] += surfaces[ielem].size();
        //_source_indices.push_back(make_pair(iproc_distant, ielem));

        for (map<int,double>::const_iterator iter = surfaces[ielem].begin();
             iter != surfaces[ielem].end();
             iter++)
          {
            //surface of the target triangle
            double surf;
            if(needTargetSurf)
              surf = target_triangle_surf->getIJ(iter->first,0);

            //locating the (iproc, itriangle) pair in the list of columns
            map<pair<int,int>,int >::iterator iter2 = _col_offsets.find(make_pair(iproc_distant,iter->first));
            int col_id;

            if (iter2 == _col_offsets.end())
              {
                //(iproc, itriangle) is not registered in the list
                //of distant elements

                col_id =_col_offsets.size();
                _col_offsets.insert(make_pair(make_pair(iproc_distant,iter->first),col_id));
                _mapping.addElementFromSource(iproc_distant,
                                              distant_elems[iter->first]);
                //target_volume.push_back(surf);
              }
            else 
              {
                col_id = iter2->second;
              }

            //the non zero coefficient is stored 
            //ielem is the row,
            //col_id is the number of the column
            //iter->second is the value of the coefficient
            if(needTargetSurf)
              _target_volume[ielem].push_back(surf);
            _coeffs[ielem].push_back(make_pair(col_id,iter->second));
          }
      }
    if(needTargetSurf)
      target_triangle_surf->decrRef();
  }

  void InterpolationMatrix::finishContributionW(GlobalizerMeshWorkingSide& globalizer)
  {
    NatureOfField nature=globalizer.getLocalNature();
    switch(nature)
      {
      case ConservativeVolumic:
        computeConservVolDenoW(globalizer);
        break;
      case Integral:
        computeIntegralDenoW(globalizer);
        break;
      case IntegralGlobConstraint:
        computeGlobConstraintDenoW(globalizer);
        break;
      default:
        throw INTERP_KERNEL::Exception("Not recognized nature of field. Change nature of Field.");
        break;
      }
    /*for(map<pair<int,int>,int>::iterator iter=_col_offsets.begin();iter!=_col_offsets.end();iter++)
          if((*iter).second==9)
            cout << (*iter).first.first << " ** " << (*iter).first.second << endl;
    cout << "--> " << _col_offsets[make_pair(4,74)] << endl;
    for(vector<vector<pair<int,double> > >::iterator iter3=_coeffs.begin();iter3!=_coeffs.end();iter3++)
      for(vector<pair<int,double> >::iterator iter4=(*iter3).begin();iter4!=(*iter3).end();iter4++)
        if((*iter4).first==569)
          cout << " __ " << iter3-_coeffs.begin() << "___" << (*iter4).second << endl; 
    ostringstream st; st << "deno_" << _deno_multiply[0][0];
    ofstream fs(st.str().c_str());
    for(vector<vector<double> >::iterator iter1=_deno_multiply.begin();iter1!=_deno_multiply.end();iter1++)
      {
        for(vector<double>::iterator iter2=(*iter1).begin();iter2!=(*iter1).end();iter2++)
          fs << *iter2 << " ";
        fs << endl;
        }*/
  }

  void InterpolationMatrix::finishContributionL(GlobalizerMeshLazySide& globalizer)
  {
    NatureOfField nature=globalizer.getLocalNature();
    switch(nature)
      {
      case ConservativeVolumic:
        computeConservVolDenoL(globalizer);
        break;
      case Integral:
        computeIntegralDenoL(globalizer);
        break;
      case IntegralGlobConstraint:
        computeGlobConstraintDenoL(globalizer);
        break;
      default:
        throw INTERP_KERNEL::Exception("Not recognized nature of field. Change nature of Field.");
        break;
      }
  }
  
  void InterpolationMatrix::computeConservVolDenoW(GlobalizerMeshWorkingSide& globalizer)
  {
    computeGlobalRowSum(globalizer,_deno_multiply);
    computeGlobalColSum(_deno_reverse_multiply);
  }
  
  void InterpolationMatrix::computeConservVolDenoL(GlobalizerMeshLazySide& globalizer)
  {
    globalizer.recvFromWorkingSide();
    globalizer.sendToWorkingSide();
  }

  void InterpolationMatrix::computeIntegralDenoW(GlobalizerMeshWorkingSide& globalizer)
  {
    MEDCouplingFieldDouble *source_triangle_surf = _source_support->getMeasureField();
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
  
  /*!
   * Nothing to do because surface computation is on working side.
   */
  void InterpolationMatrix::computeIntegralDenoL(GlobalizerMeshLazySide& globalizer)
  {
  }

  void InterpolationMatrix::computeGlobConstraintDenoW(GlobalizerMeshWorkingSide& globalizer)
  {
    computeGlobalColSum(_deno_multiply);
    computeGlobalRowSum(globalizer,_deno_reverse_multiply);
  }

  void InterpolationMatrix::computeGlobConstraintDenoL(GlobalizerMeshLazySide& globalizer)
  {
    globalizer.recvFromWorkingSide();
    globalizer.sendToWorkingSide();
  }

  void InterpolationMatrix::computeGlobalRowSum(GlobalizerMeshWorkingSide& globalizer, std::vector<std::vector<double> >& denoStrorage)
  {
    //stores id in distant procs sorted by lazy procs connected with
    vector< vector<int> > rowsPartialSumI;
    //stores the corresponding values.
    vector< vector<double> > rowsPartialSumD;
    computeLocalRowSum(globalizer.getProcIdsInInteraction(),rowsPartialSumI,rowsPartialSumD);
    globalizer.sendSumToLazySide(rowsPartialSumI,rowsPartialSumD);
    globalizer.recvSumFromLazySide(rowsPartialSumD);
    divideByGlobalRowSum(globalizer.getProcIdsInInteraction(),rowsPartialSumI,rowsPartialSumD,denoStrorage);
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
    int id;
    const vector<std::pair<int,int> >& mapping=_mapping.getSendingIds();
    for(vector<std::pair<int,int> >::const_iterator iter2=mapping.begin();iter2!=mapping.end();iter2++)
      {
        std::pair<set<int>::iterator,bool> isIns=procsSet.insert((*iter2).first);
        if(isIns.second)
          id=std::find(distantProcs.begin(),distantProcs.end(),(*iter2).first)-distantProcs.begin();
        resPerProcI[id].push_back((*iter2).second);
        resPerProcD[id].push_back(res[iter2-mapping.begin()]);
      }
    /*
    for(map<pair<int,int>, int >::const_iterator iter2=_col_offsets.begin();iter2!=_col_offsets.end();iter2++)
      {
        std::pair<set<int>::iterator,bool> isIns=procsSet.insert((*iter2).first.first);
        if(isIns.second)
          id=std::find(distantProcs.begin(),distantProcs.end(),(*iter2).first.first)-distantProcs.begin();
        resPerProcI[id].push_back((*iter2).first.second);
        resPerProcD[id].push_back(res[(*iter2).second]);
        }*/
  }

  void InterpolationMatrix::divideByGlobalRowSum(const std::vector<int>& distantProcs, const std::vector<std::vector<int> >& resPerProcI,
                                                 const std::vector<std::vector<double> >& resPerProcD, std::vector<std::vector<double> >& deno)
  {
    map<int,double> fastSums;
    int procId=0;
    const vector<std::pair<int,int> >& mapping=_mapping.getSendingIds();
    map< int, map<int,int> > distIdToLocId;
    for(map< pair<int,int>,int >::iterator iter8=_col_offsets.begin();iter8!=_col_offsets.end();iter8++)
          distIdToLocId[(*iter8).first.first][mapping[(*iter8).second].second]=(*iter8).first.second;

    for(vector<int>::const_iterator iter1=distantProcs.begin();iter1!=distantProcs.end();iter1++,procId++)
      {
        const std::vector<int>& currentProcI=resPerProcI[procId];
        const std::vector<double>& currentProcD=resPerProcD[procId];
        vector<double>::const_iterator iter3=currentProcD.begin();
        for(vector<int>::const_iterator iter2=currentProcI.begin();iter2!=currentProcI.end();iter2++,iter3++)
          fastSums[_col_offsets[std::make_pair(*iter1,distIdToLocId[*iter1][*iter2])]]=*iter3;
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

  // ==================================================================
  // The call to this method updates the arrays on the target side
  //   so that they know which amount of data from which processor they 
  //   should expect. 
  //   That call makes actual interpolations via multiply method 
  //   available.
  // ==================================================================

  void InterpolationMatrix::prepare()
  {
    int nbelems = _source_field->getField()->getNumberOfTuples();
    for (int ielem=0; ielem < nbelems; ielem++)
      {
        _row_offsets[ielem+1]+=_row_offsets[ielem];
      }
    _mapping.prepareSendRecv();
  }


  //   =======================================================================
  //   brief performs t=Ws, where t is the target field, s is the source field

  //   The call to this method must be called both on the working side 
  //   and on the idle side. On the working side, the vector  T=VT^(-1).(W.S)
  //   is computed and sent. On the idle side, no computation is done, but the 
  //   result from the working side is received and the field is updated.

  //   param field source field on processors involved on the source side,
  //   target field on processors on the target side
  //   =======================================================================

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
  

  // =========================================================================
  // brief performs s=WTt, where t is the target field, s is the source field,
  // WT is the transpose matrix from W

  //   The call to this method must be called both on the working side 
  //   and on the idle side. On the working side, the target vector T is
  //   received and the vector  S=VS^(-1).(WT.T) is computed to update
  //   the field. 
  //   On the idle side, no computation is done, but the field is sent.

  //   param field source field on processors involved on the source side,
  //   target field on processors on the target side
  // =========================================================================

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
