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

  InterpolationMatrix::InterpolationMatrix(ParaMEDMEM::ParaMESH *source_support, 
                                           const ProcessorGroup& source_group,
                                           const ProcessorGroup& target_group,
                                           const DECOptions& dec_options,
                                           const INTERP_KERNEL::InterpolationOptions& interp_options):
    _source_support(source_support->getCellMesh()),
    _mapping(source_group, target_group, dec_options),
    _source_group(source_group),
    _target_group(target_group),
    _source_volume(0),
    DECOptions(dec_options),
    INTERP_KERNEL::InterpolationOptions(interp_options)
  {
    int nbelems = _source_support->getNumberOfCells();

    _row_offsets.resize(nbelems+1);
    for (int i=0; i<nbelems+1; i++)
      {
        _row_offsets[i]=0;
      }

    _coeffs.resize(nbelems);
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
    std::string interpMethod(srcMeth);
    interpMethod+=targetMeth;
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
              && distant_support.getSpaceDimension() ==3 )
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

    MEDCouplingFieldDouble *target_triangle_surf = distant_support.getMeasureField();
    MEDCouplingFieldDouble *source_triangle_surf = _source_support->getMeasureField();

    // Storing the source volumes
    _source_volume.resize(source_size);
    for (int i=0; i<source_size; i++)
      {
        _source_volume[i] = source_triangle_surf->getIJ(i,0);
      }

    source_triangle_surf->decrRef();

    //loop over the elements to build the interpolation
    //matrix structures
    for (int ielem=0; ielem < surfaces.size(); ielem++) 
      {
        _row_offsets[ielem+1] += surfaces[ielem].size();
        //_source_indices.push_back(make_pair(iproc_distant, ielem));

        for (map<int,double>::const_iterator iter = surfaces[ielem].begin();
             iter != surfaces[ielem].end();
             iter++)
          {
            //surface of the target triangle
            double surf = target_triangle_surf->getIJ(iter->first,0);

            //locating the (iproc, itriangle) pair in the list of columns
            map<pair<int,int>,int >::iterator iter2 = _col_offsets.find(make_pair(iproc_distant,iter->first));
            int col_id;

            if (iter2 == _col_offsets.end())
              {
                //(iproc, itriangle) is not registered in the list
                //of distant elements

                col_id =_col_offsets.size()+1;
                _col_offsets.insert(make_pair(make_pair(iproc_distant,iter->first),col_id));
                _mapping.addElementFromSource(iproc_distant,
                                              distant_elems[iter->first]);
                _target_volume.push_back(surf);
              }
            else 
              {
                col_id = iter2->second;
              }

            //the non zero coefficient is stored 
            //ielem is the row,
            //col_id is the number of the column
            //iter->second is the value of the coefficient

            _coeffs[ielem].push_back(make_pair(col_id,iter->second));
          }
      }
    target_triangle_surf->decrRef();
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
    int nbelems = _source_support->getNumberOfCells();
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
        int nbrows = _source_support->getNumberOfCells();

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
                    target_value[(colid-1)*nbcomp+icomp]+=value*coeff_row;
                  }
              }
          }

        // performing VT^(-1).(W.S)
        // where VT^(-1) is the inverse of the diagonal matrix containing 
        // the volumes of target cells

        for (int i=0; i<_col_offsets.size();i++)
          {
            for (int icomp=0; icomp<nbcomp; icomp++)
              {
                target_value[i*nbcomp+icomp] /= _target_volume[i];
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
    //  int nbcomp = field.getNumberOfComponents();
    int nbcomp = field.getArray()->getNumberOfComponents();
    vector<double> source_value(_col_offsets.size()* nbcomp,0.0);
    _mapping.reverseSendRecv(&source_value[0],field);

    //treatment of the transpose matrix multiply on the source side
    if (_source_group.containsMyRank())
      {
        int nbrows    = _source_support->getNumberOfCells();
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
          
                for (int icomp=0; icomp<nbcomp; icomp++)
                  {
                    double coeff_row = source_value[(colid-1)*nbcomp+icomp];
                    array[irow*nbcomp+icomp] += value*coeff_row;
                  }
              }
          }

        //performing VS^(-1).(WT.T)
        //VS^(-1) is the inverse of the diagonal matrix storing
        //volumes of the source cells

        for (int irow=0; irow<nbrows; irow++)
          {
            for (int icomp=0; icomp<nbcomp; icomp++)
              {
                array[irow*nbcomp+icomp] /= _source_volume[irow];
              }
          }

      }
  }
}
