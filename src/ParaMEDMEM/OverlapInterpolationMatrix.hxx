// Copyright (C) 2007-2020  CEA/DEN, EDF R&D
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

#ifndef __OVERLAPINTERPOLATIONMATRIX_HXX__
#define __OVERLAPINTERPOLATIONMATRIX_HXX__

#include "MPIAccessDEC.hxx"
#include "OverlapMapping.hxx"
#include "InterpolationOptions.hxx"
#include "DECOptions.hxx"

namespace MEDCoupling
{
  class ParaFIELD;
  class MEDCouplingPointSet;

  /*!
   * Internal class, not part of the public API.
   *
   * Similar to InterpolationMatrix, but for the OverlapDEC instead of the InterpKernelDEC.
   */
  class OverlapInterpolationMatrix : public INTERP_KERNEL::InterpolationOptions,
                                     public DECOptions
  {
  public:
    
    OverlapInterpolationMatrix(ParaFIELD *source_field,
                               ParaFIELD *target_field,
                               const ProcessorGroup& group,
                               const DECOptions& dec_opt,
                               const InterpolationOptions& i_opt,
                               const OverlapElementLocator & loc);

    void keepTracksOfSourceIds(int procId, DataArrayIdType *ids);

    void keepTracksOfTargetIds(int procId, DataArrayIdType *ids);

    void computeLocalIntersection(const MEDCouplingPointSet *src, const DataArrayIdType *srcIds, const std::string& srcMeth, int srcProcId,
                         const MEDCouplingPointSet *trg, const DataArrayIdType *trgIds, const std::string& trgMeth, int trgProcId);

    void prepare(const std::vector< int > & procsToSendField);
    
    void computeSurfacesAndDeno();

    void multiply(double default_val);

    void transposeMultiply();
    
    virtual ~OverlapInterpolationMatrix();
  private:

    static void TransposeMatrix(const std::vector<SparseDoubleVec>& matIn, mcIdType nbColsMatIn,
                                std::vector<SparseDoubleVec>& matOut);
  private:
    ParaFIELD           *_source_field;
    ParaFIELD           *_target_field;
    MEDCouplingPointSet *_source_support;
    MEDCouplingPointSet *_target_support;
    OverlapMapping      _mapping;
  };
}

#endif
