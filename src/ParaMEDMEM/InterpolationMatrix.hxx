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
#ifndef __INTERPOLATIONMATRIX_HXX__
#define __INTERPOLATIONMATRIX_HXX__

#include "MPIAccessDEC.hxx"
#include "MxN_Mapping.hxx"
#include "InterpolationOptions.hxx"
#include "DECOptions.hxx"

namespace ParaMEDMEM
{
  class GlobalizerMeshWorkingSide;
  class GlobalizerMeshLazySide;

  class InterpolationMatrix : public INTERP_KERNEL::InterpolationOptions,
                              public DECOptions
  {
  public:
    
    InterpolationMatrix(const ParaMEDMEM::ParaFIELD *source_field, 
                        const ProcessorGroup& source_group,
                        const ProcessorGroup& target_group,
                        const DECOptions& dec_opt,
                        const InterpolationOptions& i_opt);

    
    virtual ~InterpolationMatrix();
    void addContribution(MEDCouplingPointSet& distant_support, int iproc_distant,
                         int* distant_elems, const std::string& srcMeth, const std::string& targetMeth);
    void finishContributionW(GlobalizerMeshWorkingSide& globalizer);
    void finishContributionL(GlobalizerMeshLazySide& globalizer);
    void multiply(MEDCouplingFieldDouble& field) const;
    void transposeMultiply(MEDCouplingFieldDouble& field)const;
    void prepare();
    int getNbRows() const { return _row_offsets.size(); }
    MPIAccessDEC* getAccessDEC() { return _mapping.getAccessDEC(); }
  private:
    void computeConservVolDenoW(GlobalizerMeshWorkingSide& globalizer);
    void computeIntegralDenoW(GlobalizerMeshWorkingSide& globalizer);
    void computeGlobConstraintDenoW(GlobalizerMeshWorkingSide& globalizer);
    void computeConservVolDenoL(GlobalizerMeshLazySide& globalizer);
    void computeIntegralDenoL(GlobalizerMeshLazySide& globalizer);
    void computeGlobConstraintDenoL(GlobalizerMeshLazySide& globalizer);
    
    void computeLocalColSum(std::vector<double>& res) const;
    void computeLocalRowSum(const std::vector<int>& distantProcs, std::vector<std::vector<int> >& resPerProcI,
                            std::vector<std::vector<double> >& resPerProcD) const;
    void computeGlobalRowSum(GlobalizerMeshWorkingSide& globalizer, std::vector<std::vector<double> >& denoStrorage);
    void computeGlobalColSum(std::vector<std::vector<double> >& denoStrorage);
    void divideByGlobalRowSum(const std::vector<int>& distantProcs, const std::vector<std::vector<int> >& resPerProcI,
                              const std::vector<std::vector<double> >& resPerProcD, std::vector<std::vector<double> >& deno);
  private:
    bool isSurfaceComputationNeeded(const std::string& method) const;
  private:
    const ParaMEDMEM::ParaFIELD *_source_field;
    std::vector<int> _row_offsets;
    std::map<std::pair<int,int>, int > _col_offsets;
    MEDCouplingPointSet *_source_support;
    MxN_Mapping _mapping;
 
    const ProcessorGroup& _source_group;
    const ProcessorGroup& _target_group;
    std::vector< std::vector<double> > _target_volume;
    std::vector<std::vector<std::pair<int,double> > > _coeffs;
    std::vector<std::vector<double> > _deno_multiply;
    std::vector<std::vector<double> > _deno_reverse_multiply;
  };
}

#endif
