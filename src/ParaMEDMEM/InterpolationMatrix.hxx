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

#ifndef __INTERPOLATIONMATRIX_HXX__
#define __INTERPOLATIONMATRIX_HXX__

#include "MPIAccessDEC.hxx"
#include "MxN_Mapping.hxx"
#include "InterpolationOptions.hxx"
#include "DECOptions.hxx"

namespace MEDCoupling
{
  class ElementLocator;

  /*!
   Internal class, not part of the public API.

   This class enables the storage of an interpolation matrix Wij mapping
   a source field Sj to a target field Ti via Ti=Vi^(-1).Wij.Sj.
   The matrix is built and stored on the processors belonging to the source
   group.
   */
  class InterpolationMatrix : public INTERP_KERNEL::InterpolationOptions,
                              public DECOptions
  {
  public:
    
    InterpolationMatrix(const MEDCoupling::ParaFIELD *source_field, 
                        const ProcessorGroup& source_group,
                        const ProcessorGroup& target_group,
                        const DECOptions& dec_opt,
                        const InterpolationOptions& i_opt);

    
    virtual ~InterpolationMatrix();
    void addContribution(MEDCouplingPointSet& distant_support, int iproc_distant,
                         const int* distant_elems, const std::string& srcMeth, const std::string& targetMeth);
    void finishContributionW(ElementLocator& elementLocator);
    void finishContributionL(ElementLocator& elementLocator);
    void multiply(MEDCouplingFieldDouble& field) const;
    void transposeMultiply(MEDCouplingFieldDouble& field)const;
    void prepare();
    int getNbRows() const { return _row_offsets.size(); }
    MPIAccessDEC* getAccessDEC() { return _mapping.getAccessDEC(); }
  private:
    void computeConservVolDenoW(ElementLocator& elementLocator);
    void computeIntegralDenoW(ElementLocator& elementLocator);
    void computeRevIntegralDenoW(ElementLocator& elementLocator);
    void computeGlobConstraintDenoW(ElementLocator& elementLocator);
    void computeConservVolDenoL(ElementLocator& elementLocator);
    void computeIntegralDenoL(ElementLocator& elementLocator);
    void computeRevIntegralDenoL(ElementLocator& elementLocator);
    
    void computeLocalColSum(std::vector<double>& res) const;
    void computeLocalRowSum(const std::vector<int>& distantProcs, std::vector<std::vector<int> >& resPerProcI,
                            std::vector<std::vector<double> >& resPerProcD) const;
    void computeGlobalRowSum(ElementLocator& elementLocator, std::vector<std::vector<double> >& denoStrorage, std::vector<std::vector<double> >& denoStrorageInv);
    void computeGlobalColSum(std::vector<std::vector<double> >& denoStrorage);
    void resizeGlobalColSum(std::vector<std::vector<double> >& denoStrorage);
    void fillDSFromVM(int iproc_distant, const int* distant_elems, const std::vector< std::map<int,double> >& values, MEDCouplingFieldDouble *surf);
    void serializeMe(std::vector< std::vector< std::map<int,double> > >& data1, std::vector<int>& data2) const;
    void initialize();
    void findAdditionnalElements(ElementLocator& elementLocator, std::vector<std::vector<int> >& elementsToAdd,
                                 const std::vector<std::vector<int> >& resPerProcI, const std::vector<std::vector<int> >& globalIdsPartial);
    void addGhostElements(const std::vector<int>& distantProcs, const std::vector<std::vector<int> >& elementsToAdd);
    int mergePolicies(const std::vector<int>& policyPartial);
    void mergeRowSum(const std::vector< std::vector<double> >& rowsPartialSumD, const std::vector< std::vector<int> >& globalIdsPartial,
                     std::vector<int>& globalIdsLazySideInteraction, std::vector<double>& sumCorresponding);
    void mergeRowSum2(const std::vector< std::vector<int> >& globalIdsPartial, std::vector< std::vector<double> >& rowsPartialSumD,
                      const std::vector<int>& globalIdsLazySideInteraction, const std::vector<double>& sumCorresponding);
    void mergeRowSum3(const std::vector< std::vector<int> >& globalIdsPartial, std::vector< std::vector<double> >& rowsPartialSumD);
    void mergeCoeffs(const std::vector<int>& procsInInteraction, const std::vector< std::vector<int> >& rowsPartialSumI,
                     const std::vector<std::vector<int> >& globalIdsPartial, std::vector<std::vector<double> >& denoStrorageInv);
    void divideByGlobalRowSum(const std::vector<int>& distantProcs, const std::vector<std::vector<int> >& resPerProcI,
                              const std::vector<std::vector<double> >& resPerProcD, std::vector<std::vector<double> >& deno);
  private:
    bool isSurfaceComputationNeeded(const std::string& method) const;
  private:
    const MEDCoupling::ParaFIELD *_source_field;
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
