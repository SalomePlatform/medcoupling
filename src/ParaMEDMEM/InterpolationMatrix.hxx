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
  class InterpolationMatrix : public INTERP_KERNEL::InterpolationOptions,
                              public DECOptions
  {
  public:
    
    InterpolationMatrix(ParaMEDMEM::ParaMESH *source_support, 
                        const ProcessorGroup& source_group,
                        const ProcessorGroup& target_group,
                        const DECOptions& dec_opt,
                        const InterpolationOptions& i_opt);

    
    virtual ~InterpolationMatrix();
    void addContribution(MEDCouplingUMesh& distant_support, int iproc_distant,
                         int* distant_elems, const std::string& srcMeth, const std::string& targetMeth);
    void multiply(MEDCouplingFieldDouble& field) const;
    void transposeMultiply(MEDCouplingFieldDouble& field)const;
    void prepare();
    int getNbRows() const {return _row_offsets.size();}
    MPIAccessDEC* getAccessDEC(){return _mapping.getAccessDEC();}

    static MEDCouplingFieldDouble *getSupportVolumes(MEDCouplingMesh *field);
    static MEDCouplingFieldDouble *getSupportUnstructuredVolumes(MEDCouplingUMesh *field);
  private:
    std::vector<int> _row_offsets;
    std::vector<std::pair<int,int> > _col_offsets;
    MEDCouplingUMesh *_source_support; 
    MxN_Mapping _mapping;
 
    const ProcessorGroup& _source_group;
    const ProcessorGroup& _target_group;
    std::vector<double> _target_volume;
    std::vector<double> _source_volume;
    std::vector<std::vector<std::pair<int,double> > > _coeffs;
  };
}

#endif
