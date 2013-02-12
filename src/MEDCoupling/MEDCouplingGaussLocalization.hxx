// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#ifndef __PARAMEDMEM_MEDCOUPLINGGAUSSLOCALIZATION_HXX__
#define __PARAMEDMEM_MEDCOUPLINGGAUSSLOCALIZATION_HXX__

#include "MEDCoupling.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpKernelException.hxx"

#include <vector>

namespace ParaMEDMEM
{
  class MEDCouplingMesh;

  class MEDCOUPLING_EXPORT MEDCouplingGaussLocalization
  {
  public:
    MEDCouplingGaussLocalization(INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                 const std::vector<double>& gsCoo, const std::vector<double>& w) throw(INTERP_KERNEL::Exception);
    INTERP_KERNEL::NormalizedCellType getType() const { return _type; }
    int getNumberOfGaussPt() const { return (int)_weight.size(); }
    int getDimension() const;
    int getNumberOfPtsInRefCell() const;
    std::string getStringRepr() const;
    std::size_t getHeapMemorySize() const;
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    bool isEqual(const MEDCouplingGaussLocalization& other, double eps) const;
    void pushTinySerializationIntInfo(std::vector<int>& tinyInfo) const;
    void pushTinySerializationDblInfo(std::vector<double>& tinyInfo) const;
    const double *fillWithValues(const double *vals);
    //
    const std::vector<double>& getRefCoords() const { return _ref_coord; }
    double getRefCoord(int ptIdInCell, int comp) const throw(INTERP_KERNEL::Exception);
    const std::vector<double>& getGaussCoords() const { return _gauss_coord; }
    double getGaussCoord(int gaussPtIdInCell, int comp) const throw(INTERP_KERNEL::Exception);
    const std::vector<double>& getWeights() const { return _weight; }
    double getWeight(int gaussPtIdInCell, double newVal) const throw(INTERP_KERNEL::Exception);
    void setRefCoord(int ptIdInCell, int comp, double newVal) throw(INTERP_KERNEL::Exception);
    void setGaussCoord(int gaussPtIdInCell, int comp, double newVal) throw(INTERP_KERNEL::Exception);
    void setWeight(int gaussPtIdInCell, double newVal) throw(INTERP_KERNEL::Exception);
    //
    static MEDCouplingGaussLocalization BuildNewInstanceFromTinyInfo(int dim, const std::vector<int>& tinyData);
    static bool AreAlmostEqual(const std::vector<double>& v1, const std::vector<double>& v2, double eps);
  private:
    int checkCoherencyOfRequest(int gaussPtIdInCell, int comp) const throw(INTERP_KERNEL::Exception);
  private:
    INTERP_KERNEL::NormalizedCellType _type;
    std::vector<double> _ref_coord;
    std::vector<double> _gauss_coord;
    std::vector<double> _weight;
  };
}

#endif
