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
// Author : Anthony Geay

#ifndef __BOXSPLITTINGOPTIONS_HXX__
#define __BOXSPLITTINGOPTIONS_HXX__

#include "INTERPKERNELDefines.hxx"

#include <string>

namespace INTERP_KERNEL
{ 
  /*!
   * \class BoxSplittingOptions
   * Class defining the options for box splitting used for AMR algorithm like creation of patches following a criterion.
   */
  class INTERPKERNEL_EXPORT BoxSplittingOptions
  {
  private:
    double _efficiency_goal;
    double _efficiency_threshold;
    int _min_patch_length;
    int _max_patch_length;
    int _max_nb_cells_in_patch;
  public:
    BoxSplittingOptions() { init(); }
    void init();
    double getEfficiencyGoal() const { return _efficiency_goal; }
    void setEfficiencyGoal(double efficiency) { _efficiency_goal=efficiency; }
    double getEfficiencyThreshold() const { return _efficiency_threshold; }
    void setEfficiencyThreshold(double efficiencyThreshold) { _efficiency_threshold=efficiencyThreshold; }
    int getMinimumPatchLength() const { return _min_patch_length; }
    void setMinimumPatchLength(int minPatchLength) { _min_patch_length=minPatchLength; }
    int getMaximumPatchLength() const { return _max_patch_length; }
    void setMaximumPatchLength(int maxNbCellPatch) { _max_patch_length=maxNbCellPatch; }
    int getMaximumNbOfCellsInPatch() const { return _max_nb_cells_in_patch; }
    void setMaximumNbOfCellsInPatch(int maxNbCellsInPatch) { _max_nb_cells_in_patch=maxNbCellsInPatch; }
    void copyOptions(const BoxSplittingOptions & other) { *this=other; }
    std::string printOptions() const;
  private:
    static const int DFT_MIN_PATCH_LENGTH=1;
    static const int DFT_MAX_PATCH_LENGTH=2147483647;
    static const int DFT_MAX_PATCH_MEASURE=2147483647;
    static const double DFT_EFFICIENCY_GOAL;
    static const double DFT_EFFICIENCY_THRESHOLD;
  public:
    static const char PRECISION_STR[];
  };
}

#endif
