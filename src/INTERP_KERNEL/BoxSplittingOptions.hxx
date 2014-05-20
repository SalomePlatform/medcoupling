// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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
#include "NormalizedUnstructuredMesh.hxx"

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
    double _effeciency;
    double _effeciency_snd;
    int _min_cell_direction;
    int _max_cells;
  public:
    BoxSplittingOptions() { init(); }
    void init();
    double getEffeciency() const { return _effeciency; }
    void setEffeciency(double effeciency) { _effeciency=effeciency; }
    double getEffeciencySnd() const { return _effeciency_snd; }
    void setEffeciencySnd(double effeciencySnd) { _effeciency_snd=effeciencySnd; }
    int getMinCellDirection() const { return _min_cell_direction; }
    void setMinCellDirection(int minCellDirection) { _min_cell_direction=minCellDirection; }
    int getMaxCells() const { return _max_cells; }
    void setMaxCells(int maxCells) { _max_cells=maxCells; }
    void copyOptions(const BoxSplittingOptions & other) { *this=other; }
    std::string printOptions() const;
  private:
    static const int DFT_MIN_CELL_DIRECTION=3;
    static const int DFT_MAX_CELLS=1000;
    static const double DFT_EFFECIENCY;
    static const double DFT_EFFECIENCY_SND;
  public:
    static const char PRECISION_STR[];
  };
}

#endif
