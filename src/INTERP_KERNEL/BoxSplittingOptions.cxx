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

#include "BoxSplittingOptions.hxx"

#include <sstream>

const double INTERP_KERNEL::BoxSplittingOptions::DFT_EFFICIENCY_GOAL=0.5;

const double INTERP_KERNEL::BoxSplittingOptions::DFT_EFFICIENCY_THRESHOLD=0.7;

void INTERP_KERNEL::BoxSplittingOptions::init()
{
  _efficiency_goal=DFT_EFFICIENCY_GOAL;
  _efficiency_threshold=DFT_EFFICIENCY_THRESHOLD;
  _min_patch_length=DFT_MIN_PATCH_LENGTH;
  _max_patch_length=DFT_MAX_PATCH_LENGTH;
  _max_nb_cells_in_patch=DFT_MAX_PATCH_MEASURE;
}

std::string INTERP_KERNEL::BoxSplittingOptions::printOptions() const
{
  std::ostringstream oss;
  oss << "Efficiency goal: " << 100*_efficiency_goal << "%" << std::endl;
  oss << "Efficiency threshold: " << 100*_efficiency_threshold << "%" << std::endl;
  oss << "Min. patch side length: " << _min_patch_length << std::endl;
  oss << "Max. patch side length: " << _max_patch_length << std::endl;
  oss << "Max. patch measure: " << _max_nb_cells_in_patch << std::endl;
  return oss.str();
}
