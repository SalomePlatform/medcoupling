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

#include "BoxSplittingOptions.hxx"

#include <sstream>

const double INTERP_KERNEL::BoxSplittingOptions::DFT_EFFECIENCY=0.5;

const double INTERP_KERNEL::BoxSplittingOptions::DFT_EFFECIENCY_SND=0.5;

void INTERP_KERNEL::BoxSplittingOptions::init()
{
  _effeciency=DFT_EFFECIENCY;
  _effeciency_snd=DFT_EFFECIENCY_SND;
  _min_cell_direction=DFT_MIN_CELL_DIRECTION;
  _max_cells=DFT_MAX_CELLS;
}

std::string INTERP_KERNEL::BoxSplittingOptions::printOptions() const
{
  std::ostringstream oss;
  oss << "Efficiency threshold: " << 100*_effeciency << "%"<< std::endl;
  oss << "Efficiency Snd threshold: " << 100*_effeciency_snd << "%"<< std::endl;
  oss << "Min. box side length: " << _min_cell_direction << std::endl;
  oss << "Max. box cells : " << _max_cells << std::endl;
  return oss.str();
}
