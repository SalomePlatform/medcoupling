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
// Author : Anthony Geay (EDF R&D)

#include "OrientationInverter.hxx"
#include "InterpKernelException.hxx"

#include <sstream>
#include <algorithm>

using namespace INTERP_KERNEL;

void OrientationInverterChecker::check(int *beginPt, int *endPt) const
{
  if(std::distance(beginPt,endPt)!=getNbNodes())
    {
      std::ostringstream oss; oss << "OrientationInverterChecker::check : length of nodal connectivity mismatches ! Expecting " << getNbNodes() << " having " << std::distance(beginPt,endPt) << " !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
}

void OrientationInverterSEG2::operateAndShutUp(int *beginPt) const
{
  std::swap(beginPt[0],beginPt[1]);
}

void OrientationInverterSEG3::operateAndShutUp(int *beginPt) const
{
  std::swap(beginPt[0],beginPt[2]);
}

void OrientationInverter2DLinear::operateAndShutUp(int *beginPt) const
{
  std::reverse(beginPt+1,beginPt+getNbNodes());
}

void OrientationInverter2DQuadratic::operateAndShutUp(int *beginPt) const
{
  int nbNodes(getNbNodes());
  std::reverse(beginPt+1,beginPt+nbNodes/2);
  std::reverse(beginPt+nbNodes/2,beginPt+nbNodes);
}

void OrientationInverterPolygon::operate(int *beginPt, int *endPt) const
{
  std::reverse(beginPt+1,endPt);
}

void OrientationInverterQPolygon::operate(int *beginPt, int *endPt) const
{
  std::size_t sz(std::distance(beginPt,endPt));
  std::reverse(beginPt+1,beginPt+sz/2);
  std::reverse(beginPt+sz/2,endPt);
}

void OrientationInverterTetra4::operateAndShutUp(int *beginPt) const
{
  std::swap(beginPt[1],beginPt[2]);
}

void OrientationInverterTetra10::operateAndShutUp(int *beginPt) const
{
  std::swap(beginPt[1],beginPt[2]);
  std::swap(beginPt[4],beginPt[6]);
  std::swap(beginPt[8],beginPt[9]);
}

void OrientationInverterPyra5::operateAndShutUp(int *beginPt) const
{
  std::reverse(beginPt+1,beginPt+4);
}

void OrientationInverterPyra13::operateAndShutUp(int *beginPt) const
{
  std::reverse(beginPt+1,beginPt+4);
  std::reverse(beginPt+5,beginPt+9);
  std::swap(beginPt[10],beginPt[12]);
}

void OrientationInverter3DExtrusionLinear::operateAndShutUp(int *beginPt) const
{
  int nbNodes(getNbNodes());
  std::reverse(beginPt+1,beginPt+nbNodes/2);
  std::reverse(beginPt+nbNodes/2+1,beginPt+nbNodes);
}

void OrientationInverter3DExtrusionQuadratic::operateAndShutUp(int *beginPt) const
{
  int nbNodes(getNbNodes()),nbNodesLinearBase(nbNodes/5);
  std::reverse(beginPt+1,beginPt+nbNodesLinearBase);
  std::reverse(beginPt+nbNodesLinearBase+1,beginPt+2*nbNodesLinearBase);
  std::reverse(beginPt+2*nbNodesLinearBase,beginPt+3*nbNodesLinearBase);
  std::reverse(beginPt+3*nbNodesLinearBase,beginPt+4*nbNodesLinearBase);
  std::reverse(beginPt+4*nbNodesLinearBase+1,beginPt+5*nbNodesLinearBase);
}
