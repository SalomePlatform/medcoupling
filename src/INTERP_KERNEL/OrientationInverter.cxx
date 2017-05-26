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
#include "CellModel.hxx"

#include <sstream>
#include <algorithm>

using namespace INTERP_KERNEL;

OrientationInverter *OrientationInverter::BuildInstanceFrom(NormalizedCellType gt)
{
  switch(gt)
    {
    case NORM_SEG2:
      return new OrientationInverterSEG2;
    case NORM_SEG3:
      return new OrientationInverterSEG3;
    case NORM_TRI3:
      return new OrientationInverter2DLinear(3u);
    case NORM_QUAD4:
      return new OrientationInverter2DLinear(4u);
    case NORM_POLYGON:
      return new OrientationInverterPolygon;
    case NORM_TRI6:
      return new OrientationInverter2DQuadratic(6u);
    case NORM_QUAD8:
      return new OrientationInverter2DQuadratic(8u);
    case NORM_QPOLYG:
      return new OrientationInverterQPolygon;
    case NORM_TETRA4:
      return new OrientationInverterTetra4;
    case NORM_PYRA5:
      return new OrientationInverterPyra5;
    case NORM_PENTA6:
      return new OrientationInverter3DExtrusionLinear(6u);
    case NORM_HEXA8:
      return new OrientationInverter3DExtrusionLinear(8u);
    case NORM_TETRA10:
      return new OrientationInverterTetra10;
    case NORM_PYRA13:
      return new OrientationInverterPyra13;
    case NORM_PENTA15:
      return new OrientationInverter3DExtrusionQuadratic(15u);
    case NORM_HEXA20:
      return new OrientationInverter3DExtrusionQuadratic(20u);
    default:
      {
        const CellModel& cm(CellModel::GetCellModel(gt));
        std::ostringstream oss; oss << "OrientationInverter::BuildInstanceFrom : Sorry no inverter for geo type " << cm.getRepr() << " !";
        throw INTERP_KERNEL::Exception(oss.str());
      }
    }
}

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
