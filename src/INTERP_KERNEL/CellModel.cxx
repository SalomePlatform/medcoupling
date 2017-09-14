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

#include "CellModel.hxx"

#include "InterpKernelException.hxx"
#include "DiameterCalculator.hxx"
#include "OrientationInverter.hxx"

#include <algorithm>
#include <sstream>
#include <vector>
#include <limits>

namespace INTERP_KERNEL
{
  const char *CellModel::CELL_TYPES_REPR[]={"NORM_POINT1", "NORM_SEG2", "NORM_SEG3", "NORM_TRI3", "NORM_QUAD4",// 0->4
                                            "NORM_POLYGON", "NORM_TRI6", "NORM_TRI7" , "NORM_QUAD8", "NORM_QUAD9",//5->9
                                            "NORM_SEG4", "", "", "", "NORM_TETRA4",//10->14
                                            "NORM_PYRA5", "NORM_PENTA6", "", "NORM_HEXA8", "",//15->19
                                            "NORM_TETRA10", "", "NORM_HEXGP12", "NORM_PYRA13", "",//20->24
                                            "NORM_PENTA15", "", "NORM_HEXA27", "NORM_PENTA18", "",//25->29
                                            "NORM_HEXA20", "NORM_POLYHED", "NORM_QPOLYG", "NORM_POLYL", "",//30->34
                                            "", "", "", "", "",//35->39
                                            "NORM_ERROR"};

  std::map<NormalizedCellType,CellModel> CellModel::_map_of_unique_instance;

  const CellModel& CellModel::GetCellModel(NormalizedCellType type)
  {
    if(_map_of_unique_instance.empty())
      buildUniqueInstance();
    const std::map<NormalizedCellType,CellModel>::iterator iter=_map_of_unique_instance.find(type);
    if(iter==_map_of_unique_instance.end())
      {
        std::ostringstream stream; stream << "no cellmodel for normalized type " << type;
        throw Exception(stream.str().c_str());
      }
    return (*iter).second;
  }

  const char *CellModel::getRepr() const
  {
    return CELL_TYPES_REPR[(int)_type];
  }

  /*!
   * This method is compatible with all types including dynamic one.
   */
  bool CellModel::isCompatibleWith(NormalizedCellType type) const
  {
    if(_type==type)
      return true;
    const CellModel& other=GetCellModel(type);
    if(_dim!=other.getDimension())
      return false;
    bool b1=isQuadratic();
    bool b2=other.isQuadratic();
    if((b1 && !b2) || (!b1 && b2))
      return false;
    b1=isDynamic();
    b2=other.isDynamic();
    return b1 || b2;
  }

  void CellModel::buildUniqueInstance()
  {
    _map_of_unique_instance.insert(std::make_pair(NORM_POINT1,CellModel(NORM_POINT1)));
    _map_of_unique_instance.insert(std::make_pair(NORM_SEG2,CellModel(NORM_SEG2)));
    _map_of_unique_instance.insert(std::make_pair(NORM_SEG3,CellModel(NORM_SEG3)));
    _map_of_unique_instance.insert(std::make_pair(NORM_SEG4,CellModel(NORM_SEG4)));
    _map_of_unique_instance.insert(std::make_pair(NORM_TRI3,CellModel(NORM_TRI3)));
    _map_of_unique_instance.insert(std::make_pair(NORM_QUAD4,CellModel(NORM_QUAD4)));
    _map_of_unique_instance.insert(std::make_pair(NORM_TRI6,CellModel(NORM_TRI6)));
    _map_of_unique_instance.insert(std::make_pair(NORM_TRI7,CellModel(NORM_TRI7)));
    _map_of_unique_instance.insert(std::make_pair(NORM_QUAD8,CellModel(NORM_QUAD8)));
    _map_of_unique_instance.insert(std::make_pair(NORM_QUAD9,CellModel(NORM_QUAD9)));
    _map_of_unique_instance.insert(std::make_pair(NORM_TETRA4,CellModel(NORM_TETRA4)));
    _map_of_unique_instance.insert(std::make_pair(NORM_HEXA8,CellModel(NORM_HEXA8)));
    _map_of_unique_instance.insert(std::make_pair(NORM_PYRA5,CellModel(NORM_PYRA5)));
    _map_of_unique_instance.insert(std::make_pair(NORM_PENTA6,CellModel(NORM_PENTA6)));
    _map_of_unique_instance.insert(std::make_pair(NORM_TETRA10,CellModel(NORM_TETRA10)));
    _map_of_unique_instance.insert(std::make_pair(NORM_HEXGP12,CellModel(NORM_HEXGP12)));
    _map_of_unique_instance.insert(std::make_pair(NORM_PYRA13,CellModel(NORM_PYRA13)));
    _map_of_unique_instance.insert(std::make_pair(NORM_PENTA15,CellModel(NORM_PENTA15)));
    _map_of_unique_instance.insert(std::make_pair(NORM_PENTA18,CellModel(NORM_PENTA18)));
    _map_of_unique_instance.insert(std::make_pair(NORM_HEXA20,CellModel(NORM_HEXA20)));
    _map_of_unique_instance.insert(std::make_pair(NORM_HEXA27,CellModel(NORM_HEXA27)));
    _map_of_unique_instance.insert(std::make_pair(NORM_POLYGON,CellModel(NORM_POLYGON)));
    _map_of_unique_instance.insert(std::make_pair(NORM_POLYHED,CellModel(NORM_POLYHED)));
    _map_of_unique_instance.insert(std::make_pair(NORM_QPOLYG,CellModel(NORM_QPOLYG)));
    _map_of_unique_instance.insert(std::make_pair(NORM_POLYL,CellModel(NORM_POLYL)));
    _map_of_unique_instance.insert(std::make_pair(NORM_ERROR,CellModel(NORM_ERROR)));
  }

  CellModel::CellModel(NormalizedCellType type):_type(type)
  {
    _is_extruded=false;
    _quadratic=false;
    _dyn=false;
    _extruded_type=NORM_ERROR;
    _reverse_extruded_type=NORM_ERROR;
    _linear_type=NORM_ERROR;
    _quadratic_type=NORM_ERROR;
    _quadratic_type2=NORM_ERROR;
    _nb_of_little_sons=std::numeric_limits<unsigned>::max();
    switch(type)
      {
      case NORM_POINT1:
        {
          _nb_of_pts=1; _nb_of_sons=0; _dim=0; _extruded_type=NORM_SEG2; _is_simplex=true;
        }
        break;
      case NORM_SEG2:
        {
          _nb_of_pts=2; _nb_of_sons=2; _dim=1; _extruded_type=NORM_QUAD4; _quadratic_type=NORM_SEG3; _quadratic_type2=NORM_SEG3; _is_simplex=true; _is_extruded=true; _reverse_extruded_type=NORM_POINT1;
          _sons_type[0]=NORM_POINT1; _sons_type[1]=NORM_POINT1;
          _sons_con[0][0]=0; _nb_of_sons_con[0]=1;
          _sons_con[1][0]=1; _nb_of_sons_con[1]=1;
        }
        break;
      case NORM_SEG3:
        {
          _nb_of_pts=3; _nb_of_sons=3; _dim=1; _extruded_type=NORM_QUAD8; _linear_type=NORM_SEG2; _quadratic=true; _is_simplex=false;
          _sons_type[0]=NORM_POINT1; _sons_type[1]=NORM_POINT1; _sons_type[2]=NORM_POINT1;
          _sons_con[0][0]=0; _nb_of_sons_con[0]=1;
          _sons_con[1][0]=1; _nb_of_sons_con[1]=1;
          _sons_con[2][0]=2; _nb_of_sons_con[2]=1;
        }
        break;
      case NORM_SEG4:
        {
          _nb_of_pts=4; _nb_of_sons=4; _dim=1; _linear_type=NORM_SEG2; _quadratic=true; _is_simplex=false; // no _extruded_type because no cubic 2D cell
          _sons_type[0]=NORM_POINT1; _sons_type[1]=NORM_POINT1; _sons_type[2]=NORM_POINT1; _sons_type[3]=NORM_POINT1;
          _sons_con[0][0]=0; _nb_of_sons_con[0]=1;
          _sons_con[1][0]=1; _nb_of_sons_con[1]=1;
          _sons_con[2][0]=2; _nb_of_sons_con[2]=1;
          _sons_con[3][0]=3; _nb_of_sons_con[3]=1;
        }
        break;
      case NORM_TETRA4:
        {
          _nb_of_pts=4; _nb_of_sons=4; _dim=3; _quadratic_type=NORM_TETRA10; _is_simplex=true;
          _sons_type[0]=NORM_TRI3; _sons_type[1]=NORM_TRI3; _sons_type[2]=NORM_TRI3; _sons_type[3]=NORM_TRI3;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _nb_of_sons_con[0]=3;
          _sons_con[1][0]=0; _sons_con[1][1]=3; _sons_con[1][2]=1; _nb_of_sons_con[1]=3;
          _sons_con[2][0]=1; _sons_con[2][1]=3; _sons_con[2][2]=2; _nb_of_sons_con[2]=3;
          _sons_con[3][0]=2; _sons_con[3][1]=3; _sons_con[3][2]=0; _nb_of_sons_con[3]=3;
          _little_sons_con[0][0]=0; _little_sons_con[0][1]=1;  _nb_of_little_sons=6;
          _little_sons_con[1][0]=1; _little_sons_con[1][1]=2;
          _little_sons_con[2][0]=2; _little_sons_con[2][1]=0;
          _little_sons_con[3][0]=0; _little_sons_con[3][1]=3;
          _little_sons_con[4][0]=1; _little_sons_con[4][1]=3;
          _little_sons_con[5][0]=2; _little_sons_con[5][1]=3;
        }
        break;
      case NORM_HEXA8:
        {
          _nb_of_pts=8; _nb_of_sons=6; _dim=3; _quadratic_type=NORM_HEXA20; _quadratic_type2=NORM_HEXA27; _is_simplex=false; _is_extruded=true; _reverse_extruded_type=NORM_QUAD4;
          _sons_type[0]=NORM_QUAD4; _sons_type[1]=NORM_QUAD4; _sons_type[2]=NORM_QUAD4; _sons_type[3]=NORM_QUAD4; _sons_type[4]=NORM_QUAD4; _sons_type[5]=NORM_QUAD4;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _sons_con[0][3]=3; _nb_of_sons_con[0]=4;
          _sons_con[1][0]=4; _sons_con[1][1]=7; _sons_con[1][2]=6; _sons_con[1][3]=5; _nb_of_sons_con[1]=4;
          _sons_con[2][0]=0; _sons_con[2][1]=4; _sons_con[2][2]=5; _sons_con[2][3]=1; _nb_of_sons_con[2]=4;
          _sons_con[3][0]=1; _sons_con[3][1]=5; _sons_con[3][2]=6; _sons_con[3][3]=2; _nb_of_sons_con[3]=4;
          _sons_con[4][0]=2; _sons_con[4][1]=6; _sons_con[4][2]=7; _sons_con[4][3]=3; _nb_of_sons_con[4]=4;
          _sons_con[5][0]=3; _sons_con[5][1]=7; _sons_con[5][2]=4; _sons_con[5][3]=0; _nb_of_sons_con[5]=4;
          _little_sons_con[0][0]=0; _little_sons_con[0][1]=1;  _nb_of_little_sons=12;
          _little_sons_con[1][0]=1; _little_sons_con[1][1]=2;
          _little_sons_con[2][0]=2; _little_sons_con[2][1]=3;
          _little_sons_con[3][0]=3; _little_sons_con[3][1]=0;
          _little_sons_con[4][0]=4; _little_sons_con[4][1]=5;
          _little_sons_con[5][0]=5; _little_sons_con[5][1]=6;
          _little_sons_con[6][0]=6; _little_sons_con[6][1]=7;
          _little_sons_con[7][0]=7; _little_sons_con[7][1]=4;
          _little_sons_con[8][0]=0; _little_sons_con[8][1]=4;
          _little_sons_con[9][0]=1; _little_sons_con[9][1]=5;
          _little_sons_con[10][0]=2; _little_sons_con[10][1]=6;
          _little_sons_con[11][0]=3; _little_sons_con[11][1]=7;
        }
        break;
      case NORM_QUAD4:
        {
          _nb_of_pts=4; _nb_of_sons=4; _dim=2; _quadratic_type=NORM_QUAD8; _quadratic_type2=NORM_QUAD9; _is_simplex=false; _is_extruded=true;
          _sons_type[0]=NORM_SEG2; _sons_type[1]=NORM_SEG2; _sons_type[2]=NORM_SEG2; _sons_type[3]=NORM_SEG2;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _nb_of_sons_con[0]=2;
          _sons_con[1][0]=1; _sons_con[1][1]=2; _nb_of_sons_con[1]=2;
          _sons_con[2][0]=2; _sons_con[2][1]=3; _nb_of_sons_con[2]=2;
          _sons_con[3][0]=3; _sons_con[3][1]=0; _nb_of_sons_con[3]=2; _extruded_type=NORM_HEXA8;
        }
        break;
      case NORM_TRI3:
        {
          _nb_of_pts=3; _nb_of_sons=3; _dim=2; _quadratic_type=NORM_TRI6; _quadratic_type2=NORM_TRI7; _is_simplex=true;
          _sons_type[0]=NORM_SEG2; _sons_type[1]=NORM_SEG2; _sons_type[2]=NORM_SEG2;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _nb_of_sons_con[0]=2;
          _sons_con[1][0]=1; _sons_con[1][1]=2; _nb_of_sons_con[1]=2;
          _sons_con[2][0]=2; _sons_con[2][1]=0; _nb_of_sons_con[2]=2; _extruded_type=NORM_PENTA6;
        }
        break;
      case NORM_TRI6:
        {
          _nb_of_pts=6; _nb_of_sons=3; _dim=2; _linear_type=NORM_TRI3; _is_simplex=false;
          _sons_type[0]=NORM_SEG3; _sons_type[1]=NORM_SEG3; _sons_type[2]=NORM_SEG3;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=3; _nb_of_sons_con[0]=3;
          _sons_con[1][0]=1; _sons_con[1][1]=2; _sons_con[1][2]=4; _nb_of_sons_con[1]=3;
          _sons_con[2][0]=2; _sons_con[2][1]=0; _sons_con[2][2]=5; _nb_of_sons_con[2]=3; _quadratic=true; _extruded_type=NORM_PENTA15;
        }
        break;
      case NORM_TRI7:
        {
          _nb_of_pts=7; _nb_of_sons=3; _dim=2; _linear_type=NORM_TRI3; _is_simplex=false;
          _sons_type[0]=NORM_SEG3; _sons_type[1]=NORM_SEG3; _sons_type[2]=NORM_SEG3;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=3; _nb_of_sons_con[0]=3;
          _sons_con[1][0]=1; _sons_con[1][1]=2; _sons_con[1][2]=4; _nb_of_sons_con[1]=3;
          _sons_con[2][0]=2; _sons_con[2][1]=0; _sons_con[2][2]=5; _nb_of_sons_con[2]=3; _quadratic=true; //no extruded type because no penta20
        }
        break;
      case NORM_QUAD8:
        {
          _nb_of_pts=8; _nb_of_sons=4; _dim=2; _linear_type=NORM_QUAD4; _is_simplex=false;
          _sons_type[0]=NORM_SEG3; _sons_type[1]=NORM_SEG3; _sons_type[2]=NORM_SEG3; _sons_type[3]=NORM_SEG3;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=4; _nb_of_sons_con[0]=3;
          _sons_con[1][0]=1; _sons_con[1][1]=2; _sons_con[1][2]=5; _nb_of_sons_con[1]=3;
          _sons_con[2][0]=2; _sons_con[2][1]=3; _sons_con[2][2]=6; _nb_of_sons_con[2]=3;
          _sons_con[3][0]=3; _sons_con[3][1]=0; _sons_con[3][2]=7; _nb_of_sons_con[3]=3; _quadratic=true; _extruded_type=NORM_HEXA20;
        }
        break;
      case NORM_QUAD9:
        {
          _nb_of_pts=9; _nb_of_sons=4; _dim=2; _linear_type=NORM_QUAD4; _is_simplex=false;
          _sons_type[0]=NORM_SEG3; _sons_type[1]=NORM_SEG3; _sons_type[2]=NORM_SEG3; _sons_type[3]=NORM_SEG3;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=4; _nb_of_sons_con[0]=3;
          _sons_con[1][0]=1; _sons_con[1][1]=2; _sons_con[1][2]=5; _nb_of_sons_con[1]=3;
          _sons_con[2][0]=2; _sons_con[2][1]=3; _sons_con[2][2]=6; _nb_of_sons_con[2]=3;
          _sons_con[3][0]=3; _sons_con[3][1]=0; _sons_con[3][2]=7; _nb_of_sons_con[3]=3; _quadratic=true; _extruded_type=NORM_HEXA27;
        }
        break;
      case NORM_PYRA5:
        {
          _nb_of_pts=5; _nb_of_sons=5; _dim=3; _quadratic_type=NORM_PYRA13; _is_simplex=false;
          _sons_type[0]=NORM_QUAD4; _sons_type[1]=NORM_TRI3; _sons_type[2]=NORM_TRI3; _sons_type[3]=NORM_TRI3; _sons_type[4]=NORM_TRI3;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _sons_con[0][3]=3; _nb_of_sons_con[0]=4;
          _sons_con[1][0]=0; _sons_con[1][1]=4; _sons_con[1][2]=1; _nb_of_sons_con[1]=3;
          _sons_con[2][0]=1; _sons_con[2][1]=4; _sons_con[2][2]=2; _nb_of_sons_con[2]=3;
          _sons_con[3][0]=2; _sons_con[3][1]=4; _sons_con[3][2]=3; _nb_of_sons_con[3]=3;
          _sons_con[4][0]=3; _sons_con[4][1]=4; _sons_con[4][2]=0; _nb_of_sons_con[4]=3;
          _little_sons_con[0][0]=0; _little_sons_con[0][1]=1;  _nb_of_little_sons=8;
          _little_sons_con[1][0]=1; _little_sons_con[1][1]=2;
          _little_sons_con[2][0]=2; _little_sons_con[2][1]=3;
          _little_sons_con[3][0]=3; _little_sons_con[3][1]=0;
          _little_sons_con[4][0]=0; _little_sons_con[4][1]=4;
          _little_sons_con[5][0]=1; _little_sons_con[5][1]=4;
          _little_sons_con[6][0]=2; _little_sons_con[6][1]=4;
          _little_sons_con[7][0]=3; _little_sons_con[7][1]=4;
        }
        break;
      case NORM_PENTA6:
        {
          _nb_of_pts=6; _nb_of_sons=5; _dim=3; _quadratic_type=NORM_PENTA15; _is_simplex=false; _is_extruded=true; _reverse_extruded_type=NORM_TRI3;
          _sons_type[0]=NORM_TRI3; _sons_type[1]=NORM_TRI3; _sons_type[2]=NORM_QUAD4; _sons_type[3]=NORM_QUAD4; _sons_type[4]=NORM_QUAD4;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _nb_of_sons_con[0]=3;
          _sons_con[1][0]=3; _sons_con[1][1]=5; _sons_con[1][2]=4; _nb_of_sons_con[1]=3;
          _sons_con[2][0]=0; _sons_con[2][1]=3; _sons_con[2][2]=4; _sons_con[2][3]=1; _nb_of_sons_con[2]=4;
          _sons_con[3][0]=1; _sons_con[3][1]=4; _sons_con[3][2]=5; _sons_con[3][3]=2; _nb_of_sons_con[3]=4;
          _sons_con[4][0]=2; _sons_con[4][1]=5; _sons_con[4][2]=3; _sons_con[4][3]=0; _nb_of_sons_con[4]=4;
          _little_sons_con[0][0]=0; _little_sons_con[0][1]=1;  _nb_of_little_sons=9;
          _little_sons_con[1][0]=1; _little_sons_con[1][1]=2;
          _little_sons_con[2][0]=2; _little_sons_con[2][1]=0;
          _little_sons_con[3][0]=3; _little_sons_con[3][1]=4;
          _little_sons_con[4][0]=4; _little_sons_con[4][1]=5;
          _little_sons_con[5][0]=5; _little_sons_con[5][1]=3;
          _little_sons_con[6][0]=0; _little_sons_con[6][1]=3;
          _little_sons_con[7][0]=1; _little_sons_con[7][1]=4;
          _little_sons_con[8][0]=2; _little_sons_con[8][1]=5;
        }
        break;
      case NORM_TETRA10:
        {
          _nb_of_pts=10; _nb_of_sons=4; _dim=3; _linear_type=NORM_TETRA4; _is_simplex=false;
          _sons_type[0]=NORM_TRI6; _sons_type[1]=NORM_TRI6; _sons_type[2]=NORM_TRI6; _sons_type[3]=NORM_TRI6;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _sons_con[0][3]=4; _sons_con[0][4]=5; _sons_con[0][5]=6; _nb_of_sons_con[0]=6;
          _sons_con[1][0]=0; _sons_con[1][1]=3; _sons_con[1][2]=1; _sons_con[1][3]=7; _sons_con[1][4]=8; _sons_con[1][5]=4; _nb_of_sons_con[1]=6;
          _sons_con[2][0]=1; _sons_con[2][1]=3; _sons_con[2][2]=2; _sons_con[2][3]=8; _sons_con[2][4]=9; _sons_con[2][5]=5; _nb_of_sons_con[2]=6;
          _sons_con[3][0]=2; _sons_con[3][1]=3; _sons_con[3][2]=0; _sons_con[3][3]=9; _sons_con[3][4]=7; _sons_con[3][5]=6; _nb_of_sons_con[3]=6;  _quadratic=true;
          _little_sons_con[0][0]=0; _little_sons_con[0][1]=1;  _little_sons_con[0][2]=4;  _nb_of_little_sons=6;
          _little_sons_con[1][0]=1; _little_sons_con[1][1]=2;  _little_sons_con[1][2]=5;
          _little_sons_con[2][0]=2; _little_sons_con[2][1]=0;  _little_sons_con[2][2]=6;
          _little_sons_con[3][0]=0; _little_sons_con[3][1]=3;  _little_sons_con[3][2]=7;
          _little_sons_con[4][0]=1; _little_sons_con[4][1]=3;  _little_sons_con[4][2]=8;
          _little_sons_con[5][0]=2; _little_sons_con[5][1]=3;  _little_sons_con[5][2]=9;
        }
        break;
      case NORM_HEXGP12:
        {
          _nb_of_pts=12; _nb_of_sons=8; _dim=3; _is_simplex=false; _is_extruded=true;
          _sons_type[0]=NORM_POLYGON; _sons_type[1]=NORM_POLYGON; _sons_type[2]=NORM_QUAD4; _sons_type[3]=NORM_QUAD4; _sons_type[4]=NORM_QUAD4; _sons_type[5]=NORM_QUAD4;
          _sons_type[6]=NORM_QUAD4; _sons_type[7]=NORM_QUAD4;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _sons_con[0][3]=3; _sons_con[0][4]=4; _sons_con[0][5]=5; _nb_of_sons_con[0]=6;
          _sons_con[1][0]=6; _sons_con[1][1]=11; _sons_con[1][2]=10; _sons_con[1][3]=9; _sons_con[1][4]=8; _sons_con[1][5]=7; _nb_of_sons_con[1]=6;
          _sons_con[2][0]=0; _sons_con[2][1]=6; _sons_con[2][2]=7; _sons_con[2][3]=1; _nb_of_sons_con[2]=4;
          _sons_con[3][0]=1; _sons_con[3][1]=7; _sons_con[3][2]=8; _sons_con[3][3]=2; _nb_of_sons_con[3]=4;
          _sons_con[4][0]=2; _sons_con[4][1]=8; _sons_con[4][2]=9; _sons_con[4][3]=3; _nb_of_sons_con[4]=4;
          _sons_con[5][0]=3; _sons_con[5][1]=9; _sons_con[5][2]=10; _sons_con[5][3]=4; _nb_of_sons_con[5]=4;
          _sons_con[6][0]=4; _sons_con[6][1]=10; _sons_con[6][2]=11; _sons_con[6][3]=5; _nb_of_sons_con[6]=4;
          _sons_con[7][0]=5; _sons_con[7][1]=11; _sons_con[7][2]=6; _sons_con[7][3]=0; _nb_of_sons_con[7]=4;
        }
        break;
      case NORM_PYRA13:
        {
          _nb_of_pts=13; _nb_of_sons=5; _dim=3; _linear_type=NORM_PYRA5; _is_simplex=false;
          _sons_type[0]=NORM_QUAD8; _sons_type[1]=NORM_TRI6; _sons_type[2]=NORM_TRI6; _sons_type[3]=NORM_TRI6; _sons_type[4]=NORM_TRI6;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _sons_con[0][3]=3; _sons_con[0][4]=5; _sons_con[0][5]=6; _sons_con[0][6]=7; _sons_con[0][7]=8; _nb_of_sons_con[0]=8;
          _sons_con[1][0]=0; _sons_con[1][1]=4; _sons_con[1][2]=1; _sons_con[1][3]=9; _sons_con[1][4]=10; _sons_con[1][5]=5; _nb_of_sons_con[1]=6;
          _sons_con[2][0]=1; _sons_con[2][1]=4; _sons_con[2][2]=2; _sons_con[2][3]=10; _sons_con[2][4]=11; _sons_con[2][5]=6; _nb_of_sons_con[2]=6;
          _sons_con[3][0]=2; _sons_con[3][1]=4; _sons_con[3][2]=3; _sons_con[3][3]=11; _sons_con[3][4]=12; _sons_con[3][5]=7;  _nb_of_sons_con[3]=6;
          _sons_con[4][0]=3; _sons_con[4][1]=4; _sons_con[4][2]=0; _sons_con[4][3]=12; _sons_con[4][4]=9; _sons_con[4][5]=8; _nb_of_sons_con[4]=6; _quadratic=true;
          _little_sons_con[0][0]=0; _little_sons_con[0][1]=1; _little_sons_con[0][2]=5;  _nb_of_little_sons=8;
          _little_sons_con[1][0]=1; _little_sons_con[1][1]=2; _little_sons_con[1][2]=6;
          _little_sons_con[2][0]=2; _little_sons_con[2][1]=3; _little_sons_con[2][2]=7;
          _little_sons_con[3][0]=3; _little_sons_con[3][1]=0; _little_sons_con[3][2]=8;
          _little_sons_con[4][0]=0; _little_sons_con[4][1]=4; _little_sons_con[4][2]=9;
          _little_sons_con[5][0]=1; _little_sons_con[5][1]=4; _little_sons_con[5][2]=10;
          _little_sons_con[6][0]=2; _little_sons_con[6][1]=4; _little_sons_con[6][2]=11;
          _little_sons_con[7][0]=3; _little_sons_con[7][1]=4; _little_sons_con[7][2]=12;
        }
        break;
      case NORM_PENTA15:
        {
          _nb_of_pts=15; _nb_of_sons=5; _dim=3; _linear_type=NORM_PENTA6; _is_simplex=false;
          _sons_type[0]=NORM_TRI6; _sons_type[1]=NORM_TRI6; _sons_type[2]=NORM_QUAD8; _sons_type[3]=NORM_QUAD8; _sons_type[4]=NORM_QUAD8;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _sons_con[0][3]=6; _sons_con[0][4]=7; _sons_con[0][5]=8; _nb_of_sons_con[0]=6;
          _sons_con[1][0]=3; _sons_con[1][1]=5; _sons_con[1][2]=4; _sons_con[1][3]=11; _sons_con[1][4]=10; _sons_con[1][5]=9; _nb_of_sons_con[1]=6;
          _sons_con[2][0]=0; _sons_con[2][1]=3; _sons_con[2][2]=4; _sons_con[2][3]=1; _sons_con[2][4]=12; _sons_con[2][5]=9; _sons_con[2][6]=13; _sons_con[2][7]=6; _nb_of_sons_con[2]=8;
          _sons_con[3][0]=1; _sons_con[3][1]=4; _sons_con[3][2]=5; _sons_con[3][3]=2; _sons_con[3][4]=13; _sons_con[3][5]=10; _sons_con[3][6]=14; _sons_con[3][7]=7; _nb_of_sons_con[3]=8;
          _sons_con[4][0]=2; _sons_con[4][1]=4; _sons_con[4][2]=5; _sons_con[4][3]=0; _sons_con[4][4]=14; _sons_con[4][5]=11; _sons_con[4][6]=12; _sons_con[4][7]=8; _nb_of_sons_con[4]=8; _quadratic=true;
          _little_sons_con[0][0]=0; _little_sons_con[0][1]=1; _little_sons_con[0][2]=6;  _nb_of_little_sons=9;
          _little_sons_con[1][0]=1; _little_sons_con[1][1]=2; _little_sons_con[1][2]=7;
          _little_sons_con[2][0]=2; _little_sons_con[2][1]=0; _little_sons_con[2][2]=8;
          _little_sons_con[3][0]=3; _little_sons_con[3][1]=4; _little_sons_con[3][2]=9;
          _little_sons_con[4][0]=4; _little_sons_con[4][1]=5; _little_sons_con[4][2]=10;
          _little_sons_con[5][0]=5; _little_sons_con[5][1]=3; _little_sons_con[5][2]=11;
          _little_sons_con[6][0]=0; _little_sons_con[6][1]=3; _little_sons_con[6][2]=12;
          _little_sons_con[7][0]=1; _little_sons_con[7][1]=4; _little_sons_con[7][2]=13;
          _little_sons_con[8][0]=2; _little_sons_con[8][1]=5; _little_sons_con[8][2]=14;
        }
        break;
      case NORM_PENTA18:
        {
          _nb_of_pts=18; _nb_of_sons=5; _dim=3; _linear_type=NORM_PENTA6; _is_simplex=false;
          _sons_type[0]=NORM_TRI6; _sons_type[1]=NORM_TRI6; _sons_type[2]=NORM_QUAD9; _sons_type[3]=NORM_QUAD9; _sons_type[4]=NORM_QUAD9;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _sons_con[0][3]=6; _sons_con[0][4]=7; _sons_con[0][5]=8; _nb_of_sons_con[0]=6;
          _sons_con[1][0]=3; _sons_con[1][1]=5; _sons_con[1][2]=4; _sons_con[1][3]=11; _sons_con[1][4]=10; _sons_con[1][5]=9; _nb_of_sons_con[1]=6;
          _sons_con[2][0]=0; _sons_con[2][1]=3; _sons_con[2][2]=4; _sons_con[2][3]=1; _sons_con[2][4]=12; _sons_con[2][5]=9; _sons_con[2][6]=13; _sons_con[2][7]=6; _sons_con[2][8]=15; _nb_of_sons_con[2]=9;
          _sons_con[3][0]=1; _sons_con[3][1]=4; _sons_con[3][2]=5; _sons_con[3][3]=2; _sons_con[3][4]=13; _sons_con[3][5]=10; _sons_con[3][6]=14; _sons_con[3][7]=7; _sons_con[3][8]=16; _nb_of_sons_con[3]=9;
          _sons_con[4][0]=2; _sons_con[4][1]=4; _sons_con[4][2]=5; _sons_con[4][3]=0; _sons_con[4][4]=14; _sons_con[4][5]=11; _sons_con[4][6]=12; _sons_con[4][7]=8; _sons_con[4][8]=17; _nb_of_sons_con[4]=9; _quadratic=true;
          _little_sons_con[0][0]=0; _little_sons_con[0][1]=1; _little_sons_con[0][2]=6;  _nb_of_little_sons=9;
          _little_sons_con[1][0]=1; _little_sons_con[1][1]=2; _little_sons_con[1][2]=7;
          _little_sons_con[2][0]=2; _little_sons_con[2][1]=0; _little_sons_con[2][2]=8;
          _little_sons_con[3][0]=3; _little_sons_con[3][1]=4; _little_sons_con[3][2]=9;
          _little_sons_con[4][0]=4; _little_sons_con[4][1]=5; _little_sons_con[4][2]=10;
          _little_sons_con[5][0]=5; _little_sons_con[5][1]=3; _little_sons_con[5][2]=11;
          _little_sons_con[6][0]=0; _little_sons_con[6][1]=3; _little_sons_con[6][2]=12;
          _little_sons_con[7][0]=1; _little_sons_con[7][1]=4; _little_sons_con[7][2]=13;
          _little_sons_con[8][0]=2; _little_sons_con[8][1]=5; _little_sons_con[8][2]=14;
        }
        break;
      case NORM_HEXA20:
        {
          _nb_of_pts=20; _nb_of_sons=6; _dim=3; _linear_type=NORM_HEXA8; _is_simplex=false;
          _sons_type[0]=NORM_QUAD8; _sons_type[1]=NORM_QUAD8; _sons_type[2]=NORM_QUAD8; _sons_type[3]=NORM_QUAD8; _sons_type[4]=NORM_QUAD8; _sons_type[5]=NORM_QUAD8;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _sons_con[0][3]=3; _sons_con[0][4]=8; _sons_con[0][5]=9; _sons_con[0][6]=10; _sons_con[0][7]=11; _nb_of_sons_con[0]=8;
          _sons_con[1][0]=4; _sons_con[1][1]=7; _sons_con[1][2]=6; _sons_con[1][3]=5; _sons_con[1][4]=15; _sons_con[1][5]=14; _sons_con[1][6]=13; _sons_con[1][7]=12; _nb_of_sons_con[1]=8;
          _sons_con[2][0]=0; _sons_con[2][1]=4; _sons_con[2][2]=5; _sons_con[2][3]=1; _sons_con[2][4]=16; _sons_con[2][5]=12; _sons_con[2][6]=17; _sons_con[2][7]=8; _nb_of_sons_con[2]=8;
          _sons_con[3][0]=1; _sons_con[3][1]=5; _sons_con[3][2]=6; _sons_con[3][3]=2; _sons_con[3][4]=17; _sons_con[3][5]=13; _sons_con[3][6]=18; _sons_con[3][7]=9;_nb_of_sons_con[3]=8;
          _sons_con[4][0]=2; _sons_con[4][1]=6; _sons_con[4][2]=7; _sons_con[4][3]=3; _sons_con[4][4]=18; _sons_con[4][5]=14; _sons_con[4][6]=19; _sons_con[4][7]=10; _nb_of_sons_con[4]=8;
          _sons_con[5][0]=3; _sons_con[5][1]=7; _sons_con[5][2]=4; _sons_con[5][3]=0; _sons_con[5][4]=19; _sons_con[5][5]=15; _sons_con[5][6]=16; _sons_con[5][7]=11; _nb_of_sons_con[5]=8; _quadratic=true;
          _little_sons_con[0][0]=0; _little_sons_con[0][1]=1;  _little_sons_con[0][2]=8; _nb_of_little_sons=12;
          _little_sons_con[1][0]=1; _little_sons_con[1][1]=2;  _little_sons_con[1][2]=9;
          _little_sons_con[2][0]=2; _little_sons_con[2][1]=3;  _little_sons_con[2][2]=10;
          _little_sons_con[3][0]=3; _little_sons_con[3][1]=0;  _little_sons_con[3][2]=11;
          _little_sons_con[4][0]=4; _little_sons_con[4][1]=5;  _little_sons_con[4][2]=12;
          _little_sons_con[5][0]=5; _little_sons_con[5][1]=6;  _little_sons_con[5][2]=13;
          _little_sons_con[6][0]=6; _little_sons_con[6][1]=7;  _little_sons_con[6][2]=14;
          _little_sons_con[7][0]=7; _little_sons_con[7][1]=4;  _little_sons_con[7][2]=15;
          _little_sons_con[8][0]=0; _little_sons_con[8][1]=4;  _little_sons_con[8][2]=16;
          _little_sons_con[9][0]=1; _little_sons_con[9][1]=5;  _little_sons_con[9][2]=17;
          _little_sons_con[10][0]=2; _little_sons_con[10][1]=6;  _little_sons_con[10][2]=18;
          _little_sons_con[11][0]=3; _little_sons_con[11][1]=7;  _little_sons_con[11][2]=19;
        }
        break;
      case NORM_HEXA27:
        {
          _nb_of_pts=27; _nb_of_sons=6; _dim=3; _linear_type=NORM_HEXA8; _is_simplex=false;
          _sons_type[0]=NORM_QUAD9; _sons_type[1]=NORM_QUAD9; _sons_type[2]=NORM_QUAD9; _sons_type[3]=NORM_QUAD9; _sons_type[4]=NORM_QUAD9; _sons_type[5]=NORM_QUAD9;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _sons_con[0][3]=3; _sons_con[0][4]=8; _sons_con[0][5]=9; _sons_con[0][6]=10; _sons_con[0][7]=11; _sons_con[0][8]=20; _nb_of_sons_con[0]=9;
          _sons_con[1][0]=4; _sons_con[1][1]=7; _sons_con[1][2]=6; _sons_con[1][3]=5; _sons_con[1][4]=15; _sons_con[1][5]=14; _sons_con[1][6]=13; _sons_con[1][7]=12; _sons_con[1][8]=25; _nb_of_sons_con[1]=9;
          _sons_con[2][0]=0; _sons_con[2][1]=4; _sons_con[2][2]=5; _sons_con[2][3]=1; _sons_con[2][4]=16; _sons_con[2][5]=12; _sons_con[2][6]=17; _sons_con[2][7]=8; _sons_con[2][8]=21; _nb_of_sons_con[2]=9;   
          _sons_con[3][0]=1; _sons_con[3][1]=5; _sons_con[3][2]=6; _sons_con[3][3]=2; _sons_con[3][4]=17; _sons_con[3][5]=13; _sons_con[3][6]=18; _sons_con[3][7]=9; _sons_con[3][8]=22; _nb_of_sons_con[3]=9;
          _sons_con[4][0]=2; _sons_con[4][1]=6; _sons_con[4][2]=7; _sons_con[4][3]=3; _sons_con[4][4]=18; _sons_con[4][5]=14; _sons_con[4][6]=19; _sons_con[4][7]=10; _sons_con[4][8]=23; _nb_of_sons_con[4]=9;
          _sons_con[5][0]=3; _sons_con[5][1]=7; _sons_con[5][2]=4; _sons_con[5][3]=0; _sons_con[5][4]=19; _sons_con[5][5]=15; _sons_con[5][6]=16; _sons_con[5][7]=11; _sons_con[5][8]=24; _nb_of_sons_con[5]=9;
          _quadratic=true;
        }
        break;
      case NORM_POLYGON:
        {
          _nb_of_pts=0; _nb_of_sons=0; _dim=2; _dyn=true; _extruded_type=NORM_POLYHED; _is_simplex=false; _quadratic_type=NORM_QPOLYG;
        }
        break;
      case NORM_POLYHED:
        {
          _nb_of_pts=0; _nb_of_sons=0; _dim=3; _dyn=true; _is_simplex=false;
        }
        break;
      case NORM_QPOLYG:
        {
          _nb_of_pts=0; _nb_of_sons=0; _dim=2; _dyn=true; _is_simplex=false; _quadratic=true; _linear_type=NORM_POLYGON;
        }
        break;
      case NORM_POLYL:
        {
          _nb_of_pts=0; _nb_of_sons=0; _dim=1; _dyn=true; _extruded_type=NORM_POLYGON; _is_simplex=false;
        }
        break;
      case NORM_ERROR:
        {
          _nb_of_pts=std::numeric_limits<unsigned>::max(); _nb_of_sons=std::numeric_limits<unsigned>::max(); _dim=std::numeric_limits<unsigned>::max();
        }
        break;
      }
  }

  /*!
   * Equivalent to getNumberOfSons except that this method deals with dynamic type.
   */
  unsigned CellModel::getNumberOfSons2(const int *conn, int lgth) const
  {
    if(!isDynamic())
      return getNumberOfSons();
    if(_dim==2)
      {
        if(_type==NORM_POLYGON)
          return lgth;
        else
          return lgth/2;
      }
    else if(_dim==1)
      return lgth;//NORM_POLYL
    else
      return std::count(conn,conn+lgth,-1)+1;
  }

  unsigned CellModel::getNumberOfEdgesIn3D(const int *conn, int lgth) const
  {
    if(!isDynamic())
      return _nb_of_little_sons;
    else//polyhedron
      return (lgth-std::count(conn,conn+lgth,-1))/2;
  }
  
  /*!
   * \sa fillMicroEdgeNodalConnectivity
   */
  unsigned CellModel::getNumberOfMicroEdges() const
  {
    unsigned mul(isQuadratic()?2:1);
    if(!isDynamic())
      {
        switch(getDimension())
          {
          case 2:
            return mul*getNumberOfSons();
          case 3:
            return mul*_nb_of_little_sons;
          default:
            throw INTERP_KERNEL::Exception("CellModel::getNumberOfMacroEdges : only 2D and 3D cells support this !");
          }
      }
    else
      throw INTERP_KERNEL::Exception("CellModel::getNumberOfMacroEdges : not supported by dynamic type !");
  }
  
  NormalizedCellType CellModel::getCorrespondingPolyType() const
  {
    switch(getDimension())
      {
      case 0:
        return NORM_POINT1;
      case 1:
        {
          if(!isQuadratic())
            return NORM_POLYL;
          throw INTERP_KERNEL::Exception("CellModel::getPolyType : no poly type for quadratic 1D !");
        }
      case 2:
        {
          if(!isQuadratic())
            return NORM_POLYGON;
          else
            return NORM_QPOLYG;
        }
      case 3:
        {
          if(!isQuadratic())
            return NORM_POLYHED;
          throw INTERP_KERNEL::Exception("CellModel::getPolyType : no poly type for quadratic 3D !");
        }
      default:
        throw INTERP_KERNEL::Exception("CellModel::getPolyType : only dimension 0, 1, 2, 3 are supported !");
      }
  }

  /*!
   * Equivalent to getSonType except that this method deals with dynamic type.
   */
  NormalizedCellType CellModel::getSonType2(unsigned sonId) const
  {
    if(!isDynamic())
      return getSonType(sonId);
    if(_dim==2)
      {
        if(_type==NORM_POLYGON)
          return NORM_SEG2;
        else
          return NORM_SEG3;
      }
    else if(_dim==1)
      return NORM_ERROR;//NORM_POLYL
    //polyedron
    return NORM_POLYGON;
  }

  /*!
   * \b WARNING this method do not manage correctly types that return true at the call of isDynamic. Use fillSonCellNodalConnectivity2 instead.
   */
  unsigned CellModel::fillSonCellNodalConnectivity(int sonId, const int *nodalConn, int *sonNodalConn) const
  {
    unsigned nbOfTurnLoop=_nb_of_sons_con[sonId];
    const unsigned *sonConn=_sons_con[sonId];
    for(unsigned i=0;i<nbOfTurnLoop;i++)
      sonNodalConn[i]=nodalConn[sonConn[i]];
    return nbOfTurnLoop;
  }

  unsigned CellModel::fillSonCellNodalConnectivity2(int sonId, const int *nodalConn, int lgth, int *sonNodalConn, NormalizedCellType& typeOfSon) const
  {
    typeOfSon=getSonType2(sonId);
    if(!isDynamic())
      return fillSonCellNodalConnectivity(sonId,nodalConn,sonNodalConn);
    else
      {
        if(_dim==2)//polygon
          {
            if(_type==NORM_POLYGON)
              {
                sonNodalConn[0]=nodalConn[sonId];
                sonNodalConn[1]=nodalConn[(sonId+1)%lgth];
                return 2;
              }
            else
              {
                sonNodalConn[0]=nodalConn[sonId];
                sonNodalConn[1]=nodalConn[(sonId+1)%(lgth/2)];
                sonNodalConn[2]=nodalConn[sonId+(lgth/2)];
                return 3;
              }
          }
        else if(_dim==3)
          {//polyedron
            const int *where=nodalConn;
            for(int i=0;i<sonId;i++)
              {
                where=std::find(where,nodalConn+lgth,-1);
                where++;
              }
            const int *where2=std::find(where,nodalConn+lgth,-1);
            std::copy(where,where2,sonNodalConn);
            return where2-where;
          }
        else
          throw INTERP_KERNEL::Exception("CellModel::fillSonCellNodalConnectivity2 : no sons on NORM_POLYL !");
      }
  }
  
  /*!
   * Equivalent to CellModel::fillSonCellNodalConnectivity2 except for HEXA8 where the order of sub faces is not has MED file numbering for transformation HEXA8->HEXA27
   */
  unsigned CellModel::fillSonCellNodalConnectivity4(int sonId, const int *nodalConn, int lgth, int *sonNodalConn, NormalizedCellType& typeOfSon) const
  {
    if(_type==NORM_HEXA8)
      {
        static const int permutation[6]={0,2,3,4,5,1};
        return fillSonCellNodalConnectivity2(permutation[sonId],nodalConn,lgth,sonNodalConn,typeOfSon);
      }
    else
      return fillSonCellNodalConnectivity2(sonId,nodalConn,lgth,sonNodalConn,typeOfSon);
  }

  unsigned CellModel::fillSonEdgesNodalConnectivity3D(int sonId, const int *nodalConn, int lgth, int *sonNodalConn, NormalizedCellType& typeOfSon) const
  {
    if(!isDynamic())
      {
        if(!isQuadratic())
          {
            typeOfSon=NORM_SEG2;
            sonNodalConn[0]=nodalConn[_little_sons_con[sonId][0]];
            sonNodalConn[1]=nodalConn[_little_sons_con[sonId][1]];
            return 2;
          }
        else
          {
            typeOfSon=NORM_SEG3;
            sonNodalConn[0]=nodalConn[_little_sons_con[sonId][0]];
            sonNodalConn[1]=nodalConn[_little_sons_con[sonId][1]];
            sonNodalConn[2]=nodalConn[_little_sons_con[sonId][2]];
            return 3;
          }
      }
    else
      throw INTERP_KERNEL::Exception("CellModel::fillSonEdgesNodalConnectivity3D : not implemented yet for NORM_POLYHED !");   
  }

  /*!
   * \sa getNumberOfMicroEdges
   */
  unsigned CellModel::fillMicroEdgeNodalConnectivity(int sonId, const int *nodalConn, int *sonNodalConn, NormalizedCellType& typeOfSon) const
  {
    if(isQuadratic())
      {
        int edgeId(sonId/2),subEdgeId(sonId%2);
        typeOfSon=NORM_SEG2;
        const unsigned *sonConn(0);
        switch(getDimension())
          {
          case 2:
            {
              sonConn=_sons_con[edgeId];
              break;
            }
          case 3:
            {
              sonConn=_little_sons_con[edgeId];
              break;
            }
          default:
            throw INTERP_KERNEL::Exception("CellModel::fillMicroEdgeNodalConnectivity : only 2D and 3D cells support this !");
          }
        const unsigned tmp[3]={sonConn[0],sonConn[2],sonConn[1]};
        sonNodalConn[0]=nodalConn[tmp[subEdgeId]];
        sonNodalConn[1]=nodalConn[tmp[subEdgeId+1]];
        return 2;
      }
    else
      {
        switch(getDimension())
          {
          case 2:
            return fillSonCellNodalConnectivity2(sonId,nodalConn,0,sonNodalConn,typeOfSon);
          case 3:
            return fillSonEdgesNodalConnectivity3D(sonId,nodalConn,0,sonNodalConn,typeOfSon);
          default:
            throw INTERP_KERNEL::Exception("CellModel::fillMicroEdgeNodalConnectivity : only 2D and 3D cells support this #2 !");
          }
      }
  }

  void CellModel::changeOrientationOf2D(int *nodalConn, unsigned int sz) const
  {
    if(sz<1)
      return ;
    if(!isQuadratic())
      {
        std::vector<int> tmp(sz-1);
        std::copy(nodalConn+1,nodalConn+sz,tmp.rbegin());
        std::copy(tmp.begin(),tmp.end(),nodalConn+1);
      }
    else
      {
        unsigned int sz2(sz/2);
        std::vector<int> tmp0(sz2-1),tmp1(sz2);
        std::copy(nodalConn+1,nodalConn+sz2,tmp0.rbegin());
        std::copy(nodalConn+sz2,nodalConn+sz,tmp1.rbegin());
        std::copy(tmp0.begin(),tmp0.end(),nodalConn+1);
        std::copy(tmp1.begin(),tmp1.end(),nodalConn+sz2);
      }
  }

  void CellModel::changeOrientationOf1D(int *nodalConn, unsigned int sz) const
  {
    if(!isDynamic())
      {
        if(sz==2 || sz==3)
          {
            std::swap(nodalConn[0],nodalConn[1]);
            return ;
          }
        else if(sz==4)
          {
            std::swap(nodalConn[0],nodalConn[1]);
            std::swap(nodalConn[2],nodalConn[3]);
          }
        else
          throw Exception("CellModel::changeOrientationOf1D : unrecognized 1D cell type !");
      }
    else
      {
        std::vector<int> tmp(sz-1);
        std::copy(nodalConn+1,nodalConn+sz,tmp.rbegin());
        std::copy(tmp.begin(),tmp.end(),nodalConn+1);
      }
  }

  //================================================================================
  /*!
   * \brief Return number of nodes in sonId-th son of a Dynamic() cell
   */
  //================================================================================

  unsigned CellModel::getNumberOfNodesConstituentTheSon2(unsigned sonId, const int *nodalConn, int lgth) const
  {
    if(!isDynamic())
      return getNumberOfNodesConstituentTheSon(sonId);

    if(_dim==2)//polygon
      {
        if(_type==NORM_POLYGON)
          return 2;
        else
          return 3;
      }
    else if(_dim==3)
      {//polyedron
        const int *where=nodalConn;
        for(unsigned int i=0;i<sonId;i++)
          {
            where=std::find(where,nodalConn+lgth,-1);
            where++;
          }
        const int *where2=std::find(where,nodalConn+lgth,-1);
        return where2-where;
      }
    else
      throw INTERP_KERNEL::Exception("CellModel::getNumberOfNodesConstituentTheSon2 : no sons on NORM_POLYL !");
  }

  /*!
   * This method retrieves if cell1 represented by 'conn1' and cell2 represented by 'conn2'
   * are equivalent by a permutation or not. This method expects to work on 1D or 2D (only mesh dimension where it is possible to have a spaceDim) strictly higher than meshDim.
   * If not an exception will be thrown.
   * @return True if two cells have same orientation, false if not.
   */
  bool CellModel::getOrientationStatus(unsigned lgth, const int *conn1, const int *conn2) const
  {
    if(_dim!=1 && _dim!=2)
      throw INTERP_KERNEL::Exception("CellModel::getOrientationStatus : invalid dimension ! Must be 1 or 2 !");
    if(!_quadratic)
      {
        std::vector<int> tmp(2*lgth);
        std::vector<int>::iterator it=std::copy(conn1,conn1+lgth,tmp.begin());
        std::copy(conn1,conn1+lgth,it);
        it=std::search(tmp.begin(),tmp.end(),conn2,conn2+lgth);
        if(it==tmp.begin())
          return true;
        if(it!=tmp.end())
          return _dim!=1;
        std::vector<int>::reverse_iterator it2=std::search(tmp.rbegin(),tmp.rend(),conn2,conn2+lgth);
        if(it2!=tmp.rend())
          return false;
        throw INTERP_KERNEL::Exception("CellModel::getOrientationStatus : Request of orientation status of non equal connectively cells !");
      }
    else
      {
        if(_dim!=1)
          {
            std::vector<int> tmp(lgth);
            std::vector<int>::iterator it=std::copy(conn1,conn1+lgth/2,tmp.begin());
            std::copy(conn1,conn1+lgth/2,it);
            it=std::search(tmp.begin(),tmp.end(),conn2,conn2+lgth/2);
            int d=std::distance(tmp.begin(),it);
            if(it==tmp.end())
              return false;
            it=std::copy(conn1+lgth/2,conn1+lgth,tmp.begin());
            std::copy(conn1+lgth/2,conn1+lgth,it);
            it=std::search(tmp.begin(),tmp.end(),conn2,conn2+lgth);
            if(it==tmp.end())
              return false;
            int d2=std::distance(tmp.begin(),it);
            return d==d2;
          }
        else
          {
            int p=(lgth+1)/2;
            std::vector<int> tmp(2*p);
            std::vector<int>::iterator it=std::copy(conn1,conn1+p,tmp.begin());
            std::copy(conn1,conn1+p,it);
            it=std::search(tmp.begin(),tmp.end(),conn2,conn2+p);
            int d=std::distance(tmp.begin(),it);
            if(it==tmp.end())
              return false;
            tmp.resize(2*p-2);
            it=std::copy(conn1+p,conn1+lgth,tmp.begin());
            std::copy(conn1+p,conn1+lgth,it);
            it=std::search(tmp.begin(),tmp.end(),conn2+p,conn2+lgth);
            if(it==tmp.end())
              return false;
            int d2=std::distance(tmp.begin(),it);
            return d==d2;
          }
      }
  }
  
  DiameterCalculator *CellModel::buildInstanceOfDiameterCalulator(int spaceDim) const
  {
    switch(_type)
      {
      case NORM_TRI3:
        {
          switch(spaceDim)
            {
            case 2:
              return new DiameterCalulatorTRI3S2;
            case 3:
              return new DiameterCalulatorTRI3S3;
            default:
              throw Exception("CellModel::buildInstanceOfDiameterCalulator : For TRI3 only space dimension 2 and 3 implemented !");
            }
          break;
        }
      case NORM_QUAD4:
        {
          switch(spaceDim)
            {
            case 2:
              return new DiameterCalulatorQUAD4S2;
            case 3:
              return new DiameterCalulatorQUAD4S3;
            default:
              throw Exception("CellModel::buildInstanceOfDiameterCalulator : For QUAD4 only space dimension 2 and 3 implemented !");
            }
          break;
        }
      case NORM_TRI6:
        {
          switch(spaceDim)
          {
            case 2:
              return new DiameterCalulatorTRI6S2;
            case 3:
              return new DiameterCalulatorTRI6S3;
            default:
              throw Exception("CellModel::buildInstanceOfDiameterCalulator : For TRI6 only space dimension 2 and 3 implemented !");
          }
          break;
        }
      case NORM_TRI7:
        {
          switch(spaceDim)
          {
            case 2:
              return new DiameterCalulatorTRI7S2;
            case 3:
              return new DiameterCalulatorTRI7S3;
            default:
              throw Exception("CellModel::buildInstanceOfDiameterCalulator : For TRI7 only space dimension 2 and 3 implemented !");
          }
          break;
        }
      case NORM_QUAD8:
        {
          switch(spaceDim)
          {
            case 2:
              return new DiameterCalulatorQUAD8S2;
            case 3:
              return new DiameterCalulatorQUAD8S3;
            default:
              throw Exception("CellModel::buildInstanceOfDiameterCalulator : For QUAD8 only space dimension 2 and 3 implemented !");
          }
          break;
        }
      case NORM_QUAD9:
        {
          switch(spaceDim)
          {
            case 2:
              return new DiameterCalulatorQUAD9S2;
            case 3:
              return new DiameterCalulatorQUAD9S3;
            default:
              throw Exception("CellModel::buildInstanceOfDiameterCalulator : For QUAD9 only space dimension 2 and 3 implemented !");
          }
          break;
        }
      case NORM_TETRA4:
        {
          if(spaceDim==3)
            return new DiameterCalulatorTETRA4;
          else
            throw Exception("CellModel::buildInstanceOfDiameterCalulator : For TETRA4 space dimension 3 expected !");
        }
      case NORM_TETRA10:
        {
          if(spaceDim==3)
            return new DiameterCalulatorTETRA10;
          else
            throw Exception("CellModel::buildInstanceOfDiameterCalulator : For TETRA10 space dimension 3 expected !");
        }
      case NORM_HEXA8:
        {
          if(spaceDim==3)
            return new DiameterCalulatorHEXA8;
          else
            throw Exception("CellModel::buildInstanceOfDiameterCalulator : For HEXA8 space dimension 3 expected !");
        }
      case NORM_HEXA20:
        {
          if(spaceDim==3)
            return new DiameterCalulatorHEXA20;
          else
            throw Exception("CellModel::buildInstanceOfDiameterCalulator : For HEXA20 space dimension 3 expected !");
        }
      case NORM_HEXA27:
        {
          if(spaceDim==3)
            return new DiameterCalulatorHEXA27;
          else
            throw Exception("CellModel::buildInstanceOfDiameterCalulator : For HEXA27 space dimension 3 expected !");
        }
      case NORM_PENTA6:
        {
          if(spaceDim==3)
            return new DiameterCalulatorPENTA6;
          else
            throw Exception("CellModel::buildInstanceOfDiameterCalulator : For PENTA6 space dimension 3 expected !");
        }
      case NORM_PENTA15:
        {
          if(spaceDim==3)
            return new DiameterCalulatorPENTA15;
          else
            throw Exception("CellModel::buildInstanceOfDiameterCalulator : For PENTA15 space dimension 3 expected !");
        }
      case NORM_PYRA5:
        {
          if(spaceDim==3)
            return new DiameterCalulatorPYRA5;
          else
            throw Exception("CellModel::buildInstanceOfDiameterCalulator : For PYRA5 space dimension 3 expected !");
        }
      case NORM_PYRA13:
        {
          if(spaceDim==3)
            return new DiameterCalulatorPYRA13;
          else
            throw Exception("CellModel::buildInstanceOfDiameterCalulator : For PYRA13 space dimension 3 expected !");
        }
      default:
        throw Exception("CellModel::buildInstanceOfDiameterCalulator : implemented only for TRI3, QUAD4, TETRA4, HEXA8, PENTA6, PYRA5 !");
      }
  }

  OrientationInverter *CellModel::buildOrientationInverter() const
  {
    switch(_type)
      {
      case NORM_SEG2:
        return new OrientationInverterSEG2;
      case NORM_SEG3:
        return new OrientationInverterSEG3;
      case NORM_TRI3:
      case NORM_QUAD4:
        return new OrientationInverter2DLinear(getNumberOfNodes());
      case NORM_TRI6:
      case NORM_QUAD8:
        return new OrientationInverter2DQuadratic(getNumberOfNodes());
      case NORM_POLYGON:
        return new OrientationInverterPolygon;
      case NORM_QPOLYG:
        return new OrientationInverterQPolygon;
      case NORM_TETRA4:
        return new OrientationInverterTetra4;
      case NORM_PYRA5:
        return new OrientationInverterPyra5;
      case NORM_TETRA10:
        return new OrientationInverterTetra10;
      case NORM_PYRA13:
        return new OrientationInverterPyra13;
      case NORM_PENTA6:
      case NORM_HEXA8:
        return new OrientationInverter3DExtrusionLinear(getNumberOfNodes());
      case NORM_PENTA15:
      case NORM_HEXA20:
        return new OrientationInverter3DExtrusionQuadratic(getNumberOfNodes());
      default:
        {
          std::ostringstream oss; oss << "CellModel::buildOrientationInverter : not managed geometric type " << getRepr() << " yet !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      }
  }
}
