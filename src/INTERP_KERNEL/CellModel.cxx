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
#include "CellModel.hxx"

#include "InterpKernelException.hxx"

#include <sstream>
#include <limits>

using namespace std;

namespace INTERP_KERNEL
{
  std::map<NormalizedCellType,CellModel> CellModel::_map_of_unique_instance;

  const CellModel& CellModel::getCellModel(NormalizedCellType type)
  {
    if(_map_of_unique_instance.empty())
      buildUniqueInstance();
    const map<NormalizedCellType,CellModel>::iterator iter=_map_of_unique_instance.find(type);
    if(iter==_map_of_unique_instance.end())
      {
        ostringstream stream; stream << "no cellmodel for normalized type " << type;
        throw Exception(stream.str().c_str());
      }
    return (*iter).second;
  }

  void CellModel::buildUniqueInstance()
  {
    _map_of_unique_instance.insert(make_pair(NORM_SEG2,CellModel(NORM_SEG2)));
    _map_of_unique_instance.insert(make_pair(NORM_SEG3,CellModel(NORM_SEG3)));
    _map_of_unique_instance.insert(make_pair(NORM_TRI3,CellModel(NORM_TRI3)));
    _map_of_unique_instance.insert(make_pair(NORM_QUAD4,CellModel(NORM_QUAD4)));
    _map_of_unique_instance.insert(make_pair(NORM_TRI6,CellModel(NORM_TRI6)));
    _map_of_unique_instance.insert(make_pair(NORM_QUAD8,CellModel(NORM_QUAD8)));
    _map_of_unique_instance.insert(make_pair(NORM_TETRA4,CellModel(NORM_TETRA4)));
    _map_of_unique_instance.insert(make_pair(NORM_HEXA8,CellModel(NORM_HEXA8)));
    _map_of_unique_instance.insert(make_pair(NORM_PYRA5,CellModel(NORM_PYRA5)));
    _map_of_unique_instance.insert(make_pair(NORM_PENTA6,CellModel(NORM_PENTA6)));
    _map_of_unique_instance.insert(make_pair(NORM_TETRA10,CellModel(NORM_TETRA10)));
    _map_of_unique_instance.insert(make_pair(NORM_PYRA13,CellModel(NORM_PYRA13)));
    _map_of_unique_instance.insert(make_pair(NORM_PENTA15,CellModel(NORM_PENTA15)));
    _map_of_unique_instance.insert(make_pair(NORM_HEXA20,CellModel(NORM_HEXA20)));
    _map_of_unique_instance.insert(make_pair(NORM_POLYGON,CellModel(NORM_POLYGON)));
    _map_of_unique_instance.insert(make_pair(NORM_POLYHED,CellModel(NORM_POLYHED)));
  }

  CellModel::CellModel(NormalizedCellType type)
  {
    _quadratic=false;
    _dyn=false;
    switch(type)
      {
      case NORM_SEG2:
        {
          _nb_of_pts=2; _nb_of_sons=0; _dim=1;
        }
        break;
      case NORM_SEG3:
        {
          _nb_of_pts=3; _nb_of_sons=0; _dim=1;
        }
        break;
      case NORM_TETRA4:
        {
          _nb_of_pts=4; _nb_of_sons=4; _dim=3;
          _sons_type[0]=NORM_TRI3; _sons_type[1]=NORM_TRI3; _sons_type[2]=NORM_TRI3; _sons_type[3]=NORM_TRI3;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _nb_of_sons_con[0]=3;
          _sons_con[1][0]=0; _sons_con[1][1]=3; _sons_con[1][2]=1; _nb_of_sons_con[1]=3;
          _sons_con[2][0]=1; _sons_con[2][1]=3; _sons_con[2][2]=2; _nb_of_sons_con[2]=3;
          _sons_con[3][0]=2; _sons_con[3][1]=3; _sons_con[3][2]=0; _nb_of_sons_con[3]=3;
        }
        break;
      case NORM_HEXA8:
        {
          _nb_of_pts=8; _nb_of_sons=6; _dim=3;
          _sons_type[0]=NORM_QUAD4; _sons_type[1]=NORM_QUAD4; _sons_type[2]=NORM_QUAD4; _sons_type[3]=NORM_QUAD4; _sons_type[4]=NORM_QUAD4; _sons_type[5]=NORM_QUAD4;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _sons_con[0][3]=3; _nb_of_sons_con[0]=4;
          _sons_con[1][0]=4; _sons_con[1][1]=7; _sons_con[1][2]=6; _sons_con[1][3]=5; _nb_of_sons_con[1]=4;
          _sons_con[2][0]=0; _sons_con[2][1]=4; _sons_con[2][2]=5; _sons_con[2][3]=1; _nb_of_sons_con[2]=4;
          _sons_con[3][0]=1; _sons_con[3][1]=5; _sons_con[3][2]=6; _sons_con[3][3]=2; _nb_of_sons_con[3]=4;
          _sons_con[4][0]=2; _sons_con[4][1]=6; _sons_con[4][2]=7; _sons_con[4][3]=3; _nb_of_sons_con[4]=4;
          _sons_con[5][0]=3; _sons_con[5][1]=7; _sons_con[5][2]=4; _sons_con[5][3]=0; _nb_of_sons_con[5]=4;
        }
        break;
      case NORM_QUAD4:
        {
          _nb_of_pts=4; _nb_of_sons=4; _dim=2;
          _sons_type[0]=NORM_SEG2; _sons_type[1]=NORM_SEG2; _sons_type[2]=NORM_SEG2; _sons_type[3]=NORM_SEG2;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _nb_of_sons_con[0]=2;
          _sons_con[1][0]=1; _sons_con[1][1]=2; _nb_of_sons_con[1]=2;
          _sons_con[2][0]=2; _sons_con[2][1]=3; _nb_of_sons_con[2]=2;
          _sons_con[3][0]=3; _sons_con[3][1]=0; _nb_of_sons_con[3]=2;
        }
        break;
      case NORM_TRI3:
        {
          _nb_of_pts=3; _nb_of_sons=3; _dim=2;
          _sons_type[0]=NORM_SEG2; _sons_type[1]=NORM_SEG2; _sons_type[2]=NORM_SEG2;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _nb_of_sons_con[0]=2;
          _sons_con[1][0]=1; _sons_con[1][1]=2; _nb_of_sons_con[1]=2;
          _sons_con[2][0]=2; _sons_con[2][1]=0; _nb_of_sons_con[2]=2;
        }
        break;
      case NORM_TRI6:
        {
          _nb_of_pts=6; _nb_of_sons=3; _dim=2;
          _sons_type[0]=NORM_SEG3; _sons_type[1]=NORM_SEG3; _sons_type[2]=NORM_SEG3;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=3; _nb_of_sons_con[0]=3;
          _sons_con[1][0]=1; _sons_con[1][1]=2; _sons_con[1][2]=4; _nb_of_sons_con[1]=3;
          _sons_con[2][0]=2; _sons_con[2][1]=0; _sons_con[2][2]=5; _nb_of_sons_con[2]=3; _quadratic=true;
        }
        break;
      case NORM_QUAD8:
        {
          _nb_of_pts=8; _nb_of_sons=4; _dim=2;
          _sons_type[0]=NORM_SEG3; _sons_type[1]=NORM_SEG3; _sons_type[2]=NORM_SEG3; _sons_type[3]=NORM_SEG3;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=4; _nb_of_sons_con[0]=3;
          _sons_con[1][0]=1; _sons_con[1][1]=2; _sons_con[1][2]=5; _nb_of_sons_con[1]=3;
          _sons_con[2][0]=2; _sons_con[2][1]=3; _sons_con[2][2]=6; _nb_of_sons_con[2]=3;
          _sons_con[3][0]=3; _sons_con[3][1]=0; _sons_con[3][2]=7; _nb_of_sons_con[3]=3; _quadratic=true;
        }
        break;
      case NORM_PYRA5:
        {
          _nb_of_pts=5; _nb_of_sons=5; _dim=3;
          _sons_type[0]=NORM_QUAD4; _sons_type[1]=NORM_TRI3; _sons_type[2]=NORM_TRI3; _sons_type[3]=NORM_TRI3; _sons_type[4]=NORM_TRI3;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _sons_con[0][3]=3; _nb_of_sons_con[0]=4;
          _sons_con[1][0]=0; _sons_con[1][1]=4; _sons_con[1][2]=1; _nb_of_sons_con[1]=3;
          _sons_con[2][0]=1; _sons_con[2][1]=4; _sons_con[2][2]=2; _nb_of_sons_con[2]=3;
          _sons_con[3][0]=2; _sons_con[3][1]=4; _sons_con[3][2]=3; _nb_of_sons_con[3]=3;
          _sons_con[4][0]=3; _sons_con[4][1]=4; _sons_con[4][2]=0; _nb_of_sons_con[4]=3;
        }
        break;
      case NORM_PENTA6:
        {
          _nb_of_pts=6; _nb_of_sons=5; _dim=3;
          _sons_type[0]=NORM_TRI3; _sons_type[1]=NORM_TRI3; _sons_type[2]=NORM_QUAD4; _sons_type[3]=NORM_QUAD4; _sons_type[4]=NORM_QUAD4;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _nb_of_sons_con[0]=3;
          _sons_con[1][0]=3; _sons_con[1][1]=5; _sons_con[1][2]=4; _nb_of_sons_con[1]=3;
          _sons_con[2][0]=0; _sons_con[2][1]=3; _sons_con[2][2]=4; _sons_con[2][3]=1; _nb_of_sons_con[2]=4;
          _sons_con[3][0]=1; _sons_con[3][1]=4; _sons_con[3][2]=5; _sons_con[3][3]=2; _nb_of_sons_con[3]=4;
          _sons_con[4][0]=2; _sons_con[4][1]=4; _sons_con[4][2]=5; _sons_con[4][3]=0; _nb_of_sons_con[4]=4;
        }
        break;
      case NORM_TETRA10:
        {
          _nb_of_pts=10; _nb_of_sons=4; _dim=3;
          _sons_type[0]=NORM_TRI6; _sons_type[1]=NORM_TRI6; _sons_type[2]=NORM_TRI6; _sons_type[3]=NORM_TRI6;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _sons_con[0][3]=4; _sons_con[0][4]=5; _sons_con[0][5]=6; _nb_of_sons_con[0]=6;
          _sons_con[1][0]=0; _sons_con[1][1]=3; _sons_con[1][2]=1; _sons_con[1][3]=7; _sons_con[1][4]=8; _sons_con[1][5]=4; _nb_of_sons_con[1]=6;
          _sons_con[2][0]=1; _sons_con[2][1]=3; _sons_con[2][2]=2; _sons_con[2][3]=8; _sons_con[2][4]=9; _sons_con[2][5]=5; _nb_of_sons_con[2]=6;
          _sons_con[3][0]=2; _sons_con[3][1]=3; _sons_con[3][2]=0; _sons_con[3][3]=9; _sons_con[3][4]=7; _sons_con[3][5]=6; _nb_of_sons_con[3]=6;  _quadratic=true;
        }
        break;
      case NORM_PYRA13:
        {
          _nb_of_pts=13; _nb_of_sons=5; _dim=3;
          _sons_type[0]=NORM_QUAD8; _sons_type[1]=NORM_TRI6; _sons_type[2]=NORM_TRI6; _sons_type[3]=NORM_TRI6; _sons_type[4]=NORM_TRI6;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _sons_con[0][3]=3; _sons_con[0][4]=5; _sons_con[0][5]=6; _sons_con[0][6]=7; _sons_con[0][7]=8; _nb_of_sons_con[0]=8;
          _sons_con[1][0]=0; _sons_con[1][1]=4; _sons_con[1][2]=1; _sons_con[1][3]=9; _sons_con[1][4]=10; _sons_con[1][5]=5; _nb_of_sons_con[1]=6;
          _sons_con[2][0]=1; _sons_con[2][1]=4; _sons_con[2][2]=2; _sons_con[2][3]=10; _sons_con[2][4]=11; _sons_con[2][5]=6; _nb_of_sons_con[2]=6;
          _sons_con[3][0]=2; _sons_con[3][1]=4; _sons_con[3][2]=3; _sons_con[3][3]=11; _sons_con[3][4]=12; _sons_con[3][5]=7;  _nb_of_sons_con[3]=6;
          _sons_con[4][0]=3; _sons_con[4][1]=4; _sons_con[4][2]=0; _sons_con[4][3]=12; _sons_con[4][4]=9; _sons_con[4][5]=8; _nb_of_sons_con[4]=6; _quadratic=true;
        }
        break;
      case NORM_PENTA15:
        {
          _nb_of_pts=15; _nb_of_sons=5; _dim=3;
          _sons_type[0]=NORM_TRI6; _sons_type[1]=NORM_TRI6; _sons_type[2]=NORM_QUAD8; _sons_type[3]=NORM_QUAD8; _sons_type[4]=NORM_QUAD8;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _sons_con[0][3]=6; _sons_con[0][4]=7; _sons_con[0][5]=8; _nb_of_sons_con[0]=6;
          _sons_con[1][0]=3; _sons_con[1][1]=5; _sons_con[1][2]=4; _sons_con[1][3]=11; _sons_con[1][4]=10; _sons_con[1][5]=9; _nb_of_sons_con[1]=6;
          _sons_con[2][0]=0; _sons_con[2][1]=3; _sons_con[2][2]=4; _sons_con[2][3]=1; _sons_con[2][4]=12; _sons_con[2][5]=9; _sons_con[2][6]=13; _sons_con[2][7]=6; _nb_of_sons_con[2]=8;
          _sons_con[3][0]=1; _sons_con[3][1]=4; _sons_con[3][2]=5; _sons_con[3][3]=2; _sons_con[3][4]=13; _sons_con[3][5]=10; _sons_con[3][6]=14; _sons_con[3][7]=7; _nb_of_sons_con[3]=8;
          _sons_con[4][0]=2; _sons_con[4][1]=4; _sons_con[4][2]=5; _sons_con[4][3]=0; _sons_con[4][4]=14; _sons_con[4][5]=11; _sons_con[4][6]=12; _sons_con[4][7]=8; _nb_of_sons_con[4]=8; _quadratic=true;
        }
        break;
      case NORM_HEXA20:
        {
          _nb_of_pts=20; _nb_of_sons=6; _dim=3;
          _sons_type[0]=NORM_QUAD8; _sons_type[1]=NORM_QUAD8; _sons_type[2]=NORM_QUAD8; _sons_type[3]=NORM_QUAD8; _sons_type[4]=NORM_QUAD8; _sons_type[5]=NORM_QUAD8;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _sons_con[0][3]=3; _sons_con[0][4]=8; _sons_con[0][5]=9; _sons_con[0][6]=10; _sons_con[0][7]=11; _nb_of_sons_con[0]=8;
          _sons_con[1][0]=4; _sons_con[1][1]=7; _sons_con[1][2]=6; _sons_con[1][3]=5; _sons_con[1][4]=15; _sons_con[1][5]=14; _sons_con[1][6]=13; _sons_con[1][7]=12; _nb_of_sons_con[1]=8;
          _sons_con[2][0]=0; _sons_con[2][1]=4; _sons_con[2][2]=5; _sons_con[2][3]=1; _sons_con[2][4]=16; _sons_con[2][5]=12; _sons_con[2][6]=17; _sons_con[2][7]=8; _nb_of_sons_con[2]=8;
          _sons_con[3][0]=1; _sons_con[3][1]=5; _sons_con[3][3]=6; _sons_con[3][3]=2; _sons_con[3][4]=17; _sons_con[3][5]=13; _sons_con[3][6]=18; _sons_con[3][7]=9;_nb_of_sons_con[3]=8;
          _sons_con[4][0]=2; _sons_con[4][1]=6; _sons_con[4][3]=7; _sons_con[4][3]=3; _sons_con[4][4]=18; _sons_con[4][5]=14; _sons_con[4][6]=19; _sons_con[4][7]=10; _nb_of_sons_con[4]=8;
          _sons_con[5][0]=3; _sons_con[5][1]=7; _sons_con[5][3]=4; _sons_con[5][3]=0; _sons_con[5][4]=19; _sons_con[5][5]=15; _sons_con[5][6]=16; _sons_con[5][7]=11; _nb_of_sons_con[5]=8; _quadratic=true;
        }
        break;
      case NORM_POLYGON:
        {
          _nb_of_pts=0; _nb_of_sons=0; _dim=2; _dyn=true;
        }
        break;
      case NORM_POLYHED:
        {
          _nb_of_pts=0; _nb_of_sons=0; _dim=3; _dyn=true;
        }
        break;
      case NORM_ERROR:
        {
          _nb_of_pts=std::numeric_limits<unsigned>::max(); _nb_of_sons=std::numeric_limits<unsigned>::max(); _dim=std::numeric_limits<unsigned>::max();
        }
        break;
      }
  }

  void CellModel::fillSonCellNodalConnectivity(int sonId, const int *nodalConn, int *sonNodalConn) const
  {
    unsigned nbOfTurnLoop=_nb_of_sons_con[sonId];
    const unsigned *sonConn=_sons_con[sonId];
    for(unsigned i=0;i<nbOfTurnLoop;i++)
      sonNodalConn[i]=nodalConn[sonConn[i]];
  }
}
