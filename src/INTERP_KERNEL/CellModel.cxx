//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

#include <algorithm>
#include <sstream>
#include <limits>

namespace INTERP_KERNEL
{
  const char *CellModel::CELL_TYPES_REPR[]={"NORM_POINT0", "NORM_SEG2", "NORM_SEG3", "NORM_TRI3", "NORM_QUAD4",// 0->4
                                            "NORM_POLYGON", "NORM_TRI6", "" , "NORM_QUAD8", "",//5->9
                                            "", "", "", "", "NORM_TETRA4",//10->14
                                            "NORM_PYRA5", "NORM_PENTA6", "", "NORM_HEXA8", "",//15->19
                                            "NORM_TETRA10", "", "", "NORM_PYRA13", "",//20->24
                                            "NORM_PENTA15", "", "", "", "",//25->29
                                            "NORM_HEXA20", "NORM_POLYHED", "", "", "",//30->34
                                            "", "", "", "", "",//35->39
                                            "NORM_ERROR"};

  std::map<NormalizedCellType,CellModel> CellModel::_map_of_unique_instance;

  const CellModel& CellModel::getCellModel(NormalizedCellType type)
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
    const CellModel& other=getCellModel(type);
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
    _map_of_unique_instance.insert(std::make_pair(NORM_POINT0,CellModel(NORM_POINT0)));
    _map_of_unique_instance.insert(std::make_pair(NORM_SEG2,CellModel(NORM_SEG2)));
    _map_of_unique_instance.insert(std::make_pair(NORM_SEG3,CellModel(NORM_SEG3)));
    _map_of_unique_instance.insert(std::make_pair(NORM_TRI3,CellModel(NORM_TRI3)));
    _map_of_unique_instance.insert(std::make_pair(NORM_QUAD4,CellModel(NORM_QUAD4)));
    _map_of_unique_instance.insert(std::make_pair(NORM_TRI6,CellModel(NORM_TRI6)));
    _map_of_unique_instance.insert(std::make_pair(NORM_QUAD8,CellModel(NORM_QUAD8)));
    _map_of_unique_instance.insert(std::make_pair(NORM_TETRA4,CellModel(NORM_TETRA4)));
    _map_of_unique_instance.insert(std::make_pair(NORM_HEXA8,CellModel(NORM_HEXA8)));
    _map_of_unique_instance.insert(std::make_pair(NORM_PYRA5,CellModel(NORM_PYRA5)));
    _map_of_unique_instance.insert(std::make_pair(NORM_PENTA6,CellModel(NORM_PENTA6)));
    _map_of_unique_instance.insert(std::make_pair(NORM_TETRA10,CellModel(NORM_TETRA10)));
    _map_of_unique_instance.insert(std::make_pair(NORM_PYRA13,CellModel(NORM_PYRA13)));
    _map_of_unique_instance.insert(std::make_pair(NORM_PENTA15,CellModel(NORM_PENTA15)));
    _map_of_unique_instance.insert(std::make_pair(NORM_HEXA20,CellModel(NORM_HEXA20)));
    _map_of_unique_instance.insert(std::make_pair(NORM_POLYGON,CellModel(NORM_POLYGON)));
    _map_of_unique_instance.insert(std::make_pair(NORM_POLYHED,CellModel(NORM_POLYHED)));
  }

  CellModel::CellModel(NormalizedCellType type):_type(type)
  {
    _quadratic=false;
    _dyn=false;
    _extruded_type=NORM_ERROR;
    _linear_type=NORM_ERROR;
    _quadratic_type=NORM_ERROR;
    switch(type)
      {
      case NORM_POINT0:
        {
          _nb_of_pts=1; _nb_of_sons=0; _dim=0; _extruded_type=NORM_SEG2; _is_simplex=true;
        }
        break;
      case NORM_SEG2:
        {
          _nb_of_pts=2; _nb_of_sons=0; _dim=1; _extruded_type=NORM_QUAD4; _quadratic_type=NORM_SEG3; _is_simplex=true;
        }
        break;
      case NORM_SEG3:
        {
          _nb_of_pts=3; _nb_of_sons=0; _dim=1; _extruded_type=NORM_QUAD8; _linear_type=NORM_SEG2; _is_simplex=false;
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
        }
        break;
      case NORM_HEXA8:
        {
          _nb_of_pts=8; _nb_of_sons=6; _dim=3; _quadratic_type=NORM_HEXA20; _is_simplex=false;
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
          _nb_of_pts=4; _nb_of_sons=4; _dim=2; _quadratic_type=NORM_QUAD8; _is_simplex=false;
          _sons_type[0]=NORM_SEG2; _sons_type[1]=NORM_SEG2; _sons_type[2]=NORM_SEG2; _sons_type[3]=NORM_SEG2;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _nb_of_sons_con[0]=2;
          _sons_con[1][0]=1; _sons_con[1][1]=2; _nb_of_sons_con[1]=2;
          _sons_con[2][0]=2; _sons_con[2][1]=3; _nb_of_sons_con[2]=2;
          _sons_con[3][0]=3; _sons_con[3][1]=0; _nb_of_sons_con[3]=2; _extruded_type=NORM_HEXA8;
        }
        break;
      case NORM_TRI3:
        {
          _nb_of_pts=3; _nb_of_sons=3; _dim=2; _quadratic_type=NORM_TRI6; _is_simplex=true;
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
      case NORM_PYRA5:
        {
          _nb_of_pts=5; _nb_of_sons=5; _dim=3; _quadratic_type=NORM_PYRA13; _is_simplex=false;
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
          _nb_of_pts=6; _nb_of_sons=5; _dim=3; _quadratic_type=NORM_PENTA15; _is_simplex=false;
          _sons_type[0]=NORM_TRI3; _sons_type[1]=NORM_TRI3; _sons_type[2]=NORM_QUAD4; _sons_type[3]=NORM_QUAD4; _sons_type[4]=NORM_QUAD4;
          _sons_con[0][0]=0; _sons_con[0][1]=1; _sons_con[0][2]=2; _nb_of_sons_con[0]=3;
          _sons_con[1][0]=3; _sons_con[1][1]=5; _sons_con[1][2]=4; _nb_of_sons_con[1]=3;
          _sons_con[2][0]=0; _sons_con[2][1]=3; _sons_con[2][2]=4; _sons_con[2][3]=1; _nb_of_sons_con[2]=4;
          _sons_con[3][0]=1; _sons_con[3][1]=4; _sons_con[3][2]=5; _sons_con[3][3]=2; _nb_of_sons_con[3]=4;
          _sons_con[4][0]=2; _sons_con[4][1]=5; _sons_con[4][2]=3; _sons_con[4][3]=0; _nb_of_sons_con[4]=4;
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
        }
        break;
      case NORM_HEXA20:
        {
          _nb_of_pts=20; _nb_of_sons=6; _dim=3; _linear_type=NORM_HEXA8; _is_simplex=false;
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
          _nb_of_pts=0; _nb_of_sons=0; _dim=2; _dyn=true; _extruded_type=NORM_POLYHED; _is_simplex=false;
        }
        break;
      case NORM_POLYHED:
        {
          _nb_of_pts=0; _nb_of_sons=0; _dim=3; _dyn=true; _is_simplex=false;
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
    if(_dim==2)//polygon
      return lgth;
    else
      return std::count(conn,conn+lgth,-1)+1;
  }

  /*!
   * Equivalent to getSonType except that this method deals with dynamic type.
   */
  NormalizedCellType CellModel::getSonType2(unsigned sonId) const
  {
    if(!isDynamic())
      return getSonType(sonId);
    if(_dim==2)//polygon
      return NORM_SEG2;
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
            sonNodalConn[0]=nodalConn[sonId];
            sonNodalConn[1]=nodalConn[(sonId+1)%lgth];
            return 2;
          }
        else
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
        return 2;
      }
    else
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
  }

}
