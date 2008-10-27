#include "CellModel.hxx"

#include "InterpolationUtils.hxx"

#include <sstream>

using namespace std;

namespace INTERP_KERNEL
{
  std::map<NormalizedCellType,CellModel> CellModel::_mapOfUniqueInstance;

  const CellModel& CellModel::getCellModel(NormalizedCellType type)
  {
    if(_mapOfUniqueInstance.empty())
      buildUniqueInstance();
    const map<NormalizedCellType,CellModel>::iterator iter=_mapOfUniqueInstance.find(type);
    if(iter==_mapOfUniqueInstance.end())
      {
        ostringstream stream; stream << "no cellmodel for normalized type " << type;
        throw Exception(stream.str().c_str());
      }
    return (*iter).second;
  }

  void CellModel::buildUniqueInstance()
  {
    _mapOfUniqueInstance.insert(make_pair(NORM_TRI3,CellModel(NORM_TRI3)));
    _mapOfUniqueInstance.insert(make_pair(NORM_QUAD4,CellModel(NORM_QUAD4)));
    _mapOfUniqueInstance.insert(make_pair(NORM_TRI6,CellModel(NORM_TRI6)));
    _mapOfUniqueInstance.insert(make_pair(NORM_QUAD8,CellModel(NORM_QUAD8)));
    _mapOfUniqueInstance.insert(make_pair(NORM_TETRA4,CellModel(NORM_TETRA4)));
    _mapOfUniqueInstance.insert(make_pair(NORM_HEXA8,CellModel(NORM_HEXA8)));
    _mapOfUniqueInstance.insert(make_pair(NORM_PYRA5,CellModel(NORM_PYRA5)));
    _mapOfUniqueInstance.insert(make_pair(NORM_PENTA6,CellModel(NORM_PENTA6)));
    _mapOfUniqueInstance.insert(make_pair(NORM_TETRA10,CellModel(NORM_TETRA10)));
    _mapOfUniqueInstance.insert(make_pair(NORM_PYRA13,CellModel(NORM_PYRA13)));
    _mapOfUniqueInstance.insert(make_pair(NORM_PENTA15,CellModel(NORM_PENTA15)));
    _mapOfUniqueInstance.insert(make_pair(NORM_HEXA20,CellModel(NORM_HEXA20)));
  }

  CellModel::CellModel(NormalizedCellType type)
  {
    switch(type)
      {
      case NORM_TETRA4:
        {
          _nbOfPts=4; _nbOfSons=4; _dim=3;
          _sonsType[0]=NORM_TRI3; _sonsType[1]=NORM_TRI3; _sonsType[2]=NORM_TRI3; _sonsType[3]=NORM_TRI3;
          _sonsCon[0][0]=0; _sonsCon[0][1]=1; _sonsCon[0][2]=2; _nbOfSonsCon[0]=3;
          _sonsCon[1][0]=0; _sonsCon[1][1]=3; _sonsCon[1][2]=1; _nbOfSonsCon[1]=3;
          _sonsCon[2][0]=1; _sonsCon[2][1]=3; _sonsCon[2][2]=2; _nbOfSonsCon[2]=3;
          _sonsCon[3][0]=2; _sonsCon[3][1]=3; _sonsCon[3][2]=0; _nbOfSonsCon[3]=3;
        }
        break;
      case NORM_HEXA8:
        {
          _nbOfPts=8; _nbOfSons=6; _dim=3;
          _sonsType[0]=NORM_QUAD4; _sonsType[1]=NORM_QUAD4; _sonsType[2]=NORM_QUAD4; _sonsType[3]=NORM_QUAD4; _sonsType[4]=NORM_QUAD4; _sonsType[5]=NORM_QUAD4;
          _sonsCon[0][0]=0; _sonsCon[0][1]=1; _sonsCon[0][2]=2; _sonsCon[0][3]=3; _nbOfSonsCon[0]=4;
          _sonsCon[1][0]=4; _sonsCon[1][1]=7; _sonsCon[1][2]=6; _sonsCon[1][3]=5; _nbOfSonsCon[1]=4;
          _sonsCon[2][0]=0; _sonsCon[2][1]=4; _sonsCon[2][2]=5; _sonsCon[2][3]=1; _nbOfSonsCon[2]=4;
          _sonsCon[3][0]=1; _sonsCon[3][1]=5; _sonsCon[3][2]=6; _sonsCon[3][3]=2; _nbOfSonsCon[3]=4;
          _sonsCon[4][0]=2; _sonsCon[4][1]=6; _sonsCon[4][2]=7; _sonsCon[4][3]=3; _nbOfSonsCon[4]=4;
          _sonsCon[5][0]=3; _sonsCon[5][1]=7; _sonsCon[5][2]=4; _sonsCon[5][3]=0; _nbOfSonsCon[5]=4;
        }
        break;
      case NORM_QUAD4:
        {
          _nbOfPts=4; _nbOfSons=4; _dim=2;
          _sonsType[0]=NORM_SEG2; _sonsType[1]=NORM_SEG2; _sonsType[2]=NORM_SEG2; _sonsType[3]=NORM_SEG2;
          _sonsCon[0][0]=0; _sonsCon[0][1]=1; _nbOfSonsCon[0]=2;
          _sonsCon[1][0]=1; _sonsCon[1][1]=2; _nbOfSonsCon[1]=2;
          _sonsCon[2][0]=2; _sonsCon[2][1]=3; _nbOfSonsCon[2]=2;
          _sonsCon[3][0]=3; _sonsCon[3][1]=0; _nbOfSonsCon[3]=2;
        }
        break;
      case NORM_TRI3:
        {
          _nbOfPts=3; _nbOfSons=3; _dim=2;
          _sonsType[0]=NORM_SEG2; _sonsType[1]=NORM_SEG2; _sonsType[2]=NORM_SEG2;
          _sonsCon[0][0]=0; _sonsCon[0][1]=1; _nbOfSonsCon[0]=2;
          _sonsCon[1][0]=1; _sonsCon[1][1]=2; _nbOfSonsCon[1]=2;
          _sonsCon[2][0]=2; _sonsCon[2][1]=0; _nbOfSonsCon[2]=2;
        }
        break;
      case NORM_TRI6:
        {
          _nbOfPts=6; _nbOfSons=3; _dim=2;
          _sonsType[0]=NORM_SEG3; _sonsType[1]=NORM_SEG3; _sonsType[2]=NORM_SEG3;
          _sonsCon[0][0]=0; _sonsCon[0][1]=1; _sonsCon[0][2]=3; _nbOfSonsCon[0]=3;
          _sonsCon[1][0]=1; _sonsCon[1][1]=2; _sonsCon[1][2]=4; _nbOfSonsCon[1]=3;
          _sonsCon[2][0]=2; _sonsCon[2][1]=0; _sonsCon[2][2]=5; _nbOfSonsCon[2]=3;
        }
        break;
      case NORM_QUAD8:
        {
          _nbOfPts=8; _nbOfSons=4; _dim=2;
          _sonsType[0]=NORM_SEG3; _sonsType[1]=NORM_SEG3; _sonsType[2]=NORM_SEG3; _sonsType[3]=NORM_SEG3;
          _sonsCon[0][0]=0; _sonsCon[0][1]=1; _sonsCon[0][2]=4; _nbOfSonsCon[0]=3;
          _sonsCon[1][0]=1; _sonsCon[1][1]=2; _sonsCon[1][2]=5; _nbOfSonsCon[1]=3;
          _sonsCon[2][0]=2; _sonsCon[2][1]=3; _sonsCon[2][2]=6; _nbOfSonsCon[2]=3;
          _sonsCon[3][0]=3; _sonsCon[3][1]=0; _sonsCon[3][2]=7; _nbOfSonsCon[3]=3;
        }
        break;
      case NORM_PYRA5:
        {
          _nbOfPts=5; _nbOfSons=5; _dim=3;
          _sonsType[0]=NORM_QUAD4; _sonsType[1]=NORM_TRI3; _sonsType[2]=NORM_TRI3; _sonsType[3]=NORM_TRI3; _sonsType[4]=NORM_TRI3;
          _sonsCon[0][0]=0; _sonsCon[0][1]=1; _sonsCon[0][2]=2; _sonsCon[0][3]=3; _nbOfSonsCon[0]=4;
          _sonsCon[1][0]=0; _sonsCon[1][1]=4; _sonsCon[1][2]=1; _nbOfSonsCon[1]=3;
          _sonsCon[2][0]=1; _sonsCon[2][1]=4; _sonsCon[2][2]=2; _nbOfSonsCon[2]=3;
          _sonsCon[3][0]=2; _sonsCon[3][1]=4; _sonsCon[3][2]=3; _nbOfSonsCon[3]=3;
          _sonsCon[4][0]=3; _sonsCon[4][1]=4; _sonsCon[4][2]=0; _nbOfSonsCon[4]=3;
        }
        break;
      case NORM_PENTA6:
        {
          _nbOfPts=6; _nbOfSons=5; _dim=3;
          _sonsType[0]=NORM_TRI3; _sonsType[1]=NORM_TRI3; _sonsType[2]=NORM_QUAD4; _sonsType[3]=NORM_QUAD4; _sonsType[4]=NORM_QUAD4;
          _sonsCon[0][0]=0; _sonsCon[0][1]=1; _sonsCon[0][2]=2; _nbOfSonsCon[0]=3;
          _sonsCon[1][0]=3; _sonsCon[1][1]=5; _sonsCon[1][2]=4; _nbOfSonsCon[1]=3;
          _sonsCon[2][0]=0; _sonsCon[2][1]=3; _sonsCon[2][2]=4; _sonsCon[2][3]=1; _nbOfSonsCon[2]=4;
          _sonsCon[3][0]=1; _sonsCon[3][1]=4; _sonsCon[3][2]=5; _sonsCon[3][3]=2; _nbOfSonsCon[3]=4;
          _sonsCon[4][0]=2; _sonsCon[4][1]=4; _sonsCon[4][2]=5; _sonsCon[4][3]=0; _nbOfSonsCon[4]=4;
        }
        break;
      case NORM_TETRA10:
        {
          _nbOfPts=10; _nbOfSons=4; _dim=3;
          _sonsType[0]=NORM_TRI6; _sonsType[1]=NORM_TRI6; _sonsType[2]=NORM_TRI6; _sonsType[3]=NORM_TRI6;
          _sonsCon[0][0]=0; _sonsCon[0][1]=1; _sonsCon[0][2]=2; _sonsCon[0][3]=4; _sonsCon[0][4]=5; _sonsCon[0][5]=6; _nbOfSonsCon[0]=6;
          _sonsCon[1][0]=0; _sonsCon[1][1]=3; _sonsCon[1][2]=1; _sonsCon[1][3]=7; _sonsCon[1][4]=8; _sonsCon[1][5]=4; _nbOfSonsCon[1]=6;
          _sonsCon[2][0]=1; _sonsCon[2][1]=3; _sonsCon[2][2]=2; _sonsCon[2][3]=8; _sonsCon[2][4]=9; _sonsCon[2][5]=5; _nbOfSonsCon[2]=6;
          _sonsCon[3][0]=2; _sonsCon[3][1]=3; _sonsCon[3][2]=0; _sonsCon[3][3]=9; _sonsCon[3][4]=7; _sonsCon[3][5]=6; _nbOfSonsCon[3]=6;
        }
        break;
      case NORM_PYRA13:
        {
          _nbOfPts=13; _nbOfSons=5; _dim=3;
          _sonsType[0]=NORM_QUAD8; _sonsType[1]=NORM_TRI6; _sonsType[2]=NORM_TRI6; _sonsType[3]=NORM_TRI6; _sonsType[4]=NORM_TRI6;
          _sonsCon[0][0]=0; _sonsCon[0][1]=1; _sonsCon[0][2]=2; _sonsCon[0][3]=3; _sonsCon[0][4]=5; _sonsCon[0][5]=6; _sonsCon[0][6]=7; _sonsCon[0][7]=8; _nbOfSonsCon[0]=8;
          _sonsCon[1][0]=0; _sonsCon[1][1]=4; _sonsCon[1][2]=1; _sonsCon[1][3]=9; _sonsCon[1][4]=10; _sonsCon[1][5]=5; _nbOfSonsCon[1]=6;
          _sonsCon[2][0]=1; _sonsCon[2][1]=4; _sonsCon[2][2]=2; _sonsCon[2][3]=10; _sonsCon[2][4]=11; _sonsCon[2][5]=6; _nbOfSonsCon[2]=6;
          _sonsCon[3][0]=2; _sonsCon[3][1]=4; _sonsCon[3][2]=3; _sonsCon[3][3]=11; _sonsCon[3][4]=12; _sonsCon[3][5]=7;  _nbOfSonsCon[3]=6;
          _sonsCon[4][0]=3; _sonsCon[4][1]=4; _sonsCon[4][2]=0; _sonsCon[4][3]=12; _sonsCon[4][4]=9; _sonsCon[4][5]=8; _nbOfSonsCon[4]=6;
        }
        break;
      case NORM_PENTA15:
        {
          _nbOfPts=15; _nbOfSons=5; _dim=3;
          _sonsType[0]=NORM_TRI6; _sonsType[1]=NORM_TRI6; _sonsType[2]=NORM_QUAD8; _sonsType[3]=NORM_QUAD8; _sonsType[4]=NORM_QUAD8;
          _sonsCon[0][0]=0; _sonsCon[0][1]=1; _sonsCon[0][2]=2; _sonsCon[0][3]=6; _sonsCon[0][4]=7; _sonsCon[0][5]=8; _nbOfSonsCon[0]=6;
          _sonsCon[1][0]=3; _sonsCon[1][1]=5; _sonsCon[1][2]=4; _sonsCon[1][3]=11; _sonsCon[1][4]=10; _sonsCon[1][5]=9; _nbOfSonsCon[1]=6;
          _sonsCon[2][0]=0; _sonsCon[2][1]=3; _sonsCon[2][2]=4; _sonsCon[2][3]=1; _sonsCon[2][4]=12; _sonsCon[2][5]=9; _sonsCon[2][6]=13; _sonsCon[2][7]=6; _nbOfSonsCon[2]=8;
          _sonsCon[3][0]=1; _sonsCon[3][1]=4; _sonsCon[3][2]=5; _sonsCon[3][3]=2; _sonsCon[3][4]=13; _sonsCon[3][5]=10; _sonsCon[3][6]=14; _sonsCon[3][7]=7; _nbOfSonsCon[3]=8;
          _sonsCon[4][0]=2; _sonsCon[4][1]=4; _sonsCon[4][2]=5; _sonsCon[4][3]=0; _sonsCon[4][4]=14; _sonsCon[4][5]=11; _sonsCon[4][6]=12; _sonsCon[4][7]=8; _nbOfSonsCon[4]=8;
        }
        break;
      case NORM_HEXA20:
        {
          _nbOfPts=20; _nbOfSons=6; _dim=3;
          _sonsType[0]=NORM_QUAD8; _sonsType[1]=NORM_QUAD8; _sonsType[2]=NORM_QUAD8; _sonsType[3]=NORM_QUAD8; _sonsType[4]=NORM_QUAD8; _sonsType[5]=NORM_QUAD8;
          _sonsCon[0][0]=0; _sonsCon[0][1]=1; _sonsCon[0][2]=2; _sonsCon[0][3]=3; _sonsCon[0][4]=8; _sonsCon[0][5]=9; _sonsCon[0][6]=10; _sonsCon[0][7]=11; _nbOfSonsCon[0]=8;
          _sonsCon[1][0]=4; _sonsCon[1][1]=7; _sonsCon[1][2]=6; _sonsCon[1][3]=5; _sonsCon[1][4]=15; _sonsCon[1][5]=14; _sonsCon[1][6]=13; _sonsCon[1][7]=12; _nbOfSonsCon[1]=8;
          _sonsCon[2][0]=0; _sonsCon[2][1]=4; _sonsCon[2][2]=5; _sonsCon[2][3]=1; _sonsCon[2][4]=16; _sonsCon[2][5]=12; _sonsCon[2][6]=17; _sonsCon[2][7]=8; _nbOfSonsCon[2]=8;
          _sonsCon[3][0]=1; _sonsCon[3][1]=5; _sonsCon[3][3]=6; _sonsCon[3][3]=2; _sonsCon[3][4]=17; _sonsCon[3][5]=13; _sonsCon[3][6]=18; _sonsCon[3][7]=9;_nbOfSonsCon[3]=8;
          _sonsCon[4][0]=2; _sonsCon[4][1]=6; _sonsCon[4][3]=7; _sonsCon[4][3]=3; _sonsCon[4][4]=18; _sonsCon[4][5]=14; _sonsCon[4][6]=19; _sonsCon[4][7]=10; _nbOfSonsCon[4]=8;
          _sonsCon[5][0]=3; _sonsCon[5][1]=7; _sonsCon[5][3]=4; _sonsCon[5][3]=0; _sonsCon[5][4]=19; _sonsCon[5][5]=15; _sonsCon[5][6]=16; _sonsCon[5][7]=11; _nbOfSonsCon[5]=8;
        }
        break;
      }
  }

  void CellModel::fillSonCellNodalConnectivity(int sonId, const int *nodalConn, int *sonNodalConn) const
  {
    unsigned nbOfTurnLoop=_nbOfSonsCon[sonId];
    const unsigned *sonConn=_sonsCon[sonId];
    for(unsigned i=0;i<nbOfTurnLoop;i++)
      sonNodalConn[i]=nodalConn[sonConn[i]];
  }
}
