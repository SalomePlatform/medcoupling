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
// Author : Anthony Geay (CEA/DEN)

#include "InterpKernelCellSimplify.hxx"
#include "CellModel.hxx"

#include <algorithm>
#include <iterator>
#include <sstream>
#include <numeric>
#include <cstring>
#include <limits>
#include <vector>
#include <list>
#include <set>

using namespace INTERP_KERNEL;

/*!
 * This method takes as input a cell with type 'type' and whose connectivity is defined by (conn,lgth)
 * It retrieves the same cell with a potentially different type (in return) whose connectivity is defined by (retConn,retLgth)
 * \b WARNING for optimization reason the arrays 'retConn' and 'conn' can overlapped !
 */
INTERP_KERNEL::NormalizedCellType CellSimplify::simplifyDegeneratedCell(INTERP_KERNEL::NormalizedCellType type, const int *conn, int lgth, int *retConn, int& retLgth)
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  std::set<int> c(conn,conn+lgth);
  c.erase(-1);
  bool isObviousNonDegeneratedCell=((int)c.size()==lgth);
  if(cm.isQuadratic() || isObviousNonDegeneratedCell)
    {//quadratic do nothing for the moment.
      retLgth=lgth;
      int *tmp=new int[lgth];//no direct std::copy ! overlapping of conn and retConn !
      std::copy(conn,conn+lgth,tmp);
      std::copy(tmp,tmp+lgth,retConn);
      delete [] tmp;
      return type;
    }
  if(cm.getDimension()==2)
    {
      int *tmp=new int[lgth];
      tmp[0]=conn[0];
      int newPos=1;
      for(int i=1;i<lgth;i++)
        if(std::find(tmp,tmp+newPos,conn[i])==tmp+newPos)
          tmp[newPos++]=conn[i];
      INTERP_KERNEL::NormalizedCellType ret=tryToUnPoly2D(cm.isQuadratic(),tmp,newPos,retConn,retLgth);
      delete [] tmp;
      return ret;
    }
  if(cm.getDimension()==3)
    {
      int nbOfFaces,lgthOfPolyhConn;
      int *zipFullReprOfPolyh=getFullPolyh3DCell(type,conn,lgth,nbOfFaces,lgthOfPolyhConn);
      INTERP_KERNEL::NormalizedCellType ret=tryToUnPoly3D(zipFullReprOfPolyh,nbOfFaces,lgthOfPolyhConn,retConn,retLgth);
      delete [] zipFullReprOfPolyh;
      return ret;
    }
  throw INTERP_KERNEL::Exception("CellSimplify::simplifyDegeneratedCell : works only with 2D and 3D cell !");
}


/*!
 * This static method tries to unpolygonize a cell whose connectivity is given by 'conn' and 'lgth'.
 * Contrary to INTERP_KERNEL::CellSimplify::simplifyDegeneratedCell method 'conn' and 'retConn' do not overlap. 
 */
INTERP_KERNEL::NormalizedCellType CellSimplify::tryToUnPoly2D(bool isQuad, const int *conn, int lgth, int *retConn, int& retLgth)
{
  retLgth=lgth;
  std::copy(conn,conn+lgth,retConn);
  if(!isQuad)
    {
      switch(lgth)
        {
        case 3:
          return INTERP_KERNEL::NORM_TRI3;
        case 4:
          return INTERP_KERNEL::NORM_QUAD4;
        default:
          return INTERP_KERNEL::NORM_POLYGON;
        }
    }
  else
    {
      switch(lgth)
        {
          case 6:
            return INTERP_KERNEL::NORM_TRI6;
          case 8:
            return INTERP_KERNEL::NORM_QUAD8;
          default:
            return INTERP_KERNEL::NORM_QPOLYG;
        }
    }
}

/*!
 * This method takes as input a 3D linear cell and put its representation in returned array. Warning the returned array has to be deallocated.
 * The length of the returned array is specified by out parameter
 * The format of output array is the following :
 * 1,2,3,-1,3,4,2,-1,3,4,1,-1,1,2,4,NORM_TRI3,NORM_TRI3,NORM_TRI3 (faces type at the end of classical polyhedron nodal description)
 */
int *CellSimplify::getFullPolyh3DCell(INTERP_KERNEL::NormalizedCellType type, const int *conn, int lgth,
                                      int& retNbOfFaces, int& retLgth)
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  unsigned nbOfFaces=cm.getNumberOfSons2(conn,lgth);
  int *tmp=new int[nbOfFaces*(lgth+1)];
  int *work=tmp;
  std::vector<int> faces;
  for(unsigned j=0;j<nbOfFaces;j++)
    {
      INTERP_KERNEL::NormalizedCellType type2;
      unsigned offset=cm.fillSonCellNodalConnectivity2(j,conn,lgth,work,type2);
      //
      int *tmp2=new int[offset];
      tmp2[0]=work[0];
      int newPos=1;
      for(unsigned k=1;k<offset;k++)
        if(std::find(tmp2,tmp2+newPos,work[k])==tmp2+newPos)
          tmp2[newPos++]=work[k];
      if(newPos<3)
        {
          delete [] tmp2;
          continue;
        }
      int tmp3;
      faces.push_back(tryToUnPoly2D(CellModel::GetCellModel(type2).isQuadratic(),tmp2,newPos,work,tmp3));
      delete [] tmp2;
      //
      work+=newPos;
      *work++=-1;
    }
  std::copy(faces.begin(),faces.end(),--work);
  retNbOfFaces=(int)faces.size();
  retLgth=(int)std::distance(tmp,work);
  return tmp;
}

/*!
 * This static method tries to unpolygonize a cell whose connectivity is given by 'conn' (format is the same as specified in
 * method INTERP_KERNEL::CellSimplify::getFullPolyh3DCell ) and 'lgth'+'nbOfFaces'.
 * Contrary to INTERP_KERNEL::CellSimplify::simplifyDegeneratedCell method 'conn' and 'retConn' do not overlap. 
 */
INTERP_KERNEL::NormalizedCellType CellSimplify::tryToUnPoly3D(const int *conn, int nbOfFaces, int lgth, int *retConn, int& retLgth)
{
  std::set<int> nodes(conn,conn+lgth);
  nodes.erase(-1);
  int nbOfNodes=(int)nodes.size();
  int magicNumber=100*nbOfNodes+nbOfFaces;
  switch(magicNumber)
    {
    case 806:
      return tryToUnPolyHex8(conn,nbOfFaces,lgth,retConn,retLgth);
    case 1208:
      return tryToUnPolyHexp12(conn,nbOfFaces,lgth,retConn,retLgth);
    case 605:
      return tryToUnPolyPenta6(conn,nbOfFaces,lgth,retConn,retLgth);
    case 505:
      return tryToUnPolyPyra5(conn,nbOfFaces,lgth,retConn,retLgth);
    case 404:
      return tryToUnPolyTetra4(conn,nbOfFaces,lgth,retConn,retLgth);
    default:
      retLgth=lgth;
      std::copy(conn,conn+lgth,retConn);
      return INTERP_KERNEL::NORM_POLYHED;
    }
}

bool CellSimplify::orientOppositeFace(const int *baseFace, int *retConn, const int *sideFace, int lgthBaseFace)
{
  std::vector<int> tmp2;
  std::set<int> bases(baseFace,baseFace+lgthBaseFace);
  std::set<int> sides(sideFace,sideFace+4);
  std::set_intersection(bases.begin(),bases.end(),sides.begin(),sides.end(),std::back_insert_iterator< std::vector<int> >(tmp2));
  if(tmp2.size()!=2)
    return false;
  std::vector< std::pair<int,int> > baseEdges(lgthBaseFace);
  std::vector< std::pair<int,int> > oppEdges(lgthBaseFace);
  std::vector< std::pair<int,int> > sideEdges(4);
  for(int i=0;i<lgthBaseFace;i++)
    {
      baseEdges[i]=std::pair<int,int>(baseFace[i],baseFace[(i+1)%lgthBaseFace]);
      oppEdges[i]=std::pair<int,int>(retConn[i],retConn[(i+1)%lgthBaseFace]);
    }
  for(int i=0;i<4;i++)
    sideEdges[i]=std::pair<int,int>(sideFace[i],sideFace[(i+1)%4]);
  std::vector< std::pair<int,int> > tmp;
  std::set< std::pair<int,int> > baseEdgesS(baseEdges.begin(),baseEdges.end());
  std::set< std::pair<int,int> > sideEdgesS(sideEdges.begin(),sideEdges.end());
  std::set_intersection(baseEdgesS.begin(),baseEdgesS.end(),sideEdgesS.begin(),sideEdgesS.end(),std::back_insert_iterator< std::vector< std::pair<int,int> > >(tmp));
  if(tmp.empty())
    {
      //reverse sideFace
      for(int i=0;i<4;i++)
        {
          std::pair<int,int> p=sideEdges[i];
          std::pair<int,int> r(p.second,p.first);
          sideEdges[i]=r;
        }
      //end reverse sideFace
      std::set< std::pair<int,int> > baseEdgesS2(baseEdges.begin(),baseEdges.end());
      std::set< std::pair<int,int> > sideEdgesS2(sideEdges.begin(),sideEdges.end());
      std::set_intersection(baseEdgesS2.begin(),baseEdgesS2.end(),sideEdgesS2.begin(),sideEdgesS2.end(),std::back_insert_iterator< std::vector< std::pair<int,int> > >(tmp));
      if(tmp.empty())
        return false;
    }
  if(tmp.size()!=1)
    return false;
  bool found=false;
  std::pair<int,int> pInOpp;
  for(int i=0;i<4 && !found;i++)
    {//finding the pair(edge) in sideFace that do not include any node of tmp[0] edge
      found=(tmp[0].first!=sideEdges[i].first && tmp[0].first!=sideEdges[i].second &&
             tmp[0].second!=sideEdges[i].first && tmp[0].second!=sideEdges[i].second);
      if(found)
        {//found ! reverse it
          pInOpp.first=sideEdges[i].second;
          pInOpp.second=sideEdges[i].first;
        }
    }
  if(!found)
    return false;
  int pos=(int)std::distance(baseEdges.begin(),std::find(baseEdges.begin(),baseEdges.end(),tmp[0]));
  std::vector< std::pair<int,int> >::iterator it=std::find(oppEdges.begin(),oppEdges.end(),pInOpp);
  if(it==oppEdges.end())//the opposite edge of side face is not found opposite face ... maybe problem of orientation of polyhedron
    return false;
  int pos2=(int)std::distance(oppEdges.begin(),it);
  int offset=pos-pos2;
  if(offset<0)
    offset+=lgthBaseFace;
  //this is the end copy the result
  int *tmp3=new int[lgthBaseFace];
  for(int i=0;i<lgthBaseFace;i++)
    tmp3[(offset+i)%lgthBaseFace]=oppEdges[i].first;
  std::copy(tmp3,tmp3+lgthBaseFace,retConn);
  delete [] tmp3;
  return true;
}

bool CellSimplify::isWellOriented(const int *baseFace, int *retConn, const int *sideFace, int lgthBaseFace)
{
  return true;
}

/*!
 * This method is trying to permute the connectivity of 'oppFace' face so that the k_th node of 'baseFace' is associated to the 
 * k_th node in retConnOfOppFace. Excluded faces 'baseFace' and 'oppFace' all the other faces in 'conn' must be QUAD4 faces.
 * If the arragement process succeeds true is returned and retConnOfOppFace is filled.
 */
bool CellSimplify::tryToArrangeOppositeFace(const int *conn, int lgth, int lgthBaseFace, const int *baseFace, const int *oppFace, int nbOfFaces, int *retConnOfOppFace)
{
  retConnOfOppFace[0]=oppFace[0];
  for(int j=1;j<lgthBaseFace;j++)
    retConnOfOppFace[j]=oppFace[lgthBaseFace-j];
  const int *curFace=conn;
  int sideFace=0;
  bool ret=true;
  for(int i=0;i<nbOfFaces && ret;i++)
    {
      if(curFace!=baseFace && curFace!=oppFace)
        {
          if(sideFace==0)
            ret=orientOppositeFace(baseFace,retConnOfOppFace,curFace,lgthBaseFace);
          else
            ret=isWellOriented(baseFace,retConnOfOppFace,curFace,lgthBaseFace);
          sideFace++;
        }
      curFace=std::find(curFace,conn+lgth,-1);
      curFace++;
    }
  return ret;
}

/*!
 * Cell with 'conn' connectivity has been detected as a good candidate. Full check of this. If yes NORM_HEXA8 is returned.
 * This method is only callable if in 'conn' there is 8 nodes and 6 faces.
 * If fails a POLYHED is returned. 
 */
INTERP_KERNEL::NormalizedCellType CellSimplify::tryToUnPolyHex8(const int *conn, int nbOfFaces, int lgth, int *retConn, int& retLgth)
{
  if(std::find_if(conn+lgth,conn+lgth+nbOfFaces,std::bind2nd(std::not_equal_to<int>(),(int)INTERP_KERNEL::NORM_QUAD4))==conn+lgth+nbOfFaces)
    {//6 faces are QUAD4.
      int oppositeFace=-1;
      std::set<int> conn1(conn,conn+4);
      for(int i=1;i<6 && oppositeFace<0;i++)
        {
          std::vector<int> tmp;
          std::set<int> conn2(conn+5*i,conn+5*i+4);
          std::set_intersection(conn1.begin(),conn1.end(),conn2.begin(),conn2.end(),std::back_insert_iterator< std::vector<int> >(tmp));
          if(tmp.empty())
            oppositeFace=i;
        }
      if(oppositeFace>=1)
        {//oppositeFace of face#0 found.
          int tmp2[4];
          if(tryToArrangeOppositeFace(conn,lgth,4,conn,conn+5*oppositeFace,6,tmp2))
            {
              std::copy(conn,conn+4,retConn);
              std::copy(tmp2,tmp2+4,retConn+4);
              retLgth=8;
              return INTERP_KERNEL::NORM_HEXA8;
            }
        }
    }
  retLgth=lgth;
  std::copy(conn,conn+lgth,retConn);
  return INTERP_KERNEL::NORM_POLYHED;
}

INTERP_KERNEL::NormalizedCellType CellSimplify::tryToUnPolyHexp12(const int *conn, int nbOfFaces, int lgth, int *retConn, int& retLgth)
{
  std::size_t nbOfHexagon=std::count(conn+lgth,conn+lgth+nbOfFaces,(int)INTERP_KERNEL::NORM_POLYGON);
  std::size_t nbOfQuad=std::count(conn+lgth,conn+lgth+nbOfFaces,(int)INTERP_KERNEL::NORM_QUAD4);
  if(nbOfQuad==6 && nbOfHexagon==2)
    {
      const int *hexag0=std::find(conn+lgth,conn+lgth+nbOfFaces,(int)INTERP_KERNEL::NORM_POLYGON);
      std::size_t hexg0Id=std::distance(conn+lgth,hexag0);
      const int *hexag1=std::find(hexag0+1,conn+lgth+nbOfFaces,(int)INTERP_KERNEL::NORM_POLYGON);
      std::size_t hexg1Id=std::distance(conn+lgth,hexag1);
      const int *connHexag0=conn+5*hexg0Id;
      std::size_t lgthH0=std::distance(connHexag0,std::find(connHexag0,conn+lgth,-1));
      if(lgthH0==6)
        {
          const int *connHexag1=conn+5*hexg0Id+7+(hexg1Id-hexg0Id-1)*5;
          std::size_t lgthH1=std::distance(connHexag1,std::find(connHexag1,conn+lgth,-1));
          if(lgthH1==6)
            {
              std::vector<int> tmp;
              std::set<int> conn1(connHexag0,connHexag0+6);
              std::set<int> conn2(connHexag1,connHexag1+6);
              std::set_intersection(conn1.begin(),conn1.end(),conn2.begin(),conn2.end(),std::back_insert_iterator< std::vector<int> >(tmp));
              if(tmp.empty())
                {
                  int tmp2[6];
                  if(tryToArrangeOppositeFace(conn,lgth,6,connHexag0,connHexag1,8,tmp2))
                    {
                      std::copy(connHexag0,connHexag0+6,retConn);
                      std::copy(tmp2,tmp2+6,retConn+6);
                      retLgth=12;
                      return INTERP_KERNEL::NORM_HEXGP12;
                    }
                }
            }
        }
    }
  retLgth=lgth;
  std::copy(conn,conn+lgth,retConn);
  return INTERP_KERNEL::NORM_POLYHED;
}

/*!
 * Cell with 'conn' connectivity has been detected as a good candidate. Full check of this. If yes NORM_PENTA6 is returned.
 * If fails a POLYHED is returned. 
 */
INTERP_KERNEL::NormalizedCellType CellSimplify::tryToUnPolyPenta6(const int *conn, int nbOfFaces, int lgth, int *retConn, int& retLgth)
{
  std::size_t nbOfTriFace=std::count(conn+lgth,conn+lgth+nbOfFaces,(int)INTERP_KERNEL::NORM_TRI3);
  std::size_t nbOfQuadFace=std::count(conn+lgth,conn+lgth+nbOfFaces,(int)INTERP_KERNEL::NORM_QUAD4);
  if(nbOfTriFace==2 && nbOfQuadFace==3)
    {
      std::size_t tri3_0=std::distance(conn+lgth,std::find(conn+lgth,conn+lgth+nbOfFaces,(int)INTERP_KERNEL::NORM_TRI3));
      std::size_t tri3_1=std::distance(conn+lgth,std::find(conn+lgth+tri3_0+1,conn+lgth+nbOfFaces,(int)INTERP_KERNEL::NORM_TRI3));
      const int *tri_0=0,*tri_1=0;
      const int *w=conn;
      for(std::size_t i=0;i<5;i++)
        {
          if(i==tri3_0)
            tri_0=w;
          if(i==tri3_1)
            tri_1=w;
          w=std::find(w,conn+lgth,-1);
          w++;
        }
      std::vector<int> tmp;
      std::set<int> conn1(tri_0,tri_0+3);
      std::set<int> conn2(tri_1,tri_1+3);
      std::set_intersection(conn1.begin(),conn1.end(),conn2.begin(),conn2.end(),std::back_insert_iterator< std::vector<int> >(tmp));
      if(tmp.empty())
        {
          int tmp2[3];
          if(tryToArrangeOppositeFace(conn,lgth,3,tri_0,tri_1,5,tmp2))
            {
              std::copy(tri_0,tri_0+3,retConn);
              std::copy(tmp2,tmp2+3,retConn+3);
              retLgth=6;
              return INTERP_KERNEL::NORM_PENTA6;
            }
        }
    }
  retLgth=lgth;
  std::copy(conn,conn+lgth,retConn);
  return INTERP_KERNEL::NORM_POLYHED;
}

/*!
 * Cell with 'conn' connectivity has been detected as a good candidate. Full check of this. If yes NORM_PYRA5 is returned.
 * If fails a POLYHED is returned. 
 */
INTERP_KERNEL::NormalizedCellType CellSimplify::tryToUnPolyPyra5(const int *conn, int nbOfFaces, int lgth, int *retConn, int& retLgth)
{
  std::size_t nbOfTriFace=std::count(conn+lgth,conn+lgth+nbOfFaces,(int)INTERP_KERNEL::NORM_TRI3);
  std::size_t nbOfQuadFace=std::count(conn+lgth,conn+lgth+nbOfFaces,(int)INTERP_KERNEL::NORM_QUAD4);
  if(nbOfTriFace==4 && nbOfQuadFace==1)
    {
      std::size_t quad4_pos=std::distance(conn+lgth,std::find(conn+lgth,conn+lgth+nbOfFaces,(int)INTERP_KERNEL::NORM_QUAD4));
      const int *quad4=0;
      const int *w=conn;
      for(std::size_t i=0;i<5 && quad4==0;i++)
        {
          if(i==quad4_pos)
            quad4=w;
          w=std::find(w,conn+lgth,-1);
          w++;
        }
      std::set<int> quad4S(quad4,quad4+4);
      w=conn;
      bool ok=true;
      int point=-1;
      for(std::size_t i=0;i<5 && ok;i++)
        {
          if(i!=quad4_pos)
            {
              std::vector<int> tmp;
              std::set<int> conn2(w,w+3);
              std::set_intersection(conn2.begin(),conn2.end(),quad4S.begin(),quad4S.end(),std::back_insert_iterator< std::vector<int> >(tmp));
              ok=tmp.size()==2;
              tmp.clear();
              std::set_difference(conn2.begin(),conn2.end(),quad4S.begin(),quad4S.end(),std::back_insert_iterator< std::vector<int> >(tmp));
              ok=ok && tmp.size()==1;
              if(ok)
                {
                  if(point>=0)
                    ok=point==tmp[0];
                  else
                    point=tmp[0];
                }
            }
          w=std::find(w,conn+lgth,-1);
          w++;
        }
      if(ok && point>=0)
        {
          std::copy(quad4,quad4+4,retConn);
          retConn[4]=point;
          retLgth=5;
          return INTERP_KERNEL::NORM_PYRA5;
        }
    }
  retLgth=lgth;
  std::copy(conn,conn+lgth,retConn);
  return INTERP_KERNEL::NORM_POLYHED;
}

/*!
 * Cell with 'conn' connectivity has been detected as a good candidate. Full check of this. If yes NORM_TETRA4 is returned.
 * If fails a POLYHED is returned. 
 */
INTERP_KERNEL::NormalizedCellType CellSimplify::tryToUnPolyTetra4(const int *conn, int nbOfFaces, int lgth, int *retConn, int& retLgth)
{
  if(std::find_if(conn+lgth,conn+lgth+nbOfFaces,std::bind2nd(std::not_equal_to<int>(),(int)INTERP_KERNEL::NORM_TRI3))==conn+lgth+nbOfFaces)
    {
      std::set<int> tribase(conn,conn+3);
      int point=-1;
      bool ok=true;
      for(int i=1;i<4 && ok;i++)
        {
          std::vector<int> tmp;
          std::set<int> conn2(conn+i*4,conn+4*i+3);
          std::set_intersection(conn2.begin(),conn2.end(),tribase.begin(),tribase.end(),std::back_insert_iterator< std::vector<int> >(tmp));
          ok=tmp.size()==2;
          tmp.clear();
          std::set_difference(conn2.begin(),conn2.end(),tribase.begin(),tribase.end(),std::back_insert_iterator< std::vector<int> >(tmp));
          ok=ok && tmp.size()==1;
          if(ok)
            {
              if(point>=0)
                ok=point==tmp[0];
              else
                point=tmp[0];
            }
        }
      if(ok && point>=0)
        {
          std::copy(conn,conn+3,retConn);
          retConn[3]=point;
          retLgth=4;
          return INTERP_KERNEL::NORM_TETRA4;
        }
    }
  retLgth=lgth;
  std::copy(conn,conn+lgth,retConn);
  return INTERP_KERNEL::NORM_POLYHED;
}
