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

#include "SplitterTetra.hxx"

namespace INTERP_KERNEL
{

  void SplitHexa8IntoTetras(SplittingPolicy policy, const int *nodalConnBg, const int *nodalConnEnd, const double *coords,
                            std::vector<int>& tetrasNodalConn, std::vector<double>& addCoords)
  {
    if(std::distance(nodalConnBg,nodalConnEnd)!=8)
      throw INTERP_KERNEL::Exception("SplitHexa8IntoTetras : input hexa do not have 8 nodes !");
    switch(policy)
      {
      case PLANAR_FACE_5:
        {
          tetrasNodalConn.resize(20);
          int *conn(&tetrasNodalConn[0]);
          conn[0]=nodalConnBg[SPLIT_NODES_5_WO[0]]; conn[1]=nodalConnBg[SPLIT_NODES_5_WO[1]]; conn[2]=nodalConnBg[SPLIT_NODES_5_WO[2]]; conn[3]=nodalConnBg[SPLIT_NODES_5_WO[3]];
          conn[4]=nodalConnBg[SPLIT_NODES_5_WO[4]]; conn[5]=nodalConnBg[SPLIT_NODES_5_WO[5]]; conn[6]=nodalConnBg[SPLIT_NODES_5_WO[6]]; conn[7]=nodalConnBg[SPLIT_NODES_5_WO[7]];
          conn[8]=nodalConnBg[SPLIT_NODES_5_WO[8]]; conn[9]=nodalConnBg[SPLIT_NODES_5_WO[9]]; conn[10]=nodalConnBg[SPLIT_NODES_5_WO[10]]; conn[11]=nodalConnBg[SPLIT_NODES_5_WO[11]];
          conn[12]=nodalConnBg[SPLIT_NODES_5_WO[12]]; conn[13]=nodalConnBg[SPLIT_NODES_5_WO[13]]; conn[14]=nodalConnBg[SPLIT_NODES_5_WO[14]]; conn[15]=nodalConnBg[SPLIT_NODES_5_WO[15]];
          conn[16]=nodalConnBg[SPLIT_NODES_5_WO[16]]; conn[17]=nodalConnBg[SPLIT_NODES_5_WO[17]]; conn[18]=nodalConnBg[SPLIT_NODES_5_WO[18]]; conn[19]=nodalConnBg[SPLIT_NODES_5_WO[19]];
          return ;
        }
      case PLANAR_FACE_6:
        {
          tetrasNodalConn.resize(24);
          int *conn(&tetrasNodalConn[0]);
          conn[0]=nodalConnBg[SPLIT_NODES_6_WO[0]]; conn[1]=nodalConnBg[SPLIT_NODES_6_WO[1]]; conn[2]=nodalConnBg[SPLIT_NODES_6_WO[2]]; conn[3]=nodalConnBg[SPLIT_NODES_6_WO[3]];
          conn[4]=nodalConnBg[SPLIT_NODES_6_WO[4]]; conn[5]=nodalConnBg[SPLIT_NODES_6_WO[5]]; conn[6]=nodalConnBg[SPLIT_NODES_6_WO[6]]; conn[7]=nodalConnBg[SPLIT_NODES_6_WO[7]];
          conn[8]=nodalConnBg[SPLIT_NODES_6_WO[8]]; conn[9]=nodalConnBg[SPLIT_NODES_6_WO[9]]; conn[10]=nodalConnBg[SPLIT_NODES_6_WO[10]]; conn[11]=nodalConnBg[SPLIT_NODES_6_WO[11]];
          conn[12]=nodalConnBg[SPLIT_NODES_6_WO[12]]; conn[13]=nodalConnBg[SPLIT_NODES_6_WO[13]]; conn[14]=nodalConnBg[SPLIT_NODES_6_WO[14]]; conn[15]=nodalConnBg[SPLIT_NODES_6_WO[15]];
          conn[16]=nodalConnBg[SPLIT_NODES_6_WO[16]]; conn[17]=nodalConnBg[SPLIT_NODES_6_WO[17]]; conn[18]=nodalConnBg[SPLIT_NODES_6_WO[18]]; conn[19]=nodalConnBg[SPLIT_NODES_6_WO[19]];
          conn[20]=nodalConnBg[SPLIT_NODES_6_WO[20]]; conn[21]=nodalConnBg[SPLIT_NODES_6_WO[21]]; conn[22]=nodalConnBg[SPLIT_NODES_6_WO[22]]; conn[23]=nodalConnBg[SPLIT_NODES_6_WO[23]];
          return ;
        }
      case GENERAL_24:
        {
          addCoords.resize(7*3);
          tetrasNodalConn.resize(24*4);
          int *conn(&tetrasNodalConn[0]);
          double *tmp(&addCoords[18]);
          tmp[0]=0.; tmp[1]=0.; tmp[2]=0.;
          double *tmp2(&addCoords[0]);
          for(int i=0;i<6;i++,tmp2+=3)
            {
              tmp2[0]=0.; tmp2[1]=0.; tmp2[2]=0.;
              for(int j=0;j<4;j++,conn+=4)
                {
                  int tmp3(nodalConnBg[GENERAL_24_SUB_NODES_WO[4*i+j]]);
                  tmp2[0]+=coords[3*tmp3+0];
                  tmp2[1]+=coords[3*tmp3+1];
                  tmp2[2]+=coords[3*tmp3+2];
                  conn[0]=tmp3;
                  conn[1]=nodalConnBg[GENERAL_24_SUB_NODES_WO[4*i+(j+1)%4]];
                  conn[2]=-(i+1); conn[3]=-7;
                }
              tmp2[0]/=4.; tmp2[1]/=4.; tmp2[2]/=4.;
              tmp[0]+=tmp2[0]; tmp[1]+=tmp2[1]; tmp[2]+=tmp2[2];
            }
          tmp[0]/=6.; tmp[1]/=6.; tmp[2]/=6.;
          return ;
        }
      case GENERAL_48:
        {
          addCoords.resize(19*3);
          tetrasNodalConn.resize(48*4);
          double *tmp2(&addCoords[0]),*tmp(&addCoords[0]);
          for(int i=0;i<12;i++,tmp2+=3)
            {
              tmp2[0]=(coords[3*nodalConnBg[GENERAL_48_SUB_NODES[2*i]]+0]+coords[3*nodalConnBg[GENERAL_48_SUB_NODES[2*i+1]]+0])/2.;
              tmp2[1]=(coords[3*nodalConnBg[GENERAL_48_SUB_NODES[2*i]]+1]+coords[3*nodalConnBg[GENERAL_48_SUB_NODES[2*i+1]]+1])/2.;
              tmp2[2]=(coords[3*nodalConnBg[GENERAL_48_SUB_NODES[2*i]]+2]+coords[3*nodalConnBg[GENERAL_48_SUB_NODES[2*i+1]]+2])/2.;
            }
          for(int i=0;i<7;i++,tmp2+=3)
            {
              tmp2[0]=(tmp[3*(GENERAL_48_SUB_NODES[2*i+24]-8)+0]+tmp[3*(GENERAL_48_SUB_NODES[2*i+25]-8)+0])/2.;
              tmp2[1]=(tmp[3*(GENERAL_48_SUB_NODES[2*i+24]-8)+1]+tmp[3*(GENERAL_48_SUB_NODES[2*i+25]-8)+1])/2.;
              tmp2[2]=(tmp[3*(GENERAL_48_SUB_NODES[2*i+24]-8)+2]+tmp[3*(GENERAL_48_SUB_NODES[2*i+25]-8)+2])/2.;
            }
          int *conn(&tetrasNodalConn[0]);
          std::vector<double> dummy;
          for(int i=0;i<8;i++)
            {
              std::vector<int> c;
              SplitHexa8IntoTetras(PLANAR_FACE_6,GENERAL_48_SUBZONES_2+i*8,GENERAL_48_SUBZONES_2+(i+1)*8,coords,c,dummy);
              int *conn2(&c[0]);
              for(int j=0;j<6;j++,conn+=4,conn2+=4)
                {
                  conn[0]=conn2[0]>=0?nodalConnBg[conn2[0]]:conn2[0];
                  conn[1]=conn2[1]>=0?nodalConnBg[conn2[1]]:conn2[1];
                  conn[2]=conn2[2]>=0?nodalConnBg[conn2[2]]:conn2[2];
                  conn[3]=conn2[3]>=0?nodalConnBg[conn2[3]]:conn2[3];
                }
            }
          return ;
        }
      default:
        throw INTERP_KERNEL::Exception("SplitHexa8IntoTetras : invalid input policy ! Should be in [PLANAR_FACE_5,PLANAR_FACE_6,GENERAL_24,GENERAL_48] !");
      }
  }

  void SplitIntoTetras(SplittingPolicy policy, NormalizedCellType gt, const int *nodalConnBg, const int *nodalConnEnd, const double *coords,
                       std::vector<int>& tetrasNodalConn, std::vector<double>& addCoords)
  {
    switch(gt)
      {
      case NORM_TETRA4:
        {
          std::size_t sz(std::distance(nodalConnBg,nodalConnEnd));
          if(sz!=4)
            throw INTERP_KERNEL::Exception("SplitIntoTetras : input tetra do not have 4 nodes !");
          tetrasNodalConn.insert(tetrasNodalConn.end(),nodalConnBg,nodalConnEnd);
          return ;
        }
      case NORM_HEXA8:
        {
          SplitHexa8IntoTetras(policy,nodalConnBg,nodalConnEnd,coords,tetrasNodalConn,addCoords);
          return ;
        }
      case NORM_PYRA5:
        {
          std::size_t sz(std::distance(nodalConnBg,nodalConnEnd));
          if(sz!=5)
            throw INTERP_KERNEL::Exception("SplitIntoTetras : input pyra5 do not have 5 nodes !");
          tetrasNodalConn.resize(8);
          int *conn(&tetrasNodalConn[0]);
          conn[0]=nodalConnBg[0]; conn[1]=nodalConnBg[1]; conn[2]=nodalConnBg[2]; conn[3]=nodalConnBg[4];
          conn[4]=nodalConnBg[0]; conn[5]=nodalConnBg[2]; conn[6]=nodalConnBg[3]; conn[7]=nodalConnBg[4];
          return ;
        }
      case NORM_PENTA6:
        {
          std::size_t sz(std::distance(nodalConnBg,nodalConnEnd));
          if(sz!=6)
            throw INTERP_KERNEL::Exception("SplitIntoTetras : input penta6 do not have 6 nodes !");
          tetrasNodalConn.resize(12);
          int *conn(&tetrasNodalConn[0]);
          conn[0]=nodalConnBg[0]; conn[1]=nodalConnBg[1]; conn[2]=nodalConnBg[2]; conn[3]=nodalConnBg[3];
          conn[4]=nodalConnBg[3]; conn[5]=nodalConnBg[5]; conn[6]=nodalConnBg[4]; conn[7]=nodalConnBg[2];
          conn[8]=nodalConnBg[4]; conn[9]=nodalConnBg[2]; conn[10]=nodalConnBg[1]; conn[11]=nodalConnBg[3];
          return ;
        }
      case NORM_HEXGP12:
        {
          std::size_t sz(std::distance(nodalConnBg,nodalConnEnd));
          if(sz!=12)
            throw INTERP_KERNEL::Exception("SplitIntoTetras : input octa12 (hexagone prism) do not have 12 nodes !");
          tetrasNodalConn.resize(48);
          int *conn(&tetrasNodalConn[0]);
          conn[0]=nodalConnBg[0]; conn[1]=nodalConnBg[1]; conn[2]=nodalConnBg[5]; conn[3]=nodalConnBg[6];
          conn[4]=nodalConnBg[6]; conn[5]=nodalConnBg[11]; conn[6]=nodalConnBg[7]; conn[7]=nodalConnBg[5];
          conn[8]=nodalConnBg[7]; conn[9]=nodalConnBg[5]; conn[10]=nodalConnBg[1]; conn[11]=nodalConnBg[6];
          //
          conn[12]=nodalConnBg[1]; conn[13]=nodalConnBg[4]; conn[14]=nodalConnBg[5]; conn[15]=nodalConnBg[7];
          conn[16]=nodalConnBg[7]; conn[17]=nodalConnBg[11]; conn[18]=nodalConnBg[10]; conn[19]=nodalConnBg[5];
          conn[20]=nodalConnBg[10]; conn[21]=nodalConnBg[5]; conn[22]=nodalConnBg[4]; conn[23]=nodalConnBg[7];
          //
          conn[24]=nodalConnBg[1]; conn[25]=nodalConnBg[2]; conn[26]=nodalConnBg[4]; conn[27]=nodalConnBg[7];
          conn[28]=nodalConnBg[7]; conn[29]=nodalConnBg[10]; conn[30]=nodalConnBg[8]; conn[31]=nodalConnBg[4];
          conn[32]=nodalConnBg[8]; conn[33]=nodalConnBg[4]; conn[34]=nodalConnBg[2]; conn[35]=nodalConnBg[7];
          //
          conn[36]=nodalConnBg[2]; conn[37]=nodalConnBg[3]; conn[38]=nodalConnBg[4]; conn[39]=nodalConnBg[8];
          conn[40]=nodalConnBg[8]; conn[41]=nodalConnBg[10]; conn[42]=nodalConnBg[9]; conn[43]=nodalConnBg[4];
          conn[44]=nodalConnBg[9]; conn[45]=nodalConnBg[4]; conn[46]=nodalConnBg[3]; conn[47]=nodalConnBg[8];
          return ;
        }
      case NORM_POLYHED:
        {
          std::size_t nbOfFaces(std::count(nodalConnBg,nodalConnEnd,-1)+1);
          std::size_t nbOfTetra(std::distance(nodalConnBg,nodalConnEnd)-nbOfFaces+1);
          addCoords.resize((nbOfFaces+1)*3);
          tetrasNodalConn.resize(nbOfTetra*4);
          int *conn(&tetrasNodalConn[0]);
          const int *work(nodalConnBg);
          double *tmp(&addCoords[0]),*tmp2(&addCoords[3*nbOfFaces]);
          tmp2[0]=0.; tmp2[1]=0.; tmp2[2]=0.;
          for(std::size_t i=0;i<nbOfFaces;i++,tmp+=3)
            {
              tmp[0]=0.; tmp[1]=0.; tmp[2]=0.;
              std::size_t nbOfNodesOfFace(std::distance(work,std::find(work,nodalConnEnd,-1)));
              for(std::size_t j=0;j<nbOfNodesOfFace;j++,conn+=4)
                {
                  conn[0]=work[j]; conn[1]=work[(j+1)%nbOfNodesOfFace]; conn[2]=-((int)i+1); conn[3]=-((int)nbOfFaces+1);
                  tmp[0]+=coords[3*work[j]+0]; tmp[1]+=coords[3*work[j]+1]; tmp[2]+=coords[3*work[j]+2];
                }
              tmp[0]/=(int)nbOfNodesOfFace; tmp[1]/=(int)nbOfNodesOfFace; tmp[2]/=(int)nbOfNodesOfFace;
              tmp2[0]+=tmp[0]; tmp2[1]+=tmp[1]; tmp2[2]+=tmp[2];
              work+=nbOfNodesOfFace+1;
            }
          tmp2[0]/=(int)nbOfFaces; tmp2[1]/=(int)nbOfFaces; tmp2[2]/=(int)nbOfFaces;
          return ;
        }
      default:
        throw INTERP_KERNEL::Exception("SplitIntoTetras : not managed such Geometric type ! Available geometric types are all 3D linear cells !");
      }
  }
}
