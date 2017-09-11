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
#ifndef __VOLSURFUSER_TXX__
#define __VOLSURFUSER_TXX__

#include "VolSurfUser.hxx"
#include "VolSurfFormulae.hxx"
#include "InterpolationUtils.hxx"

#include <algorithm>

namespace INTERP_KERNEL
{
  template<class ConnType, NumberingPolicy numPol, int SPACEDIM>
  double computeVolSurfOfCell(NormalizedCellType type, const ConnType *connec, int lgth, const double *coords)
  {
    switch(type)
      {
      case INTERP_KERNEL::NORM_SEG2 :
      case INTERP_KERNEL::NORM_SEG4 :
        {
          int N1 = OTT<ConnType,numPol>::coo2C(connec[0]);
          int N2 = OTT<ConnType,numPol>::coo2C(connec[1]);
          return INTERP_KERNEL::calculateLgthForSeg2(coords+(SPACEDIM*N1),coords+(SPACEDIM*N2),SPACEDIM);
        }
      case INTERP_KERNEL::NORM_SEG3 :
        {
          int beginNode = OTT<ConnType,numPol>::coo2C(connec[0]);
          int endNode = OTT<ConnType,numPol>::coo2C(connec[1]);
          int middleNode = OTT<ConnType,numPol>::coo2C(connec[2]);
          return INTERP_KERNEL::calculateLgthForSeg3(coords+(SPACEDIM*beginNode),coords+(SPACEDIM*endNode),coords+(SPACEDIM*middleNode),SPACEDIM);
        }
      case INTERP_KERNEL::NORM_TRI3 :
        {
          int N1 = OTT<ConnType,numPol>::coo2C(connec[0]);
          int N2 = OTT<ConnType,numPol>::coo2C(connec[1]);
          int N3 = OTT<ConnType,numPol>::coo2C(connec[2]);
              
          return INTERP_KERNEL::calculateAreaForTria(coords+(SPACEDIM*N1),
                                                     coords+(SPACEDIM*N2),
                                                     coords+(SPACEDIM*N3),
                                                     SPACEDIM);
        }
        break;
            
      case INTERP_KERNEL::NORM_TRI6 :
      case INTERP_KERNEL::NORM_TRI7 :
        {
          const double *pts[6];
          pts[0] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[0]);
          pts[1] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[1]);
          pts[2] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[2]);
          pts[3] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[3]);
          pts[4] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[4]);
          pts[5] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[5]);
          return INTERP_KERNEL::calculateAreaForQPolyg(pts,6,SPACEDIM);
        }
        break;
      case INTERP_KERNEL::NORM_QUAD4 :
        {
          int N1 = OTT<ConnType,numPol>::coo2C(connec[0]);
          int N2 = OTT<ConnType,numPol>::coo2C(connec[1]);
          int N3 = OTT<ConnType,numPol>::coo2C(connec[2]);
          int N4 = OTT<ConnType,numPol>::coo2C(connec[3]);
              
          return INTERP_KERNEL::calculateAreaForQuad(coords+SPACEDIM*N1,
                                                     coords+SPACEDIM*N2,
                                                     coords+SPACEDIM*N3,
                                                     coords+SPACEDIM*N4,
                                                     SPACEDIM);
        }
        break;
      case INTERP_KERNEL::NORM_QUAD8 :
      case INTERP_KERNEL::NORM_QUAD9 :  
        {
          const double *pts[8];
          pts[0] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[0]);
          pts[1] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[1]);
          pts[2] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[2]);
          pts[3] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[3]);
          pts[4] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[4]);
          pts[5] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[5]);
          pts[6] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[6]);
          pts[7] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[7]);
          return INTERP_KERNEL::calculateAreaForQPolyg(pts,8,SPACEDIM);
        }
        break;
      case INTERP_KERNEL::NORM_POLYGON :
        {          
          const double **pts=new const double *[lgth];
          for(int inod=0;inod<lgth;inod++)
            pts[inod] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[inod]);
          double val=INTERP_KERNEL::calculateAreaForPolyg(pts,lgth,SPACEDIM);
          delete [] pts;
          return val;
        }
        break;
      case INTERP_KERNEL::NORM_QPOLYG :
        {
          const double **pts=new const double *[lgth];
          for(int inod=0;inod<lgth;inod++)
            pts[inod] = coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[inod]);
          double val=INTERP_KERNEL::calculateAreaForQPolyg(pts,lgth,SPACEDIM);
          delete [] pts;
          return val;
        }
        break;
      case INTERP_KERNEL::NORM_TETRA4 :
      case INTERP_KERNEL::NORM_TETRA10 :
        {
          int N1 = OTT<ConnType,numPol>::coo2C(connec[0]);
          int N2 = OTT<ConnType,numPol>::coo2C(connec[1]);
          int N3 = OTT<ConnType,numPol>::coo2C(connec[2]);
          int N4 = OTT<ConnType,numPol>::coo2C(connec[3]);
              
          return INTERP_KERNEL::calculateVolumeForTetra(coords+SPACEDIM*N1,
                                                        coords+SPACEDIM*N2,
                                                        coords+SPACEDIM*N3,
                                                        coords+SPACEDIM*N4);
        }
        break;
            
      case INTERP_KERNEL::NORM_PYRA5 :
      case INTERP_KERNEL::NORM_PYRA13 :
        {
          int N1 = OTT<ConnType,numPol>::coo2C(connec[0]);
          int N2 = OTT<ConnType,numPol>::coo2C(connec[1]);
          int N3 = OTT<ConnType,numPol>::coo2C(connec[2]);
          int N4 = OTT<ConnType,numPol>::coo2C(connec[3]);
          int N5 = OTT<ConnType,numPol>::coo2C(connec[4]);
              
          return INTERP_KERNEL::calculateVolumeForPyra(coords+SPACEDIM*N1,
                                                       coords+SPACEDIM*N2,
                                                       coords+SPACEDIM*N3,
                                                       coords+SPACEDIM*N4,
                                                       coords+SPACEDIM*N5);
        }
        break;
            
      case INTERP_KERNEL::NORM_PENTA6 :
      case INTERP_KERNEL::NORM_PENTA15 :
      case INTERP_KERNEL::NORM_PENTA18 :
        {
          int N1 = OTT<ConnType,numPol>::coo2C(connec[0]);
          int N2 = OTT<ConnType,numPol>::coo2C(connec[1]);
          int N3 = OTT<ConnType,numPol>::coo2C(connec[2]);
          int N4 = OTT<ConnType,numPol>::coo2C(connec[3]);
          int N5 = OTT<ConnType,numPol>::coo2C(connec[4]);
          int N6 = OTT<ConnType,numPol>::coo2C(connec[5]);
              
          return INTERP_KERNEL::calculateVolumeForPenta(coords+SPACEDIM*N1,
                                                        coords+SPACEDIM*N2,
                                                        coords+SPACEDIM*N3,
                                                        coords+SPACEDIM*N4,
                                                        coords+SPACEDIM*N5,
                                                        coords+SPACEDIM*N6);
        }
        break;
            
      case INTERP_KERNEL::NORM_HEXA8 :
      case INTERP_KERNEL::NORM_HEXA20 :
      case INTERP_KERNEL::NORM_HEXA27 :
        {
          int N1 = OTT<ConnType,numPol>::coo2C(connec[0]);
          int N2 = OTT<ConnType,numPol>::coo2C(connec[1]);
          int N3 = OTT<ConnType,numPol>::coo2C(connec[2]);
          int N4 = OTT<ConnType,numPol>::coo2C(connec[3]);
          int N5 = OTT<ConnType,numPol>::coo2C(connec[4]);
          int N6 = OTT<ConnType,numPol>::coo2C(connec[5]);
          int N7 = OTT<ConnType,numPol>::coo2C(connec[6]);
          int N8 = OTT<ConnType,numPol>::coo2C(connec[7]);
              
          return INTERP_KERNEL::calculateVolumeForHexa(coords+SPACEDIM*N1,
                                                       coords+SPACEDIM*N2,
                                                       coords+SPACEDIM*N3,
                                                       coords+SPACEDIM*N4,
                                                       coords+SPACEDIM*N5,
                                                       coords+SPACEDIM*N6,
                                                       coords+SPACEDIM*N7,
                                                       coords+SPACEDIM*N8);
        }
        break;
      case INTERP_KERNEL::NORM_HEXGP12:
        {
          const int connecHexa12[43]={
            OTT<ConnType,numPol>::coo2C(connec[0]),OTT<ConnType,numPol>::coo2C(connec[1]),OTT<ConnType,numPol>::coo2C(connec[2]),OTT<ConnType,numPol>::coo2C(connec[3]),OTT<ConnType,numPol>::coo2C(connec[4]),OTT<ConnType,numPol>::coo2C(connec[5]),-1,
            OTT<ConnType,numPol>::coo2C(connec[6]),OTT<ConnType,numPol>::coo2C(connec[11]),OTT<ConnType,numPol>::coo2C(connec[10]),OTT<ConnType,numPol>::coo2C(connec[9]),OTT<ConnType,numPol>::coo2C(connec[8]),OTT<ConnType,numPol>::coo2C(connec[7]),-1,
            OTT<ConnType,numPol>::coo2C(connec[0]),OTT<ConnType,numPol>::coo2C(connec[6]),OTT<ConnType,numPol>::coo2C(connec[7]),OTT<ConnType,numPol>::coo2C(connec[1]),-1,
            OTT<ConnType,numPol>::coo2C(connec[1]),OTT<ConnType,numPol>::coo2C(connec[7]),OTT<ConnType,numPol>::coo2C(connec[8]),OTT<ConnType,numPol>::coo2C(connec[2]),-1,
            OTT<ConnType,numPol>::coo2C(connec[2]),OTT<ConnType,numPol>::coo2C(connec[8]),OTT<ConnType,numPol>::coo2C(connec[9]),OTT<ConnType,numPol>::coo2C(connec[3]),-1,
            OTT<ConnType,numPol>::coo2C(connec[3]),OTT<ConnType,numPol>::coo2C(connec[9]),OTT<ConnType,numPol>::coo2C(connec[10]),OTT<ConnType,numPol>::coo2C(connec[4]),-1,
            OTT<ConnType,numPol>::coo2C(connec[4]),OTT<ConnType,numPol>::coo2C(connec[10]),OTT<ConnType,numPol>::coo2C(connec[11]),OTT<ConnType,numPol>::coo2C(connec[5]),-1,
            OTT<ConnType,numPol>::coo2C(connec[5]),OTT<ConnType,numPol>::coo2C(connec[11]),OTT<ConnType,numPol>::coo2C(connec[6]),OTT<ConnType,numPol>::coo2C(connec[0])};
          return calculateVolumeForPolyh2<ConnType,numPol>(connecHexa12,43,coords);
        }
      case INTERP_KERNEL::NORM_POLYHED :
        {
          return calculateVolumeForPolyh2<ConnType,numPol>(connec,lgth,coords);
        }
        break;
      default:
        throw INTERP_KERNEL::Exception("Not recognized cell type to get Length/Area/Volume on it !");
      }
  }

  template<class ConnType, NumberingPolicy numPolConn>
  double computeVolSurfOfCell2(NormalizedCellType type, const ConnType *connec, int lgth, const double *coords, int spaceDim)
  {
    if(spaceDim==3)
      return computeVolSurfOfCell<ConnType,numPolConn,3>(type,connec,lgth,coords);
    if(spaceDim==2)
      return computeVolSurfOfCell<ConnType,numPolConn,2>(type,connec,lgth,coords);
    if(spaceDim==1)
      return computeVolSurfOfCell<ConnType,numPolConn,1>(type,connec,lgth,coords);
    throw INTERP_KERNEL::Exception("Invalid spaceDim specified : must be 1, 2 or 3");
  }


  template<class ConnType, NumberingPolicy numPol,int SPACEDIM>
  void computeBarycenter(NormalizedCellType type, const ConnType *connec, int lgth, const double *coords, double *res)
  {
    switch(type)
      {
      case NORM_SEG2:
      case NORM_SEG4:
        {
          std::copy(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[0]),
                    coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[0]+1),res);
          std::transform(res,res+SPACEDIM,coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[1]),res,std::plus<double>());
          std::transform(res,res+SPACEDIM,res,std::bind2nd(std::multiplies<double>(),0.5));
          break;
        }
      case NORM_SEG3:
        {
          if(SPACEDIM==2)
            {
              Edge *ed=Edge::BuildEdgeFrom3Points(coords+2*OTT<ConnType,numPol>::coo2C(connec[0]),coords+2*OTT<ConnType,numPol>::coo2C(connec[2]),coords+2*OTT<ConnType,numPol>::coo2C(connec[1]));
              ed->getBarycenter(res);
              ed->decrRef();
            }
          else if(SPACEDIM==1)
            {
              *res=(coords[OTT<ConnType,numPol>::coo2C(connec[0])]+coords[OTT<ConnType,numPol>::coo2C(connec[1])])/2.;
            }
          else if(SPACEDIM==3)
            {
              std::copy(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[0]),
                        coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[0]+1),res);
              std::transform(res,res+SPACEDIM,coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[1]),res,std::plus<double>());
              std::transform(res,res+SPACEDIM,res,std::bind2nd(std::multiplies<double>(),0.5));
            }
          else
            throw INTERP_KERNEL::Exception("computeBarycenter for SEG3 only SPACEDIM 1,2 or 3 supported !");
          break;
        }
      case NORM_TRI3:
      case NORM_TRI7:
        {
          std::copy(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[0]),
                    coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[0]+1),res);
          std::transform(res,res+SPACEDIM,coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[1]),res,std::plus<double>());
          std::transform(res,res+SPACEDIM,coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[2]),res,std::plus<double>());
          std::transform(res,res+SPACEDIM,res,std::bind2nd(std::multiplies<double>(),1./3.));
          break;
        }
      case NORM_TRI6:
        {
          if(SPACEDIM==2)
            {
              double *pts[6];
              pts[0] = const_cast<double *>(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[0]));
              pts[1] = const_cast<double *>(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[1]));
              pts[2] = const_cast<double *>(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[2]));
              pts[3] = const_cast<double *>(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[3]));
              pts[4] = const_cast<double *>(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[4]));
              pts[5] = const_cast<double *>(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[5]));
              computeQPolygonBarycenter2D(pts,6,2,res);
            }
          else if(SPACEDIM==3)
            computePolygonBarycenter3D<ConnType,numPol>(connec,lgth/2,coords,res);
          else
            throw INTERP_KERNEL::Exception("Impossible spacedim linked to cell 2D Cell !");
          break;
        }
      case NORM_QUAD4:
      case NORM_POLYGON:
        {
          if(SPACEDIM==2)
            computePolygonBarycenter2D<ConnType,numPol>(connec,lgth,coords,res);
          else if(SPACEDIM==3)
            computePolygonBarycenter3D<ConnType,numPol>(connec,lgth,coords,res);
          else
            throw INTERP_KERNEL::Exception("Impossible spacedim linked to cell 2D Cell !");
          break;
        }
      case NORM_QUAD8:
        {
          if(SPACEDIM==2)
            {
              double *pts[8];
              pts[0] = const_cast<double *>(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[0]));
              pts[1] = const_cast<double *>(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[1]));
              pts[2] = const_cast<double *>(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[2]));
              pts[3] = const_cast<double *>(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[3]));
              pts[4] = const_cast<double *>(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[4]));
              pts[5] = const_cast<double *>(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[5]));
              pts[6] = const_cast<double *>(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[6]));
              pts[7] = const_cast<double *>(coords+SPACEDIM*OTT<ConnType,numPol>::coo2C(connec[7]));
              computeQPolygonBarycenter2D(pts,8,2,res);
            }
          else if(SPACEDIM==3)
            computePolygonBarycenter3D<ConnType,numPol>(connec,lgth/2,coords,res);
          else
            throw INTERP_KERNEL::Exception("Impossible spacedim linked to cell 2D Cell !");
          break;
        }
      case INTERP_KERNEL::NORM_QPOLYG :
        {
          if(SPACEDIM==2)
            {
              double **pts=new double *[lgth];
              for(int i=0;i<lgth;i++)
                pts[i]=const_cast<double *>(coords+2*OTT<ConnType,numPol>::coo2C(connec[i]));
              computeQPolygonBarycenter2D(pts,lgth,2,res);
              delete [] pts;
            }
          else if(SPACEDIM==3)
            computePolygonBarycenter3D<ConnType,numPol>(connec,lgth/2,coords,res);
          else
            throw INTERP_KERNEL::Exception("Impossible spacedim linked to cell 2D Cell !");
          break;
        }
        break;
      case NORM_TETRA4:
        {
          res[0]=coords[3*OTT<ConnType,numPol>::coo2C(connec[0])]; 
          res[1]=coords[3*OTT<ConnType,numPol>::coo2C(connec[0])+1];
          res[2]=coords[3*OTT<ConnType,numPol>::coo2C(connec[0])+2];
          res[0]+=coords[3*OTT<ConnType,numPol>::coo2C(connec[1])]; 
          res[1]+=coords[3*OTT<ConnType,numPol>::coo2C(connec[1])+1];
          res[2]+=coords[3*OTT<ConnType,numPol>::coo2C(connec[1])+2];
          res[0]+=coords[3*OTT<ConnType,numPol>::coo2C(connec[2])]; 
          res[1]+=coords[3*OTT<ConnType,numPol>::coo2C(connec[2])+1];
          res[2]+=coords[3*OTT<ConnType,numPol>::coo2C(connec[2])+2];
          res[0]+=coords[3*OTT<ConnType,numPol>::coo2C(connec[3])]; 
          res[1]+=coords[3*OTT<ConnType,numPol>::coo2C(connec[3])+1];
          res[2]+=coords[3*OTT<ConnType,numPol>::coo2C(connec[3])+2];
          res[0]*=0.25; res[1]*=0.25; res[2]*=0.25;
          break;
        }
      case NORM_PYRA5:
        {
          double tmp[3];
          computePolygonBarycenter3D<ConnType,numPol>(connec,lgth-1,coords,tmp);
          res[0]=(coords[3*OTT<ConnType,numPol>::coo2C(connec[4])]+3.*tmp[0])/4.;
          res[1]=(coords[3*OTT<ConnType,numPol>::coo2C(connec[4])+1]+3.*tmp[1])/4.;
          res[2]=(coords[3*OTT<ConnType,numPol>::coo2C(connec[4])+2]+3.*tmp[2])/4.;
          break;
        }
      case NORM_HEXA8:
        {
          const int conn[29]={
            OTT<ConnType,numPol>::coo2C(connec[0]),OTT<ConnType,numPol>::coo2C(connec[1]),OTT<ConnType,numPol>::coo2C(connec[2]),OTT<ConnType,numPol>::coo2C(connec[3]),-1,
            OTT<ConnType,numPol>::coo2C(connec[4]),OTT<ConnType,numPol>::coo2C(connec[7]),OTT<ConnType,numPol>::coo2C(connec[6]),OTT<ConnType,numPol>::coo2C(connec[5]),-1,
            OTT<ConnType,numPol>::coo2C(connec[0]),OTT<ConnType,numPol>::coo2C(connec[3]),OTT<ConnType,numPol>::coo2C(connec[7]),OTT<ConnType,numPol>::coo2C(connec[4]),-1,
            OTT<ConnType,numPol>::coo2C(connec[3]),OTT<ConnType,numPol>::coo2C(connec[2]),OTT<ConnType,numPol>::coo2C(connec[6]),OTT<ConnType,numPol>::coo2C(connec[7]),-1,
            OTT<ConnType,numPol>::coo2C(connec[2]),OTT<ConnType,numPol>::coo2C(connec[1]),OTT<ConnType,numPol>::coo2C(connec[5]),OTT<ConnType,numPol>::coo2C(connec[6]),-1,
            OTT<ConnType,numPol>::coo2C(connec[0]),OTT<ConnType,numPol>::coo2C(connec[4]),OTT<ConnType,numPol>::coo2C(connec[5]),OTT<ConnType,numPol>::coo2C(connec[1]),
            };
          barycenterOfPolyhedron<ConnType,numPol>(conn,29,coords,res);
          break;
        }
      case NORM_PENTA6:
        {
          const int conn[22]={
            OTT<ConnType,numPol>::coo2C(connec[0]),OTT<ConnType,numPol>::coo2C(connec[1]),OTT<ConnType,numPol>::coo2C(connec[2]),-1,
            OTT<ConnType,numPol>::coo2C(connec[3]),OTT<ConnType,numPol>::coo2C(connec[5]),OTT<ConnType,numPol>::coo2C(connec[4]),-1,
            OTT<ConnType,numPol>::coo2C(connec[0]),OTT<ConnType,numPol>::coo2C(connec[2]),OTT<ConnType,numPol>::coo2C(connec[5]),OTT<ConnType,numPol>::coo2C(connec[3]),-1,
            OTT<ConnType,numPol>::coo2C(connec[2]),OTT<ConnType,numPol>::coo2C(connec[1]),OTT<ConnType,numPol>::coo2C(connec[4]),OTT<ConnType,numPol>::coo2C(connec[5]),-1,
            OTT<ConnType,numPol>::coo2C(connec[1]),OTT<ConnType,numPol>::coo2C(connec[0]),OTT<ConnType,numPol>::coo2C(connec[3]),OTT<ConnType,numPol>::coo2C(connec[4])
          };
          barycenterOfPolyhedron<ConnType,numPol>(conn,22,coords,res);
          break;
        }
      case INTERP_KERNEL::NORM_HEXGP12:
        {
          const int connecHexa12[43]={
            OTT<ConnType,numPol>::coo2C(connec[0]),OTT<ConnType,numPol>::coo2C(connec[1]),OTT<ConnType,numPol>::coo2C(connec[2]),OTT<ConnType,numPol>::coo2C(connec[3]),OTT<ConnType,numPol>::coo2C(connec[4]),OTT<ConnType,numPol>::coo2C(connec[5]),-1,
            OTT<ConnType,numPol>::coo2C(connec[6]),OTT<ConnType,numPol>::coo2C(connec[11]),OTT<ConnType,numPol>::coo2C(connec[10]),OTT<ConnType,numPol>::coo2C(connec[9]),OTT<ConnType,numPol>::coo2C(connec[8]),OTT<ConnType,numPol>::coo2C(connec[7]),-1,
            OTT<ConnType,numPol>::coo2C(connec[0]),OTT<ConnType,numPol>::coo2C(connec[6]),OTT<ConnType,numPol>::coo2C(connec[7]),OTT<ConnType,numPol>::coo2C(connec[1]),-1,
            OTT<ConnType,numPol>::coo2C(connec[1]),OTT<ConnType,numPol>::coo2C(connec[7]),OTT<ConnType,numPol>::coo2C(connec[8]),OTT<ConnType,numPol>::coo2C(connec[2]),-1,
            OTT<ConnType,numPol>::coo2C(connec[2]),OTT<ConnType,numPol>::coo2C(connec[8]),OTT<ConnType,numPol>::coo2C(connec[9]),OTT<ConnType,numPol>::coo2C(connec[3]),-1,
            OTT<ConnType,numPol>::coo2C(connec[3]),OTT<ConnType,numPol>::coo2C(connec[9]),OTT<ConnType,numPol>::coo2C(connec[10]),OTT<ConnType,numPol>::coo2C(connec[4]),-1,
            OTT<ConnType,numPol>::coo2C(connec[4]),OTT<ConnType,numPol>::coo2C(connec[10]),OTT<ConnType,numPol>::coo2C(connec[11]),OTT<ConnType,numPol>::coo2C(connec[5]),-1,
            OTT<ConnType,numPol>::coo2C(connec[5]),OTT<ConnType,numPol>::coo2C(connec[11]),OTT<ConnType,numPol>::coo2C(connec[6]),OTT<ConnType,numPol>::coo2C(connec[0])};
          barycenterOfPolyhedron<ConnType,numPol>(connecHexa12,43,coords,res);
          break;
        }
      case NORM_POLYHED:
        {
          barycenterOfPolyhedron<ConnType,numPol>(connec,lgth,coords,res);
          break;
        }
      default:
        throw INTERP_KERNEL::Exception("Not recognized cell type to get Barycenter on it !");
      }
  }

  template<class ConnType, NumberingPolicy numPolConn>
  void computeBarycenter2(NormalizedCellType type, const ConnType *connec, int lgth, const double *coords, int spaceDim, double *res)
  {
    if(spaceDim==3)
      return computeBarycenter<ConnType,numPolConn,3>(type,connec,lgth,coords,res);
    if(spaceDim==2)
      return computeBarycenter<ConnType,numPolConn,2>(type,connec,lgth,coords,res);
    if(spaceDim==1)
      return computeBarycenter<ConnType,numPolConn,1>(type,connec,lgth,coords,res);
    throw INTERP_KERNEL::Exception("Invalid spaceDim specified for compute barycenter : must be 1, 2 or 3");
  }
}

#endif
