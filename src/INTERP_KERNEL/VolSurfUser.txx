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
      case INTERP_KERNEL::NORM_SEG3 :
        {
          int N1 = OTT<ConnType,numPol>::coo2C(connec[0]);
          int N2 = OTT<ConnType,numPol>::coo2C(connec[1]);
          return INTERP_KERNEL::calculateLgthForSeg2(coords+(SPACEDIM*N1),coords+(SPACEDIM*N2),SPACEDIM);
        }
      case INTERP_KERNEL::NORM_TRI3 :
      case INTERP_KERNEL::NORM_TRI6 :
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
            
      case INTERP_KERNEL::NORM_QUAD4 :
      case INTERP_KERNEL::NORM_QUAD8 :
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
            
      case INTERP_KERNEL::NORM_POLYGON :
        {          
          const double **pts=new const double *[lgth];
          for(int inod=0;inod<lgth;inod++)
            pts[inod] = coords+3*OTT<ConnType,numPol>::coo2C(connec[inod]);
          double val=INTERP_KERNEL::calculateAreaForPolyg(pts,lgth,SPACEDIM);
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

  template<class ConnType, NumberingPolicy numPolConn,int SPACEDIM>
  void computeBarycenter(NormalizedCellType type, const ConnType *connec, int lgth, const double *coords, double *res)
  {
    switch(type)
      {
      case NORM_TRI3:
      case NORM_QUAD4:
      case NORM_POLYGON:
        {
          if(SPACEDIM==2)
            computePolygonBarycenter2D<ConnType,numPolConn>(connec,lgth,coords,res);
          else if(SPACEDIM==3)
            computePolygonBarycenter3D<ConnType,numPolConn>(connec,lgth,coords,res);
          else
            throw INTERP_KERNEL::Exception("Impossible spacedim linked to cell 2D Cell !");
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
    //if(spaceDim==1)
    //  return computeBarycenter<ConnType,numPolConn,1>(type,connec,lgth,coords,res);
    throw INTERP_KERNEL::Exception("Invalid spaceDim specified for compute barycenter : must be 1, 2 or 3");
  }
}

#endif
