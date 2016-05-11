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

#include "VolSurfUser.hxx"
#include "InterpKernelAutoPtr.hxx"

#include <cmath>
#include <limits>
#include <algorithm>
#include <functional>

namespace INTERP_KERNEL
{
  /* Orthogonal distance from a point to a plane defined by three points p1, p2, p3.
   * Returns a signed distance, the normal of the plane being defined by (p1-p2)x(p3-p2)
   */
  double OrthoDistanceFromPtToPlaneInSpaceDim3(const double *p, const double *p1, const double *p2, const double *p3)
  {
    double prec = 1.0e-14;
    double T[2][3] = {{p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]},
                      {p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2]}};
    double N[3] = {T[0][1]*T[1][2]-T[0][2]*T[1][1],
                   T[0][2]*T[1][0]-T[0][0]*T[1][2],
                   T[0][0]*T[1][1]-T[0][1]*T[1][0]};

    double norm2 = N[0]*N[0] + N[1]*N[1] + N[2]*N[2];
    if (norm2 < prec)
      throw INTERP_KERNEL::Exception("OrthoDistanceFromPtToPlaneInSpaceDim3: degenerated normal vector!");
    double num = N[0]*(p[0]-p1[0]) + N[1]*(p[1]-p1[1]) + N[2]*(p[2]-p1[2]);
    return num/sqrt(norm2);
  }

  double SquareDistanceFromPtToSegInSpaceDim2(const double *pt, const double *pt0Seg2, const double *pt1Seg2, std::size_t &nbOfHint)
  {
    double dx=pt1Seg2[0]-pt0Seg2[0],dy=pt1Seg2[1]-pt0Seg2[1];
    double norm=sqrt(dx*dx+dy*dy);
    if(norm==0.)
      return (pt[0]-pt0Seg2[0])*(pt[0]-pt0Seg2[0])+(pt[1]-pt0Seg2[1])*(pt[1]-pt0Seg2[1]);//return std::numeric_limits<double>::max();
    dx/=norm; dy/=norm;
    double dx2=pt[0]-pt0Seg2[0],dy2=pt[1]-pt0Seg2[1];
    double dotP=(dx2*dx+dy2*dy);
    if(dotP<0. || dotP>norm)
      return dotP<0.?(pt[0]-pt0Seg2[0])*(pt[0]-pt0Seg2[0])+(pt[1]-pt0Seg2[1])*(pt[1]-pt0Seg2[1]):(pt[0]-pt1Seg2[0])*(pt[0]-pt1Seg2[0])+(pt[1]-pt1Seg2[1])*(pt[1]-pt1Seg2[1]);
    nbOfHint++;
    double x=pt0Seg2[0]+dotP*dx,y=pt0Seg2[1]+dotP*dy;
    return (x-pt[0])*(x-pt[0])+(y-pt[1])*(y-pt[1]);
  }

  double DistanceFromPtToTriInSpaceDim3(const double *pt, const double *pt0Tri3, const double *pt1Tri3, const double *pt2Tri3)
  {
    double matrix[12];
    if(!ComputeRotTranslationMatrixToPut3PointsOnOXY(pt0Tri3,pt1Tri3,pt2Tri3,matrix))
      return std::numeric_limits<double>::max();
    double xy0[2],xy1[2],xy2[2],xy[2]; xy0[0]=0.; xy0[1]=0.;
    xy1[0]=matrix[0]*pt1Tri3[0]+matrix[1]*pt1Tri3[1]+matrix[2]*pt1Tri3[2]+matrix[3]; xy1[1]=0.;
    xy2[0]=matrix[0]*pt2Tri3[0]+matrix[1]*pt2Tri3[1]+matrix[2]*pt2Tri3[2]+matrix[3];
    xy2[1]=matrix[4]*pt2Tri3[0]+matrix[5]*pt2Tri3[1]+matrix[6]*pt2Tri3[2]+matrix[7];
    xy[0]=matrix[0]*pt[0]+matrix[1]*pt[1]+matrix[2]*pt[2]+matrix[3];
    xy[1]=matrix[4]*pt[0]+matrix[5]*pt[1]+matrix[6]*pt[2]+matrix[7];
    double z=matrix[8]*pt[0]+matrix[9]*pt[1]+matrix[10]*pt[2]+matrix[11];
    double ret=std::numeric_limits<double>::max();
    std::size_t nbOfHint=0;
    if(xy[0]>0. && xy[0]<xy1[0])
      { ret=std::min(ret,z*z+xy[1]*xy[1]); nbOfHint++; } //distance pt to edge [pt0Tri3,pt1Tri3]
    double tmp=SquareDistanceFromPtToSegInSpaceDim2(xy,xy1,xy2,nbOfHint); //distance pt to edge [pt1Tri3,pt2Tri3]
    ret=std::min(ret,z*z+tmp);
    tmp=SquareDistanceFromPtToSegInSpaceDim2(xy,xy2,xy0,nbOfHint);//distance pt to edge [pt2Tri3,pt0Tri3]
    ret=std::min(ret,z*z+tmp);
    if(nbOfHint==3)
      ret=std::min(ret,z*z);
  return sqrt(ret);
  }

  double DistanceFromPtToPolygonInSpaceDim3(const double *pt, const int *connOfPolygonBg, const int *connOfPolygonEnd, const double *coords)
  {
    std::size_t nbOfEdges=std::distance(connOfPolygonBg,connOfPolygonEnd);
    if(nbOfEdges<3)
      throw INTERP_KERNEL::Exception("DistanceFromPtToPolygonInSpaceDim3 : trying to compute a distance to a polygon containing less than 3 edges !");
    double baryOfNodes[3]={0.,0.,0.};
    for(std::size_t i=0;i<nbOfEdges;i++)
      { baryOfNodes[0]+=coords[3*connOfPolygonBg[i]]; baryOfNodes[1]+=coords[3*connOfPolygonBg[i]+1]; baryOfNodes[2]+=coords[3*connOfPolygonBg[i]+2]; }
    std::transform(baryOfNodes,baryOfNodes+3,baryOfNodes,std::bind2nd(std::multiplies<double>(),1./((double)nbOfEdges)));
    double matrix[12];
    if(!ComputeRotTranslationMatrixToPut3PointsOnOXY(coords+3*connOfPolygonBg[0],coords+3*connOfPolygonBg[1],baryOfNodes,matrix))
      return std::numeric_limits<double>::max();
    INTERP_KERNEL::AutoPtr<double> ptXY=new double[2*nbOfEdges]; ptXY[0]=0.; ptXY[1]=0.;
    ptXY[2]=matrix[0]*coords[3*connOfPolygonBg[1]]+matrix[1]*coords[3*connOfPolygonBg[1]+1]+matrix[2]*coords[3*connOfPolygonBg[1]+2]+matrix[3]; ptXY[3]=0.;
    for(std::size_t i=2;i<nbOfEdges;i++)
      {
        ptXY[2*i]=matrix[0]*coords[3*connOfPolygonBg[i]]+matrix[1]*coords[3*connOfPolygonBg[i]+1]+matrix[2]*coords[3*connOfPolygonBg[i]+2]+matrix[3];
        ptXY[2*i+1]=matrix[4]*coords[3*connOfPolygonBg[i]]+matrix[5]*coords[3*connOfPolygonBg[i]+1]+matrix[6]*coords[3*connOfPolygonBg[i]+2]+matrix[7];
      }
    double xy[2]={matrix[0]*pt[0]+matrix[1]*pt[1]+matrix[2]*pt[2]+matrix[3],matrix[4]*pt[0]+matrix[5]*pt[1]+matrix[6]*pt[2]+matrix[7]};
    double z=matrix[8]*pt[0]+matrix[9]*pt[1]+matrix[10]*pt[2]+matrix[11];
    double ret=std::numeric_limits<double>::max();
    std::size_t nbOfHint=0;
    for(std::size_t i=0;i<nbOfEdges;i++)
      {
        double tmp=SquareDistanceFromPtToSegInSpaceDim2(xy,((double *)ptXY)+2*i,((double *)ptXY)+2*((i+1)%nbOfEdges),nbOfHint);
        ret=std::min(ret,z*z+tmp);
      }
    if(nbOfHint==nbOfEdges)
      ret=std::min(ret,z*z);
    return sqrt(ret);
  }

  /*!
   * \param [out] matrix contain a dense matrix of size 12 with 3 rows containing each 4 colums. This matrix is the reduction of 4x4 matrix but the last
   *              line containing [0,0,0,1] is omitted.
   */
  bool ComputeRotTranslationMatrixToPut3PointsOnOXY(const double *p0, const double *p1, const double *p2, double *matrix)
  {
    double norm=sqrt((p1[0]-p0[0])*(p1[0]-p0[0])+(p1[1]-p0[1])*(p1[1]-p0[1])+(p1[2]-p0[2])*(p1[2]-p0[2]));
    double c=(p1[0]-p0[0])/norm;
    double s=sqrt(1-c*c);
    double y=p1[2]-p0[2],z=p0[1]-p1[1];
    norm=sqrt(y*y+z*z);
    if(norm!=0.)
      { y/=norm; z/=norm; }
    double r0[9]={c,-z*s,y*s,
                  z*s,y*y*(1-c)+c,y*z*(1-c),
                  -y*s,z*y*(1-c),z*z*(1-c)+c};
    // 2nd rotation matrix
    double x=p2[0]-p0[0];
    y=p2[1]-p0[1]; z=p2[2]-p0[2];
    double y1=x*r0[3]+y*r0[4]+z*r0[5],z1=x*r0[6]+y*r0[7]+z*r0[8];
    c=y1/sqrt(y1*y1+z1*z1);
    s=sqrt(1.-c*c);
    //
    std::copy(r0,r0+3,matrix);
    matrix[4]=c*r0[3]-s*r0[6]; matrix[5]=c*r0[4]-s*r0[7]; matrix[6]=c*r0[5]-s*r0[8];
    matrix[8]=s*r0[3]+c*r0[6]; matrix[9]=s*r0[4]+c*r0[7]; matrix[10]=s*r0[5]+c*r0[8];
    matrix[3]=-p0[0]*matrix[0]-p0[1]*matrix[1]-p0[2]*matrix[2];
    matrix[7]=-p0[0]*matrix[4]-p0[1]*matrix[5]-p0[2]*matrix[6];
    matrix[11]=-p0[0]*matrix[8]-p0[1]*matrix[9]-p0[2]*matrix[10];
    return true;
  }
}

