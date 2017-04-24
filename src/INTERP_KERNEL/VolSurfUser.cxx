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
#include "InterpolationUtils.hxx"

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

  /**
   * See http://geomalgorithms.com/a02-_lines.html#Distance-to-Ray-or-Segment
   */
  double DistanceFromPtToSegInSpaceDim3(const double *pt, const double *pt0Seg2, const double *pt1Seg2)
  {
    double v[3], w[3];
    for(int i=0; i < 3; i++) {
        v[i]=pt1Seg2[i]-pt0Seg2[i];
        w[i] = pt[i] - pt0Seg2[i];
    }

    double c1 = dotprod<3>(w,v);
    if ( c1 <= 0 )
      return norm<3>(w);
    double c2 = dotprod<3>(v,v);
    if ( c2 <= c1 )
      {
        for(int i=0; i < 3; i++)
          w[i] = pt[i] - pt1Seg2[i];
        return norm<3>(w);
      }
    double b = c1 / c2;
    for(int i=0; i < 3; i++)
      w[i] = pt0Seg2[i] + b * v[i] - pt[i];
    return norm<3>(w);
  }

  /**
     Helper for DistanceFromPtToTriInSpaceDim3
   */
  inline double _HelperDistancePtToTri3D_1(const double aXX, const double bX, const double c)
  {
    if (bX >= 0)
      return c;
    if (-bX >= aXX)
      return aXX + 2*bX + c;
    return bX*(-bX / aXX) + c;
  }

  /**
    Helper for DistanceFromPtToTriInSpaceDim3
   */
  inline double _HelperDistancePtToTri3D_2(const double a01, const double aXX, const double aYY,
                                           const double bX, const double bY, const double c)
  {
    double tmp0 = a01 + bX, tmp1 = aXX + bY;
    if (tmp1 > tmp0) {
        double numer = tmp1 - tmp0, denom = aXX - 2*a01 + aYY;
        if (numer >= denom)
          return aXX + 2*bX + c;
        else {
            double s, t;
            s = numer / denom; t = 1 - s;
            return s*(aXX*s + a01*t + 2*bX) + t*(a01*s + aYY*t + 2*bY) + c;
        }
    }
    else
      {
        if (tmp1 <= 0)   return aYY + 2*bY + c;
        else {
            if (bY >= 0) return c;
            else         return bY*(-bY / aYY) + c;
        }
      }
  }

  /**
   * From https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
   */
  double DistanceFromPtToTriInSpaceDim3(const double *pt, const double *pt0Tri3, const double *pt1Tri3, const double *pt2Tri3)
  {
    double diff[3], edge0[3], edge1[3];
    for(int i=0; i < 3; i++) diff[i]=pt0Tri3[i]-pt[i];
    for(int i=0; i < 3; i++) edge0[i]=pt1Tri3[i]-pt0Tri3[i];
    for(int i=0; i < 3; i++) edge1[i]=pt2Tri3[i]-pt0Tri3[i];

    double a00=dotprod<3>(edge0, edge0), a01=dotprod<3>(edge0,edge1), a11=dotprod<3>(edge1,edge1);
    double b0=dotprod<3>(diff, edge0), b1=dotprod<3>(diff, edge1), c=dotprod<3>(diff, diff);
    double det = fabs(a00*a11 - a01*a01);
    double s = a01*b1 - a11*b0, t = a01*b0 - a00*b1;
    double sDist;

    if (s + t <= det)
      {
        if (s < 0)  {
            if (t < 0) { // region 4
                if (b0 < 0) {
                    if (-b0 >= a00)  sDist = a00 + 2*b0 + c;
                    else             sDist = b0*(-b0 / a00) + c;
                  }
                else
                  sDist = _HelperDistancePtToTri3D_1(a11, b1, c);
              }
            else  // region 3
              sDist = _HelperDistancePtToTri3D_1(a11, b1, c);
          }
        else       {
            if (t < 0)  // region 5
              sDist = _HelperDistancePtToTri3D_1(a00, b0, c);
            else  // region 0
              {
                // minimum at interior point
                if (fabs(det) < 1.0e-12)
                  {
                    // points are colinear (degenerated triangle)
                    // => Compute distance between segments
                     double distance = std::min(DistanceFromPtToSegInSpaceDim3(pt, pt0Tri3, pt1Tri3),
                                                DistanceFromPtToSegInSpaceDim3(pt, pt1Tri3, pt2Tri3));
                     return distance;
                  }

                // else we can divide by non-zero
                double invDet = 1 / det;
                s *= invDet;    t *= invDet;
                sDist = s*(a00*s + a01*t + 2*b0) + t*(a01*s + a11*t + 2*b1) + c;
              }
          }
      }
    else  // s+t > det
      {
        if (s < 0.0)  // region 2
          sDist = _HelperDistancePtToTri3D_2(a01, a00, a11, b0, b1, c);
        else {
            if (t < 0.0)  // region 6
              sDist = _HelperDistancePtToTri3D_2(a01, a11, a00, b1, b0, c);
            else {  // region 1
                double numer = a11 + b1 - a01 - b0;
                if (numer <= 0.0)
                  sDist = a11 + 2*b1 + c;
                else {
                    double denom = a00 - 2*a01 + a11;
                    if (numer >= denom)
                      sDist = a00 + 2*b0 + c;
                    else {
                        s = numer / denom; t = 1 - s;
                        sDist = s*(a00*s + a01*t + 2*b0) + t*(a01*s + a11*t + 2*b1) + c;
                    }
                }
            }
        }
      }
    // Account for numerical round-off error.
    if (sDist < 0)
      sDist = 0.0;

    return sqrt(sDist);
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

