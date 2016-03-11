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

#include "InterpKernelMeshQuality.hxx"

#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>

double INTERP_KERNEL::quadSkew(const double *coo)
{
  double pa0[3]={
    coo[3]+coo[6]-coo[0]-coo[9],
    coo[4]+coo[7]-coo[1]-coo[10],
    coo[5]+coo[8]-coo[2]-coo[11]
  };
  double pa1[3]={
    coo[6]+coo[9]-coo[0]-coo[3],
    coo[7]+coo[10]-coo[1]-coo[4],
    coo[8]+coo[11]-coo[2]-coo[5],
  };
  double l0=sqrt(pa0[0]*pa0[0]+pa0[1]*pa0[1]+pa0[2]*pa0[2]);
  double l1=sqrt(pa1[0]*pa1[0]+pa1[1]*pa1[1]+pa1[2]*pa1[2]);
  if(l0<1.e-15)
    return 0.;
  if(l1<1.e-15)
    return 0.;
  pa0[0]/=l0; pa0[1]/=l0; pa0[2]/=l0;
  pa1[0]/=l1; pa1[1]/=l1; pa1[2]/=l1;
  return pa0[0]*pa1[0]+pa0[1]*pa1[1]+pa0[2]*pa1[2];
}

double INTERP_KERNEL::quadEdgeRatio(const double *coo)
{
  double a2=(coo[3]-coo[0])*(coo[3]-coo[0])+(coo[4]-coo[1])*(coo[4]-coo[1])+(coo[5]-coo[2])*(coo[5]-coo[2]);
  double b2=(coo[6]-coo[3])*(coo[6]-coo[3])+(coo[7]-coo[4])*(coo[7]-coo[4])+(coo[8]-coo[5])*(coo[8]-coo[5]);
  double c2=(coo[9]-coo[6])*(coo[9]-coo[6])+(coo[10]-coo[7])*(coo[10]-coo[7])+(coo[11]-coo[8])*(coo[11]-coo[8]);
  double d2=(coo[0]-coo[9])*(coo[0]-coo[9])+(coo[1]-coo[10])*(coo[1]-coo[10])+(coo[2]-coo[11])*(coo[2]-coo[11]);
  double mab=a2<b2?a2:b2;
  double Mab=a2<b2?b2:a2;
  double mcd=c2<d2?c2:d2;
  double Mcd=c2<d2?d2:c2;
  double m2=mab<mcd?mab:mcd;
  double M2=Mab>Mcd?Mab:Mcd;
  if(m2>1.e-15)
    return sqrt(M2/m2);
  else
    return std::numeric_limits<double>::max();
}

double INTERP_KERNEL::quadAspectRatio(const double *coo)
{
  double a=sqrt((coo[3]-coo[0])*(coo[3]-coo[0])+(coo[4]-coo[1])*(coo[4]-coo[1])+(coo[5]-coo[2])*(coo[5]-coo[2]));
  double b=sqrt((coo[6]-coo[3])*(coo[6]-coo[3])+(coo[7]-coo[4])*(coo[7]-coo[4])+(coo[8]-coo[5])*(coo[8]-coo[5]));
  double c=sqrt((coo[9]-coo[6])*(coo[9]-coo[6])+(coo[10]-coo[7])*(coo[10]-coo[7])+(coo[11]-coo[8])*(coo[11]-coo[8]));
  double d=sqrt((coo[0]-coo[9])*(coo[0]-coo[9])+(coo[1]-coo[10])*(coo[1]-coo[10])+(coo[2]-coo[11])*(coo[2]-coo[11]));
  double ma=a>b?a:b;
  double mb=c>d?c:d;
  double hm=ma>mb?ma:mb;
  double ab[3]={(coo[4]-coo[1])*(coo[8]-coo[5])-(coo[7]-coo[4])*(coo[5]-coo[2]),
                (coo[5]-coo[2])*(coo[6]-coo[3])-(coo[3]-coo[0])*(coo[8]-coo[5]),
                (coo[3]-coo[0])*(coo[7]-coo[4])-(coo[4]-coo[1])*(coo[6]-coo[3])};
  double cd[3]={(coo[10]-coo[7])*(coo[2]-coo[11])-(coo[1]-coo[10])*(coo[11]-coo[8]),
                (coo[11]-coo[8])*(coo[0]-coo[9])-(coo[9]-coo[6])*(coo[2]-coo[11]),
                (coo[9]-coo[6])*(coo[1]-coo[10])-(coo[10]-coo[7])*(coo[0]-coo[9])};
  double e=sqrt(ab[0]*ab[0]+ab[1]*ab[1]+ab[2]*ab[2])+sqrt(cd[0]*cd[0]+cd[1]*cd[1]+cd[2]*cd[2]);
  if(d>1e-15)
    return 0.5*(a+b+c+d)*hm/e;
  else
    return std::numeric_limits<double>::max();
}

double INTERP_KERNEL::quadWarp(const double *coo)
{
  double e0[3]={coo[3]-coo[0],coo[4]-coo[1],coo[5]-coo[2]};
  double e1[3]={coo[6]-coo[3],coo[7]-coo[4],coo[8]-coo[5]};
  double e2[3]={coo[9]-coo[6],coo[10]-coo[7],coo[11]-coo[8]};
  double e3[3]={coo[0]-coo[9],coo[1]-coo[10],coo[2]-coo[11]};
  
  double n0[3]={e3[1]*e0[2]-e3[2]*e0[1],e3[2]*e0[0]-e3[0]*e0[2],e3[0]*e0[1]-e3[1]*e0[0]};
  double n1[3]={e0[1]*e1[2]-e0[2]*e1[1],e0[2]*e1[0]-e0[0]*e1[2],e0[0]*e1[1]-e0[1]*e1[0]};
  double n2[3]={e1[1]*e2[2]-e1[2]*e2[1],e1[2]*e2[0]-e1[0]*e2[2],e1[0]*e2[1]-e1[1]*e2[0]};
  double n3[3]={e2[1]*e3[2]-e2[2]*e3[1],e2[2]*e3[0]-e2[0]*e3[2],e2[0]*e3[1]-e2[1]*e3[0]};

  double l0=sqrt(n0[0]*n0[0]+n0[1]*n0[1]+n0[2]*n0[2]);
  double l1=sqrt(n1[0]*n1[0]+n1[1]*n1[1]+n1[2]*n1[2]);
  double l2=sqrt(n2[0]*n2[0]+n2[1]*n2[1]+n2[2]*n2[2]);
  double l3=sqrt(n3[0]*n3[0]+n3[1]*n3[1]+n3[2]*n3[2]);

  if(l0<1.e-15 || l1<1.e-15 || l2<1.e-15 || l3<1e-15)
    return std::numeric_limits<double>::min();

  double warp=std::min(n0[0]/l0*n2[0]/l2+n0[1]/l0*n2[1]/l2+n0[2]/l0*n2[2]/l2,n1[0]/l1*n3[0]/l3+n1[1]/l1*n3[1]/l3+n1[2]/l1*n3[2]/l3);
  return warp*warp*warp;
}

double INTERP_KERNEL::triEdgeRatio(const double *coo)
{
  double a2=(coo[3]-coo[0])*(coo[3]-coo[0])+(coo[4]-coo[1])*(coo[4]-coo[1])+(coo[5]-coo[2])*(coo[5]-coo[2]);
  double b2=(coo[6]-coo[3])*(coo[6]-coo[3])+(coo[7]-coo[4])*(coo[7]-coo[4])+(coo[8]-coo[5])*(coo[8]-coo[5]);
  double c2=(coo[0]-coo[6])*(coo[0]-coo[6])+(coo[1]-coo[7])*(coo[1]-coo[7])+(coo[2]-coo[8])*(coo[2]-coo[8]);
  double mab=a2<b2?a2:b2;
  double Mab=a2<b2?b2:a2;
  double m2=c2>mab?mab:c2;
  double M2=c2>Mab?c2:Mab;
  if(m2>1.e-15)
    return sqrt(M2/m2);
  else
    return std::numeric_limits<double>::max();
}

double INTERP_KERNEL::triAspectRatio(const double *coo)
{
  double a=sqrt((coo[3]-coo[0])*(coo[3]-coo[0])+(coo[4]-coo[1])*(coo[4]-coo[1])+(coo[5]-coo[2])*(coo[5]-coo[2]));
  double b=sqrt((coo[6]-coo[3])*(coo[6]-coo[3])+(coo[7]-coo[4])*(coo[7]-coo[4])+(coo[8]-coo[5])*(coo[8]-coo[5]));
  double c=sqrt((coo[0]-coo[6])*(coo[0]-coo[6])+(coo[1]-coo[7])*(coo[1]-coo[7])+(coo[2]-coo[8])*(coo[2]-coo[8]));
 
  double hm=a>b?a:b;
  hm=hm>c?hm:c;

  double ab[3]={(coo[4]-coo[1])*(coo[8]-coo[5])-(coo[7]-coo[4])*(coo[5]-coo[2]),
                (coo[5]-coo[2])*(coo[6]-coo[3])-(coo[3]-coo[0])*(coo[8]-coo[5]),
                (coo[3]-coo[0])*(coo[7]-coo[4])-(coo[4]-coo[1])*(coo[6]-coo[3])};
  double d=sqrt(ab[0]*ab[0]+ab[1]*ab[1]+ab[2]*ab[2]);
  static const double normalizeCoeff=sqrt(3.)/6.;
  if(d>1.e-15) 
    return normalizeCoeff*hm*(a+b+c)/d;
  else
    return std::numeric_limits<double>::max();
}

double INTERP_KERNEL::tetraEdgeRatio(const double *coo)
{
  double a[3]={coo[3]-coo[0],coo[4]-coo[1],coo[5]-coo[2]};
  double b[3]={coo[6]-coo[3],coo[7]-coo[4],coo[8]-coo[5]};
  double c[3]={coo[0]-coo[6],coo[1]-coo[7],coo[2]-coo[8]};
  double d[3]={coo[9]-coo[0],coo[10]-coo[1],coo[11]-coo[2]};
  double e[3]={coo[9]-coo[3],coo[10]-coo[4],coo[11]-coo[5]};
  double f[3]={coo[9]-coo[6],coo[10]-coo[7],coo[11]-coo[8]};
  
  double l2[6]=
    {a[0]*a[0]+a[1]*a[1]+a[2]*a[2],
     b[0]*b[0]+b[1]*b[1]+b[2]*b[2],
     c[0]*c[0]+c[1]*c[1]+c[2]*c[2],
     d[0]*d[0]+d[1]*d[1]+d[2]*d[2],
     e[0]*e[0]+e[1]*e[1]+e[2]*e[2],
     f[0]*f[0]+f[1]*f[1]+f[2]*f[2]};

  double M2=*std::max_element(l2,l2+6);
  double m2=*std::min_element(l2,l2+6);
  if(m2>1e-15)
    return sqrt(M2/m2);
  else
    return std::numeric_limits<double>::max();
}

double INTERP_KERNEL::tetraAspectRatio(const double *coo)
{
  static const double normalizeCoeff=sqrt(6.)/12.;
  double ab[3]={coo[3]-coo[0],coo[4]-coo[1],coo[5]-coo[2]};
  double ac[3]={coo[6]-coo[0],coo[7]-coo[1],coo[8]-coo[2]};
  double ad[3]={coo[9]-coo[0],coo[10]-coo[1],coo[11]-coo[2]};
  double detTet=(ab[0]*(ac[1]*ad[2]-ac[2]*ad[1]))+(ab[1]*(ac[2]*ad[0]-ac[0]*ad[2]))+(ab[2]*(ac[0]*ad[1]-ac[1]*ad[0]));
  //if(detTet<1.e-15)
  //  return std::numeric_limits<double>::max();
  double bc[3]={coo[6]-coo[3],coo[7]-coo[4],coo[8]-coo[5]};
  double bd[3]={coo[9]-coo[3],coo[10]-coo[4],coo[11]-coo[5]};
  double cd[3]={coo[9]-coo[6],coo[10]-coo[7],coo[11]-coo[8]};

  double ab2=ab[0]*ab[0]+ab[1]*ab[1]+ab[2]*ab[2];
  double bc2=bc[0]*bc[0]+bc[1]*bc[1]+bc[2]*bc[2];
  double ac2=ac[0]*ac[0]+ac[1]*ac[1]+ac[2]*ac[2];
  double ad2=ad[0]*ad[0]+ad[1]*ad[1]+ad[2]*ad[2];
  double bd2=bd[0]*bd[0]+bd[1]*bd[1]+bd[2]*bd[2];
  double cd2=cd[0]*cd[0]+cd[1]*cd[1]+cd[2]*cd[2];

  double A=ab2>bc2?ab2:bc2;
  double B=ac2>ad2?ac2:ad2;
  double C=bd2>cd2?bd2:cd2;
  double D=A>B?A:B;
  double hm=D>C?sqrt(D):sqrt(C);

  bd[0]=ab[1]*bc[2]-ab[2]*bc[1]; bd[1]=ab[2]*bc[0]-ab[0]*bc[2]; bd[2]=ab[0]*bc[1]-ab[1]*bc[0];
  A=sqrt(bd[0]*bd[0]+bd[1]*bd[1]+bd[2]*bd[2]);
  bd[0]=ab[1]*ad[2]-ab[2]*ad[1]; bd[1]=ab[2]*ad[0]-ab[0]*ad[2]; bd[2]=ab[0]*ad[1]-ab[1]*ad[0];
  B=sqrt(bd[0]*bd[0]+bd[1]*bd[1]+bd[2]*bd[2]);
  bd[0]=ac[1]*ad[2]-ac[2]*ad[1]; bd[1]=ac[2]*ad[0]-ac[0]*ad[2]; bd[2]=ac[0]*ad[1]-ac[1]*ad[0];
  C=sqrt(bd[0]*bd[0]+bd[1]*bd[1]+bd[2]*bd[2]);
  bd[0]=bc[1]*cd[2]-bc[2]*cd[1]; bd[1]=bc[2]*cd[0]-bc[0]*cd[2]; bd[2]=bc[0]*cd[1]-bc[1]*cd[0];
  D=sqrt(bd[0]*bd[0]+bd[1]*bd[1]+bd[2]*bd[2]);
  return normalizeCoeff*hm*(A+B+C+D)/fabs(detTet);
}
