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

#ifndef __GENMATHFORMULAE_HXX__
#define __GENMATHFORMULAE_HXX__

#include "InterpKernelException.hxx"

#include <cmath>

namespace INTERP_KERNEL
{
  /*!
   * This method computes eigenvalues of a 3x3 symetric matrix stored with 6 values in 'matrix'. The convension chosen for 'matrix' is described here:
   * matrix[0]=m_xx, matrix[1]=m_yy, matrix[2]=m_zz,
   * matrix[3]=m_xy, matrix[4]=m_yz, matrix[5]=m_xz
   * This method returns the 3 eigenvalues in 'eigenVals'.
   */
  void computeEigenValues6(const double *matrix, double *eigenVals)
  {
    double tr=(matrix[0]+matrix[1]+matrix[2])/3.;
    double K[6]={matrix[0]-tr,matrix[1]-tr,matrix[2]-tr,matrix[3],matrix[4],matrix[5]};
    double q=(K[0]*K[1]*K[2]+2.*K[4]*K[5]*K[3]-K[0]*K[4]*K[4]-K[2]*K[3]*K[3]-K[1]*K[5]*K[5])/2.;
    double p=K[0]*K[0]+K[1]*K[1]+K[2]*K[2]+2*(K[3]*K[3]+K[4]*K[4]+K[5]*K[5]);
    p/=6.;
    double sqp=sqrt(p);
    double tmp=p*sqp;
    double phi;
    if(fabs(q)<=fabs(tmp))
      phi=1./3.*acos(q/tmp);
    else
      phi=0.;
    if(phi<0.)
      phi+=M_PI/3.;
    eigenVals[0]=tr+2.*sqp*cos(phi);
    eigenVals[1]=tr-sqp*(cos(phi)+sqrt(3.)*sin(phi));
    eigenVals[2]=tr-sqp*(cos(phi)-sqrt(3.)*sin(phi));
  }
  
  /*!
   * This method computes one eigenvector of a 3x3 symetric matrix stored with 6 values in 'matrix'. The convension chosen for 'matrix' is described here:
   * matrix[0]=m_xx, matrix[1]=m_yy, matrix[2]=m_zz,
   * matrix[3]=m_xy, matrix[4]=m_yz, matrix[5]=m_xz
   * This method returns the eigenvector of the corresponding eigenvalue in 'eigenVal'. The returned eigenValue is normalized.
   */
  void computeEigenVectorForEigenValue6(const double *matrix, double eigenVal, double eps, double *eigenVector)
  {
    //if(fabs(eigenVal)>eps)
      {
        const double m9[9]={matrix[0]-eigenVal,matrix[3],matrix[5],matrix[3],matrix[1]-eigenVal,matrix[4],matrix[5],matrix[4],matrix[2]-eigenVal};
        for(int i=0;i<3;i++)
          {
            double w[9]={m9[0+3*i],m9[1+3*i],m9[2+3*i],m9[0+(3*(i+1))%6],m9[1+(3*(i+1))%6],m9[2+(3*(i+1))%6],1.,1.,1.};
            double det=w[0]*w[4]*w[8]+w[1]*w[5]*w[6]+w[2]*w[3]*w[7]-w[0]*w[5]*w[7]-w[1]*w[3]*w[8]-w[2]*w[4]*w[6];
            if(fabs(det)>eps)
              {
                eigenVector[0]=(w[1]*w[5]-w[4]*w[2])/det;
                eigenVector[1]=(w[2]*w[3]-w[0]*w[5])/det;
                eigenVector[2]=(w[0]*w[4]-w[1]*w[3])/det;
                double norm=sqrt(eigenVector[0]*eigenVector[0]+eigenVector[1]*eigenVector[1]+eigenVector[2]*eigenVector[2]);
                eigenVector[0]/=norm;
                eigenVector[1]/=norm;
                eigenVector[2]/=norm;
                return;
              }
          }
      }
      //else
      {
        eigenVector[0]=0.;
        eigenVector[1]=0.;
        eigenVector[2]=0.;
        return;
      }
      //throw INTERP_KERNEL::Exception("computeEigenVector : Do not succed in finding eigen vector !");
  }
}

#endif
