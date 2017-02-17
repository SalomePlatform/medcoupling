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

#ifndef __TRANSLATIONROTATIONMATRIX_HXX__
#define __TRANSLATIONROTATIONMATRIX_HXX__

#include "INTERPKERNELDefines.hxx"

#include <cmath>

namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT TranslationRotationMatrix
  {

  public:

    TranslationRotationMatrix()
    { 
      unsigned i;
      for(i=0;i<TRANSL_SIZE;i++)
        _translation_coeffs[i]=0.;
      for(i=0;i<ROT_SIZE;i++)
        _rotation_coeffs[i]=i%4?0.:1.;
    }

    void multiply(const TranslationRotationMatrix& A)
    {
      TranslationRotationMatrix result;
      //setting output matrix to zero
      for (int i=0; i<3; i++)
        result._rotation_coeffs[i*4]=0.0;
      //multiplying
      for (int i=0; i<3;i++)
        for (int j=0; j<3; j++)
          for (int k=0; k<3; k++)
            result._rotation_coeffs[j+i*3]+=A._rotation_coeffs[3*i+k]*_rotation_coeffs[j+k*3];
    
      for (int i=0;i<9; i++)
        _rotation_coeffs[i]=result._rotation_coeffs[i];
    }
  
    void rotate_vector(double* P)
    {
      double temp[3]={0.0, 0.0, 0.0};
    
      for (int i=0; i<3;i++)
        for (int j=0; j<3; j++)
          temp[i] +=_rotation_coeffs[3*i+j]*P[j];
       
      P[0]=temp[0];P[1]=temp[1];P[2]=temp[2];
    }
 
    void transform_vector(double*P)
    {
      P[0]+=_translation_coeffs[0];
      P[1]+=_translation_coeffs[1];
      P[2]+=_translation_coeffs[2];
      rotate_vector(P);
    }

    void translate(const double* P)
    {
      _translation_coeffs[0]=P[0];
      _translation_coeffs[1]=P[1];
      _translation_coeffs[2]=P[2];
    }
  
    void  rotate_x (double* P)
    {
      _rotation_coeffs[0]=1.0;
      double r_sqr = P[1]*P[1]+P[2]*P[2];
      if (r_sqr < EPS)
        {_rotation_coeffs[4]=1.0; _rotation_coeffs[8]=1.0; return;}
      double r = sqrt(r_sqr);
      double cos =P[1]/r;
      double sin =P[2]/r;

      _rotation_coeffs[4]=cos;
      _rotation_coeffs[5]=sin;
      _rotation_coeffs[7]=-sin;
      _rotation_coeffs[8]=cos;


      rotate_vector(P);
    }

    void  rotate_z (double* P)
    {
      _rotation_coeffs[8]=1.0;
      double r_sqr = P[0]*P[0]+P[1]*P[1];
      if (r_sqr < EPS)
        {_rotation_coeffs[4]=1.0; _rotation_coeffs[0]=1.0; return;}
      double r = sqrt(r_sqr);
      double cos =P[0]/r;
      double sin =P[1]/r;
    
      _rotation_coeffs[0]=cos;
      _rotation_coeffs[1]=sin; 
      _rotation_coeffs[3]=-sin;
      _rotation_coeffs[4]=cos;
    
      rotate_vector(P);
    }
                     
    /**
     * Fills in rotation_matrix so that the triangle PP1, PP2, PP3 is in the plane Oxy, with PP1 at the origin and PP2 on Ox.
     */
    static void Rotate3DTriangle(const double* PP1,const double*PP2,const double*PP3, TranslationRotationMatrix& rotation_matrix)
    {
      double P1w[3];
      double P2w[3];
      double P3w[3];
      P2w[0]=PP2[0]; P2w[1]=PP2[1];P2w[2]=PP2[2];
      P3w[0]=PP3[0]; P3w[1]=PP3[1];P3w[2]=PP3[2];

      // translating to set P1 at the origin
      for (int i=0; i<3; i++)
        {
          P1w[i] = -PP1[i];
          P2w[i] -= PP1[i];
          P3w[i] -= PP1[i];
        }

      //initializes matrix
      rotation_matrix.translate(P1w);

      // rotating to set P2 on the Oxy plane
      TranslationRotationMatrix A;
      A.rotate_x(P2w);
      A.rotate_vector(P3w);
      rotation_matrix.multiply(A);

      //rotating to set P2 on the Ox axis
      TranslationRotationMatrix B;
      B.rotate_z(P2w);
      B.rotate_vector(P3w);
      rotation_matrix.multiply(B);

      //rotating to set P3 on the Oxy plane
      TranslationRotationMatrix C;
      C.rotate_x(P3w);
      rotation_matrix.multiply(C);
    }

    /**
     * Fills in rotation_matrix so that the bipoint PP1, PP2 is on the Ox axis.
     */
    static void Rotate3DBipoint(const double* PP1,const double*PP2, TranslationRotationMatrix& rotation_matrix)
    {
      double P1w[3];
      double P2w[3];
      P2w[0]=PP2[0]; P2w[1]=PP2[1];P2w[2]=PP2[2];

      // translating to set P1 at the origin
      for (int i=0; i<3; i++)
        {
          P1w[i] = -PP1[i];
          P2w[i] -= PP1[i];
        }

      //initializes
      rotation_matrix.translate(P1w);

      // rotating to set P2 on the Oxy plane
      TranslationRotationMatrix A;
      A.rotate_x(P2w);
      rotation_matrix.multiply(A);

      //rotating to set P2 on the Ox axis
      TranslationRotationMatrix B;
      B.rotate_z(P2w);
      rotation_matrix.multiply(B);
    }



  private:
    static const double EPS;
    static const unsigned ROT_SIZE=9;
    static const unsigned TRANSL_SIZE=3;
    double _rotation_coeffs[ROT_SIZE];
    double _translation_coeffs[TRANSL_SIZE];
  };
}

#endif
