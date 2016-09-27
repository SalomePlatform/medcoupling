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

#ifndef __INTERPOLATIONUTILS_HXX__
#define __INTERPOLATIONUTILS_HXX__

#include "INTERPKERNELDefines.hxx"
#include "InterpKernelException.hxx"
#include "VolSurfUser.hxx"

#include "NormalizedUnstructuredMesh.hxx"

#include <deque>
#include <map>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <functional>

namespace INTERP_KERNEL
{
  template<class ConnType, NumberingPolicy numPol>
  class OTT//OffsetToolTrait
  {
  };
  
  template<class ConnType>
  class OTT<ConnType,ALL_C_MODE>
  {
  public:
    static ConnType indFC(ConnType i) { return i; }
    static ConnType ind2C(ConnType i) { return i; }
    static ConnType conn2C(ConnType i) { return i; }
    static ConnType coo2C(ConnType i) { return i; }
  };
  
  template<class ConnType>
  class OTT<ConnType,ALL_FORTRAN_MODE>
  {
  public:
    static ConnType indFC(ConnType i) { return i+1; }
    static ConnType ind2C(ConnType i) { return i-1; }
    static ConnType conn2C(ConnType i) { return i-1; }
    static ConnType coo2C(ConnType i) { return i-1; }
  };

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */
  /*   calcul la surface d'un triangle                  */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */

  inline double Surf_Tri(const double* P_1,const double* P_2,const double* P_3)
  {
    double A=(P_3[1]-P_1[1])*(P_2[0]-P_1[0])-(P_2[1]-P_1[1])*(P_3[0]-P_1[0]);
    double Surface = 0.5*fabs(A);
    return Surface;
  }

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */
  /*     fonction qui calcul le determinant            */
  /*      de deux vecteur(cf doc CGAL).                */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/

  //fonction qui calcul le determinant des vecteurs: P3P1 et P3P2
  //(cf doc CGAL).

  inline double mon_determinant(const double* P_1,
                                const double*  P_2,
                                const double* P_3)
  {
    double mon_det=(P_1[0]-P_3[0])*(P_2[1]-P_3[1])-(P_2[0]-P_3[0])*(P_1[1]-P_3[1]);
    return mon_det;
  }

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  //calcul la norme du vecteur P1P2

  inline double norme_vecteur(const double* P_1,const double* P_2)
  {
    double X=P_1[0]-P_2[0];
    double Y=P_1[1]-P_2[1];
    double norme=sqrt(X*X+Y*Y);
    return norme;
  }

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */
  /*         calcul le cos et le sin de l'angle P1P2,P1P3     */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */

  inline std::vector<double> calcul_cos_et_sin(const double* P_1,
                                               const double* P_2,
                                               const double* P_3)
  {
       
    std::vector<double> Vect;
    double P1_P2=norme_vecteur(P_1,P_2);
    double P2_P3=norme_vecteur(P_2,P_3);
    double P3_P1=norme_vecteur(P_3,P_1);
       
    double N=P1_P2*P1_P2+P3_P1*P3_P1-P2_P3*P2_P3;
    double D=2.0*P1_P2*P3_P1;
    double COS=N/D;
    if (COS>1.0) COS=1.0;
    if (COS<-1.0) COS=-1.0;
    Vect.push_back(COS);
    double V=mon_determinant(P_2,P_3,P_1);
    double D_1=P1_P2*P3_P1;
    double SIN=V/D_1;
    if (SIN>1.0) SIN=1.0;
    if (SIN<-1.0) SIN=-1.0;
    Vect.push_back(SIN);
       
    return Vect;
       
  }

  /*!
   * This method builds a quadrangle built with the first point of 'triIn' the barycenter of two edges starting or ending with
   * the first point of 'triIn' and the barycenter of 'triIn'.
   *
   * @param triIn is a 6 doubles array in full interlace mode, that represents a triangle.
   * @param quadOut is a 8 doubles array filled after the following call.
   */
  template<int SPACEDIM>
  inline void fillDualCellOfTri(const double *triIn, double *quadOut)
  {
    //1st point
    std::copy(triIn,triIn+SPACEDIM,quadOut);
    double tmp[SPACEDIM];
    std::transform(triIn,triIn+SPACEDIM,triIn+SPACEDIM,tmp,std::plus<double>());
    //2nd point
    std::transform(tmp,tmp+SPACEDIM,quadOut+SPACEDIM,std::bind2nd(std::multiplies<double>(),0.5));
    std::transform(tmp,tmp+SPACEDIM,triIn+2*SPACEDIM,tmp,std::plus<double>());
    //3rd point
    std::transform(tmp,tmp+SPACEDIM,quadOut+2*SPACEDIM,std::bind2nd(std::multiplies<double>(),1/3.));
    //4th point
    std::transform(triIn,triIn+SPACEDIM,triIn+2*SPACEDIM,tmp,std::plus<double>());
    std::transform(tmp,tmp+SPACEDIM,quadOut+3*SPACEDIM,std::bind2nd(std::multiplies<double>(),0.5));
  }

  /*!
   * This method builds a potentially non-convex polygon cell built with the first point of 'triIn' the barycenter of two edges starting or ending with
   * the first point of 'triIn' and the barycenter of 'triIn'.
   *
   * @param triIn is a 6 doubles array in full interlace mode, that represents a triangle.
   * @param quadOut is a 8 doubles array filled after the following call.
   */
  template<int SPACEDIM>
  inline void fillDualCellOfPolyg(const double *polygIn, int nPtsPolygonIn, double *polygOut)
  {
    //1st point
    std::copy(polygIn,polygIn+SPACEDIM,polygOut);
    std::transform(polygIn,polygIn+SPACEDIM,polygIn+SPACEDIM,polygOut+SPACEDIM,std::plus<double>());
    //2nd point
    std::transform(polygOut+SPACEDIM,polygOut+2*SPACEDIM,polygOut+SPACEDIM,std::bind2nd(std::multiplies<double>(),0.5));
    double tmp[SPACEDIM];
    //
    for(int i=0;i<nPtsPolygonIn-2;i++)
      {
        std::transform(polygIn,polygIn+SPACEDIM,polygIn+(i+2)*SPACEDIM,tmp,std::plus<double>());
        std::transform(tmp,tmp+SPACEDIM,polygOut+(2*i+3)*SPACEDIM,std::bind2nd(std::multiplies<double>(),0.5));
        std::transform(polygIn+(i+1)*SPACEDIM,polygIn+(i+2)*SPACEDIM,tmp,tmp,std::plus<double>());
        std::transform(tmp,tmp+SPACEDIM,polygOut+(2*i+2)*SPACEDIM,std::bind2nd(std::multiplies<double>(),1./3.));
      }
  }

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */
  /*     calcul les coordonnees du barycentre d'un polygone   */ 
  /*     le vecteur en entree est constitue des coordonnees   */
  /*     des sommets du polygone                              */                             
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */

  inline std::vector<double> bary_poly(const std::vector<double>& V)
  {
    std::vector<double> Bary;
    long taille=V.size();
    double x=0;
    double y=0;

    for(long i=0;i<taille/2;i++)
      {
        x=x+V[2*i];
        y=y+V[2*i+1];
      }
    double A=2*x/((double)taille);
    double B=2*y/((double)taille);
    Bary.push_back(A);//taille vecteur=2*nb de points.
    Bary.push_back(B);


    return Bary;
  }
  
  /*!
   * Given 6 coeffs of a Tria6 returns the corresponding value of a given pos
   */
  inline double computeTria6RefBase(const double *coeffs, const double *pos)
  {
    return coeffs[0]+coeffs[1]*pos[0]+coeffs[2]*pos[1]+coeffs[3]*pos[0]*pos[0]+coeffs[4]*pos[0]*pos[1]+coeffs[5]*pos[1]*pos[1];
  }
  
  /*!
   * Given xsi,eta in refCoo (length==2) return 6 coeffs in weightedPos.
   */
  inline void computeWeightedCoeffsInTria6FromRefBase(const double *refCoo, double *weightedPos)
  {
    weightedPos[0]=(1.-refCoo[0]-refCoo[1])*(1.-2*refCoo[0]-2.*refCoo[1]);
    weightedPos[1]=refCoo[0]*(2.*refCoo[0]-1.);
    weightedPos[2]=refCoo[1]*(2.*refCoo[1]-1.);
    weightedPos[3]=4.*refCoo[0]*(1.-refCoo[0]-refCoo[1]);
    weightedPos[4]=4.*refCoo[0]*refCoo[1];
    weightedPos[5]=4.*refCoo[1]*(1.-refCoo[0]-refCoo[1]);
  }

  /*!
   * Given 10 coeffs of a Tetra10 returns the corresponding value of a given pos
   */
  inline double computeTetra10RefBase(const double *coeffs, const double *pos)
  {
    return coeffs[0]+coeffs[1]*pos[0]+coeffs[2]*pos[1]+coeffs[3]*pos[2]+
      coeffs[4]*pos[0]*pos[0]+coeffs[5]*pos[0]*pos[1]+coeffs[6]*pos[0]*pos[2]+
      coeffs[7]*pos[1]*pos[1]+coeffs[8]*pos[1]*pos[2]+coeffs[9]*pos[2]*pos[2];
  }

  /*!
   * Given xsi,eta,z in refCoo (length==3) return 10 coeffs in weightedPos.
   */
  inline void computeWeightedCoeffsInTetra10FromRefBase(const double *refCoo, double *weightedPos)
  {
    //http://www.cadfamily.com/download/CAE/ABAQUS/The%20Finite%20Element%20Method%20-%20A%20practical%20course%20abaqus.pdf page 217
    //L1=1-refCoo[0]-refCoo[1]-refCoo[2]
    //L2=refCoo[0] L3=refCoo[1] L4=refCoo[2]
    weightedPos[0]=(-2.*(refCoo[0]+refCoo[1]+refCoo[2])+1)*(1-refCoo[0]-refCoo[1]-refCoo[2]);//(2*L1-1)*L1
    weightedPos[1]=(2.*refCoo[0]-1.)*refCoo[0];//(2*L2-1)*L2
    weightedPos[2]=(2.*refCoo[1]-1.)*refCoo[1];//(2*L3-1)*L3
    weightedPos[3]=(2.*refCoo[2]-1.)*refCoo[2];//(2*L4-1)*L4
    weightedPos[4]=4.*(1-refCoo[0]-refCoo[1]-refCoo[2])*refCoo[0];//4*L1*L2
    weightedPos[5]=4.*refCoo[0]*refCoo[1];//4*L2*L3
    weightedPos[6]=4.*(1-refCoo[0]-refCoo[1]-refCoo[2])*refCoo[1];//4*L1*L3
    weightedPos[7]=4.*(1-refCoo[0]-refCoo[1]-refCoo[2])*refCoo[2];//4*L1*L4
    weightedPos[8]=4.*refCoo[0]*refCoo[2];//4*L2*L4
    weightedPos[9]=4.*refCoo[1]*refCoo[2];//4*L3*L4
  }

  /*!
   * \brief Solve system equation in matrix form using Gaussian elimination algorithm
   *  \param M - N x N+1 matrix
   *  \param sol - vector of N solutions
   *  \retval bool - true if succeeded
   */
  template<unsigned nbRow>
  bool solveSystemOfEquations(double M[nbRow][nbRow+1], double* sol)
  {
    const int nbCol=nbRow+1;

    // make upper triangular matrix (forward elimination)

    int iR[nbRow];// = { 0, 1, 2 };
    for ( int i = 0; i < (int) nbRow; ++i )
      iR[i] = i;
    for ( int i = 0; i < (int)(nbRow-1); ++i ) // nullify nbRow-1 rows
      {
        // swap rows to have max value of i-th column in i-th row
        double max = std::fabs( M[ iR[i] ][i] );
        for ( int r = i+1; r < (int)nbRow; ++r )
          {
            double m = std::fabs( M[ iR[r] ][i] );
            if ( m > max )
              {
                max = m;
                std::swap( iR[r], iR[i] );
              }
          }
        if ( max < std::numeric_limits<double>::min() )
          {
            //sol[0]=1; sol[1]=sol[2]=sol[3]=0;
            return false; // no solution
          }
        // make 0 below M[i][i] (actually we do not modify i-th column)
        double* tUpRow = M[ iR[i] ];
        for ( int r = i+1; r < (int)nbRow; ++r )
          {
            double* mRow = M[ iR[r] ];
            double coef = mRow[ i ] / tUpRow[ i ];
            for ( int c = i+1; c < nbCol; ++c )
              mRow[ c ] -= tUpRow[ c ] * coef;
          }
      }
    double* mRow = M[ iR[nbRow-1] ];
    if ( std::fabs( mRow[ nbRow-1 ] ) < std::numeric_limits<double>::min() )
      {
        //sol[0]=1; sol[1]=sol[2]=sol[3]=0;
        return false; // no solution
      }
    mRow[ nbRow ] /= mRow[ nbRow-1 ];

    // calculate solution (back substitution)

    sol[ nbRow-1 ] = mRow[ nbRow ];

    for ( int i = nbRow-2; i+1; --i )
      {
        mRow = M[ iR[i] ];
        sol[ i ] = mRow[ nbRow ];
        for ( int j = nbRow-1; j > i; --j )
          sol[ i ] -= sol[j]*mRow[ j ];
        sol[ i ] /= mRow[ i ];
      }

    return true;
  }

  
  /*!
   * \brief Solve system equation in matrix form using Gaussian elimination algorithm
   *  \param M - N x N+NB_OF_VARS matrix
   *  \param sol - vector of N solutions
   *  \retval bool - true if succeeded
   */
  template<unsigned SZ, unsigned NB_OF_RES>
  bool solveSystemOfEquations2(const double *matrix, double *solutions, double eps)
  {
    unsigned k,j;
    int nr,n,m,np;
    double s,g;
    int mb;
    //
    double B[SZ*(SZ+NB_OF_RES)];
    std::copy(matrix,matrix+SZ*(SZ+NB_OF_RES),B);
    //
    nr=SZ+NB_OF_RES;
    for(k=0;k<SZ;k++)
      {
        np=nr*k+k;
        if(fabs(B[np])<eps)
          {
            n=k;
            do
              {
                n++;
                if(fabs(B[nr*k+n])>eps)
                  {/* Rows permutation */
                    for(m=0;m<nr;m++)
                      std::swap(B[nr*k+m],B[nr*n+m]);
                  }
              }
            while (n<(int)SZ);
          }
        s=B[np];//s is the Pivot
        std::transform(B+k*nr,B+(k+1)*nr,B+k*nr,std::bind2nd(std::divides<double>(),s));
        for(j=0;j<SZ;j++)
          {
            if(j!=k)
              {
                g=B[j*nr+k];
                for(mb=k;mb<nr;mb++)
                  B[j*nr+mb]-=B[k*nr+mb]*g;
              }
          }
      }
    for(j=0;j<NB_OF_RES;j++)
      for(k=0;k<SZ;k++)
        solutions[j*SZ+k]=B[nr*k+SZ+j];
    //
    return true;
  }

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  /*     Calculate barycentric coordinates of a 2D point p */ 
  /*     with respect to the triangle verices.             */
  /*     triaCoords are in full interlace                  */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/

  template<int SPACEDIM>
  inline void barycentric_coords(const double* triaCoords, const double* p, double* bc)
  {
    // matrix 2x2
    double
      T11 = triaCoords[0]-triaCoords[2*SPACEDIM], T12 = triaCoords[SPACEDIM]-triaCoords[2*SPACEDIM],
      T21 = triaCoords[1]-triaCoords[2*SPACEDIM+1], T22 = triaCoords[SPACEDIM+1]-triaCoords[2*SPACEDIM+1];
    // matrix determinant
    double Tdet = T11*T22 - T12*T21;
    if ( fabs( Tdet ) < std::numeric_limits<double>::min() ) {
      bc[0]=1; bc[1]=0; bc[2]=0;
      return;
    }
    // matrix inverse
    double t11 = T22, t12 = -T12, t21 = -T21, t22 = T11;
    // vector
    double r11 = p[0]-triaCoords[2*SPACEDIM], r12 = p[1]-triaCoords[2*SPACEDIM+1];
    // barycentric coordinates: mutiply matrix by vector
    bc[0] = (t11 * r11 + t12 * r12)/Tdet;
    bc[1] = (t21 * r11 + t22 * r12)/Tdet;
    bc[2] = 1. - bc[0] - bc[1];
  }

  /*!
   * Calculate barycentric coordinates of a point p with respect to triangle or tetra vertices.
   * This method makes 2 assumptions :
   *    - this is a simplex
   *    - spacedim == meshdim. For TRI3 and TRI6 spaceDim is expected to be equal to 2 and for TETRA4 spaceDim is expected to be equal to 3.
   *      If not the case (3D surf for example) a previous projection should be done before.
   */
  inline void barycentric_coords(const std::vector<const double*>& n, const double *p, double *bc)
  {
    enum { _XX=0, _YY, _ZZ };
    switch(n.size())
      {
      case 2:
        {// SEG 2
          double delta=n[0][0]-n[1][0];
          bc[0]=fabs((*p-n[1][0])/delta);
          bc[1]=fabs((*p-n[0][0])/delta);
          break;
        }
      case 3:
        { // TRIA3
          // matrix 2x2
          double
            T11 = n[0][_XX]-n[2][_XX], T12 = n[1][_XX]-n[2][_XX],
            T21 = n[0][_YY]-n[2][_YY], T22 = n[1][_YY]-n[2][_YY];
          // matrix determinant
          double Tdet = T11*T22 - T12*T21;
          if ( (std::fabs( Tdet) ) < (std::numeric_limits<double>::min()) )
            {
              bc[0]=1; bc[1]=bc[2]=0; // no solution
              return;
            }
          // matrix inverse
          double t11 = T22, t12 = -T12, t21 = -T21, t22 = T11;
          // vector
          double r11 = p[_XX]-n[2][_XX], r12 = p[_YY]-n[2][_YY];
          // barycentric coordinates: mutiply matrix by vector
          bc[0] = (t11 * r11 + t12 * r12)/Tdet;
          bc[1] = (t21 * r11 + t22 * r12)/Tdet;
          bc[2] = 1. - bc[0] - bc[1];
          break;
        }
      case 4:
        { // TETRA4
          // Find bc by solving system of 3 equations using Gaussian elimination algorithm
          // bc1*( x1 - x4 ) + bc2*( x2 - x4 ) + bc3*( x3 - x4 ) = px - x4
          // bc1*( y1 - y4 ) + bc2*( y2 - y4 ) + bc3*( y3 - y4 ) = px - y4
          // bc1*( z1 - z4 ) + bc2*( z2 - z4 ) + bc3*( z3 - z4 ) = px - z4
          
          double T[3][4]=
            {{ n[0][_XX]-n[3][_XX], n[1][_XX]-n[3][_XX], n[2][_XX]-n[3][_XX], p[_XX]-n[3][_XX] },
             { n[0][_YY]-n[3][_YY], n[1][_YY]-n[3][_YY], n[2][_YY]-n[3][_YY], p[_YY]-n[3][_YY] },
             { n[0][_ZZ]-n[3][_ZZ], n[1][_ZZ]-n[3][_ZZ], n[2][_ZZ]-n[3][_ZZ], p[_ZZ]-n[3][_ZZ] }};
          
          if ( !solveSystemOfEquations<3>( T, bc ) )
            bc[0]=1., bc[1] = bc[2] = bc[3] = 0;
          else
            bc[ 3 ] = 1. - bc[0] - bc[1] - bc[2];
          break;
        }
      case 6:
        {
          // TRIA6
          double matrix2[48]={1., 0., 0., 0., 0., 0., 0., 0.,
                              1., 0., 0., 0., 0., 0., 1., 0., 
                              1., 0., 0., 0., 0., 0., 0., 1.,
                              1., 0., 0., 0., 0., 0., 0.5, 0.,
                              1., 0., 0., 0., 0., 0., 0.5, 0.5,
                              1., 0., 0., 0., 0., 0., 0.,0.5};
          for(int i=0;i<6;i++)
            {
              matrix2[8*i+1]=n[i][0];
              matrix2[8*i+2]=n[i][1];
              matrix2[8*i+3]=n[i][0]*n[i][0];
              matrix2[8*i+4]=n[i][0]*n[i][1];
              matrix2[8*i+5]=n[i][1]*n[i][1];
            }
          double res[12];
          solveSystemOfEquations2<6,2>(matrix2,res,std::numeric_limits<double>::min());
          double refCoo[2];
          refCoo[0]=computeTria6RefBase(res,p);
          refCoo[1]=computeTria6RefBase(res+6,p);
          computeWeightedCoeffsInTria6FromRefBase(refCoo,bc);
          break;
        }
      case 10:
        {//TETRA10
          double matrix2[130]={1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                               1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,
                               1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,
                               1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,
                               1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0., 0.,
                               1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0.5, 0.,
                               1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,0.5, 0.,
                               1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5,
                               1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0., 0.5,
                               1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.5, 0.5};
          for(int i=0;i<10;i++)
            {
              matrix2[13*i+1]=n[i][0];
              matrix2[13*i+2]=n[i][1];
              matrix2[13*i+3]=n[i][2];
              matrix2[13*i+4]=n[i][0]*n[i][0];
              matrix2[13*i+5]=n[i][0]*n[i][1];
              matrix2[13*i+6]=n[i][0]*n[i][2];
              matrix2[13*i+7]=n[i][1]*n[i][1];
              matrix2[13*i+8]=n[i][1]*n[i][2];
              matrix2[13*i+9]=n[i][2]*n[i][2];
            }
          double res[30];
          solveSystemOfEquations2<10,3>(matrix2,res,std::numeric_limits<double>::min());
          double refCoo[3];
          refCoo[0]=computeTetra10RefBase(res,p);
          refCoo[1]=computeTetra10RefBase(res+10,p);
          refCoo[2]=computeTetra10RefBase(res+20,p);
          computeWeightedCoeffsInTetra10FromRefBase(refCoo,bc);
          break;
        }
      default:
        throw INTERP_KERNEL::Exception("INTERP_KERNEL::barycentric_coords : unrecognized simplex !");
      }
  }

  /*!
   * Calculate pseudo barycentric coordinates of a point p with respect to the quadrangle vertices.
   * This method makes the assumption that:
   *  - spacedim == meshdim (2 here).
   *  - the point is within the quad
   *  Quadratic elements are not supported yet.
   *
   *  A quadrangle can be described as 3 vectors, one point being taken as the origin.
   *  Denoting A, B, C the three other points, any point P within the quad is written as
   *    P = xA+ yC + xy(B-A-C)
   *  This method solve those 2 equations (one per component) for x and y.
   *

          A------B
          |      |
          |      |
          0------C
   */
  inline void quad_mapped_coords(const std::vector<const double*>& n, const double *p, double *bc)
  {
    double prec = 1.0e-14;
    enum { _XX=0, _YY, _ZZ };

    if(n.size() != 4)
      throw INTERP_KERNEL::Exception("INTERP_KERNEL::quad_mapped_coords : unrecognized geometric type! Only QUAD4 supported.");

    double A[2] = {n[1][_XX] - n[0][_XX],  n[1][_YY] - n[0][_YY]};
    double B[2] = {n[2][_XX] - n[0][_XX],  n[2][_YY] - n[0][_YY]};
    double C[2] = {n[3][_XX] - n[0][_XX],  n[3][_YY] - n[0][_YY]};
    double N[2] = {B[_XX] - A[_XX] - C[_XX], B[_YY] - A[_YY] - C[_YY]};
    double P[2] = {p[_XX] - n[0][_XX], p[_YY] - n[0][_YY]};

    // degenerated case: a rectangle:
    if (fabs(N[0]) < prec && fabs(N[1]) < prec)
      {
        double det = C[0]*A[1] -C[1]*A[0];
        if (fabs(det) < prec)
          throw INTERP_KERNEL::Exception("MappedBarycentric intersection type: quad_mapped_coords() has a degenerated 2x2 system!");
        bc[0] = (P[0]*A[1]-P[1]*A[0])/det;
        bc[1] = (P[1]*C[0]-P[0]*C[1])/det;
        return;
      }
    double b,c ,a = A[1]*N[0]-A[0]*N[1];
    bool cas1;
    if (fabs(a) > 1.0e-14)
      {
        b = A[1]*C[0]+N[1]*P[0]-N[0]*P[1]-A[0]*C[1];
        c = P[0]*C[1] - P[1]*C[0];
        cas1 = true;
      }
    else
      {
        a = -C[1]*N[0]+C[0]*N[1];
        b = A[1]*C[0]-N[1]*P[0]+N[0]*P[1]-A[0]*C[1];
        c = -P[0]*A[1] + P[1]*A[0];
        cas1 = false;
      }
    double delta = b*b - 4.0*a*c;
    if (delta < 0.0)
      throw INTERP_KERNEL::Exception("MappedBarycentric intersection type: quad_mapped_coords(): imaginary solutions!");
    bc[1] = 0.5*(-b+sqrt(delta))/a;
    if (bc[1] < -prec || bc[1] > (1.0+prec))
      bc[1] = 0.5*(-b-sqrt(delta))/a;
    if (bc[1] < -prec || bc[1] > (1.0+prec))
      throw INTERP_KERNEL::Exception("MappedBarycentric intersection type: quad_mapped_coords(): point doesn't seem to be in quad4!");
    if (cas1)
      {
        double denom = C[0]+bc[1]*N[0];
        if (fabs(denom) < prec)
          throw INTERP_KERNEL::Exception("MappedBarycentric intersection type: quad_mapped_coords(): point doesn't seem to be in quad4!");
        bc[0] = (P[0]-bc[1]*A[0])/denom;
        if (bc[0] < -prec || bc[0] > (1.0+prec))
          throw INTERP_KERNEL::Exception("MappedBarycentric intersection type: quad_mapped_coords(): point doesn't seem to be in quad4!");
      }
    else
      {
        bc[0] = bc[1];
        double denom = A[1]+bc[0]*N[1];
        if (fabs(denom) < prec)
          throw INTERP_KERNEL::Exception("MappedBarycentric intersection type: cuboid_mapped_coord(): point doesn't seem to be in quad4!");
        bc[1] = (P[1]-bc[0]*C[1])/denom;
        if (bc[1] < -prec || bc[1] > (1.0+prec))
          throw INTERP_KERNEL::Exception("MappedBarycentric intersection type: cuboid_mapped_coord(): point doesn't seem to be in quad4!");
      }
  }

  /*!
   * Doing as in quad_mapped_coords() would lead to a 4th order equation ... So go simpler here:
   * orthogonal distance to each pair of parallel faces is computed. The ratio gives a number in [0,1]
   *
   * Conventions:
   *   - for HEXA8, point F (5) is taken to be the origin (see med file ref connec):
   *          0 ------ 3
             /|       /|
            / |      / |
           1 ------ 2  |
           |  |     |  |
           |  |     |  |
           |  4-----|- 7
           | /      | /
           5 ------ 6

   *
   */

  inline void cuboid_mapped_coords(const std::vector<const double*>& n, const double *p, double *bc)
  {
    double prec = 1.0e-14;
    enum { _XX=0, _YY };
    if (n.size() != 8)
      throw INTERP_KERNEL::Exception("INTERP_KERNEL::cuboid_mapped_coords: unrecognized geometric type! Only HEXA8 supported.");

    double dx1, dx2, dy1, dy2, dz1, dz2;
    dx1 = OrthoDistanceFromPtToPlaneInSpaceDim3(p, n[4],n[5],n[1]);
    dx2 = OrthoDistanceFromPtToPlaneInSpaceDim3(p, n[7],n[3],n[2]);

    dy1 = OrthoDistanceFromPtToPlaneInSpaceDim3(p, n[5],n[6],n[2]);
    dy2 = OrthoDistanceFromPtToPlaneInSpaceDim3(p, n[4],n[0],n[3]);

    dz1 = OrthoDistanceFromPtToPlaneInSpaceDim3(p, n[5],n[4],n[7]);
    dz2 = OrthoDistanceFromPtToPlaneInSpaceDim3(p, n[1],n[2],n[3]);

    if (dx1 < -prec || dx2 < -prec || dy1 < -prec || dy2 < -prec || dz1 < -prec || dz2 < -prec)
      throw INTERP_KERNEL::Exception("INTERP_KERNEL::cuboid_mapped_coords: point outside HEXA8");

    bc[0] = dx1+dx2 < prec ? 0.5 : dx1/(dx1+dx2);
    bc[1] = dy1+dy2 < prec ? 0.5 : dy1/(dy1+dy2);
    bc[2] = dz1+dz2 < prec ? 0.5 : dz1/(dz1+dz2);
  }

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */
  /*         calcul la surface d'un polygone.                 */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */

  inline double  Surf_Poly(const std::vector<double>& Poly)
  { 

    double Surface=0;
    for(unsigned long i=0; i<(Poly.size())/2-2; i++)
      {
        double Surf=Surf_Tri( &Poly[0],&Poly[2*(i+1)],&Poly[2*(i+2)] ); 
        Surface=Surface + Surf ;
      }
    return Surface ;
  }

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */
  /*   fonction qui teste si un point est dans une maille     */
  /*   point: P_0                                             */
  /*   P_1, P_2, P_3 sommet des mailles                       */   
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */

  inline bool point_dans_triangle(const double* P_0,const double* P_1,
                                  const double* P_2,const double* P_3,
                                  double eps)
  {

    bool A=false;
    double det_1=mon_determinant(P_1,P_3,P_0);
    double det_2=mon_determinant(P_3,P_2,P_0);
    double det_3=mon_determinant(P_2,P_1,P_0);
    if( (det_1>=-eps && det_2>=-eps && det_3>=-eps) || (det_1<=eps && det_2<=eps && det_3<=eps) )
      {
        A=true;
      }

    return A;
  }

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */
  /*fonction pour verifier qu'un point n'a pas deja ete considerer dans   */ 
  /*      le vecteur et le rajouter au vecteur sinon.                     */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */

  inline void verif_point_dans_vect(const double* P, std::vector<double>& V, double absolute_precision )
  {
    long taille=V.size();
    bool isPresent=false;
    for(long i=0;i<taille/2;i++) 
      {
        if (sqrt(((P[0]-V[2*i])*(P[0]-V[2*i])+(P[1]-V[2*i+1])*(P[1]-V[2*i+1])))<absolute_precision)
          isPresent=true;
      
      }
    if(!isPresent)
      {
      
        V.push_back(P[0]);
        V.push_back(P[1]);    
      }
  }

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */
  /* fonction qui rajoute les sommet du triangle P dans le vecteur V        */ 
  /* si ceux-ci sont compris dans le triangle S et ne sont pas deja dans    */
  /* V.                                                                     */
  /*sommets de P: P_1, P_2, P_3                                             */
  /*sommets de S: P_4, P_5, P_6                                             */  
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */ 

  inline void rajou_sommet_triangl(const double* P_1,const double* P_2,const double* P_3,
                                   const double* P_4,const double* P_5,const double* P_6,
                                   std::vector<double>& V, double dim_caracteristic, double precision)
  {

    double absolute_precision = precision*dim_caracteristic;
    bool A_1=INTERP_KERNEL::point_dans_triangle(P_1,P_4,P_5,P_6,absolute_precision);
    if(A_1)
      verif_point_dans_vect(P_1,V,absolute_precision);
    bool A_2=INTERP_KERNEL::point_dans_triangle(P_2,P_4,P_5,P_6,absolute_precision);
    if(A_2)
      verif_point_dans_vect(P_2,V,absolute_precision);
    bool A_3=INTERP_KERNEL::point_dans_triangle(P_3,P_4,P_5,P_6,absolute_precision);
    if(A_3)
      verif_point_dans_vect(P_3,V,absolute_precision);
  }


  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  _ _ _ _ _ _ _ _*/
  /*  calcul de l'intersection de deux segments: segments P1P2 avec P3P4      */
  /*  . Si l'intersection est non nulle et si celle-ci n'est                  */
  /*  n'est pas deja contenue dans Vect on la rajoute a Vect                  */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  _ _ _ _ _ _ _ _*/ 

  inline void inters_de_segment(const double * P_1,const double * P_2,
                                const double * P_3,const double * P_4,
                                std::vector<double>& Vect, 
                                double dim_caracteristic, double precision)
  {
    // calcul du determinant de P_1P_2 et P_3P_4.
    double det=(P_2[0]-P_1[0])*(P_4[1]-P_3[1])-(P_4[0]-P_3[0])*(P_2[1]-P_1[1]);

    double absolute_precision = dim_caracteristic*precision;
    if(fabs(det)>absolute_precision)
      {
        double k_1=-((P_3[1]-P_4[1])*(P_3[0]-P_1[0])+(P_4[0]-P_3[0])*(P_3[1]-P_1[1]))/det;

        if (k_1 >= -absolute_precision && k_1 <= 1+absolute_precision)
          //if( k_1 >= -precision && k_1 <= 1+precision)
          {
            double k_2= ((P_1[1]-P_2[1])*(P_1[0]-P_3[0])+(P_2[0]-P_1[0])*(P_1[1]-P_3[1]))/det;

            if (k_2 >= -absolute_precision && k_2 <= 1+absolute_precision)
              //if( k_2 >= -precision && k_2 <= 1+precision)
              {
                double P_0[2];
                P_0[0]=P_1[0]+k_1*(P_2[0]-P_1[0]);
                P_0[1]=P_1[1]+k_1*(P_2[1]-P_1[1]);
                verif_point_dans_vect(P_0,Vect,absolute_precision);
              }
          }
      }
  }



  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/  
  /*      calcul l'intersection de deux triangles            */
  /* P_1, P_2, P_3: sommets du premier triangle              */
  /* P_4, P_5, P_6: sommets du deuxi�me triangle             */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/ 

  inline void intersec_de_triangle(const double* P_1,const double* P_2, const double* P_3,
                                   const double* P_4,const double* P_5,const double* P_6, 
                                   std::vector<double>& Vect, double dim_caracteristic, double precision)
  {
    inters_de_segment(P_1,P_2,P_4,P_5,Vect, dim_caracteristic, precision);
    inters_de_segment(P_1,P_2,P_5,P_6,Vect, dim_caracteristic, precision);
    inters_de_segment(P_1,P_2,P_6,P_4,Vect, dim_caracteristic, precision);
    inters_de_segment(P_2,P_3,P_4,P_5,Vect, dim_caracteristic, precision);
    inters_de_segment(P_2,P_3,P_5,P_6,Vect, dim_caracteristic, precision);
    inters_de_segment(P_2,P_3,P_6,P_4,Vect, dim_caracteristic, precision);
    inters_de_segment(P_3,P_1,P_4,P_5,Vect, dim_caracteristic, precision);
    inters_de_segment(P_3,P_1,P_5,P_6,Vect, dim_caracteristic, precision);
    inters_de_segment(P_3,P_1,P_6,P_4,Vect, dim_caracteristic, precision);
    rajou_sommet_triangl(P_1,P_2,P_3,P_4,P_5,P_6,Vect, dim_caracteristic, precision);
    rajou_sommet_triangl(P_4,P_5,P_6,P_1,P_2,P_3,Vect, dim_caracteristic, precision);
  }

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  /* fonction pour verifier qu'un node maille n'a pas deja ete considerer  */
  /*  dans le vecteur et le rajouter au vecteur sinon.                     */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/

  inline void verif_maill_dans_vect(int Num, std::vector<int>& V)
  {
    long taille=V.size();
    int A=0;
    for(long i=0;i<taille;i++)
      {
        if(Num==V[i])
          {
            A=1;
            break;
          } 
      }
    if(A==0)
      {V.push_back(Num); }
  }

  /*! Function that compares two angles from the values of the pairs (sin,cos)*/
  /*! Angles are considered in [0, 2Pi] bt are not computed explicitely */
  class AngleLess
  {
  public:
    bool operator()(std::pair<double,double>theta1, std::pair<double,double> theta2) const
    {
      double norm1 = sqrt(theta1.first*theta1.first +theta1.second*theta1.second);
      double norm2 = sqrt(theta2.first*theta2.first +theta2.second*theta2.second);
      
      double epsilon = 1.e-12;
      
      if( norm1 < epsilon || norm2 < epsilon  ) 
        std::cout << "Warning InterpolationUtils.hxx: AngleLess : Vector with zero norm, cannot define the angle !!!! " << std::endl;
      
      return theta1.second*(norm2 + theta2.first) < theta2.second*(norm1 + theta1.first);
    
    }
  };


  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */  
  /* fonction pour reconstituer un polygone convexe a partir  */
  /*              d'un nuage de point.                        */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */  

  inline std::vector<double> reconstruct_polygon(const std::vector<double>& V)
  {

    int taille((int)V.size());

    //VB : why 6 ?

    if(taille<=6)
      {return V;}
    else
      {
        double *COS=new double[taille/2];
        double *SIN=new double[taille/2];
        //double *angle=new double[taille/2];
        std::vector<double> Bary=bary_poly(V);
        COS[0]=1.0;
        SIN[0]=0.0;
        //angle[0]=0.0;
        for(int i=0; i<taille/2-1;i++)
          {
            std::vector<double> Trigo=calcul_cos_et_sin(&Bary[0],&V[0],&V[2*(i+1)]);
            COS[i+1]=Trigo[0];
            SIN[i+1]=Trigo[1];
            //if(SIN[i+1]>=0)
            //    {angle[i+1]=atan2(SIN[i+1],COS[i+1]);}
            //             else
            //               {angle[i+1]=-atan2(SIN[i+1],COS[i+1]);}
          }
                     
        //ensuite on ordonne les angles.
        std::vector<double> Pt_ordonne;
        Pt_ordonne.reserve(taille);
        //        std::multimap<double,int> Ordre;
        std::multimap<std::pair<double,double>,int, AngleLess> CosSin;
        for(int i=0;i<taille/2;i++)       
          {
            //  Ordre.insert(std::make_pair(angle[i],i));
            CosSin.insert(std::make_pair(std::make_pair(SIN[i],COS[i]),i));
          }
        //        std::multimap <double,int>::iterator mi;
        std::multimap<std::pair<double,double>,int, AngleLess>::iterator   micossin;
        //         for(mi=Ordre.begin();mi!=Ordre.end();mi++)
        //           {
        //             int j=(*mi).second;
        //             Pt_ordonne.push_back(V[2*j]);
        //             Pt_ordonne.push_back(V[2*j+1]);
        //           }
        for(micossin=CosSin.begin();micossin!=CosSin.end();micossin++)
          {
            int j=(*micossin).second;
            Pt_ordonne.push_back(V[2*j]);
            Pt_ordonne.push_back(V[2*j+1]);
          }
        delete [] COS;
        delete [] SIN;
        //        delete [] angle;
        return Pt_ordonne;
      }
  }

  template<int DIM, NumberingPolicy numPol, class MyMeshType>
  inline void getElemBB(double* bb, const double *coordsOfMesh, int iP, int nb_nodes)
  {
    bb[0]=std::numeric_limits<double>::max();
    bb[1]=-std::numeric_limits<double>::max();
    bb[2]=std::numeric_limits<double>::max();
    bb[3]=-std::numeric_limits<double>::max();
    bb[4]=std::numeric_limits<double>::max();
    bb[5]=-std::numeric_limits<double>::max();
    
    for (int i=0; i<nb_nodes; i++)
      {
        double x = coordsOfMesh[3*(iP+i)];
        double y = coordsOfMesh[3*(iP+i)+1];
        double z = coordsOfMesh[3*(iP+i)+2];
        bb[0]=(x<bb[0])?x:bb[0];
        bb[1]=(x>bb[1])?x:bb[1];
        bb[2]=(y<bb[2])?y:bb[2];
        bb[3]=(y>bb[3])?y:bb[3];
        bb[4]=(z<bb[4])?z:bb[4];
        bb[5]=(z>bb[5])?z:bb[5];
      }              
  }

  /*!
   * Find a vector orthogonal to the input vector
   */
  inline void orthogonalVect3(const double inpVect[3], double outVect[3])
  {
    std::vector<bool> sw(3,false);
    double inpVect2[3];
    std::transform(inpVect,inpVect+3,inpVect2,std::ptr_fun<double,double>(fabs));
    std::size_t posMin(std::distance(inpVect2,std::min_element(inpVect2,inpVect2+3)));
    sw[posMin]=true;
    std::size_t posMax(std::distance(inpVect2,std::max_element(inpVect2,inpVect2+3)));
    if(posMax==posMin)
      { posMax=(posMin+1)%3; }
    sw[posMax]=true;
    std::size_t posMid(std::distance(sw.begin(),std::find(sw.begin(),sw.end(),false)));
    outVect[posMin]=0.; outVect[posMid]=1.; outVect[posMax]=-inpVect[posMid]/inpVect[posMax];
  }
  
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  /* Computes the dot product of a and b */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  template<int dim> 
  inline double dotprod( const double * a, const double * b)
  {
    double result=0;
    for(int idim = 0; idim < dim ; idim++) result += a[idim]*b[idim];
    return result;
  }
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  /* Computes the norm of vector v */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/  
  template<int dim> 
  inline double norm(const double * v)
  {   
    double result =0;
    for(int idim =0; idim<dim; idim++) result+=v[idim]*v[idim];
    return sqrt(result);
  }
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  /* Computes the square norm of vector a-b */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/  
  template<int dim> 
  inline double distance2( const double * a, const double * b)
  {   
    double result =0;
    for(int idim =0; idim<dim; idim++) result+=(a[idim]-b[idim])*(a[idim]-b[idim]);
    return result;
  }
  template<class T, int dim> 
  inline double distance2(  T * a, int inda, T * b, int indb)
  {   
    double result =0;
    for(int idim =0; idim<dim; idim++) result += ((*a)[inda+idim] - (*b)[indb+idim])* ((*a)[inda+idim] - (*b)[indb+idim]);
    return result;
  }
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  /* Computes the determinant of a and b */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  inline double determinant ( double *  a, double * b)
  {
    return a[0]*b[1]-a[1]*b[0];
  }
  inline double determinant ( double *  a, double * b, double * c)
  {
    return a[0]*determinant(b+1,c+1)-b[0]*determinant(a+1,c+1)+c[0]*determinant(a+1,b+1);
  }
  
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  /* Computes the cross product of AB and AC */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/  

  template<int dim> inline void crossprod( const double * A, const double * B, const double * C, double * V);
  
  template<> inline
  void crossprod<2>( const double * A, const double * B, const double * C, double * V)
  {   
    double AB[2];
    double AC[2];
    for(int idim =0; idim<2; idim++) AB[idim] = B[idim]-A[idim];//B-A
    for(int idim =0; idim<2; idim++) AC[idim] = C[idim]-A[idim];//C-A;

    V[0]=determinant(AB,AC);
    V[1]=0;
  }
  template<> inline
  void crossprod<3>( const double * A, const double * B, const double * C, double * V)
  {   
    double AB[3];
    double AC[3];
    for(int idim =0; idim<3; idim++) AB[idim] = B[idim]-A[idim];//B-A
    for(int idim =0; idim<3; idim++) AC[idim] = C[idim]-A[idim];//C-A;

    V[0]=AB[1]*AC[2]-AB[2]*AC[1];
    V[1]=-AB[0]*AC[2]+AB[2]*AC[0];
    V[2]=AB[0]*AC[1]-AB[1]*AC[0];    
  }
  template<> inline
  void crossprod<1>( const double * /*A*/, const double * /*B*/, const double * /*C*/, double * /*V*/)
  {
    // just to be able to compile
  }
  
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  /* Checks wether point A is inside the quadrangle BCDE */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/  

  template<int dim> inline double check_inside(const double* A,const double* B,const double* C,const double* D,
                                               const double* E,double* ABC, double* ADE)
  {
    crossprod<dim>(A,B,C,ABC);
    crossprod<dim>(A,D,E,ADE);
    return dotprod<dim>(ABC,ADE);
  }   

  
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  /* Computes the geometric angle (in [0,Pi]) between two non zero vectors AB and AC */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/  
  template<int dim> inline double angle(const double * A, const double * B, const double * C, double * n)
  {   
    double AB[dim];
    double AC[dim];
    double orthAB[dim];

    for(int idim =0; idim<dim; idim++) AB[idim] = B[idim]-A[idim];//B-A;
    for(int idim =0; idim<dim; idim++) AC[idim] = C[idim]-A[idim];//C-A;

    double normAB= norm<dim>(AB);
    for(int idim =0; idim<dim; idim++) AB[idim]/=normAB;

    double normAC= norm<dim>(AC);
    double AB_dot_AC=dotprod<dim>(AB,AC);
    for(int idim =0; idim<dim; idim++) orthAB[idim] = AC[idim]-AB_dot_AC*AB[idim];

    double denom= normAC+AB_dot_AC;
    double numer=norm<dim>(orthAB);
    
    return 2*atan2(numer,denom);
  }    
  
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  /* Tells whether the frame constituted of vectors AB, AC and n is direct */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/  
  template<int dim> inline double direct_frame(const double * A, const double * B, const double * C, double * n);
  template<> inline
  double direct_frame<2>(const double * A, const double * B, const double * C, double * n)
  {
    double AB[2];
    double AC[2];
    for(int idim =0; idim<2; idim++) AB[idim] = B[idim]-A[idim];//B-A;
    for(int idim =0; idim<2; idim++) AC[idim] = C[idim]-A[idim];//C-A;
    
    return  determinant(AB,AC)*n[0];
  }
  template<> inline
  double direct_frame<3>(const double * A, const double * B, const double * C, double * n)
  {
    double AB[3];
    double AC[3];
    for(int idim =0; idim<3; idim++) AB[idim] = B[idim]-A[idim];//B-A;
    for(int idim =0; idim<3; idim++) AC[idim] = C[idim]-A[idim];//C-A;
    
    return determinant(AB,AC,n)>0;
  }

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/  
  /*      calcul l'intersection de deux polygones COPLANAIRES */
  /* en dimension DIM (2 ou 3). Si DIM=3 l'algorithme ne considere*/
  /* que les deux premieres coordonnees de chaque point */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/ 
  template<int DIM> inline void intersec_de_polygone(const double * Coords_A, const double * Coords_B, 
                                                     int nb_NodesA, int nb_NodesB,
                                                     std::vector<double>& inter,
                                                     double dim_caracteristic, double precision)
  {
    for(int i_A = 1; i_A<nb_NodesA-1; i_A++)
      {
        for(int i_B = 1; i_B<nb_NodesB-1; i_B++)
          {
            INTERP_KERNEL::intersec_de_triangle(&Coords_A[0],&Coords_A[DIM*i_A],&Coords_A[DIM*(i_A+1)],
                                                &Coords_B[0],&Coords_B[DIM*i_B],&Coords_B[DIM*(i_B+1)],
                                                inter, dim_caracteristic, precision);
          }
      }
    int nb_inter=((int)inter.size())/DIM;
    if(nb_inter >3) inter=INTERP_KERNEL::reconstruct_polygon(inter);
  }

  /*_ _ _ _ _ _ _ _ _
   *_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
   *  fonctions qui calcule l'aire d'un polygone en dimension 2 ou 3    
   *_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */
  template<int DIM> inline double polygon_area(std::vector<double>& inter)
  {
    double result=0.;
    double area[DIM];
                  
    for(int i = 1; i<(int)inter.size()/DIM-1; i++)
      {
        INTERP_KERNEL::crossprod<DIM>(&inter[0],&inter[DIM*i],&inter[DIM*(i+1)],area);
        result +=0.5*norm<DIM>(area);
      }
    return result;
  }
         
  template<int DIM> inline double polygon_area(std::deque<double>& inter)
  {
    double result=0.;
    double area[DIM];
                  
    for(int i = 1; i<(int)inter.size()/DIM-1; i++)
      {
        INTERP_KERNEL::crossprod<DIM>(&inter[0],&inter[DIM*i],&inter[DIM*(i+1)],area);
        result +=0.5*norm<DIM>(area);
      }
    return result;
  }
  
  /*! Computes the triple product (XA^XB).XC (in 3D)*/
  inline double triple_product(const double* A, const double*B, const double*C, const double*X)
  {
    double XA[3];
    XA[0]=A[0]-X[0];
    XA[1]=A[1]-X[1];
    XA[2]=A[2]-X[2];
    double XB[3];
    XB[0]=B[0]-X[0];
    XB[1]=B[1]-X[1];
    XB[2]=B[2]-X[2];
    double XC[3];
    XC[0]=C[0]-X[0];
    XC[1]=C[1]-X[1];
    XC[2]=C[2]-X[2];
    
    return 
      (XA[1]*XB[2]-XA[2]*XB[1])*XC[0]+
      (XA[2]*XB[0]-XA[0]*XB[2])*XC[1]+
      (XA[0]*XB[1]-XA[1]*XB[0])*XC[2];
  }
  
  /*! Subroutine of checkEqualPolygins that tests if two list of nodes (not necessarily distincts) describe the same polygon, assuming they share a comon point.*/
  /*! Indexes istart1 and istart2 designate two points P1 in L1 and P2 in L2 that have identical coordinates. Generally called with istart1=0.*/
  /*! Integer sign ( 1 or -1) indicate the direction used in going all over L2. */
  template<class T, int dim> 
  bool checkEqualPolygonsOneDirection(T * L1, T * L2, int size1, int size2, int istart1, int istart2, double epsilon, int sign)
  {
    int i1 = istart1;
    int i2 = istart2;
    int i1next = ( i1 + 1 ) % size1;
    int i2next = ( i2 + sign +size2) % size2;
    
    while(true)
      {
        while( i1next!=istart1 && distance2<T,dim>(L1,i1*dim, L1,i1next*dim) < epsilon ) i1next = (  i1next + 1 ) % size1;  
        while( i2next!=istart2 && distance2<T,dim>(L2,i2*dim, L2,i2next*dim) < epsilon ) i2next = (  i2next + sign +size2 ) % size2;  
        
        if(i1next == istart1)
          {
            if(i2next == istart2)
              return true;
            else return false;
          }
        else
          if(i2next == istart2)
            return false;
          else 
            {
              if(distance2<T,dim>(L1,i1next*dim, L2,i2next*dim) > epsilon )
                return false;
              else
                {
                  i1 = i1next;
                  i2 = i2next;
                  i1next = ( i1 + 1 ) % size1;
                  i2next = ( i2 + sign + size2 ) % size2;
                }
            }
      }
  }

  /*! Tests if two list of nodes (not necessarily distincts) describe the same polygon.*/
  /*! Existence of multiple points in the list is considered.*/
  template<class T, int dim> 
  bool checkEqualPolygons(T * L1, T * L2, double epsilon)
  {
    if(L1==NULL || L2==NULL) 
      {
        std::cout << "Warning InterpolationUtils.hxx:checkEqualPolygonsPointer: Null pointer " << std::endl;
        throw(Exception("big error: not closed polygon..."));
      }
    
    int size1 = (*L1).size()/dim;
    int size2 = (*L2).size()/dim;
    int istart1 = 0;
    int istart2 = 0;
    
    while( istart2 < size2  && distance2<T,dim>(L1,istart1*dim, L2,istart2*dim) > epsilon ) istart2++;
  
    if(istart2 == size2)
      {  
        return (size1 == 0) && (size2 == 0);
      }
    else 
      return   checkEqualPolygonsOneDirection<T,dim>( L1, L2, size1, size2, istart1, istart2, epsilon,  1)
        || checkEqualPolygonsOneDirection<T,dim>( L1, L2, size1, size2, istart1, istart2, epsilon, -1);

  }
}


#endif
