//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#ifndef __INTERPOLATIONUTILS_HXX__
#define __INTERPOLATIONUTILS_HXX__

#include "INTERPKERNELDefines.hxx"
#include "InterpKernelException.hxx"

#include "NormalizedUnstructuredMesh.hxx"

#include <deque>
#include <map>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>

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
  /*     fonction qui calcul le d�terminant            */
  /*      de deux vecteur(cf doc CGAL).                */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/

  //fonction qui calcul le d�terminant des vecteurs: P3P1 et P3P2
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

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */
  /*     calcul les coordonn�es du barycentre d'un polygone   */ 
  /*     le vecteur en entr�e est constitu� des coordonn�es   */
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
    double A=2*x/taille;
    double B=2*y/taille;
    Bary.push_back(A);//taille vecteur=2*nb de points.
    Bary.push_back(B);


    return Bary;
  }


  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  /*     Calculate barycentric coordinates of a 2D point p */ 
  /*     with respect to the triangle verices.             */
  /*     triaCoords are in full interlace                  */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/

  inline void barycentric_coords(const double* triaCoords, const double* p, double* bc)
  {
    // matrix 2x2
    double
      T11 = triaCoords[0]-triaCoords[4], T12 = triaCoords[2]-triaCoords[4],
      T21 = triaCoords[1]-triaCoords[5], T22 = triaCoords[3]-triaCoords[5];
    // matrix determinant
    double Tdet = T11*T22 - T12*T21;
    if ( fabs( Tdet ) < std::numeric_limits<double>::min() ) {
      bc[0]=1; bc[1]=0; bc[2]=0;
      return;
    }
    // matrix inverse
    double t11 = T22, t12 = -T12, t21 = -T21, t22 = T11;
    // vector
    double r11 = p[0]-triaCoords[4], r12 = p[1]-triaCoords[5];
    // barycentric coordinates: mutiply matrix by vector
    bc[0] = (t11 * r11 + t12 * r12)/Tdet;
    bc[1] = (t21 * r11 + t22 * r12)/Tdet;
    bc[2] = 1. - bc[0] - bc[1];
  }

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  /*     Calculate barycentric coordinates of a point p    */ 
  /*     with respect to triangle or tetra verices.        */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/

  inline void barycentric_coords(const std::vector<const double*>& n, const double* p, double* bc)
  {
    enum { _X, _Y, _Z };
    if ( n.size() == 3 ) // TRIA3
      {
        // matrix 2x2
        double
          T11 = n[0][_X]-n[2][_X], T12 = n[1][_X]-n[2][_X],
          T21 = n[0][_Y]-n[2][_Y], T22 = n[1][_Y]-n[2][_Y];
        // matrix determinant
        double Tdet = T11*T22 - T12*T21;
        if ( std::fabs( Tdet ) < std::numeric_limits<double>::min() ) {
          bc[0]=1; bc[1]=bc[2]=0; // no solution
          return;
        }
        // matrix inverse
        double t11 = T22, t12 = -T12, t21 = -T21, t22 = T11;
        // vector
        double r11 = p[_X]-n[2][_X], r12 = p[_Y]-n[2][_Y];
        // barycentric coordinates: mutiply matrix by vector
        bc[0] = (t11 * r11 + t12 * r12)/Tdet;
        bc[1] = (t21 * r11 + t22 * r12)/Tdet;
        bc[2] = 1. - bc[0] - bc[1];
      }
    else // TETRA4
      {
        bc[3]=0; // for no solution

        // Find bc by solving system of 3 equations using Gaussian elimination algorithm
        // bc1*( x1 - x4 ) + bc2*( x2 - x4 ) + bc3*( x3 - x4 ) = px - x4
        // bc1*( y1 - y4 ) + bc2*( y2 - y4 ) + bc3*( y3 - y4 ) = px - y4
        // bc1*( z1 - z4 ) + bc2*( z2 - z4 ) + bc3*( z3 - z4 ) = px - z4
        const int nbCol=4, nbRow=3;

        double T[nbRow][nbCol]=
          {{ n[0][_X]-n[3][_X], n[1][_X]-n[3][_X], n[2][_X]-n[3][_X], p[_X]-n[3][_X] },
           { n[0][_Y]-n[3][_Y], n[1][_Y]-n[3][_Y], n[2][_Y]-n[3][_Y], p[_Y]-n[3][_Y] },
           { n[0][_Z]-n[3][_Z], n[1][_Z]-n[3][_Z], n[2][_Z]-n[3][_Z], p[_Z]-n[3][_Z] }};

        // make upper triangular matrix (forward elimination)

        int iR[nbRow] = { 0, 1, 2 };

        for ( int i = 0; i < 2; ++i ) // nullify 2 rows
          {
            // swap rows to have max value of i-th column in i-th row
            double max = std::fabs( T[ iR[i] ][i] );
            for ( int r = i+1; r < nbRow; ++r ) {
              double t = std::fabs( T[ iR[r] ][i] );
              if ( t > max ) {
                max = t;
                std::swap( iR[r], iR[i] );
              }
            }
            if ( max < std::numeric_limits<double>::min() ) {
              bc[0]=1; bc[1]=bc[2]=bc[3]=0;
              return; // no solution
            }
            // make 0 below T[i][i] (actually we do not modify i-th column)
            double* tUpRow = T[ iR[i] ];
            for ( int r = i+1; r < nbRow; ++r ) {
              double* tRow = T[ iR[r] ];
              double coef = tRow[ i ] / tUpRow[ i ];
              for ( int c = i+1; c < nbCol; ++c )
                tRow[ c ] -= tUpRow[ c ] * coef;
            }
          }
        double* tRow = T[ iR[2] ];
        if ( std::fabs( tRow[ 2 ] ) < std::numeric_limits<double>::min() ) {
          bc[0]=1; bc[1]=bc[2]=bc[3]=0;
          return; // no solution
        }
        tRow[ 3 ] /= tRow[ 2 ];

        // calculate solution (back substitution)

        bc[ 2 ] = tRow[ 3 ];

        tRow = T[ iR[1] ];
        bc[ 1 ] = (tRow[ 3 ] - bc[2]*tRow[ 2 ]) / tRow[ 1 ];

        tRow = T[ iR[0] ];
        bc[ 0 ] = (tRow[ 3 ] - bc[2]*tRow[ 2 ] - bc[1]*tRow[ 1 ]) / tRow[ 0 ];

        bc[ 3 ] = 1. - bc[0] - bc[1] - bc[2];
      }
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
  /*fonction pour v�rifier qu'un point n'a pas d�ja �t� consid�rer dans   */ 
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
  /* si ceux-ci sont compris dans le triangle S et ne sont pas d�j� dans    */
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
  /*  n'est pas d�j� contenue dans Vect on la rajoute � Vect                  */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  _ _ _ _ _ _ _ _*/ 

  inline void inters_de_segment(const double * P_1,const double * P_2,
                                const double * P_3,const double * P_4,
                                std::vector<double>& Vect, 
                                double dim_caracteristic, double precision)
  {
    // calcul du d�terminant de P_1P_2 et P_3P_4.
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
  /* fonction pour v�rifier qu'un n�de maille n'a pas d�ja �t� consid�rer  */
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
    bool operator()(std::pair<double,double>theta1, std::pair<double,double> theta2) 
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
  /* fonction pour reconstituer un polygone convexe � partir  */
  /*              d'un nuage de point.                        */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */  

  inline std::vector<double> reconstruct_polygon(const std::vector<double>& V)
  {

    int taille=V.size();

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

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  /* Computes the dot product of a and b */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  template<int dim> 
  inline double dotprod( double * a, double * b)
  {
    double result=0;
    for(int idim = 0; idim < dim ; idim++) result += a[idim]*b[idim];
    return result;
  }
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/
  /* Computes the norm of vector v */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/  
  template<int dim> 
  inline double norm( double * v)
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
  /* en dimension DIM (2 ou 3). Si DIM=3 l'algorithme ne consid�re*/
  /* que les deux premi�res coordonn�es de chaque point */
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
