//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D, OPEN CASCADE
//
//  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
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
#ifndef _INTERPOLATIONUTILS_HXX_
#define _INTERPOLATIONUTILS_HXX_

#include <INTERPKERNEL_defines.hxx>

#include "NormalizedUnstructuredMesh.hxx"

#include <deque>
#include <map>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

namespace INTERP_KERNEL
{

  class Exception : std::exception
  {
  public:
    Exception(const char *what):_reason(what) { }
    ~Exception() throw () { }
    const char *what() { return _reason.c_str(); }
  protected:
    std::string _reason;
  };

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
  /*     fonction qui calcul le déterminant            */
  /*      de deux vecteur(cf doc CGAL).                */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _*/

  //fonction qui calcul le déterminant des vecteurs: P3P1 et P3P2
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

  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ */
  /*     calcul les coordonnées du barycentre d'un polygone   */ 
  /*     le vecteur en entrée est constitué des coordonnées   */
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
  /*fonction pour vérifier qu'un point n'a pas déja été considérer dans   */ 
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
  /* si ceux-ci sont compris dans le triangle S et ne sont pas déjà dans    */
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
  /*  n'est pas déjà contenue dans Vect on la rajoute à Vect                  */
  /*_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _  _ _ _ _ _ _ _ _*/ 

  inline void inters_de_segment(const double * P_1,const double * P_2,
                                const double * P_3,const double * P_4,
                                std::vector<double>& Vect, 
                                double dim_caracteristic, double precision)
  {
    // calcul du déterminant de P_1P_2 et P_3P_4.
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
  /* P_4, P_5, P_6: sommets du deuxième triangle             */
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
  /* fonction pour vérifier qu'un n°de maille n'a pas déja été considérer  */
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
  /* fonction pour reconstituer un polygone convexe à partir  */
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
						//	Ordre.insert(std::make_pair(angle[i],i));
						CosSin.insert(std::make_pair(std::make_pair(SIN[i],COS[i]),i));
					}
				//        std::multimap <double,int>::iterator mi;
				std::multimap<std::pair<double,double>,int, AngleLess>::iterator   micossin;
// 				for(mi=Ordre.begin();mi!=Ordre.end();mi++)
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
  /* en dimension DIM (2 ou 3). Si DIM=3 l'algorithme ne considère*/
  /* que les deux premières coordonnées de chaque point */
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
				while( i1next!=istart1 && distance2<T,dim>(L1,i1*dim, L1,i1next*dim) < epsilon ) i1next = (	i1next + 1 ) % size1;	
				while( i2next!=istart2 && distance2<T,dim>(L2,i2*dim, L2,i2next*dim) < epsilon ) i2next = (	i2next + sign +size2 ) % size2;	
				
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
