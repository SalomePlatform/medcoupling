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

#ifndef __POLYGONALGORITHMS_HXX__
#define __POLYGONALGORITHMS_HXX__

#include <vector>
#include <deque>
#include <map>

namespace INTERP_KERNEL
{
  template<int DIM>
  class VertexLess
  {
  public:
    bool operator()(const double * P_1, const double * P_2) const
    {
      for(int idim=0; idim<DIM; idim++)
        {        
          if(P_1[idim] < P_2[idim] )  return true;
          else if( P_1[idim] > P_2[idim]) return false;
        }
      return false;
    }
  };
  
  template<int DIM>
  class PolygonAlgorithms
  {
  public:
    PolygonAlgorithms(double epsilon, double precision);
    std::deque<double> intersectConvexPolygons(const double* P_1,const double* P_2, int N1, int N2);

    //Not yet tested
    int convexDecomposition(const double * P, int N, std::vector< std::map< int,int > >& components,
                            std::vector< int >& components_index, const double epsilon);
  private:
    void defineIndices(int& i_loc, int& i_next, int& i_prev, 
                       const double *& Poly1, const double *& Poly2,
                       int& j1, int& j1_glob, int& j2, int& j2_glob,
                       int& j3, int& j3_glob, int& j4, int& j4_glob, 
                       int& i_glob, int& i_next_glob, int& i_prev_glob, 
                       const double * P_1, const double * P_2, 
                       int N1, int N2, int sign);
    void addCrossings( const double * A, const double * B, int i , int i_next,
                       const double * C, const double * D, int j1, int j2,
                       const double * E, const double * F, int j3, int j4,
                       const double * G);
    void addCrossing0(const double * A, const double * B, int i, int i_next,
                      const double * C, const double * D, int j, int j_next);
    void addCrossing( double * ABCD, std::pair< int,int > i_i_next, std::pair< int,int > j_j_next);
    void addNewVertex( int i, int i_glob, int i_next_glob, int i_prev_glob, const double * P);
    bool intersectSegmentSegment(const double * A,  const double * B, const double * C,
                                 const double * D,  const double * E, double * V);


    //Not yet tested
    void convexDecomposition(const double* P, int N, double* n,  std::vector< int > subP, int NsubP, 
                             std::vector< std::map< int,int > >& components, std::vector< int >& components_index,
                             int& Ncomp, int sign, const double epsilon);
    void convHull(const double *P, int N, double * n,  std::map< int,int >& subP,
                  std::map< int,int >& not_in_hull, int& NsubP, const double epsilon);
  private:
    std::deque< double > _Inter;/* vertices of the intersection  P1^P2 */
    std::vector< std::pair< int,int > > _End_segments; /* segments containing inter final edges */   
    /* status list of segments (ending point, starting point) intersected by the sweeping line */
    /* and a boolean true if the ending point is in the intersection */
    std::multimap< int, std::pair< int,bool> > _Status;
    bool _is_in_intersection;
    bool _terminus;
    double _vdouble[DIM];
    double _epsilon;
    double _precision;
  };
}

#endif
