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
#ifndef __POLYGONALGORITHMS_TXX__
#define __POLYGONALGORITHMS_TXX__

#include "PolygonAlgorithms.hxx"
#include "InterpolationUtils.hxx"
#include <list>
#include <map>
#include <iostream>

namespace INTERP_KERNEL
{
  template<int DIM>
  PolygonAlgorithms<DIM>::PolygonAlgorithms(double epsilon, double precision)//: (0)
  {
    _is_in_intersection = false;
    _epsilon = epsilon;
    _precision = precision;
  }
  /*************************************************************/
  /* Computes the 3D intersection between two COPLANAR         */
  /* Segments [A,B] and [C,D], stores the result in V.         */
  /* If A belongs to [CD] then the vertex E (preceeding A)     */
  /* is used to decide if the crossing is real. If A coincides */
  /* with C or D, a special treatment is performed             */
  /*************************************************************/
  template<int DIM>
  bool PolygonAlgorithms<DIM>::intersectSegmentSegment(const double * A,  const double * B, const double * C,
                                                       const double * D, const double * E, double * V)
  {    
    double AB[DIM], DC[DIM], AC[DIM], det, t1, t2, inv_det;
    
    /******* Initialisation of the linear system  t1*AB+t2*DC=AC ***********/
    for(int idim=0;idim<DIM;idim++)
      {
        AB[idim] = B[idim]-A[idim];//B-A
        DC[idim] = C[idim]-D[idim];//C-D
        AC[idim] = C[idim]-A[idim];//C-A
      }
    
    /******* Resolution of the linear system  t1*AB+t2*DC=AC ***********/    
    det = determinant(AB,DC);//determinant of the first two coordinates
    if(fabs(det) >_epsilon)
      {   
        inv_det = 1/det;
        t1 = determinant(AC,DC)*inv_det;//solves the linear system t1*AB+t2*DC=AC
        t2 = determinant(AB,AC)*inv_det;//solves the linear system t1*AB+t2*DC=AC
      }
    else 
      {
        switch(DIM)
          {
          case 2:
            {
              if(distance2<DIM>(A,D)<_epsilon)
                crossprod<DIM>(A,C,E,_vdouble);//store the crossprod between vectors AC and AE (E=vertex preceding A)                     
              return false;//case of paralell segments
            }
          case 3://beware AB and CD may belong to a vertical plane
            det = determinant(&AB[1],&DC[1]);//determinant of the last two coefficients
            if(fabs(det) > _epsilon) 
              {
                inv_det = 1/det;
                t1=(AC[1]*DC[DIM-1]-AC[DIM-1]*DC[1])*inv_det;
                t2=(AB[1]*AC[DIM-1]-AB[DIM-1]*AC[1])*inv_det;
              }
            else //beware AB and CD may belong to a plane y = constant
              {
                det = AB[0]*DC[DIM-1]-AB[DIM-1]*DC[0];
                if(fabs(det) > _epsilon) 
                  {
                    inv_det = 1/det;
                    t1=(AC[0]*DC[DIM-1]-AC[DIM-1]*DC[0])*inv_det;
                    t2=(AB[0]*AC[DIM-1]-AB[DIM-1]*AC[0])*inv_det;
                  }
                else
                  {
                    if(distance2<DIM>(A,D)<_epsilon)
                      crossprod<DIM>(A,C,E,_vdouble);//store the crossprod between vectors AC and AE (E=vertex preceding A)                     
                    return false;//case of paralell segments
                  }
              }
          }
      }
    
    if(t1>_precision && t1<1-_precision)
      {
        if( t2>_precision && t2<1-_precision)
          {         
            for(int idim=0;idim<DIM;idim++) V[idim]=t1*AB[idim] + A[idim];
            return true;
          }
      }
    else if(fabs(t1) <= _precision)
      {
        if( t2>_precision && t2<1-_precision)//vertex on an edge
          {
            double V12[DIM];
            double V34[DIM];
            crossprod<DIM>(A,D,B,V12);
            crossprod<DIM>(A,D,E,V34);
            double same_side =dotprod<DIM>(V12, V34);        
            if( same_side < -_epsilon ) // <= epsilon or 0 ?//crossing
              {
                for(int idim=0;idim<DIM;idim++) V[idim]=A[idim]; 
                return true;
              }
            else if( same_side > _epsilon ) _terminus= !_is_in_intersection;//reflexion
            else //separation of overlaping edges
              {
                if(_Inter.empty() ) _terminus=true;
                else if(!_is_in_intersection)
                  {
                    for(int idim=0;idim<DIM;idim++) V[idim]=A[idim];
                    return true;
                  }
              }         
          }
        else if(fabs(t2-1) <= _precision)//vertex on a vertex (A=D), first run
          crossprod<DIM>(A,C,E,_vdouble);//store the angle between vectors AC and AE (E=vertex preceding A)                     
        else if(fabs(t2) <= _precision)//vertex on a vertex (A=C), second run
          {
            double Vdoublebis[DIM];
            //crossprod<DIM>(A,C,E,_vdouble);
            crossprod<DIM>(A,B,D,Vdoublebis);
            double in_between =dotprod<DIM>(Vdoublebis,_vdouble);
            if(in_between>_epsilon)//crossing
              {
                for(int idim=0;idim<DIM;idim++) V[idim]=A[idim]; 
                return true;
              }
            else if(fabs(in_between)<=_epsilon && dotprod<DIM>(Vdoublebis,Vdoublebis) > _epsilon)
              //ie _vdouble=0, separation of overlaping edges at a double point
              {
                //crossprod<DIM>(A,E,B,_vdouble); 
                if(dotprod<DIM>(_vdouble,Vdoublebis) >=_epsilon )//crossing
                  {
                    if(_Inter.empty()) _terminus=true;
                    else if(!_is_in_intersection)
                      {
                        for(int idim=0;idim<DIM;idim++) V[idim]=A[idim];
                        return true;
                      }
                  }
              } 
          }
      }
    return false;
  }
  
  /*************************************************************/  
  /* adds vertex i  to the list inter and updates _End_segments */
  /* i is the local index of the current vertex                */
  /*************************************************************/ 
  template<int DIM>
  inline void PolygonAlgorithms<DIM>::addNewVertex( int i, int i_glob, int i_next_glob, int i_prev_glob,
                                                    const double * P)
  {    
    /* Question:Should we add vertex i to the front or back ? */
    if( _End_segments[1].second == i_glob)
      {
        for(int idim=0;idim<DIM;idim++) _Inter.push_back(P[DIM*i+idim]);
        _End_segments[1] = std::make_pair(i_glob, i_next_glob);
      }
    else
      {
        for(int idim=DIM-1;idim>-1;idim--) _Inter.push_front(P[DIM*i+idim]);
        _End_segments[0] = std::make_pair(i_glob, i_next_glob);
      }
  }  
  
  /************************************************************/  
  /* adds a crossing between two segments starting at i and j */
  /* to the double ended list inter in the correct order      */
  /* according to endsegments, updates _End_segments           */
  /************************************************************/ 
  template<int DIM>
  inline void PolygonAlgorithms<DIM>::addCrossing( double * ABCD, std::pair< int,int > i_i_next, 
                                                   std::pair< int,int > j_j_next)
  {    
    if(!_Inter.empty() )
      {
        if(_End_segments[0] ==i_i_next)
          {
            for(int idim=DIM-1;idim>=0;idim--) _Inter.push_front(ABCD[idim]);
            _terminus= (_End_segments[1]== j_j_next);
            _End_segments[0] = j_j_next;
          }
        else
          {
            if( _End_segments[0]== j_j_next)
              {        
                for(int idim=DIM-1;idim>=0;idim--) _Inter.push_front(ABCD[idim]);
                _terminus= (_End_segments[1]== i_i_next);
                _End_segments[0] = i_i_next;
              }
            else
              {
                for(int idim=0;idim<DIM;idim++) _Inter.push_back(ABCD[idim]);
                _End_segments[1] = (_End_segments[1]== i_i_next) ? j_j_next : i_i_next;
              }
          }
      }      
    else
      {
        for(int i=0;i<DIM;i++) _Inter.push_back(ABCD[i]);
        _End_segments.push_back(i_i_next);
        _End_segments.push_back(j_j_next);
      }
  }
  
  /*******************************************************/
  /* checks the possible crossing between segments [A,B] */
  /* (with end-point global indices i and i_next)        */
  /* and  [C,D] (end-point global indices j and j_next). */
  /* If no intersection is detected, checks whether B is */
  /* inside the quadrangle AEDC.                         */
  /* Updates the lists inter and _End_segments            */
  /*******************************************************/
 
  template<int DIM>
  void PolygonAlgorithms<DIM>::addCrossing0(const double * A, const double * B, int i, int i_next,
                                            const double * C, const double * D, int j, int j_next)
  {
    double ABCD[DIM];                
    if(intersectSegmentSegment(A,B,C,D,ABCD, ABCD))
      //fifth and sixth arguments are useless here
      {
        /* Updating _End_segments */
        std::pair< int,int > i_i_next = std::make_pair(i, i_next);
        std::pair< int,int > j_j_next = std::make_pair(j, j_next); 
        if( _End_segments[0] == i_i_next)
          {         
            for(int idim=DIM-1;idim>-1;idim--) _Inter.push_front(ABCD[idim]);
            _End_segments[0] = j_j_next;
          }
        else 
          {
            for(int idim=0;idim<DIM;idim++) _Inter.push_back(ABCD[idim]);
            _End_segments[1] = j_j_next;
            _terminus = _End_segments[0]== j_j_next;
          }
         
        /* Updating _Status */
        _Status.insert(make_pair(i_next,std::make_pair(i, false)));
        std::multimap< int, std::pair< int,bool> >::iterator mi =_Status.find(j_next);
        ((* mi).second).second= !((* mi).second).second;
      }
    else        _Status.insert(std::make_pair(i_next,std::make_pair(i,true)));
  }  

  /*******************************************************/
  /* adds the possible crossings between segments [A,B] (with end-point global indices i and i_next) */
  /*and segments [C,D] and [E,F] to the list inter and updates _End_segments */
  /* In cases of ambiguity, the vertex G is used to decide wether the crossing should be accepted */
  /*******************************************************/
  template<int DIM>
  inline void PolygonAlgorithms<DIM>::addCrossings( const double * A, const double * B, int i , int i_next,
                                                    const double * C, const double * D, int j1, int j2,
                                                    const double * E, const double * F, int j3, int j4,
                                                    const double * G)
  {
    double ABCD[DIM];
    double ABEF[DIM];
    std::multimap< int, std::pair< int,bool> >::iterator mi;
    
    if(intersectSegmentSegment(A,B,C,D,G,ABCD))
      {
        if(intersectSegmentSegment(A,B,E,F,G,ABEF))
          {
            VertexLess<DIM> vl;
            if (vl(ABCD,ABEF))
              {
                addCrossing(ABCD,  std::make_pair(i, i_next),  std::make_pair(j1, j2));
                addCrossing(ABEF,  std::make_pair(i, i_next),  std::make_pair(j3, j4));
              }
            else
              {
                addCrossing(ABEF,  std::make_pair(i, i_next),  std::make_pair(j3, j4));
                addCrossing(ABCD,  std::make_pair(i, i_next),  std::make_pair(j1, j2));
              }
            _Status.insert(std::make_pair(i_next,std::make_pair(i, _is_in_intersection)));
            mi=_Status.find(j2);
            ((* mi).second).second= !((* mi).second).second;
            mi=_Status.find(j4); 
            ((* mi).second).second= !((* mi).second).second;
          }
        else
          {
            addCrossing(ABCD, std::make_pair( i, i_next),  std::make_pair(j1,j2));
            _Status.insert(std::make_pair(i_next,std::make_pair(i, !_is_in_intersection)));
            mi=_Status.find(j2); 
            ((* mi).second).second= !((* mi).second).second;
          }
      }
    else
      {
        if(intersectSegmentSegment(A,B,E,F,G, ABEF))
          {
            addCrossing(ABEF, std::make_pair( i, i_next), std::make_pair( j3, j4));  
            _Status.insert(std::make_pair(i_next,std::make_pair(i, !_is_in_intersection)));
            mi=_Status.find(j4);
            ((* mi).second).second= !((* mi).second).second;
          }
        else           _Status.insert(std::make_pair(i_next,std::make_pair(i, _is_in_intersection)));      
      }
  }


  /* define various indices required in the function intersect_conv_polygon */
  /* vertices from the both polygons are supposed to be present in the status */
  template<int DIM>
  inline void PolygonAlgorithms<DIM>::defineIndices(int& i_loc, int& i_next, int& i_prev, 
                                                    const double *& Poly1, const double *& Poly2,
                                                    int& j1, int& j1_glob, int& j2, int& j2_glob,
                                                    int& j3, int& j3_glob, int& j4,  int& j4_glob, 
                                                    int& i_glob, int& i_next_glob, int& i_prev_glob, 
                                                    const double * P_1, const double * P_2, 
                                                    int N1, int N2, int sign)
  {
    int N0, shift;
    if(i_glob < N1)
      { 
        N0 = N1;
        shift = 0;
        Poly1 = P_1;
        Poly2 = P_2;
       
        std::multimap< int, std::pair< int,bool> >::reverse_iterator mi1=_Status.rbegin();
        j1_glob=((*mi1).second).first;
        j1=j1_glob-N1;
        j2_glob=(*mi1).first;
        j2=j2_glob-N1;
        mi1++;
        j3_glob=((*mi1).second).first;
        j3=j3_glob-N1;
        j4_glob=(*mi1).first;
        j4=j4_glob-N1;
      }
    else
      { 
        N0 = N2;
        shift = N1;
        Poly1 = P_2;
        Poly2 = P_1;
       
        std::multimap< int, std::pair< int,bool> >::iterator mi2= _Status.begin();
        j1_glob=((*mi2).second).first;
        j1=j1_glob;
        j2_glob=(*mi2).first;
        j2=j2_glob;
        mi2++;
        j3_glob=((*mi2).second).first;
        j3=j3_glob;
        j4_glob=(*mi2).first;
        j4=j4_glob;
      }
    i_loc = i_glob-shift;
    i_next = (i_next_glob-shift+N0)%N0;//end-point of segment starting at i
    i_prev = (i_prev_glob-shift+N0)%N0;
    i_next_glob = i_next+shift;
    i_prev_glob = i_prev+shift;
    //warning: sign is either 1 or -1;
    //To do: test and remove from Convex_intersecor.cxx
    //        while(distance2<DIM>(&Poly1[DIM*i_loc],&Poly1[DIM*i_next])< _epsilon && i_next != i_loc)
    //          i_next =(i_next+sign+N0)%N0; 
    //        while(distance2<DIM>(&Poly1[DIM*i_loc],&Poly1[DIM*i_prev])< _epsilon && i_prev != i_loc) 
    //          i_prev =(i_prev+sign+N0)%N0; 
  }
  /*******************************************************/
  /* computes the vertices of the intersection of two COPLANAR */
  /* simple (no dble points)convex polygons using line sweep algorithm */
  /* P1 and P2 contain the 3D coordinates of the successive vertices */
  /*******************************************************/
  template<int DIM>
  std::deque< double > PolygonAlgorithms<DIM>::intersectConvexPolygons(const double* P_1,const double* P_2,
                                                                       int N1, int N2)
  {    
    int i_loc, i_glob, j1, j1_glob, j2,j2_glob, j3, j3_glob, j4,j4_glob,
      i_prev, i_prev_glob, i_next, i_next_glob, nb_prev, sign, idim;
    const double * Poly1, * Poly2;
    bool four_neighbours=false;
    _terminus = N1 < 3 || N2<3;

    /* list of future events ordered according to their coordinates (x,y,z) (lexicographical order) */
    std::multimap< const double *, int, VertexLess<DIM> > mmap_events;
    typename std::list< std::pair< const double *, int > >::iterator mi1,mi2;

    std::multimap< int, std::pair< int,bool> >::iterator mi;

    /********** Initalisation of events with P1 and P2 vertices ************/
    for(i_loc=0;i_loc<N1;i_loc++)
      mmap_events.insert(std::make_pair(&P_1[DIM*i_loc],i_loc));
    for(i_loc=0;i_loc<N2;i_loc++)
      mmap_events.insert(std::make_pair(&P_2[DIM*i_loc],i_loc+N1));
                
    std::list< std::pair< const double *, int > > events(mmap_events.begin(),mmap_events.end());
                
    if(!_terminus)
      {
        /******** Treatment of the first vertex ********/
        mi1=events.begin();
        i_glob = (* mi1).second;
        bool which_start = i_glob < N1;
        if(i_glob < N1){ i_next_glob = (i_glob   +1)%N1;     i_prev_glob = (i_glob   -1+N1)%N1;}
        else{            i_next_glob = (i_glob-N1+1)%N2 + N1;i_prev_glob = (i_glob-N1-1+N2)%N2 + N1;}
        _Status.insert(std::make_pair(i_next_glob,std::make_pair(i_glob, false)));
        _Status.insert(std::make_pair(i_prev_glob,std::make_pair(i_glob, false))); 
        mi1++;
        //std::cout<< "nb_prev= "<< 0 << " i_glob= " << i_glob << std::endl;
                                
        /******* Loop until the second polygon is reached *******/  
        while( !four_neighbours)
          {
            i_glob=(* mi1).second;//global index of vertex i
            nb_prev = _Status.count(i_glob);//counts the number of segments ending at i
                                                
            //std::cout<< "nb_prev= "<< nb_prev << " i_glob= " << i_glob << std::endl;
            switch (nb_prev)
              {
              case 1 :           
                mi=_Status.find(i_glob);// pointer to the segment ending at i
                i_prev_glob = ((*mi).second).first;//starting point of the segment ending at i
                i_next= (i_prev_glob - i_glob > 0) == (abs(i_prev_glob - i_glob) == 1)  ? i_glob - 1 : i_glob + 1;
                if(i_glob < N1) i_next_glob = (i_next   +N1)%N1;
                else            i_next_glob = (i_next-N1+N2)%N2 + N1;
                _Status.erase(mi);
                _Status.insert(std::make_pair(i_next_glob,std::make_pair(i_glob, false))); 
                mi1++;
                break;
              case 2 :
                return _Inter;
              case 0 :
                if( (i_glob < N1) != which_start)
                  {
                    mi2=mi1;
                    mi2++;
                    /* detection of double points */
                    if(distance2<DIM>((* mi1).first, (*mi2).first) > _epsilon)
                      four_neighbours = true;
                    else         /* Rare pothological case:  */
                      {
                        const std::pair< const double *, int > next_pt= *mi2;
                        events.erase(mi2);
                        mi1=events.insert(mi1,next_pt);
                      }
                  }
                break;
              default:
                throw Exception("intersectConvexPolygon: sequence of nodes does not describe a simple polygon (1)"); 
              }
          }
        /******** Loop until a terminal point or crossing is reached ************/
        while( !_terminus)  
          {
            //std::cout<< "nb_prev= "<< nb_prev<< " nb_inter= " << _Inter.size()/DIM << std::endl;
            switch (nb_prev)
              {
              case 1 :           
                mi=_Status.find(i_glob);// pointer to the segment ending at i
                i_prev_glob = ((*mi).second).first;//starting point of the segment ending at i
                sign = (i_prev_glob - i_glob > 0) == (abs(i_prev_glob - i_glob) == 1)  ? - 1 : + 1;
                i_next_glob = i_glob+sign;
                _is_in_intersection = ((*mi).second).second;//boolean that tells if i is in the intersection
                _Status.erase(mi);
                defineIndices(i_loc,i_next,i_prev, Poly1,Poly2,
                              j1,j1_glob,j2,j2_glob,j3,j3_glob,j4,j4_glob,
                              i_glob,i_next_glob,i_prev_glob, P_1,P_2, N1, N2, sign);
                if( _is_in_intersection ) addNewVertex(i_loc, i_glob, i_next_glob, i_prev_glob, Poly1);
                addCrossings(&Poly1[DIM*i_loc], &Poly1[DIM*i_next], i_glob, i_next_glob,
                             &Poly2[DIM*j1]   , &Poly2[DIM*j2]    , j1_glob,j2_glob,
                             &Poly2[DIM*j3]   , &Poly2[DIM*j4]    , j3_glob,j4_glob, &Poly1[DIM*i_prev]); 
                break;
              case 2 :
                if(!_Inter.empty())
                  {
                    if(i_glob < N1)  for(idim=0;idim<DIM;idim++) _Inter.push_back(P_1[DIM*i_glob+idim]);
                    else for(idim=0;idim<DIM;idim++) _Inter.push_back(P_2[DIM*(i_glob-N1)+idim]);
                  }
                return _Inter;
              case 0 ://To do if possible : remove this case from here
                if(_Inter.empty() && (i_glob < N1) != which_start){
                  i_next_glob = i_glob+1;
                  i_prev_glob = i_glob-1;
                  defineIndices(i_loc,i_next,i_prev, Poly1,Poly2,
                                j1,j1_glob,j2,j2_glob,j3,j3_glob,j4,j4_glob,
                                i_glob,i_next_glob,i_prev_glob, P_1,P_2, N1, N2, 1);
                  double V12[DIM], V34[DIM];
                  double inside = check_inside<DIM>(&Poly1[DIM*i_loc],&Poly2[DIM*j1],&Poly2[DIM*j2],
                                                    &Poly2[DIM*j3],   &Poly2[DIM*j4],V12, V34);      
                  _is_in_intersection=( inside < _epsilon ); // <= epsilon or 0 ?                
              
                  if(fabs(inside) > _epsilon)//vertex clearly inside or outside
                    {                                                                                        
                      //std::cout<<"coucou1" << std::endl;
                      if( _is_in_intersection)
                        {
                          for(int iidim=0;iidim<DIM;iidim++)
                            _Inter.push_back(Poly1[DIM*i_loc+iidim]);
                          _End_segments.push_back(std::make_pair(i_glob,i_next_glob));
                          _End_segments.push_back(std::make_pair(i_glob,i_prev_glob));
                        }
                      addCrossings(&Poly1[DIM*i_loc], &Poly1[DIM*i_next], i_glob, i_next_glob,
                                   &Poly2[DIM*j1]   , &Poly2[DIM*j2]    , j1_glob,j2_glob,
                                   &Poly2[DIM*j3]   , &Poly2[DIM*j4]    , j3_glob,j4_glob, &Poly1[DIM*i_prev]);
                      addCrossings(&Poly1[DIM*i_loc], &Poly1[DIM*i_prev], i_glob, i_prev_glob,
                                   &Poly2[DIM*j1]   , &Poly2[DIM*j2]    , j1_glob,j2_glob,
                                   &Poly2[DIM*j3]   , &Poly2[DIM*j4]    , j3_glob,j4_glob, &Poly1[DIM*i_next]); 
                    }
                  else //vertex on an edge
                    {
                      //std::cout<<"coucou2" << std::endl;
                      bool is_inside_next, is_inside_prev;
                      double Vnext[DIM], Vprev[DIM];
                      for(idim=0;idim<DIM;idim++) _Inter.push_back(Poly1[DIM*i_loc+idim]); 
                  
                      if(dotprod<DIM>(V34,V34) > _epsilon)//vertex i on edge (j1,j2), not on (j3,j4)
                        {
                          crossprod<DIM>(&Poly1[DIM*i_loc], &Poly2[DIM*j2], &Poly1[DIM*i_next],Vnext);
                          crossprod<DIM>(&Poly1[DIM*i_loc], &Poly2[DIM*j2], &Poly1[DIM*i_prev],Vprev);
                          is_inside_next= (dotprod<DIM>(Vnext,V34)<0);
                          is_inside_prev= (dotprod<DIM>(Vprev,V34)<0);
                     
                          if(!(is_inside_next || is_inside_prev)) return std::deque< double >();
                     
                          if(is_inside_next)
                            {
                              _End_segments.push_back(std::make_pair(i_glob,i_next_glob));
                              addCrossing0(&Poly1[DIM*i_loc], &Poly1[DIM*i_next], i_glob, i_next_glob,
                                           &Poly2[DIM*j3]   , &Poly2[DIM*j4]    , j3_glob,j4_glob); 
                            }
                          else
                            {
                              _End_segments.push_back(std::make_pair(j1_glob,j2_glob));
                              _Status.insert(std::make_pair(i_next_glob,std::make_pair(i_glob, false)));
                              mi=_Status.find(j2_glob); 
                              ((* mi).second).second= !((* mi).second).second;
                            }
                          if(is_inside_prev)
                            {
                              _End_segments.push_back(std::make_pair(i_glob,i_prev_glob));
                              addCrossing0(&Poly1[DIM*i_loc], &Poly1[DIM*i_prev], i_glob, i_prev_glob,
                                           &Poly2[DIM*j3]   , &Poly2[DIM*j4]    , j3_glob,j4_glob); 
                            }
                          else
                            {
                              _End_segments.push_back(std::make_pair(j1_glob,j2_glob));
                              _Status.insert(std::make_pair(i_prev_glob,std::make_pair(i_glob, false)));
                              mi=_Status.find(j2_glob);
                              ((* mi).second).second= !((* mi).second).second;
                            }
                        }
                      else if(dotprod<DIM>(V12,V12) > _epsilon)//vertex i on a edge (j3,j4), not on (j1,j2)
                        {
                          crossprod<DIM>(&Poly1[DIM*i_loc], &Poly2[DIM*j4], &Poly1[DIM*i_next],Vnext);
                          crossprod<DIM>(&Poly1[DIM*i_loc], &Poly2[DIM*j4], &Poly1[DIM*i_prev],Vprev);
                          is_inside_next= dotprod<DIM>(Vnext,V12)<0;
                          is_inside_prev= dotprod<DIM>(Vprev,V12)<0;
                     
                          if(!(is_inside_next || is_inside_prev)) return std::deque< double >();
                     
                          if(is_inside_next)
                            {
                              _End_segments.push_back(std::make_pair(i_glob,i_next_glob));
                              addCrossing0(&Poly1[DIM*i_loc], &Poly1[DIM*i_next], i_glob, i_next_glob,
                                           &Poly2[DIM*j1]   , &Poly2[DIM*j2]    , j1_glob,j2_glob); 
                            }
                          else
                            {
                              _End_segments.push_back(std::make_pair(j3_glob,j4_glob));
                              _Status.insert(std::make_pair(i_next_glob,std::make_pair(i_glob, false)));
                              mi=_Status.find(j4_glob); 
                              ((* mi).second).second= ! ((* mi).second).second;
                            }
                          if(is_inside_prev)
                            {
                              _End_segments.push_back(std::make_pair(i_glob,i_prev_glob));
                              addCrossing0(&Poly1[DIM*i_loc], &Poly1[DIM*i_prev], i_glob, i_prev_glob,
                                           &Poly2[DIM*j1]   , &Poly2[DIM*j2]    , j1_glob,j2_glob); 
                            }
                          else
                            {
                              _End_segments.push_back(std::make_pair(j3_glob,j4_glob));
                              _Status.insert(std::make_pair(i_prev_glob,std::make_pair(i_glob, false)));
                              mi=_Status.find(j4_glob); 
                              ((* mi).second).second= !((* mi).second).second;
                            }
                        }
                      else //vertices i, j1 and j3 share the same coordinates
                        {
                          crossprod<DIM>(&Poly1[DIM*i_loc], &Poly2[DIM*j2], &Poly1[DIM*i_next],Vnext);
                          crossprod<DIM>(&Poly1[DIM*i_loc], &Poly2[DIM*j2], &Poly1[DIM*i_prev],Vprev);
                          crossprod<DIM>(&Poly1[DIM*i_loc], &Poly2[DIM*j4], &Poly1[DIM*i_next],V12);
                          crossprod<DIM>(&Poly1[DIM*i_loc], &Poly2[DIM*j4], &Poly1[DIM*i_prev],V34);
                                                                                                        
                          double inside_next= dotprod<DIM>(Vnext,V12);
                          double inside_prev= dotprod<DIM>(Vprev,V34);
                          double inside_j2  = dotprod<DIM>(Vnext,Vprev);
                          double inside_j4  = dotprod<DIM>(V12,V34);
                                                                                                        
                          std::map<double, std::pair<int,int> > which_is_inside;
                          which_is_inside[inside_next] = std::make_pair(i_glob,i_next_glob);
                          which_is_inside[inside_prev] = std::make_pair(i_glob,i_prev_glob);
                          which_is_inside[inside_j2] =  std::make_pair(j1_glob,j2_glob);
                          which_is_inside[inside_j4] =  std::make_pair(j3_glob,j4_glob);

                          std::map<double, std::pair<int,int> >::iterator min = which_is_inside.begin();
                          std::map<double, std::pair<int,int> >::iterator minext = min;
                          minext++;
                          std::map<double, std::pair<int,int> >::reverse_iterator max = which_is_inside.rbegin();
                          std::multimap< int, std::pair< int,bool> >::iterator j2_in_status = _Status.find(((*min).second).second);
                          std::multimap< int, std::pair< int,bool> >::iterator j4_in_status = _Status.find(((*minext).second).second);

                          if((*min).first < -_epsilon) //there is someone clearly inside
                            {
                              _End_segments.push_back( (*min).second );
                              _End_segments.push_back((* minext).second);
                              if(j2_in_status != _Status.end())
                                ((*j2_in_status).second).second        = !        ((*j2_in_status).second).second;                                                                                        
                              if(j4_in_status != _Status.end())
                                ((*j4_in_status).second).second        = !        ((*j4_in_status).second).second;                                                                                        
                              is_inside_next = ((*min).second).second == i_next_glob || ((*minext).second).second == i_next_glob;
                              is_inside_prev = ((*min).second).second == i_prev_glob || ((*minext).second).second == i_prev_glob;
                            }
                          else
                            if(fabs((*min).first) <= _epsilon) //nobody is clearly inside but two segments are superposed
                              {
                                if(fabs((*max).first) > _epsilon) 
                                  return std::deque< double >();
                                else //all four segments are superposed
                                  {
                                    _End_segments.push_back(std::make_pair(i_glob,i_next_glob));
                                    _End_segments.push_back(std::make_pair(i_glob,i_prev_glob));
                                    is_inside_next= true;
                                    is_inside_prev= true;
                                  }
                              } 
                            else                //there is nobody inside
                              return std::deque< double >();
                                                                                                                         
                          _Status.insert(std::make_pair(i_prev_glob,std::make_pair(i_glob,is_inside_prev)));
                          _Status.insert(std::make_pair(i_next_glob,std::make_pair(i_glob,is_inside_next)));
                        }
                    }
                }
                break;
              default:
                std::cout << "Problem: nbprev= " << nb_prev << " ; i_glob = " << i_glob << std::endl;
                throw Exception("intersectConvexPolygon: sequence of nodes does not describe a simple polygon (2)"); 
              } 
            mi1++;
            i_glob=(* mi1).second;//global index of vertex i
            nb_prev = _Status.count(i_glob);
          }
      }
    return _Inter;
  }
  
  /**************************************************************************/
  /* computes the convex hull of a polygon subP which is a sub polygon of P */
  /* P is the array of coordinates, subP is a map containing initially the indices of a subpolygon of P */
  /* in the end, subP contains only the elements belonging to the convex hull, and  not_in_hull the others */
  /**************************************************************************/
  template<int DIM>
  inline void PolygonAlgorithms<DIM>::convHull(const double *P, int N, double * normal,  
                                               std::map< int,int >& subP, std::map< int,int >& not_in_hull,
                                               int& NsubP, const double epsilon)
  {
    if(NsubP>3)
      {
        std::map< int,int >::iterator mi_prev = subP.begin();
        std::map< int,int >::iterator  mi = mi_prev;
        mi++;
        std::map< int,int >::iterator  mi_next = mi;
        mi_next++;
        double directframe=0.;
       
        /* Check if the polygon subP is positively oriented */
        std::map< int,int >::iterator mi1=mi;
        while(mi1 != subP.end() && distance2<DIM>(&P[DIM*(*subP.begin()).second],&P[DIM*(*mi1).second])< epsilon) 
          mi1++;
        std::map< int,int >::iterator mi2=mi1;
        while(mi2 != subP.end() && fabs(directframe)<epsilon)
          {
            directframe =direct_frame<DIM>(&P[DIM* (*mi1).second],
                                           &P[DIM* (*subP.begin()).second],
                                           &P[DIM* (*mi2).second], normal);
            mi2++;
          }
        if(directframe < 0) for(int idim=0; idim< DIM; idim++) normal[idim] *= -1;
       
        /* Core of the algorithm */
        while(mi_next != subP.end())
          {
            directframe = direct_frame<DIM>(&P[DIM* (*mi).second],
                                            &P[DIM* (*mi_prev).second],
                                            &P[DIM* (*mi_next).second], normal);
            if(directframe > -epsilon){
              mi ++;
              mi_prev++;
              mi_next++;
            }
            else
              {
                not_in_hull.insert(*mi);
                subP.erase(mi); 
                NsubP--;
                mi--;
              }
          }
        directframe = direct_frame<DIM>(&P[DIM*(*mi).second],
                                        &P[DIM*(*mi_prev).second],
                                        &P[DIM*(*subP.begin()).second], normal);
        if(directframe < -epsilon)
          {
            not_in_hull.insert(*mi);
            subP.erase(mi); 
            NsubP--;
          }
      }
  }
  
  template<int DIM>
  void PolygonAlgorithms<DIM>::convexDecomposition(const double * P, int N, double *normal, std::vector< int > subP, int NsubP,
                                                   std::vector< std::map< int,int > >& components, std::vector< int >& components_index,
                                                   int& Ncomp, int sign, const double epsilon)
  {
    int i;
    std::map< int, int > hull;
    std::map< int, int > not_in_hull;
    std::map< int, int >::iterator mi, mj;
    std::vector< int > reflex_region;
    int Nreflex;
    int i_xmax=0;
    const double * xmax=&P[DIM*subP[0]];
    /* checking an extremal point of subP */
    for(i=0; i<NsubP; i++)
      {
        if(&P[DIM*subP[i]]> xmax)
          {
            i_xmax=i;
            xmax=&P[DIM*subP[i]];
          }
      }
    /* renumbering of SubP elements for the convex hull*/
    for(i=0; i<NsubP; i++) hull.insert(hull.end(),std::make_pair(i,subP[(i+i_xmax)%NsubP]));
    /* compute the convex hull */
    convHull(P, N, normal, hull, not_in_hull, NsubP,epsilon);
    /* convex hull is the next component */
    components.push_back(hull);
    components_index.push_back(sign*NsubP);
    Ncomp++;
    /* searching for reflex regions */
    for(mi=not_in_hull.begin(); mi!=not_in_hull.end(); mi++)
      {
        reflex_region.clear();
        reflex_region.push_back(hull[(*mi).first-1]);
        reflex_region.push_back( (*mi).second );
        Nreflex=2;
        mj=mi;
        mj++;
        while((mj != not_in_hull.end()) && ((*mj).first == (*mi).first+1))
          {
            reflex_region.push_back((*mj).second);
            Nreflex++;       
            mi++;
            mj++;
          }
        reflex_region.push_back(hull[(*mi).first+1]);
        Nreflex++;       
        convexDecomposition( P, N,normal, reflex_region, Nreflex, components, components_index, Ncomp, -sign, epsilon);
      }
  }

  /**************************************************************************/
  /* decomposes a non convex polygon P with N vertices contained in a plane */
  /* into a sequence of convex polygons */
  /* the input vectors 'components' and 'components_index' should be empty */
  /* returns the number of convex components */
  /* if P is composed of a single point, then an empty polygon is returned */
  /**************************************************************************/
  template<int DIM>
  int PolygonAlgorithms<DIM>::convexDecomposition(const double * P, int N, std::vector< std::map< int,int > >& components,
                                                  std::vector< int >& components_index, const double epsilon)
  {
    int Ncomp=0;
    std::vector< int > subP(N);
    double normal[3]={0,0,0};

    for(int i = 0; i<N; i++) subP[i]=i;

    //Build the normal of polygon P
    int i1=1;
    while(i1<N && distance2<DIM>(&P[0],&P[i1])< epsilon) i1++;
    int i2=i1+1;
    while(i2<N && fabs(dotprod<DIM>(normal,normal))<epsilon)
      {
        crossprod<DIM>(&P[i1], &P[0], &P[i2],normal);
        i2++;
      }
    
    convexDecomposition(P, N, normal, subP, N, components, components_index, Ncomp, 1, epsilon);
    return Ncomp;
  }
}

#endif
