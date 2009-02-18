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
// File      : UnitTetraIntersectionBary.cxx
// Created   : Tue Dec  9 16:48:49 2008
// Author    : Edward AGAPOV (eap)

#include "UnitTetraIntersectionBary.hxx"

#include "VectorUtils.hxx"
#include "InterpolationUtils.hxx"
#include "VolSurfFormulae.hxx"

#define NB_TETRA_SIDES 4
#define NB_TETRA_NODES 4

using namespace std;

namespace INTERP_KERNEL
{
  enum { _X=0, _Y, _Z };

  inline bool samePoint( const double* p1, const double* p2 )
  {
    return ( p1[0]==p2[0] && p1[1]==p2[1] && p1[2]==p2[2]);
  }

  //================================================================================
  /*!
   * \brief Creates a ready-to-use tool
   */
  //================================================================================

  UnitTetraIntersectionBary::UnitTetraIntersectionBary():TransformedTriangle()
  {
    _int_volume = 0;
    //init();
  }
  //================================================================================
  /*!
   * \brief Initializes fields
   */
  //================================================================================

  void UnitTetraIntersectionBary::init()
  {
    _int_volume = 0;
    _faces.clear();
  }
  
  //================================================================================
  /*!
   * \brief Stores a part of triangle common with the unit tetrahedron
   *  \param triangle - triangle side of other cell
   */
  //================================================================================

  void UnitTetraIntersectionBary::addSide(const TransformedTriangle& triangle)
  {
    _int_volume += triangle.getVolume();

    double triNormal[3], polyNormal[3];
    crossprod<3>( triangle.getCorner(P),triangle.getCorner(Q),triangle.getCorner(R), triNormal);

    const std::vector<double*> * pPolygonA = &triangle.getPolygonA();
    if ( pPolygonA->size() < 3 )
      {
        if ( !epsilonEqual( triNormal[_Z], 0 ))
          return; // not vertical triangle does not intersect the unit tetra

        // Vertical triangle. Use inherited methods of TransformedTriangle to
        // calculate intesection polygon
        *((TransformedTriangle*)this) = triangle; // copy triangle fields
        _polygonA.clear();
        _polygonB.clear();
        calculateIntersectionPolygons();
        if (this->_polygonA.size() < 3)
          return;
        calculatePolygonBarycenter(A, _barycenterA);
        sortIntersectionPolygon(A, _barycenterA);
        pPolygonA = & _polygonA;
      }

    // check if polygon orientation is same as the one of triangle
    vector<double*>::const_iterator p = pPolygonA->begin(), pEnd = pPolygonA->end();
    double* p1 = *p++;
    double* p2 = *p;
    while ( samePoint( p1, p2 ) && ++p != pEnd )
      p2 = *p;
    if ( p == pEnd )
      {
        clearPolygons();
        return;
      }
    double* p3 = *p;
    while ( samePoint( p2, p3 ) && ++p != pEnd )
      p3 = *p;
    if ( p == pEnd )
      {
        clearPolygons();
        return ;
      }
    crossprod<3>( p1, p2, p3, polyNormal );
    bool reverse = ( dotprod<3>( triNormal, polyNormal ) < 0.0 );

    // store polygon
    _faces.push_back( vector< double* > () );
    vector< double* >& faceCorner = _faces.back();
    faceCorner.resize( pPolygonA->size()/* + 1*/ );

    int i = 0;
    if ( reverse )
      {
        vector<double*>::const_reverse_iterator polyF = pPolygonA->rbegin(), polyEnd;
        for ( polyEnd = pPolygonA->rend(); polyF != polyEnd; ++i, ++polyF )
          if ( i==0 || !samePoint( *polyF, faceCorner[i-1] ))
            copyVector3( *polyF, faceCorner[i] = new double[3] );
          else
            --i;
      }
    else
      {
        vector<double*>::const_iterator polyF = pPolygonA->begin(), polyEnd;
        for ( polyEnd = pPolygonA->end(); polyF != polyEnd; ++i, ++polyF )
          if ( i==0 || !samePoint( *polyF, faceCorner[i-1] ))
            copyVector3( *polyF, faceCorner[i] = new double[3] );
          else
            --i;
      }
    if ( i < 3 )
      {
        clearPolygons(); // free memory of _polygonA
        _polygonA = faceCorner;
        _faces.pop_back();
      }
    else
      {
        if ( i < (int)pPolygonA->size() )
          faceCorner.resize( i );
      }

#ifdef DMP_ADDSIDE
    cout << "**** addSide() " << _faces.size() << endl;
    for ( int i = 0; i < faceCorner.size(); ++i )
      {
        double* p = faceCorner[i];
        cout << i << ": ( " << p[0] << ", " << p[1] << ", " << p[2] << " )" << endl;
      }
#endif
    clearPolygons();
  }

  //================================================================================
  /*!
   * \brief Computes and returns coordinates of barycentre
   */
  //================================================================================

  bool UnitTetraIntersectionBary::getBary(double* baryCenter)
  {
    baryCenter[0] = baryCenter[1] = baryCenter[2] = -1.0;
    if ( addSideFaces() < NB_TETRA_SIDES )
      {
        // tetra is not intersected
        if ( fabs(_int_volume) > 1e-10 )
          {
            // tetra is fully inside the other cell
            baryCenter[0] = baryCenter[1] = baryCenter[2] = 0.25;
            _int_volume = 0.16666666666666666;
            return true;
          }
        return false;
      }
    // Algo:
    // - pick up one point P among the summits of the polyhedron
    // - for each face of the polyhedron which does not contain the point :
    //   - compute the barycenter of the volume obtained by forming the "pyramid" with
    //     the face as a base and point P as a summit
    //   - compute the volume of the "pyramid"
    // - Add up all barycenter positions weighting them with the volumes.

    baryCenter[0] = baryCenter[1] = baryCenter[2] = 0.;

    list< vector< double* > >::iterator f = _faces.begin(), fEnd = _faces.end();
    double * P = f->at(0);

    for ( ++f; f != fEnd; ++f )
      {
        vector< double* >& polygon = *f;
        if ( polygon.empty() )
          continue;

        bool pBelongsToPoly = false;
        vector<double*>::iterator v = polygon.begin(), vEnd = polygon.end();
        for ( ; !pBelongsToPoly && v != vEnd; ++v )
          pBelongsToPoly = samePoint( P, *v );
        if ( pBelongsToPoly )
          continue;

        // Compute the barycenter of the volume. Barycenter of pyramid is on line
        // ( barycenter of polygon -> P ) with 1/4 of pyramid height from polygon.

        double bary[] = { 0, 0, 0 };

        // base polygon bary
        for ( v = polygon.begin(); v != vEnd ; ++v )
          {
            double* p = *v;
            bary[0] += p[0];
            bary[1] += p[1];
            bary[2] += p[2];
          }
        bary[0] /= polygon.size();
        bary[1] /= polygon.size();
        bary[2] /= polygon.size();

        // pyramid volume
        double vol = 0;
        for ( int i = 0; i < (int)polygon.size(); ++i )
          {
            double* p1 = polygon[i];
            double* p2 = polygon[(i+1)%polygon.size()];
            vol += std::fabs( calculateVolumeForTetra( p1, p2, bary, P ));
          }

        // put bary on the line ( barycenter of polygon -> P ) and multiply by volume
        baryCenter[0] += ( bary[0] * 0.75 + P[0] * 0.25 ) * vol;
        baryCenter[1] += ( bary[1] * 0.75 + P[1] * 0.25 ) * vol;
        baryCenter[2] += ( bary[2] * 0.75 + P[2] * 0.25 ) * vol;
      }
    if ( _int_volume < 0. )
      _int_volume = -_int_volume;
    baryCenter[0] /= _int_volume;
    baryCenter[1] /= _int_volume;
    baryCenter[2] /= _int_volume;

    return true;
  }

  //================================================================================
  /*!
   * \brief Add faces of the intersection polyhedron formed on faces of the
   *        unit tetrahedron by sides of already added faces
   *  \retval int - number of faces of intersection polyhedron
   */
  //================================================================================

  int UnitTetraIntersectionBary::addSideFaces()
  {
    int nbPolyhedraFaces = 0;

    if ( _faces.empty() )
      return nbPolyhedraFaces;

    // -------------------------------------------
    // Detect polygons laying on sides of a tetra
    // -------------------------------------------

    bool sideAdded[NB_TETRA_SIDES] = { false, false, false, false };
    int nbAddedSides = 0;
    list< vector< double* > >::iterator f = _faces.begin(), fEnd = _faces.end();
    for ( ; f != fEnd; ++f )
      {
        vector< double* >& polygon = *f;
        double coordSum[3] = {0,0,0};
        for ( int i = 0; i < (int)polygon.size(); ++i )
          {
            double* p = polygon[i];
            coordSum[0] += p[0];
            coordSum[1] += p[1];
            coordSum[2] += p[2];
          }
        for ( int j = 0; j < 3 && !sideAdded[j]; ++j )
          {
            if ( coordSum[j] == 0 )
              sideAdded[j] = bool( ++nbAddedSides );
          }
        if ( !sideAdded[3] &&
             ( epsilonEqual( (coordSum[0]+coordSum[1]+coordSum[2]) / polygon.size(), 1. )))
          sideAdded[3] = bool( ++nbAddedSides );
      }
    if ( nbAddedSides == NB_TETRA_SIDES )
      return nbAddedSides;

    // ---------------------------------------------------------------------------------
    // Add segments of already added polygons to future polygonal faces on sides of tetra
    // ---------------------------------------------------------------------------------

    int nbIntersectPolygs = _faces.size();

    vector< double* > * sideFaces[ 4 ]; // future polygons on sides of tetra
    for ( int i = 0; i < NB_TETRA_SIDES; ++i )
      {
        sideFaces[ i ]=0;
        if ( !sideAdded[ i ] )
          {
            _faces.push_back( vector< double* > () );
            sideFaces[ i ] = &_faces.back();
          }
      }
    f = _faces.begin(), fEnd = _faces.end();
    for ( int iF = 0; iF < nbIntersectPolygs; ++f, ++iF ) // loop on added intersection polygons
      {
        vector< double* >& polygon = *f;
        for ( int i = 0; i < (int)polygon.size(); ++i )
          {
            // segment ends
            double* p1 = polygon[i];
            double* p2 = polygon[(i+1)%polygon.size()];
            bool p1OnSide, p2OnSide;//, onZeroSide = false;
            for ( int j = 0; j < 3; ++j )
              {
                if ( !sideFaces[ j ] )
                  continue;
                p1OnSide = ( p1[j] == 0 );
                p2OnSide = ( p2[j] == 0 );
                if ( p1OnSide && p2OnSide )
                  {
                    // segment p1-p2 is on j-th orthogonal side of tetra
                    sideFaces[j]->push_back( new double[3] );
                    copyVector3( p1, sideFaces[j]->back() );
                    sideFaces[j]->push_back( new double[3] );
                    copyVector3( p2, sideFaces[j]->back() );
                    //break;
                  }
              }
            // check if the segment p1-p2 is on the inclined side
            if ( sideFaces[3] &&
                 epsilonEqual( p1[_X] + p1[_Y] + p1[_Z], 1.0 ) &&
                 epsilonEqual( p2[_X] + p2[_Y] + p2[_Z], 1.0 ))
              {
                sideFaces[3]->push_back( new double[3] );
                copyVector3( p1, sideFaces[3]->back() );
                sideFaces[3]->push_back( new double[3] );
                copyVector3( p2, sideFaces[3]->back() );
              }
          }
      }
#ifdef DMP_ADDSIDEFACES
    cout << "**** after Add segments to sides " << endl;
    for ( int i = 0; i < NB_TETRA_SIDES; ++i )
      {
        cout << "\t Side " << i << endl;
        if ( !sideFaces[i] )
          {
            cout << "\t cut by triagle" << endl;
          }
        else
          {
            vector< double* >& sideFace = *sideFaces[i];
            for ( int i = 0; i < sideFace.size(); ++i )
              {
                double* p = sideFace[i];
                cout << "\t" << i << ": ( " << p[0] << ", " << p[1] << ", " << p[2] << " )" << endl;
              }
          }
      }
#endif

    // ---------------------------------------------------------------------------
    // Make closed polygons on tetra sides by adding not cut off corners of tetra
    // ---------------------------------------------------------------------------

    double origin[3] = { 0,0,0 };

    // find corners of tetra cut off by triangles of other tetra
    // ---------------------------------------------------------

    // corners are coded like this: index = 1*X + 2*Y + 3*Z
    // (0,0,0) -> index == 0; (0,0,1) -> index == 3
    vector< int > cutOffCorners(NB_TETRA_NODES, false), passedCorners(NB_TETRA_NODES, false);

    int nbCutOffCorners = 0;
    for ( int i = 0; i < 3; ++i ) // loop on orthogonal faces of the unit tetra
      {
        if ( !sideFaces[i] ) continue;
        vector< double* >& sideFace = *sideFaces[i];

        int nbPoints = sideFace.size();
        if ( nbPoints == 0 )
          continue; // not intersected face at all - no cut off corners can be detected

        int ind1 = (i+1)%3, ind2 = (i+2)%3; // indices of coords on i-th tetra side

        int nbCutOnSide = 0;
        bool isSegmentOnEdge;
        for ( int ip = 0; ip < nbPoints; ++ip )
          {
            int isSegmentEnd = ( ip % 2 );

            double* p = sideFace[ ip ];
            double* p2 = isSegmentEnd ? 0 : sideFace[ip+1];

            if ( !isSegmentEnd )
              isSegmentOnEdge = false; // initialize

            int cutOffIndex = -1, passThIndex = -1;// no cut off neither pass through
            int pCut[] = { 0,0,0 }, pPass[] = { 0,0,0 };

            if ( p[ind1]==0.)
              {
                // point is in on orthogonal edge
                if ( !isSegmentEnd && p2[ind1]==0 )
                  isSegmentOnEdge = true;

                if ( !isSegmentOnEdge )
                  { // segment ends are on different edges
                    pCut[ind2] = isSegmentEnd; // believe that cutting triangles are well oriented
                    cutOffIndex = pCut[0] + 2*pCut[1] + 3*pCut[2];
                  }
                if ( p[ind2]==0. || p[ind2]==1.)
                  {
                    pPass[ind2] = int(p[ind2]);
                    passThIndex = pPass[0] + 2*pPass[1] + 3*pPass[2];
                  }
              }
            else if ( p[ind2]==0.)
              {
                // point is on orthogonal edge
                if ( !isSegmentEnd && p2[ind2]==0 )
                  isSegmentOnEdge = true;
                if ( !isSegmentEnd )
                  {// segment ends are on different edges
                    pCut[ind1] = 1-isSegmentEnd;
                    cutOffIndex = pCut[0] + 2*pCut[1] + 3*pCut[2];
                  }
                if ( p[ind1]==0. || p[ind1]==1.)
                  {
                    pPass[ind1] = int(p[ind1]);
                    passThIndex = pPass[0] + 2*pPass[1] + 3*pPass[2];
                  }
              }
            else if ( epsilonEqual(p[ind1] + p[ind2], 1.0 ))
              {
                // point is on inclined edge
                if ( !isSegmentEnd && epsilonEqual(p2[ind1] + p2[ind2], 1.0 ))
                  isSegmentOnEdge = true;
                if ( !isSegmentOnEdge )
                  { //segment ends are on different edges
                    pCut[ind1] = isSegmentEnd;
                    pCut[ind2] = 1-isSegmentEnd;
                    cutOffIndex = pCut[0] + 2*pCut[1] + 3*pCut[2];
                  }
              }
            else
              {
                continue;
              }
            // remember cut off and passed through points
            if ( passThIndex >= 0. )
              {
                passedCorners[ passThIndex ] = true;
                if ( cutOffCorners[ passThIndex ] )
                  {
                    nbCutOffCorners--;
                    cutOffCorners[ passThIndex ] = false;
                  }
              }
            if ( cutOffIndex >= 0. )
              {
                nbCutOnSide++;
                if ( !passedCorners[ cutOffIndex ] && !cutOffCorners[ cutOffIndex ] )
                  {
                    nbCutOffCorners++;
                    cutOffCorners[ cutOffIndex ] = true;
                  }
              }
          } // loop on points on a unit tetra side

        if ( nbCutOnSide == 0 && nbPoints <= 2 )
          continue; // one segment laying on edge at most

        if ( nbCutOffCorners == NB_TETRA_NODES )
          break; // all tetra corners are cut off

        if ( /*nbCutOnSide <= 2 &&*/ nbPoints >= 6 )
          {
            // at least 3 segments - all corners of a side are cut off
            for (int cutIndex = 0; cutIndex < NB_TETRA_NODES; ++cutIndex )
              if ( cutIndex != i+1 && !passedCorners[ cutIndex ] && !cutOffCorners[ cutIndex ])
                cutOffCorners[ cutIndex ] = bool ( ++nbCutOffCorners );
          }

      } // loop on orthogonal faces of tetra

    // check if all corners are cut off on the inclined tetra side
    if ( sideFaces[ XYZ ] && sideFaces[ XYZ ]->size() >= 6 )
      {
        for (int cutIndex = 1; cutIndex < NB_TETRA_NODES; ++cutIndex )
          if ( !passedCorners[ cutIndex ] && !cutOffCorners[ cutIndex ])
            cutOffCorners[ cutIndex ] = bool ( ++nbCutOffCorners );
      }

    // Add to faces on tetra sides the corners not cut off by segments of intersection polygons
    // ----------------------------------------------------------------------------------
    if ( nbCutOffCorners > 0 )
      {
        for ( int i = 0; i < NB_TETRA_SIDES; ++i )
          {
            if ( !sideFaces[ i ] ) continue;
            vector< double* >& sideFace = *sideFaces[i];

            int excludeCorner = (i + 1) % NB_TETRA_NODES;
            for ( int ic = 0; ic < NB_TETRA_NODES; ++ic )
              {
                if ( !cutOffCorners[ ic ] && ic != excludeCorner )
                  {
                    sideFace.push_back( new double[3] );
                    copyVector3( origin, sideFace.back() );
                    if ( ic )
                      sideFace.back()[ ic-1 ] = 1.0;
                  }
              }
          }
      }

#ifdef DMP_ADDSIDEFACES
    cout << "**** after Add corners to sides " << endl;
    for ( int i = 0; i < NB_TETRA_SIDES; ++i )
      {
        cout << "\t Side " << i << endl;
        if ( !sideFaces[i] ) {
          cout << "\t cut by triagle" << endl;
        }
        else 
          {
            vector< double* >& sideFace = *sideFaces[i];
            for ( int i = 0; i < sideFace.size(); ++i )
              {
                double* p = sideFace[i];
                cout << "\t" << i << ": ( " << p[0] << ", " << p[1] << ", " << p[2] << " )" << endl;
              }
          }
      }
    cout << "Cut off corners: ";
    if ( nbCutOffCorners == 0 )
      cout << "NO";
    else 
      for ( int ic = 0; ic < NB_TETRA_NODES; ++ic )
        cout << cutOffCorners[ ic ];
    cout << endl;
#endif
    // ------------------------------------------------------------------------
    // Sort corners of filled up faces on tetra sides and exclude equal points
    // ------------------------------------------------------------------------

    int iF = 0;
    for ( f = _faces.begin(); f != fEnd; ++f, ++iF )
      {
        vector< double* >&  face = *f;
        if ( face.size() >= 3 )
          {
            clearPolygons(); // free memory of _polygonA
            _polygonA = face;
            face.clear();
            face.reserve( _polygonA.size() );
            if ( iF >= nbIntersectPolygs )
              { // sort points of side faces
                calculatePolygonBarycenter( A, _barycenterA );
                setTriangleOnSide( iF - nbIntersectPolygs );
                sortIntersectionPolygon( A, _barycenterA );
              }
            // exclude equal points
            vector< double* >::iterator v = _polygonA.begin(), vEnd = _polygonA.end();
            face.push_back( *v );
            *v = 0;
            for ( ++v; v != vEnd; ++v )
              {
                double* pPrev = face.back();
                double* p     = *v;
                if ( !samePoint( p, pPrev ))
                  {
                    face.push_back( p );
                    *v = 0;
                  }
              }
          }
        if ( face.size() < 3 )
          { // size could decrease
            clearPolygons(); // free memory of _polygonA
            _polygonA = face;
            face.clear();
          }
        else
          {
            nbPolyhedraFaces++;
          }
      }
#ifdef DMP_ADDSIDEFACES
    cout << "**** after HEALING all faces " << endl;
    for (iF=0, f = _faces.begin(); f != fEnd; ++f, ++iF )
      {
        cout << "\t Side " << iF << endl;
        vector< double* >& sideFace = *f;
        for ( int i = 0; i < sideFace.size(); ++i )
          {
            double* p = sideFace[i];
            cout << "\t" << i << ": ( " << p[0] << ", " << p[1] << ", " << p[2] << " )" << endl;
          }
      }
#endif
    return nbPolyhedraFaces;
  }

  //================================================================================
  /*!
   * \brief set corners of inherited TransformedTriangle as corners of i-th side of
   * the Unit tetra. It is necessary to sort points of faces on sides of the unit
   * tetrahedron using sortIntersectionPolygon(A)
   */
  //================================================================================

  void UnitTetraIntersectionBary::setTriangleOnSide(int iSide)
  {
    if ( iSide >= 3 )
      iSide = 0;
    for(int i = 0 ; i < 3 ; ++i) 
      {
        _coords[5*i] = _coords[5*i + 1] = _coords[5*i + 2] = 0.;
        if ( i != iSide )
          _coords[5*i + i] = 1.;
      }
  }

  //================================================================================
  /*!
   * \brief Free memory of polygons
   */
  //================================================================================

  void UnitTetraIntersectionBary::clearPolygons(bool andFaces)
  {
    for(vector<double*>::iterator it = _polygonA.begin() ; it != _polygonA.end() ; ++it)
      {  delete[] *it;
        *it = 0; 
      }
    for(vector<double*>::iterator it = _polygonB.begin() ; it != _polygonB.end() ; ++it)
      { 
        delete[] *it; 
        *it = 0; 
      }

    _polygonA.clear();
    _polygonB.clear();

    if ( andFaces )
      {
        list< vector< double* > >::iterator f = this->_faces.begin(), fEnd = this->_faces.end();
        for ( ; f != fEnd; ++f )
          {
            vector< double* >& polygon = *f;
            for(vector<double*>::iterator it = polygon.begin() ; it != polygon.end() ; ++it)
              { 
                delete[] *it;
                *it = 0;
              }
          }
        this->_faces.clear();
      }
  }

  //================================================================================
  /*!
   * \brief Destructor clears coordinates of faces
   */
  //================================================================================

  UnitTetraIntersectionBary::~UnitTetraIntersectionBary()
  {
    clearPolygons(/*andFaces=*/true );
  }

}
