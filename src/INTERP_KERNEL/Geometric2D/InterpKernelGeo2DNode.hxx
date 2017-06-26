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

#ifndef __INTERPKERNELGEO2DNODE_HXX__
#define __INTERPKERNELGEO2DNODE_HXX__

#include "InterpKernelGeo2DPrecision.hxx"
#include "INTERPKERNELDefines.hxx"

#include <map>
#include <cmath>
#include <vector>
#include <iostream>

namespace INTERP_KERNEL
{
  typedef enum
  {
    IN_1      =  7,
    ON_1      =  8,
    ON_LIM_1  = 12,
    ON_TANG_1 =  9,
    OUT_1     = 10,
    UNKNOWN   = 11
  } TypeOfLocInPolygon;

  class Bounds;

  /*!
   * Representation of a 2D point, and potentially its location relative to a polygon.
   * As nodes can be shared between edges it is handled with ref counting.
   */
  class INTERPKERNEL_EXPORT Node
  {
  public:
    Node(double x, double y);
    Node(const double *coords);
    Node(std::istream& stream);
    void incrRef() const { _cnt++; }
    bool decrRef();
    void initHitStatus() const { _hit=0; }
    char getHitStatus() const { return _hit; }
    void hitMeAlone(double xBary, double yBary, double dimChar) { if(_hit==0) { _hit=1; applySimilarity(xBary,yBary,dimChar); } }
    void unHitMeAlone(double xBary, double yBary, double dimChar) { if(_hit==0) { _hit=1; unApplySimilarity(xBary,yBary,dimChar); } }
    void hitMeAfter(double xBary, double yBary, double dimChar) { if(_hit==0) { hitMeAlone(xBary,yBary,dimChar); _hit=2; } else if(_hit==1) declareOn(); }
    void unHitMeAfter(double xBary, double yBary, double dimChar) { if(_hit==0) { unHitMeAlone(xBary,yBary,dimChar); _hit=2; } }
    void initLocs() const { _loc=UNKNOWN; }
    void setLoc(TypeOfLocInPolygon loc) const { _loc=loc; }
    TypeOfLocInPolygon getLoc() const { return _loc; }
    void declareIn() const { if(_loc==UNKNOWN) _loc=IN_1; }
    void declareOn() const { if(_loc==UNKNOWN) _loc=ON_1; }
    void declareOnLim() const { if(_loc==UNKNOWN || _loc==ON_1) _loc=ON_LIM_1; }
    void declareOut() { if(_loc==UNKNOWN) _loc=OUT_1; }
    void declareOnTangent() { _loc=ON_TANG_1; }
    operator const double*() const { return _coords; }
    bool isEqual(const Node& other) const;
    //returns an angle in -Pi/2;Pi/2.
    double getSlope(const Node& other) const;
    bool isEqualAndKeepTrack(const Node& other, std::vector<Node *>& track) const;
    void dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const;
    double distanceWithSq(const Node& other) const;
    double operator[](int i) const { return _coords[i]; }
    //! use with caution
    void setNewCoords(double x, double y) { _coords[0]=x; _coords[1]=y; }
    //returns an angle in -Pi/2;Pi/2.
    static double computeSlope(const double *pt1, const double *pt2);
    //returns an angle in -Pi;Pi
    static double computeAngle(const double *pt1, const double *pt2);
    void applySimilarity(double xBary, double yBary, double dimChar);
    void unApplySimilarity(double xBary, double yBary, double dimChar);
    static double dot(const double *vect1, const double *vect2) { return vect1[0]*vect2[0]+vect1[1]*vect2[1]; }
    static double sign(double val) { if(val>=0) return 1.; else return -1.; }
    static double norm(const double *vect) { return sqrt(vect[0]*vect[0]+vect[1]*vect[1]); }
    static bool areDoubleEquals(double a, double b) { return fabs(a-b) < QuadraticPlanarPrecision::getPrecision(); }
    //! idem areDoubleEquals except that precision of comparison is modified.
    static bool areDoubleEqualsWP(double a, double b, double k) { return fabs(a-b) < k*QuadraticPlanarPrecision::getPrecision(); }
    static double distanceBtw2Pt(const double *a, const double *b) { return sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])); }
    static double distanceBtw2PtSq(const double *a, const double *b) { return (a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1]); }
    //
    void fillGlobalInfoAbs(const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                           std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int>& mapAddCoo, int *nodeId) const;
    void fillGlobalInfoAbs2(const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                            std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int>& mapAddCoo, std::vector<int>& pointsOther) const;
  protected:
    ~Node();
  protected:
    mutable char _hit;
    mutable unsigned char _cnt;
    mutable TypeOfLocInPolygon _loc;
    double _coords[2];
  };
}

#endif
