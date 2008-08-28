#ifndef __NODE_HXX__
#define __NODE_HXX__

#include "Geometric2D_defines.hxx"

#include "Precision.hxx"

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
   * As nodes can be shared between edges it is dealed with ref counting.
   */
  class GEOMETRIC2D_EXPORT Node
  {
  public:
    Node(double x, double y);
    Node(const double *coords);
    Node(std::istream& stream);
    void incrRef() const { _cnt++; }
    bool decrRef();
    TypeOfLocInPolygon getLoc() const { return _loc; }
    void declareIn() const { if(_loc==UNKNOWN) _loc=IN_1; }
    void declareOn() const { if(_loc==UNKNOWN) _loc=ON_1; }
    void declareOnLim() const { if(_loc==UNKNOWN || _loc==ON_1) _loc=ON_LIM_1; }
    void declareOut() { if(_loc==UNKNOWN) _loc=OUT_1; }
    void declareOnTangent() { _loc=ON_TANG_1; }
    operator const double*() const { return _coords; }
    bool isEqual(const Node& other) const;
    double getSlope(const Node& other) const;
    bool isEqualAndKeepTrack(const Node& other, std::vector<Node *>& track) const;
    void dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const;
    double distanceWithSq(const Node& other) const;
    double operator[]( int i) const { return _coords[i]; }
    //!for tests only !
    void setNewCoords(double x, double y) { _coords[0]=x; _coords[1]=y; }
    static double sign(double val) { if(val>=0) return 1.; else return -1.; }
    static bool areDoubleEquals(double a, double b) { return fabs(a-b) < QUADRATIC_PLANAR::_precision; }
    //! idem areDoubleEquals except that precision of comparison is modified.
    static bool areDoubleEqualsWP(double a, double b, double k) { return fabs(a-b) < k*QUADRATIC_PLANAR::_precision; }
    static double distanceBtw2Pt(const double *a, const double *b) { return sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])); }
  protected:
    ~Node();
  protected:
    bool _isToDel;
    mutable unsigned char _cnt;
    mutable TypeOfLocInPolygon _loc;
    double *_coords;
  };
}

#endif
