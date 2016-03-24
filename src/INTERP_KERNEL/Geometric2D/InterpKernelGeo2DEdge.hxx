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

#ifndef __INTERPKERNELGEO2DEDGE_HXX__
#define __INTERPKERNELGEO2DEDGE_HXX__

#include "INTERPKERNELDefines.hxx"
#include "InterpKernelException.hxx"
#include "InterpKernelGeo2DBounds.hxx"
#include "InterpKernelGeo2DNode.hxx"

#include <iostream>
#include <vector>
#include <list>
#include <map>

namespace INTERP_KERNEL
{
  typedef enum
  {
    SEG         = 1,
    ARC_CIRCLE  = 4,
    ARC_PARABOL = 8
  } TypeOfFunction;

  typedef enum
  {
    CIRCLE  = 0 ,
    PARABOL = 1
  } TypeOfMod4QuadEdge;

  typedef enum
  {
    START       = 5,
    END         = 1,
    INSIDE      = 2,
    OUT_BEFORE  = 3,
    OUT_AFTER   = 4
  } TypeOfLocInEdge; //see Edge::OFFSET_FOR_TYPEOFLOCINEDGE

  typedef enum
  {
    FULL_IN_1    = 1,
    FULL_ON_1    = 4,
    FULL_OUT_1   = 2,
    FULL_UNKNOWN = 3
  } TypeOfEdgeLocInPolygon;

  class INTERPKERNEL_EXPORT MergePoints
  {
  public:
    MergePoints();

    //methods called during intersection edge-edge
    void start1Replaced();
    void end1Replaced();
    void start1OnStart2();
    void start1OnEnd2();
    void end1OnStart2();
    void end1OnEnd2();
    //methods to be called during aggregation
    bool isStart1(unsigned rk) const;
    bool isEnd1(unsigned rk) const;
    bool isStart2(unsigned rk) const;
    bool isEnd2(unsigned rk) const;
    void clear();
    unsigned getNumberOfAssociations() const;
    void updateMergedNodes(int e1Start, int e1End, int e2Start, int e2End, std::map<int,int>& mergedNodes);
  private:
    static void PushInMap(int key, int value, std::map<int,int>& mergedNodes);
  private:
    unsigned _ass1Start1  : 1;
    unsigned _ass1End1    : 1;
    unsigned _ass1Start2  : 1;
    unsigned _ass1End2    : 1;
    unsigned _ass2Start1  : 1;
    unsigned _ass2End1    : 1;
    unsigned _ass2Start2  : 1;
    unsigned _ass2End2    : 1;
  };

  class Edge;
  class ComposedEdge;
  /*!
   * This class is in charge to store an intersection point as result of \b non oververlapping edge intersection.
   * This class manages the cases when intersect element is one of the extrimities of edge1 and/or edge2.
   */
  class INTERPKERNEL_EXPORT IntersectElement
  {
  public:
    IntersectElement(double val1, double val2, bool start1, bool end1, bool start2, bool end2, Node *node, const Edge& e1, const Edge& e2, bool keepOrder);
    IntersectElement(const IntersectElement& other);
    //! The sort operator is done on the edge 1 \b not edge 2.
    bool operator<(const IntersectElement& other) const;
    IntersectElement& operator=(const IntersectElement& other);
    double getVal1() const { return _chararct_val_for_e1; }
    double getVal2() const { return _chararct_val_for_e2; }
    //! idem operator< method except that the orientation is done on edge 2 \b not edge 1.
    bool isLowerOnOther(const IntersectElement& other) const;
    unsigned isOnExtrForAnEdgeAndInForOtherEdge() const;
    void attachLoc() { _node->setLoc(_loc_of_node); }
    bool isOnMergedExtremity() const;
    bool isIncludedByBoth() const;
    void setNode(Node *node) const;
    void performMerging(MergePoints& commonNode) const;
    Node *getNodeOnly() const { return _node; }
    Node *getNodeAndReleaseIt() { Node *tmp=_node; _node=0; return tmp; }
    ~IntersectElement();
  private:
    bool _1S;
    bool _1E;
    bool _2S;
    bool _2E;
    double _chararct_val_for_e1;
    double _chararct_val_for_e2;
    Node *_node;
    TypeOfLocInPolygon _loc_of_node;
    const Edge& _e1;
    const Edge& _e2;
  public:
    static const unsigned LIMIT_ALONE = 22;
    static const unsigned LIMIT_ON = 73;
    static const unsigned NO_LIMIT = 19;
  };

  /*!
   * This abstract interface specifies all the methods to be overloaded of all possibilities edge-intersection.
   */
  class INTERPKERNEL_EXPORT EdgeIntersector
  {
  protected:
    //! All non symetric methods are relative to 'e1'.
    EdgeIntersector(const Edge& e1, const Edge& e2):_e1(e1),_e2(e2) { }
  public:
    virtual ~EdgeIntersector() { }
    virtual bool keepOrder() const = 0;
    virtual bool areColinears() const = 0;
    //!to call only if 'areOverlapped' have been set to true when areOverlappedOrOnlyColinears was called
    virtual bool haveTheySameDirection() const = 0;
    //!to call only if 'areOverlapped' have been set to true when areOverlappedOrOnlyColinears was called
    virtual void getPlacements(Node *start, Node *end, TypeOfLocInEdge& whereStart, TypeOfLocInEdge& whereEnd, MergePoints& commonNode) const = 0;
    //! When true is returned, newNodes should contains at least 1 element. All merging nodes betw _e1 and _e2 extremities must be done.
    bool intersect(const Bounds *whereToFind, std::vector<Node *>& newNodes, bool& order, MergePoints& commonNode);
    //! Should be called only once per association.
    virtual void areOverlappedOrOnlyColinears(const Bounds *whereToFind, bool& obviousNoIntersection, bool& areOverlapped) = 0;
    //! The size of returned vector is equal to number of potential intersections point. The values are so that their are interpretable by virtual Edge::isIn method.
    virtual std::list< IntersectElement > getIntersectionsCharacteristicVal() const = 0;
  protected:
    void obviousCaseForCurvAbscisse(Node *node, TypeOfLocInEdge& where, MergePoints& commonNode, bool& obvious) const;
  protected:
    const Edge& _e1;
    const Edge& _e2;
  };

  class INTERPKERNEL_EXPORT SameTypeEdgeIntersector : public EdgeIntersector
  {
  protected:
    SameTypeEdgeIntersector(const Edge& e1, const Edge& e2):EdgeIntersector(e1,e2) { }
    bool keepOrder() const { return true; }
  };

  class INTERPKERNEL_EXPORT CrossTypeEdgeIntersector : public EdgeIntersector
  {
  protected:
    CrossTypeEdgeIntersector(const Edge& e1, const Edge& e2, bool reverse):EdgeIntersector(e1,e2),_reverse(reverse) { }
    bool keepOrder() const { return _reverse; }
    bool haveTheySameDirection() const { throw Exception("Cross type intersector is not supposed to deal with overlapped in cross type."); }
    const Edge *myE1() { if(_reverse) return &_e1; else return &_e2; }
    const Edge *myE2() { if(_reverse) return &_e2; else return &_e1; }
  protected:
    //! boolean to inform intersector that unsymetrics treatments reverse of e1 and e2 should be done.
    bool _reverse;
  };

  class EdgeLin;
  class EdgeInfLin;
  class EdgeArcCircle;

  /*!
   * Deal with an oriented edge of a polygon.
   * An Edge is defined with a start node, an end node and an equation of 1D curve.
   * All other attributes are mutable because they don't impact these 3 invariant attributes.
   * To be exact start and end nodes can change (address) but their location remain
   * the same (at precision).
   */
  class INTERPKERNEL_EXPORT Edge
  {
  public:
    Edge(Node *start, Node *end, bool direction=true):_cnt(1),_loc(FULL_UNKNOWN) { if(direction) { _start=start; _end=end; } else { _start=end; _end=start; } _start->incrRef(); _end->incrRef(); }
    Edge(double sX, double sY, double eX, double eY);
    TypeOfEdgeLocInPolygon getLoc() const { return _loc; }
    void incrRef() const { _cnt++; }
    bool decrRef();
    void initLocs() const { _loc=FULL_UNKNOWN; _start->initLocs(); _end->initLocs(); }
    void declareOn() const;
    void declareIn() const;
    void declareOut() const;
    void initHitStatus() const { _hit=false; }
    bool getHitStatus() const { return _hit; }
    void hitMeAlone(double xBary, double yBary, double dimChar) { _hit=true; applySimilarity(xBary,yBary,dimChar); }
    void unHitMeAlone(double xBary, double yBary, double dimChar) { _hit=true; unApplySimilarity(xBary,yBary,dimChar); }
    void hitMeAfter(double xBary, double yBary, double dimChar) { if(!_hit) hitMeAlone(xBary,yBary,dimChar); }
    void unHitMeAfter(double xBary, double yBary, double dimChar) { if(!_hit) unHitMeAlone(xBary,yBary,dimChar); }
    const Bounds& getBounds() const { return _bounds; }
    void fillXfigStreamForLoc(std::ostream& stream) const;
    Node *getNode(TypeOfLocInEdge where) const { if(where==START) return _start; else if(where==END) return _end; else return 0; }
    Node *getStartNode() const { return _start; }
    Node *getEndNode() const { return _end; }
    void setEndNodeWithoutChange(Node *newEnd);
    void setStartNodeWithoutChange(Node *newStart);
    bool changeStartNodeWith(Node *otherStartNode) const;
    bool changeStartNodeWithAndKeepTrack(Node *otherStartNode, std::vector<Node *>& track) const;
    bool changeEndNodeWith(Node *otherEndNode) const;
    bool changeEndNodeWithAndKeepTrack(Node *otherEndNode, std::vector<Node *>& track) const;
    void addSubEdgeInVector(Node *start, Node *end, ComposedEdge& vec) const;
    void getNormalVector(double *vectOutput) const;
    static EdgeIntersector *BuildIntersectorWith(const Edge *e1, const Edge *e2);
    static Edge *BuildFromXfigLine(std::istream& str);
    static Edge *BuildEdgeFrom(Node *start, Node *end);
    template<TypeOfMod4QuadEdge type>
    static Edge *BuildEdgeFrom(Node *start, Node *middle, Node *end);
    static Edge *BuildEdgeFrom3Points(const double *start, const double *middle, const double *end);
    virtual void update(Node *m) = 0;
    //! returns area between this and axe Ox delimited along Ox by _start and _end.
    virtual double getAreaOfZone() const = 0;
    //! apply a similiraty transformation on 'this'
    virtual void applySimilarity(double xBary, double yBary, double dimChar);
    //! apply the inverse similiraty transformation on 'this'
    virtual void unApplySimilarity(double xBary, double yBary, double dimChar);
    //! return the length of arc. Value is always > 0. !
    virtual double getCurveLength() const = 0;
    virtual void getBarycenter(double *bary) const = 0;
    virtual void getBarycenterOfZone(double *bary) const = 0;
    //! return the middle of two points
    virtual void getMiddleOfPoints(const double *p1, const double *p2, double *mid) const = 0;
    //! return the middle of two points respecting the orientation defined by this (relevant for arc of circle). By default same as getMiddleOfPoints()
    virtual void getMiddleOfPointsOriented(const double *p1, const double *p2, double *mid) const;
    //! Retrieves a point that is owning to this, well placed for IN/OUT detection of this. Typically midlle of this is returned.
    virtual Node *buildRepresentantOfMySelf() const = 0;
    //! Given a magnitude specified by sub-type returns if in or not. See getCharactValue method.
    virtual bool isIn(double characterVal) const = 0;
    //! With the same magnitude as defined in 'isIn' method perform a compararison. Precondition : val1 and val2 are different and exactly INSIDE this.
    virtual bool isLower(double val1, double val2) const = 0;
    //! node is expected to lay on 'this'. It returns a characteristic magnitude usable by isIn method.
    virtual double getCharactValue(const Node& node) const = 0;
    //! node is expected to lay on 'this'. It returns a characteristic magnitude between 0 and 1.
    virtual double getCharactValueBtw0And1(const Node& node) const = 0;
    //! retrieves the distance to this : The min distance from pt and any point of this.
    virtual double getDistanceToPoint(const double *pt) const = 0;
    //! return if node with coords 'coordOfNode' is on this (with precision).
    virtual bool isNodeLyingOn(const double *coordOfNode) const = 0;
    virtual TypeOfFunction getTypeOfFunc() const = 0;
    virtual void dynCastFunction(const EdgeLin * &seg,
                                 const EdgeArcCircle * &arcSeg) const = 0;
    bool intersectWith(const Edge *other, MergePoints& commonNode,
                       ComposedEdge& outVal1, ComposedEdge& outVal2) const;
    static bool IntersectOverlapped(const Edge *f1, const Edge *f2, EdgeIntersector *intersector, MergePoints& commonNode,
                                    ComposedEdge& outValForF1, ComposedEdge& outValForF2);
    static void Interpolate1DLin(const std::vector<double>& distrib1, const std::vector<double>& distrib2,
                                 std::map<int, std::map<int,double> >& result);
    virtual void dumpInXfigFile(std::ostream& stream, bool direction, int resolution, const Bounds& box) const = 0;
    bool isEqual(const Edge& other) const;
  public:
    bool sortSubNodesAbs(const double *coo, std::vector<int>& subNodes);
    void sortIdsAbs(const std::vector<INTERP_KERNEL::Node *>& addNodes, const std::map<INTERP_KERNEL::Node *, int>& mapp1, const std::map<INTERP_KERNEL::Node *, int>& mapp2, std::vector<int>& edgesThis);
    virtual void fillGlobalInfoAbs(bool direction, const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                                   std::vector<int>& edgesThis, std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int> mapAddCoo) const = 0;
    virtual void fillGlobalInfoAbs2(const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                                    std::vector<int>& edgesOther, std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int>& mapAddCoo) const = 0;
    virtual Edge *buildEdgeLyingOnMe(Node *start, Node *end, bool direction=true) const = 0;
  protected:
    Edge():_cnt(1),_loc(FULL_UNKNOWN),_start(0),_end(0) { }
    virtual ~Edge();
    static int CombineCodes(TypeOfLocInEdge code1, TypeOfLocInEdge code2);
    static bool Intersect(const Edge *f1, const Edge *f2, EdgeIntersector *intersector, const Bounds *whereToFind, MergePoints& commonNode,
                          ComposedEdge& outValForF1, ComposedEdge& outValForF2);
    //! The code 'code' is built by method combineCodes
    static bool SplitOverlappedEdges(const Edge *e1, const Edge *e2, Node *nS, Node *nE, bool direction, int code,
                                     ComposedEdge& outVal1, ComposedEdge& outVal2);
  protected:
    mutable bool _hit;
    mutable unsigned char _cnt;
    mutable TypeOfEdgeLocInPolygon _loc;
    Bounds _bounds;
    Node *_start;
    Node *_end;
  protected:
    //In relation with max possible value of TypeOfLocInEdge.
    static const int OFFSET_FOR_TYPEOFLOCINEDGE = 8;
  };
}

#endif
