#ifndef __EDGE_HXX__
#define __EDGE_HXX__

#include "InterpolationUtils.hxx"
#include "ComposedEdge.hxx"
#include "Bounds.hxx"
#include "Node.hxx"

#include <iostream>
#include <vector>
#include <list>

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

  class MergePoints
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

  class IntersectElement
  {
  public:
    IntersectElement(double val1, double val2, bool start1, bool end1, bool start2, bool end2, Node *node, const Edge& e1, const Edge& e2, bool keepOrder);
    IntersectElement(const IntersectElement& other);
    bool operator<(const IntersectElement& other) const;
    IntersectElement& operator=(const IntersectElement& other);
    double getVal1() const { return _chararctValForE1; }
    double getVal2() const { return _chararctValForE2; }
    bool isLowerOnOther(const IntersectElement& other) const;
    unsigned isOnExtrForAnEdgeAndInForOtherEdge() const;
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
    double _chararctValForE1;
    double _chararctValForE2;
    Node *_node;
    const Edge& _e1;
    const Edge& _e2;
  public:
    static const unsigned LIMIT_ALONE = 22;
    static const unsigned LIMIT_ON = 73;
    static const unsigned NO_LIMIT = 19;
  };

  class Intersector
  {
  protected:
    //! All non symetric methods are relative to 'e1'.
    Intersector(const Edge& e1, const Edge& e2):_e1(e1),_e2(e2) { }
  public:
    virtual ~Intersector() { }
    virtual bool keepOrder() const = 0;
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

  class SameTypeIntersector : public Intersector
  {
  protected:
    SameTypeIntersector(const Edge& e1, const Edge& e2):Intersector(e1,e2) { }
    bool keepOrder() const { return true; }
  };

  class CrossTypeIntersector : public Intersector
  {
  protected:
    CrossTypeIntersector(const Edge& e1, const Edge& e2, bool reverse):Intersector(e1,e2),_reverse(reverse) { }
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
   */
  class Edge
  {
  public:
    Edge(Node *start, Node *end, bool direction=true):_cnt(1),_loc(FULL_UNKNOWN) { if(direction) { _start=start; _end=end; } else { _start=end; _end=start; } _start->incrRef(); _end->incrRef(); }
    Edge(double sX, double sY, double eX, double eY);
    TypeOfEdgeLocInPolygon getLoc() const { return _loc; }
    void incrRef() const { _cnt++; }
    bool decrRef();
    void declareOn() const;
    void declareIn() const;
    void declareOut() const;
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
    static Intersector *buildIntersectorWith(const Edge *e1, const Edge *e2);
    static Edge *buildFromXfigLine(std::istream& str);
    static Edge *buildEdgeFrom(Node *start, Node *end);
    template<TypeOfMod4QuadEdge type>
    static Edge *buildEdgeFrom(Node *start, Node *middle, Node *end);
    virtual void update(Node *m) = 0;
    //! returns area between this and axe Ox delimited along Ox by _start and _end.
    virtual double getAreaOfZone() const = 0;
    virtual double getCurveLength() const = 0;
    virtual void getBarycenter(double *bary) const = 0;
    //! Retrieves a point that is owning to this, well placed for IN/OUT detection of this. Typically midlle of this is returned.
    virtual Node *buildRepresentantOfMySelf() const = 0;
    //! Given a magnitude specified by sub-type returns if in or not. See getCharactValue method.
    virtual bool isIn(double characterVal) const = 0;
    //! With the same magnitude as defined in 'isIn' method perform a compararison. Precondition : val1 and val2 are different and exactly INSIDE this.
    virtual bool isLower(double val1, double val2) const = 0;
    //! node is expected to lay on 'this'. It returns a characteristic magnitude usable by isIn method.
    virtual double getCharactValue(const Node& node) const = 0;
    virtual TypeOfFunction getTypeOfFunc() const = 0;
    virtual Edge *buildEdgeLyingOnMe(Node *start, Node *end, bool direction=true) const = 0;
    virtual void dynCastFunction(const EdgeLin * &seg,
                                 const EdgeArcCircle * &arcSeg) const = 0;
    bool intersectWith(const Edge *other, MergePoints& commonNode,
                       ComposedEdge& outVal1, ComposedEdge& outVal2) const;
    virtual void dumpInXfigFile(std::ostream& stream, bool direction, int resolution, const Bounds& box) const = 0;
  protected:
    Edge():_cnt(1),_loc(FULL_UNKNOWN),_start(0),_end(0) { }
    virtual ~Edge();
    static int combineCodes(TypeOfLocInEdge code1, TypeOfLocInEdge code2);
    static bool intersect(const Edge *f1, const Edge *f2, Intersector *intersector, const Bounds *whereToFind, MergePoints& commonNode,
                          ComposedEdge& outValForF1, ComposedEdge& outValForF2);
    //! The code 'code' is built by method combineCodes
    static bool splitOverlappedEdges(const Edge *e1, const Edge *e2, Node *nS, Node *nE, bool direction, int code,
                                     ComposedEdge& outVal1, ComposedEdge& outVal2);
  protected:
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
