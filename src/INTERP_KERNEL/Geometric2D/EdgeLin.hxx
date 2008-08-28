#ifndef __EDGELIN_HXX__
#define __EDGELIN_HXX__

#include "Geometric2D_defines.hxx"

#include "Edge.hxx"

namespace INTERP_KERNEL
{
  class GEOMETRIC2D_EXPORT SegSegIntersector : SameTypeIntersector
    {
      friend class Edge;
    public:
      SegSegIntersector(const EdgeLin& e1, const EdgeLin& e2);
      bool areColinears() const;
      bool haveTheySameDirection() const;
      void getPlacements(Node *start, Node *end, TypeOfLocInEdge& whereStart, TypeOfLocInEdge& whereEnd, MergePoints& commonNode) const;
      void areOverlappedOrOnlyColinears(const Bounds *whereToFind, bool& obviousNoIntersection, bool& areOverlapped);
      std::list< IntersectElement > getIntersectionsCharacteristicVal() const;
    private:
      void getCurveAbscisse(Node *node, TypeOfLocInEdge& where, MergePoints& commonNode) const;
    private:
      //! index on which all single index op will be performed. Filled in case colinearity is equal to true.
      int _ind;
      double _col[2];
      double _matrix[4];//SPACEDIM*SPACEDIM
    };

  class GEOMETRIC2D_EXPORT EdgeLin : public Edge
  {
    friend class SegSegIntersector;
  public:
    EdgeLin(std::istream& lineInXfig);
    EdgeLin(Node *start, Node *end, bool direction=true);
    EdgeLin(double sX, double sY, double eX, double eY);
    ~EdgeLin();
    TypeOfFunction getTypeOfFunc() const { return SEG; }
    void dumpInXfigFile(std::ostream& stream, bool direction, int resolution, const Bounds& box) const;
    void update(Node *m);
    double getNormSq() const;
    double getAreaOfZone() const;
    double getCurveLength() const;
    void getBarycenter(double *bary) const;
    bool isIn(double characterVal) const;
    Node *buildRepresentantOfMySelf() const;
    double getCharactValue(const Node& node) const;
    bool isLower(double val1, double val2) const { return val1<val2; }
    bool doIHaveSameDirectionAs(const Edge& other) const;
    Edge *buildEdgeLyingOnMe(Node *start, Node *end, bool direction=true) const;
    void dynCastFunction(const EdgeLin * &seg,
                         const EdgeArcCircle * &arcSeg) const { seg=this; }
  protected:
    EdgeLin() { }
    void updateBounds();
  };
}

#endif
