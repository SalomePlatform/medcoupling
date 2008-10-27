#ifndef __QUADRATICPOLYGON_HXX__
#define __QUADRATICPOLYGON_HXX__

#include "Geometric2D_defines.hxx"

#include "ComposedEdge.hxx"
#include "AbstractEdge.hxx"
#include "ElementaryEdge.hxx"

#include <list>
#include <vector>

namespace INTERP_KERNEL
{
  class Edge;
  class MergePoints;

  class GEOMETRIC2D_EXPORT QuadraticPolygon : public ComposedEdge
  {
  public:
    QuadraticPolygon() { }
    QuadraticPolygon(const QuadraticPolygon& other):ComposedEdge(other) { }
    QuadraticPolygon(const char *fileName);
    static QuadraticPolygon *buildLinearPolygon(std::vector<Node *>& nodes);
    static QuadraticPolygon *buildArcCirclePolygon(std::vector<Node *>& nodes);
    static void buildDbgFile(const std::vector<Node *>& nodes, const char *fileName);
    ~QuadraticPolygon();
    void closeMe() const;
    void circularPermute();
    void dumpInXfigFile(const char *fileName) const;
    void dumpInXfigFileWithOther(const ComposedEdge& other, const char *fileName) const;
    //! Before intersecting as intersectWith a normalization is done.
    double intersectWithAbs(QuadraticPolygon& other);
    double intersectWith(const QuadraticPolygon& other) const;
    std::vector<QuadraticPolygon *> intersectMySelfWith(const QuadraticPolygon& other) const;
    void intersectForPerimeter(const QuadraticPolygon& other, double& perimeterThisPart, double& perimeterOtherPart, double& perimeterCommonPart, double& area) const;
    double intersectForPerimeterAdvanced(const QuadraticPolygon& other, std::vector< double >& polThis, std::vector< double >& polOther, double& area) const;
    void intersectForPoint(const QuadraticPolygon& other, std::vector< int >& numberOfCreatedPointsPerEdge) const;
  public://Only public for tests reasons
    void performLocatingOperation(QuadraticPolygon& pol2) const;
    void splitPolygonsEachOther(QuadraticPolygon& pol1, QuadraticPolygon& pol2, int& nbOfSplits) const;
    std::vector<QuadraticPolygon *> buildIntersectionPolygons(const QuadraticPolygon& pol1, const QuadraticPolygon& pol2) const;
    bool amIAChanceToBeCompletedBy(const QuadraticPolygon& pol1Splitted, const QuadraticPolygon& pol2NotSplitted, bool& direction);
  protected:
    std::list<QuadraticPolygon *> zipConsecutiveInSegments() const;
    void dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const;
    void closePolygons(std::list<QuadraticPolygon *>& pol2Zip, const QuadraticPolygon& pol1, std::vector<QuadraticPolygon *>& results) const;
    template<class EDGES>
    static void updateNeighbours(const MergePoints& merger, IteratorOnComposedEdge it1, IteratorOnComposedEdge it2,
				 const EDGES *e1, const EDGES *e2);
    std::list<QuadraticPolygon *>::iterator fillAsMuchAsPossibleWith(const QuadraticPolygon& pol1Splitted,
                                                                     std::list<QuadraticPolygon *>::iterator iStart,
                                                                     std::list<QuadraticPolygon *>::iterator iEnd,
                                                                     bool direction);
    static std::list<QuadraticPolygon *>::iterator checkInList(Node *n, std::list<QuadraticPolygon *>::iterator iStart,
                                                               std::list<QuadraticPolygon *>::iterator iEnd);
  };
}

namespace INTERP_KERNEL
{
  template<class EDGES>
  void QuadraticPolygon::updateNeighbours(const MergePoints& merger, IteratorOnComposedEdge it1, IteratorOnComposedEdge it2,
                                          const EDGES *e1, const EDGES *e2)
  {
    it1.previousLoop(); it2.previousLoop();
    ElementaryEdge *curE1=it1.current(); ElementaryEdge *curE2=it2.current();
    curE1->changeEndNodeWith(e1->getStartNode()); curE2->changeEndNodeWith(e2->getStartNode());
    it1.nextLoop(); it1.nextLoop(); it2.nextLoop(); it2.nextLoop();
    curE1->changeStartNodeWith(e1->getEndNode()); curE2->changeStartNodeWith(e2->getEndNode());
  }
}

#endif
