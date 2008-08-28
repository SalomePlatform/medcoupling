#ifndef __QUADRATICPOLYGON_HXX__
#define __QUADRATICPOLYGON_HXX__

#include "Geometric2D_defines.hxx"

#include "ComposedEdge.hxx"

#include <list>

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
    ~QuadraticPolygon();
    void closeMe() const;
    void circularPermute();
    //! warning : use it if and only if this is composed of ElementaryEdges only : typical case.
    double getAreaFast() const;
    double getPerimeterFast() const;
    double getHydraulicDiameter() const;
    void dumpInXfigFile(const char *fileName) const;
    void dumpInXfigFileWithOther(const ComposedEdge& other, const char *fileName) const;
    double intersectWith(const QuadraticPolygon& other) const;
    std::vector<QuadraticPolygon *> intersectMySelfWith(const QuadraticPolygon& other) const;
  public://Only public for tests reasons
    void performLocatingOperation(QuadraticPolygon& pol2) const;
    void splitPolygonsEachOther(QuadraticPolygon& pol1, QuadraticPolygon& pol2, int& nbOfSplits) const;
    std::vector<QuadraticPolygon *> buildIntersectionPolygons(const QuadraticPolygon& pol1, const QuadraticPolygon& pol2) const;
    bool amIAChanceToBeCompletedBy(const QuadraticPolygon& pol1Splitted, const QuadraticPolygon& pol2NotSplitted, bool& direction);
  protected:
    std::list<QuadraticPolygon *> zipConsecutiveInSegments() const;
    void dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const;
    void closePolygons(std::list<QuadraticPolygon *>& pol2Zip, const QuadraticPolygon& pol1, std::vector<QuadraticPolygon *>& results) const;
    static void updateNeighbours(const MergePoints& merger, IteratorOnComposedEdge it1, IteratorOnComposedEdge it2,
				 const AbstractEdge *e1, const AbstractEdge *e2);
    std::list<QuadraticPolygon *>::iterator fillAsMuchAsPossibleWith(const QuadraticPolygon& pol1Splitted,
                                                                     std::list<QuadraticPolygon *>::iterator iStart,
                                                                     std::list<QuadraticPolygon *>::iterator iEnd,
                                                                     bool direction);
    static std::list<QuadraticPolygon *>::iterator checkInList(Node *n, std::list<QuadraticPolygon *>::iterator iStart,
                                                               std::list<QuadraticPolygon *>::iterator iEnd);
  };
}

#endif
