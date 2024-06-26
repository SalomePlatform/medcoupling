// Copyright (C) 2007-2024  CEA, EDF
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

#ifndef __INTERPKERNELGEO2DQUADRATICPOLYGON_HXX__
#define __INTERPKERNELGEO2DQUADRATICPOLYGON_HXX__

#include "INTERPKERNELDefines.hxx"

#include "InterpKernelGeo2DComposedEdge.hxx"
#include "InterpKernelGeo2DAbstractEdge.hxx"
#include "InterpKernelGeo2DElementaryEdge.hxx"

#include <list>
#include <vector>

namespace INTERP_KERNEL
{
  class Edge;
  class MergePoints;

  enum NodeUsage { USAGE_UNKNOWN, USAGE_LINEAR, USAGE_QUADRATIC_ONLY };
  typedef std::pair<INTERP_KERNEL::Node *,NodeUsage> NodeWithUsage;

  /**
   * A set of quadratic or linear edges, not necessarily connected to form a closed polygon.
   * Some methods however requires a closed form.
   * Class ComposedEdge focuses more on connectivity aspect.
   */
  class QuadraticPolygon : public ComposedEdge
  {
  public:
    INTERPKERNEL_EXPORT QuadraticPolygon() { }
    INTERPKERNEL_EXPORT QuadraticPolygon(const QuadraticPolygon& other):ComposedEdge(other) { }
    INTERPKERNEL_EXPORT QuadraticPolygon(const char *fileName);
    INTERPKERNEL_EXPORT static QuadraticPolygon *BuildLinearPolygon(std::vector<Node *>& nodes);
    INTERPKERNEL_EXPORT static QuadraticPolygon *BuildArcCirclePolygon(std::vector<Node *>& nodes);
    INTERPKERNEL_EXPORT static Edge *BuildLinearEdge(std::vector<Node *>& nodes);
    INTERPKERNEL_EXPORT static Edge *BuildArcCircleEdge(std::vector<Node *>& nodes);
    INTERPKERNEL_EXPORT static void BuildDbgFile(const std::vector<Node *>& nodes, const char *fileName);
    INTERPKERNEL_EXPORT ~QuadraticPolygon();
    INTERPKERNEL_EXPORT void closeMe() const;
    INTERPKERNEL_EXPORT void circularPermute();
    INTERPKERNEL_EXPORT bool isButterflyAbs();
    INTERPKERNEL_EXPORT bool isButterfly() const;
    INTERPKERNEL_EXPORT void dumpInXfigFile(const char *fileName) const;
    INTERPKERNEL_EXPORT void dumpInXfigFileWithOther(const ComposedEdge& other, const char *fileName) const;
    //! Before intersecting as intersectWith a normalization is done.
    INTERPKERNEL_EXPORT double intersectWithAbs(QuadraticPolygon& other);
    INTERPKERNEL_EXPORT double intersectWithAbs1D(QuadraticPolygon& other, bool& isColinear);
    //! Before intersecting as intersectWith a normalization is done.
    INTERPKERNEL_EXPORT double intersectWithAbs(QuadraticPolygon& other, double* barycenter);
    INTERPKERNEL_EXPORT void splitAbs(QuadraticPolygon& other, const std::map<INTERP_KERNEL::Node *,mcIdType>& mapThis, const std::map<INTERP_KERNEL::Node *,mcIdType>& mapOther, mcIdType offset1, mcIdType offset2, const std::vector<mcIdType>& otherEdgeIds,
                                      std::vector<mcIdType>& edgesThis, mcIdType cellIdThis, std::vector< std::vector<mcIdType> >& edgesInOtherColinearWithThis, std::vector< std::vector<mcIdType> >& subDivOther, std::vector<double>& addCoo, std::map<mcIdType,mcIdType>& mergedNodes);
    INTERPKERNEL_EXPORT void buildFromCrudeDataArray(const std::map<mcIdType,INTERP_KERNEL::Node *>& mapp, bool isQuad, const mcIdType *nodalBg, const double *coords,
                                                     const mcIdType *descBg, const mcIdType *descEnd, const std::vector<std::vector<mcIdType> >& intersectEdges);
    INTERPKERNEL_EXPORT void buildFromCrudeDataArray2(const std::map<mcIdType,INTERP_KERNEL::Node *>& mapp, bool isQuad, const mcIdType *nodalBg, const double *coords, const mcIdType *descBg, const mcIdType *descEnd, const std::vector<std::vector<mcIdType> >& intersectEdges,
                                                      const INTERP_KERNEL::QuadraticPolygon& pol1, const mcIdType *descBg1, const mcIdType *descEnd1, const std::vector<std::vector<mcIdType> >& intersectEdges1,
                                                      const std::vector< std::vector<mcIdType> >& colinear1,
                                                      std::map<mcIdType,std::vector<INTERP_KERNEL::ElementaryEdge *> >& alreadyExistingIn2);
    INTERPKERNEL_EXPORT void updateLocOfEdgeFromCrudeDataArray2(const mcIdType *descBg, const mcIdType *descEnd, const std::vector<std::vector<mcIdType> >& intersectEdges, const INTERP_KERNEL::QuadraticPolygon& pol1, const mcIdType *descBg1, const mcIdType *descEnd1, const std::vector<std::vector<mcIdType> >& intersectEdges1, const std::vector< std::vector<mcIdType> >& colinear1) const;
    INTERPKERNEL_EXPORT void appendEdgeFromCrudeDataArray(std::size_t edgeId, const std::map<mcIdType,INTERP_KERNEL::Node *>& mapp, bool isQuad, const mcIdType *nodalBg, const double *coords,
                                                          const mcIdType *descBg,  const mcIdType *descEnd, const std::vector<std::vector<mcIdType> >& intersectEdges);
    INTERPKERNEL_EXPORT void appendSubEdgeFromCrudeDataArray(Edge *baseEdge, std::size_t j, bool direct, mcIdType edgeId, const std::vector<mcIdType>& subEdge, const std::map<mcIdType,INTERP_KERNEL::Node *>& mapp);
    INTERPKERNEL_EXPORT void appendCrudeData(const std::map<INTERP_KERNEL::Node *,mcIdType>& mapp, double xBary, double yBary, double fact, mcIdType offset, std::vector<double>& addCoordsQuadratic, std::vector<mcIdType>& conn, std::vector<mcIdType>& connI) const;
    INTERPKERNEL_EXPORT void buildPartitionsAbs(QuadraticPolygon& other, std::set<INTERP_KERNEL::Edge *>& edgesThis, std::set<INTERP_KERNEL::Edge *>& edgesBoundaryOther, const std::map<INTERP_KERNEL::Node *,mcIdType>& mapp, mcIdType idThis, mcIdType idOther, mcIdType offset,
                                                std::vector<double>& addCoordsQuadratic, std::vector<mcIdType>& conn, std::vector<mcIdType>& connI, std::vector<mcIdType>& nb1, std::vector<mcIdType>& nb2);
    //
    INTERPKERNEL_EXPORT double intersectWith(const QuadraticPolygon& other) const;
    INTERPKERNEL_EXPORT double intersectWith(const QuadraticPolygon& other, double* barycenter) const;
    INTERPKERNEL_EXPORT std::vector<QuadraticPolygon *> intersectMySelfWith(const QuadraticPolygon& other) const;
    INTERPKERNEL_EXPORT void intersectForPerimeter(const QuadraticPolygon& other, double& perimeterThisPart, double& perimeterOtherPart, double& perimeterCommonPart) const;
    INTERPKERNEL_EXPORT void intersectForPerimeterAdvanced(const QuadraticPolygon& other, std::vector< double >& polThis, std::vector< double >& polOther) const;
    INTERPKERNEL_EXPORT void intersectForPoint(const QuadraticPolygon& other, std::vector< int >& numberOfCreatedPointsPerEdge) const;
  public://Only public for tests reasons
    INTERPKERNEL_EXPORT void performLocatingOperation(QuadraticPolygon& pol2) const;
    INTERPKERNEL_EXPORT void performLocatingOperationSlow(QuadraticPolygon& pol2) const;
    INTERPKERNEL_EXPORT static void SplitPolygonsEachOther(QuadraticPolygon& pol1, QuadraticPolygon& pol2, int& nbOfSplits);
    INTERPKERNEL_EXPORT std::vector<QuadraticPolygon *> buildIntersectionPolygons(const QuadraticPolygon& pol1, const QuadraticPolygon& pol2) const;
    INTERPKERNEL_EXPORT bool haveIAChanceToBeCompletedBy(const QuadraticPolygon& pol1NotSplitted, const QuadraticPolygon& pol2Splitted,
                                                         bool& direction, bool& needCleaning) const;
    INTERPKERNEL_EXPORT static void ComputeResidual(const QuadraticPolygon& pol1, const std::set<Edge *>& notUsedInPol1, const std::set<Edge *>& edgesInPol2OnBoundary, const std::map<INTERP_KERNEL::Node *,mcIdType>& mapp, mcIdType offset, mcIdType idThis,
                                                    std::vector<double>& addCoordsQuadratic, std::vector<mcIdType>& conn, std::vector<mcIdType>& connI, std::vector<mcIdType>& nb1, std::vector<mcIdType>& nb2);
    INTERPKERNEL_EXPORT void cleanDegeneratedConsecutiveEdges();
  protected:
    std::list<QuadraticPolygon *> zipConsecutiveInSegments() const;
    void dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const;
    static void ClosePolygons(std::list<QuadraticPolygon *>& pol1Zip, const QuadraticPolygon& pol1, const QuadraticPolygon& pol2,
                              std::vector<QuadraticPolygon *>& results);
    template<class EDGES>
    static void UpdateNeighbours(const MergePoints& merger, IteratorOnComposedEdge it1, IteratorOnComposedEdge it2,
                                 const EDGES *e1, const EDGES *e2);
    std::list<QuadraticPolygon *>::iterator fillAsMuchAsPossibleWith(const QuadraticPolygon& pol2Splitted,
                                                                     std::list<QuadraticPolygon *>::iterator iStart,
                                                                     std::list<QuadraticPolygon *>::iterator iEnd,
                                                                     bool direction);
    static std::list<QuadraticPolygon *>::iterator CheckInList(Node *n, std::list<QuadraticPolygon *>::iterator iStart,
                                                               std::list<QuadraticPolygon *>::iterator iEnd);
  };
}

namespace INTERP_KERNEL
{
  template<class EDGES>
  void QuadraticPolygon::UpdateNeighbours(const MergePoints& merger, IteratorOnComposedEdge it1, IteratorOnComposedEdge it2,
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
