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

#ifndef __INTERPKERNELGEO2DEDGELIN_HXX__
#define __INTERPKERNELGEO2DEDGELIN_HXX__

#include "INTERPKERNELDefines.hxx"
#include "InterpKernelGeo2DEdge.hxx"

namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT SegSegIntersector : SameTypeEdgeIntersector
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

  class INTERPKERNEL_EXPORT EdgeLin : public Edge
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
    void getBarycenterOfZone(double *bary) const;
    void getMiddleOfPoints(const double *p1, const double *p2, double *mid) const;
    bool isIn(double characterVal) const;
    Node *buildRepresentantOfMySelf() const;
    double getCharactValue(const Node& node) const;
    double getCharactValueBtw0And1(const Node& node) const;
    double getDistanceToPoint(const double *pt) const;
    bool isNodeLyingOn(const double *coordOfNode) const;
    bool isLower(double val1, double val2) const { return val1<val2; }
    double getCharactValueEng(const double *node) const;
    bool doIHaveSameDirectionAs(const Edge& other) const;
    void dynCastFunction(const EdgeLin * &seg,
                         const EdgeArcCircle * &arcSeg) const { seg=this; }
  protected:
    EdgeLin() { }
    void updateBounds();
    Edge *buildEdgeLyingOnMe(Node *start, Node *end, bool direction) const;
    void fillGlobalInfoAbs(bool direction, const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                           std::vector<int>& edgesThis, std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int> mapAddCoo) const;
    void fillGlobalInfoAbs2(const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                            std::vector<int>& edgesOther, std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int>& mapAddCoo) const;
  };
}

#endif
