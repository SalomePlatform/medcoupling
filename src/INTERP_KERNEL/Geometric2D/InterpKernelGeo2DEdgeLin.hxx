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

#ifndef __INTERPKERNELGEO2DEDGELIN_HXX__
#define __INTERPKERNELGEO2DEDGELIN_HXX__

#include "INTERPKERNELDefines.hxx"
#include "InterpKernelGeo2DBounds.hxx"
#include "InterpKernelGeo2DEdge.hxx"
#include "InterpKernelGeo2DNode.hxx"
#include <list>
#include <istream>
#include <ostream>

namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT SegSegIntersector : SameTypeEdgeIntersector
  {
    friend class Edge;
  public:
    SegSegIntersector(const EdgeLin& e1, const EdgeLin& e2);
    bool areColinears() const override;
    bool haveTheySameDirection() const override;
    void getPlacements(Node *start, Node *end, TypeOfLocInEdge& whereStart, TypeOfLocInEdge& whereEnd, MergePoints& commonNode) const override;
    void areOverlappedOrOnlyColinears(bool& obviousNoIntersection, bool& areOverlapped) override;
    std::list< IntersectElement > getIntersectionsCharacteristicVal() const override;
  private:
    void getCurveAbscisse(Node *node, TypeOfLocInEdge& where, MergePoints& commonNode) const;
  private:
    //! index on which all single index op will be performed. Filled in case colinearity is equal to true.
    int _ind;
    double _col[2];
    double _matrix[4];               //SPACEDIM*SPACEDIM  = [e1_x, e1_y, e2_x, e2_y]
    double _determinant;
  };

  class INTERPKERNEL_EXPORT EdgeLin : public Edge
  {
    friend class SegSegIntersector;
  public:
    EdgeLin(std::istream& lineInXfig);
    EdgeLin(Node *start, Node *end, bool direction=true);
    EdgeLin(double sX, double sY, double eX, double eY);
    ~EdgeLin() override;
    TypeOfFunction getTypeOfFunc() const override { return SEG; }
    void dumpInXfigFile(std::ostream& stream, bool direction, int resolution, const Bounds& box) const override;
    void update(Node *m) override;
    double getNormSq() const;
    double getAreaOfZone() const override;
    double getCurveLength() const override;
    void getBarycenter(double *bary) const override;
    void getBarycenterOfZone(double *bary) const override;
    void getMiddleOfPoints(const double *p1, const double *p2, double *mid) const override;
    bool isIn(double characterVal) const override;
    Node *buildRepresentantOfMySelf() const override;
    double getCharactValue(const Node& node) const override;
    double getCharactValueBtw0And1(const Node& node) const override;
    double getDistanceToPoint(const double *pt) const override;
    bool isNodeLyingOn(const double *coordOfNode) const override;
    bool isLower(double val1, double val2) const override { return val1<val2; }
    double getCharactValueEng(const double *node) const;
    bool doIHaveSameDirectionAs(const Edge& other) const;
    void dynCastFunction(const EdgeLin * &seg,
                         const EdgeArcCircle * & /*arcSeg*/) const override { seg=this; }
  protected:
    EdgeLin() { }
    void updateBounds();
    Edge *buildEdgeLyingOnMe(Node *start, Node *end, bool direction) const override;
  };
}

#endif
