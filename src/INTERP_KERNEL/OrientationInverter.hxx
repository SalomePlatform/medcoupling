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
// Author : Anthony Geay (EDF R&D)

#ifndef __ORIENTATIONINVERTER_HXX__
#define __ORIENTATIONINVERTER_HXX__

#include "INTERPKERNELDefines.hxx"
#include "NormalizedGeometricTypes"

namespace INTERP_KERNEL
{
  class OrientationInverter
  {
  public:
    INTERPKERNEL_EXPORT static OrientationInverter *BuildInstanceFrom(NormalizedCellType gt);
    INTERPKERNEL_EXPORT virtual ~OrientationInverter() { }
    INTERPKERNEL_EXPORT virtual void operate(int *beginPt, int *endPt) const = 0;
  };

  class OrientationInverterChecker : public OrientationInverter
  {
  public:
    OrientationInverterChecker(unsigned nbNodes):_nb_nodes(nbNodes) { }
    void operate(int *beginPt, int *endPt) const { check(beginPt,endPt); operateAndShutUp(beginPt); }
    virtual void operateAndShutUp(int *beginPt) const = 0;
  protected:
    unsigned getNbNodes() const { return _nb_nodes; }
  private:
    void check(int *beginPt, int *endPt) const;
  private:
    unsigned _nb_nodes;
  };

  class OrientationInverterSEG2 : public OrientationInverterChecker
  {
  public:
    OrientationInverterSEG2():OrientationInverterChecker(2u) { }
    void operateAndShutUp(int *beginPt) const;
  };

  class OrientationInverterSEG3 : public OrientationInverterChecker
  {
  public:
    OrientationInverterSEG3():OrientationInverterChecker(3u) { }
    void operateAndShutUp(int *beginPt) const;
  };

  class OrientationInverter2DLinear : public OrientationInverterChecker
  {
  public:
    OrientationInverter2DLinear(unsigned nbNodes):OrientationInverterChecker(nbNodes) { }
    void operateAndShutUp(int *beginPt) const;
  };

  class OrientationInverter2DQuadratic : public OrientationInverterChecker
  {
  public:
    OrientationInverter2DQuadratic(unsigned nbNodes):OrientationInverterChecker(nbNodes) { }
    void operateAndShutUp(int *beginPt) const;
  };

  class OrientationInverterPolygon : public OrientationInverter
  {
  public:
    void operate(int *beginPt, int *endPt) const;
  };

  class OrientationInverterQPolygon : public OrientationInverter
  {
  public:
    void operate(int *beginPt, int *endPt) const;
  };

  class OrientationInverterTetra4 : public OrientationInverterChecker
  {
  public:
    OrientationInverterTetra4():OrientationInverterChecker(4u) { }
    void operateAndShutUp(int *beginPt) const;
  };

  class OrientationInverterTetra10 : public OrientationInverterChecker
  {
  public:
    OrientationInverterTetra10():OrientationInverterChecker(10u) { }
    void operateAndShutUp(int *beginPt) const;
  };

  class OrientationInverterPyra5 : public OrientationInverterChecker
  {
  public:
    OrientationInverterPyra5():OrientationInverterChecker(5u) { }
    void operateAndShutUp(int *beginPt) const;
  };

  class OrientationInverterPyra13 : public OrientationInverterChecker
  {
  public:
    OrientationInverterPyra13():OrientationInverterChecker(13u) { }
    void operateAndShutUp(int *beginPt) const;
  };

  class OrientationInverter3DExtrusionLinear : public OrientationInverterChecker
  {
  public:
    OrientationInverter3DExtrusionLinear(unsigned nbNodes):OrientationInverterChecker(nbNodes) { }
    void operateAndShutUp(int *beginPt) const;
  };

  class OrientationInverter3DExtrusionQuadratic : public OrientationInverterChecker
  {
  public:
    OrientationInverter3DExtrusionQuadratic(unsigned nbNodes):OrientationInverterChecker(nbNodes) { }
    void operateAndShutUp(int *beginPt) const;
  };
}

#endif
