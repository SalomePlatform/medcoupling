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
// Author : Anthony Geay (EDF R&D)

#ifndef __DIAMETERCALCULATOR_HXX__
#define __DIAMETERCALCULATOR_HXX__

#include "INTERPKERNELDefines.hxx"

#include "NormalizedGeometricTypes"
#include "MCIdType.hxx"

namespace INTERP_KERNEL
{
  class DiameterCalculator
  {
  public:
    INTERPKERNEL_EXPORT virtual ~DiameterCalculator() { }
    INTERPKERNEL_EXPORT virtual NormalizedCellType getType() const = 0;
    INTERPKERNEL_EXPORT virtual double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const = 0;
    INTERPKERNEL_EXPORT virtual void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const = 0;
    INTERPKERNEL_EXPORT virtual void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const = 0;
    INTERPKERNEL_EXPORT virtual void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const = 0;
  };

  class DiameterCalulatorTRI3S2 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorTRI3S3 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorTRI6S2 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorTRI6S3 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorTRI7S2 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorTRI7S3 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorQUAD4S2 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorQUAD4S3 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorQUAD8S2 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorQUAD8S3 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorQUAD9S2 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorQUAD9S3 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorTETRA4 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorTETRA10 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorHEXA8 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorHEXA20 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorHEXA27 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorPENTA6 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorPENTA15 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorPYRA5 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };

  class DiameterCalulatorPYRA13 : public DiameterCalculator
  {
  public:
    NormalizedCellType getType() const { return TYPE; }
    double computeForOneCell(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr) const { return ComputeForOneCellInternal(bg,endd,coordsPtr); }
    static double ComputeForOneCellInternal(const mcIdType *bg, const mcIdType *endd, const double *coordsPtr);
    void computeForListOfCellIdsUMeshFrmt(const mcIdType *bgIds, const mcIdType *endIds, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeForRangeOfCellIdsUMeshFrmt(mcIdType bgId, mcIdType endId, const mcIdType *indPtr, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
    void computeFor1SGTUMeshFrmt(mcIdType nbOfCells, const mcIdType *connPtr, const double *coordsPtr, double *resPtr) const;
  public:
    static NormalizedCellType TYPE;
  };
}

#endif
