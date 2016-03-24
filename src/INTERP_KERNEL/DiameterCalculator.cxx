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

#include "DiameterCalculator.hxx"
#include "InterpKernelException.hxx"
#include "CellModel.hxx"

#include <algorithm>
#include <sstream>
#include <cmath>

using namespace INTERP_KERNEL;

NormalizedCellType DiameterCalulatorTRI3S2::TYPE=NORM_TRI3;

NormalizedCellType DiameterCalulatorTRI3S3::TYPE=NORM_TRI3;

NormalizedCellType DiameterCalulatorTRI6S2::TYPE=NORM_TRI6;

NormalizedCellType DiameterCalulatorTRI6S3::TYPE=NORM_TRI6;

NormalizedCellType DiameterCalulatorTRI7S2::TYPE=NORM_TRI7;

NormalizedCellType DiameterCalulatorTRI7S3::TYPE=NORM_TRI7;

NormalizedCellType DiameterCalulatorQUAD4S2::TYPE=NORM_QUAD4;

NormalizedCellType DiameterCalulatorQUAD4S3::TYPE=NORM_QUAD4;

NormalizedCellType DiameterCalulatorQUAD8S2::TYPE=NORM_QUAD8;

NormalizedCellType DiameterCalulatorQUAD8S3::TYPE=NORM_QUAD8;

NormalizedCellType DiameterCalulatorQUAD9S2::TYPE=NORM_QUAD9;

NormalizedCellType DiameterCalulatorQUAD9S3::TYPE=NORM_QUAD9;

NormalizedCellType DiameterCalulatorTETRA4::TYPE=NORM_TETRA4;

NormalizedCellType DiameterCalulatorTETRA10::TYPE=NORM_TETRA10;

NormalizedCellType DiameterCalulatorHEXA8::TYPE=NORM_HEXA8;

NormalizedCellType DiameterCalulatorHEXA20::TYPE=NORM_HEXA20;

NormalizedCellType DiameterCalulatorHEXA27::TYPE=NORM_HEXA27;

NormalizedCellType DiameterCalulatorPENTA6::TYPE=NORM_PENTA6;

NormalizedCellType DiameterCalulatorPENTA15::TYPE=NORM_PENTA15;

NormalizedCellType DiameterCalulatorPYRA5::TYPE=NORM_PYRA5;

NormalizedCellType DiameterCalulatorPYRA13::TYPE=NORM_PYRA13;

inline double SqNormV2(const double tab[2])
{
  return tab[0]*tab[0]+tab[1]*tab[1];
}

inline void DiffV2(const double a[2], const double b[2], double c[2])
{
  c[0]=a[0]-b[0];
  c[1]=a[1]-b[1];
}

inline double SqNormV3(const double tab[3])
{
  return tab[0]*tab[0]+tab[1]*tab[1]+tab[2]*tab[2];
}

inline void DiffV3(const double a[3], const double b[3], double c[3])
{
  c[0]=a[0]-b[0];
  c[1]=a[1]-b[1];
  c[2]=a[2]-b[2];
}

template<class Evaluator>
void ComputeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr)
{
  Evaluator evtor;
  NormalizedCellType ct(Evaluator::TYPE);
  int cti((int) ct);
  for(const int *it=bgIds;it!=endIds;it++)
    {
      int offset(indPtr[*it]);
      if(connPtr[offset]==cti)
        resPtr[*it]=evtor.ComputeForOneCellInternal(connPtr+offset+1,connPtr+indPtr[(*it)+1],coordsPtr);
      else
        {
          std::ostringstream oss; oss << "DiameterCalculator::computeForListOfCellIdsUMeshFrmt : invalid nodal connectivity format at cell # " << *it <<  " !";
          throw Exception(oss.str().c_str());
        }
    }
}

template<class Evaluator>
void ComputeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr)
{
  Evaluator evtor;
  NormalizedCellType ct(Evaluator::TYPE);
  int cti((int) ct);
  for(int it=bgId;it<endId;it++)
    {
      int offset(indPtr[it]);
      if(connPtr[offset]==cti)
        resPtr[it]=evtor.ComputeForOneCellInternal(connPtr+offset+1,connPtr+indPtr[it+1],coordsPtr);
      else
        {
          std::ostringstream oss; oss << "DiameterCalculator::computeForListOfCellIdsUMeshFrmt : invalid nodal connectivity format at cell # " << it <<  " !";
          throw Exception(oss.str().c_str());
        }
    }
}

template<class Evaluator>
void ComputeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr)
{
  Evaluator evtor;
  NormalizedCellType ct(Evaluator::TYPE);
  const CellModel& cm(CellModel::GetCellModel(ct));
  unsigned nbNodes(cm.getNumberOfNodes());
  const int *ptr(connPtr);
  for(int i=0;i<nbOfCells;i++,ptr+=nbNodes,resPtr++)
    *resPtr=evtor.ComputeForOneCellInternal(ptr,ptr+nbNodes,coordsPtr);
}

//=================================================================

double DiameterCalulatorTRI3S2::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==3)
    {
      const double *a(coordsPtr+2*bg[0]),*b(coordsPtr+2*bg[1]),*c(coordsPtr+2*bg[2]);
      double l0[2],l1[2],l2[2];
      DiffV2(a,b,l0);
      DiffV2(a,c,l1);
      DiffV2(b,c,l2);
      double res(std::max(SqNormV2(l0),SqNormV2(l1)));
      res=std::max(res,SqNormV2(l2));
      return std::sqrt(res);
    }
  else
    throw Exception("DiameterCalulatorTRI3S2::ComputeForOneCellInternal : input connectivity must be of size 3 !");
}

void DiameterCalulatorTRI3S2::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorTRI3S2>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTRI3S2::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorTRI3S2>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTRI3S2::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorTRI3S2>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorTRI3S3::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==3)
    {
      const double *a(coordsPtr+3*bg[0]),*b(coordsPtr+3*bg[1]),*c(coordsPtr+3*bg[2]);
      double l0[3],l1[3],l2[3];
      DiffV3(a,b,l0);
      DiffV3(a,c,l1);
      DiffV3(b,c,l2);
      double res(std::max(SqNormV3(l0),SqNormV3(l1)));
      res=std::max(res,SqNormV3(l2));
      return std::sqrt(res);
    }
  else
    throw Exception("DiameterCalulatorTRI3S2::ComputeForOneCellInternal : input connectivity must be of size 3 !");
}

void DiameterCalulatorTRI3S3::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorTRI3S3>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTRI3S3::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorTRI3S3>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTRI3S3::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorTRI3S3>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorTRI6S2::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==6)
    return DiameterCalulatorTRI3S2::ComputeForOneCellInternal(bg,bg+3,coordsPtr);
  else
    throw Exception("DiameterCalulatorTRI6S2::ComputeForOneCellInternal : input connectivity must be of size 6 !");
}

void DiameterCalulatorTRI6S2::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorTRI6S2>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTRI6S2::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorTRI6S2>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTRI6S2::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorTRI6S2>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorTRI6S3::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==6)
    return DiameterCalulatorTRI3S3::ComputeForOneCellInternal(bg,bg+3,coordsPtr);
  else
    throw Exception("DiameterCalulatorTRI6S3::ComputeForOneCellInternal : input connectivity must be of size 6 !");
}

void DiameterCalulatorTRI6S3::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorTRI6S3>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTRI6S3::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorTRI6S3>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTRI6S3::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorTRI6S3>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorTRI7S2::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==7)
    return DiameterCalulatorTRI3S2::ComputeForOneCellInternal(bg,bg+3,coordsPtr);
  else
    throw Exception("DiameterCalulatorTRI7S2::ComputeForOneCellInternal : input connectivity must be of size 7 !");
}

void DiameterCalulatorTRI7S2::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorTRI7S2>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTRI7S2::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorTRI7S2>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTRI7S2::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorTRI7S2>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorTRI7S3::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==7)
    return DiameterCalulatorTRI3S3::ComputeForOneCellInternal(bg,bg+3,coordsPtr);
  else
    throw Exception("DiameterCalulatorTRI7S3::ComputeForOneCellInternal : input connectivity must be of size 7 !");
}

void DiameterCalulatorTRI7S3::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorTRI7S3>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTRI7S3::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorTRI7S3>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTRI7S3::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorTRI7S3>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorQUAD4S2::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==4)
    {
      const double *a(coordsPtr+2*bg[0]),*b(coordsPtr+2*bg[1]),*c(coordsPtr+2*bg[2]),*d(coordsPtr+2*bg[3]);
      double l0[2],l1[2];
      DiffV2(a,c,l0);
      DiffV2(b,d,l1);
      return std::sqrt(std::max(SqNormV2(l0),SqNormV2(l1)));
    }
  else
    throw Exception("DiameterCalulatorQUAD4S2::ComputeForOneCellInternal : input connectivity must be of size 4 !");
}

void DiameterCalulatorQUAD4S2::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorQUAD4S2>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorQUAD4S2::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorQUAD4S2>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorQUAD4S2::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorQUAD4S2>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorQUAD4S3::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==4)
    {
      const double *a(coordsPtr+3*bg[0]),*b(coordsPtr+3*bg[1]),*c(coordsPtr+3*bg[2]),*d(coordsPtr+3*bg[3]);
      double l0[3],l1[3];
      DiffV3(a,c,l0);
      DiffV3(b,d,l1);
      return std::sqrt(std::max(SqNormV3(l0),SqNormV3(l1)));
    }
  else
    throw Exception("DiameterCalulatorQUAD4S3::ComputeForOneCellInternal : input connectivity must be of size 4 !");
}

void DiameterCalulatorQUAD4S3::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorQUAD4S3>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorQUAD4S3::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorQUAD4S3>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorQUAD4S3::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorQUAD4S3>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorQUAD8S2::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==8)
    return DiameterCalulatorQUAD4S2::ComputeForOneCellInternal(bg,bg+4,coordsPtr);
  else
    throw Exception("DiameterCalulatorQUAD8S2::ComputeForOneCellInternal : input connectivity must be of size 8 !");
}

void DiameterCalulatorQUAD8S2::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorQUAD8S2>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorQUAD8S2::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorQUAD8S2>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorQUAD8S2::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorQUAD8S2>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorQUAD8S3::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==8)
    return DiameterCalulatorQUAD4S3::ComputeForOneCellInternal(bg,bg+4,coordsPtr);
  else
    throw Exception("DiameterCalulatorQUAD8S3::ComputeForOneCellInternal : input connectivity must be of size 8 !");
}

void DiameterCalulatorQUAD8S3::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorQUAD8S3>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorQUAD8S3::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorQUAD8S3>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorQUAD8S3::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorQUAD8S3>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorQUAD9S2::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==9)
    return DiameterCalulatorQUAD4S2::ComputeForOneCellInternal(bg,bg+4,coordsPtr);
  else
    throw Exception("DiameterCalulatorQUAD9S2::ComputeForOneCellInternal : input connectivity must be of size 9 !");
}

void DiameterCalulatorQUAD9S2::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorQUAD9S2>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorQUAD9S2::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorQUAD9S2>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorQUAD9S2::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorQUAD9S2>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorQUAD9S3::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==9)
    return DiameterCalulatorQUAD4S3::ComputeForOneCellInternal(bg,bg+4,coordsPtr);
  else
    throw Exception("DiameterCalulatorQUAD8S3::ComputeForOneCellInternal : input connectivity must be of size 9 !");
}

void DiameterCalulatorQUAD9S3::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorQUAD9S3>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorQUAD9S3::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorQUAD9S3>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorQUAD9S3::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorQUAD9S3>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorTETRA4::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==4)
    {
      const double *a(coordsPtr+3*bg[0]),*b(coordsPtr+3*bg[1]),*c(coordsPtr+3*bg[2]),*d(coordsPtr+3*bg[3]);
      double l0[3],l1[3],l2[3],l3[3],l4[3],l5[3];
      DiffV3(a,b,l0);
      DiffV3(a,c,l1);
      DiffV3(b,c,l2);
      DiffV3(a,d,l3);
      DiffV3(b,d,l4);
      DiffV3(c,d,l5);
      double tmp[6];
      tmp[0]=SqNormV3(l0); tmp[1]=SqNormV3(l1); tmp[2]=SqNormV3(l2); tmp[3]=SqNormV3(l3); tmp[4]=SqNormV3(l4); tmp[5]=SqNormV3(l5);
      return std::sqrt(*std::max_element(tmp,tmp+6));
    }
  else
    throw Exception("DiameterCalulatorTETRA4::ComputeForOneCellInternal : input connectivity must be of size 4 !");
}

void DiameterCalulatorTETRA4::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorTETRA4>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTETRA4::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorTETRA4>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTETRA4::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorTETRA4>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorTETRA10::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==10)
    return DiameterCalulatorTETRA4::ComputeForOneCellInternal(bg,bg+4,coordsPtr);
  else
    throw Exception("DiameterCalulatorTETRA10::ComputeForOneCellInternal : input connectivity must be of size 10 !");
}

void DiameterCalulatorTETRA10::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorTETRA10>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTETRA10::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorTETRA10>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorTETRA10::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorTETRA10>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorHEXA8::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==8)
    {
      const double *p0(coordsPtr+3*bg[0]),*p1(coordsPtr+3*bg[1]),*p2(coordsPtr+3*bg[2]),*p3(coordsPtr+3*bg[3]),*p4(coordsPtr+3*bg[4]),*p5(coordsPtr+3*bg[5]),*p6(coordsPtr+3*bg[6]),*p7(coordsPtr+3*bg[7]);
      double l0[3],l1[3],l2[3],l3[3];
      DiffV3(p0,p6,l0);
      DiffV3(p1,p7,l1);
      DiffV3(p2,p4,l2);
      DiffV3(p3,p5,l3);
      double tmp[4];
      tmp[0]=SqNormV3(l0); tmp[1]=SqNormV3(l1); tmp[2]=SqNormV3(l2); tmp[3]=SqNormV3(l3);
      return std::sqrt(*std::max_element(tmp,tmp+4));
    }
  else
    throw Exception("DiameterCalulatorHEXA8::ComputeForOneCellInternal : input connectivity must be of size 8 !");
}

void DiameterCalulatorHEXA8::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorHEXA8>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorHEXA8::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorHEXA8>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorHEXA8::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorHEXA8>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorHEXA20::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==20)
    return DiameterCalulatorHEXA8::ComputeForOneCellInternal(bg,bg+8,coordsPtr);
  else
    throw Exception("DiameterCalulatorHEXA20::ComputeForOneCellInternal : input connectivity must be of size 20 !");
}

void DiameterCalulatorHEXA20::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorHEXA20>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorHEXA20::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorHEXA20>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorHEXA20::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorHEXA20>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorHEXA27::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==27)
    return DiameterCalulatorHEXA8::ComputeForOneCellInternal(bg,bg+8,coordsPtr);
  else
    throw Exception("DiameterCalulatorHEXA27::ComputeForOneCellInternal : input connectivity must be of size 27 !");
}

void DiameterCalulatorHEXA27::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorHEXA27>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorHEXA27::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorHEXA27>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorHEXA27::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorHEXA27>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorPENTA6::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==6)
    {
      const double *p0(coordsPtr+3*bg[0]),*p1(coordsPtr+3*bg[1]),*p2(coordsPtr+3*bg[2]),*p3(coordsPtr+3*bg[3]),*p4(coordsPtr+3*bg[4]),*p5(coordsPtr+3*bg[5]);
      double l0[3],l1[3],l2[3],l3[3],l4[3],l5[3];
      DiffV3(p0,p4,l0);
      DiffV3(p1,p3,l1);
      DiffV3(p1,p5,l2);
      DiffV3(p2,p4,l3);
      DiffV3(p0,p5,l4);
      DiffV3(p2,p3,l5);
      double tmp[6];
      tmp[0]=SqNormV3(l0); tmp[1]=SqNormV3(l1); tmp[2]=SqNormV3(l2); tmp[3]=SqNormV3(l3); tmp[4]=SqNormV3(l4); tmp[5]=SqNormV3(l5);
      return std::sqrt(*std::max_element(tmp,tmp+6));
    }
  else
    throw Exception("DiameterCalulatorPENTA6::ComputeForOneCellInternal : input connectivity must be of size 6 !");
}

void DiameterCalulatorPENTA6::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorPENTA6>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorPENTA6::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorPENTA6>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorPENTA6::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorPENTA6>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorPENTA15::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==15)
    return DiameterCalulatorPENTA6::ComputeForOneCellInternal(bg,bg+6,coordsPtr);
  else
    throw Exception("DiameterCalulatorPENTA15::ComputeForOneCellInternal : input connectivity must be of size 15 !");
}

void DiameterCalulatorPENTA15::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorPENTA15>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorPENTA15::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorPENTA15>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorPENTA15::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorPENTA15>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorPYRA5::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==5)
    {
      const double *p0(coordsPtr+3*bg[0]),*p1(coordsPtr+3*bg[1]),*p2(coordsPtr+3*bg[2]),*p3(coordsPtr+3*bg[3]),*p4(coordsPtr+3*bg[4]);
      double l0[3],l1[3],l2[3],l3[3],l4[3],l5[3];
      DiffV3(p0,p2,l0);
      DiffV3(p1,p3,l1);
      DiffV3(p0,p4,l2);
      DiffV3(p1,p4,l3);
      DiffV3(p2,p4,l4);
      DiffV3(p3,p4,l5);
      double tmp[6];
      tmp[0]=SqNormV3(l0); tmp[1]=SqNormV3(l1); tmp[2]=SqNormV3(l2); tmp[3]=SqNormV3(l3); tmp[4]=SqNormV3(l4); tmp[5]=SqNormV3(l5);
      return std::sqrt(*std::max_element(tmp,tmp+6));
    }
  else
    throw Exception("DiameterCalulatorPYRA5::ComputeForOneCellInternal : input connectivity must be of size 5 !");
}

void DiameterCalulatorPYRA5::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorPYRA5>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorPYRA5::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorPYRA5>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorPYRA5::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorPYRA5>(nbOfCells,connPtr,coordsPtr,resPtr);
}

//=================================================================

double DiameterCalulatorPYRA13::ComputeForOneCellInternal(const int *bg, const int *endd, const double *coordsPtr)
{
  if(std::distance(bg,endd)==13)
    return DiameterCalulatorPYRA5::ComputeForOneCellInternal(bg,bg+5,coordsPtr);
  else
    throw Exception("DiameterCalulatorPYRA13::ComputeForOneCellInternal : input connectivity must be of size 13 !");
}

void DiameterCalulatorPYRA13::computeForListOfCellIdsUMeshFrmt(const int *bgIds, const int *endIds, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForListOfCellIdsUMeshFrmt<DiameterCalulatorPYRA13>(bgIds,endIds,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorPYRA13::computeForRangeOfCellIdsUMeshFrmt(int bgId, int endId, const int *indPtr, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeForRangeOfCellIdsUMeshFrmt<DiameterCalulatorPYRA13>(bgId,endId,indPtr,connPtr,coordsPtr,resPtr);
}

void DiameterCalulatorPYRA13::computeFor1SGTUMeshFrmt(int nbOfCells, const int *connPtr, const double *coordsPtr, double *resPtr) const
{
  ComputeFor1SGTUMeshFrmt<DiameterCalulatorPYRA13>(nbOfCells,connPtr,coordsPtr,resPtr);
}
