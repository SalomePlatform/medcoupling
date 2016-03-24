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

#ifndef __PARAMEDMEM_MEDCOUPLINGREMAPPER_HXX__
#define __PARAMEDMEM_MEDCOUPLINGREMAPPER_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "InterpolationOptions.hxx"
#include "MEDCouplingNatureOfField.hxx"
#include "MCAuto.hxx"

#include "InterpKernelException.hxx"

#include <map>
#include <vector>

namespace MEDCoupling
{
  class MEDCouplingMesh;
  class MEDCouplingFieldDouble;
  class MEDCouplingFieldTemplate;
}

namespace MEDCoupling
{
  typedef enum
  {
    IK_ONLY_PREFERED = 0,
    NOT_IK_ONLY_PREFERED = 1,
    IK_ONLY_FORCED = 2,
    NOT_IK_ONLY_FORCED =3
  } InterpolationMatrixPolicy;

  class MEDCouplingRemapper : public TimeLabel, public INTERP_KERNEL::InterpolationOptions
  {
  public:
    MEDCOUPLINGREMAPPER_EXPORT MEDCouplingRemapper();
    MEDCOUPLINGREMAPPER_EXPORT ~MEDCouplingRemapper();
    MEDCOUPLINGREMAPPER_EXPORT int prepare(const MEDCouplingMesh *srcMesh, const MEDCouplingMesh *targetMesh, const std::string& method);
    MEDCOUPLINGREMAPPER_EXPORT int prepareEx(const MEDCouplingFieldTemplate *src, const MEDCouplingFieldTemplate *target);
    MEDCOUPLINGREMAPPER_EXPORT void transfer(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField, double dftValue);
    MEDCOUPLINGREMAPPER_EXPORT void partialTransfer(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField);
    MEDCOUPLINGREMAPPER_EXPORT void reverseTransfer(MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *targetField, double dftValue);
    MEDCOUPLINGREMAPPER_EXPORT MEDCouplingFieldDouble *transferField(const MEDCouplingFieldDouble *srcField, double dftValue);
    MEDCOUPLINGREMAPPER_EXPORT MEDCouplingFieldDouble *reverseTransferField(const MEDCouplingFieldDouble *targetField, double dftValue);
    MEDCOUPLINGREMAPPER_EXPORT bool setOptionInt(const std::string& key, int value);
    MEDCOUPLINGREMAPPER_EXPORT bool setOptionDouble(const std::string& key, double value);
    MEDCOUPLINGREMAPPER_EXPORT bool setOptionString(const std::string& key, const std::string& value);
    MEDCOUPLINGREMAPPER_EXPORT int getInterpolationMatrixPolicy() const;
    MEDCOUPLINGREMAPPER_EXPORT void setInterpolationMatrixPolicy(int newInterpMatPol);
    //
    MEDCOUPLINGREMAPPER_EXPORT int nullifiedTinyCoeffInCrudeMatrixAbs(double maxValAbs);
    MEDCOUPLINGREMAPPER_EXPORT int nullifiedTinyCoeffInCrudeMatrix(double scaleFactor);
    MEDCOUPLINGREMAPPER_EXPORT double getMaxValueInCrudeMatrix() const;
  public:
    MEDCOUPLINGREMAPPER_EXPORT const std::vector<std::map<int,double> >& getCrudeMatrix() const;
    MEDCOUPLINGREMAPPER_EXPORT int getNumberOfColsOfMatrix() const;
    MEDCOUPLINGREMAPPER_EXPORT static void PrintMatrix(const std::vector<std::map<int,double> >& m);
    MEDCOUPLINGREMAPPER_EXPORT static std::string BuildMethodFrom(const std::string& meth1, const std::string& meth2);
  private:
    int prepareInterpKernelOnly();
    int prepareInterpKernelOnlyUU();
    int prepareInterpKernelOnlyEE();
    int prepareInterpKernelOnlyUC();
    int prepareInterpKernelOnlyCU();
    int prepareInterpKernelOnlyCC();
    //
    int prepareNotInterpKernelOnly();
    int prepareNotInterpKernelOnlyGaussGauss();
    //
    static int CheckInterpolationMethodManageableByNotOnlyInterpKernel(const std::string& method);
    //
    bool isInterpKernelOnlyOrNotOnly() const;
    void updateTime() const;
    void checkPrepare() const;
    std::string checkAndGiveInterpolationMethodStr(std::string& srcMeth, std::string& trgMeth) const;
    void releaseData(bool matrixSuppression);
    void transferUnderground(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField, bool isDftVal, double dftValue);
    void computeDeno(NatureOfField nat, const MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *trgField);
    void computeDenoFromScratch(NatureOfField nat, const MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *trgField);
    void computeProduct(const double *inputPointer, int inputNbOfCompo, bool isDftVal, double dftValue, double *resPointer);
    void computeReverseProduct(const double *inputPointer, int inputNbOfCompo, double dftValue, double *resPointer);
    void buildFinalInterpolationMatrixByConvolution(const std::vector< std::map<int,double> >& m1D,
                                                    const std::vector< std::map<int,double> >& m2D,
                                                    const int *corrCellIdSrc, int nbOf2DCellsSrc, int nbOf1DCellsSrc,
                                                    const int *corrCellIdTrg);
    static void ReverseMatrix(const std::vector<std::map<int,double> >& matIn, int nbColsMatIn,
                              std::vector<std::map<int,double> >& matOut);
    static void ComputeRowSumAndColSum(const std::vector<std::map<int,double> >& matrixDeno,
                                       std::vector<std::map<int,double> >& deno, std::vector<std::map<int,double> >& denoReverse);
    static void ComputeColSumAndRowSum(const std::vector<std::map<int,double> >& matrixDeno,
                                       std::vector<std::map<int,double> >& deno, std::vector<std::map<int,double> >& denoReverse);
  private:
    MCAuto<MEDCouplingFieldTemplate> _src_ft;
    MCAuto<MEDCouplingFieldTemplate> _target_ft;
    InterpolationMatrixPolicy _interp_matrix_pol;
    NatureOfField _nature_of_deno;
    unsigned int _time_deno_update;
    std::vector<std::map<int,double> > _matrix;
    std::vector<std::map<int,double> > _deno_multiply;
    std::vector<std::map<int,double> > _deno_reverse_multiply;
  };
}

#endif
