//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#ifndef __PARAMEDMEM_MEDCOUPLINGREMAPPER_HXX__
#define __PARAMEDMEM_MEDCOUPLINGREMAPPER_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "InterpolationOptions.hxx"
#include "MEDCouplingNatureOfField.hxx"
#include "InterpKernelException.hxx"

#include <map>
#include <vector>

namespace ParaMEDMEM
{
  class MEDCouplingMesh;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class MEDCouplingFieldTemplate;
}

namespace ParaMEDMEM
{
  class MEDCouplingRemapper : public TimeLabel, public INTERP_KERNEL::InterpolationOptions
  {
  public:
    MEDCOUPLINGREMAPPER_EXPORT MEDCouplingRemapper();
    MEDCOUPLINGREMAPPER_EXPORT ~MEDCouplingRemapper();
    MEDCOUPLINGREMAPPER_EXPORT int prepare(const MEDCouplingMesh *srcMesh, const MEDCouplingMesh *targetMesh, const char *method) throw(INTERP_KERNEL::Exception);
    MEDCOUPLINGREMAPPER_EXPORT int prepareEx(const MEDCouplingFieldTemplate *src, const MEDCouplingFieldTemplate *target) throw(INTERP_KERNEL::Exception);
    MEDCOUPLINGREMAPPER_EXPORT void transfer(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField, double dftValue) throw(INTERP_KERNEL::Exception);
    MEDCOUPLINGREMAPPER_EXPORT void reverseTransfer(MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *targetField, double dftValue) throw(INTERP_KERNEL::Exception);
    MEDCOUPLINGREMAPPER_EXPORT MEDCouplingFieldDouble *transferField(const MEDCouplingFieldDouble *srcField, double dftValue) throw(INTERP_KERNEL::Exception);
    MEDCOUPLINGREMAPPER_EXPORT MEDCouplingFieldDouble *reverseTransferField(const MEDCouplingFieldDouble *targetField, double dftValue) throw(INTERP_KERNEL::Exception);
    MEDCOUPLINGREMAPPER_EXPORT bool setOptionInt(const std::string& key, int value);
    MEDCOUPLINGREMAPPER_EXPORT bool setOptionDouble(const std::string& key, double value);
    MEDCOUPLINGREMAPPER_EXPORT bool setOptionString(const std::string& key, const std::string& value);
  public:
    MEDCOUPLINGREMAPPER_EXPORT static void printMatrix(const std::vector<std::map<int,double> >& m);
  private:
    int prepareUU(const char *method) throw(INTERP_KERNEL::Exception);
    int prepareEE(const char *method) throw(INTERP_KERNEL::Exception);
    void updateTime() const;
    void releaseData(bool matrixSuppression);
    void computeDeno(NatureOfField nat, const MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *trgField);
    void computeDenoFromScratch(NatureOfField nat, const MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *trgField) throw(INTERP_KERNEL::Exception);
    void computeProduct(const double *inputPointer, int inputNbOfCompo, double dftValue, double *resPointer);
    void computeReverseProduct(const double *inputPointer, int inputNbOfCompo, double dftValue, double *resPointer);
    void buildFinalInterpolationMatrixByConvolution(const std::vector< std::map<int,double> >& m1D,
                                                    const std::vector< std::map<int,double> >& m2D,
                                                    const int *corrCellIdSrc, int nbOf2DCellsSrc, int nbOf1DCellsSrc,
                                                    const int *corrCellIdTrg);
    static void reverseMatrix(const std::vector<std::map<int,double> >& matIn, int nbColsMatIn,
                              std::vector<std::map<int,double> >& matOut);
    static void computeRowSumAndColSum(const std::vector<std::map<int,double> >& matrixDeno,
                                       std::vector<std::map<int,double> >& deno, std::vector<std::map<int,double> >& denoReverse);
    static void computeColSumAndRowSum(const std::vector<std::map<int,double> >& matrixDeno,
                                       std::vector<std::map<int,double> >& deno, std::vector<std::map<int,double> >& denoReverse);
  private:
    MEDCouplingMesh *_src_mesh;
    MEDCouplingMesh *_target_mesh;
    std::string _src_method;
    std::string _target_method;
    NatureOfField _nature_of_deno;
    unsigned int _time_deno_update;
    std::vector<std::map<int,double> > _matrix;
    std::vector<std::map<int,double> > _deno_multiply;
    std::vector<std::map<int,double> > _deno_reverse_multiply;
  };
}

#endif
