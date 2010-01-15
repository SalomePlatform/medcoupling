//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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

#include <map>
#include <vector>

namespace ParaMEDMEM
{
  class MEDCouplingMesh;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
}

namespace ParaMEDMEM
{
  class MEDCOUPLING_EXPORT MEDCouplingRemapper : public TimeLabel, public INTERP_KERNEL::InterpolationOptions
  {
  public:
    MEDCouplingRemapper();
    ~MEDCouplingRemapper();
    int prepare(const MEDCouplingMesh *srcMesh, const MEDCouplingMesh *targetMesh, const char *method);
    void transfer(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField, double dftValue);
    void reverseTransfer(MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *targetField, double dftValue);
    MEDCouplingFieldDouble *transferField(const MEDCouplingFieldDouble *srcField, double dftValue);
    MEDCouplingFieldDouble *reverseTransferField(const MEDCouplingFieldDouble *targetField, double dftValue);
    bool setOptionInt(const std::string& key, int value);
    bool setOptionDouble(const std::string& key, double value);
    bool setOptionString(const std::string& key, std::string& value);
  private:
    int prepareUU(const char *method);
    int prepareEE(const char *method);
    void updateTime();
    void releaseData(bool matrixSuppression);
    void computeDeno(NatureOfField nat, const MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *trgField);
    void computeDenoFromScratch(NatureOfField nat, const MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *trgField);
    void computeProduct(const double *inputPointer, int inputNbOfCompo, double dftValue, double *resPointer);
    void computeReverseProduct(const double *inputPointer, int inputNbOfCompo, double dftValue, double *resPointer);
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
