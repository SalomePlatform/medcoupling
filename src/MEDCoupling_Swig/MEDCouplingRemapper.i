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

%module MEDCouplingRemapper

#define MEDCOUPLING_EXPORT
#define INTERPKERNEL_EXPORT
#define MEDCOUPLINGREMAPPER_EXPORT

%{
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingField.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingRemapper.hxx"

using namespace MEDCoupling;
using namespace INTERP_KERNEL;
%}

%newobject MEDCoupling::MEDCouplingRemapper::transferField;
%newobject MEDCoupling::MEDCouplingRemapper::reverseTransferField;

%include "MEDCouplingCommon.i"
%include "InterpolationOptions.hxx"

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
    private:
      void updateTime() const;
    public:
      MEDCouplingRemapper();
      ~MEDCouplingRemapper();
      int prepare(const MEDCouplingMesh *srcMesh, const MEDCouplingMesh *targetMesh, const std::string& method) throw(INTERP_KERNEL::Exception);
      int prepareEx(const MEDCouplingFieldTemplate *src, const MEDCouplingFieldTemplate *target) throw(INTERP_KERNEL::Exception);
      void transfer(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField, double dftValue) throw(INTERP_KERNEL::Exception);
      void partialTransfer(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField) throw(INTERP_KERNEL::Exception);
      void reverseTransfer(MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *targetField, double dftValue) throw(INTERP_KERNEL::Exception);
      MEDCouplingFieldDouble *transferField(const MEDCouplingFieldDouble *srcField, double dftValue) throw(INTERP_KERNEL::Exception);
      MEDCouplingFieldDouble *reverseTransferField(const MEDCouplingFieldDouble *targetField, double dftValue) throw(INTERP_KERNEL::Exception);
      bool setOptionInt(const std::string& key, int value) throw(INTERP_KERNEL::Exception);
      bool setOptionDouble(const std::string& key, double value) throw(INTERP_KERNEL::Exception);
      bool setOptionString(const std::string& key, const std::string& value) throw(INTERP_KERNEL::Exception);
      int getInterpolationMatrixPolicy() const throw(INTERP_KERNEL::Exception);
      void setInterpolationMatrixPolicy(int newInterpMatPol) throw(INTERP_KERNEL::Exception);
      //
      int nullifiedTinyCoeffInCrudeMatrixAbs(double maxValAbs) throw(INTERP_KERNEL::Exception);
      int nullifiedTinyCoeffInCrudeMatrix(double scaleFactor) throw(INTERP_KERNEL::Exception);
      double getMaxValueInCrudeMatrix() const throw(INTERP_KERNEL::Exception);
      int getNumberOfColsOfMatrix() const throw(INTERP_KERNEL::Exception);
      static std::string BuildMethodFrom(const std::string& meth1, const std::string& meth2) throw(INTERP_KERNEL::Exception);
      %extend
         {
           PyObject *getCrudeMatrix() const throw(INTERP_KERNEL::Exception)
           {
             const std::vector<std::map<int,double> >& m=self->getCrudeMatrix();
             std::size_t sz=m.size();
             PyObject *ret=PyList_New(sz);
             for(std::size_t i=0;i<sz;i++)
               {
                 const std::map<int,double>& row=m[i];
                 PyObject *ret0=PyDict_New();
                 for(std::map<int,double>::const_iterator it=row.begin();it!=row.end();it++)
                   PyDict_SetItem(ret0,PyInt_FromLong((*it).first),PyFloat_FromDouble((*it).second));
                 PyList_SetItem(ret,i,ret0);
               }
             return ret;
           }
#if defined(WITH_NUMPY) && defined(WITH_SCIPY)
           PyObject *getCrudeCSRMatrix() const throw(INTERP_KERNEL::Exception)
           {
             return ToCSRMatrix(self->getCrudeMatrix(),self->getNumberOfColsOfMatrix());
           }
#endif
         }
    };
}

%pythoncode %{
def MEDCouplingDataArrayDoublenew(cls,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayDouble____new___(cls,args)
def MEDCouplingDataArrayDoubleIadd(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayDouble____iadd___(self, self, *args)
def MEDCouplingDataArrayDoubleIsub(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayDouble____isub___(self, self, *args)
def MEDCouplingDataArrayDoubleImul(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayDouble____imul___(self, self, *args)
def MEDCouplingDataArrayDoubleIdiv(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayDouble____idiv___(self, self, *args)
def MEDCouplingDataArrayDoubleIpow(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayDouble____ipow___(self, self, *args)
def MEDCouplingFieldDoublenew(cls,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.MEDCouplingFieldDouble____new___(cls,args)
def MEDCouplingFieldDoubleIadd(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.MEDCouplingFieldDouble____iadd___(self, self, *args)
def MEDCouplingFieldDoubleIsub(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.MEDCouplingFieldDouble____isub___(self, self, *args)
def MEDCouplingFieldDoubleImul(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.MEDCouplingFieldDouble____imul___(self, self, *args)
def MEDCouplingFieldDoubleIdiv(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.MEDCouplingFieldDouble____idiv___(self, self, *args)
def MEDCouplingFieldDoubleIpow(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.MEDCouplingFieldDouble____ipow___(self, self, *args)
def MEDCouplingFieldIntnew(cls,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.MEDCouplingFieldInt____new___(cls,args)
def MEDCouplingFieldFloatnew(cls,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.MEDCouplingFieldFloat____new___(cls,args)
def MEDCouplingDataArrayBytenew(cls,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayByte____new___(cls,args)
def MEDCouplingDataArrayFloatnew(cls,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayFloat____new___(cls,args)
def MEDCouplingDataArrayFloatIadd(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayFloat____iadd___(self, self, *args)
def MEDCouplingDataArrayFloatIsub(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayFloat____isub___(self, self, *args)
def MEDCouplingDataArrayFloatImul(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayFloat____imul___(self, self, *args)
def MEDCouplingDataArrayFloatIdiv(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayFloat____idiv___(self, self, *args)
def MEDCouplingDataArrayIntnew(cls,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayInt____new___(cls,args)
def MEDCouplingDataArrayIntIadd(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayInt____iadd___(self, self, *args)
def MEDCouplingDataArrayIntIsub(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayInt____isub___(self, self, *args)
def MEDCouplingDataArrayIntImul(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayInt____imul___(self, self, *args)
def MEDCouplingDataArrayIntIdiv(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayInt____idiv___(self, self, *args)
def MEDCouplingDataArrayIntImod(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayInt____imod___(self, self, *args)
def MEDCouplingDataArrayIntIpow(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayInt____ipow___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIadd(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayDoubleTuple____iadd___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIsub(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayDoubleTuple____isub___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleImul(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayDoubleTuple____imul___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIdiv(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayDoubleTuple____idiv___(self, self, *args)
def MEDCouplingDataArrayIntTupleIadd(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayIntTuple____iadd___(self, self, *args)
def MEDCouplingDataArrayIntTupleIsub(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayIntTuple____isub___(self, self, *args)
def MEDCouplingDataArrayIntTupleImul(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayIntTuple____imul___(self, self, *args)
def MEDCouplingDataArrayIntTupleIdiv(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayIntTuple____idiv___(self, self, *args)
def MEDCouplingDataArrayIntTupleImod(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DataArrayIntTuple____imod___(self, self, *args)
def ParaMEDMEMDenseMatrixIadd(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DenseMatrix____iadd___(self, self, *args)
def ParaMEDMEMDenseMatrixIsub(self,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.DenseMatrix____isub___(self, self, *args)
def MEDCouplingUMeshnew(cls,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.MEDCouplingUMesh____new___(cls,args)
def MEDCoupling1DGTUMeshnew(cls,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.MEDCoupling1DGTUMesh____new___(cls,args)
def MEDCoupling1SGTUMeshnew(cls,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.MEDCoupling1SGTUMesh____new___(cls,args)
def MEDCouplingCurveLinearMeshnew(cls,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.MEDCouplingCurveLinearMesh____new___(cls,args)
def MEDCouplingCMeshnew(cls,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.MEDCouplingCMesh____new___(cls,args)
def MEDCouplingIMeshnew(cls,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.MEDCouplingIMesh____new___(cls,args)
def MEDCouplingExtrudedMeshnew(cls,*args):
    import _MEDCouplingRemapper
    return _MEDCouplingRemapper.MEDCouplingMappedExtrudedMesh____new___(cls,args)
%}

%include "MEDCouplingFinalize.i"
