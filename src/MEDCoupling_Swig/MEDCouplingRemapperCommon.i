// Copyright (C) 2017-2019  CEA/DEN, EDF R&D
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

#define MEDCOUPLING_EXPORT
#define INTERPKERNEL_EXPORT
#define MEDCOUPLINGREMAPPER_EXPORT

%newobject MEDCoupling::MEDCouplingRemapper::transferField;
%newobject MEDCoupling::MEDCouplingRemapper::reverseTransferField;

%{
#include "MEDCouplingRemapper.hxx"
%}

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
      int prepare(const MEDCouplingMesh *srcMesh, const MEDCouplingMesh *targetMesh, const std::string& method);
      int prepareEx(const MEDCouplingFieldTemplate *src, const MEDCouplingFieldTemplate *target);
      void transfer(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField, double dftValue);
      void partialTransfer(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField);
      void reverseTransfer(MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *targetField, double dftValue);
      MEDCouplingFieldDouble *transferField(const MEDCouplingFieldDouble *srcField, double dftValue);
      MEDCouplingFieldDouble *reverseTransferField(const MEDCouplingFieldDouble *targetField, double dftValue);
      bool setOptionInt(const std::string& key, int value);
      bool setOptionDouble(const std::string& key, double value);
      bool setOptionString(const std::string& key, const std::string& value);
      int getInterpolationMatrixPolicy() const;
      void setInterpolationMatrixPolicy(int newInterpMatPol);
      //
      int nullifiedTinyCoeffInCrudeMatrixAbs(double maxValAbs);
      int nullifiedTinyCoeffInCrudeMatrix(double scaleFactor);
      double getMaxValueInCrudeMatrix() const;
      int getNumberOfColsOfMatrix() const;
      static std::string BuildMethodFrom(const std::string& meth1, const std::string& meth2);
      %extend
         {
           PyObject *getCrudeMatrix() const
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
           PyObject *getCrudeCSRMatrix() const
           {
             return ToCSRMatrix(self->getCrudeMatrix(),self->getNumberOfColsOfMatrix());
           }
#endif
           void setCrudeMatrix(const MEDCouplingMesh *srcMesh, const MEDCouplingMesh *targetMesh, const std::string& method, PyObject *m)
           {
             std::vector<std::map<int,double> > mCpp;
             if(isCSRMatrix(m))
               {
#if defined(WITH_NUMPY) && defined(WITH_SCIPY)
                 PyObject *indptr(PyObject_GetAttrString(m,"indptr"));
                 PyObject *indices(PyObject_GetAttrString(m,"indices"));
                 PyObject *data(PyObject_GetAttrString(m,"data"));
                 MCAuto<DataArrayInt> indptrPtr(MEDCoupling_DataArrayInt_New__SWIG_1(indptr,NULL,NULL));
                 MCAuto<DataArrayInt> indicesPtr(MEDCoupling_DataArrayInt_New__SWIG_1(indices,NULL,NULL));
                 MCAuto<DataArrayDouble> dataPtr(MEDCoupling_DataArrayDouble_New__SWIG_1(data,NULL,NULL));
                 convertCSR_MCDataToVectMapIntDouble(indptrPtr,indicesPtr,dataPtr,mCpp);
                 Py_XDECREF(data); Py_XDECREF(indptr); Py_XDECREF(indices);
#else
                 throw INTERP_KERNEL::Exception("pywrap of MEDCouplingRemapper::setCrudeMatrix : unexpected situation regarding numpy/scipy !");
#endif
               }
             else
               convertToVectMapIntDouble(m,mCpp);
             self->setCrudeMatrix(srcMesh,targetMesh,method,mCpp);
           }

           void setCrudeMatrixEx(const MEDCouplingFieldTemplate *src, const MEDCouplingFieldTemplate *target, PyObject *m)
           {
             std::vector<std::map<int,double> > mCpp;
             convertToVectMapIntDouble(m,mCpp);
             self->setCrudeMatrixEx(src,target,mCpp);
           }
         }
    };
}

