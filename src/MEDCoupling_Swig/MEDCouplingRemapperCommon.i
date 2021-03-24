// Copyright (C) 2017-2020  CEA/DEN, EDF R&D
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
#include <memory>
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
             const std::vector<std::map<mcIdType,double> >& m=self->getCrudeMatrix();
             std::size_t sz=m.size();
             PyObject *ret=PyList_New(sz);
             for(std::size_t i=0;i<sz;i++)
               {
                 const std::map<mcIdType,double>& row=m[i];
                 PyObject *ret0=PyDict_New();
                 for(std::map<mcIdType,double>::const_iterator it=row.begin();it!=row.end();it++)
                   {
                     std::unique_ptr<PyObject,std::function<void(PyObject*)>> k(PyInt_FromLong((*it).first),[](PyObject *obj) { Py_XDECREF(obj); } ),v(PyFloat_FromDouble((*it).second),[](PyObject *obj) { Py_XDECREF(obj); } );
                     PyDict_SetItem(ret0,k.get(),v.get());
                   }
                 PyList_SetItem(ret,i,ret0);
               }
             return ret;
           }
#if defined(WITH_NUMPY) && defined(WITH_SCIPY)
           PyObject *getCrudeCSRMatrix() const
           {
             return ToCSRMatrix(self->getCrudeMatrix(),self->getNumberOfColsOfMatrix());
           }
           static PyObject *ToCSRMatrix(PyObject *m, mcIdType nbOfCols)
           {
              std::vector<std::map<mcIdType,double> > mCpp;
              convertToVectMapIntDouble(m,mCpp);
              return ToCSRMatrix(mCpp,nbOfCols);
           }
#endif
           void setCrudeMatrix(const MEDCouplingMesh *srcMesh, const MEDCouplingMesh *targetMesh, const std::string& method, PyObject *m)
           {
             std::vector<std::map<mcIdType,double> > mCpp;
             if(isCSRMatrix(m))
               {
#if defined(WITH_NUMPY) && defined(WITH_SCIPY)
                 PyObject *indptr(PyObject_GetAttrString(m,"indptr"));
                 PyObject *indices(PyObject_GetAttrString(m,"indices"));
                 PyObject *data(PyObject_GetAttrString(m,"data"));
                 MCAuto<DataArrayInt32> indptrPtr, indicesPtr;
                 // csr_matrix.indptr and csr_matrix.indices are always dtype==int32
// #if defined(MEDCOUPLING_USE_64BIT_IDS)
//                  indptrPtr = MEDCoupling_DataArrayInt64_New__SWIG_1(indptr,NULL,NULL);
//                  indicesPtr = MEDCoupling_DataArrayInt64_New__SWIG_1(indices,NULL,NULL);
// #else
                 indptrPtr = MEDCoupling_DataArrayInt32_New__SWIG_1(indptr,NULL,NULL);
                 indicesPtr = MEDCoupling_DataArrayInt32_New__SWIG_1(indices,NULL,NULL);
//#endif
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
             std::vector<std::map<mcIdType,double> > mCpp;
             convertToVectMapIntDouble(m,mCpp);
             self->setCrudeMatrixEx(src,target,mCpp);
           }
         }
    };
}

