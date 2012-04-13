// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

using namespace ParaMEDMEM;
using namespace INTERP_KERNEL;
%}

%newobject ParaMEDMEM::MEDCouplingRemapper::transferField;
%newobject ParaMEDMEM::MEDCouplingRemapper::reverseTransferField;

%include "MEDCoupling.i"
%include "InterpolationOptions.hxx"

namespace ParaMEDMEM
{
  class MEDCouplingRemapper : public TimeLabel, public INTERP_KERNEL::InterpolationOptions
    {
    private:
      void updateTime() const;
    public:
      MEDCouplingRemapper();
      ~MEDCouplingRemapper();
      int prepare(const MEDCouplingMesh *srcMesh, const MEDCouplingMesh *targetMesh, const char *method) throw(INTERP_KERNEL::Exception);
      int prepareEx(const MEDCouplingFieldTemplate *src, const MEDCouplingFieldTemplate *target) throw(INTERP_KERNEL::Exception);
      void transfer(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField, double dftValue) throw(INTERP_KERNEL::Exception);
      void partialTransfer(const MEDCouplingFieldDouble *srcField, MEDCouplingFieldDouble *targetField) throw(INTERP_KERNEL::Exception);
      void reverseTransfer(MEDCouplingFieldDouble *srcField, const MEDCouplingFieldDouble *targetField, double dftValue) throw(INTERP_KERNEL::Exception);
      MEDCouplingFieldDouble *transferField(const MEDCouplingFieldDouble *srcField, double dftValue) throw(INTERP_KERNEL::Exception);
      MEDCouplingFieldDouble *reverseTransferField(const MEDCouplingFieldDouble *targetField, double dftValue) throw(INTERP_KERNEL::Exception);
      bool setOptionInt(const std::string& key, int value);
      bool setOptionDouble(const std::string& key, double value);
      bool setOptionString(const std::string& key, const std::string& value);
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
         }
    };
}
