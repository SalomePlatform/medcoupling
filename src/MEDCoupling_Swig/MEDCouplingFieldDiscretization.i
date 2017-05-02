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

%newobject MEDCoupling::MEDCouplingFieldDiscretization::New;
%newobject MEDCoupling::MEDCouplingFieldDiscretization::deepCopy;
%newobject MEDCoupling::MEDCouplingFieldDiscretization::clone;
%newobject MEDCoupling::MEDCouplingFieldDiscretization::clonePartRange;
%newobject MEDCoupling::MEDCouplingFieldDiscretization::getOffsetArr;
%newobject MEDCoupling::MEDCouplingFieldDiscretization::getLocalizationOfDiscValues;
%newobject MEDCoupling::MEDCouplingFieldDiscretization::getMeasureField;
%newobject MEDCoupling::MEDCouplingFieldDiscretization::clonePart;
%newobject MEDCoupling::MEDCouplingFieldDiscretization::getValueOnMulti;
%newobject MEDCoupling::MEDCouplingFieldDiscretization::computeTupleIdsToSelectFromCellIds;
%newobject MEDCoupling::MEDCouplingFieldDiscretizationKriging::PerformDriftOfVec;

namespace MEDCoupling
{
  class MEDCouplingFieldDiscretization : public RefCountObject, public TimeLabel
  {
  public:
    static MEDCouplingFieldDiscretization *New(TypeOfField type) throw(INTERP_KERNEL::Exception);
    double getPrecision() const throw(INTERP_KERNEL::Exception);
    void setPrecision(double val) throw(INTERP_KERNEL::Exception);
    static TypeOfField GetTypeOfFieldFromStringRepr(const std::string& repr) throw(INTERP_KERNEL::Exception);
    virtual TypeOfField getEnum() const throw(INTERP_KERNEL::Exception);
    virtual bool isEqual(const MEDCouplingFieldDiscretization *other, double eps) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingFieldDiscretization *other, double eps) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDiscretization *deepCopy() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDiscretization *clone() const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDiscretization *clonePartRange(int beginCellIds, int endCellIds, int stepCellIds) const throw(INTERP_KERNEL::Exception);
    virtual std::string getStringRepr() const throw(INTERP_KERNEL::Exception);
    virtual const char *getRepr() const throw(INTERP_KERNEL::Exception);
    virtual int getNumberOfTuples(const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception);
    virtual int getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getOffsetArr(const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayDouble *getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception);
    virtual void checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception);
    virtual double getIJK(const MEDCouplingMesh *mesh, const DataArrayDouble *da, int cellId, int nodeIdInCell, int compoId) const throw(INTERP_KERNEL::Exception);
    virtual void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArray *da) const throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingFieldDouble *getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const throw(INTERP_KERNEL::Exception);
    virtual void setGaussLocalizationOnType(const MEDCouplingMesh *m, INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                            const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception);
    virtual void clearGaussLocalizations() throw(INTERP_KERNEL::Exception);
    virtual MEDCouplingGaussLocalization& getGaussLocalization(int locId) throw(INTERP_KERNEL::Exception);
    virtual int getNbOfGaussLocalization() const throw(INTERP_KERNEL::Exception);
    virtual int getGaussLocalizationIdOfOneCell(int cellId) const throw(INTERP_KERNEL::Exception);
    virtual int getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception);
    %extend
    {
      virtual MEDCouplingFieldDiscretization *clonePart(PyObject *li)
      {
        int sz=0,sw=-1,val1=-1;
        std::vector<int> val2;
        const int *inp=convertIntStarLikePyObjToCppIntStar(li,sw,sz,val1,val2);
        return self->clonePart(inp,inp+sz);
      }
      
      virtual PyObject *buildSubMeshDataRange(const MEDCouplingMesh *mesh, int beginCellIds, int endCellIds, int stepCellIds, int& beginOut, int& endOut, int& stepOut, DataArrayInt *&di) const throw(INTERP_KERNEL::Exception)
      {
        DataArrayInt *ret1=0;
        int bb,ee,ss;
        MEDCouplingMesh *ret0=self->buildSubMeshDataRange(mesh,beginCellIds,endCellIds,stepCellIds,bb,ee,ss,ret1);
        PyObject *res=PyTuple_New(2);
        PyTuple_SetItem(res,0,convertMesh(ret0, SWIG_POINTER_OWN | 0 ));
        if(ret1)
          PyTuple_SetItem(res,1,SWIG_NewPointerObj((void*)ret1,SWIGTYPE_p_MEDCoupling__DataArrayInt,SWIG_POINTER_OWN | 0));
        else
          {
            PyObject *res1=PySlice_New(PyInt_FromLong(bb),PyInt_FromLong(ee),PyInt_FromLong(ss));
            PyTuple_SetItem(res,1,res1);
          }
        return res;
      }

      virtual int getNumberOfTuplesExpectedRegardingCode(PyObject *code, PyObject *idsPerType) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> inp0;
        convertPyToNewIntArr4(code,1,3,inp0);
        std::vector<const DataArrayInt *> inp1;
        convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayInt *>(idsPerType,SWIGTYPE_p_MEDCoupling__DataArrayInt,"DataArrayInt",inp1);
        return self->getNumberOfTuplesExpectedRegardingCode(inp0,inp1);
      }

      virtual PyObject *computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, PyObject *tupleIds) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> vVal; int iVal=-1;
        int sz=-1,sw=0;
        const int *tupleIdsBg=convertIntStarLikePyObjToCppIntStar(tupleIds,sw,sz,iVal,vVal);
        if(sw==0)
          throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretization::computeMeshRestrictionFromTupleIds : none parameter in input !");
        DataArrayInt *ret0=0,*ret1=0;
        self->computeMeshRestrictionFromTupleIds(mesh,tupleIdsBg,tupleIdsBg+sz,ret0,ret1);
        PyObject *pyRet=PyTuple_New(2);
        PyTuple_SetItem(pyRet,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(pyRet,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return pyRet;
      }

      virtual PyObject *normL1(const MEDCouplingMesh *mesh, const DataArrayDouble *arr) const throw(INTERP_KERNEL::Exception)
      {
        if(!arr)
          throw INTERP_KERNEL::Exception("wrap of MEDCouplingFieldDiscretization::normL1 : input array is null !");
        int sz(arr->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->normL1(mesh,arr,tmp);
        return convertDblArrToPyList<double>(tmp,sz);
      }

      virtual PyObject *normL2(const MEDCouplingMesh *mesh, const DataArrayDouble *arr) const throw(INTERP_KERNEL::Exception)
      {
        if(!arr)
          throw INTERP_KERNEL::Exception("wrap of MEDCouplingFieldDiscretization::normL2 : input array is null !");
        int sz(arr->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->normL2(mesh,arr,tmp);
        return convertDblArrToPyList<double>(tmp,sz);
      }

      virtual PyObject *integral(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, bool isWAbs) const throw(INTERP_KERNEL::Exception)
      {
        if(!arr)
          throw INTERP_KERNEL::Exception("wrap of MEDCouplingFieldDiscretization::integral : input array is null !");
        int sz(arr->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->integral(mesh,arr,isWAbs,tmp);
        return convertDblArrToPyList<double>(tmp,sz);
      }

      virtual PyObject *getCellIdsHavingGaussLocalization(int locId) const throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> tmp;
        self->getCellIdsHavingGaussLocalization(locId,tmp);
        DataArrayInt *ret=DataArrayInt::New();
        ret->alloc((int)tmp.size(),1);
        std::copy(tmp.begin(),tmp.end(),ret->getPointer());
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
      }

      virtual void setGaussLocalizationOnCells(const MEDCouplingMesh *m, PyObject *li, const std::vector<double>& refCoo,
                                               const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTYPE_p_MEDCoupling__DataArrayInt, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            int size;
            INTERP_KERNEL::AutoPtr<int> tmp=convertPyToNewIntArr2(li,&size);
            self->setGaussLocalizationOnCells(m,tmp,((int *)tmp)+size,refCoo,gsCoo,wg);
          }
        else
          {
            DataArrayInt *da2=reinterpret_cast< DataArrayInt * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
            da2->checkAllocated();
            self->setGaussLocalizationOnCells(m,da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems(),refCoo,gsCoo,wg);
          }
      }

      virtual PyObject *getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception)
      {
        std::set<int> ret=self->getGaussLocalizationIdsOfOneType(type);
        return convertIntArrToPyList3(ret);
      }

      virtual PyObject *getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, PyObject *sl) const throw(INTERP_KERNEL::Exception)
      {
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        int sw;
        if(!mesh)
          throw INTERP_KERNEL::Exception("Python wrap of MEDCouplingFieldDiscretization::getValueOn : no underlying mesh !");
        int spaceDim=mesh->getSpaceDimension();
        const char msg[]="Python wrap of MEDCouplingFieldDiscretization::getValueOn : ";
        const double *spaceLoc=convertObjToPossibleCpp5_Safe(sl,sw,val,a,aa,bb,msg,1,spaceDim,true);
        //
        INTERP_KERNEL::AutoPtr<double> res(new double[spaceDim]);
        self->getValueOn(arr,mesh,spaceLoc,res);
        return convertDblArrToPyList<double>(res,spaceDim);
      }

      virtual PyObject *getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, int i, int j, int k) const throw(INTERP_KERNEL::Exception)
      {
        if(!arr)
          throw INTERP_KERNEL::Exception("wrap of MEDCouplingFieldDiscretization::getValueOnPos : input array is null !");
        int sz(arr->getNumberOfComponents());
         INTERP_KERNEL::AutoPtr<double> res=new double[sz];
         self->getValueOnPos(arr,mesh,i,j,k,res);
         return convertDblArrToPyList<double>(res,sz);
       }
      
      virtual DataArrayDouble *getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, PyObject *loc) const throw(INTERP_KERNEL::Exception)
      {
        if(!mesh)
          throw INTERP_KERNEL::Exception("Python wrap MEDCouplingFieldDiscretization::getValueOnMulti : null input mesh !");
        //
        int sw,nbPts;
        double v0; MEDCoupling::DataArrayDouble *v1(0); MEDCoupling::DataArrayDoubleTuple *v2(0); std::vector<double> v3;
        const double *inp=convertObjToPossibleCpp5_Safe2(loc,sw,v0,v1,v2,v3,"wrap of MEDCouplingFieldDouble::getValueOnMulti",
                                                         mesh->getSpaceDimension(),true,nbPts);
        return self->getValueOnMulti(arr,mesh,inp,nbPts);
      }

      virtual void renumberCells(PyObject *li, bool check=true) throw(INTERP_KERNEL::Exception)
      {
        int sw,sz(-1);
        int v0; std::vector<int> v1;
        const int *ids(convertIntStarLikePyObjToCppIntStar(li,sw,sz,v0,v1));
        self->renumberCells(ids,check);
      }

      virtual void renumberArraysForCell(const MEDCouplingMesh *mesh, PyObject *arrays,
                                         PyObject *old2New, bool check) throw(INTERP_KERNEL::Exception)
      {
        std::vector<DataArray *> input1;
        convertFromPyObjVectorOfObj<MEDCoupling::DataArray *>(arrays,SWIGTYPE_p_MEDCoupling__DataArray,"DataArray",input1);
        //
        int sw,sz(-1);
        int v0; std::vector<int> v1;
        const int *old2NewBg(convertIntStarLikePyObjToCppIntStar(old2New,sw,sz,v0,v1));
        //
        self->renumberArraysForCell(mesh,input1,old2NewBg,check);
      }
      
      virtual DataArrayInt *computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, PyObject *cellIds) const throw(INTERP_KERNEL::Exception)
      {
        int sw,sz(-1);
        int v0; std::vector<int> v1;
        const int *cellIdsBg(convertIntStarLikePyObjToCppIntStar(cellIds,sw,sz,v0,v1));
        return self->computeTupleIdsToSelectFromCellIds(mesh,cellIdsBg,cellIdsBg+sz);
      }

      virtual PyObject *buildSubMeshData(const MEDCouplingMesh *mesh, PyObject *ids) const throw(INTERP_KERNEL::Exception)
      {
        int sw,sz(-1);
        int v0; std::vector<int> v1;
        const int *idsBg(convertIntStarLikePyObjToCppIntStar(ids,sw,sz,v0,v1));
        DataArrayInt *di(0);
        MEDCouplingMesh *ret0=self->buildSubMeshData(mesh,idsBg,idsBg+sz,di);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,convertMesh(ret0, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(di),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      virtual void renumberValuesOnNodes(double epsOnVals, PyObject *old2New, int newNbOfNodes, DataArrayDouble *arr) const throw(INTERP_KERNEL::Exception)
      {
        int sw,sz(-1);
        int v0; std::vector<int> v1;
        const int *old2NewBg(convertIntStarLikePyObjToCppIntStar(old2New,sw,sz,v0,v1));
        self->renumberValuesOnNodes(epsOnVals,old2NewBg,newNbOfNodes,arr);
      }

      virtual void renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, PyObject *old2New, int newSz, DataArrayDouble *arr) const throw(INTERP_KERNEL::Exception)
      {
        int sw,sz(-1);
        int v0; std::vector<int> v1;
        const int *old2NewBg(convertIntStarLikePyObjToCppIntStar(old2New,sw,sz,v0,v1));
        self->renumberValuesOnCells(epsOnVals,mesh,old2NewBg,newSz,arr);
      }

      virtual void renumberValuesOnCellsR(const MEDCouplingMesh *mesh, PyObject *new2old, int newSz, DataArrayDouble *arr) const throw(INTERP_KERNEL::Exception)
      {
        int sw,sz(-1);
        int v0; std::vector<int> v1;
        const int *new2oldBg(convertIntStarLikePyObjToCppIntStar(new2old,sw,sz,v0,v1));
        self->renumberValuesOnCellsR(mesh,new2oldBg,newSz,arr);
      }
    }
  };

  class MEDCouplingFieldDiscretizationP0 : public MEDCouplingFieldDiscretization
  {
  };

  class MEDCouplingFieldDiscretizationOnNodes : public MEDCouplingFieldDiscretization
  {
  };

  class MEDCouplingFieldDiscretizationP1 : public MEDCouplingFieldDiscretizationOnNodes
  {
  };
  
  class MEDCouplingFieldDiscretizationPerCell : public MEDCouplingFieldDiscretization
  {
  public:
    void setArrayOfDiscIds(const DataArrayInt *adids) throw(INTERP_KERNEL::Exception);
    void checkNoOrphanCells() const throw(INTERP_KERNEL::Exception);
    %extend
    {
      PyObject *getArrayOfDiscIds() const
      {
        DataArrayInt *ret=const_cast<DataArrayInt *>(self->getArrayOfDiscIds());
        if(ret)
          ret->incrRef();
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 );
      }

      PyObject *splitIntoSingleGaussDicrPerCellType() const throw(INTERP_KERNEL::Exception)
      {
        std::vector<int> ret1;
        std::vector<DataArrayInt *> ret0=self->splitIntoSingleGaussDicrPerCellType(ret1);
        std::size_t sz=ret0.size();
        PyObject *pyRet=PyTuple_New(2);
        PyObject *pyRet0=PyList_New((int)sz);
        PyObject *pyRet1=PyList_New((int)sz);
        for(std::size_t i=0;i<sz;i++)
          {
            PyList_SetItem(pyRet0,i,SWIG_NewPointerObj(SWIG_as_voidptr(ret0[i]),SWIGTYPE_p_MEDCoupling__DataArrayInt, SWIG_POINTER_OWN | 0 ));
            PyList_SetItem(pyRet1,i,PyInt_FromLong(ret1[i]));
          }
        PyTuple_SetItem(pyRet,0,pyRet0);
        PyTuple_SetItem(pyRet,1,pyRet1);
        return pyRet;
      }
    }
  protected:
    ~MEDCouplingFieldDiscretizationPerCell();
  };

  class MEDCouplingFieldDiscretizationGauss : public MEDCouplingFieldDiscretizationPerCell
  {
  public:
    MEDCouplingFieldDiscretizationGauss();
  };

  class MEDCouplingFieldDiscretizationGaussNE : public MEDCouplingFieldDiscretization
  {
  public:
    %extend
    {
      static PyObject *GetWeightArrayFromGeometricType(INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception)
      {
        std::size_t sz(0);
        const double *ret(MEDCouplingFieldDiscretizationGaussNE::GetWeightArrayFromGeometricType(geoType,sz));
        return convertDblArrToPyList<double>(ret,sz);
      }
      
      static PyObject *GetRefCoordsFromGeometricType(INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception)
      {
        std::size_t sz(0);
        const double *ret(MEDCouplingFieldDiscretizationGaussNE::GetRefCoordsFromGeometricType(geoType,sz));
        return convertDblArrToPyList<double>(ret,sz);
      }
      
      static PyObject *GetLocsFromGeometricType(INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception)
      {
        std::size_t sz(0);
        const double *ret(MEDCouplingFieldDiscretizationGaussNE::GetLocsFromGeometricType(geoType,sz));
        return convertDblArrToPyList<double>(ret,sz);
      }
    }
  };

  class MEDCouplingFieldDiscretizationKriging : public MEDCouplingFieldDiscretizationOnNodes
  {
  public:
    static DataArrayDouble *PerformDriftOfVec(const DataArrayDouble *arr, int isDrift) throw(INTERP_KERNEL::Exception);
    %extend
    {
      PyObject *computeVectorOfCoefficients(const MEDCouplingMesh *mesh, const DataArrayDouble *arr) const throw(INTERP_KERNEL::Exception)
      {
        int ret1;
        DataArrayDouble *ret0=self->computeVectorOfCoefficients(mesh,arr,ret1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,PyInt_FromLong(ret1));
        return ret;
      }
      
      PyObject *computeInverseMatrix(const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception)
      {
        int ret1(-1),ret2(-1);
        DataArrayDouble *ret0=self->computeInverseMatrix(mesh,ret1,ret2);
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,PyInt_FromLong(ret1));
        PyTuple_SetItem(ret,2,PyInt_FromLong(ret2));
        return ret;
      }

      PyObject *computeMatrix(const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception)
      {
        int ret1(-1),ret2(-1);
        DataArrayDouble *ret0=self->computeMatrix(mesh,ret1,ret2);
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,PyInt_FromLong(ret1));
        PyTuple_SetItem(ret,2,PyInt_FromLong(ret2));
        return ret;
      }
      
      PyObject *computeEvaluationMatrixOnGivenPts(const MEDCouplingMesh *mesh, PyObject *locs) const throw(INTERP_KERNEL::Exception)
      {
        if(!mesh)
          throw INTERP_KERNEL::Exception("wrap of MEDCouplingFieldDiscretizationKriging::computeEvaluationMatrixOnGivenPts : input mesh is empty !");
        int sw,nbPts;
        double v0; MEDCoupling::DataArrayDouble *v1(0); MEDCoupling::DataArrayDoubleTuple *v2(0); std::vector<double> v3;
        const double *inp=convertObjToPossibleCpp5_Safe2(locs,sw,v0,v1,v2,v3,"wrap of MEDCouplingFieldDiscretizationKriging::computeEvaluationMatrixOnGivenPts",
                                                         mesh->getSpaceDimension(),true,nbPts);
        //
        int ret1(-1);
        DataArrayDouble *ret0=self->computeEvaluationMatrixOnGivenPts(mesh,inp,nbPts,ret1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,PyInt_FromLong(ret1));
        return ret;
      }

      void operateOnDenseMatrix(int spaceDimension, DataArrayDouble *myMatrix) const throw(INTERP_KERNEL::Exception)
      {
        if(!myMatrix || !myMatrix->isAllocated() || myMatrix->getNumberOfComponents()!=1)
          throw INTERP_KERNEL::Exception("Wrap of MEDCouplingFieldDiscretizationKriging::operateOnDenseMatrix : invalid input matrix as DataArrayDouble ! Must be allocated with one component !");
        self->operateOnDenseMatrix(spaceDimension,myMatrix->getNumberOfTuples(),myMatrix->getPointer());
      }

      PyObject *performDrift(const DataArrayDouble *matr, const DataArrayDouble *arr) const throw(INTERP_KERNEL::Exception)
      {
        int ret1(-1);
        DataArrayDouble *ret0(self->performDrift(matr,arr,ret1));
        PyObject *res(PyTuple_New(2));
        PyTuple_SetItem(res,0,SWIG_NewPointerObj((void*)ret0,SWIGTYPE_p_MEDCoupling__DataArrayDouble,SWIG_POINTER_OWN | 0));
        PyTuple_SetItem(res,1,PyInt_FromLong(ret1));
        return res;
      }

      static PyObject *PerformDriftRect(const DataArrayDouble *matr, const DataArrayDouble *arr) throw(INTERP_KERNEL::Exception)
      {
        int ret1(-1);
        DataArrayDouble *ret0(MEDCouplingFieldDiscretizationKriging::PerformDriftRect(matr,arr,ret1));
        PyObject *res(PyTuple_New(2));
        PyTuple_SetItem(res,0,SWIG_NewPointerObj((void*)ret0,SWIGTYPE_p_MEDCoupling__DataArrayDouble,SWIG_POINTER_OWN | 0));
        PyTuple_SetItem(res,1,PyInt_FromLong(ret1));
        return res;
      }

      static void OperateOnDenseMatrixH3(DataArrayDouble *myMatrix) throw(INTERP_KERNEL::Exception)
      {
        if(!myMatrix || !myMatrix->isAllocated() || myMatrix->getNumberOfComponents()!=1)
          throw INTERP_KERNEL::Exception("Wrap of MEDCouplingFieldDiscretizationKriging::OperateOnDenseMatrixH3 : invalid input matrix as DataArrayDouble ! Must be allocated with one component !");
        MEDCouplingFieldDiscretizationKriging::OperateOnDenseMatrixH3(myMatrix->getNumberOfTuples(),myMatrix->getPointer());
      }
      
      static void OperateOnDenseMatrixH2Ln(DataArrayDouble *myMatrix) throw(INTERP_KERNEL::Exception)
      {
        if(!myMatrix || !myMatrix->isAllocated() || myMatrix->getNumberOfComponents()!=1)
          throw INTERP_KERNEL::Exception("Wrap of MEDCouplingFieldDiscretizationKriging::OperateOnDenseMatrixH2Ln : invalid input matrix as DataArrayDouble ! Must be allocated with one component !");
        MEDCouplingFieldDiscretizationKriging::OperateOnDenseMatrixH2Ln(myMatrix->getNumberOfTuples(),myMatrix->getPointer());
      }
    }
  };
}
