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
    static MEDCouplingFieldDiscretization *New(TypeOfField type);
    double getPrecision() const;
    void setPrecision(double val);
    static TypeOfField GetTypeOfFieldFromStringRepr(const std::string& repr);
    virtual TypeOfField getEnum() const;
    virtual bool isEqual(const MEDCouplingFieldDiscretization *other, double eps) const;
    virtual bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const;
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingFieldDiscretization *other, double eps) const;
    virtual MEDCouplingFieldDiscretization *deepCopy() const;
    virtual MEDCouplingFieldDiscretization *clone() const;
    virtual MEDCouplingFieldDiscretization *clonePartRange(int beginCellIds, int endCellIds, int stepCellIds) const;
    virtual std::string getStringRepr() const;
    virtual const char *getRepr() const;
    virtual int getNumberOfTuples(const MEDCouplingMesh *mesh) const;
    virtual int getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const;
    virtual DataArrayIdType *getOffsetArr(const MEDCouplingMesh *mesh) const;
    virtual DataArrayDouble *getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const;
    virtual void checkCompatibilityWithNature(NatureOfField nat) const;
    virtual double getIJK(const MEDCouplingMesh *mesh, const DataArrayDouble *da, int cellId, int nodeIdInCell, int compoId) const;
    virtual void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArray *da) const;
    virtual MEDCouplingFieldDouble *getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const;
    virtual void setGaussLocalizationOnType(const MEDCouplingMesh *m, INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                            const std::vector<double>& gsCoo, const std::vector<double>& wg);
    virtual void clearGaussLocalizations();
    virtual MEDCouplingGaussLocalization& getGaussLocalization(int locId);
    virtual int getNbOfGaussLocalization() const;
    virtual int getGaussLocalizationIdOfOneCell(int cellId) const;
    virtual int getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const;
    %extend
    {
      virtual MEDCouplingFieldDiscretization *clonePart(PyObject *li)
      {
        mcIdType sz=0,sw=-1,val1=-1;
        std::vector<mcIdType> val2;
        const mcIdType *inp=convertIntStarLikePyObjToCppIntStar(li,sw,sz,val1,val2);
        return self->clonePart(inp,inp+sz);
      }
      
      virtual PyObject *buildSubMeshDataRange(const MEDCouplingMesh *mesh, mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds, mcIdType& beginOut, mcIdType& endOut, mcIdType& stepOut, DataArrayIdType *&di) const
      {
        DataArrayIdType *ret1=0;
        mcIdType bb,ee,ss;
        MEDCouplingMesh *ret0=self->buildSubMeshDataRange(mesh,beginCellIds,endCellIds,stepCellIds,bb,ee,ss,ret1);
        PyObject *res=PyTuple_New(2);
        PyTuple_SetItem(res,0,convertMesh(ret0, SWIG_POINTER_OWN | 0 ));
        if(ret1)
          PyTuple_SetItem(res,1,SWIG_NewPointerObj((void*)ret1,SWIGTITraits<mcIdType>::TI,SWIG_POINTER_OWN | 0));
        else
          {
            PyObject *res1=PySlice_New(PyInt_FromLong(bb),PyInt_FromLong(ee),PyInt_FromLong(ss));
            PyTuple_SetItem(res,1,res1);
          }
        return res;
      }

      virtual mcIdType getNumberOfTuplesExpectedRegardingCode(PyObject *code, PyObject *idsPerType) const
      {
        std::vector<mcIdType> inp0;
        convertPyToNewIntArr4(code,1,3,inp0);
        std::vector<const DataArrayIdType *> inp1;
        convertFromPyObjVectorOfObj<const MEDCoupling::DataArrayIdType *>(idsPerType,SWIGTITraits<mcIdType>::TI,"DataArrayIdType",inp1);
        return self->getNumberOfTuplesExpectedRegardingCode(inp0,inp1);
      }

      virtual PyObject *computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, PyObject *tupleIds) const
      {
        std::vector<mcIdType> vVal; mcIdType iVal=-1;
        mcIdType sz=-1,sw=0;
        const mcIdType *tupleIdsBg=convertIntStarLikePyObjToCppIntStar(tupleIds,sw,sz,iVal,vVal);
        if(sw==0)
          throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretization::computeMeshRestrictionFromTupleIds : none parameter in input !");
        DataArrayIdType *ret0=0,*ret1=0;
        self->computeMeshRestrictionFromTupleIds(mesh,tupleIdsBg,tupleIdsBg+sz,ret0,ret1);
        PyObject *pyRet=PyTuple_New(2);
        PyTuple_SetItem(pyRet,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(pyRet,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return pyRet;
      }

      virtual PyObject *normL1(const MEDCouplingMesh *mesh, const DataArrayDouble *arr) const
      {
        if(!arr)
          throw INTERP_KERNEL::Exception("wrap of MEDCouplingFieldDiscretization::normL1 : input array is null !");
        std::size_t sz(arr->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->normL1(mesh,arr,tmp);
        return convertDblArrToPyList<double>(tmp,sz);
      }

      virtual PyObject *normL2(const MEDCouplingMesh *mesh, const DataArrayDouble *arr) const
      {
        if(!arr)
          throw INTERP_KERNEL::Exception("wrap of MEDCouplingFieldDiscretization::normL2 : input array is null !");
        std::size_t sz(arr->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->normL2(mesh,arr,tmp);
        return convertDblArrToPyList<double>(tmp,sz);
      }

      virtual PyObject *integral(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, bool isWAbs) const
      {
        if(!arr)
          throw INTERP_KERNEL::Exception("wrap of MEDCouplingFieldDiscretization::integral : input array is null !");
        std::size_t sz(arr->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<double> tmp=new double[sz];
        self->integral(mesh,arr,isWAbs,tmp);
        return convertDblArrToPyList<double>(tmp,sz);
      }

      virtual PyObject *getCellIdsHavingGaussLocalization(mcIdType locId) const
      {
        std::vector<mcIdType> tmp;
        self->getCellIdsHavingGaussLocalization(locId,tmp);
        DataArrayIdType *ret=DataArrayIdType::New();
        ret->alloc(tmp.size(),1);
        std::copy(tmp.begin(),tmp.end(),ret->getPointer());
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 );
      }

      virtual void setGaussLocalizationOnCells(const MEDCouplingMesh *m, PyObject *li, const std::vector<double>& refCoo,
                                               const std::vector<double>& gsCoo, const std::vector<double>& wg)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTITraits<mcIdType>::TI, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            mcIdType size;
            INTERP_KERNEL::AutoPtr<mcIdType> tmp=convertPyToNewIntArr2(li,&size);
            self->setGaussLocalizationOnCells(m,tmp,((mcIdType *)tmp)+size,refCoo,gsCoo,wg);
          }
        else
          {
            DataArrayIdType *da2=reinterpret_cast< DataArrayIdType * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayIdType instance expected !");
            da2->checkAllocated();
            self->setGaussLocalizationOnCells(m,da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems(),refCoo,gsCoo,wg);
          }
      }

      virtual PyObject *getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const
      {
        std::set<mcIdType> ret=self->getGaussLocalizationIdsOfOneType(type);
        return convertIntArrToPyList3(ret);
      }

      virtual PyObject *getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, PyObject *sl) const
      {
        double val;
        DataArrayDouble *a;
        DataArrayDoubleTuple *aa;
        std::vector<double> bb;
        mcIdType sw;
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

      virtual PyObject *getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, mcIdType i, mcIdType j, mcIdType k) const
      {
        if(!arr)
          throw INTERP_KERNEL::Exception("wrap of MEDCouplingFieldDiscretization::getValueOnPos : input array is null !");
        std::size_t sz(arr->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<double> res=new double[sz];
        self->getValueOnPos(arr,mesh,i,j,k,res);
        return convertDblArrToPyList<double>(res,sz);
      }
      
      virtual DataArrayDouble *getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, PyObject *loc) const
      {
        if(!mesh)
          throw INTERP_KERNEL::Exception("Python wrap MEDCouplingFieldDiscretization::getValueOnMulti : null input mesh !");
        //
        mcIdType sw,nbPts;
        double v0; MEDCoupling::DataArrayDouble *v1(0); MEDCoupling::DataArrayDoubleTuple *v2(0); std::vector<double> v3;
        const double *inp=convertObjToPossibleCpp5_Safe2(loc,sw,v0,v1,v2,v3,"wrap of MEDCouplingFieldDouble::getValueOnMulti",
                                                         mesh->getSpaceDimension(),true,nbPts);
        return self->getValueOnMulti(arr,mesh,inp,nbPts);
      }

      virtual void renumberCells(PyObject *li, bool check=true)
      {
        mcIdType sw,sz(-1);
        mcIdType v0; std::vector<mcIdType> v1;
        const mcIdType *ids(convertIntStarLikePyObjToCppIntStar(li,sw,sz,v0,v1));
        self->renumberCells(ids,check);
      }

      virtual void renumberArraysForCell(const MEDCouplingMesh *mesh, PyObject *arrays,
                                         PyObject *old2New, bool check)
      {
        std::vector<DataArray *> input1;
        convertFromPyObjVectorOfObj<MEDCoupling::DataArray *>(arrays,SWIGTYPE_p_MEDCoupling__DataArray,"DataArray",input1);
        //
        mcIdType sw,sz(-1);
        mcIdType v0; std::vector<mcIdType> v1;
        const mcIdType *old2NewBg(convertIntStarLikePyObjToCppIntStar(old2New,sw,sz,v0,v1));
        //
        self->renumberArraysForCell(mesh,input1,old2NewBg,check);
      }
      
      virtual DataArrayIdType *computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, PyObject *cellIds) const
      {
        mcIdType sw,sz(-1);
        mcIdType v0; std::vector<mcIdType> v1;
        const mcIdType *cellIdsBg(convertIntStarLikePyObjToCppIntStar(cellIds,sw,sz,v0,v1));
        return self->computeTupleIdsToSelectFromCellIds(mesh,cellIdsBg,cellIdsBg+sz);
      }

      virtual PyObject *buildSubMeshData(const MEDCouplingMesh *mesh, PyObject *ids)
      {
        mcIdType sw,sz(-1);
        mcIdType v0; std::vector<mcIdType> v1;
        const mcIdType *idsBg(convertIntStarLikePyObjToCppIntStar(ids,sw,sz,v0,v1));
        DataArrayIdType *di(0);
        MEDCouplingMesh *ret0=self->buildSubMeshData(mesh,idsBg,idsBg+sz,di);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,convertMesh(ret0, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(di),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      virtual void renumberValuesOnNodes(double epsOnVals, PyObject *old2New, mcIdType newNbOfNodes, DataArrayDouble *arr) const
      {
        mcIdType sw,sz(-1);
        mcIdType v0; std::vector<mcIdType> v1;
        const mcIdType *old2NewBg(convertIntStarLikePyObjToCppIntStar(old2New,sw,sz,v0,v1));
        self->renumberValuesOnNodes(epsOnVals,old2NewBg,newNbOfNodes,arr);
      }

      virtual void renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, PyObject *old2New, mcIdType newSz, DataArrayDouble *arr) const
      {
        mcIdType sw,sz(-1);
        mcIdType v0; std::vector<mcIdType> v1;
        const mcIdType *old2NewBg(convertIntStarLikePyObjToCppIntStar(old2New,sw,sz,v0,v1));
        self->renumberValuesOnCells(epsOnVals,mesh,old2NewBg,newSz,arr);
      }

      virtual void renumberValuesOnCellsR(const MEDCouplingMesh *mesh, PyObject *new2old, mcIdType newSz, DataArrayDouble *arr) const
      {
        mcIdType sw,sz(-1);
        mcIdType v0; std::vector<mcIdType> v1;
        const mcIdType *new2oldBg(convertIntStarLikePyObjToCppIntStar(new2old,sw,sz,v0,v1));
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
    void setArrayOfDiscIds(const DataArrayIdType *adids);
    void checkNoOrphanCells() const;
    %extend
    {
      PyObject *getArrayOfDiscIds() const
      {
        DataArrayIdType *ret=const_cast<DataArrayIdType *>(self->getArrayOfDiscIds());
        if(ret)
          ret->incrRef();
        return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 );
      }

      PyObject *splitIntoSingleGaussDicrPerCellType() const
      {
        std::vector<mcIdType> ret1;
        std::vector<DataArrayIdType *> ret0=self->splitIntoSingleGaussDicrPerCellType(ret1);
        std::size_t sz=ret0.size();
        PyObject *pyRet=PyTuple_New(2);
        PyObject *pyRet0=PyList_New((mcIdType)sz);
        PyObject *pyRet1=PyList_New((mcIdType)sz);
        for(std::size_t i=0;i<sz;i++)
          {
            PyList_SetItem(pyRet0,i,SWIG_NewPointerObj(SWIG_as_voidptr(ret0[i]),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
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
      static PyObject *GetWeightArrayFromGeometricType(INTERP_KERNEL::NormalizedCellType geoType)
      {
        std::size_t sz(0);
        const double *ret(MEDCouplingFieldDiscretizationGaussNE::GetWeightArrayFromGeometricType(geoType,sz));
        return convertDblArrToPyList<double>(ret,ToIdType(sz));
      }
      
      static PyObject *GetRefCoordsFromGeometricType(INTERP_KERNEL::NormalizedCellType geoType)
      {
        std::size_t sz(0);
        const double *ret(MEDCouplingFieldDiscretizationGaussNE::GetRefCoordsFromGeometricType(geoType,sz));
        return convertDblArrToPyList<double>(ret,ToIdType(sz));
      }
      
      static PyObject *GetLocsFromGeometricType(INTERP_KERNEL::NormalizedCellType geoType)
      {
        std::size_t sz(0);
        const double *ret(MEDCouplingFieldDiscretizationGaussNE::GetLocsFromGeometricType(geoType,sz));
        return convertDblArrToPyList<double>(ret,ToIdType(sz));
      }
    }
  };

  class MEDCouplingFieldDiscretizationKriging : public MEDCouplingFieldDiscretizationOnNodes
  {
  public:
    static DataArrayDouble *PerformDriftOfVec(const DataArrayDouble *arr, mcIdType isDrift);
    %extend
    {
      PyObject *computeVectorOfCoefficients(const MEDCouplingMesh *mesh, const DataArrayDouble *arr) const
      {
        mcIdType ret1;
        DataArrayDouble *ret0=self->computeVectorOfCoefficients(mesh,arr,ret1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,PyInt_FromLong(ret1));
        return ret;
      }
      
      PyObject *computeInverseMatrix(const MEDCouplingMesh *mesh) const
      {
        mcIdType ret1(-1),ret2(-1);
        DataArrayDouble *ret0=self->computeInverseMatrix(mesh,ret1,ret2);
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,PyInt_FromLong(ret1));
        PyTuple_SetItem(ret,2,PyInt_FromLong(ret2));
        return ret;
      }

      PyObject *computeMatrix(const MEDCouplingMesh *mesh) const
      {
        mcIdType ret1(-1),ret2(-1);
        DataArrayDouble *ret0=self->computeMatrix(mesh,ret1,ret2);
        PyObject *ret=PyTuple_New(3);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,PyInt_FromLong(ret1));
        PyTuple_SetItem(ret,2,PyInt_FromLong(ret2));
        return ret;
      }
      
      PyObject *computeEvaluationMatrixOnGivenPts(const MEDCouplingMesh *mesh, PyObject *locs) const
      {
        if(!mesh)
          throw INTERP_KERNEL::Exception("wrap of MEDCouplingFieldDiscretizationKriging::computeEvaluationMatrixOnGivenPts : input mesh is empty !");
        mcIdType sw,nbPts;
        double v0; MEDCoupling::DataArrayDouble *v1(0); MEDCoupling::DataArrayDoubleTuple *v2(0); std::vector<double> v3;
        const double *inp=convertObjToPossibleCpp5_Safe2(locs,sw,v0,v1,v2,v3,"wrap of MEDCouplingFieldDiscretizationKriging::computeEvaluationMatrixOnGivenPts",
                                                         mesh->getSpaceDimension(),true,nbPts);
        //
        mcIdType ret1(-1);
        DataArrayDouble *ret0=self->computeEvaluationMatrixOnGivenPts(mesh,inp,nbPts,ret1);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTYPE_p_MEDCoupling__DataArrayDouble, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,PyInt_FromLong(ret1));
        return ret;
      }

      void operateOnDenseMatrix(int spaceDimension, DataArrayDouble *myMatrix) const
      {
        if(!myMatrix || !myMatrix->isAllocated() || myMatrix->getNumberOfComponents()!=1)
          throw INTERP_KERNEL::Exception("Wrap of MEDCouplingFieldDiscretizationKriging::operateOnDenseMatrix : invalid input matrix as DataArrayDouble ! Must be allocated with one component !");
        self->operateOnDenseMatrix(spaceDimension,myMatrix->getNumberOfTuples(),myMatrix->getPointer());
      }

      PyObject *performDrift(const DataArrayDouble *matr, const DataArrayDouble *arr) const
      {
        mcIdType ret1(-1);
        DataArrayDouble *ret0(self->performDrift(matr,arr,ret1));
        PyObject *res(PyTuple_New(2));
        PyTuple_SetItem(res,0,SWIG_NewPointerObj((void*)ret0,SWIGTYPE_p_MEDCoupling__DataArrayDouble,SWIG_POINTER_OWN | 0));
        PyTuple_SetItem(res,1,PyInt_FromLong(ret1));
        return res;
      }

      static PyObject *PerformDriftRect(const DataArrayDouble *matr, const DataArrayDouble *arr)
      {
        mcIdType ret1(-1);
        DataArrayDouble *ret0(MEDCouplingFieldDiscretizationKriging::PerformDriftRect(matr,arr,ret1));
        PyObject *res(PyTuple_New(2));
        PyTuple_SetItem(res,0,SWIG_NewPointerObj((void*)ret0,SWIGTYPE_p_MEDCoupling__DataArrayDouble,SWIG_POINTER_OWN | 0));
        PyTuple_SetItem(res,1,PyInt_FromLong(ret1));
        return res;
      }

      static void OperateOnDenseMatrixH3(DataArrayDouble *myMatrix)
      {
        if(!myMatrix || !myMatrix->isAllocated() || myMatrix->getNumberOfComponents()!=1)
          throw INTERP_KERNEL::Exception("Wrap of MEDCouplingFieldDiscretizationKriging::OperateOnDenseMatrixH3 : invalid input matrix as DataArrayDouble ! Must be allocated with one component !");
        MEDCouplingFieldDiscretizationKriging::OperateOnDenseMatrixH3(myMatrix->getNumberOfTuples(),myMatrix->getPointer());
      }
      
      static void OperateOnDenseMatrixH2Ln(DataArrayDouble *myMatrix) //throw(INTERP_KERNEL::Exception)
      {
        if(!myMatrix || !myMatrix->isAllocated() || myMatrix->getNumberOfComponents()!=1)
          throw INTERP_KERNEL::Exception("Wrap of MEDCouplingFieldDiscretizationKriging::OperateOnDenseMatrixH2Ln : invalid input matrix as DataArrayDouble ! Must be allocated with one component !");
        MEDCouplingFieldDiscretizationKriging::OperateOnDenseMatrixH2Ln(myMatrix->getNumberOfTuples(),myMatrix->getPointer());
      }
    }
  };

  class MEDCouplingFieldDiscretizationOnNodesFE : public MEDCouplingFieldDiscretizationOnNodes
  {
    public:
    %extend
    {
      DataArrayDouble *getCooInRefElement(const MEDCouplingMesh *mesh, PyObject *locs)
      {
        mcIdType sw,nbPts;
        double v0; MEDCoupling::DataArrayDouble *v1(0); MEDCoupling::DataArrayDoubleTuple *v2(0); std::vector<double> v3;
        const double *inp=convertObjToPossibleCpp5_Safe2(locs,sw,v0,v1,v2,v3,"wrap of MEDCouplingFieldDouble::getValueOnMulti",
                                                         mesh->getSpaceDimension(),true,nbPts);
        MCAuto<DataArrayDouble> ret(self->getCooInRefElement(mesh,inp,nbPts));
        return ret.retn();
      }
    }
  };
}
