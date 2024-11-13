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
// Author : Anthony Geay (EDF R&D)

%include "MEDCouplingMemArray.i"

%define ARRAYDEF( ARRAY, INT )

// namespace MEDCoupling
// {
//   class ARRAY ## Iterator;

//   class ARRAY : public DataArray -- #ifdef doesn't work inside
//   {
  public:
    static ARRAY *New();
    INT intValue() const;
    INT getHashCode() const;
    bool empty() const;
    void aggregate(const ARRAY *other);
    ARRAY *performCopyOrIncrRef(bool deepCopy) const;
    void deepCopyFrom(const ARRAY& other);
    void reserve(std::size_t nbOfElems);
    void pushBackSilent(INT val);
    INT popBackSilent();
    void pack() const;
    void allocIfNecessary(INT nbOfTuple, INT nbOfCompo);
    bool isEqual(const ARRAY& other) const;
    bool isEqualWithoutConsideringStr(const ARRAY& other) const;
    bool isEqualWithoutConsideringStrAndOrder(const ARRAY& other) const;
    DataArrayIdType *occurenceRankInThis() const;
    DataArrayIdType *buildPermutationArr(const ARRAY& other) const;
    ARRAY *sumPerTuple() const;
    void sort(bool asc=true);
    void reverse();
    void checkMonotonic(bool increasing) const;
    bool isMonotonic(bool increasing) const;
    void checkStrictlyMonotonic(bool increasing) const;
    bool isStrictlyMonotonic(bool increasing) const;
    void fillWithZero();
    void fillWithValue(INT val);
    void iota(INT init=0);
    std::string repr() const;
    std::string reprZip() const;
    std::string reprNotTooLong() const;
    ARRAY *invertArrayO2N2N2O(mcIdType newNbOfElem) const;
    ARRAY *invertArrayN2O2O2N(mcIdType oldNbOfElem) const;
    ARRAY *invertArrayO2N2N2OBis(mcIdType newNbOfElem) const;
    DataArrayIdType *indicesOfSubPart(const ARRAY& partOfThis) const;
    ARRAY *fromNoInterlace() const;
    ARRAY *toNoInterlace() const;
    ARRAY *selectByTupleIdSafeSlice(mcIdType bg, mcIdType end, mcIdType step) const;
    DataArrayIdType *checkAndPreparePermutation() const;
    DataArrayIdType *buildPermArrPerLevel() const;
    bool isIota(mcIdType sizeExpected) const;
    bool isUniform(INT val) const;
    INT checkUniformAndGuess() const;
    bool hasUniqueValues() const;
    ARRAY *subArray(mcIdType tupleIdBg, mcIdType tupleIdEnd=-1) const;
    void transpose();
    ARRAY *changeNbOfComponents(std::size_t newNbOfComp, INT dftValue) const;
    void meldWith(const ARRAY *other);
    void setPartOfValues1(const ARRAY *a, mcIdType bgTuples, mcIdType endTuples, mcIdType stepTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp, bool strictCompoCompare=true);
    void setPartOfValuesSimple1(INT a, mcIdType bgTuples, mcIdType endTuples, mcIdType stepTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp);
    void setPartOfValuesAdv(const ARRAY *a, const DataArrayIdType *tuplesSelec);
    void getTuple(mcIdType tupleId, INT *res) const;
    INT getIJ(std::size_t tupleId, std::size_t compoId) const;
    INT getIJSafe(std::size_t tupleId, std::size_t compoId) const;
    INT front() const;
    INT back() const;
    void setIJ(mcIdType tupleId, mcIdType compoId, INT newVal);
    void setIJSilent(mcIdType tupleId, mcIdType compoId, INT newVal);
    INT *getPointer();
    const INT *getConstPointer() const;
    ARRAY ## Iterator *iterator();
    const INT *begin() const;
    const INT *end() const;
    DataArrayIdType *findIdsEqual(INT val) const;
    DataArrayIdType *findIdsNotEqual(INT val) const;
    mcIdType changeValue(INT oldValue, INT newValue);
    mcIdType findIdFirstEqualTuple(const std::vector<INT>& tupl) const;
    mcIdType findIdFirstEqual(INT value) const;
    mcIdType findIdFirstEqual(const std::vector<INT>& vals) const;
    mcIdType findIdSequence(const std::vector<INT>& vals) const;
    bool presenceOfTuple(const std::vector<INT>& tupl) const;
    bool presenceOfValue(INT value) const;
    bool presenceOfValue(const std::vector<INT>& vals) const;
    INT count(INT value) const;
    INT accumulate(INT compId) const;
    INT getMaxValueInArray() const;
    INT getMaxAbsValueInArray() const;
    INT getMinValueInArray() const;
    void abs();
    void sortPerTuple(bool asc);
    ARRAY *computeAbs() const;
    void applyLin(INT a, INT b, INT compoId);
    void applyLin(INT a, INT b);
    void applyInv(INT numerator);
    ARRAY *negate() const;
    void applyDivideBy(INT val);
    void applyModulus(INT val);
    void applyRModulus(INT val);
    void applyPow(INT val);
    void applyRPow(INT val);
    ARRAY *findIdsInRange(INT vmin, INT vmax) const;
    ARRAY *findIdsNotInRange(INT vmin, INT vmax) const;
    ARRAY *findIdsStrictlyNegative() const;
    bool checkAllIdsInRange(INT vmin, INT vmax) const;
    static ARRAY *Aggregate(const ARRAY *a1, const ARRAY *a2, INT offsetA2);
    static ARRAY *Meld(const ARRAY *a1, const ARRAY *a2);
    static DataArrayIdType *MakePartition(const std::vector<const ARRAY *>& groups, mcIdType newNb, std::vector< std::vector<mcIdType> >& fidsOfGroups);
    static ARRAY *BuildUnion(const std::vector<const ARRAY *>& arr);
    static ARRAY *BuildIntersection(const std::vector<const ARRAY *>& arr);
    static DataArrayIdType *FindPermutationFromFirstToSecond(const ARRAY *ids1, const ARRAY *ids2);
    static DataArrayIdType *FindPermutationFromFirstToSecondDuplicate(const ARRAY *ids1, const ARRAY *ids2);
    DataArrayIdType *buildComplement(mcIdType nbOfElement) const;
    ARRAY *buildSubstraction(const ARRAY *other) const;
    ARRAY *buildSubstractionOptimized(const ARRAY *other) const;
    ARRAY *buildUnion(const ARRAY *other) const;
    ARRAY *buildIntersection(const ARRAY *other) const;
    DataArrayIdType *indexOfSameConsecutiveValueGroups() const;
    ARRAY *buildUnique() const;
    ARRAY *buildUniqueNotSorted() const;
    ARRAY *deltaShiftIndex() const;
    void computeOffsets();
    void computeOffsetsFull();
    ARRAY *buildExplicitArrByRanges(const ARRAY *offsets) const;
    DataArrayIdType *findRangeIdForEachTuple(const ARRAY *ranges) const;
    ARRAY *findIdInRangeForEachTuple(const ARRAY *ranges) const;
    void sortEachPairToMakeALinkedList();
    void sortToHaveConsecutivePairs();
    ARRAY *duplicateEachTupleNTimes(mcIdType nbTimes) const;
    ARRAY *getDifferentValues() const;
    static ARRAY *Add(const ARRAY *a1, const ARRAY *a2);
    void addEqual(const ARRAY *other);
    static ARRAY *Substract(const ARRAY *a1, const ARRAY *a2);
    void substractEqual(const ARRAY *other);
    static ARRAY *Multiply(const ARRAY *a1, const ARRAY *a2);
    void multiplyEqual(const ARRAY *other);
    static ARRAY *Divide(const ARRAY *a1, const ARRAY *a2);
    void divideEqual(const ARRAY *other);
    static ARRAY *Modulus(const ARRAY *a1, const ARRAY *a2);
    void modulusEqual(const ARRAY *other);
    static ARRAY *Pow(const ARRAY *a1, const ARRAY *a2);
    void powEqual(const ARRAY *other);
    MCAuto<ARRAY> fromLinkedListOfPairToList() const;
    MCAuto<DataArrayIdType> findIdsGreaterOrEqualTo(INT val) const;
    MCAuto<DataArrayIdType> findIdsGreaterThan(INT val) const;
    MCAuto<DataArrayIdType> findIdsLowerOrEqualTo(INT val) const;
    MCAuto<DataArrayIdType> findIdsLowerThan(INT val) const;
    MCAuto<ARRAY> selectPartDef(const PartDefinition* pd) const;
    MCAuto<DataArrayDouble> convertToDblArr() const;
    MCAuto<DataArrayFloat> convertToFloatArr() const;
  public:
    static ARRAY *Range(INT begin, INT end, INT step);
    %extend
    {
      ARRAY()
        {
          return ARRAY::New();
        }

      static ARRAY *New(PyObject *elt0, PyObject *nbOfTuples=0, PyObject *nbOfComp=0)
      {
        const char *msgBase="MEDCoupling::ARRAY::New : Available API are : \n-ARRAY.New()\n-ARRAY.New([1,3,4])\n-ARRAY.New([1,3,4],3)\n-ARRAY.New([1,3,4,5],2,2)\n-ARRAY.New([1,3,4,5,7,8],3,2)\n-ARRAY.New([(1,3),(4,5),(7,8)])\n-ARRAY.New(5)\n-ARRAY.New(5,2)";
        std::string msg(msgBase);
        if ( MEDCouplingHasNumPyBindings() )
          msg+="\n-ARRAY.New(numpy array with dtype=int32)";

        msg+=" !";
        if(PyList_Check(elt0) || PyTuple_Check(elt0))
          {
            if(nbOfTuples)
              {
                if(PyInt_Check(nbOfTuples))
                  {
                    mcIdType nbOfTuples1=ToIdType(PyInt_AS_LONG(nbOfTuples));
                    if(nbOfTuples1<0)
                      throw INTERP_KERNEL::Exception("ARRAY::New : should be a positive set of allocated memory !");
                    if(nbOfComp)
                      {
                        if(PyInt_Check(nbOfComp))
                          {//ARRAY.New([1,3,4,5],2,2)
                            mcIdType nbOfCompo=ToIdType(PyInt_AS_LONG(nbOfComp));
                            if(nbOfCompo<0)
                              throw INTERP_KERNEL::Exception("ARRAY::New : should be a positive number of components !");
                            MCAuto<ARRAY> ret=ARRAY::New();
                            std::vector<INT> tmp=fillArrayWithPyListInt2<INT>(elt0,nbOfTuples1,nbOfCompo);
                            ret->alloc(nbOfTuples1,nbOfCompo); std::copy(tmp.begin(),tmp.end(),ret->getPointer());
                            return ret.retn();
                          }
                        else
                          throw INTERP_KERNEL::Exception(msg.c_str());
                      }
                    else
                      {//ARRAY.New([1,3,4],3)
                        MCAuto<ARRAY> ret=ARRAY::New();
                        mcIdType tmpp1=-1;
                        std::vector<INT> tmp=fillArrayWithPyListInt2<INT>(elt0,nbOfTuples1,tmpp1);
                        ret->alloc(nbOfTuples1,tmpp1); std::copy(tmp.begin(),tmp.end(),ret->getPointer());
                        return ret.retn();
                      }
                  }
                else
                  throw INTERP_KERNEL::Exception(msg.c_str());
              }
            else
              {// ARRAY.New([1,3,4])
                MCAuto<ARRAY> ret=ARRAY::New();
                mcIdType tmpp1=-1,tmpp2=-1;
                std::vector<INT> tmp=fillArrayWithPyListInt2<INT>(elt0,tmpp1,tmpp2);
                ret->alloc(tmpp1,tmpp2); std::copy(tmp.begin(),tmp.end(),ret->getPointer());
                return ret.retn();
              }
          }
        else if(PyInt_Check(elt0))
          {
            INT nbOfTuples1=(INT)PyInt_AS_LONG(elt0);
            if(nbOfTuples1<0)
              throw INTERP_KERNEL::Exception("ARRAY::New : should be a positive set of allocated memory !");
            if(nbOfTuples)
              {
                if(!nbOfComp)
                  {
                    if(PyInt_Check(nbOfTuples))
                      {//ARRAY.New(5,2)
                        INT nbOfCompo=(INT)PyInt_AS_LONG(nbOfTuples);
                        if(nbOfCompo<0)
                          throw INTERP_KERNEL::Exception("ARRAY::New : should be a positive number of components !");
                        MCAuto<ARRAY> ret=ARRAY::New();
                        ret->alloc(nbOfTuples1,nbOfCompo);
                        return ret.retn();
                      }
                    else
                      throw INTERP_KERNEL::Exception(msg.c_str());
                  }
                else
                  throw INTERP_KERNEL::Exception(msg.c_str());
              }
            else
              {//ARRAY.New(5)
                MCAuto<ARRAY> ret=ARRAY::New();
                ret->alloc(nbOfTuples1,1);
                return ret.retn();
              }
          }
#if defined(WITH_NUMPY)
        else if(MEDCouplingHasNumPyBindings() && PyArray_Check(elt0) && nbOfTuples==NULL && nbOfComp==NULL)
          {//ARRAY.New(numpyArray)
            return BuildNewInstance<ARRAY,INT>(elt0,NPYTraits<INT>::NPYObjectType,NPYTraits<INT>::NPYFunc,MEDCoupling::Traits<INT>::NPYStr);
          }
#endif
        else
          throw INTERP_KERNEL::Exception(msg.c_str());
        throw INTERP_KERNEL::Exception(msg.c_str());//to make g++ happy
      }

      ARRAY(PyObject *elt0, PyObject *nbOfTuples=0, PyObject *nbOfComp=0)
        {
          return MEDCoupling_ ## ARRAY ## _New__SWIG_1(elt0,nbOfTuples,nbOfComp);
        }
      
      std::string __str__() const
      {
        return self->reprNotTooLong();
      }

      mcIdType __len__() const
      {
        if(self->isAllocated())
          {
            return self->getNumberOfTuples();
          }
        else
          {
            throw INTERP_KERNEL::Exception("ARRAY::__len__ : Instance is NOT allocated !");
          }
      }

      INT __int__() const
      {
        return self->intValue();
      }

      ARRAY ## Iterator *__iter__()
      {
        return self->iterator();
      }
   
      PyObject *accumulate() const
      {
        mcIdType sz=ToIdType(self->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<INT> tmp=new INT[sz];
        self->accumulate((INT *)tmp);
        return convertIntArrToPyList((const INT *)tmp,sz);
      }

      ARRAY *accumulatePerChunck(PyObject *indexArr) const
      {
        mcIdType sw,sz,val;
        std::vector<mcIdType> val2;
        const mcIdType *bg=convertIntStarLikePyObjToCppIntStar(indexArr,sw,sz,val,val2);
        return self->accumulatePerChunck(bg,bg+sz);
      }

      DataArrayIdType *findIdsEqualTuple(PyObject *inputTuple) const
      {
        mcIdType sw,sz;
        INT val;
        std::vector<INT> val2;
        const INT *bg(convertIntStarLikePyObjToCppIntStar(inputTuple,sw,sz,val,val2));
        return self->findIdsEqualTuple(bg,bg+sz);
      }

      DataArrayIdType *findIdForEach(PyObject *vals) const
      {
        mcIdType sw,sz;
        INT val;
        std::vector<INT> val2;
        const INT *bg(convertIntStarLikePyObjToCppIntStar(vals,sw,sz,val,val2));
        MCAuto<DataArrayIdType> ret(self->findIdForEach(bg,bg+sz));
        return ret.retn();
      }

      PyObject *splitInBalancedSlices(mcIdType nbOfSlices) const
      {
        std::vector< std::pair<mcIdType,mcIdType> > slcs(self->splitInBalancedSlices(nbOfSlices));
        PyObject *ret=PyList_New(slcs.size());
        for(std::size_t i=0;i<slcs.size();i++)
          PyList_SetItem(ret,i,PySlice_New(PyInt_FromLong(slcs[i].first),PyInt_FromLong(slcs[i].second),PyInt_FromLong(1)));
        return ret;
      }

      ARRAY *buildExplicitArrOfSliceOnScaledArr(PyObject *slic) const
      {
        if(!PySlice_Check(slic))
          throw INTERP_KERNEL::Exception("ARRAY::buildExplicitArrOfSliceOnScaledArr (wrap) : expecting a pyslice as second (first) parameter !");
        Py_ssize_t strt=2,stp=2,step=2;
        GetIndicesOfSliceExplicitely(slic,&strt,&stp,&step,"ARRAY::buildExplicitArrOfSliceOnScaledArr (wrap) : the input slice is invalid !");
        if(strt==std::numeric_limits<INT>::max() || stp==std::numeric_limits<INT>::max())
          throw INTERP_KERNEL::Exception("ARRAY::buildExplicitArrOfSliceOnScaledArr (wrap) : the input slice contains some unknowns that can't be determined in static method ! Call DataArray::getSlice (non static) instead !");
        return self->buildExplicitArrOfSliceOnScaledArr((INT)strt,(INT)stp,(INT)step);
      }

      PyObject *getMinMaxValues() const
      {
        INT a,b;
        self->getMinMaxValues(a,b);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,PyInt_FromLong(a));
        PyTuple_SetItem(ret,1,PyInt_FromLong(b));
        return ret;
      }
   
      static PyObject *ConvertIndexArrayToO2N(mcIdType nbOfOldTuples, PyObject *arr, PyObject *arrI)
      {
        mcIdType newNbOfTuples=-1;
        mcIdType szArr,szArrI,sw,iTypppArr,iTypppArrI;
        std::vector<mcIdType> stdvecTyyppArr;
        std::vector<mcIdType> stdvecTyyppArrI;
        const mcIdType *arrPtr=convertIntStarLikePyObjToCppIntStar(arr,sw,szArr,iTypppArr,stdvecTyyppArr);
        const mcIdType *arrIPtr=convertIntStarLikePyObjToCppIntStar(arrI,sw,szArrI,iTypppArrI,stdvecTyyppArrI);
        DataArrayIdType *ret0=MEDCoupling::ARRAY::ConvertIndexArrayToO2N(nbOfOldTuples,arrPtr,arrIPtr,arrIPtr+szArrI,newNbOfTuples);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj((void*)ret0,SWIGTITraits<mcIdType>::TI,SWIG_POINTER_OWN | 0));
        PyTuple_SetItem(ret,1,PyInt_FromLong(newNbOfTuples));
        return ret;
      }

      static DataArrayIdType *CheckAndPreparePermutation(PyObject *arr)
      {
        MCAuto<DataArrayIdType> ret(DataArrayIdType::New());
        mcIdType szArr,sw;
        INT iTypppArr;
        std::vector<INT> stdvecTyyppArr;
        const INT *arrPtr(convertIntStarLikePyObjToCppIntStar(arr,sw,szArr,iTypppArr,stdvecTyyppArr));
        mcIdType *pt(MEDCoupling::ARRAY::CheckAndPreparePermutation(arrPtr,arrPtr+szArr));
        ret->useArray(pt,true,MEDCoupling::DeallocType::C_DEALLOC,szArr,1);
        return ret.retn();
      }

      void setValues(PyObject *li, PyObject *nbOfTuples=0, PyObject *nbOfComp=0)
      {
        const char *msg="MEDCoupling::ARRAY::setValues : Available API are : \n-ARRAY.setValues([1,3,4])\n-ARRAY.setValues([1,3,4],3)\n-ARRAY.setValues([1,3,4,5],2,2)\n-ARRAY.New(5)\n !";
        if(PyList_Check(li) || PyTuple_Check(li))
          {
            if(nbOfTuples && nbOfTuples != Py_None)
              {
                if(PyInt_Check(nbOfTuples))
                  {
                    mcIdType nbOfTuples1=ToIdType(PyInt_AS_LONG(nbOfTuples));
                    if(nbOfTuples1<0)
                      throw INTERP_KERNEL::Exception("ARRAY::setValue : should be a positive set of allocated memory !");
                    if(nbOfComp && nbOfComp != Py_None)
                      {
                        if(PyInt_Check(nbOfComp))
                          {//ARRAY.setValues([1,3,4,5],2,2)
                            mcIdType nbOfCompo=ToIdType(PyInt_AS_LONG(nbOfComp));
                            if(nbOfCompo<0)
                              throw INTERP_KERNEL::Exception("ARRAY::setValue : should be a positive number of components !");
                            std::vector<INT> tmp=fillArrayWithPyListInt2<INT>(li,nbOfTuples1,nbOfCompo);
                            self->alloc(nbOfTuples1,nbOfCompo); std::copy(tmp.begin(),tmp.end(),self->getPointer());
                          }
                        else
                          throw INTERP_KERNEL::Exception(msg);
                      }
                    else
                      {//ARRAY.setValues([1,3,4],3)
                        mcIdType tmpp1=-1;
                        std::vector<INT> tmp=fillArrayWithPyListInt2<INT>(li,nbOfTuples1,tmpp1);
                        self->alloc(nbOfTuples1,tmpp1); std::copy(tmp.begin(),tmp.end(),self->getPointer());
                      }
                  }
                else
                  throw INTERP_KERNEL::Exception(msg);
              }
            else
              {// ARRAY.setValues([1,3,4])
                mcIdType tmpp1=-1,tmpp2=-1;
                std::vector<INT> tmp=fillArrayWithPyListInt2<INT>(li,tmpp1,tmpp2);
                self->alloc(tmpp1,tmpp2); std::copy(tmp.begin(),tmp.end(),self->getPointer());
              }
          }
        else
          throw INTERP_KERNEL::Exception(msg);
      }

      PyObject *getValues() const
      {
        const INT *vals=self->getConstPointer();
        return convertIntArrToPyList(vals,self->getNbOfElems());
      }

      PyObject *isEqualIfNotWhy(const ARRAY& other) const
      {
        std::string ret1;
        bool ret0=self->isEqualIfNotWhy(other,ret1);
        PyObject *ret=PyTuple_New(2);
        PyObject *ret0Py=ret0?Py_True:Py_False;
        Py_XINCREF(ret0Py);
        PyTuple_SetItem(ret,0,ret0Py);
        PyTuple_SetItem(ret,1,PyString_FromString(ret1.c_str()));
        return ret;
      }

      PyObject *getValuesAsTuple() const
      {
        const INT *vals=self->getConstPointer();
        mcIdType nbOfComp=ToIdType(self->getNumberOfComponents());
        mcIdType nbOfTuples=self->getNumberOfTuples();
        return convertIntArrToPyListOfTuple(vals,nbOfComp,nbOfTuples);
      }

      static PyObject *MakePartition(PyObject *gps, mcIdType newNb)
      {
        std::vector<const ARRAY *> groups;
        std::vector< std::vector<mcIdType> > fidsOfGroups;
        convertFromPyObjVectorOfObj(gps,SWIGTITraits<INT>::TI,"ARRAY",groups);
        DataArrayIdType *ret0=MEDCoupling::ARRAY::MakePartition(groups,newNb,fidsOfGroups);
        PyObject *ret = PyList_New(2);
        PyList_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        std::size_t sz=fidsOfGroups.size();
        PyObject *ret1 = PyList_New(sz);
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(ret1,i,convertIntArrToPyList2(fidsOfGroups[i]));
        PyList_SetItem(ret,1,ret1);
        return ret;
      }

      DataArrayIdType *findIdsEqualList(PyObject *obj)
      {
        mcIdType sw;
        INT singleVal;
        std::vector<INT> multiVal;
        std::pair<mcIdType, std::pair<mcIdType,mcIdType> > slic;
        ARRAY *daIntTyypp=0;
        convertIntStarOrSliceLikePyObjToCpp(obj,self->getNumberOfTuples(),sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            return self->findIdsEqualList(&singleVal,&singleVal+1);
          case 2:
            return self->findIdsEqualList(&multiVal[0],&multiVal[0]+multiVal.size());
          case 4:
            return self->findIdsEqualList(daIntTyypp->begin(),daIntTyypp->end());
          default:
            throw INTERP_KERNEL::Exception("ARRAY::findIdsEqualList : unrecognized type entered, expected list of int, tuple of int or ARRAY !");
          }
      }

      DataArrayIdType *findIdsNotEqualList(PyObject *obj)
      {
        mcIdType sw;
        INT singleVal;
        std::vector<INT> multiVal;
        std::pair<mcIdType, std::pair<mcIdType,mcIdType> > slic;
        ARRAY *daIntTyypp=0;
        convertIntStarOrSliceLikePyObjToCpp(obj,self->getNumberOfTuples(),sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            return self->findIdsNotEqualList(&singleVal,&singleVal+1);
          case 2:
            return self->findIdsNotEqualList(&multiVal[0],&multiVal[0]+multiVal.size());
          case 4:
            return self->findIdsNotEqualList(daIntTyypp->begin(),daIntTyypp->end());
          default:
            throw INTERP_KERNEL::Exception("ARRAY::findIdsNotEqualList : unrecognized type entered, expected list of int, tuple of int or ARRAY !");
          }
      }

      PyObject *splitByValueRange(PyObject *li) const
      {
        ARRAY *ret0=0,*ret1=0,*ret2=0;
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTITraits<INT>::TI, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            mcIdType size;
            INTERP_KERNEL::AutoPtr<INT> tmp=convertPyToNewIntArr2<INT>(li,&size);
            self->splitByValueRange(tmp,(INT *)tmp+size,ret0,ret1,ret2);
          }
        else
          {
            ARRAY *da2=reinterpret_cast< ARRAY * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null ARRAY instance expected !");
            da2->checkAllocated();
            self->splitByValueRange(da2->begin(),da2->end(),ret0,ret1,ret2);
          }
        PyObject *ret = PyList_New(3);
        PyList_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(ret,2,SWIG_NewPointerObj(SWIG_as_voidptr(ret2),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      DataArrayIdType *transformWithIndArrR(PyObject *li) const
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTITraits<INT>::TI, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            mcIdType size;
            INTERP_KERNEL::AutoPtr<INT> tmp=convertPyToNewIntArr2<INT>(li,&size);
            return self->transformWithIndArrR(tmp,tmp+size);
          }
        else
          {
            ARRAY *da2=reinterpret_cast< ARRAY * >(da);
            return self->transformWithIndArrR(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems());
          }
      }

      ARRAY *renumberAndReduce(PyObject *li, mcIdType newNbOfTuple)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTITraits<mcIdType>::TI, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            mcIdType size;
            INTERP_KERNEL::AutoPtr<mcIdType> tmp=convertPyToNewIntArr2(li,&size);
            if(size!=self->getNumberOfTuples())
              {
                throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
              }
            return self->renumberAndReduce(tmp,newNbOfTuple);
          }
        else
          {
            DataArrayIdType *da2=reinterpret_cast< DataArrayIdType * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
            da2->checkAllocated();
            mcIdType size=self->getNumberOfTuples();
            if(size!=self->getNumberOfTuples())
              {
                throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
              }
            return self->renumberAndReduce(da2->getConstPointer(),newNbOfTuple);
          }
      }

      ARRAY *renumber(PyObject *li)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTITraits<mcIdType>::TI, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            mcIdType size;
            INTERP_KERNEL::AutoPtr<mcIdType> tmp=convertPyToNewIntArr2(li,&size);
            if(size!=self->getNumberOfTuples())
              {
                throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
              }
            return self->renumber(tmp);
          }
        else
          {
            DataArrayIdType *da2=reinterpret_cast< DataArrayIdType * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
            da2->checkAllocated();
            mcIdType size=self->getNumberOfTuples();
            if(size!=self->getNumberOfTuples())
              {
                throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
              }
            return self->renumber(da2->getConstPointer());
          }
      }

      ARRAY *renumberR(PyObject *li)
      {
        void *da=0;
        int res1=SWIG_ConvertPtr(li,&da,SWIGTITraits<mcIdType>::TI, 0 |  0 );
        if (!SWIG_IsOK(res1))
          {
            mcIdType size;
            INTERP_KERNEL::AutoPtr<mcIdType> tmp=convertPyToNewIntArr2(li,&size);
            if(size!=self->getNumberOfTuples())
              {
                throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
              }
            return self->renumberR(tmp);
          }
        else
          {
            DataArrayIdType *da2=reinterpret_cast< DataArrayIdType * >(da);
            if(!da2)
              throw INTERP_KERNEL::Exception("Not null DataArrayInt instance expected !");
            da2->checkAllocated();
            mcIdType size=self->getNumberOfTuples();
            if(size!=self->getNumberOfTuples())
              {
                throw INTERP_KERNEL::Exception("Invalid list length ! Must be equal to number of tuples !");
              }
            return self->renumberR(da2->getConstPointer());
          }
      }

      void setSelectedComponents(const ARRAY *a, PyObject *li)
      {
        std::vector<std::size_t> tmp;
        convertPyToNewIntArr3(li,tmp);
        self->setSelectedComponents(a,tmp);
      }

      PyObject *explodeComponents() const
      {
        std::vector< MCAuto<ARRAY> > retCpp(self->explodeComponents());
        std::size_t sz(retCpp.size());
        PyObject *res(PyList_New(sz));
        for(std::size_t i=0;i<sz;i++)
          PyList_SetItem(res,i,SWIG_NewPointerObj(SWIG_as_voidptr(retCpp[i].retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return res;
      }

      PyObject *getTuple(mcIdType tupleId)
      {
        mcIdType sz=ToIdType(self->getNumberOfComponents());
        INTERP_KERNEL::AutoPtr<INT> tmp=new INT[sz];
        self->getTuple(tupleId,tmp);
        return convertIntArrToPyList((const INT*)tmp,sz);
      }

      PyObject *changeSurjectiveFormat(INT targetNb) const
      {
        DataArrayIdType *arr=0;
        DataArrayIdType *arrI=0;
        self->changeSurjectiveFormat(targetNb,arr,arrI);
        PyObject *res = PyList_New(2);
        PyList_SetItem(res,0,SWIG_NewPointerObj((void*)arr,SWIGTITraits<mcIdType>::TI,SWIG_POINTER_OWN | 0));
        PyList_SetItem(res,1,SWIG_NewPointerObj((void*)arrI,SWIGTITraits<mcIdType>::TI,SWIG_POINTER_OWN | 0));
        return res;
      }

      static ARRAY *Meld(PyObject *li)
      {
        std::vector<const ARRAY *> tmp;
        convertFromPyObjVectorOfObj(li,SWIGTITraits<INT>::TI,"ARRAY",tmp);
        return ARRAY::Meld(tmp);
      }

      static ARRAY *Aggregate(PyObject *li)
      {
        std::vector<const ARRAY *> tmp;
        convertFromPyObjVectorOfObj(li,SWIGTITraits<INT>::TI,"ARRAY",tmp);
        return ARRAY::Aggregate(tmp);
      }

      static ARRAY *AggregateIndexes(PyObject *li)
      {
        std::vector<const ARRAY *> tmp;
        convertFromPyObjVectorOfObj(li,SWIGTITraits<INT>::TI,"ARRAY",tmp);
        return ARRAY::AggregateIndexes(tmp);
      }

      static ARRAY *BuildUnion(PyObject *li)
      {
        std::vector<const ARRAY *> tmp;
        convertFromPyObjVectorOfObj(li,SWIGTITraits<INT>::TI,"ARRAY",tmp);
        return ARRAY::BuildUnion(tmp);
      }

      static ARRAY *BuildIntersection(PyObject *li)
      {
        std::vector<const ARRAY *> tmp;
        convertFromPyObjVectorOfObj(li,SWIGTITraits<INT>::TI,"ARRAY",tmp);
        return ARRAY::BuildIntersection(tmp);
      }

      PyObject *getMaxValue() const
      {
        mcIdType tmp;
        INT r1=self->getMaxValue(tmp);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,PyInt_FromLong(r1));
        PyTuple_SetItem(ret,1,PyInt_FromLong(tmp));
        return ret;
      }
    
      PyObject *getMaxAbsValue(std::size_t& tupleId) const
      {
      	std::size_t tmp;
        INT r1=self->getMaxAbsValue(tmp);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,PyInt_FromLong(r1));
        PyTuple_SetItem(ret,1,PyInt_FromLong(tmp));
        return ret;
      }

      PyObject *getMinValue() const
      {
        mcIdType tmp;
        INT r1=self->getMinValue(tmp);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,PyInt_FromLong(r1));
        PyTuple_SetItem(ret,1,PyInt_FromLong(tmp));
        return ret;
      }

      mcIdType index(PyObject *obj) const
      {
        std::size_t nbOfCompo=self->getNumberOfComponents();
        switch(nbOfCompo)
          {
          case 1:
            {
              if(PyInt_Check(obj))
                {
                  INT val=(INT)PyInt_AS_LONG(obj);
                  return self->findIdFirstEqual(val);
                }
              else
                throw INTERP_KERNEL::Exception("ARRAY::index : 'this' contains one component and trying to find an element which is not an integer !");
            }
          default:
            {
              std::vector<INT> arr;
              convertPyToNewIntArr3(obj,arr);
              return self->findIdFirstEqualTuple(arr);
            }
          }
      }

      bool __contains__(PyObject *obj) const
      {
        std::size_t nbOfCompo=self->getNumberOfComponents();
        switch(nbOfCompo)
          {
          case 0:
            return false;
          case 1:
            {
              if(PyInt_Check(obj))
                {
                  INT val=(INT)PyInt_AS_LONG(obj);
                  return self->presenceOfValue(val);
                }
              else
                throw INTERP_KERNEL::Exception("ARRAY::__contains__ : 'this' contains one component and trying to find an element which is not an integer !");
            }
          default:
            {
              std::vector<INT> arr;
              convertPyToNewIntArr3(obj,arr);
              return self->presenceOfTuple(arr);
            }
          }
      }

      PyObject *__getitem__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in ARRAY::__getitem__ !";
        const char msg2[]="ARRAY::__getitem__ : Mismatch of slice values in 2nd parameter (components) !";
        self->checkAllocated();
        mcIdType nbOfTuples=self->getNumberOfTuples();
        std::size_t nbOfComponents=self->getNumberOfComponents();
        mcIdType it1;
        std::size_t ic1;
        std::vector<mcIdType> vt1;
        std::vector<std::size_t> vc1;
        std::pair<mcIdType, std::pair<mcIdType,mcIdType> > pt1,pc1;
        DataArrayIdType *dt1=0,*dc1=0;
        mcIdType sw;
        convertObjToPossibleCpp3(obj,nbOfTuples,(int)nbOfComponents,sw,it1,ic1,vt1,vc1,pt1,pc1,dt1,dc1);
        MCAuto<ARRAY> ret;
        switch(sw)
          {
          case 1:
            {
              if(nbOfComponents==1)
                return PyInt_FromLong(self->getIJSafe(it1,0));
              return SWIG_NewPointerObj(SWIG_as_voidptr(self->selectByTupleIdSafe(&it1,&it1+1)),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 );
            }
          case 2:
            return SWIG_NewPointerObj(SWIG_as_voidptr(self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size())),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 );
          case 3:
            return SWIG_NewPointerObj(SWIG_as_voidptr(self->selectByTupleIdSafeSlice(pt1.first,pt1.second.first,pt1.second.second)),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 );
          case 4:
            return SWIG_NewPointerObj(SWIG_as_voidptr(self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems())),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 );
          case 5:
            return PyInt_FromLong(self->getIJSafe(it1,ic1));
          case 6:
            {
              ret=self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size());
              std::vector<std::size_t> v2(1,ic1);
              return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(v2)),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 );
            }
          case 7:
            {
              ret=self->selectByTupleIdSafeSlice(pt1.first,pt1.second.first,pt1.second.second);
              std::vector<std::size_t> v2(1,ic1);
              return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(v2)),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 );
            }
          case 8:
            {
              ret=self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems());
              std::vector<std::size_t> v2(1,ic1);
              return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(v2)),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 );
            }
          case 9:
            {
              ret=self->selectByTupleIdSafe(&it1,&it1+1);
              return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(vc1)),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 );
            }
          case 10:
            {
              ret=self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size());
              return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(vc1)),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 );
            }
          case 11:
            {
              ret=self->selectByTupleIdSafeSlice(pt1.first,pt1.second.first,pt1.second.second);
              return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(vc1)),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 );
            }
          case 12:
            {
              ret=self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems());
              return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(vc1)),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 );
            }
          case 13:
            {
              ret=self->selectByTupleIdSafe(&it1,&it1+1);
              mcIdType nbOfComp=DataArray::GetNumberOfItemGivenBESRelative(pc1.first,pc1.second.first,pc1.second.second,msg2);
              std::vector<std::size_t> v2(nbOfComp);
              for(INT i=0;i<nbOfComp;i++)
                v2[i]=pc1.first+i*pc1.second.second;
              return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(v2)),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 );
            }
          case 14:
            {
              ret=self->selectByTupleIdSafe(&vt1[0],&vt1[0]+vt1.size());
              mcIdType nbOfComp=DataArray::GetNumberOfItemGivenBESRelative(pc1.first,pc1.second.first,pc1.second.second,msg2);
              std::vector<std::size_t> v2(nbOfComp);
              for(INT i=0;i<nbOfComp;i++)
                v2[i]=pc1.first+i*pc1.second.second;
              return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(v2)),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 );
            }
          case 15:
            {
              ret=self->selectByTupleIdSafeSlice(pt1.first,pt1.second.first,pt1.second.second);
              mcIdType nbOfComp=DataArray::GetNumberOfItemGivenBESRelative(pc1.first,pc1.second.first,pc1.second.second,msg2);
              std::vector<std::size_t> v2(nbOfComp);
              for(mcIdType i=0;i<nbOfComp;i++)
                v2[i]=pc1.first+i*pc1.second.second;
              return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(v2)),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 );
            }
          case 16:
            {
              ret=self->selectByTupleIdSafe(dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems());
              mcIdType nbOfComp=DataArray::GetNumberOfItemGivenBESRelative(pc1.first,pc1.second.first,pc1.second.second,msg2);
              std::vector<std::size_t> v2(nbOfComp);
              for(INT i=0;i<nbOfComp;i++)
                v2[i]=pc1.first+i*pc1.second.second;
              return SWIG_NewPointerObj(SWIG_as_voidptr(ret->keepSelectedComponents(v2)),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 );
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      ARRAY *__setitem__(PyObject *obj, PyObject *value)
      {
        self->checkAllocated();
        const char msg[]="Unexpected situation in __setitem__ !";
        mcIdType nbOfTuples=self->getNumberOfTuples();
        int nbOfComponents=(int)self->getNumberOfComponents();
        mcIdType sw1,sw2;
        INT i1;
        std::vector<INT> v1;
        ARRAY *d1=0;
        ARRAY ## Tuple *dd1=0;
        convertIntStarLikePyObjToCpp(value,sw1,i1,v1,d1,dd1);
        mcIdType it1,ic1;
        std::vector<mcIdType> vt1,vc1;
        std::pair<mcIdType, std::pair<mcIdType,mcIdType> > pt1,pc1;
        DataArrayIdType *dt1=0,*dc1=0;
        convertObjToPossibleCpp3(obj,nbOfTuples,nbOfComponents,sw2,it1,ic1,vt1,vc1,pt1,pc1,dt1,dc1);
        MCAuto<ARRAY> tmp;
        switch(sw2)
          {
          case 1:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple1(i1,it1,it1+1,1,0,nbOfComponents,1);
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues1(tmp,it1,it1+1,1,0,nbOfComponents,1,false);
                  return self;
                case 3:
                  self->setPartOfValues1(d1,it1,it1+1,1,0,nbOfComponents,1);
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues1(tmp,it1,it1+1,1,0,nbOfComponents,1);
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          case 2:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple3(i1,&vt1[0],&vt1[0]+vt1.size(),0,nbOfComponents,1);
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues3(tmp,&vt1[0],&vt1[0]+vt1.size(),0,nbOfComponents,1,false);
                  return self;
                case 3:
                  self->setPartOfValues3(d1,&vt1[0],&vt1[0]+vt1.size(),0,nbOfComponents,1);
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues3(tmp,&vt1[0],&vt1[0]+vt1.size(),0,nbOfComponents,1);
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          case 3:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple1(i1,pt1.first,pt1.second.first,pt1.second.second,0,nbOfComponents,1);
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues1(tmp,pt1.first,pt1.second.first,pt1.second.second,0,nbOfComponents,1,false);
                  return self;
                case 3:
                  self->setPartOfValues1(d1,pt1.first,pt1.second.first,pt1.second.second,0,nbOfComponents,1);
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues1(tmp,pt1.first,pt1.second.first,pt1.second.second,0,nbOfComponents,1);
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          case 4:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple3(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),0,nbOfComponents,1);
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues3(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),0,nbOfComponents,1,false);
                  return self;
                case 3:
                  self->setPartOfValues3(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),0,nbOfComponents,1);
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues3(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),0,nbOfComponents,1);
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          case 5:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple1(i1,it1,it1+1,1,ic1,ic1+1,1);
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues1(tmp,it1,it1+1,1,ic1,ic1+1,1,false);
                  return self;
                case 3:
                  self->setPartOfValues1(d1,it1,it1+1,1,ic1,ic1+1,1);
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues1(tmp,it1,it1+1,1,ic1,ic1+1,1);
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          case 6:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple3(i1,&vt1[0],&vt1[0]+vt1.size(),ic1,ic1+1,1);
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues3(tmp,&vt1[0],&vt1[0]+vt1.size(),ic1,ic1+1,1,false);
                  return self;
                case 3:
                  self->setPartOfValues3(d1,&vt1[0],&vt1[0]+vt1.size(),ic1,ic1+1,1);
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues3(tmp,&vt1[0],&vt1[0]+vt1.size(),ic1,ic1+1,1);
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          case 7:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple1(i1,pt1.first,pt1.second.first,pt1.second.second,ic1,ic1+1,1);
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues1(tmp,pt1.first,pt1.second.first,pt1.second.second,ic1,ic1+1,1,false);
                  return self;
                case 3:
                  self->setPartOfValues1(d1,pt1.first,pt1.second.first,pt1.second.second,ic1,ic1+1,1);
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues1(tmp,pt1.first,pt1.second.first,pt1.second.second,ic1,ic1+1,1);
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          case 8:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple3(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),ic1,ic1+1,1);
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues3(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),ic1,ic1+1,1,false);
                  return self;
                case 3:
                  self->setPartOfValues3(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),ic1,ic1+1,1);
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues3(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),ic1,ic1+1,1);
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          case 9:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple2(i1,&it1,&it1+1,&vc1[0],&vc1[0]+vc1.size());
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues2(tmp,&it1,&it1+1,&vc1[0],&vc1[0]+vc1.size(),false);
                  return self;
                case 3:
                  self->setPartOfValues2(d1,&it1,&it1+1,&vc1[0],&vc1[0]+vc1.size());
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues2(tmp,&it1,&it1+1,&vc1[0],&vc1[0]+vc1.size());
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          case 10:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple2(i1,&vt1[0],&vt1[0]+vt1.size(),&vc1[0],&vc1[0]+vc1.size());
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues2(tmp,&vt1[0],&vt1[0]+vt1.size(),&vc1[0],&vc1[0]+vc1.size(),false);
                  return self;
                case 3:
                  self->setPartOfValues2(d1,&vt1[0],&vt1[0]+vt1.size(),&vc1[0],&vc1[0]+vc1.size());
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues2(tmp,&vt1[0],&vt1[0]+vt1.size(),&vc1[0],&vc1[0]+vc1.size());
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          case 11:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple4(i1,pt1.first,pt1.second.first,pt1.second.second,&vc1[0],&vc1[0]+vc1.size());
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues4(tmp,pt1.first,pt1.second.first,pt1.second.second,&vc1[0],&vc1[0]+vc1.size(),false);
                  return self;
                case 3:
                  self->setPartOfValues4(d1,pt1.first,pt1.second.first,pt1.second.second,&vc1[0],&vc1[0]+vc1.size());
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues4(tmp,pt1.first,pt1.second.first,pt1.second.second,&vc1[0],&vc1[0]+vc1.size());
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          case 12:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple2(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),&vc1[0],&vc1[0]+vc1.size());
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues2(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),&vc1[0],&vc1[0]+vc1.size(),false);
                  return self;
                case 3:
                  self->setPartOfValues2(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),&vc1[0],&vc1[0]+vc1.size());
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues2(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),&vc1[0],&vc1[0]+vc1.size());
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          case 13:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple1(i1,it1,it1+1,1,pc1.first,pc1.second.first,pc1.second.second);
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues1(tmp,it1,it1+1,1,pc1.first,pc1.second.first,pc1.second.second,false);
                  return self;
                case 3:
                  self->setPartOfValues1(d1,it1,it1+1,1,pc1.first,pc1.second.first,pc1.second.second);
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues1(tmp,it1,it1+1,1,pc1.first,pc1.second.first,pc1.second.second);
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          case 14:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple3(i1,&vt1[0],&vt1[0]+vt1.size(),pc1.first,pc1.second.first,pc1.second.second);
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues3(tmp,&vt1[0],&vt1[0]+vt1.size(),pc1.first,pc1.second.first,pc1.second.second,false);
                  return self;
                case 3:
                  self->setPartOfValues3(d1,&vt1[0],&vt1[0]+vt1.size(),pc1.first,pc1.second.first,pc1.second.second);
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues3(tmp,&vt1[0],&vt1[0]+vt1.size(),pc1.first,pc1.second.first,pc1.second.second);
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          case 15:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple1(i1,pt1.first,pt1.second.first,pt1.second.second,pc1.first,pc1.second.first,pc1.second.second);
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues1(tmp,pt1.first,pt1.second.first,pt1.second.second,pc1.first,pc1.second.first,pc1.second.second,false);
                  return self;
                case 3:
                  self->setPartOfValues1(d1,pt1.first,pt1.second.first,pt1.second.second,pc1.first,pc1.second.first,pc1.second.second);
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues1(tmp,pt1.first,pt1.second.first,pt1.second.second,pc1.first,pc1.second.first,pc1.second.second);
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          case 16:
            {
              switch(sw1)
                {
                case 1:
                  self->setPartOfValuesSimple3(i1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),pc1.first,pc1.second.first,pc1.second.second);
                  return self;
                case 2:
                  tmp=ARRAY::New();
                  tmp->useArray(&v1[0],false,DeallocType::CPP_DEALLOC,1,v1.size());
                  self->setPartOfValues3(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),pc1.first,pc1.second.first,pc1.second.second,false);
                  return self;
                case 3:
                  self->setPartOfValues3(d1,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),pc1.first,pc1.second.first,pc1.second.second);
                  return self;
                case 4:
                  tmp=dd1->buildDAInt(1,self->getNumberOfComponents());
                  self->setPartOfValues3(tmp,dt1->getConstPointer(),dt1->getConstPointer()+dt1->getNbOfElems(),pc1.first,pc1.second.first,pc1.second.second);
                  return self;
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
              break;
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
        return self;
      }

      ARRAY *__neg__() const
      {
        return self->negate();
      }
 
      ARRAY *__add__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in __add__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              MCAuto<ARRAY> ret=self->deepCopy();
              ret->applyLin(1,val);
              return ret.retn();
            }
          case 2:
            {
              MCAuto<ARRAY> aaaa=ARRAY::New(); aaaa->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              return ARRAY::Add(self,aaaa);
            }
          case 3:
            {
              return ARRAY::Add(self,a);
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              return ARRAY::Add(self,aaaa);
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      ARRAY *__radd__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in __radd__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              MCAuto<ARRAY> ret=self->deepCopy();
              ret->applyLin(1,val);
              return ret.retn();
            }
          case 2:
            {
              MCAuto<ARRAY> aaaa=ARRAY::New(); aaaa->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              return ARRAY::Add(self,aaaa);
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              return ARRAY::Add(self,aaaa);
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      PyObject *___iadd___(PyObject *trueSelf, PyObject *obj)
      {
        const char msg[]="Unexpected situation in __iadd__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              self->applyLin(1,val);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 2:
            {
              MCAuto<ARRAY> bb=ARRAY::New(); bb->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              self->addEqual(bb);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              self->addEqual(a);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              self->addEqual(aaaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      ARRAY *__sub__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in __sub__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              MCAuto<ARRAY> ret=self->deepCopy();
              ret->applyLin(1,-val);
              return ret.retn();
            }
          case 2:
            {
              MCAuto<ARRAY> aaaa=ARRAY::New(); aaaa->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              return ARRAY::Substract(self,aaaa);
            }
          case 3:
            {
              return ARRAY::Substract(self,a);
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              return ARRAY::Substract(self,aaaa);
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      ARRAY *__rsub__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in __rsub__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              MCAuto<ARRAY> ret=self->deepCopy();
              ret->applyLin(-1,val);
              return ret.retn();
            }
          case 2:
            {
              MCAuto<ARRAY> aaaa=ARRAY::New(); aaaa->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              return ARRAY::Substract(aaaa,self);
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              return ARRAY::Substract(aaaa,self);
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      PyObject *___isub___(PyObject *trueSelf, PyObject *obj)
      {
        const char msg[]="Unexpected situation in __isub__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              self->applyLin(1,-val);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 2:
            {
              MCAuto<ARRAY> bb=ARRAY::New(); bb->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              self->substractEqual(bb);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              self->substractEqual(a);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              self->substractEqual(aaaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      ARRAY *__mul__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in __mul__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              MCAuto<ARRAY> ret=self->deepCopy();
              ret->applyLin(val,0);
              return ret.retn();
            }
          case 2:
            {
              MCAuto<ARRAY> aaaa=ARRAY::New(); aaaa->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              return ARRAY::Multiply(self,aaaa);
            }
          case 3:
            {
              return ARRAY::Multiply(self,a);
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              return ARRAY::Multiply(self,aaaa);
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      ARRAY *__rmul__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in __rmul__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              MCAuto<ARRAY> ret=self->deepCopy();
              ret->applyLin(val,0);
              return ret.retn();
            }
          case 2:
            {
              MCAuto<ARRAY> aaaa=ARRAY::New(); aaaa->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              return ARRAY::Multiply(self,aaaa);
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              return ARRAY::Multiply(self,aaaa);
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      PyObject *___imul___(PyObject *trueSelf, PyObject *obj)
      {
        const char msg[]="Unexpected situation in __imul__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              self->applyLin(val,0);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 2:
            {
              MCAuto<ARRAY> bb=ARRAY::New(); bb->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              self->multiplyEqual(bb);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              self->multiplyEqual(a);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              self->multiplyEqual(aaaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      ARRAY *__div__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in __div__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              MCAuto<ARRAY> ret=self->deepCopy();
              ret->applyDivideBy(val);
              return ret.retn();
            }
          case 2:
            {
              MCAuto<ARRAY> aaaa=ARRAY::New(); aaaa->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              return ARRAY::Divide(self,aaaa);
            }
          case 3:
            {
              return ARRAY::Divide(self,a);
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              return ARRAY::Divide(self,aaaa);
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      ARRAY *__rdiv__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in __rdiv__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              MCAuto<ARRAY> ret=self->deepCopy();
              ret->applyInv(val);
              return ret.retn();
            }
          case 2:
            {
              MCAuto<ARRAY> aaaa=ARRAY::New(); aaaa->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              return ARRAY::Divide(aaaa,self);
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              return ARRAY::Divide(aaaa,self);
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      PyObject *___idiv___(PyObject *trueSelf, PyObject *obj)
      {
        const char msg[]="Unexpected situation in __idiv__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              self->applyDivideBy(val);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 2:
            {
              MCAuto<ARRAY> bb=ARRAY::New(); bb->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              self->divideEqual(bb);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              self->divideEqual(a);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              self->divideEqual(aaaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      ARRAY *__mod__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in __mod__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              MCAuto<ARRAY> ret=self->deepCopy();
              ret->applyModulus(val);
              return ret.retn();
            }
          case 2:
            {
              MCAuto<ARRAY> aaaa=ARRAY::New(); aaaa->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              return ARRAY::Modulus(self,aaaa);
            }
          case 3:
            {
              return ARRAY::Modulus(self,a);
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              return ARRAY::Modulus(self,aaaa);
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      ARRAY *__rmod__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in __rmod__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              MCAuto<ARRAY> ret=self->deepCopy();
              ret->applyRModulus(val);
              return ret.retn();
            }
          case 2:
            {
              MCAuto<ARRAY> aaaa=ARRAY::New(); aaaa->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              return ARRAY::Modulus(aaaa,self);
            }
          case 3:
            {
              return ARRAY::Modulus(a,self);
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              return ARRAY::Modulus(aaaa,self);
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      PyObject *___imod___(PyObject *trueSelf, PyObject *obj)
      {
        const char msg[]="Unexpected situation in __imod__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              self->applyModulus(val);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              self->modulusEqual(a);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              self->modulusEqual(aaaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      ARRAY *__pow__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in __pow__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              MCAuto<ARRAY> ret=self->deepCopy();
              ret->applyPow(val);
              return ret.retn();
            }
          case 2:
            {
              MCAuto<ARRAY> aaaa=ARRAY::New(); aaaa->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              return ARRAY::Pow(self,aaaa);
            }
          case 3:
            {
              return ARRAY::Pow(self,a);
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              return ARRAY::Pow(self,aaaa);
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      ARRAY *__rpow__(PyObject *obj)
      {
        const char msg[]="Unexpected situation in __rpow__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              MCAuto<ARRAY> ret=self->deepCopy();
              ret->applyRPow(val);
              return ret.retn();
            }
          case 2:
            {
              MCAuto<ARRAY> aaaa=ARRAY::New(); aaaa->useArray(&aa[0],false,DeallocType::CPP_DEALLOC,1,aa.size());
              return ARRAY::Pow(aaaa,self);
            }
          case 3:
            {
              return ARRAY::Pow(a,self);
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              return ARRAY::Pow(aaaa,self);
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }
   
      PyObject *___ipow___(PyObject *trueSelf, PyObject *obj)
      {
        const char msg[]="Unexpected situation in __ipow__ !";
        INT val;
        ARRAY *a;
        std::vector<INT> aa;
        ARRAY ## Tuple *aaa;
        mcIdType sw;
        convertIntStarLikePyObjToCpp(obj,sw,val,aa,a,aaa);
        switch(sw)
          {
          case 1:
            {
              self->applyPow(val);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 3:
            {
              self->powEqual(a);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          case 4:
            {
              MCAuto<ARRAY> aaaa=aaa->buildDAInt(1,self->getNumberOfComponents());
              self->powEqual(aaaa);
              Py_XINCREF(trueSelf);
              return trueSelf;
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }

      std::string __repr__() const
      {
        std::ostringstream oss;
        self->reprQuickOverview(oss);
        return oss.str();
      }
      
      void pushBackValsSilent(PyObject *li)
      {
        mcIdType szArr,sw;
        INT iTypppArr;
        std::vector<INT> stdvecTyyppArr;
        const INT *tmp=convertIntStarLikePyObjToCppIntStar(li,sw,szArr,iTypppArr,stdvecTyyppArr);
        self->pushBackValsSilent(tmp,tmp+szArr);
      }
      
      PyObject *partitionByDifferentValues() const
      {
        std::vector<INT> ret1;
        std::vector<DataArrayIdType *> ret0=self->partitionByDifferentValues(ret1);
        std::size_t sz=ret0.size();
        PyObject *pyRet=PyTuple_New(2);
        PyObject *pyRet0=PyList_New((INT)sz);
        PyObject *pyRet1=PyList_New((INT)sz);
        for(std::size_t i=0;i<sz;i++)
          {
            PyList_SetItem(pyRet0,i,SWIG_NewPointerObj(SWIG_as_voidptr(ret0[i]),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
            PyList_SetItem(pyRet1,i,PyInt_FromLong(ret1[i]));
          }
        PyTuple_SetItem(pyRet,0,pyRet0);
        PyTuple_SetItem(pyRet,1,pyRet1);
        return pyRet;
      }
      
      DataArrayIdType *locateComponentId(const ARRAY *valToSearchIntoTuples, const DataArrayIdType *tupleIdHint) const
      {
        return self->locateComponentId(valToSearchIntoTuples,tupleIdHint);
      }
      
      PyObject *findIdsRangesInListOfIds(const ARRAY *listOfIds) const
      {
        DataArrayIdType *ret0=0;
        ARRAY *ret1=0;
        self->findIdsRangesInListOfIds(listOfIds,ret0,ret1);
        PyObject *pyRet=PyTuple_New(2);
        PyTuple_SetItem(pyRet,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(pyRet,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 ));
        return pyRet;
      }
      
      PyObject *forThisAsPartitionBuildReduction(DataArrayIdType *commonEntities, DataArrayIdType *commonEntitiesIndex) const
      {
        MCAuto<DataArrayIdType> ceCpp( MCAuto<DataArrayIdType>::TakeRef(commonEntities) );
        MCAuto<DataArrayIdType> ceiCpp( MCAuto<DataArrayIdType>::TakeRef(commonEntitiesIndex) );
        MCAuto<ARRAY> ret1;
        MCAuto<DataArrayIdType> ret2;
        MCAuto<ARRAY> ret0 = self->forThisAsPartitionBuildReduction(ceCpp,ceiCpp,ret1,ret2);
        PyObject *pyRet( PyTuple_New(3) );
        PyTuple_SetItem(pyRet,0,SWIG_NewPointerObj(SWIG_as_voidptr(ret0.retn()),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(pyRet,1,SWIG_NewPointerObj(SWIG_as_voidptr(ret1.retn()),SWIGTITraits<INT>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(pyRet,2,SWIG_NewPointerObj(SWIG_as_voidptr(ret2.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return pyRet;
      }

      PyObject *isRange() const
      {
        INT a(0),b(0),c(0);
        bool ret(self->isRange(a,b,c));
        PyObject *pyRet=PyTuple_New(2);
        PyObject *ret0Py=ret?Py_True:Py_False,*ret1Py(0);
        Py_XINCREF(ret0Py);
        PyTuple_SetItem(pyRet,0,ret0Py);
        if(ret)
          ret1Py=PySlice_New(PyInt_FromLong(a),PyInt_FromLong(b),PyInt_FromLong(c));
        else
          {
            ret1Py=Py_None;
            Py_XINCREF(ret1Py);
          }
        PyTuple_SetItem(pyRet,1,ret1Py);
        return pyRet;
      }

      static bool RemoveIdsFromIndexedArrays(PyObject *li, ARRAY *arr, DataArrayIdType *arrIndx, mcIdType offsetForRemoval=0) throw(INTERP_KERNEL::Exception)
      {
        mcIdType sw;
        INT singleVal;
        std::vector<INT> multiVal;
        std::pair<mcIdType, std::pair<mcIdType,mcIdType> > slic;
        MEDCoupling::ARRAY *daIntTyypp=0;
        if(!arrIndx)
          throw INTERP_KERNEL::Exception("ARRAY::RemoveIdsFromIndexedArrays : null pointer as arrIndex !");
        convertIntStarOrSliceLikePyObjToCpp(li,arrIndx->getNumberOfTuples()-1,sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            return ARRAY::RemoveIdsFromIndexedArrays(&singleVal,&singleVal+1,arr,arrIndx,offsetForRemoval);
          case 2:
            return ARRAY::RemoveIdsFromIndexedArrays(&multiVal[0],&multiVal[0]+multiVal.size(),arr,arrIndx,offsetForRemoval);
          case 4:
            return ARRAY::RemoveIdsFromIndexedArrays(daIntTyypp->begin(),daIntTyypp->end(),arr,arrIndx,offsetForRemoval);
          default:
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::RemoveIdsFromIndexedArrays : unrecognized type entered, expected list of int, tuple of int or ARRAY !");
          }
      }

      PyObject *findCommonTuples(mcIdType limitNodeId=-1) const
      {
        MCAuto<DataArrayIdType> comm,commIndex;
        self->findCommonTuples(limitNodeId,comm,commIndex);
        PyObject *res = PyList_New(2);
        PyList_SetItem(res,0,SWIG_NewPointerObj(SWIG_as_voidptr(comm.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyList_SetItem(res,1,SWIG_NewPointerObj(SWIG_as_voidptr(commIndex.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return res;
      }

      static PyObject *ExtractFromIndexedArrays(PyObject *li, const ARRAY *arrIn, const DataArrayIdType *arrIndxIn) throw(INTERP_KERNEL::Exception)
      {
        ARRAY *arrOut=0;
        DataArrayIdType *arrIndexOut=0;
        mcIdType sw;
        mcIdType singleVal;
        std::vector<mcIdType> multiVal;
        std::pair<mcIdType, std::pair<mcIdType,mcIdType> > slic;
        MEDCoupling::DataArrayIdType *daIntTyypp=0;
        if(!arrIndxIn)
          throw INTERP_KERNEL::Exception("ARRAY::ExtractFromIndexedArrays : null pointer as arrIndxIn !");
        convertIntStarOrSliceLikePyObjToCpp(li,arrIndxIn->getNumberOfTuples()-1,sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            {
              ARRAY::ExtractFromIndexedArrays(&singleVal,&singleVal+1,arrIn,arrIndxIn,arrOut,arrIndexOut);
              break;
            }
          case 2:
            {
              ARRAY::ExtractFromIndexedArrays(&multiVal[0],&multiVal[0]+multiVal.size(),arrIn,arrIndxIn,arrOut,arrIndexOut);
              break;
            }
          case 4:
            {
              ARRAY::ExtractFromIndexedArrays(daIntTyypp->begin(),daIntTyypp->end(),arrIn,arrIndxIn,arrOut,arrIndexOut);
              break;
            }
          default:
            throw INTERP_KERNEL::Exception("ARRAY::ExtractFromIndexedArrays : unrecognized type entered, expected list of int, tuple of int or ARRAY !");
          }
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(arrOut),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(arrIndexOut),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *ExtractFromIndexedArraysSlice(mcIdType strt, mcIdType stp, mcIdType step, const ARRAY *arrIn, const DataArrayIdType *arrIndxIn) throw(INTERP_KERNEL::Exception)
      {
        ARRAY *arrOut=0;
        DataArrayIdType *arrIndexOut=0;
        ARRAY::ExtractFromIndexedArraysSlice(strt,stp,step,arrIn,arrIndxIn,arrOut,arrIndexOut);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(arrOut),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(arrIndexOut),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *ExtractFromIndexedArraysSlice(PyObject *slic, const ARRAY *arrIn, const DataArrayIdType *arrIndxIn) throw(INTERP_KERNEL::Exception)
      {
        if(!PySlice_Check(slic))
          throw INTERP_KERNEL::Exception("ExtractFromIndexedArraysSlice (wrap) : the first param is not a pyslice !");
        Py_ssize_t strt=2,stp=2,step=2;
        if(!arrIndxIn)
          throw INTERP_KERNEL::Exception("ExtractFromIndexedArraysSlice (wrap) : last array is null !");
        arrIndxIn->checkAllocated();
        if(arrIndxIn->getNumberOfComponents()!=1)
          throw INTERP_KERNEL::Exception("ExtractFromIndexedArraysSlice (wrap) : number of components of last argument must be equal to one !");
        GetIndicesOfSlice(slic,arrIndxIn->getNumberOfTuples(),&strt,&stp,&step,"ExtractFromIndexedArraysSlice (wrap) : Invalid slice regarding nb of elements !");
        ARRAY *arrOut=0;
        DataArrayIdType *arrIndexOut=0;
        ARRAY::ExtractFromIndexedArraysSlice(ToIdType(strt),ToIdType(stp),ToIdType(step),arrIn,arrIndxIn,arrOut,arrIndexOut);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(arrOut),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(arrIndexOut),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static PyObject *SetPartOfIndexedArrays(PyObject *li,
                                              const ARRAY *arrIn, const DataArrayIdType *arrIndxIn,
                                              const ARRAY *srcArr, const DataArrayIdType *srcArrIndex) throw(INTERP_KERNEL::Exception)
      {
        ARRAY *arrOut=0;
        DataArrayIdType *arrIndexOut=0;
        mcIdType sw;
        mcIdType singleVal;
        std::vector<mcIdType> multiVal;
        std::pair<mcIdType, std::pair<mcIdType,mcIdType> > slic;
        MEDCoupling::DataArrayIdType *daIntTyypp=0;
        if(!arrIndxIn)
          throw INTERP_KERNEL::Exception("ARRAY::SetPartOfIndexedArrays : null pointer as arrIndex !");
        convertIntStarOrSliceLikePyObjToCpp(li,arrIndxIn->getNumberOfTuples()-1,sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            {
              ARRAY::SetPartOfIndexedArrays(&singleVal,&singleVal+1,arrIn,arrIndxIn,srcArr,srcArrIndex,arrOut,arrIndexOut);
              break;
            }
          case 2:
            {
              ARRAY::SetPartOfIndexedArrays(&multiVal[0],&multiVal[0]+multiVal.size(),arrIn,arrIndxIn,srcArr,srcArrIndex,arrOut,arrIndexOut);
              break;
            }
          case 4:
            {
              ARRAY::SetPartOfIndexedArrays(daIntTyypp->begin(),daIntTyypp->end(),arrIn,arrIndxIn,srcArr,srcArrIndex,arrOut,arrIndexOut);
              break;
            }
          default:
            throw INTERP_KERNEL::Exception("ARRAY::SetPartOfIndexedArrays : unrecognized type entered, expected list of int, tuple of int or ARRAY !");
          }
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(arrOut),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(arrIndexOut),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }
      
      static PyObject *FromVTKInternalReprOfPolyedra(const ARRAY *arrIn, const DataArrayIdType *arrIndxIn) throw(INTERP_KERNEL::Exception)
      {
        MCAuto<ARRAY> arrOut;
        MCAuto<DataArrayIdType> arrIndexOut;
        ARRAY::FromVTKInternalReprOfPolyedra(arrIn,arrIndxIn,arrOut,arrIndexOut);
        PyObject *ret=PyTuple_New(2);
        PyTuple_SetItem(ret,0,SWIG_NewPointerObj(SWIG_as_voidptr(arrOut.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        PyTuple_SetItem(ret,1,SWIG_NewPointerObj(SWIG_as_voidptr(arrIndexOut.retn()),SWIGTITraits<mcIdType>::TI, SWIG_POINTER_OWN | 0 ));
        return ret;
      }

      static void SetPartOfIndexedArraysSameIdx(PyObject *li, ARRAY *arrIn, const DataArrayIdType *arrIndxIn,
                                                const ARRAY *srcArr, const DataArrayIdType *srcArrIndex) throw(INTERP_KERNEL::Exception)
      {
        mcIdType sw;
        mcIdType singleVal;
        std::vector<mcIdType> multiVal;
        std::pair<mcIdType, std::pair<mcIdType,mcIdType> > slic;
        MEDCoupling::DataArrayIdType *daIntTyypp=0;
        if(!arrIndxIn)
          throw INTERP_KERNEL::Exception("ARRAY::SetPartOfIndexedArraysSameIdx : null pointer as arrIndex !");
        convertIntStarOrSliceLikePyObjToCpp(li,arrIndxIn->getNumberOfTuples()-1,sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            {
              ARRAY::SetPartOfIndexedArraysSameIdx(&singleVal,&singleVal+1,arrIn,arrIndxIn,srcArr,srcArrIndex);
              break;
            }
          case 2:
            {
              ARRAY::SetPartOfIndexedArraysSameIdx(&multiVal[0],&multiVal[0]+multiVal.size(),arrIn,arrIndxIn,srcArr,srcArrIndex);
              break;
            }
          case 4:
            {
              ARRAY::SetPartOfIndexedArraysSameIdx(daIntTyypp->begin(),daIntTyypp->end(),arrIn,arrIndxIn,srcArr,srcArrIndex);
              break;
            }
          default:
            throw INTERP_KERNEL::Exception("ARRAY::SetPartOfIndexedArraysSameIdx : unrecognized type entered, expected list of int, tuple of int or ARRAY !");
          }
      }

    } // end extent
  };

  class ARRAY ## Tuple;

  class ARRAY ## Iterator
  {
  public:
    ARRAY ## Iterator(ARRAY *da);
    ~ARRAY ## Iterator();
    %extend
    {
      PyObject *next()
      {
        ARRAY ## Tuple *ret=self->nextt();
        if(ret)
          return SWIG_NewPointerObj(SWIG_as_voidptr(ret),SWIGTYPE_p_MEDCoupling__ ## ARRAY ## Tuple,SWIG_POINTER_OWN | 0);
        else
          {
            PyErr_SetString(PyExc_StopIteration,"No more data.");
            return 0;
          }
      }
    }
  };

  class ARRAY ## Tuple
  {
  public:
    std::size_t getNumberOfCompo() const;
    ARRAY *buildDAInt(INT nbOfTuples, INT nbOfCompo) const;
    %extend
    {
      std::string __str__() const
      {
        return self->repr();
      }

      INT __int__() const
      {
        return self->intValue();
      }

      ARRAY *buildDAInt()
      {
        return self->buildDAInt(1,self->getNumberOfCompo());
      }

      PyObject *___iadd___(PyObject *trueSelf, PyObject *obj)
      {
        MCAuto<ARRAY> ret=self->buildDAInt(1,self->getNumberOfCompo());
        MEDCoupling_ ## ARRAY ## ____iadd___(ret,0,obj);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
  
      PyObject *___isub___(PyObject *trueSelf, PyObject *obj)
      {
        MCAuto<ARRAY> ret=self->buildDAInt(1,self->getNumberOfCompo());
        MEDCoupling_ ## ARRAY ## ____isub___(ret,0,obj);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
  
      PyObject *___imul___(PyObject *trueSelf, PyObject *obj)
      {
        MCAuto<ARRAY> ret=self->buildDAInt(1,self->getNumberOfCompo());
        MEDCoupling_ ## ARRAY ## ____imul___(ret,0,obj);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }
      PyObject *___idiv___(PyObject *trueSelf, PyObject *obj)
      {
        MCAuto<ARRAY> ret=self->buildDAInt(1,self->getNumberOfCompo());
        MEDCoupling_ ## ARRAY ## ____idiv___(ret,0,obj);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }

      PyObject *___imod___(PyObject *trueSelf, PyObject *obj)
      {
        MCAuto<ARRAY> ret=self->buildDAInt(1,self->getNumberOfCompo());
        MEDCoupling_ ## ARRAY ## ____imod___(ret,0,obj);
        Py_XINCREF(trueSelf);
        return trueSelf;
      }

      PyObject *__len__()
      {
        return PyInt_FromLong(self->getNumberOfCompo());
      }
  
      PyObject *__getitem__(PyObject *obj)
      {
        const char msg2[]="ARRAY ## Tuple::__getitem__ : Mismatch of slice values in 2nd parameter (components) !";
        mcIdType sw;
        INT singleVal;
        std::vector<INT> multiVal;
        std::pair<mcIdType, std::pair<mcIdType,mcIdType> > slic;
        MEDCoupling::DataArrayIdType *daIntTyypp=0;
        const INT *pt=self->getConstPointer();
        INT nbc=(INT)self->getNumberOfCompo();
        convertIntStarOrSliceLikePyObjToCppWithNegIntInterp(obj,ToIdType(nbc),sw,singleVal,multiVal,slic,daIntTyypp);
        switch(sw)
          {
          case 1:
            {
              if(singleVal>=(INT)nbc)
                {
                  std::ostringstream oss;
                  oss << "Requesting for id " << singleVal << " having only " << nbc << " components !";
                  PyErr_SetString(PyExc_StopIteration,oss.str().c_str());
                  return 0;
                }
              if(singleVal>=0)
                return PyInt_FromLong(pt[singleVal]);
              else
                {
                  if(nbc+singleVal>0)
                    return PyInt_FromLong(pt[nbc+singleVal]);
                  else
                    {
                      std::ostringstream oss;
                      oss << "Requesting for id " << singleVal << " having only " << nbc << " components !";
                      throw INTERP_KERNEL::Exception(oss.str().c_str());
                    }
                }
            }
          case 2:
            {
              PyObject *t=PyTuple_New(multiVal.size());
              for(std::size_t j=0;j<multiVal.size();j++)
                {
                  INT cid=multiVal[j];
                  if(cid>=(INT)nbc)
                    {
                      std::ostringstream oss;
                      oss << "Requesting for id #" << cid << " having only " << nbc << " components !";
                      throw INTERP_KERNEL::Exception(oss.str().c_str());
                    }
                  PyTuple_SetItem(t,j,PyInt_FromLong(pt[cid]));
                }
              return t;
            }
          case 3:
            {
              mcIdType sz=DataArray::GetNumberOfItemGivenBES(slic.first,slic.second.first,slic.second.second,msg2);
              PyObject *t=PyTuple_New(sz);
              for(INT j=0;j<sz;j++)
                PyTuple_SetItem(t,j,PyInt_FromLong(pt[slic.first+j*slic.second.second]));
              return t;
            }
          default:
            throw INTERP_KERNEL::Exception("ARRAY ## Tuple::__getitem__ : unrecognized type entered !");
          }
      }

      ARRAY ## Tuple *__setitem__(PyObject *obj, PyObject *value)
      {
        const char msg[]="DataArrayIntTuple::__setitem__ : unrecognized type entered, int, slice, list<int>, tuple<int> !";
        const char msg2[]="DataArrayIntTuple::__setitem__ : Mismatch of slice values in 2nd parameter (components) !";
        mcIdType sw1,sw2;
        mcIdType singleValV;
        std::vector<mcIdType> multiValV;
        std::pair<mcIdType, std::pair<mcIdType,mcIdType> > slicV;
        MEDCoupling::ARRAY ## Tuple *daIntTyyppV=0;
        mcIdType nbc=ToIdType(self->getNumberOfCompo());
        convertObjToPossibleCpp22<INT>(value,nbc,sw1,singleValV,multiValV,slicV,daIntTyyppV);
        INT singleVal;
        std::vector<INT> multiVal;
        std::pair<mcIdType, std::pair<mcIdType,mcIdType> > slic;
        MEDCoupling::ARRAY *daIntTyypp=0;
        INT *pt=self->getPointer();
        convertIntStarOrSliceLikePyObjToCppWithNegIntInterp(obj,nbc,sw2,singleVal,multiVal,slic,daIntTyypp);
        switch(sw2)
          {
          case 1:
            {
              if(singleVal>=nbc)
                {
                  std::ostringstream oss;
                  oss << "Requesting for setting id # " << singleVal << " having only " << nbc << " components !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
              switch(sw1)
                {
                case 1:
                  {
                    pt[singleVal]=(INT)singleValV;
                    return self;
                  }
                case 2:
                  {
                    if(multiValV.size()!=1)
                      {
                        std::ostringstream oss;
                        oss << "Requesting for setting id # " << singleVal << " with a list or tuple with size != 1 ! ";
                        throw INTERP_KERNEL::Exception(oss.str().c_str());
                      }
                    pt[singleVal]=(INT)multiValV[0];
                    return self;
                  }
                case 4:
                  {
                    pt[singleVal]=daIntTyyppV->getConstPointer()[0];
                    return self;
                  }
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
            }
          case 2:
            {
              switch(sw1)
                {
                case 1:
                  {
                    for(std::vector<INT>::const_iterator it=multiVal.begin();it!=multiVal.end();it++)
                      {
                        if(*it>=nbc)
                          {
                            std::ostringstream oss;
                            oss << "Requesting for setting id # " << *it << " having only " << nbc << " components !";
                            throw INTERP_KERNEL::Exception(oss.str().c_str());
                          }
                        pt[*it]=(INT)singleValV;
                      }
                    return self;
                  }
                case 2:
                  {
                    if(multiVal.size()!=multiValV.size())
                      {
                        std::ostringstream oss;
                        oss << "Mismatch length of during assignment : " << multiValV.size() << " != " << multiVal.size() << " !";
                        throw INTERP_KERNEL::Exception(oss.str().c_str());
                      }
                    for(INT i=0;i<(INT)multiVal.size();i++)
                      {
                        INT pos=multiVal[i];
                        if(pos>=nbc)
                          {
                            std::ostringstream oss;
                            oss << "Requesting for setting id # " << pos << " having only " << nbc << " components !";
                            throw INTERP_KERNEL::Exception(oss.str().c_str());
                          }
                        pt[multiVal[i]]=(INT)multiValV[i];
                      }
                    return self;
                  }
                case 4:
                  {
                    const INT *ptV=daIntTyyppV->getConstPointer();
                    if(nbc>(INT)daIntTyyppV->getNumberOfCompo())
                      {
                        std::ostringstream oss;
                        oss << "Mismatch length of during assignment : " << nbc << " != " << daIntTyyppV->getNumberOfCompo() << " !";
                        throw INTERP_KERNEL::Exception(oss.str().c_str());
                      }
                    std::copy(ptV,ptV+nbc,pt);
                    return self;
                  }
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
            }
          case 3:
            {
              std::size_t sz=DataArray::GetNumberOfItemGivenBES(slic.first,slic.second.first,slic.second.second,msg2);
              switch(sw1)
                {
                case 1:
                  {
                    for(std::size_t j=0;j<sz;j++)
                      pt[slic.first+j*slic.second.second]=(INT)singleValV;
                    return self;
                  }
                case 2:
                  {
                    if(sz!=multiValV.size())
                      {
                        std::ostringstream oss;
                        oss << "Mismatch length of during assignment : " << multiValV.size() << " != " << sz << " !";
                        throw INTERP_KERNEL::Exception(oss.str().c_str());
                      }
                    for(std::size_t j=0;j<sz;j++)
                      pt[slic.first+j*slic.second.second]=(INT)multiValV[j];
                    return self;
                  }
                case 4:
                  {
                    const INT *ptV=daIntTyyppV->getConstPointer();
                    if(sz>daIntTyyppV->getNumberOfCompo())
                      {
                        std::ostringstream oss;
                        oss << "Mismatch length of during assignment : " << nbc << " != " << daIntTyyppV->getNumberOfCompo() << " !";
                        throw INTERP_KERNEL::Exception(oss.str().c_str());
                      }
                    for(std::size_t j=0;j<sz;j++)
                      pt[slic.first+j*slic.second.second]=ptV[j];
                    return self;
                  }
                default:
                  throw INTERP_KERNEL::Exception(msg);
                }
            }
          default:
            throw INTERP_KERNEL::Exception(msg);
          }
      }
    }
//   };
// }
%enddef


%define TRANSFORMWITHINDARR( ARRAY, INT )
void transformWithIndArr(PyObject *li)
{
  void *da=0;
  int res1(SWIG_ConvertPtr(li,&da,SWIGTITraits<INT>::TI, 0 |  0 ));
  if (!SWIG_IsOK(res1))
    {
      int res2(SWIG_ConvertPtr(li,&da,SWIGTYPE_p_MEDCoupling__MapII, 0 |  0 ));
      if(SWIG_IsOK(res2))
        {
          MapII *m=reinterpret_cast<MapII *>(da);
          self->transformWithIndArr(*m);
        }
      else
        {
          mcIdType size;
          INTERP_KERNEL::AutoPtr<INT> tmp=convertPyToNewIntArr2<INT>(li,&size);
          self->transformWithIndArr(tmp,tmp+size);
        }
    }
  else
    {
      ARRAY *da2=reinterpret_cast< ARRAY * >(da);
      self->transformWithIndArr(da2->getConstPointer(),da2->getConstPointer()+da2->getNbOfElems());
    }
}
%enddef

namespace MEDCoupling
{
  typedef int      Int32;
#ifdef WIN32
  typedef long long Int64;
#else
  typedef long int  Int64;
#endif

  class DataArrayInt32Iterator;

  class DataArrayInt32 : public DataArray
  {
    ARRAYDEF( DataArrayInt32, Int32 )
  };
  %extend DataArrayInt32 {
    DataArrayInt64 *convertToInt64Arr() const
    {
      MCAuto<DataArrayInt64> ret(self->convertToInt64Arr());
      return ret.retn();
    }
#ifdef WITH_NUMPY
    PyObject *toNumPyArray() // not const. It is not a bug !
    {
      return ToNumPyArray<DataArrayInt32,MEDCoupling::Int32>(self,NPY_INT32,"DataArrayInt32");
    }
#endif
#ifndef MEDCOUPLING_USE_64BIT_IDS
    MCAuto< MapII > invertArrayN2O2O2NOptimized()
    {
      return self->invertArrayN2O2O2NOptimized();
    }
    MCAuto< MapII > giveN2OOptimized()
    {
      return self->giveN2OOptimized();
    }
    TRANSFORMWITHINDARR( DataArrayInt32, Int32 )
#endif
  };
  
  class DataArrayInt64Iterator;

  class DataArrayInt64 : public DataArray
  {
    ARRAYDEF( DataArrayInt64, Int64 )
  };
  %extend DataArrayInt64 {
    DataArrayInt32 *convertToInt32Arr() const
    {
      MCAuto<DataArrayInt32> ret(self->convertToInt32Arr());
      return ret.retn();
    }
#ifdef WITH_NUMPY
    PyObject *toNumPyArray() // not const. It is not a bug !
    {
      return ToNumPyArray<DataArrayInt64,MEDCoupling::Int64>(self,NPY_INT64,"DataArrayInt64");
    }
#endif
#ifdef MEDCOUPLING_USE_64BIT_IDS
    MCAuto< MapII > invertArrayN2O2O2NOptimized()
    {
      return self->invertArrayN2O2O2NOptimized();
    }
    MCAuto< MapII > giveN2OOptimized()
    {
      return self->giveN2OOptimized();
    }
    TRANSFORMWITHINDARR( DataArrayInt64, Int64 )
#endif
  };

}
