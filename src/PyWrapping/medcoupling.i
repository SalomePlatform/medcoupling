// Copyright (C) 2017-2024  CEA, EDF
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

%module medcoupling

#define MEDCOUPLINGREMAPPER_EXPORT
#define INTERPKERNEL_EXPORT
#define MEDCOUPLING_EXPORT
#define MEDLOADER_EXPORT

%include "MEDCouplingCommon.i"

%include "MEDCouplingRemapperCommon.i"

%include "ICoCoMEDField.i"

#ifdef WITH_MED_FILE
%include "MEDLoaderCommon.i"
#endif

#ifdef WITH_RENUMBER
%include "MEDRenumberImpl.i"
#endif

#ifdef WITH_PARTITIONER
%include "MEDPartitionerCommon.i"
#endif

#ifdef WITH_PARALLEL_INTERPOLATOR
%include "ParaMEDMEMCommon.i"
#endif

%constant const char __version__[]=MEDCOUPLING_GIT_SHA1;
%constant const char __config_datetime__[]=MEDCOUPLING_CONFIG_DT;

%{
  static const char SEQ_INTERPOL_EXT[]="Sequential interpolator";
  static const char MEDFILEIO_EXT[]="MED file I/O";
  static const char RENUM_EXT[]="Renumberer";
  static const char PART_EXT[]="Partitioner";
  static const char PAR_INTERPOL_EXT[]="Parallel interpolator (SPMD paradigm)";
  static const char IT_STATS_EXT[] = "Iterative statistics";
  
  static const char *EXTENSIONS[]={SEQ_INTERPOL_EXT,MEDFILEIO_EXT,RENUM_EXT,PART_EXT,PAR_INTERPOL_EXT,IT_STATS_EXT};
  static const int NB_OF_EXTENSIONS=sizeof(EXTENSIONS)/sizeof(const char *);
%}

%inline
{
  std::vector<std::string> AllPossibleExtensions()
  {
    std::vector<std::string> ret(EXTENSIONS,EXTENSIONS+NB_OF_EXTENSIONS);
    return ret;
  }

  bool HasMEDFileExt()
  {
#ifdef WITH_MED_FILE
    return true;
#else
    return false;
#endif
  }

  bool HasRenumberExt()
  {
#ifdef WITH_RENUMBER
    return true;
#else
    return false;
#endif
  }

  bool HasPartitionerExt()
  {
#ifdef WITH_PARTITIONER
    return true;
#else
    return false;
#endif
  }

  bool HasScotchPartitionerAlg()
  {
#ifdef WITH_PARTITIONER
    return MEDPartitioner::HasScotchAlg();
#else
    return false;
#endif    
  }

  bool HasPTScotchPartitionerAlg()
  {
#ifdef WITH_PARTITIONER
    return MEDPartitioner::HasPTScotchAlg();
#else
    return false;
#endif    
  }

  bool HasMetisPartitionerAlg()
  {
#ifdef WITH_PARTITIONER
    return MEDPartitioner::HasMetisAlg();
#else
    return false;
#endif    
  }
  
  bool HasParallelInterpolatorExt()
  {
#ifdef WITH_PARALLEL_INTERPOLATOR
    return true;
#else
    return false;
#endif
  }
  
  std::vector<std::string> ActiveExtensions()
  {
    std::vector<std::string> ret;
    ret.push_back(std::string(SEQ_INTERPOL_EXT));
#ifdef WITH_MED_FILE
    ret.push_back(std::string(MEDFILEIO_EXT));
#endif
#ifdef WITH_RENUMBER
    ret.push_back(std::string(RENUM_EXT));
#endif
#ifdef WITH_PARTITIONER
    ret.push_back(std::string(PART_EXT));
#endif
#ifdef WITH_PARALLEL_INTERPOLATOR
    ret.push_back(std::string(PAR_INTERPOL_EXT));
#endif
#ifdef WITH_ITERATIVE_STATISTICS
    ret.push_back(std::string(IT_STATS_EXT));
#endif
    return ret;
  }
}

%pythoncode %{
def MEDCouplingDataArrayDoubleIadd(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayDouble____iadd___(self, self, *args)
def MEDCouplingDataArrayDoubleIsub(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayDouble____isub___(self, self, *args)
def MEDCouplingDataArrayDoubleImul(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayDouble____imul___(self, self, *args)
def MEDCouplingDataArrayDoubleIdiv(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayDouble____idiv___(self, self, *args)
def MEDCouplingDataArrayDoubleIpow(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayDouble____ipow___(self, self, *args)
def MEDCouplingFieldDoubleIadd(self,*args):
    import _medcoupling
    return _medcoupling.MEDCouplingFieldDouble____iadd___(self, self, *args)
def MEDCouplingFieldDoubleIsub(self,*args):
    import _medcoupling
    return _medcoupling.MEDCouplingFieldDouble____isub___(self, self, *args)
def MEDCouplingFieldDoubleImul(self,*args):
    import _medcoupling
    return _medcoupling.MEDCouplingFieldDouble____imul___(self, self, *args)
def MEDCouplingFieldDoubleIdiv(self,*args):
    import _medcoupling
    return _medcoupling.MEDCouplingFieldDouble____idiv___(self, self, *args)
def MEDCouplingFieldDoubleIpow(self,*args):
    import _medcoupling
    return _medcoupling.MEDCouplingFieldDouble____ipow___(self, self, *args)
def MEDCouplingDataArrayInt32Iadd(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt32____iadd___(self, self, *args)
def MEDCouplingDataArrayInt32Isub(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt32____isub___(self, self, *args)
def MEDCouplingDataArrayInt32Imul(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt32____imul___(self, self, *args)
def MEDCouplingDataArrayInt32Idiv(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt32____idiv___(self, self, *args)
def MEDCouplingDataArrayInt32Imod(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt32____imod___(self, self, *args)
def MEDCouplingDataArrayInt32Ipow(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt32____ipow___(self, self, *args)
def MEDCouplingDataArrayInt64Iadd(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt64____iadd___(self, self, *args)
def MEDCouplingDataArrayInt64Isub(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt64____isub___(self, self, *args)
def MEDCouplingDataArrayInt64Imul(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt64____imul___(self, self, *args)
def MEDCouplingDataArrayInt64Idiv(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt64____idiv___(self, self, *args)
def MEDCouplingDataArrayInt64Imod(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt64____imod___(self, self, *args)
def MEDCouplingDataArrayInt64Ipow(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt64____ipow___(self, self, *args)
def MEDCouplingDataArrayFloatIadd(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayFloat____iadd___(self, self, *args)
def MEDCouplingDataArrayFloatIsub(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayFloat____isub___(self, self, *args)
def MEDCouplingDataArrayFloatImul(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayFloat____imul___(self, self, *args)
def MEDCouplingDataArrayFloatIdiv(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayFloat____idiv___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIadd(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayDoubleTuple____iadd___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIsub(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayDoubleTuple____isub___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleImul(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayDoubleTuple____imul___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIdiv(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayDoubleTuple____idiv___(self, self, *args)
def MEDCouplingDataArrayInt32TupleIadd(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt32Tuple____iadd___(self, self, *args)
def MEDCouplingDataArrayInt32TupleIsub(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt32Tuple____isub___(self, self, *args)
def MEDCouplingDataArrayInt32TupleImul(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt32Tuple____imul___(self, self, *args)
def MEDCouplingDataArrayInt32TupleIdiv(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt32Tuple____idiv___(self, self, *args)
def MEDCouplingDataArrayInt32TupleImod(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt32Tuple____imod___(self, self, *args)
def MEDCouplingDataArrayInt64TupleIadd(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt64Tuple____iadd___(self, self, *args)
def MEDCouplingDataArrayInt64TupleIsub(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt64Tuple____isub___(self, self, *args)
def MEDCouplingDataArrayInt64TupleImul(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt64Tuple____imul___(self, self, *args)
def MEDCouplingDataArrayInt64TupleIdiv(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt64Tuple____idiv___(self, self, *args)
def MEDCouplingDataArrayInt64TupleImod(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt64Tuple____imod___(self, self, *args)
def MEDCouplingDenseMatrixIadd(self,*args):
    import _medcoupling
    return _medcoupling.DenseMatrix____iadd___(self, self, *args)
def MEDCouplingDenseMatrixIsub(self,*args):
    import _medcoupling
    return _medcoupling.DenseMatrix____isub___(self, self, *args)
%}

%include "MEDCouplingFinalize.i"

#ifdef WITH_MED_FILE
%include "MEDLoaderFinalize.i"
#endif

%include "medcoupling_pycode"

