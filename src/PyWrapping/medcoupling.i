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

// Author : Anthony Geay (EDF R&D)

%module medcoupling

#define MEDCOUPLINGREMAPPER_EXPORT
#define INTERPKERNEL_EXPORT
#define MEDCOUPLING_EXPORT
#define MEDLOADER_EXPORT

%include "MEDCouplingCommon.i"

%include "MEDCouplingRemapperCommon.i"

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

%{
  static const char SEQ_INTERPOL_EXT[]="Sequential interpolator";
  static const char MEDFILEIO_EXT[]="MED file I/O";
  static const char RENUM_EXT[]="Renumberer";
  static const char PART_EXT[]="Partitioner";
  static const char PAR_INTERPOL_EXT[]="Parallel interpolator (SPMD paradigm)";
  
  static const char *EXTENSIONS[]={SEQ_INTERPOL_EXT,MEDFILEIO_EXT,RENUM_EXT,PART_EXT,PAR_INTERPOL_EXT};
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
def MEDCouplingDataArrayIntIadd(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt____iadd___(self, self, *args)
def MEDCouplingDataArrayIntIsub(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt____isub___(self, self, *args)
def MEDCouplingDataArrayIntImul(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt____imul___(self, self, *args)
def MEDCouplingDataArrayIntIdiv(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt____idiv___(self, self, *args)
def MEDCouplingDataArrayIntImod(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt____imod___(self, self, *args)
def MEDCouplingDataArrayIntIpow(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayInt____ipow___(self, self, *args)
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
def MEDCouplingDataArrayIntTupleIadd(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayIntTuple____iadd___(self, self, *args)
def MEDCouplingDataArrayIntTupleIsub(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayIntTuple____isub___(self, self, *args)
def MEDCouplingDataArrayIntTupleImul(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayIntTuple____imul___(self, self, *args)
def MEDCouplingDataArrayIntTupleIdiv(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayIntTuple____idiv___(self, self, *args)
def MEDCouplingDataArrayIntTupleImod(self,*args):
    import _medcoupling
    return _medcoupling.DataArrayIntTuple____imod___(self, self, *args)
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

