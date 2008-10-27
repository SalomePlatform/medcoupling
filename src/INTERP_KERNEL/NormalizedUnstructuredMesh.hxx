#ifndef __NORMALIZEDUNSTRUCTUREDMESH_HXX__
#define __NORMALIZEDUNSTRUCTUREDMESH_HXX__

namespace INTERP_KERNEL
{
  typedef enum
    {
      ALL_C_MODE       ,
      ALL_FORTRAN_MODE
    } NumberingPolicy;


  typedef enum
    {
      NORM_SEG2    =  1,
      NORM_SEG3    =  2,
      NORM_TRI3    =  3,
      NORM_QUAD4   =  4,
      NORM_POLYGON =  5,
      NORM_TRI6    =  6,
      NORM_QUAD8   =  8,
      //
      NORM_TETRA4  = 14,
      NORM_PYRA5   = 15,
      NORM_PENTA6  = 16,
      NORM_HEXA8   = 18,
      NORM_TETRA10 = 20,
      NORM_PYRA13  = 23,
      NORM_PENTA15 = 25,
      NORM_HEXA20  = 30,
      NORM_POLYHED = 31
    } NormalizedCellType;

  class GenericMesh
  {};
}

#endif
