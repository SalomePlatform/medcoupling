//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D, OPEN CASCADE
//
//  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#ifndef REMAPPER_HXX_
#define REMAPPER_HXX_

#include "InterpKernelMatrix.hxx"
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Support.hxx"
#include "MEDMEM_Field.hxx"

namespace INTERP_KERNEL
{

  class INTERPKERNEL_EXPORT Remapper
  {
  public:
    Remapper();
    virtual ~Remapper();
    void prepare(const MEDMEM::MESH& mesh_source, const MEDMEM::MESH& mesh_target);
    void transfer(const MEDMEM::FIELD<double>& field_source, MEDMEM::FIELD<double>& field_target);
    void setOptionDouble(const std::string& key, double value);
    void setOptionInt(const std::string& key, int value);
  private :
    Matrix<double,ALL_FORTRAN_MODE>* _matrix;
    MEDMEM::FIELD<double>* getSupportVolumes(const MEDMEM::SUPPORT& support);
  };

}

#endif /*REMAPPER_HXX_*/
