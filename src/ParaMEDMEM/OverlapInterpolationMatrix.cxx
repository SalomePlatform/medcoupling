//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

#include "OverlapInterpolationMatrix.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "ProcessorGroup.hxx"
#include "MxN_Mapping.hxx"
#include "TranslationRotationMatrix.hxx"
#include "Interpolation.hxx"
#include "Interpolation1D.txx"
#include "Interpolation2DCurve.hxx"
#include "Interpolation2D.txx"
#include "Interpolation3DSurf.hxx"
#include "Interpolation3D.txx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingNormalizedUnstructuredMesh.txx"
#include "InterpolationOptions.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "ElementLocator.hxx"

#include <algorithm>

using namespace std;

namespace ParaMEDMEM
{
  OverlapInterpolationMatrix::OverlapInterpolationMatrix(const ParaFIELD *source_field,
                                                         const ParaFIELD *target_field,
                                                         const ProcessorGroup& group,
                                                         const DECOptions& dec_options,
                                                         const INTERP_KERNEL::InterpolationOptions& i_opt):
    INTERP_KERNEL::InterpolationOptions(i_opt),
    DECOptions(dec_options),
    _source_field(source_field),
    _target_field(target_field),
    _source_support(source_field->getSupport()->getCellMesh()),
    _target_support(target_field->getSupport()->getCellMesh()),
    //_mapping(source_group, target_group, dec_options),
    _group(group)
  {
    int nbelems = source_field->getField()->getNumberOfTuples();
    _row_offsets.resize(nbelems+1);
    _coeffs.resize(nbelems);
    _target_volume.resize(nbelems);
  }

  OverlapInterpolationMatrix::~OverlapInterpolationMatrix()
  {
  }
}
