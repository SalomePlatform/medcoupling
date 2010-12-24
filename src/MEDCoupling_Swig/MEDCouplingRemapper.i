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
%include "MEDCouplingRemapper.hxx"
