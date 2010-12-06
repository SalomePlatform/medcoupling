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

#ifndef __MEDFILEUTILITIES_HXX__
#define __MEDFILEUTILITIES_HXX__

#include "InterpKernelException.hxx"

extern "C"
{
#include "med.h"
}

namespace MEDFileUtilities
{
  med_mode_acces TraduceWriteMode(int medloaderwritemode) throw(INTERP_KERNEL::Exception);
  void CheckMEDCode(int code, med_idt fid, const char *msg) throw(INTERP_KERNEL::Exception);
  void CheckFileForRead(const char *fileName) throw(INTERP_KERNEL::Exception);
}
  
#endif
