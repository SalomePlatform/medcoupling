// Copyright (C) 2007-2017  CEA/DEN, EDF R&D
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

#ifndef __MEDFILEFIELDVISITOR_HXX__
#define __MEDFILEFIELDVISITOR_HXX__

#include "MEDLoaderDefines.hxx"
#include "MEDFileField.hxx"

namespace MEDCoupling
{
  class MEDFileFieldPerMesh;
  class MEDFileAnyTypeField1TSWithoutSDA;
  class MEDFileAnyTypeFieldMultiTSWithoutSDA;
  class MEDFileFieldPerMeshPerTypeCommon;
  class MEDFileFieldPerMeshPerTypePerDisc;
  
  class MEDFileFieldVisitor
  {
  public:
    virtual void newFieldEntry(const MEDFileAnyTypeFieldMultiTSWithoutSDA *field) = 0;
    virtual void endFieldEntry(const MEDFileAnyTypeFieldMultiTSWithoutSDA *field) = 0;
    //
    virtual void newTimeStepEntry(const MEDFileAnyTypeField1TSWithoutSDA *ts) = 0;
    virtual void endTimeStepEntry(const MEDFileAnyTypeField1TSWithoutSDA *ts) = 0;
    //
    virtual void newMeshEntry(const MEDFileFieldPerMesh *fpm) = 0;
    virtual void endMeshEntry(const MEDFileFieldPerMesh *fpm) = 0;
    //
    virtual void newPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt) = 0;
    virtual void endPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt) = 0;
    //
    virtual void newPerMeshPerTypePerDisc(const MEDFileFieldPerMeshPerTypePerDisc *pmptpd) = 0;
  };
}

#endif
