// Copyright (C) 2016  CEA/DEN, EDF R&D
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

#ifndef __MEDLOADERNS_HXX__
#define __MEDLOADERNS_HXX__

#include <vector>
#include <string>

namespace MEDCoupling
{
  class MEDCouplingFieldDouble;
}

#include "med.h"

namespace MEDLoaderNS
{
  int readUMeshDimFromFile(const std::string& fileName, const std::string& meshName, std::vector<int>& possibilities);
  void dispatchElems(int nbOfElemCell, int nbOfElemFace, int& nbOfElem, med_entity_type& whichEntity);
  template<class T>
  void writeFieldWithoutReadingAndMappingOfMeshInFile(const std::string& fileName, const typename MEDCoupling::Traits<T>::FieldType *f, bool writeFromScratch);
  med_int getIdFromMeshName(med_idt fid, const std::string& meshName, std::string& trueMeshName);
  std::vector<std::string> getMeshNamesFid(med_idt fid);
}

#endif
