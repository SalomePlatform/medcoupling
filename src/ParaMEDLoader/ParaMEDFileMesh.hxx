// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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

#ifndef __PARAMEDFILEMESH_HXX__
#define __PARAMEDFILEMESH_HXX__

#include "mpi.h"

#include <string>

namespace ParaMEDMEM
{
  class MEDFileMesh;
  class MEDFileUMesh;
  class MEDFileMeshes;
  class MEDFileMeshReadSelector;

  class ParaMEDFileMesh
  {
  public:
    static MEDFileMesh *New(const MPI_Comm com, const MPI_Info nfo, const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
  };

  class ParaMEDFileUMesh
  {
  public:
    static MEDFileUMesh *New(const MPI_Comm com, const MPI_Info nfo, const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
  };

  class ParaMEDFileMeshes
  {
  public:
    static MEDFileMeshes *New(const MPI_Comm com, const MPI_Info nfo, const std::string& fileName);
  };
}

#endif
