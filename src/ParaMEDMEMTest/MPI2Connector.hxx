// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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

#ifndef __MPI2CONNECTOR_HXX__
#define __MPI2CONNECTOR_HXX__

#include <mpi.h>
#include <string>
#include <sstream>

class MPI2Connector
{
public:
  MPI2Connector();
  ~MPI2Connector();
  // MPI2 connection
  MPI_Comm remoteMPI2Connect(const std::string& service);
  // MPI2 disconnection
  void remoteMPI2Disconnect(const std::string& service);
private:
  // Processus id
  int _num_proc;
  // Processus size
  int _nb_proc;
  MPI_Comm _gcom;
  bool _srv;
  std::string _port_name;
private:
  static const int TIMEOUT=5;
};

#endif
