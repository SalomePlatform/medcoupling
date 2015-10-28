// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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

#include "CommInterface.hxx"

namespace ParaMEDMEM
{
  /*! \anchor CommInterface-det
     \class CommInterface

    The class \a CommInterface is the gateway to the MPI library.

    It is a helper class that gathers the calls to the MPI
    library that are made in the %ParaMEDMEM library. This gathering
    allows easier gathering of information about the communication
    in the library.

    It is typically called after the MPI_Init() call in a program. It is afterwards passed as a parameter to the constructors of %ParaMEDMEM objects so that they access the MPI library via the CommInterface.

    As an example, the following code excerpt initializes a processor group made of the zero processor.

    \verbatim
    #include "CommInterface.hxx"
    #include "ProcessorGroup.hxx"

    int main(int argc, char** argv)
    {
    //initialization
    MPI_Init(&argc, &argv);
    ParaMEDMEM::CommInterface comm_interface;

    //setting up a processor group with proc 0
    set<int> procs;
    procs.insert(0);
    ParaMEDMEM::ProcessorGroup group(procs, comm_interface);

    //cleanup
    MPI_Finalize();
    }
    \endverbatim
  */

  CommInterface::CommInterface()
  {
  }

  CommInterface::~CommInterface()
  {
  }
}
