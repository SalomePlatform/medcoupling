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

#include "CommInterface.hxx"

namespace MEDCoupling
{
  /*! \anchor CommInterface-det
     \class CommInterface

    The class \a CommInterface is the gateway to the MPI library.
    It is a wrapper around all MPI calls, thus trying to abstract the rest of the code from using the direct MPI API
    (but this is not strictly respected overall in practice ...). It is used in all
    the \ref parallel "DEC related classes".

    It is typically instanciated after the MPI_Init() call in a program and is afterwards passed as a
    parameter to the constructors of various \ref parallel "parallel objects" so that they access the
    MPI library via this common interface.

    As an example, the following code excerpt initializes a processor group made of the zero processor.

    \verbatim
    #include "CommInterface.hxx"
    #include "ProcessorGroup.hxx"

    int main(int argc, char** argv)
    {
    //initialization
    MPI_Init(&argc, &argv);
    MEDCoupling::CommInterface comm_interface;

    //setting up a processor group with proc 0
    set<int> procs;
    procs.insert(0);
    MEDCoupling::ProcessorGroup group(procs, comm_interface);

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
