#include "CommInterface.hxx"

namespace ParaMEDMEM
{
/*! \defgroup comm_interface CommInterface
	Class \a CommInterface is the gateway to the MPI library.
	It is a helper class that gathers the calls to the MPI
	library that are made in the ParaMEDMEM library. This gathering
	allows easier gathering of information about the communication
	in the library.

It is typically called after the MPI_Init() call in a program. It is afterwards passed as a parameter to the constructors of ParaMEDMEM objects so that they access the MPI library via the CommInterface.

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
