//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#include "MEDMEM_Exception.hxx"
#include "Topology.hxx"
#include "BlockTopology.hxx"
#include "ComponentTopology.hxx"
#include "ParaSUPPORT.hxx"
#include "StructuredParaSUPPORT.hxx"
#include "UnstructuredParaSUPPORT.hxx"
#include "ExplicitCoincidentDEC.hxx"
#include "StructuredCoincidentDEC.hxx"
#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "ParaFIELD.hxx"
#include "ParaMESH.hxx"

using namespace MEDMEM;


namespace ParaMEDMEM
{



/*!
\defgroup parafield ParaFIELD
This class encapsulates parallel fields. It basically encapsulates
a MEDMEM::FIELD<double> with extra information related to parallel 
topology.
It is most conveniently created by giving a pointer to a MEDMEM::FIELD<double>
object and a \c ProcessorGroup.
By default, a ParaFIELD object will be constructed with all field components
located on the same processors. In some specific cases, it might be necessary to scatter components over several processors. In this case, the constructor
using a ComponentTopology is required.

@{ */

/*!

\brief  Constructing a \c ParaFIELD from a \c ParaSUPPORT and a \c ComponentTopology.

This constructor creates an empty field based on the ParaSUPPORT description 
and the partitioning of components described in \a component_topology.
It takes ownership over the \c _field object that it creates.

Here come the three ComponentTopology constructors :
\verbatim
ComponentTopology c; // one component in the field
ComponentTopology c(6); //six components, all of them on the same processor
ComponentTopology c(6, proc_group); // six components, evenly distributed over the processors of procgroup
\endverbatim

*/
ParaFIELD::ParaFIELD(const ParaSUPPORT* para_support, const ComponentTopology& component_topology)
:_support(para_support),
 _component_topology(component_topology),
 _field(0),
 _has_field_ownership(true),
 _has_support_ownership(false)
{
	if (dynamic_cast<const StructuredParaSUPPORT*>(para_support)!=0
			|| (para_support->getTopology()->getProcGroup()->size()==1 && component_topology.nbBlocks()!=1))
	{const BlockTopology* source_topo = dynamic_cast<const BlockTopology*>(para_support->getTopology());
		_topology=new BlockTopology(*source_topo,component_topology);
	}
	else
	{
	  if (component_topology.nbBlocks()!=1 &&  para_support->getTopology()->getProcGroup()->size()!=1)
	    throw MEDEXCEPTION(LOCALIZED(
		"ParaFIELD constructor : Unstructured Support not taken into account with component topology yet"));
	  else 
	    {
	      const BlockTopology* source_topo=
		dynamic_cast<const BlockTopology*> (para_support->getTopology());
		int nb_local_comp=component_topology.nbLocalComponents();
	      _topology=new BlockTopology(*source_topo,nb_local_comp);
					     
	    }
	}

	int nb_components = component_topology.nbLocalComponents();
	if (nb_components!=0)
		{
			_field=new FIELD<double> (para_support->getSupport(), nb_components);
		}
	else return;
	
	_field->setName("Default ParaFIELD name");
	_field->setDescription("Default ParaFIELD description");
	_field->setNumberOfValues(para_support->getSupport()->getNumberOfElements(MED_EN::MED_ALL_ELEMENTS));
	string* compnames=new string[nb_components];
	string* compdesc=new string[nb_components];
	string* compunit=new string[nb_components];
	for (int i=0; i<nb_components; i++)
	{
		ostringstream stream(compnames[i]);
		ostringstream stream2(compdesc[i]);
		ostringstream stream3(compunit[i]);
		stream<<"component "<<i;
		stream2<<"component description "<<i;
				stream3<<"component unit "<<i;

		compnames[i]=stream.str();
		compdesc[i]=stream2.str();
				compunit[i]=stream3.str();

	}
	_field->setComponentsNames(compnames);
	_field->setComponentsDescriptions(compdesc);
	_field->setMEDComponentsUnits(compunit);
	_field->setIterationNumber(0);
	_field->setOrderNumber(0);
	_field->setTime(0.0);
	delete[] compnames;
	delete[] compdesc;
  delete[] compunit;
} 

/*! \brief Constructor creating the ParaFIELD
 from a given FIELD and a processor group. 

This constructor supposes that support underlying \a subdomain_field has no ParaSUPPORT 
attached and it therefore recreates one. It therefore takes ownership over _support. The component topology associated with the field is a basic one (all components on the same processor). 
 */

ParaFIELD::ParaFIELD(MEDMEM::FIELD<double>* subdomain_field, const ProcessorGroup& proc_group):
_field(subdomain_field),
_support(new UnstructuredParaSUPPORT(subdomain_field->getSupport(), proc_group)),
_component_topology(ComponentTopology(_field->getNumberOfComponents())),
_has_field_ownership(false),
_has_support_ownership(true)
{
  const ExplicitTopology* source_topo=
    dynamic_cast<const ExplicitTopology*> (_support->getTopology());
    _topology=new ExplicitTopology(*source_topo,_component_topology.nbLocalComponents());
}

ParaFIELD::ParaFIELD(MEDMEM::driverTypes driver_type, const string& file_name, 
	const string& driver_name, const ComponentTopology& component_topology) 
	throw (MEDEXCEPTION):_component_topology(component_topology),_field(0){}

ParaFIELD::~ParaFIELD()
{
  if (_topology!=0)
    delete _topology;
	if (_has_support_ownership)
		{
			delete _support;
			_support=0;
		}

	if (_has_field_ownership)
		{
			delete _field; _field=0;
		}
			
}

	/*! writes the contents of the field in MED v2.3 file \a filename :
- \c filename is the name of the master file
-  files containing the subdomains are named \c filename1.med, \c filename2.med, etc...
	*/
void ParaFIELD::write(MEDMEM::driverTypes driverType, const string& fileName, const string& meshName){
  //	Topology* topo = dynamic_cast<BlockTopology*> (_topology);
	int myrank = _topology->getProcGroup()->myRank();
	ostringstream name;
	name <<fileName<<myrank+1<<".med";
	cout << name <<endl;
	int driver = _field->addDriver(driverType, name.str().c_str(), meshName);
	_field->write(driver);
}

void ParaFIELD::synchronizeTarget(ParaFIELD* source_field){
	DEC* data_channel;
	if (dynamic_cast<BlockTopology*>(_topology)!=0)
	{
		data_channel=new StructuredCoincidentDEC();
	}
	else
	{
		data_channel=new ExplicitCoincidentDEC();
	}
	data_channel->attachLocalField(this);
	data_channel->synchronize();
	data_channel->prepareTargetDE();
	data_channel->recvData();
	
	delete data_channel;
}

void ParaFIELD::synchronizeSource(ParaFIELD* target_field){
	DEC* data_channel;
	if (dynamic_cast<BlockTopology*>(_topology)!=0)
	{
		data_channel=new StructuredCoincidentDEC();
	}
	else
	{
		data_channel=new ExplicitCoincidentDEC();
	}
	data_channel->attachLocalField(this);
	data_channel->synchronize();
	data_channel->prepareSourceDE();
	data_channel->sendData();
	
	delete data_channel;
}

	/*! This method retrieves the integral of component \a icomp
		over the all domain. */
double ParaFIELD::getVolumeIntegral(int icomp) const
{
	CommInterface comm_interface = _topology->getProcGroup()->getCommInterface();
	const MEDMEM::SUPPORT* support = _field->getSupport();
	FIELD<double>* volume= getSupportVolumes(*support);
	double integral=0.0;
	for (int i=0; i<support->getNumberOfElements(MED_EN::MED_ALL_ELEMENTS); i++)
		integral+=_field->getValueIJ(i+1,icomp)*volume->getValueIJ(i+1,1);
	delete volume;

	double total=0.0;
	const MPI_Comm* comm = (dynamic_cast<const MPIProcessorGroup*>(_topology->getProcGroup()))->getComm();
	comm_interface.allReduce(&integral, &total, 1, MPI_DOUBLE, MPI_SUM, *comm);
	
	return total;
}

/*!
\brief returns the volumes of the cells underlying the field \a field

For 2D geometries, the returned field contains the areas.
For 3D geometries, the returned field contains the volumes.

\param field field on which cells the volumes are required
\return field containing the volumes
*/
 MEDMEM::FIELD<double>* ParaFIELD::getSupportVolumes(const MEDMEM::SUPPORT& support) const
	 {
		 const MEDMEM::MESH* mesh=support.getMesh();
		 int dim = mesh->getMeshDimension();
		 switch (dim)
			 {
			 case 2:
				 return mesh->getArea(&support);
			 case 3:
				 return mesh->getVolume(&support);
			 default:
				 throw MEDMEM::MEDEXCEPTION("interpolation is not available for this dimension");
			 }
	 }

/*! @} */
}

