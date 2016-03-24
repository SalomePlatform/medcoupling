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

#include "Topology.hxx"
#include "BlockTopology.hxx"
#include "ComponentTopology.hxx"
#include "ExplicitCoincidentDEC.hxx"
#include "StructuredCoincidentDEC.hxx"
#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "ParaFIELD.hxx"
#include "ParaMESH.hxx"
#include "InterpKernelUtilities.hxx"
#include "InterpolationMatrix.hxx"

#include <numeric>

namespace MEDCoupling
{
  /*!
    \anchor ParaFIELD-det
    \class ParaFIELD

    This class encapsulates parallel fields.

    It gathers a \ref fields "MEDCouplingField" with some extra information related to the parallel
    topology.

    It is most conveniently created by giving a pointer to a MEDCouplingFieldDouble
    object and a ProcessorGroup.
    By default, a ParaFIELD object will be constructed with all field components
    located on the same processors. In some specific cases, it might be necessary to scatter components over
    several processors. In this case, the constructor using a ComponentTopology is required.

    */

  /*!

    \brief  Constructing a \c ParaFIELD from a \c ParaMESH and a \c ComponentTopology.

    This constructor creates an empty field based on the ParaMESH description
    and the partitioning of components described in \a component_topology.
    It takes ownership over the \c _field object that it creates.

    Here come the three ComponentTopology constructors :
    \verbatim
    ComponentTopology c; // one component in the field
    ComponentTopology c(6); //six components, all of them on the same processor
    ComponentTopology c(6, proc_group); // six components, evenly distributed over the processors of procgroup
    \endverbatim

  */
  ParaFIELD::ParaFIELD(TypeOfField type, TypeOfTimeDiscretization td, ParaMESH* para_support, const ComponentTopology& component_topology)
    :_field(0),
     _component_topology(component_topology),_topology(0),_own_support(false),
     _support(para_support)
  {
    if (para_support->isStructured() || (para_support->getTopology()->getProcGroup()->size()==1 && component_topology.nbBlocks()!=1))
      {
        const BlockTopology* source_topo = dynamic_cast<const BlockTopology*>(para_support->getTopology());
        _topology=new BlockTopology(*source_topo,component_topology);
      }
    else
      {
        if (component_topology.nbBlocks()!=1 &&  para_support->getTopology()->getProcGroup()->size()!=1)
          throw INTERP_KERNEL::Exception(LOCALIZED("ParaFIELD constructor : Unstructured Support not taken into account with component topology yet"));
        else 
          {
            const BlockTopology* source_topo=dynamic_cast<const BlockTopology*> (para_support->getTopology());
            int nb_local_comp=component_topology.nbLocalComponents();
            _topology=new BlockTopology(*source_topo,nb_local_comp);
          }
      }
    int nb_components = component_topology.nbLocalComponents();
    if (nb_components!=0)
      {
        _field=MEDCouplingFieldDouble::New(type,td);
        _field->setMesh(_support->getCellMesh());
        DataArrayDouble *array=DataArrayDouble::New();
        array->alloc(_field->getNumberOfTuples(),nb_components);
        _field->setArray(array);
        array->decrRef();
      }
    else return;
  
    _field->setName("Default ParaFIELD name");
    _field->setDescription("Default ParaFIELD description");
  } 

  /*! \brief Constructor creating the ParaFIELD
    from a given FIELD and a processor group. 

    This constructor supposes that support underlying \a subdomain_field has no ParaMESH
    attached and it therefore recreates one. It therefore takes ownership over _support. The component topology associated with the field is a basic one (all components on the same processor). 
  */
  ParaFIELD::ParaFIELD(MEDCouplingFieldDouble* subdomain_field, ParaMESH *sup, const ProcessorGroup& proc_group):
    _field(subdomain_field),
    _component_topology(ComponentTopology(_field->getNumberOfComponents())),_topology(0),_own_support(false),
    _support(sup)
  {
    if(_field)
      _field->incrRef();
    const BlockTopology* source_topo=dynamic_cast<const BlockTopology*> (_support->getTopology());
    _topology=new BlockTopology(*source_topo,_component_topology.nbLocalComponents());
  }

  ParaFIELD::~ParaFIELD()
  {
    if(_field)
      _field->decrRef();
    if(_own_support)
      delete _support;
    delete _topology;
  }

  void ParaFIELD::synchronizeTarget(ParaFIELD* source_field)
  {
    DisjointDEC* data_channel;
    if (dynamic_cast<BlockTopology*>(_topology)!=0)
      {
        data_channel=new StructuredCoincidentDEC;
      }
    else
      {
        data_channel=new ExplicitCoincidentDEC;
      }
    data_channel->attachLocalField(this);
    data_channel->synchronize();
    data_channel->prepareTargetDE();
    data_channel->recvData();
  
    delete data_channel;
  }

  void ParaFIELD::synchronizeSource(ParaFIELD* target_field)
  {
    DisjointDEC* data_channel;
    if (dynamic_cast<BlockTopology*>(_topology)!=0)
      {
        data_channel=new StructuredCoincidentDEC;
      }
    else
      {
        data_channel=new ExplicitCoincidentDEC;
      }
    data_channel->attachLocalField(this);
    data_channel->synchronize();
    data_channel->prepareSourceDE();
    data_channel->sendData();
  
    delete data_channel;
  }

  /*!
   * This method returns, if it exists, an array with only one component and as many as tuples as _field has.
   * This array gives for every element on which this->_field lies, its global number, if this->_field is nodal.
   * For example if _field is a nodal field : returned array will be the nodal global numbers.
   * The content of this method is used to inform Working side to accumulate data recieved by lazy side.
   */
  DataArrayInt* ParaFIELD::returnCumulativeGlobalNumbering() const
  {
    if(!_field)
      return 0;
    TypeOfField type=_field->getTypeOfField();
    switch(type)
      {
      case ON_CELLS:
        return 0;
      case ON_NODES:
        return _support->getGlobalNumberingNodeDA();
      default:
        return 0;
      }
  }

  DataArrayInt* ParaFIELD::returnGlobalNumbering() const
  {
    if(!_field)
      return 0;
    TypeOfField type=_field->getTypeOfField();
    switch(type)
      {
      case ON_CELLS:
        return _support->getGlobalNumberingCellDA();
      case ON_NODES:
        return _support->getGlobalNumberingNodeDA();
      default:
        return 0;
      }
  }
  
  int ParaFIELD::nbComponents() const
  {
    return _component_topology.nbComponents();
  }


  /*! This method retrieves the integral of component \a icomp
    over the all domain. */
  double ParaFIELD::getVolumeIntegral(int icomp, bool isWAbs) const
  {
    CommInterface comm_interface = _topology->getProcGroup()->getCommInterface();
    double integral=_field->integral(icomp,isWAbs);
    double total=0.;
    const MPI_Comm* comm = (dynamic_cast<const MPIProcessorGroup*>(_topology->getProcGroup()))->getComm();
    comm_interface.allReduce(&integral, &total, 1, MPI_DOUBLE, MPI_SUM, *comm);
  
    return total;
  }
}
