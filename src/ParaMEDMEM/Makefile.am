#  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

include $(top_srcdir)/adm_local/unix/make_common_starter.am

lib_LTLIBRARIES= libparamedmem.la

salomeinclude_HEADERS= \
CommInterface.hxx\
MPIProcessorGroup.hxx\
ProcessorGroup.hxx\
BlockTopology.hxx\
Topology.hxx\
ParaGRID.hxx\
ParaMESH.hxx\
ComponentTopology.hxx\
ExplicitTopology.hxx\
ParaFIELD.hxx\
DEC.hxx\
DECOptions.hxx\
MxN_Mapping.hxx\
StructuredCoincidentDEC.hxx\
InterpolationMatrix.hxx\
InterpKernelDEC.hxx\
ExplicitCoincidentDEC.hxx\
ElementLocator.hxx\
ExplicitMapping.hxx\
ICoCoField.h \
ICoCoMEDField.hxx \
ICoCoTrioField.h \
MPIAccess.hxx \
MPIAccessDEC.hxx \
TimeInterpolator.hxx \
LinearTimeInterpolator.hxx

dist_libparamedmem_la_SOURCES= \
ProcessorGroup.cxx \
MPIProcessorGroup.cxx\
ParaMESH.cxx\
ComponentTopology.cxx\
MPIAccess.cxx \
InterpolationMatrix.cxx\
StructuredCoincidentDEC.cxx\
ExplicitCoincidentDEC.cxx\
InterpKernelDEC.cxx\
ElementLocator.cxx\
MPIAccessDEC.cxx \
TimeInterpolator.cxx \
LinearTimeInterpolator.cxx\
DEC.cxx\
ExplicitTopology.cxx\
MxN_Mapping.cxx\
ICoCoMEDField.cxx\
ICoCoField.cxx\
ICoCoTrioField.cxx\
ParaFIELD.cxx\
ParaGRID.cxx\
BlockTopology.cxx

#libmedmem_la_LDFLAGS= -L$(top_builddir)/lib@LIB_LOCATION_SUFFIX@/salome 
libparamedmem_la_CPPFLAGS= $(MPI_INCLUDES) @CXXTMPDPTHFLAGS@ \
	-I$(srcdir)/../INTERP_KERNEL \
	-I$(srcdir)/../INTERP_KERNEL/Bases \
	-I$(srcdir)/../INTERP_KERNEL/Geometric2D \
	-I$(srcdir)/../MEDCoupling

# change motivated by the bug KERNEL4778.
libparamedmem_la_LDFLAGS= ../MEDCoupling/libmedcoupling.la \
../INTERP_KERNEL/libinterpkernel.la \
$(MPI_LIBS)

if MED_ENABLE_FVM
	salomeinclude_HEADERS += NonCoincidentDEC.hxx
	libparamedmem_la_SOURCES+= NonCoincidentDEC.cxx
	libparamedmem_la_LDFLAGS+= $(FVM_LIBS)
	libparamedmem_la_CPPFLAGS+= $(FVM_INCLUDES)
endif

LDADD= -lm

EXTRA_DIST += \
	BASICS_JR \
	CommInterface.cxx \
	NonCoincidentDEC.cxx \
	NonCoincidentDEC.hxx \
	DECOptions.hxx \
	README_JR \
	TODO_JR \
	Topology.cxx
