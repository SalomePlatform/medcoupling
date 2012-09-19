# Copyright (C) 2007-2012  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#
# Author : Anthony Geay (CEA/DEN)

include $(top_srcdir)/adm_local/unix/make_common_starter.am

lib_LTLIBRARIES = libparamedloader.la

salomeinclude_HEADERS= \
ParaMEDLoader.hxx

dist_libparamedloader_la_SOURCES= \
ParaMEDLoader.cxx

#libmedmem_la_LDFLAGS= -L$(top_builddir)/lib@LIB_LOCATION_SUFFIX@/salome 
libparamedloader_la_CPPFLAGS= $(MPI_INCLUDES) $(MED3_INCLUDES) $(HDF5_INCLUDES) @CXXTMPDPTHFLAGS@ \
	-I$(srcdir)/../INTERP_KERNEL \
	-I$(srcdir)/../INTERP_KERNEL/Geometric2D \
	-I$(srcdir)/../INTERP_KERNEL/Bases \
	-I$(srcdir)/../MEDCoupling \
	-I$(srcdir)/../MEDLoader \
	-I$(srcdir)/../ParaMEDMEM

libparamedloader_la_LDFLAGS= ../ParaMEDMEM/libparamedmem.la \
	../MEDLoader/libmedloader.la \
	../MEDCoupling/libmedcoupling.la \
	../INTERP_KERNEL/libinterpkernel.la $(MPI_LIBS) $(MED3_LIBS_C_ONLY) $(HDF5_LIBS)
