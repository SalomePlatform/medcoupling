#  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#  File   : Makefile.am
#  Author : Anthony GEAY (CEA/DEN/DANS/DM2S/SFME/LGLS)
#  Module : MED
#
include $(top_srcdir)/adm_local/unix/make_common_starter.am

noinst_LTLIBRARIES = libinterpkernelexpreval.la

salomeinclude_HEADERS = 		\
INTERPKERNELEXPREVALDefines.hxx         \
InterpKernelExprParser.hxx              \
InterpKernelFunction.hxx                \
InterpKernelUnit.hxx                    \
InterpKernelValue.hxx


EXTRA_DIST +=		        	\
INTERPKERNELEXPREVALDefines.hxx         \
InterpKernelExprParser.hxx              \
InterpKernelFunction.hxx                \
InterpKernelUnit.hxx                    \
InterpKernelValue.hxx

dist_libinterpkernelexpreval_la_SOURCES = \
	InterpKernelExprParser.cxx        \
	InterpKernelFunction.cxx          \
	InterpKernelUnit.cxx              \
	InterpKernelValue.cxx

libinterpkernelexpreval_la_CPPFLAGS=-I$(srcdir)/../Bases

libinterpkernelexpreval_la_LDFLAGS=

AM_CPPFLAGS += $(libinterpkernelexpreval_la_CPPFLAGS)

LDADD= $(libinterpkernelexpreval_la_LDFLAGS)
