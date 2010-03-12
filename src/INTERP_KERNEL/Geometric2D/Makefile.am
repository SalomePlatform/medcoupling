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

noinst_LTLIBRARIES = libInterpGeometric2DAlg.la

dist_libInterpGeometric2DAlg_la_SOURCES = \
InterpKernelGeo2DAbstractEdge.cxx         \
InterpKernelGeo2DBounds.cxx               \
InterpKernelGeo2DPrecision.cxx            \
InterpKernelGeo2DComposedEdge.cxx         \
InterpKernelGeo2DEdgeArcCircle.cxx        \
InterpKernelGeo2DEdge.cxx                 \
InterpKernelGeo2DEdgeInfLin.cxx           \
InterpKernelGeo2DEdgeLin.cxx              \
InterpKernelGeo2DElementaryEdge.cxx       \
InterpKernelGeo2DNode.cxx                 \
InterpKernelGeo2DQuadraticPolygon.cxx

salomeinclude_HEADERS =                \
INTERPKERNELGEOMETRIC2DDefines.hxx     \
InterpKernelGeo2DAbstractEdge.hxx      \
InterpKernelGeo2DBounds.hxx            \
InterpKernelGeo2DPrecision.hxx         \
InterpKernelGeo2DComposedEdge.hxx      \
InterpKernelGeo2DEdgeArcCircle.hxx     \
InterpKernelGeo2DEdge.hxx              \
InterpKernelGeo2DEdgeInfLin.hxx        \
InterpKernelGeo2DEdgeLin.hxx           \
InterpKernelGeo2DElementaryEdge.hxx    \
InterpKernelGeo2DNode.hxx              \
InterpKernelGeo2DQuadraticPolygon.hxx  \
InterpKernelGeo2DEdge.txx

libInterpGeometric2DAlg_la_CPPFLAGS = -I$(srcdir)/../Bases
