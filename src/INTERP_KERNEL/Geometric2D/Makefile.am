#  MED MEDMEM : MED files in memory
#
#  Copyright (C) 2003  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
#  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS 
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
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#
#
#
#  File   : Makefile.am
#  Author : Anthony GEAY (CEA/DEN/DANS/DM2S/SFME/LGLS)
#  Module : MED

include $(top_srcdir)/adm_local/unix/make_common_starter.am

noinst_LTLIBRARIES = libInterpGeometric2DAlg.la

dist_libInterpGeometric2DAlg_la_SOURCES = \
AbstractEdge.hxx     \
Bounds.hxx           \
Precision.hxx        \
ComposedEdge.hxx     \
ComposedEdgeWithIt.hxx \
EdgeArcCircle.hxx    \
Edge.hxx             \
Edge.txx             \
EdgeInfLin.hxx       \
EdgeLin.hxx          \
ElementaryEdge.hxx   \
Node.hxx             \
QuadraticPolygon.hxx \
AbstractEdge.cxx     \
Bounds.cxx           \
Precision.cxx        \
ComposedEdge.cxx     \
EdgeArcCircle.cxx    \
Edge.cxx             \
EdgeInfLin.cxx       \
EdgeLin.cxx          \
ElementaryEdge.cxx   \
Node.cxx             \
QuadraticPolygon.cxx \
Geometric2D_defines.hxx

libInterpGeometric2DAlg_la_CPPFLAGS = -I$(srcdir)/..

EXTRA_DIST += \
	BasicMainTest.hxx \
	Test

