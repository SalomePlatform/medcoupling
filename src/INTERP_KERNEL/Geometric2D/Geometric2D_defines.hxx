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
#ifndef __GEOMETRIC2D_DEFINES_HXX__
#define __GEOMETRIC2D_DEFINES_HXX__

//export symbols
#ifdef WIN32
# ifdef GEOMETRIC2D_EXPORTS
#  define GEOMETRIC2D_EXPORT __declspec(dllexport)
# else
#  define GEOMETRIC2D_EXPORT __declspec(dllimport)
# endif
#else
# define GEOMETRIC2D_EXPORT
#endif 

#ifdef WIN32
# include <stdlib.h>
# define fmax __max
# define fmin __min
#endif

#endif //__GEOMETRIC2D_DEFINES_HXX__
