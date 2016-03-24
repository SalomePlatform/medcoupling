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

#ifndef _LOG_H_
#define _LOG_H_

/** 
 * \file Log.hxx
 * \brief Simple pre-processor logging utility.
 *
 * Replaces LOG( lvl, x ) with "if(lvl <= LOG_LEVEL) std::cout << x << std::endl" when logging is active 
 * (LOG_LEVEL > 0 is defined). x is the level at which the message should be logged - if it is smaller or equal to
 * LOG_LEVEL (which can be defined at compile-time for each file by passing option -DLOG_LEVEL=x to gcc)
 * than the message is logged.
 *
 *
 *
 */

/// define LOG_LEVEL here if it is not already defined
#ifndef LOG_LEVEL
#define LOG_LEVEL 0
#endif

#if LOG_LEVEL > 0

#include <iostream>

/// write message msg to std::cout if x <= LOG_LEVEL
#define LOG(x, msg) if(x <= LOG_LEVEL) std::cout << msg << std::endl; 
#define LOG3( x , msg1 , msg2 ) if(x <= LOG_LEVEL) std::cout << msg1, msg2 << std::endl; 

#else

#define LOG( x , msg )
#define LOG3( x , msg1 , msg2 )

#endif











#endif
