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
