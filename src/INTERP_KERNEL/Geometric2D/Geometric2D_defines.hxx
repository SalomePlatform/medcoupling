
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
# include <math.h>
# define fmax __max
# define fmin __min
#endif

#endif //__GEOMETRIC2D_DEFINES_HXX__