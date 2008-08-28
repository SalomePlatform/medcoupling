#ifndef __INTERPKERNEL_DEFINES_HXX__
#define __INTERPKERNEL_DEFINES_HXX__

//export symbols
#ifdef WIN32
# ifdef INTERPKERNEL_EXPORTS
#  define INTERPKERNEL_EXPORT __declspec(dllexport)
# else
#  define INTERPKERNEL_EXPORT __declspec(dllimport)
# endif
#else
# define INTERPKERNEL_EXPORT
#endif 

#endif //__INTERPKERNEL_DEFINES_HXX__
