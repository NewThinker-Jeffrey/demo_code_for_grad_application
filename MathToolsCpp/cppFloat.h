#ifndef __MathtoolsCpp__CppFloat__
#define __MathtoolsCpp__CppFloat__


//#define USE_FLOAT64_IN_CPP


#define START_JK_Math_NAMESPACE       namespace JK { namespace Math {
#define END_JK_Math_NAMESPACE         }}
#define IMPORT_JK_Math                using namespace JK::Math;
#define JKMath JK::Math


#ifdef USE_FLOAT64_IN_CPP
    typedef double Float ;


#else
    typedef float Float ;


#endif


#ifndef NAN
#define NAN (sqrt(-1.0))
#endif



#endif
