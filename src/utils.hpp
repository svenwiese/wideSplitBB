#ifndef __UTIL_H
#define __UTIL_H

#include <cmath>

#define BUFFERSIZE  516

/* Floating point constants */
#define BIGREAL  1E+75

/* Output level */
#define OUTPUTLEVEL 1

/* Math utility */
#define XMAX(x,y)   ((x)<(y)?(y):(x))
#define XMIN(x,y)   ((x)>(y)?(y):(x))
#define XABS        fabs

/* Tolerances to be tuned */
#define EPSZERO  1E-10
#define EPSRHS   1E-06
#define EPSINT   1E-05
#define EPSVIOL  1E-05

/* Boolean constants  */
#define TRUE   1
#define FALSE  0

/* Assertions */
#define USE_ASSERT
#ifdef USE_ASSERT
#include <cassert>
#define myassert(x) assert (x)
#else 
#define assert(x,msg)
#endif

/* Error codes */
#define ERR_NOMEMORY       -1
#define ERR_BADARGUMENT    -2
#define ERR_BADPROBLEM     -3
#define ERR_BADLPSTAT      -4
#define ERR_BADFILEFORMAT  -5
#define ERR_BADPARAMVALUE  -6
#define ERR_OPENFILE       -7
#define ERR_BADMIPSTAT     -8

/* Free and set to NULL */
#define FREEN(pptr) do {                                        \
      if ( (*(pptr)) != NULL ) {                                \
         free (*pptr);                                          \
         *pptr = NULL;                                          \
      }                                                         \
   } while (0)

/* Free and set to NULL for matrix */
#define FREEN_mat(pptr,nrow)    {                           \
    if ( (*(pptr)) != NULL ) {                              \
        for(int i=0; i<nrow; i++){                          \
            free ( (*pptr)[i] );                            \
        }                                                   \
        free (*pptr);                                       \
    *pptr = NULL;                                           \
    }                                                       \
}


#endif
