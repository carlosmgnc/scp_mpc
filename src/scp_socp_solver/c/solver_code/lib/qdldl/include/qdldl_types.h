#ifndef QDLDL_TYPES_H
# define QDLDL_TYPES_H

# ifdef __cplusplus
extern "C" {
# endif /* ifdef __cplusplus */

#include <limits.h> //for the QDLDL_INT_TYPE_MAX

// QDLDL integer and float types

typedef int    QDLDL_int;   /* for indices */
typedef double  QDLDL_float; /* for numerical values  */
typedef unsigned char   QDLDL_bool;  /* for boolean values  */

//Maximum value of the signed type QDLDL_int.
#define QDLDL_INT_MAX INT_MAX

/*
 * Configuration options
 */

/* When defined, QDLDL is using floats instead of doubles */
/* #undef QDLDL_FLOAT */

/* When defined, QDLDL is using long long instead of int types */
/* #undef QDLDL_LONG */

# ifdef __cplusplus
}
# endif /* ifdef __cplusplus */

#endif /* ifndef QDLDL_TYPES_H */
