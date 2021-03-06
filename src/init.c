#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(lgd1v)(double*,double*,double*,double*,int*,int*,double*,double*);
extern void F77_NAME(lgdvii)(double*,double*,double*,double*,int*,int*,int*,double*,double*);
extern void F77_NAME(lgdvvi)(double*,double*,double*,double*,double*,int*,int*,int*,double*,double*);

static const R_FortranMethodDef FortranEntries[] = {
  {"lgd1v",     (DL_FUNC) &F77_NAME(lgd1v),        8},
  {"lgdvii",    (DL_FUNC) &F77_NAME(lgdvii),       9},
  {"lgdvvi",    (DL_FUNC) &F77_NAME(lgdvvi),      10}, //
    {NULL, NULL, 0}
};

void R_init_probout(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

