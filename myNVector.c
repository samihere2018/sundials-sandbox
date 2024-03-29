/* -----------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2024, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This is the implementation file for a custom N_Vector
 * (note that this template just implements the serial N_Vector).
 * -----------------------------------------------------------------*/

#include "myNVector.h"
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_core.h>
#include "sundials/sundials_errors.h"

#define ZERO   SUN_RCONST(0.0) // check every occurance!
#define HALF   SUN_RCONST(0.5) // check every occurance!
#define ONE    SUN_RCONST(1.0) // check every occurance!
#define ONEPT5 SUN_RCONST(1.5) // check every occurance!

#define CZERO       (ZERO + ZERO*I)
#define HALFiHALF   (HALF + HALF*I) // check every occurance!
#define CONE        (ONE + ONE*I) // check every occurance!

/* Private functions for special cases of vector operations */
static void VCopy_Mine(N_Vector x, N_Vector z);             /* z=x       */
static void VSum_Mine(N_Vector x, N_Vector y, N_Vector z);  /* z=x+y     */
static void VDiff_Mine(N_Vector x, N_Vector y, N_Vector z); /* z=x-y     */
static void VNeg_Mine(N_Vector x, N_Vector z);              /* z=-x      */
static void VScaleSum_Mine(suncomplextype c, N_Vector x, N_Vector y,
                           N_Vector z); /* z=c(x+y)  */
static void VScaleDiff_Mine(suncomplextype c, N_Vector x, N_Vector y,
                            N_Vector z); /* z=c(x-y)  */
static void VLin1_Mine(suncomplextype a, N_Vector x, N_Vector y,
                       N_Vector z); /* z=ax+y    */
static void VLin2_Mine(suncomplextype a, N_Vector x, N_Vector y,
                       N_Vector z);                            /* z=ax-y    */
static void Vaxpy_Mine(suncomplextype a, N_Vector x, N_Vector y); /* y <- ax+y */
static void VScaleBy_Mine(suncomplextype a, N_Vector x);          /* x <- ax   */

/* Private functions for special cases of vector array operations */
static void VSumVectorArray_Mine(int nvec, N_Vector* X, N_Vector* Y,
                                 N_Vector* Z); /* Z=X+Y     */
static void VDiffVectorArray_Mine(int nvec, N_Vector* X, N_Vector* Y,
                                  N_Vector* Z); /* Z=X-Y     */
static void VScaleSumVectorArray_Mine(int nvec, suncomplextype c, N_Vector* X,
                                      N_Vector* Y, N_Vector* Z); /* Z=c(X+Y)  */
static void VScaleDiffVectorArray_Mine(int nvec, suncomplextype c, N_Vector* X,
                                       N_Vector* Y, N_Vector* Z); /* Z=c(X-Y)  */
static void VLin1VectorArray_Mine(int nvec, suncomplextype a, N_Vector* X,
                                  N_Vector* Y, N_Vector* Z); /* Z=aX+Y    */
static void VLin2VectorArray_Mine(int nvec, suncomplextype a, N_Vector* X,
                                  N_Vector* Y, N_Vector* Z); /* Z=aX-Y    */
static void VaxpyVectorArray_Mine(int nvec, suncomplextype a, N_Vector* X,
                                  N_Vector* Y); /* Y <- aX+Y */

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new empty serial vector
 */

N_Vector N_VNewEmpty_Mine(sunindextype length, SUNContext sunctx)
{

  N_Vector v;
  My_NVectorContent content;

  if (length < 0) { return NULL; }

  /* Create an empty vector object */
  v = NULL;
  v = N_VNewEmpty(sunctx);

  /* Attach operations */

  /* constructors, destructors, and utility operations */
  v->ops->nvgetvectorid     = N_VGetVectorID_Mine;
  v->ops->nvclone           = N_VClone_Mine;
  v->ops->nvcloneempty      = N_VCloneEmpty_Mine;
  v->ops->nvdestroy         = N_VDestroy_Mine;
  v->ops->nvspace           = N_VSpace_Mine;

  // v->ops->nvgetarraypointer = N_VGetArrayPointer_Mine;
  v->ops->nvgetarraypointer = N_VGetArrayPointer_Real;
  
  // v->ops->nvsetarraypointer = N_VSetArrayPointer_Mine;
  v->ops->nvsetarraypointer = N_VSetArrayPointer_Real;
  
  v->ops->nvgetlength       = N_VGetLength_Mine;
  v->ops->nvgetlocallength  = N_VGetLength_Mine;

  /* standard vector operations */
  // v->ops->nvlinearsum    = N_VLinearSum_Mine;
  v->ops->nvlinearsum    = N_VLinearSum_Real;
  
  // v->ops->nvconst        = N_VConst_Mine;
  v->ops->nvconst        = N_VConst_Real;

  v->ops->nvprod         = N_VProd_Mine;
  v->ops->nvdiv          = N_VDiv_Mine;
  
  // v->ops->nvscale        = N_VScale_Mine;
  v->ops->nvscale        = N_VScale_Real;
  
  v->ops->nvabs          = N_VAbs_Mine;
  v->ops->nvinv          = N_VInv_Mine;

  // v->ops->nvaddconst     = N_VAddConst_Mine;
  v->ops->nvaddconst     = N_VAddConst_Real;

  // v->ops->nvdotprod      = N_VDotProd_Mine;
  v->ops->nvdotprod      = N_VDotProd_Real;

  v->ops->nvmaxnorm      = N_VMaxNorm_Mine;
  v->ops->nvwrmsnormmask = N_VWrmsNormMask_Mine;
  v->ops->nvwrmsnorm     = N_VWrmsNorm_Mine;
  v->ops->nvmin          = N_VMin_Mine;
  v->ops->nvwl2norm      = N_VWL2Norm_Mine;
  v->ops->nvl1norm       = N_VL1Norm_Mine;
  v->ops->nvcompare      = N_VCompare_Mine;
  v->ops->nvinvtest      = N_VInvTest_Mine;
  v->ops->nvconstrmask   = N_VConstrMask_Mine;
  v->ops->nvminquotient  = N_VMinQuotient_Mine;

  /* fused and vector array operations are disabled (NULL) by default */

  /* local reduction operations */
  // v->ops->nvdotprodlocal     = N_VDotProd_Mine;
  v->ops->nvdotprodlocal     = N_VDotProd_Real;

  v->ops->nvmaxnormlocal     = N_VMaxNorm_Mine;
  v->ops->nvminlocal         = N_VMin_Mine;
  v->ops->nvl1normlocal      = N_VL1Norm_Mine;
  v->ops->nvinvtestlocal     = N_VInvTest_Mine;
  v->ops->nvconstrmasklocal  = N_VConstrMask_Mine;
  v->ops->nvminquotientlocal = N_VMinQuotient_Mine;
  v->ops->nvwsqrsumlocal     = N_VWSqrSumLocal_Mine;
  v->ops->nvwsqrsummasklocal = N_VWSqrSumMaskLocal_Mine;

  /* single buffer reduction operations */
  // v->ops->nvdotprodmultilocal = N_VDotProdMulti_Mine;
  v->ops->nvdotprodmultilocal = N_VDotProdMulti_Real;

  /* XBraid interface operations */
  v->ops->nvbufsize   = N_VBufSize_Mine;
  v->ops->nvbufpack   = N_VBufPack_Mine;
  v->ops->nvbufunpack = N_VBufUnpack_Mine;

  /* debugging functions */
  v->ops->nvprint     = N_VPrint_Mine;
  v->ops->nvprintfile = N_VPrintFile_Mine;

  /* Create content */
  content = NULL;
  content = (My_NVectorContent)malloc(sizeof *content);
  if (content == NULL) { return NULL; }

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->length   = length;
  content->own_data = SUNFALSE;
  content->real_data = NULL;
  content->complex_data = NULL;

  return (v);
}

/* ----------------------------------------------------------------------------
 * Function to create a new serial vector
 */

N_Vector N_VNew_Mine(sunindextype length, SUNContext sunctx)
{
  N_Vector v;
  suncomplextype* complex_data;

  if (length < 0) { return NULL; }

  v = NULL;
  v = N_VNewEmpty_Mine(length, sunctx);
  if (v == NULL) { return NULL; }

  /* Create data */
  complex_data = NULL;
  if (length > 0)
  {
    complex_data = (suncomplextype*)malloc(length * sizeof(suncomplextype));
    if (complex_data == NULL) { return NULL; }
  }

  /* Attach data */
  MYNV_OWN_DATA(v) = SUNTRUE;
  MYNV_COMPLEX_DATA(v) = complex_data;

  return (v);
}

/* ----------------------------------------------------------------------------
 * Function to create a serial N_Vector with user data component
 */

N_Vector N_VMake_Mine(sunindextype length, suncomplextype* v_data,
                        SUNContext sunctx)
{
  N_Vector v;

  if (length < 0) { return NULL; }

  v = NULL;
  v = N_VNewEmpty_Mine(length, sunctx);
  if (v == NULL) { return NULL; }

  if (length > 0)
  {
    /* Attach data */
    MYNV_OWN_DATA(v) = SUNFALSE;
    MYNV_COMPLEX_DATA(v) = v_data;
  }

  return (v);
}

/* ----------------------------------------------------------------------------
 * Function to return number of vector elements
 */
sunindextype N_VGetLength_Mine(N_Vector v) { return MYNV_LENGTH(v); }

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to stdout
 */

void N_VPrint_Mine(N_Vector x)
{
  N_VPrintFile_Mine(x, stdout);
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector to outfile
 */

void N_VPrintFile_Mine(N_Vector x, FILE* outfile)
{
  sunindextype i, N;
  suncomplextype* xd;

  xd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);

  for (i = 0; i < N; i++)
  {
#if defined(SUNDIALS_EXTENDED_PRECISION)
    fprintf(outfile, "%35.32Le + i%35.32Le\n", creal(xd[i]), cimag(xd[i]));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    fprintf(outfile, "%19.16e + i%19.16e\n", creal(xd[i]), cimag(xd[i]));
#else
    fprintf(outfile, "%11.8e + i%11.8e\n", creal(xd[i]), cimag(xd[i]));
#endif
  }
  fprintf(outfile, "\n");

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_Mine(N_Vector w)
{
  N_Vector v;
  My_NVectorContent content;

  /* Create vector */
  v = NULL;
  v = N_VNewEmpty(w->sunctx);
  if (v == NULL) { return NULL; }

  /* Attach operations */
  N_VCopyOps(w, v);

  /* Create content */
  content = NULL;
  content = (My_NVectorContent)malloc(sizeof *content);
  if (content == NULL) { return NULL; }

  /* Attach content */
  v->content = content;

  /* Initialize content */
  content->length   = MYNV_LENGTH(w);
  content->own_data = SUNFALSE;
  content->real_data = NULL;
  content->complex_data = NULL;

  return (v);
}

N_Vector N_VClone_Mine(N_Vector w)
{
  N_Vector v;
  sunrealtype* real_data;
  suncomplextype* complex_data;
  sunindextype length;

  v = NULL;
  v = N_VCloneEmpty_Mine(w);
  if (v == NULL) { return NULL; }

  length = MYNV_LENGTH(w);

  /* Create complex_data */
  complex_data = MYNV_COMPLEX_DATA(v) = NULL;
  real_data = MYNV_REAL_DATA(v) = NULL;
  if (length > 0)
  {
    complex_data = (suncomplextype*)malloc(length * sizeof(suncomplextype));
    if (complex_data == NULL) { return NULL; }
    
    real_data = (sunrealtype*)malloc(length * sizeof(sunrealtype));
    if (real_data == NULL) { return NULL; }

    /* Attach real_data */
    MYNV_OWN_DATA(v) = SUNTRUE;
    MYNV_REAL_DATA(v)     = real_data;

    /* Attach complex_data */
    MYNV_COMPLEX_DATA(v)     = complex_data;
  }

  return (v);
}

void N_VDestroy_Mine(N_Vector v)
{
  if (v == NULL) { return; }

  /* free content */
  if (v->content != NULL)
  {
    /* free data array if it's owned by the vector */
    if (MYNV_OWN_DATA(v) && (MYNV_COMPLEX_DATA(v) != NULL))
    {
      free(MYNV_COMPLEX_DATA(v));
      MYNV_COMPLEX_DATA(v) = NULL;
    }
    if (MYNV_OWN_DATA(v) && (MYNV_REAL_DATA(v) != NULL))
    {
      free(MYNV_REAL_DATA(v));
      MYNV_REAL_DATA(v) = NULL;
    }
    free(v->content);
    v->content = NULL;
  }

  /* free ops and vector */
  if (v->ops != NULL)
  {
    free(v->ops);
    v->ops = NULL;
  }
  free(v);
  v = NULL;

  return;
}

void N_VSpace_Mine(N_Vector v, sunindextype* lrw, sunindextype* liw)
{
  if (lrw == NULL) { return; }
  if (liw == NULL) { return; }

  *lrw = MYNV_LENGTH(v);
  *liw = 1;

  return;
}

suncomplextype* N_VGetArrayPointer_Mine(N_Vector v)
{
    
  return ((suncomplextype*)MYNV_COMPLEX_DATA(v));
}

sunrealtype* N_VGetArrayPointer_Real(N_Vector v)
{
    
  return ((sunrealtype*)MYNV_REAL_DATA(v));
}

void N_VSetArrayPointer_Mine(suncomplextype* v_data, N_Vector v)
{
  if (MYNV_LENGTH(v) > 0) { MYNV_COMPLEX_DATA(v) = v_data; }

  return;
}

void N_VSetArrayPointer_Real(sunrealtype* v_data, N_Vector v)
{
  if (MYNV_LENGTH(v) > 0) { MYNV_REAL_DATA(v) = v_data; }

  return;
}

void N_VLinearSum_Mine(suncomplextype a, N_Vector x, suncomplextype b, N_Vector y,
                       N_Vector z)
{
  sunindextype i, N;
  suncomplextype c, *xd, *yd, *zd;
  N_Vector v1, v2;
  sunbooleantype test;

  xd = yd = zd = NULL;

  if (((creal(b) == ONE) && (cimag(b) == ZERO)) && (z == y))
  { /* BLAS usage: axpy y <- ax+y */
    Vaxpy_Mine(a, x, y);
    return;
  }

  if (((creal(a) == ONE) && (cimag(a) == ZERO)) && (z == x))
  { /* BLAS usage: axpy x <- by+x */
    Vaxpy_Mine(b, y, x);
    return;
  }

  /* Case: a == b == 1.0 */

  if (((creal(a) == ONE) && (cimag(a) == ZERO)) && ((creal(b) == ONE) && (cimag(b) == ZERO)))
  {
    VSum_Mine(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = (((creal(a) == ONE) && (cimag(a) == ZERO)) && ((creal(b) == -ONE) && (cimag(b) == ZERO)))) || (((creal(a) == -ONE) && (cimag(a) == ZERO)) && ((creal(b) == ONE) && (cimag(b) == ZERO))))
  {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_Mine(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = ((creal(a) == ONE) && (cimag(a) == ZERO))) || ((creal(b) == ONE) && (cimag(b) == ZERO)))
  {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_Mine(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = ((creal(a) == -ONE) && (cimag(a) == ZERO))) || ((creal(b) == -ONE) && (cimag(b) == ZERO)))
  {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_Mine(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b)
  {
    VScaleSum_Mine(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b)
  {
    VScaleDiff_Mine(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  yd = MYNV_COMPLEX_DATA(y);
  zd = MYNV_COMPLEX_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = (a * xd[i]) + (b * yd[i]); }

  return;
}

void N_VLinearSum_Real(sunrealtype a, N_Vector x, sunrealtype b, N_Vector y,
                       N_Vector z)
{
  sunindextype i, N;
  sunrealtype c, *xd, *yd, *zd;
  N_Vector v1, v2;
  sunbooleantype test;

  xd = yd = zd = NULL;

  if ((b == ONE) && (z == y))
  { /* BLAS usage: axpy y <- ax+y */
    Vaxpy_Mine(a, x, y);
    return;
  }

  if ((a == ONE) && (z == x))
  { /* BLAS usage: axpy x <- by+x */
    Vaxpy_Mine(b, y, x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE))
  {
    VSum_Mine(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE)))
  {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff_Mine(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE))
  {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1_Mine(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE))
  {
    c  = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2_Mine(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b)
  {
    VScaleSum_Mine(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b)
  {
    VScaleDiff_Mine(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  N  = MYNV_LENGTH(x);
  xd = MYNV_REAL_DATA(x);
  yd = MYNV_REAL_DATA(y);
  zd = MYNV_REAL_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = (a * xd[i]) + (b * yd[i]); }

  return;
}

void N_VConst_Mine(suncomplextype c, N_Vector z)
{
  sunindextype i, N;
  suncomplextype* zd;

  zd = NULL;

  N  = MYNV_LENGTH(z);
  zd = MYNV_COMPLEX_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = c; }

  return;
}

void N_VConst_Real(sunrealtype c, N_Vector z)
{
  sunindextype i, N;
  sunrealtype* zd;

  zd = NULL;

  N  = MYNV_LENGTH(z);
  zd = MYNV_REAL_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = c; }

  return;
}

void N_VProd_Mine(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  suncomplextype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  yd = MYNV_COMPLEX_DATA(y);
  zd = MYNV_COMPLEX_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] * yd[i]; }

  return;
}

void N_VDiv_Mine(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  suncomplextype *xd, *yd, *zd;

  xd = yd = zd = NULL;
  
  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  yd = MYNV_COMPLEX_DATA(y);
  zd = MYNV_COMPLEX_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] / yd[i]; }

  return;
}

void N_VScale_Mine(suncomplextype c, N_Vector x, N_Vector z)
{
  
  sunindextype i, N;
  suncomplextype *xd, *zd;

  xd = zd = NULL;

  if (z == x)
  { /* BLAS usage: scale x <- cx */
    VScaleBy_Mine(c, x);
    return;
  }

  if ((creal(c) == ONE) && (cimag(c) == ZERO)) { VCopy_Mine(x, z); }
  else if ((creal(c) == -ONE) && (cimag(c) == ZERO)) { VNeg_Mine(x, z); }
  else
  {
    N  = MYNV_LENGTH(x);
    xd = MYNV_COMPLEX_DATA(x);
    zd = MYNV_COMPLEX_DATA(z);
    for (i = 0; i < N; i++) { zd[i] = c * xd[i]; }
  }

  return;
}

void N_VScale_Real(sunrealtype c, N_Vector x, N_Vector z)
{
  
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  if (z == x)
  { /* BLAS usage: scale x <- cx */
    VScaleBy_Mine(c, x);
    return;
  }

  if ((creal(c) == ONE) && (cimag(c) == ZERO)) { VCopy_Mine(x, z); }
  else if ((creal(c) == -ONE) && (cimag(c) == ZERO)) { VNeg_Mine(x, z); }
  else
  {
    N  = MYNV_LENGTH(x);
    xd = MYNV_REAL_DATA(x);
    zd = MYNV_REAL_DATA(z);
    for (i = 0; i < N; i++) { zd[i] = c * xd[i]; }
  }

  return;
}


void N_VAbs_Mine(N_Vector x, N_Vector z)
{

  sunindextype i, N;
  suncomplextype *xd, *zd;

  xd = zd = NULL;
  
  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  zd = MYNV_COMPLEX_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = SUNCabs(xd[i]); }

  return;
}

void N_VInv_Mine(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  suncomplextype *xd, *zd;

  xd = zd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  zd = MYNV_COMPLEX_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = (suncomplextype)ONE / xd[i]; }

  return;
}

void N_VAddConst_Mine(N_Vector x, suncomplextype b, N_Vector z)
{
  sunindextype i, N;
  suncomplextype *xd, *zd;

  xd = zd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  zd = MYNV_COMPLEX_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] + b; }

  return;
}

void N_VAddConst_Real(N_Vector x, sunrealtype b, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  xd = zd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_REAL_DATA(x);
  zd = MYNV_REAL_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] + b; }

  return;
}

suncomplextype N_VDotProd_Mine(N_Vector x, N_Vector y)
{
  sunindextype i, N;
  suncomplextype sum, *xd, *yd;

  sum = (suncomplextype)ZERO;
  xd = yd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  yd = MYNV_COMPLEX_DATA(y);

  for (i = 0; i < N; i++) { sum += xd[i] * conj(yd[i]); }

  return (sum);
}

sunrealtype N_VDotProd_Real(N_Vector x, N_Vector y)
{
  sunindextype i, N;
  sunrealtype sum, *xd, *yd;

  sum = ZERO;
  xd = yd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_REAL_DATA(x);
  yd = MYNV_REAL_DATA(y);

  for (i = 0; i < N; i++) { sum += xd[i] * yd[i]; }

  return (sum);
}

sunrealtype N_VMaxNorm_Mine(N_Vector x)
{
  sunindextype i, N;
  sunrealtype max;
  suncomplextype *xd;

  max = ZERO;
  xd  = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);

  for (i = 0; i < N; i++)
  {
    if (SUNCabs(xd[i]) > max) { max = SUNCabs(xd[i]); }
  }

  return (max);
}

sunrealtype N_VWrmsNorm_Mine(N_Vector x, N_Vector w)
{

printf("sunrealtype N_VWrmsNorm_Mine(N_Vector x, N_Vector w) is not ready yet! Check myNvector.c\n");

  // sunrealtype norm = N_VWSqrSumLocal_Mine(x, w);
  // norm = SUNRsqrt(norm / MYNV_LENGTH(x));
  // return norm;
}

sunrealtype N_VWSqrSumLocal_Mine(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  sunrealtype sum, prodi, *xd, *wd;

printf("sunrealtype N_VWSqrSumLocal_Mine(N_Vector x, N_Vector w) is not ready yet! Check myNvector.c\n");

  sum = ZERO;
  xd = wd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_REAL_DATA(x);
  wd = MYNV_REAL_DATA(w); //change to comlex might be required

  // for (i = 0; i < N; i++)
  // {
  //   prodi = xd[i] * wd[i];
  //   sum += SUNSQR(prodi);
  // }

  return (sum);
}

sunrealtype N_VWrmsNormMask_Mine(N_Vector x, N_Vector w, N_Vector id)
{
printf("N_VWrmsNormMask_Mine(N_Vector x, N_Vector w, N_Vector id) is not ready yet! Check myNvector.c\n");

  // sunrealtype norm = N_VWSqrSumMaskLocal_Mine(x, w, id);
  // norm = SUNRsqrt(norm / MYNV_LENGTH(x));
  // return norm;
}

sunrealtype N_VWSqrSumMaskLocal_Mine(N_Vector x, N_Vector w, N_Vector id)
{
  sunindextype i, N;
  suncomplextype sum, prodi, *xd, *wd, *idd;

  printf("sunrealtype N_VWSqrSumMaskLocal_Mine(N_Vector x, N_Vector w, N_Vector id) is not ready yet! Check myNvector.c\n");

  sum = 0.0; //change to a complex ZERO!
  xd = wd = idd = NULL;

  N   = MYNV_LENGTH(x);
  xd  = MYNV_COMPLEX_DATA(x);
  wd  = MYNV_COMPLEX_DATA(w);
  idd = MYNV_COMPLEX_DATA(id);

  // for (i = 0; i < N; i++)
  // {
  //   if (creal(idd[i]) > creal(ZERO))  //it used to be "if (idd[i] > ZERO)". Change this accordingly.
  //   {
  //     prodi = xd[i] * wd[i];
  //     sum += SUNSQR(prodi);
  //   }
  // }

  return (sum);
}

sunrealtype N_VMin_Mine(N_Vector x)
{
  sunindextype i, N;
  sunrealtype min, *xd;

  printf("sunrealtype N_VMin_Mine(N_Vector x) is not ready yet! Check myNvector.c\n");

  // xd = NULL;

  // N  = MYNV_LENGTH(x);
  // xd = MYNV_REAL_DATA(x); //change to the complex data

  // min = xd[0];

  // for (i = 1; i < N; i++)
  // {
  //   if (xd[i] < min) { min = xd[i]; }
  // }

  return (min);
}

sunrealtype N_VWL2Norm_Mine(N_Vector x, N_Vector w)
{
  sunindextype i, N;
  sunrealtype sum, prodi, *xd, *wd;

  sum = ZERO;
  xd = wd = NULL;

  printf("sunrealtype N_VWL2Norm_Mine(N_Vector x, N_Vector w) is not ready yet! Check myNvector.c\n");

  N  = MYNV_LENGTH(x);
  xd = MYNV_REAL_DATA(x);
  wd = MYNV_REAL_DATA(w);

  // for (i = 0; i < N; i++)
  // {
  //   prodi = xd[i] * wd[i];
  //   sum += SUNSQR(prodi);
  // }

  return (SUNRsqrt(sum));
}

sunrealtype N_VL1Norm_Mine(N_Vector x)
{
  sunindextype i, N;
  sunrealtype sum; 
  suncomplextype *xd;

  sum = ZERO;
  xd  = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);

  for (i = 0; i < N; i++) { sum += SUNCabs(xd[i]); }

  return (sum);
}

void N_VCompare_Mine(sunrealtype c, N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;

  printf("void N_VCompare_Mine(sunrealtype c, N_Vector x, N_Vector z) is not ready yet! Check myNvector.c\n");

  xd = zd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_REAL_DATA(x);
  zd = MYNV_REAL_DATA(z);

  // for (i = 0; i < N; i++) { zd[i] = (SUNRabs(xd[i]) >= c) ? ONE : ZERO; }

  return;
}

sunbooleantype N_VInvTest_Mine(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  sunrealtype *xd, *zd;
  sunbooleantype no_zero_found;

  xd = zd = NULL;

  printf("sunbooleantype N_VInvTest_Mine(N_Vector x, N_Vector z) is not ready yet! Check myNvector.c\n");


  N  = MYNV_LENGTH(x);
  xd = MYNV_REAL_DATA(x);
  zd = MYNV_REAL_DATA(z);

  // no_zero_found = SUNTRUE;
  // for (i = 0; i < N; i++)
  // {
  //   if (xd[i] == ZERO) { no_zero_found = SUNFALSE; }
  //   else { zd[i] = ONE / xd[i]; }
  // }

  return no_zero_found;
}

sunbooleantype N_VConstrMask_Mine(N_Vector c, N_Vector x, N_Vector m)
{
  sunindextype i, N;
  sunrealtype temp;
  sunrealtype *cd, *xd, *md;
  sunbooleantype test;

  printf("sunbooleantype N_VConstrMask_Mine(N_Vector c, N_Vector x, N_Vector m) is not ready yet! Check myNvector.c\n");

  cd = xd = md = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_REAL_DATA(x);
  cd = MYNV_REAL_DATA(c);
  md = MYNV_REAL_DATA(m);

  temp = ZERO;

  for (i = 0; i < N; i++)
  {
    md[i] = ZERO;

    /* Continue if no constraints were set for the variable */
    if (cd[i] == ZERO) { continue; }

    /* Check if a set constraint has been violated */
    test = (SUNRabs(cd[i]) > ONEPT5 && xd[i] * cd[i] <= ZERO) ||
           (SUNRabs(cd[i]) > HALF && xd[i] * cd[i] < ZERO);
    if (test) { temp = md[i] = ONE; }
  }

  /* Return false if any constraint was violated */
  return (temp == ONE) ? SUNFALSE : SUNTRUE;
}

sunrealtype N_VMinQuotient_Mine(N_Vector num, N_Vector denom)
{
  sunbooleantype notEvenOnce;
  sunindextype i, N;
  sunrealtype *nd, *dd, min;

  nd = dd = NULL;

  printf("sunrealtype N_VMinQuotient_Mine(N_Vector num, N_Vector denom) is not ready yet! Check myNvector.c\n");

  N  = MYNV_LENGTH(num);
  nd = MYNV_REAL_DATA(num);
  dd = MYNV_REAL_DATA(denom);

  notEvenOnce = SUNTRUE;
  min         = SUN_BIG_REAL;

  // for (i = 0; i < N; i++)
  // {
  //   if (dd[i] == ZERO) { continue; }
  //   else
  //   {
  //     if (!notEvenOnce) { min = SUNMIN(min, nd[i] / dd[i]); }
  //     else
  //     {
  //       min         = nd[i] / dd[i];
  //       notEvenOnce = SUNFALSE;
  //     }
  //   }
  // }

  return (min);
}

/*
 * -----------------------------------------------------------------
 * fused vector operations
 * -----------------------------------------------------------------
 */

SUNErrCode N_VLinearCombination_Mine(int nvec, suncomplextype* c, N_Vector* X,
                                     N_Vector z)
{
  int i;
  sunindextype j, N;
  suncomplextype* zd = NULL;
  suncomplextype* xd = NULL;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* should have called N_VScale */
  if (nvec == 1)
  {
    N_VScale_Mine(c[0], X[0], z);
    return SUN_SUCCESS;
  }

  /* should have called N_VLinearSum */
  if (nvec == 2)
  {
    N_VLinearSum_Mine(c[0], X[0], c[1], X[1], z);
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N  = MYNV_LENGTH(z);
  zd = MYNV_COMPLEX_DATA(z);

  /*
   * X[0] += c[i]*X[i], i = 1,...,nvec-1
   */
  if ((X[0] == z) && ((creal(c[0]) == ONE) && (cimag(c[0]) == ZERO)))
  {
    for (i = 1; i < nvec; i++)
    {
      xd = MYNV_COMPLEX_DATA(X[i]);
      for (j = 0; j < N; j++) { zd[j] += c[i] * xd[j]; }
    }
    return SUN_SUCCESS;
  }

  /*
   * X[0] = c[0] * X[0] + sum{ c[i] * X[i] }, i = 1,...,nvec-1
   */
  if (X[0] == z)
  {
    for (j = 0; j < N; j++) { zd[j] *= c[0]; }
    for (i = 1; i < nvec; i++)
    {
      xd = MYNV_COMPLEX_DATA(X[i]);
      for (j = 0; j < N; j++) { zd[j] += c[i] * xd[j]; }
    }
    return SUN_SUCCESS;
  }

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */
  xd = MYNV_COMPLEX_DATA(X[0]);
  for (j = 0; j < N; j++) { zd[j] = c[0] * xd[j]; }
  for (i = 1; i < nvec; i++)
  {
    xd = MYNV_COMPLEX_DATA(X[i]);
    for (j = 0; j < N; j++) { zd[j] += c[i] * xd[j]; }
  }

  return SUN_SUCCESS;
}

SUNErrCode N_VLinearCombination_Real(int nvec, sunrealtype* c, N_Vector* X,
                                     N_Vector z)
{
  int i;
  sunindextype j, N;
  sunrealtype* zd = NULL;
  sunrealtype* xd = NULL;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* should have called N_VScale */
  if (nvec == 1)
  {
    N_VScale_Mine(c[0], X[0], z);
    return SUN_SUCCESS;
  }

  /* should have called N_VLinearSum */
  if (nvec == 2)
  {
    N_VLinearSum_Mine(c[0], X[0], c[1], X[1], z);
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N  = MYNV_LENGTH(z);
  zd = MYNV_REAL_DATA(z);

  /*
   * X[0] += c[i]*X[i], i = 1,...,nvec-1
   */
  if ((X[0] == z) && (c[0] == ONE))
  {
    for (i = 1; i < nvec; i++)
    {
      xd = MYNV_REAL_DATA(X[i]);
      for (j = 0; j < N; j++) { zd[j] += c[i] * xd[j]; }
    }
    return SUN_SUCCESS;
  }

  /*
   * X[0] = c[0] * X[0] + sum{ c[i] * X[i] }, i = 1,...,nvec-1
   */
  if (X[0] == z)
  {
    for (j = 0; j < N; j++) { zd[j] *= c[0]; }
    for (i = 1; i < nvec; i++)
    {
      xd = MYNV_REAL_DATA(X[i]);
      for (j = 0; j < N; j++) { zd[j] += c[i] * xd[j]; }
    }
    return SUN_SUCCESS;
  }

  /*
   * z = sum{ c[i] * X[i] }, i = 0,...,nvec-1
   */
  xd = MYNV_REAL_DATA(X[0]);
  for (j = 0; j < N; j++) { zd[j] = c[0] * xd[j]; }
  for (i = 1; i < nvec; i++)
  {
    xd = MYNV_REAL_DATA(X[i]);
    for (j = 0; j < N; j++) { zd[j] += c[i] * xd[j]; }
  }

  return SUN_SUCCESS;
}

SUNErrCode N_VScaleAddMulti_Mine(int nvec, suncomplextype* a, N_Vector x,
                                 N_Vector* Y, N_Vector* Z)
{  
  int i;
  sunindextype j, N;
  suncomplextype* xd = NULL;
  suncomplextype* yd = NULL;
  suncomplextype* zd = NULL;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* should have called N_VLinearSum */
  if (nvec == 1)
  {
    N_VLinearSum_Mine(a[0], x, (suncomplextype)ONE, Y[0], Z[0]);
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z)
  {
    for (i = 0; i < nvec; i++)
    {
      yd = MYNV_COMPLEX_DATA(Y[i]);
      for (j = 0; j < N; j++) { yd[j] += a[i] * xd[j]; }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i = 0; i < nvec; i++)
  {
    yd = MYNV_COMPLEX_DATA(Y[i]);
    zd = MYNV_COMPLEX_DATA(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = a[i] * xd[j] + yd[j]; }
  }
  return SUN_SUCCESS;
}

SUNErrCode N_VScaleAddMulti_Real(int nvec, sunrealtype* a, N_Vector x,
                                 N_Vector* Y, N_Vector* Z)
{  
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* should have called N_VLinearSum */
  if (nvec == 1)
  {
    N_VLinearSum_Mine(a[0], x, ONE, Y[0], Z[0]);
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N  = MYNV_LENGTH(x);
  xd = MYNV_REAL_DATA(x);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z)
  {
    for (i = 0; i < nvec; i++)
    {
      yd = MYNV_REAL_DATA(Y[i]);
      for (j = 0; j < N; j++) { yd[j] += a[i] * xd[j]; }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i = 0; i < nvec; i++)
  {
    yd = MYNV_REAL_DATA(Y[i]);
    zd = MYNV_REAL_DATA(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = a[i] * xd[j] + yd[j]; }
  }
  return SUN_SUCCESS;
}

SUNErrCode N_VDotProdMulti_Mine(int nvec, N_Vector x, N_Vector* Y,
                                suncomplextype* dotprods)
{
  int i;
  sunindextype j, N;
  suncomplextype* xd = NULL;
  suncomplextype* yd = NULL;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* should have called N_VDotProd */
  if (nvec == 1)
  {
    dotprods[0] = N_VDotProd_Mine(x, Y[0]);
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);

  /* compute multiple dot products */
  for (i = 0; i < nvec; i++)
  {
    yd          = MYNV_COMPLEX_DATA(Y[i]);
    dotprods[i] = (suncomplextype)ZERO;
    for (j = 0; j < N; j++) { dotprods[i] += xd[j] * conj(yd[j]); }
  }

  return SUN_SUCCESS;
}

SUNErrCode N_VDotProdMulti_Real(int nvec, N_Vector x, N_Vector* Y,
                                sunrealtype* dotprods)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* should have called N_VDotProd */
  if (nvec == 1)
  {
    dotprods[0] = N_VDotProd_Mine(x, Y[0]);
    return SUN_SUCCESS;
  }

  /* get vector length and data array */
  N  = MYNV_LENGTH(x);
  xd = MYNV_REAL_DATA(x);

  /* compute multiple dot products */
  for (i = 0; i < nvec; i++)
  {
    yd          = MYNV_REAL_DATA(Y[i]);
    dotprods[i] = ZERO;
    for (j = 0; j < N; j++) { dotprods[i] += xd[j] * yd[j]; }
  }

  return SUN_SUCCESS;
}


/*
 * -----------------------------------------------------------------
 * vector array operations
 * -----------------------------------------------------------------
 */

SUNErrCode N_VLinearSumVectorArray_Mine(int nvec, suncomplextype a, N_Vector* X,
                                        suncomplextype b, N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  suncomplextype* xd = NULL;
  suncomplextype* yd = NULL;
  suncomplextype* zd = NULL;
  suncomplextype c;
  N_Vector* V1;
  N_Vector* V2;
  sunbooleantype test;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* should have called N_VLinearSum */
  if (nvec == 1)
  {
    N_VLinearSum_Mine(a, X[0], b, Y[0], Z[0]);
    return SUN_SUCCESS;
  }

  /* BLAS usage: axpy y <- ax+y */
  if (((creal(b) == ONE) && (cimag(b)==ZERO)) && (Z == Y))
  {
    VaxpyVectorArray_Mine(nvec, a, X, Y);
    return SUN_SUCCESS;
  }

  /* BLAS usage: axpy x <- by+x */
  if (((creal(a) == ONE) && (cimag(a)==ZERO)) && (Z == X))
  {
    VaxpyVectorArray_Mine(nvec, b, Y, X);
    return SUN_SUCCESS;
  }

  /* Case: a == b == 1.0 */
  if (((creal(a) == ONE) && (cimag(a)==ZERO)) && ((creal(b) == ONE) && (cimag(b)==ZERO)))
  {
    VSumVectorArray_Mine(nvec, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Cases:                    */
  /*   (1) a == 1.0, b = -1.0, */
  /*   (2) a == -1.0, b == 1.0 */
  if ((test = (((creal(a) == ONE) && (cimag(a)==ZERO)) && ((creal(b) == -ONE) && (cimag(b)==ZERO)))) || (((creal(a) == -ONE) && (cimag(a)==ZERO)) && ((creal(b) == ONE) && (cimag(b)==ZERO))))
  {
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    VDiffVectorArray_Mine(nvec, V2, V1, Z);
    return SUN_SUCCESS;
  }

  /* Cases:                                                  */
  /*   (1) a == 1.0, b == other or 0.0,                      */
  /*   (2) a == other or 0.0, b == 1.0                       */
  /* if a or b is 0.0, then user should have called N_VScale */
  if ((test = ((creal(a) == ONE) && (cimag(a)==ZERO))) || ((creal(b) == ONE) && (cimag(b)==ZERO)))
  {
    c  = test ? b : a;
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    VLin1VectorArray_Mine(nvec, c, V1, V2, Z);
    return SUN_SUCCESS;
  }

  /* Cases:                     */
  /*   (1) a == -1.0, b != 1.0, */
  /*   (2) a != 1.0, b == -1.0  */
  if ((test = ((creal(a) == -ONE) && (cimag(a)==ZERO))) || ((creal(b) == -ONE) && (cimag(b)==ZERO)))
  {
    c  = test ? b : a;
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    VLin2VectorArray_Mine(nvec, c, V1, V2, Z);
    return SUN_SUCCESS;
  }

  /* Case: a == b                                                         */
  /* catches case both a and b are 0.0 - user should have called N_VConst */
  if (a == b)
  {
    VScaleSumVectorArray_Mine(nvec, a, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Case: a == -b */
  if (a == -b)
  {
    VScaleDiffVectorArray_Mine(nvec, a, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Do all cases not handled above:                               */
  /*   (1) a == other, b == 0.0 - user should have called N_VScale */
  /*   (2) a == 0.0, b == other - user should have called N_VScale */
  /*   (3) a,b == other, a !=b, a != -b                            */

  /* get vector length */
  N = MYNV_LENGTH(Z[0]);

  /* compute linear sum for each vector pair in vector arrays */
  for (i = 0; i < nvec; i++)
  {
    xd = MYNV_COMPLEX_DATA(X[i]);
    yd = MYNV_COMPLEX_DATA(Y[i]);
    zd = MYNV_COMPLEX_DATA(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = a * xd[j] + b * yd[j]; }
  }

  return SUN_SUCCESS;
}

SUNErrCode N_VLinearSumVectorArray_Real(int nvec, sunrealtype a, N_Vector* X,
                                        sunrealtype b, N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;
  sunrealtype c;
  N_Vector* V1;
  N_Vector* V2;
  sunbooleantype test;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* should have called N_VLinearSum */
  if (nvec == 1)
  {
    N_VLinearSum_Mine(a, X[0], b, Y[0], Z[0]);
    return SUN_SUCCESS;
  }

  /* BLAS usage: axpy y <- ax+y */
  if (((creal(b) == ONE) && (cimag(b)==ZERO)) && (Z == Y))
  {
    VaxpyVectorArray_Mine(nvec, a, X, Y);
    return SUN_SUCCESS;
  }

  /* BLAS usage: axpy x <- by+x */
  if (((creal(a) == ONE) && (cimag(a)==ZERO)) && (Z == X))
  {
    VaxpyVectorArray_Mine(nvec, b, Y, X);
    return SUN_SUCCESS;
  }

  /* Case: a == b == 1.0 */
  if (((creal(a) == ONE) && (cimag(a)==ZERO)) && ((creal(b) == ONE) && (cimag(b)==ZERO)))
  {
    VSumVectorArray_Mine(nvec, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Cases:                    */
  /*   (1) a == 1.0, b = -1.0, */
  /*   (2) a == -1.0, b == 1.0 */
  if ((test = (((creal(a) == ONE) && (cimag(a)==ZERO)) && ((creal(b) == -ONE) && (cimag(b)==ZERO)))) || (((creal(a) == -ONE) && (cimag(a)==ZERO)) && ((creal(b) == ONE) && (cimag(b)==ZERO))))
  {
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    VDiffVectorArray_Mine(nvec, V2, V1, Z);
    return SUN_SUCCESS;
  }

  /* Cases:                                                  */
  /*   (1) a == 1.0, b == other or 0.0,                      */
  /*   (2) a == other or 0.0, b == 1.0                       */
  /* if a or b is 0.0, then user should have called N_VScale */
  if ((test = ((creal(a) == ONE) && (cimag(a)==ZERO))) || ((creal(b) == ONE) && (cimag(b)==ZERO)))
  {
    c  = test ? b : a;
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    VLin1VectorArray_Mine(nvec, c, V1, V2, Z);
    return SUN_SUCCESS;
  }

  /* Cases:                     */
  /*   (1) a == -1.0, b != 1.0, */
  /*   (2) a != 1.0, b == -1.0  */
  if ((test = ((creal(a) == -ONE) && (cimag(a)==ZERO))) || ((creal(b) == -ONE) && (cimag(b)==ZERO)))
  {
    c  = test ? b : a;
    V1 = test ? Y : X;
    V2 = test ? X : Y;
    VLin2VectorArray_Mine(nvec, c, V1, V2, Z);
    return SUN_SUCCESS;
  }

  /* Case: a == b                                                         */
  /* catches case both a and b are 0.0 - user should have called N_VConst */
  if (a == b)
  {
    VScaleSumVectorArray_Mine(nvec, a, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Case: a == -b */
  if (a == -b)
  {
    VScaleDiffVectorArray_Mine(nvec, a, X, Y, Z);
    return SUN_SUCCESS;
  }

  /* Do all cases not handled above:                               */
  /*   (1) a == other, b == 0.0 - user should have called N_VScale */
  /*   (2) a == 0.0, b == other - user should have called N_VScale */
  /*   (3) a,b == other, a !=b, a != -b                            */

  /* get vector length */
  N = MYNV_LENGTH(Z[0]);

  /* compute linear sum for each vector pair in vector arrays */
  for (i = 0; i < nvec; i++)
  {
    xd = MYNV_REAL_DATA(X[i]);
    yd = MYNV_REAL_DATA(Y[i]);
    zd = MYNV_REAL_DATA(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = a * xd[j] + b * yd[j]; }
  }

  return SUN_SUCCESS;
}

SUNErrCode N_VScaleVectorArray_Mine(int nvec, suncomplextype* c, N_Vector* X,
                                    N_Vector* Z)
{
  int i;
  sunindextype j, N;
  suncomplextype* xd = NULL;
  suncomplextype* zd = NULL;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* should have called N_VScale */
  if (nvec == 1)
  {
    N_VScale_Mine(c[0], X[0], Z[0]);
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = MYNV_LENGTH(Z[0]);

  /*
   * X[i] *= c[i]
   */
  if (X == Z)
  {
    for (i = 0; i < nvec; i++)
    {
      xd = MYNV_COMPLEX_DATA(X[i]);
      for (j = 0; j < N; j++) { xd[j] *= c[i]; }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i] = c[i] * X[i]
   */
  for (i = 0; i < nvec; i++)
  {
    xd = MYNV_COMPLEX_DATA(X[i]);
    zd = MYNV_COMPLEX_DATA(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = c[i] * xd[j]; }
  }
  return SUN_SUCCESS;
}

SUNErrCode N_VScaleVectorArray_Real(int nvec, sunrealtype* c, N_Vector* X,
                                    N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* xd = NULL;
  sunrealtype* zd = NULL;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* should have called N_VScale */
  if (nvec == 1)
  {
    N_VScale_Mine(c[0], X[0], Z[0]);
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = MYNV_LENGTH(Z[0]);

  /*
   * X[i] *= c[i]
   */
  if (X == Z)
  {
    for (i = 0; i < nvec; i++)
    {
      xd = MYNV_REAL_DATA(X[i]);
      for (j = 0; j < N; j++) { xd[j] *= c[i]; }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i] = c[i] * X[i]
   */
  for (i = 0; i < nvec; i++)
  {
    xd = MYNV_REAL_DATA(X[i]);
    zd = MYNV_REAL_DATA(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = c[i] * xd[j]; }
  }
  return SUN_SUCCESS;
}

SUNErrCode N_VConstVectorArray_Mine(int nvec, suncomplextype c, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  suncomplextype* zd = NULL;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* should have called N_VConst */
  if (nvec == 1)
  {
    N_VConst_Mine((suncomplextype)ONE, Z[0]);
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = MYNV_LENGTH(Z[0]);

  /* set each vector in the vector array to a constant */
  for (i = 0; i < nvec; i++)
  {
    zd = MYNV_COMPLEX_DATA(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = c; }
  }

  return SUN_SUCCESS;
}

SUNErrCode N_VConstVectorArray_Real(int nvec, sunrealtype c, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  sunrealtype* zd = NULL;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* should have called N_VConst */
  if (nvec == 1)
  {
    N_VConst_Mine(ONE, Z[0]);
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = MYNV_LENGTH(Z[0]);

  /* set each vector in the vector array to a constant */
  for (i = 0; i < nvec; i++)
  {
    zd = MYNV_REAL_DATA(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = c; }
  }

  return SUN_SUCCESS;
}

SUNErrCode N_VWrmsNormVectorArray_Mine(int nvec, N_Vector* X, N_Vector* W,
                                       sunrealtype* nrm)
{
  int i;
  sunindextype j, N;
  sunrealtype* wd = NULL;
  sunrealtype* xd = NULL;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* should have called N_VWrmsNorm */
  if (nvec == 1)
  {
    nrm[0] = N_VWrmsNorm_Mine(X[0], W[0]); //this function is not ready yet!
    return SUN_SUCCESS;
  }

  /* get vector length */
  N = MYNV_LENGTH(X[0]);

  printf("SUNErrCode N_VWrmsNormVectorArray_Mine(int nvec, N_Vector* X, N_Vector* W, sunrealtype* nrm) is not ready yet! Check myNvector.c\n");

  /* compute the WRMS norm for each vector in the vector array */
  for (i = 0; i < nvec; i++)
  {
    xd     = MYNV_REAL_DATA(X[i]); //change to complex data
    wd     = MYNV_REAL_DATA(W[i]); //change to complex data
    nrm[i] = ZERO;
    for (j = 0; j < N; j++) { nrm[i] += SUNSQR(xd[j] * wd[j]); }
    nrm[i] = SUNRsqrt(nrm[i] / N);
  }

  return SUN_SUCCESS;
}

SUNErrCode N_VWrmsNormMaskVectorArray_Mine(int nvec, N_Vector* X, N_Vector* W,
                                           N_Vector id, sunrealtype* nrm)
{
  int i;
  sunindextype j, N;
  sunrealtype* wd  = NULL;
  sunrealtype* xd  = NULL;
  sunrealtype* idd = NULL;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* should have called N_VWrmsNorm */
  if (nvec == 1)
  {
    nrm[0] = N_VWrmsNormMask_Mine(X[0], W[0], id);
    return SUN_SUCCESS;
  }

  printf("SUNErrCode N_VWrmsNormMaskVectorArray_Mine(int nvec, N_Vector* X, N_Vector* W, N_Vector id, sunrealtype* nrm) is not ready yet! Check myNvector.c\n");


  /* get vector length and mask data array */
  N   = MYNV_LENGTH(X[0]);
  idd = MYNV_REAL_DATA(id); //change to complex data

  /* compute the WRMS norm for each vector in the vector array */
  for (i = 0; i < nvec; i++)
  {
    xd     = MYNV_REAL_DATA(X[i]); //change to complex data
    wd     = MYNV_REAL_DATA(W[i]); //change to complex data
    nrm[i] = ZERO;
    for (j = 0; j < N; j++)
    {
      if (idd[j] > ZERO) { nrm[i] += SUNSQR(xd[j] * wd[j]); }
    }
    nrm[i] = SUNRsqrt(nrm[i] / N);
  }

  return SUN_SUCCESS;
}

SUNErrCode N_VScaleAddMultiVectorArray_Mine(int nvec, int nsum,
                                            suncomplextype* a, N_Vector* X,
                                            N_Vector** Y, N_Vector** Z)
{
  int i, j;
  sunindextype k, N;
  suncomplextype* xd = NULL;
  suncomplextype* yd = NULL;
  suncomplextype* zd = NULL;
  N_Vector* YY;
  N_Vector* ZZ;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1)
  {
    /* should have called N_VLinearSum */
    if (nsum == 1)
    {
      N_VLinearSum_Mine(a[0], X[0], (suncomplextype)ONE, Y[0][0], Z[0][0]);

      return SUN_SUCCESS;
    }

    /* should have called N_VScaleAddMulti */
    YY = (N_Vector*)malloc(nsum * sizeof(N_Vector));
    if (YY == NULL) { return SUN_ERR_MALLOC_FAIL; }
    ZZ = (N_Vector*)malloc(nsum * sizeof(N_Vector));
    if (ZZ == NULL) { return SUN_ERR_MALLOC_FAIL; }

    for (j = 0; j < nsum; j++)
    {
      YY[j] = Y[j][0];
      ZZ[j] = Z[j][0];
    }

    N_VScaleAddMulti_Mine(nsum, a, X[0], YY, ZZ);

    free(YY);
    free(ZZ);

    return SUN_SUCCESS;
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 1)
  {
    N_VLinearSumVectorArray_Mine(nvec, a[0], X, (suncomplextype)ONE, Y[0], Z[0]);
    return SUN_SUCCESS;
  }

  /* ----------------------------
   * Compute multiple linear sums
   * ---------------------------- */

  /* get vector length */
  N = MYNV_LENGTH(X[0]);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z)
  {
    for (i = 0; i < nvec; i++)
    {
      xd = MYNV_COMPLEX_DATA(X[i]);
      for (j = 0; j < nsum; j++)
      {
        yd = MYNV_COMPLEX_DATA(Y[j][i]);
        for (k = 0; k < N; k++) { yd[k] += a[j] * xd[k]; }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i = 0; i < nvec; i++)
  {
    xd = MYNV_COMPLEX_DATA(X[i]);
    for (j = 0; j < nsum; j++)
    {
      yd = MYNV_COMPLEX_DATA(Y[j][i]);
      zd = MYNV_COMPLEX_DATA(Z[j][i]);
      for (k = 0; k < N; k++) { zd[k] = a[j] * xd[k] + yd[k]; }
    }
  }
  return SUN_SUCCESS;
}

SUNErrCode N_VScaleAddMultiVectorArray_Real(int nvec, int nsum,
                                            sunrealtype* a, N_Vector* X,
                                            N_Vector** Y, N_Vector** Z)
{
  int i, j;
  sunindextype k, N;
  sunrealtype* xd = NULL;
  sunrealtype* yd = NULL;
  sunrealtype* zd = NULL;
  N_Vector* YY;
  N_Vector* ZZ;

  /* invalid number of vectors */
  if (nvec < 0) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1)
  {
    /* should have called N_VLinearSum */
    if (nsum == 1)
    {
      N_VLinearSum_Mine(a[0], X[0], ONE, Y[0][0], Z[0][0]);

      return SUN_SUCCESS;
    }

    /* should have called N_VScaleAddMulti */
    YY = (N_Vector*)malloc(nsum * sizeof(N_Vector));
    if (YY == NULL) { return SUN_ERR_MALLOC_FAIL; }
    ZZ = (N_Vector*)malloc(nsum * sizeof(N_Vector));
    if (ZZ == NULL) { return SUN_ERR_MALLOC_FAIL; }

    for (j = 0; j < nsum; j++)
    {
      YY[j] = Y[j][0];
      ZZ[j] = Z[j][0];
    }

    N_VScaleAddMulti_Mine(nsum, (suncomplextype*)a, X[0], YY, ZZ);

    free(YY);
    free(ZZ);

    return SUN_SUCCESS;
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 1)
  {
    N_VLinearSumVectorArray_Mine(nvec, a[0], X, ONE, Y[0], Z[0]);
    return SUN_SUCCESS;
  }

  /* ----------------------------
   * Compute multiple linear sums
   * ---------------------------- */

  /* get vector length */
  N = MYNV_LENGTH(X[0]);

  /*
   * Y[i][j] += a[i] * x[j]
   */
  if (Y == Z)
  {
    for (i = 0; i < nvec; i++)
    {
      xd = MYNV_REAL_DATA(X[i]);
      for (j = 0; j < nsum; j++)
      {
        yd = MYNV_REAL_DATA(Y[j][i]);
        for (k = 0; k < N; k++) { yd[k] += a[j] * xd[k]; }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[i][j] = Y[i][j] + a[i] * x[j]
   */
  for (i = 0; i < nvec; i++)
  {
    xd = MYNV_REAL_DATA(X[i]);
    for (j = 0; j < nsum; j++)
    {
      yd = MYNV_REAL_DATA(Y[j][i]);
      zd = MYNV_REAL_DATA(Z[j][i]);
      for (k = 0; k < N; k++) { zd[k] = a[j] * xd[k] + yd[k]; }
    }
  }
  return SUN_SUCCESS;
}























SUNErrCode N_VLinearCombinationVectorArray_Mine(int nvec, int nsum,
                                                suncomplextype* c, N_Vector** X,
                                                N_Vector* Z)
{
  int i;          /* vector arrays index in summation [0,nsum) */
  int j;          /* vector index in vector array     [0,nvec) */
  sunindextype k; /* element index in vector          [0,N)    */
  sunindextype N;
  suncomplextype* zd = NULL;
  suncomplextype* xd = NULL;
  suncomplextype* ctmp;
  N_Vector* Y;

  /* invalid number of vectors */
  if ((nvec < 1) || (nsum < 1)) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1)
  {
    /* should have called N_VScale */
    if (nsum == 1)
    {
      N_VScale_Mine(c[0], X[0][0], Z[0]);
      return SUN_SUCCESS;
    }

    /* should have called N_VLinearSum */
    if (nsum == 2)
    {
      N_VLinearSum_Mine(c[0], X[0][0], c[1], X[1][0], Z[0]);
      return SUN_SUCCESS;
    }

    /* should have called N_VLinearCombination */
    Y = (N_Vector*)malloc(nsum * sizeof(N_Vector));
    if (Y == NULL) { return SUN_ERR_MALLOC_FAIL; }

    for (i = 0; i < nsum; i++) { Y[i] = X[i][0]; }

    N_VLinearCombination_Mine(nsum, c, Y, Z[0]);

    free(Y);

    return SUN_SUCCESS;
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VScaleVectorArray */
  if (nsum == 1)
  {
    ctmp = (suncomplextype*)malloc(nvec * sizeof(suncomplextype));
    if (ctmp == NULL) { return SUN_ERR_MALLOC_FAIL; }

    for (j = 0; j < nvec; j++) { ctmp[j] = c[0]; }

    N_VScaleVectorArray_Mine(nvec, ctmp, X[0], Z);

    free(ctmp);
    return SUN_SUCCESS;
  }

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 2)
  {
    N_VLinearSumVectorArray_Mine(nvec, c[0], X[0], c[1], X[1], Z);
    return SUN_SUCCESS;
  }

  /* --------------------------
   * Compute linear combination
   * -------------------------- */

  /* get vector length */
  N = MYNV_LENGTH(Z[0]);

  /*
   * X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
   */
  if ((X[0] == Z) && ((creal(c[0]) == ONE) && cimag(c[0]) == ZERO))
  {
    for (j = 0; j < nvec; j++)
    {
      zd = MYNV_COMPLEX_DATA(Z[j]);
      for (i = 1; i < nsum; i++)
      {
        xd = MYNV_COMPLEX_DATA(X[i][j]);
        for (k = 0; k < N; k++) { zd[k] += c[i] * xd[k]; }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * X[0][j] = c[0] * X[0][j] + sum{ c[i] * X[i][j] }, i = 1,...,nvec-1
   */
  if (X[0] == Z)
  {
    for (j = 0; j < nvec; j++)
    {
      zd = MYNV_COMPLEX_DATA(Z[j]);
      for (k = 0; k < N; k++) { zd[k] *= c[0]; }
      for (i = 1; i < nsum; i++)
      {
        xd = MYNV_COMPLEX_DATA(X[i][j]);
        for (k = 0; k < N; k++) { zd[k] += c[i] * xd[k]; }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[j] = sum{ c[i] * X[i][j] }, i = 0,...,nvec-1
   */
  for (j = 0; j < nvec; j++)
  {
    xd = MYNV_COMPLEX_DATA(X[0][j]);
    zd = MYNV_COMPLEX_DATA(Z[j]);
    for (k = 0; k < N; k++) { zd[k] = c[0] * xd[k]; }
    for (i = 1; i < nsum; i++)
    {
      xd = MYNV_COMPLEX_DATA(X[i][j]);
      for (k = 0; k < N; k++) { zd[k] += c[i] * xd[k]; }
    }
  }
  return SUN_SUCCESS;
}

SUNErrCode N_VLinearCombinationVectorArray_Real(int nvec, int nsum,
                                                sunrealtype* c, N_Vector** X,
                                                N_Vector* Z)
{
  int i;          /* vector arrays index in summation [0,nsum) */
  int j;          /* vector index in vector array     [0,nvec) */
  sunindextype k; /* element index in vector          [0,N)    */
  sunindextype N;
  sunrealtype* zd = NULL;
  sunrealtype* xd = NULL;
  sunrealtype* ctmp;
  N_Vector* Y;

  /* invalid number of vectors */
  if ((nvec < 1) || (nsum < 1)) { return SUN_ERR_ARG_OUTOFRANGE; }

  /* ---------------------------
   * Special cases for nvec == 1
   * --------------------------- */

  if (nvec == 1)
  {
    /* should have called N_VScale */
    if (nsum == 1)
    {
      N_VScale_Mine(c[0], X[0][0], Z[0]);
      return SUN_SUCCESS;
    }

    /* should have called N_VLinearSum */
    if (nsum == 2)
    {
      N_VLinearSum_Mine(c[0], X[0][0], c[1], X[1][0], Z[0]);
      return SUN_SUCCESS;
    }

    /* should have called N_VLinearCombination */
    Y = (N_Vector*)malloc(nsum * sizeof(N_Vector));
    if (Y == NULL) { return SUN_ERR_MALLOC_FAIL; }

    for (i = 0; i < nsum; i++) { Y[i] = X[i][0]; }

    N_VLinearCombination_Mine(nsum, (suncomplextype*)c, Y, Z[0]);

    free(Y);

    return SUN_SUCCESS;
  }

  /* --------------------------
   * Special cases for nvec > 1
   * -------------------------- */

  /* should have called N_VScaleVectorArray */
  if (nsum == 1)
  {
    ctmp = (sunrealtype*)malloc(nvec * sizeof(sunrealtype));
    if (ctmp == NULL) { return SUN_ERR_MALLOC_FAIL; }

    for (j = 0; j < nvec; j++) { ctmp[j] = c[0]; }

    N_VScaleVectorArray_Mine(nvec, (suncomplextype*)ctmp, X[0], Z);

    free(ctmp);
    return SUN_SUCCESS;
  }

  /* should have called N_VLinearSumVectorArray */
  if (nsum == 2)
  {
    N_VLinearSumVectorArray_Mine(nvec, c[0], X[0], c[1], X[1], Z);
    return SUN_SUCCESS;
  }

  /* --------------------------
   * Compute linear combination
   * -------------------------- */

  /* get vector length */
  N = MYNV_LENGTH(Z[0]);

  /*
   * X[0][j] += c[i]*X[i][j], i = 1,...,nvec-1
   */
  if ((X[0] == Z) && ((creal(c[0]) == ONE) && cimag(c[0]) == ZERO))
  {
    for (j = 0; j < nvec; j++)
    {
      zd = MYNV_REAL_DATA(Z[j]);
      for (i = 1; i < nsum; i++)
      {
        xd = MYNV_REAL_DATA(X[i][j]);
        for (k = 0; k < N; k++) { zd[k] += c[i] * xd[k]; }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * X[0][j] = c[0] * X[0][j] + sum{ c[i] * X[i][j] }, i = 1,...,nvec-1
   */
  if (X[0] == Z)
  {
    for (j = 0; j < nvec; j++)
    {
      zd = MYNV_REAL_DATA(Z[j]);
      for (k = 0; k < N; k++) { zd[k] *= c[0]; }
      for (i = 1; i < nsum; i++)
      {
        xd = MYNV_REAL_DATA(X[i][j]);
        for (k = 0; k < N; k++) { zd[k] += c[i] * xd[k]; }
      }
    }
    return SUN_SUCCESS;
  }

  /*
   * Z[j] = sum{ c[i] * X[i][j] }, i = 0,...,nvec-1
   */
  for (j = 0; j < nvec; j++)
  {
    xd = MYNV_REAL_DATA(X[0][j]);
    zd = MYNV_REAL_DATA(Z[j]);
    for (k = 0; k < N; k++) { zd[k] = c[0] * xd[k]; }
    for (i = 1; i < nsum; i++)
    {
      xd = MYNV_REAL_DATA(X[i][j]);
      for (k = 0; k < N; k++) { zd[k] += c[i] * xd[k]; }
    }
  }
  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * OPTIONAL XBraid interface operations
 * -----------------------------------------------------------------
 */

SUNErrCode N_VBufSize_Mine(N_Vector x, sunindextype* size)
{
  *size = MYNV_LENGTH(x) * ((sunindextype)sizeof(suncomplextype));
  return SUN_SUCCESS;
}

SUNErrCode N_VBufPack_Mine(N_Vector x, void* buf)
{
  sunindextype i, N;
  suncomplextype* xd = NULL;
  suncomplextype* bd = NULL;

  if (buf == NULL) { return SUN_ERR_ARG_CORRUPT; }
  
    printf("SUNErrCode N_VBufPack_Mine(N_Vector x, void* buf) is not ready yet! Check myNvector.c\n");

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  bd = (suncomplextype*)buf;

  for (i = 0; i < N; i++) { bd[i] = xd[i]; }

  return SUN_SUCCESS;
}

SUNErrCode N_VBufUnpack_Mine(N_Vector x, void* buf)
{
  sunindextype i, N;
  suncomplextype* xd = NULL;
  suncomplextype* bd = NULL;

  printf("SUNErrCode N_VBufUnpack_Mine(N_Vector x, void* buf) is not ready yet! Check myNvector.c\n");

  if (buf == NULL) { return SUN_ERR_ARG_CORRUPT; }

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  bd = (suncomplextype*)buf;

  for (i = 0; i < N; i++) { xd[i] = bd[i]; }

  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector operations
 * -----------------------------------------------------------------
 */

static void VCopy_Mine(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  suncomplextype *xd, *zd;

  xd = zd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  zd = MYNV_COMPLEX_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i]; }

  return;
}

static void VSum_Mine(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  suncomplextype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  yd = MYNV_COMPLEX_DATA(y);
  zd = MYNV_COMPLEX_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] + yd[i]; }

  return;
}

static void VDiff_Mine(N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  suncomplextype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  yd = MYNV_COMPLEX_DATA(y);
  zd = MYNV_COMPLEX_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = xd[i] - yd[i]; }

  return;
}

static void VNeg_Mine(N_Vector x, N_Vector z)
{
  sunindextype i, N;
  suncomplextype *xd, *zd;

  xd = zd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  zd = MYNV_COMPLEX_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = -xd[i]; }

  return;
}

static void VScaleSum_Mine(suncomplextype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  suncomplextype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  yd = MYNV_COMPLEX_DATA(y);
  zd = MYNV_COMPLEX_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = c * (xd[i] + yd[i]); }

  return;
}

static void VScaleDiff_Mine(suncomplextype c, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  suncomplextype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  yd = MYNV_COMPLEX_DATA(y);
  zd = MYNV_COMPLEX_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = c * (xd[i] - yd[i]); }

  return;
}

static void VLin1_Mine(suncomplextype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  suncomplextype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  yd = MYNV_COMPLEX_DATA(y);
  zd = MYNV_COMPLEX_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = (a * xd[i]) + yd[i]; }

  return;
}

static void VLin2_Mine(suncomplextype a, N_Vector x, N_Vector y, N_Vector z)
{
  sunindextype i, N;
  suncomplextype *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  yd = MYNV_COMPLEX_DATA(y);
  zd = MYNV_COMPLEX_DATA(z);

  for (i = 0; i < N; i++) { zd[i] = (a * xd[i]) - yd[i]; }

  return;
}

static void Vaxpy_Mine(suncomplextype a, N_Vector x, N_Vector y)
{
  sunindextype i, N;
  suncomplextype *xd, *yd;

  xd = yd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);
  yd = MYNV_COMPLEX_DATA(y);

  if ((creal(a) == ONE) && (cimag(a) == ZERO))
  {
    for (i = 0; i < N; i++) { yd[i] += xd[i]; }
    return;
  }

  if ((creal(a) == -ONE) && (cimag(a) == ZERO))
  {
    for (i = 0; i < N; i++) { yd[i] -= xd[i]; }
    return;
  }

  for (i = 0; i < N; i++) { yd[i] += a * xd[i]; }

  return;
}

static void VScaleBy_Mine(suncomplextype a, N_Vector x)
{
  sunindextype i, N;
  suncomplextype* xd;

  xd = NULL;

  N  = MYNV_LENGTH(x);
  xd = MYNV_COMPLEX_DATA(x);

  for (i = 0; i < N; i++) { xd[i] *= a; }

  return;
}

/*
 * -----------------------------------------------------------------
 * private functions for special cases of vector array operations
 * -----------------------------------------------------------------
 */

static void VSumVectorArray_Mine(int nvec, N_Vector* X, N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  suncomplextype* xd = NULL;
  suncomplextype* yd = NULL;
  suncomplextype* zd = NULL;

  N = MYNV_LENGTH(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = MYNV_COMPLEX_DATA(X[i]);
    yd = MYNV_COMPLEX_DATA(Y[i]);
    zd = MYNV_COMPLEX_DATA(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = xd[j] + yd[j]; }
  }
}

static void VDiffVectorArray_Mine(int nvec, N_Vector* X, N_Vector* Y,
                                  N_Vector* Z)
{
  int i;
  sunindextype j, N;
  suncomplextype* xd = NULL;
  suncomplextype* yd = NULL;
  suncomplextype* zd = NULL;

  N = MYNV_LENGTH(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = MYNV_COMPLEX_DATA(X[i]);
    yd = MYNV_COMPLEX_DATA(Y[i]);
    zd = MYNV_COMPLEX_DATA(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = xd[j] - yd[j]; }
  }
}

static void VScaleSumVectorArray_Mine(int nvec, suncomplextype c, N_Vector* X,
                                      N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  suncomplextype* xd = NULL;
  suncomplextype* yd = NULL;
  suncomplextype* zd = NULL;

  N = MYNV_LENGTH(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = MYNV_COMPLEX_DATA(X[i]);
    yd = MYNV_COMPLEX_DATA(Y[i]);
    zd = MYNV_COMPLEX_DATA(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = c * (xd[j] + yd[j]); }
  }
}

static void VScaleDiffVectorArray_Mine(int nvec, suncomplextype c, N_Vector* X,
                                       N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  suncomplextype* xd = NULL;
  suncomplextype* yd = NULL;
  suncomplextype* zd = NULL;

  N = MYNV_LENGTH(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = MYNV_COMPLEX_DATA(X[i]);
    yd = MYNV_COMPLEX_DATA(Y[i]);
    zd = MYNV_COMPLEX_DATA(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = c * (xd[j] - yd[j]); }
  }
}

static void VLin1VectorArray_Mine(int nvec, suncomplextype a, N_Vector* X,
                                  N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  suncomplextype* xd = NULL;
  suncomplextype* yd = NULL;
  suncomplextype* zd = NULL;

  N = MYNV_LENGTH(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = MYNV_COMPLEX_DATA(X[i]);
    yd = MYNV_COMPLEX_DATA(Y[i]);
    zd = MYNV_COMPLEX_DATA(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = (a * xd[j]) + yd[j]; }
  }
}

static void VLin2VectorArray_Mine(int nvec, suncomplextype a, N_Vector* X,
                                  N_Vector* Y, N_Vector* Z)
{
  int i;
  sunindextype j, N;
  suncomplextype* xd = NULL;
  suncomplextype* yd = NULL;
  suncomplextype* zd = NULL;

  N = MYNV_LENGTH(X[0]);

  for (i = 0; i < nvec; i++)
  {
    xd = MYNV_COMPLEX_DATA(X[i]);
    yd = MYNV_COMPLEX_DATA(Y[i]);
    zd = MYNV_COMPLEX_DATA(Z[i]);
    for (j = 0; j < N; j++) { zd[j] = (a * xd[j]) - yd[j]; }
  }
}

static void VaxpyVectorArray_Mine(int nvec, suncomplextype a, N_Vector* X,
                                  N_Vector* Y)
{
  int i;
  sunindextype j, N;
  suncomplextype* xd = NULL;
  suncomplextype* yd = NULL;

  N = MYNV_LENGTH(X[0]);

  if ((creal(a) == ONE) && (cimag(a) == ZERO))
  {
    for (i = 0; i < nvec; i++)
    {
      xd = MYNV_COMPLEX_DATA(X[i]);
      yd = MYNV_COMPLEX_DATA(Y[i]);
      for (j = 0; j < N; j++) { yd[j] += xd[j]; }
    }
    return;
  }

  if ((creal(a) == -ONE) && (cimag(a) == ZERO))
  {
    for (i = 0; i < nvec; i++)
    {
      xd = MYNV_COMPLEX_DATA(X[i]);
      yd = MYNV_COMPLEX_DATA(Y[i]);
      for (j = 0; j < N; j++) { yd[j] -= xd[j]; }
    }
    return;
  }

  for (i = 0; i < nvec; i++)
  {
      xd = MYNV_COMPLEX_DATA(X[i]);
      yd = MYNV_COMPLEX_DATA(Y[i]);
    for (j = 0; j < N; j++) { yd[j] += a * xd[j]; }
  }
}

/*
 * -----------------------------------------------------------------
 * Enable / Disable fused and vector array operations
 * -----------------------------------------------------------------
 */

SUNErrCode N_VEnableFusedOps_Mine(N_Vector v, sunbooleantype tf)
{
  if (tf)
  {
    /* enable all fused vector operations */
    // v->ops->nvlinearcombination = N_VLinearCombination_Mine;
    v->ops->nvlinearcombination = N_VLinearCombination_Real;

    // v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Mine;
    v->ops->nvscaleaddmulti     = N_VScaleAddMulti_Real;
    
    // v->ops->nvdotprodmulti      = N_VDotProdMulti_Mine;
    v->ops->nvdotprodmulti      = N_VDotProdMulti_Real;


    /* enable all vector array operations */
    // v->ops->nvlinearsumvectorarray     = N_VLinearSumVectorArray_Mine;
    v->ops->nvlinearsumvectorarray     = N_VLinearSumVectorArray_Real;

    // v->ops->nvscalevectorarray         = N_VScaleVectorArray_Mine;
    v->ops->nvscalevectorarray         = N_VScaleVectorArray_Real;


    // v->ops->nvconstvectorarray         = N_VConstVectorArray_Mine;
    v->ops->nvconstvectorarray         = N_VConstVectorArray_Real;

    v->ops->nvwrmsnormvectorarray      = N_VWrmsNormVectorArray_Mine;
    v->ops->nvwrmsnormmaskvectorarray  = N_VWrmsNormMaskVectorArray_Mine;
    
    // v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Mine;
    v->ops->nvscaleaddmultivectorarray = N_VScaleAddMultiVectorArray_Real;
    
    // v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Mine;
    v->ops->nvlinearcombinationvectorarray = N_VLinearCombinationVectorArray_Real;

    /* enable single buffer reduction operations */
    
    // v->ops->nvdotprodmultilocal = N_VDotProdMulti_Mine;
    v->ops->nvdotprodmultilocal = N_VDotProdMulti_Real;

  }
  else
  {
    /* disable all fused vector operations */
    v->ops->nvlinearcombination = NULL;
    v->ops->nvscaleaddmulti     = NULL;
    v->ops->nvdotprodmulti      = NULL;
    /* disable all vector array operations */
    v->ops->nvlinearsumvectorarray         = NULL;
    v->ops->nvscalevectorarray             = NULL;
    v->ops->nvconstvectorarray             = NULL;
    v->ops->nvwrmsnormvectorarray          = NULL;
    v->ops->nvwrmsnormmaskvectorarray      = NULL;
    v->ops->nvscaleaddmultivectorarray     = NULL;
    v->ops->nvlinearcombinationvectorarray = NULL;
    /* disable single buffer reduction operations */
    v->ops->nvdotprodmultilocal = NULL;
  }

  /* return success */
  return SUN_SUCCESS;
}