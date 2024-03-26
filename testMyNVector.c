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
 * These test functions check some components of a custom NVECTOR
 * module implementation (for more thorough tests, see the main
 * SUNDIALS repository, inside examples/nvector/).
 * -----------------------------------------------------------------*/

#include <sundials/sundials_config.h>
#include "myNVector.h"
#include <stdio.h>
#include <stdlib.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_nvector.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM ".17f"
#endif

#define ZERO      SUN_RCONST(0.0)
#define HALF      SUN_RCONST(0.5)
#define ONE       SUN_RCONST(1.0)
#define TWO       SUN_RCONST(2.0)
#define NEG_HALF  SUN_RCONST(-0.5)
#define NEG_ONE   SUN_RCONST(-1.0)
#define NEG_TWO   SUN_RCONST(-2.0)

/* ----------------------------------------------------------------------
 * Implementation specific utility functions for vector tests
 * --------------------------------------------------------------------*/
int check_ans(sunrealtype ans, N_Vector X, sunindextype local_length)
{
  int failure = 0;
  sunindextype i;
  sunrealtype* Xdata;

  Xdata = N_VGetArrayPointer(X);

  /* check vector data */
  for (i = 0; i < local_length; i++) { failure += SUNRCompare(Xdata[i], ans); }

  return (failure > ZERO) ? (1) : (0);
}

sunbooleantype has_data(N_Vector X)
{
  /* check if data array is non-null */
  return (N_VGetArrayPointer(X) == NULL) ? SUNFALSE : SUNTRUE;
}

void set_element_range(N_Vector X, sunindextype is, sunindextype ie,
                       sunrealtype val)
{
  sunindextype i;

  /* set elements [is,ie] of the data array */
  sunrealtype* xd = N_VGetArrayPointer(X);
  for (i = is; i <= ie; i++) { xd[i] = val; }
}

void set_element(N_Vector X, sunindextype i, sunrealtype val)
{
  /* set i-th element of data array */
  set_element_range(X, i, i, val);
}

sunrealtype get_element(N_Vector X, sunindextype i)
{
  /* get i-th element of data array */
  return MYNV_Ith(X, i);
}



/* ----------------------------------------------------------------------
 * Main NVector Testing Routine
 * --------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
  int fails = 0;             /* counter for test failures */
  int retval;                /* function return value     */
  sunindextype length;       /* vector length             */
  N_Vector V, W, X, Y, Z;    /* test vectors              */
  SUNContext sunctx = NULL;

  if (SUNContext_Create(SUN_COMM_NULL, &sunctx))
  {
    printf("ERROR: SUNContext_Create failed\n");
    return -1;
  }

  /* check input and set vector length */
  if (argc < 2)
  {
    printf("ERROR: ONE (1) Input required: vector length \n");
    return (-1);
  }

  length = (sunindextype)atol(argv[1]);
  if (length <= 0)
  {
    printf("ERROR: length of vector must be a positive integer \n");
    return (-1);
  }

  printf("Testing custom N_Vector \n");
  printf("Vector length %ld \n", (long int)length);

  /* Create new vectors */
  V = N_VNewEmpty_Mine(length, sunctx);
  if (V == NULL)
  {
    printf("FAIL: Unable to create a new empty vector \n\n");
    return (1);
  }
  N_VDestroy(V);

  X = N_VNew_Mine(length, sunctx);
  if (X == NULL)
  {
    printf("FAIL: Unable to create a new vector \n\n");
    return (1);
  }

  /* Check vector ID */
  if (N_VGetVectorID(X) != SUNDIALS_NVEC_CUSTOM)
  {
    printf(">>> FAILED test -- N_VGetVectorID \n");
    printf("    Unrecognized vector type %d \n \n", N_VGetVectorID(X));
    fails += 1;
  }
  else { printf("PASSED test -- N_VGetVectorID \n"); }

  /* Check vector length */
  sunindextype Xlength = N_VGetLength(X);
  N_VConst(SUN_RCONST(1.0), X);
  sunindextype Xlength2 = (sunindextype)N_VDotProd(X, X);
  if (Xlength != Xlength2)
  {
    printf(">>> FAILED test -- N_VGetLength (%li != %li)\n",
           (long int)Xlength, (long int)Xlength2);
    fails += 1;
  }
  else { printf("PASSED test -- N_VGetLength\n"); }

  /* Test N_VClone */
  W = N_VClone(X);

  /*   check cloned vector */
  if (W == NULL)
  {
    printf(">>> FAILED test -- N_VClone \n");
    printf("    After N_VClone, X == NULL \n\n");
    fails += 1;
  }

  /*   check cloned vector data */
  if (!has_data(W))
  {
    printf(">>> FAILED test -- N_VClone \n");
    printf("    Vector data == NULL \n\n");
    N_VDestroy(W);
    fails += 1;
  }
  else
  {
    N_VConst(ONE, W);
    if (check_ans(ONE, W, length))
    {
      printf(">>> FAILED test -- N_VClone \n");
      printf("    Failed N_VClone check \n\n");
      N_VDestroy(W);
      fails += 1;
    }
    else
    { printf("PASSED test -- N_VClone \n"); }
    N_VDestroy(W);
  }

  /* Clone additional vectors for testing */
  Y = N_VClone(X);
  if (Y == NULL)
  {
    N_VDestroy(W);
    N_VDestroy(X);
    printf("FAIL: Unable to create a new vector \n\n");
    return (1);
  }

  Z = N_VClone(X);
  if (Z == NULL)
  {
    N_VDestroy(W);
    N_VDestroy(X);
    N_VDestroy(Y);
    printf("FAIL: Unable to create a new vector \n\n");
    return (1);
  }

  /* Standard vector operation tests */
  printf("\nTesting standard vector operations:\n\n");

  /* Test N_VConst: fill vector data with zeros to prevent passing
     in the case where the input vector is a vector of ones */
  set_element_range(X, 0, length - 1, ZERO);
  N_VConst(ONE, X);
  if (check_ans(ONE, X, length))
  {
    printf(">>> FAILED test -- N_VConst \n");
    fails += 1;
  }
  else { printf("PASSED test -- N_VConst \n"); }

  /* Test N_VLinearSum */

  /*   Case 1a: y = x + y, (Vaxpy Case 1) */
  N_VConst(ONE, X);
  N_VConst(NEG_TWO, Y);
  N_VLinearSum(ONE, X, ONE, Y, Y);
  if (check_ans(NEG_ONE, Y, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 1a \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 1a \n"); }

  /*   Case 1b: y = -x + y, (Vaxpy Case 2) */
  N_VConst(ONE, X);
  N_VConst(TWO, Y);
  N_VLinearSum(NEG_ONE, X, ONE, Y, Y);
  if (check_ans(ONE, Y, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 1b \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 1b \n"); }

  /*   Case 1c: y = ax + y, (Vaxpy Case 3) */
  N_VConst(TWO, X);
  N_VConst(NEG_TWO, Y);
  N_VLinearSum(HALF, X, ONE, Y, Y);
  if (check_ans(NEG_ONE, Y, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 1c \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 1c \n"); }

  /*   Case 2a: x = x + y, (Vaxpy Case 1) */
  N_VConst(TWO, X);
  N_VConst(NEG_ONE, Y);
  N_VLinearSum(ONE, X, ONE, Y, X);
  if (check_ans(ONE, X, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 2a \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 2a \n"); }

  /*   Case 2b: x = x - y, (Vaxpy Case 2)*/
  N_VConst(ONE, X);
  N_VConst(TWO, Y);
  N_VLinearSum(ONE, X, NEG_ONE, Y, X);
  if (check_ans(NEG_ONE, X, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 2b \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 2b \n"); }

  /*   Case 2c: x = x + by, (Vaxpy Case 3) */
  N_VConst(TWO, X);
  N_VConst(NEG_HALF, Y);
  N_VLinearSum(ONE, X, TWO, Y, X);
  if (check_ans(ONE, X, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 2c \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 2c \n"); }

  /*   Case 3: z = x + y, (VSum) */
  N_VConst(NEG_TWO, X);
  N_VConst(ONE, Y);
  N_VConst(ZERO, Z);
  N_VLinearSum(ONE, X, ONE, Y, Z);
  if (check_ans(NEG_ONE, Z, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 3 \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 3 \n"); }

  /*   Case 4a: z = x - y, (VDiff) */
  N_VConst(TWO, X);
  N_VConst(ONE, Y);
  N_VConst(ZERO, Z);
  N_VLinearSum(ONE, X, NEG_ONE, Y, Z);
  if (check_ans(ONE, Z, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 4a \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 4a \n"); }

  /*   Case 4b: z = -x + y, (VDiff) */
  N_VConst(TWO, X);
  N_VConst(ONE, Y);
  N_VConst(ZERO, Z);
  N_VLinearSum(NEG_ONE, X, ONE, Y, Z);
  if (check_ans(NEG_ONE, Z, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 4b \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 4b \n"); }

  /*   Case 5a: z = x + by, (VLin1) */
  N_VConst(TWO, X);
  N_VConst(NEG_HALF, Y);
  N_VConst(ZERO, Z);
  N_VLinearSum(ONE, X, TWO, Y, Z);
  if (check_ans(ONE, Z, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 5a \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 5a \n"); }

  /*   Case 5b: z = ax + y, (VLin1) */
  N_VConst(HALF, X);
  N_VConst(NEG_TWO, Y);
  N_VConst(ZERO, Z);
  N_VLinearSum(TWO, X, ONE, Y, Z);
  if (check_ans(NEG_ONE, Z, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 5b \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 5b \n"); }

  /*   Case 6a: z = -x + by, (VLin2) */
  N_VConst(NEG_TWO, X);
  N_VConst(NEG_HALF, Y);
  N_VConst(ZERO, Z);
  N_VLinearSum(NEG_ONE, X, TWO, Y, Z);
  if (check_ans(ONE, Z, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 6a \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 6a \n"); }

  /*   Case 6b: z = ax - y, (VLin2) */
  N_VConst(HALF, X);
  N_VConst(TWO, Y);
  N_VConst(ZERO, Z);
  N_VLinearSum(TWO, X, NEG_ONE, Y, Z);
  if (check_ans(NEG_ONE, Z, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 6b \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 6b \n"); }

  /*   Case 7: z = a(x + y), (VScaleSum) */
  N_VConst(ONE, X);
  N_VConst(NEG_HALF, Y);
  N_VConst(ZERO, Z);
  N_VLinearSum(TWO, X, TWO, Y, Z);
  if (check_ans(ONE, Z, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 7 \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 7 \n"); }

  /*   Case 8: z = a(x - y), (VScaleDiff) */
  N_VConst(HALF, X);
  N_VConst(ONE, Y);
  N_VConst(ZERO, Z);
  N_VLinearSum(TWO, X, NEG_TWO, Y, Z);
  if (check_ans(NEG_ONE, Z, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 8 \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 8 \n"); }

  /*   Case 9: z = ax + by, All Other Cases */
  N_VConst(ONE, X);
  N_VConst(NEG_TWO, Y);
  N_VConst(ZERO, Z);
  N_VLinearSum(TWO, X, HALF, Y, Z);
  if (check_ans(ONE, Z, length))
  {
    printf(">>> FAILED test -- N_VLinearSum Case 9 \n");
    fails++;
  }
  else { printf("PASSED test -- N_VLinearSum Case 9 \n"); }


  /* Test N_VScale */

  /*   Case 1: x = cx, VScaleBy */
  N_VConst(HALF, X);
  N_VScale(TWO, X, X);
  if (check_ans(ONE, X, length))
  {
    printf(">>> FAILED test -- N_VScale Case 1 \n");
    fails++;
  }
  else { printf("PASSED test -- N_VScale Case 1 \n"); }

  /*   Case 2: z = x, VCopy */
  N_VConst(NEG_ONE, X);
  N_VConst(ZERO, Z);
  N_VScale(ONE, X, Z);
  if (check_ans(NEG_ONE, Z, length))
  {
    printf(">>> FAILED test -- N_VScale Case 2\n");
    fails++;
  }
  else { printf("PASSED test -- N_VScale Case 2 \n"); }

  /*   Case 3: z = -x, VNeg */
  N_VConst(NEG_ONE, X);
  N_VConst(ZERO, Z);
  N_VScale(NEG_ONE, X, Z);
  if (check_ans(ONE, Z, length))
  {
    printf(">>> FAILED test -- N_VScale Case 3 \n");
    fails++;
  }
  else { printf("PASSED test -- N_VScale Case 3 \n"); }

  /*   Case 4: z = cx, All other cases */
  N_VConst(NEG_HALF, X);
  N_VConst(ZERO, Z);
  N_VScale(TWO, X, Z);
  if (check_ans(NEG_ONE, Z, length))
  {
    printf(">>> FAILED test -- N_VScale Case 4 \n");
    fails++;
  }
  else { printf("PASSED test -- N_VScale Case 4 \n"); }

  /* Test N_VDotProd */
  N_VConst(TWO, X);
  N_VConst(HALF, Y);
  sunindextype global_length = N_VGetLength(X);
  sunindextype ans = N_VDotProd(X, Y);

  /* ans should equal global vector length */
  if (SUNRCompare(ans, (sunrealtype)global_length))
  {
    printf(">>> FAILED test -- N_VDotProd \n");
    fails++;
  }
  else { printf("PASSED test -- N_VDotProd \n"); }

  /* Print result */
  if (fails) { printf("FAIL: NVector module failed %i tests \n\n", fails); }
  else { printf("SUCCESS: NVector module passed all tests \n\n"); }

  return (fails);
}
