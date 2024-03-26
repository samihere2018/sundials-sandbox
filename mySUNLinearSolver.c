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
 * This is the implementation file for a custom SUNLinearSolver
 * (note that this is just SUNDIALS' SPTFQMR solver).
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include "mySUNLinearSolver.h"

#define ZERO SUN_RCONST(0.0)
#define ONE  SUN_RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

#define MY_CONTENT(S) ((MySUNLinearSolverContent)(S->content))
#define LASTFLAG(S)   (MY_CONTENT(S)->last_flag)

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new linear solver
 */

SUNLinearSolver SUNLinSol_Mine(N_Vector y, int pretype, int maxl,
                               SUNContext sunctx)
{
  SUNLinearSolver S;
  MySUNLinearSolverContent content;

  /* check for legal pretype and maxl values; if illegal use defaults */
  if ((pretype != SUN_PREC_NONE) && (pretype != SUN_PREC_LEFT) &&
      (pretype != SUN_PREC_RIGHT) && (pretype != SUN_PREC_BOTH))
  {
    pretype = SUN_PREC_NONE;
  }
  if (maxl <= 0) { maxl = MY_MAXL_DEFAULT; }

  /* Create linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty(sunctx);
  if (S == NULL) { return NULL; }

  /* Attach operations */
  S->ops->gettype           = SUNLinSolGetType_Mine;
  S->ops->getid             = SUNLinSolGetID_Mine;
  S->ops->setatimes         = SUNLinSolSetATimes_Mine;
  S->ops->setpreconditioner = SUNLinSolSetPreconditioner_Mine;
  S->ops->setscalingvectors = SUNLinSolSetScalingVectors_Mine;
  S->ops->setzeroguess      = SUNLinSolSetZeroGuess_Mine;
  S->ops->initialize        = SUNLinSolInitialize_Mine;
  S->ops->setup             = SUNLinSolSetup_Mine;
  S->ops->solve             = SUNLinSolSolve_Mine;
  S->ops->numiters          = SUNLinSolNumIters_Mine;
  S->ops->resnorm           = SUNLinSolResNorm_Mine;
  S->ops->resid             = SUNLinSolResid_Mine;
  S->ops->lastflag          = SUNLinSolLastFlag_Mine;
  S->ops->space             = SUNLinSolSpace_Mine;
  S->ops->free              = SUNLinSolFree_Mine;

  /* Create content */
  content = NULL;
  content = (MySUNLinearSolverContent)malloc(sizeof *content);
  if (content == NULL) { return NULL; }

  /* Attach content */
  S->content = content;

  /* Fill content */
  content->last_flag = 0;
  content->maxl      = maxl;
  content->pretype   = pretype;
  content->zeroguess = SUNFALSE;
  content->numiters  = 0;
  content->resnorm   = ZERO;
  content->r_star    = NULL;
  content->q         = NULL;
  content->d         = NULL;
  content->v         = NULL;
  content->p         = NULL;
  content->r         = NULL;
  content->u         = NULL;
  content->vtemp1    = NULL;
  content->vtemp2    = NULL;
  content->vtemp3    = NULL;
  content->s1        = NULL;
  content->s2        = NULL;
  content->ATimes    = NULL;
  content->ATData    = NULL;
  content->Psetup    = NULL;
  content->Psolve    = NULL;
  content->PData     = NULL;

  /* Allocate content */
  content->r_star = N_VClone(y);
  if (content->r_star == NULL) { return NULL; }
  content->q = N_VClone(y);
  if (content->q == NULL) { return NULL; }
  content->d = N_VClone(y);
  if (content->d == NULL) { return NULL; }
  content->v = N_VClone(y);
  if (content->v == NULL) { return NULL; }
  content->p = N_VClone(y);
  if (content->p == NULL) { return NULL; }
  content->r = N_VCloneVectorArray(2, y);
  if (content->r == NULL) { return NULL; }
  content->u = N_VClone(y);
  if (content->u == NULL) { return NULL; }
  content->vtemp1 = N_VClone(y);
  if (content->vtemp1 == NULL) { return NULL; }
  content->vtemp2 = N_VClone(y);
  if (content->vtemp2 == NULL) { return NULL; }
  content->vtemp3 = N_VClone(y);
  if (content->vtemp3 == NULL) { return NULL; }

  return (S);
}

/* ----------------------------------------------------------------------------
 * Function to set the type of preconditioning for SPTFQMR to use
 */

SUNErrCode SUNLinSolSetPrecType_Mine(SUNLinearSolver S, int pretype)
{
  /* Set pretype */
  MY_CONTENT(S)->pretype = pretype;
  return SUN_SUCCESS;
}

/* ----------------------------------------------------------------------------
 * Function to set the maximum number of iterations for SPTFQMR to use
 */

SUNErrCode SUNLinSolSetMaxl_Mine(SUNLinearSolver S, int maxl)
{
  /* Check for legal pretype */
  if (maxl <= 0) { maxl = MY_MAXL_DEFAULT; }

  /* Set pretype */
  MY_CONTENT(S)->maxl = maxl;
  return SUN_SUCCESS;
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_Mine(SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_ITERATIVE);
}

SUNLinearSolver_ID SUNLinSolGetID_Mine(SUNLinearSolver S)
{
  return (SUNLINEARSOLVER_CUSTOM);
}

SUNErrCode SUNLinSolInitialize_Mine(SUNLinearSolver S)
{
  MySUNLinearSolverContent content;

  /* set shortcut to SPTFQMR memory structure */
  content = MY_CONTENT(S);

  /* ensure valid options */
  if (content->maxl <= 0) { content->maxl = MY_MAXL_DEFAULT; }
  if (content->ATimes == NULL) { return SUN_ERR_ARG_CORRUPT; }

  if ((content->pretype != SUN_PREC_LEFT) &&
      (content->pretype != SUN_PREC_RIGHT) && (content->pretype != SUN_PREC_BOTH))
  {
    content->pretype = SUN_PREC_NONE;
  }

  /* no additional memory to allocate */

  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetATimes_Mine(SUNLinearSolver S, void* ATData,
                                   SUNATimesFn ATimes)
{
  /* set function pointers to integrator-supplied ATimes routine
     and data, and return with success */
  MY_CONTENT(S)->ATimes = ATimes;
  MY_CONTENT(S)->ATData = ATData;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetPreconditioner_Mine(SUNLinearSolver S, void* PData,
                                           SUNPSetupFn Psetup,
                                           SUNPSolveFn Psolve)
{
  /* set function pointers to integrator-supplied Psetup and PSolve
     routines and data, and return with success */
  MY_CONTENT(S)->Psetup = Psetup;
  MY_CONTENT(S)->Psolve = Psolve;
  MY_CONTENT(S)->PData  = PData;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetScalingVectors_Mine(SUNLinearSolver S, N_Vector s1,
                                           N_Vector s2)
{
  /* set N_Vector pointers to integrator-supplied scaling vectors,
     and return with success */
  MY_CONTENT(S)->s1 = s1;
  MY_CONTENT(S)->s2 = s2;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolSetZeroGuess_Mine(SUNLinearSolver S, sunbooleantype onoff)
{
  /* set flag indicating a zero initial guess */
  MY_CONTENT(S)->zeroguess = onoff;
  return SUN_SUCCESS;
}

int SUNLinSolSetup_Mine(SUNLinearSolver S, SUNMatrix A)
{
  int status = SUN_SUCCESS;
  SUNPSetupFn Psetup;
  void* PData;

  /* Set shortcuts to SPTFQMR memory structures */
  Psetup = MY_CONTENT(S)->Psetup;
  PData  = MY_CONTENT(S)->PData;

  /* no solver-specific setup is required, but if user-supplied
     Psetup routine exists, call that here */
  if (Psetup != NULL)
  {
    status = Psetup(PData);
    if (status != 0)
    {
      LASTFLAG(S) = (status < 0) ? SUNLS_PSET_FAIL_UNREC : SUNLS_PSET_FAIL_REC;
      return (LASTFLAG(S));
    }
  }

  /* return with success */
  LASTFLAG(S) = SUN_SUCCESS;
  return SUN_SUCCESS;
}

int SUNLinSolSolve_Mine(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                        N_Vector b, sunrealtype delta)
{
  /* local data and shortcut variables */
  sunrealtype alpha, tau, eta, beta, c, sigma, v_bar, omega;
  sunrealtype rho[2];
  sunrealtype r_init_norm, r_curr_norm;
  sunrealtype temp_val;
  sunbooleantype preOnLeft, preOnRight, scale_x, scale_b, converged, b_ok;
  sunbooleantype* zeroguess;
  int n, m, l_max;
  void *A_data, *P_data;
  SUNATimesFn atimes;
  SUNPSolveFn psolve;
  sunrealtype* res_norm;
  int* nli;
  N_Vector sx, sb, r_star, q, d, v, p, *r, u, vtemp1, vtemp2, vtemp3;
  sunrealtype cv[3];
  N_Vector Xv[3];
  int status = SUN_SUCCESS;

  /* Make local shorcuts to solver variables. */
  l_max     = MY_CONTENT(S)->maxl;
  r_star    = MY_CONTENT(S)->r_star;
  q         = MY_CONTENT(S)->q;
  d         = MY_CONTENT(S)->d;
  v         = MY_CONTENT(S)->v;
  p         = MY_CONTENT(S)->p;
  r         = MY_CONTENT(S)->r;
  u         = MY_CONTENT(S)->u;
  vtemp1    = MY_CONTENT(S)->vtemp1;
  vtemp2    = MY_CONTENT(S)->vtemp2;
  vtemp3    = MY_CONTENT(S)->vtemp3;
  sb        = MY_CONTENT(S)->s1;
  sx        = MY_CONTENT(S)->s2;
  A_data    = MY_CONTENT(S)->ATData;
  P_data    = MY_CONTENT(S)->PData;
  atimes    = MY_CONTENT(S)->ATimes;
  psolve    = MY_CONTENT(S)->Psolve;
  zeroguess = &(MY_CONTENT(S)->zeroguess);
  nli       = &(MY_CONTENT(S)->numiters);
  res_norm  = &(MY_CONTENT(S)->resnorm);

  /* Initialize counters and convergence flag */
  temp_val = r_curr_norm = -ONE;
  *nli                   = 0;
  converged              = SUNFALSE;
  b_ok                   = SUNFALSE;

  /* set sunbooleantype flags for internal solver options */
  preOnLeft  = ((MY_CONTENT(S)->pretype == SUN_PREC_LEFT) ||
               (MY_CONTENT(S)->pretype == SUN_PREC_BOTH));
  preOnRight = ((MY_CONTENT(S)->pretype == SUN_PREC_RIGHT) ||
                (MY_CONTENT(S)->pretype == SUN_PREC_BOTH));
  scale_x    = (sx != NULL);
  scale_b    = (sb != NULL);

  /* Check for unsupported use case */
  if (preOnRight && !(*zeroguess))
  {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUN_ERR_ARG_INCOMPATIBLE;
    return SUN_ERR_ARG_INCOMPATIBLE;
  }

  /* Check if Atimes function has been set */
  if (atimes == NULL) { return SUN_ERR_ARG_CORRUPT; }

  /* If preconditioning, check if psolve has been set */
  if ((preOnLeft || preOnRight) && (psolve == NULL))
  { return SUN_ERR_ARG_CORRUPT; }

  /* Set r_star to initial (unscaled) residual r_star = r_0 = b - A*x_0 */
  /* NOTE: if x == 0 then just set residual to b and continue */
  if (*zeroguess) { N_VScale(ONE, b, r_star); }
  else
  {
    status = atimes(A_data, x, r_star);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC
                                 : SUNLS_ATIMES_FAIL_REC;
      return (LASTFLAG(S));
    }
    N_VLinearSum(ONE, b, -ONE, r_star, r_star);
  }

  /* Apply left preconditioner and b-scaling to r_star (or really just r_0) */
  if (preOnLeft)
  {
    status = psolve(P_data, r_star, vtemp1, delta, SUN_PREC_LEFT);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                 : SUNLS_PSOLVE_FAIL_REC;
      return (LASTFLAG(S));
    }
  }
  else { N_VScale(ONE, r_star, vtemp1); }

  if (scale_b) { N_VProd(sb, vtemp1, r_star); }
  else { N_VScale(ONE, vtemp1, r_star); }

  /* Initialize rho[0] */
  /* NOTE: initialized here to reduce number of computations - avoid need
           to compute r_star^T*r_star twice, and avoid needlessly squaring
           values */
  rho[0] = N_VDotProd(r_star, r_star);

  /* Compute norm of initial residual (r_0) to see if we really need
     to do anything */
  *res_norm = r_init_norm = SUNRsqrt(rho[0]);

  if (r_init_norm <= delta)
  {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUN_SUCCESS;
    return (LASTFLAG(S));
  }

  /* Set v = A*r_0 (preconditioned and scaled) */
  if (scale_x) { N_VDiv(r_star, sx, vtemp1); }
  else { N_VScale(ONE, r_star, vtemp1); }

  if (preOnRight)
  {
    N_VScale(ONE, vtemp1, v);
    status = psolve(P_data, v, vtemp1, delta, SUN_PREC_RIGHT);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                 : SUNLS_PSOLVE_FAIL_REC;
      return (LASTFLAG(S));
    }
  }

  status = atimes(A_data, vtemp1, v);
  if (status != 0)
  {
    *zeroguess = SUNFALSE;
    LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC : SUNLS_ATIMES_FAIL_REC;
    return (LASTFLAG(S));
  }

  if (preOnLeft)
  {
    status = psolve(P_data, v, vtemp1, delta, SUN_PREC_LEFT);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                 : SUNLS_PSOLVE_FAIL_REC;
      return (LASTFLAG(S));
    }
  }
  else { N_VScale(ONE, v, vtemp1); }

  if (scale_b) { N_VProd(sb, vtemp1, v); }
  else { N_VScale(ONE, vtemp1, v); }

  /* Initialize remaining variables */
  N_VScale(ONE, r_star, r[0]);
  N_VScale(ONE, r_star, u);
  N_VScale(ONE, r_star, p);
  N_VConst(ZERO, d);

  /* Set x = sx x if non-zero guess */
  if (scale_x && !(*zeroguess)) { N_VProd(sx, x, x); }

  tau   = r_init_norm;
  v_bar = eta = ZERO;

  /* START outer loop */
  for (n = 0; n < l_max; ++n)
  {
    /* Increment linear iteration counter */
    (*nli)++;

    /* sigma = r_star^T*v */
    sigma = N_VDotProd(r_star, v);

    /* alpha = rho[0]/sigma */
    alpha = rho[0] / sigma;

    /* q = u-alpha*v */
    N_VLinearSum(ONE, u, -alpha, v, q);

    /* r[1] = r[0]-alpha*A*(u+q) */
    N_VLinearSum(ONE, u, ONE, q, r[1]);
    if (scale_x) { N_VDiv(r[1], sx, r[1]); }

    if (preOnRight)
    {
      N_VScale(ONE, r[1], vtemp1);
      status = psolve(P_data, vtemp1, r[1], delta, SUN_PREC_RIGHT);
      if (status != 0)
      {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                   : SUNLS_PSOLVE_FAIL_REC;
        return (LASTFLAG(S));
      }
    }

    status = atimes(A_data, r[1], vtemp1);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC
                                 : SUNLS_ATIMES_FAIL_REC;
      return (LASTFLAG(S));
    }

    if (preOnLeft)
    {
      status = psolve(P_data, vtemp1, r[1], delta, SUN_PREC_LEFT);
      if (status != 0)
      {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                   : SUNLS_PSOLVE_FAIL_REC;
        return (LASTFLAG(S));
      }
    }
    else { N_VScale(ONE, vtemp1, r[1]); }

    if (scale_b) { N_VProd(sb, r[1], vtemp1); }
    else { N_VScale(ONE, r[1], vtemp1); }
    N_VLinearSum(ONE, r[0], -alpha, vtemp1, r[1]);

    /* START inner loop */
    for (m = 0; m < 2; ++m)
    {
      /* d = [*]+(v_bar^2*eta/alpha)*d */
      /* NOTES:
       *   (1) [*] = u if m == 0, and q if m == 1
       *   (2) using temp_val reduces the number of required computations
       *       if the inner loop is executed twice
       */
      if (m == 0)
      {
        temp_val = N_VDotProd(r[1], r[1]);
        temp_val = SUNRsqrt(temp_val);
        omega    = N_VDotProd(r[0], r[0]);
        omega = SUNRsqrt(SUNRsqrt(omega) * temp_val);
        N_VLinearSum(ONE, u, SUNSQR(v_bar) * eta / alpha, d, d);
      }
      else
      {
        omega = temp_val;
        N_VLinearSum(ONE, q, SUNSQR(v_bar) * eta / alpha, d, d);
      }

      /* v_bar = omega/tau */
      v_bar = omega / tau;

      /* c = (1+v_bar^2)^(-1/2) */
      c = ONE / SUNRsqrt(ONE + SUNSQR(v_bar));

      /* tau = tau*v_bar*c */
      tau = tau * v_bar * c;

      /* eta = c^2*alpha */
      eta = SUNSQR(c) * alpha;

      /* x = x+eta*d */
      if (n == 0 && m == 0 && *zeroguess) { N_VScale(eta, d, x); }
      else { N_VLinearSum(ONE, x, eta, d, x); }

      /* Check for convergence... */
      /* NOTE: just use approximation to norm of residual, if possible */
      *res_norm = r_curr_norm = tau * SUNRsqrt(m + 1);

      /* Exit inner loop if iteration has converged based upon approximation
         to norm of current residual */
      if (r_curr_norm <= delta)
      {
        converged = SUNTRUE;
        break;
      }

      /* Decide if actual norm of residual vector should be computed */
      /* NOTES:
       *   (1) if r_curr_norm > delta, then check if actual residual norm
       *       is OK (recall we first compute an approximation)
       *   (2) if r_curr_norm >= r_init_norm and m == 1 and n == l_max, then
       *       compute actual residual norm to see if the iteration can be
       *       saved
       *   (3) the scaled and preconditioned right-hand side of the given
       *       linear system (denoted by b) is only computed once, and the
       *       result is stored in vtemp3 so it can be reused - reduces the
       *       number of psovles if using left preconditioning
       */
      if ((r_curr_norm > delta) ||
          (r_curr_norm >= r_init_norm && m == 1 && n == l_max))
      {
        /* Compute norm of residual ||b-A*x||_2 (preconditioned and scaled) */
        if (scale_x) { N_VDiv(x, sx, vtemp1); }
        else { N_VScale(ONE, x, vtemp1); }

        if (preOnRight)
        {
          status = psolve(P_data, vtemp1, vtemp2, delta, SUN_PREC_RIGHT);
          if (status != 0)
          {
            *zeroguess  = SUNFALSE;
            LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                       : SUNLS_PSOLVE_FAIL_UNREC;
            return (LASTFLAG(S));
          }
          N_VScale(ONE, vtemp2, vtemp1);
        }

        status = atimes(A_data, vtemp1, vtemp2);
        if (status != 0)
        {
          *zeroguess  = SUNFALSE;
          LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC
                                     : SUNLS_ATIMES_FAIL_REC;
          return (LASTFLAG(S));
        }

        if (preOnLeft)
        {
          status = psolve(P_data, vtemp2, vtemp1, delta, SUN_PREC_LEFT);
          if (status != 0)
          {
            *zeroguess  = SUNFALSE;
            LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                       : SUNLS_PSOLVE_FAIL_REC;
            return (LASTFLAG(S));
          }
        }
        else { N_VScale(ONE, vtemp2, vtemp1); }

        if (scale_b) { N_VProd(sb, vtemp1, vtemp2); }
        else { N_VScale(ONE, vtemp1, vtemp2); }

        /* Only precondition and scale b once (result saved for reuse) */
        if (!b_ok)
        {
          b_ok = SUNTRUE;
          if (preOnLeft)
          {
            status = psolve(P_data, b, vtemp3, delta, SUN_PREC_LEFT);
            if (status != 0)
            {
              *zeroguess  = SUNFALSE;
              LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                         : SUNLS_PSOLVE_FAIL_REC;
              return (LASTFLAG(S));
            }
          }
          else { N_VScale(ONE, b, vtemp3); }

          if (scale_b) { N_VProd(sb, vtemp3, vtemp3); }
        }
        N_VLinearSum(ONE, vtemp3, -ONE, vtemp2, vtemp1);
        r_curr_norm = N_VDotProd(vtemp1, vtemp1);
        *res_norm = r_curr_norm = SUNRsqrt(r_curr_norm);

        /* Exit inner loop if inequality condition is satisfied
           (meaning exit if we have converged) */
        if (r_curr_norm <= delta)
        {
          converged = SUNTRUE;
          break;
        }
      }

    } /* END inner loop */

    /* If converged, then exit outer loop as well */
    if (converged == SUNTRUE) { break; }

    /* rho[1] = r_star^T*r_[1] */
    rho[1] = N_VDotProd(r_star, r[1]);

    /* beta = rho[1]/rho[0] */
    beta = rho[1] / rho[0];

    /* u = r[1]+beta*q */
    N_VLinearSum(ONE, r[1], beta, q, u);

    /* p = u+beta*(q+beta*p) = beta*beta*p + beta*q + u */
    cv[0] = SUNSQR(beta);
    Xv[0] = p;

    cv[1] = beta;
    Xv[1] = q;

    cv[2] = ONE;
    Xv[2] = u;

    N_VLinearCombination(3, cv, Xv, p);

    /* v = A*p */
    if (scale_x) { N_VDiv(p, sx, vtemp1); }
    else { N_VScale(ONE, p, vtemp1); }

    if (preOnRight)
    {
      N_VScale(ONE, vtemp1, v);
      status = psolve(P_data, v, vtemp1, delta, SUN_PREC_RIGHT);
      if (status != 0)
      {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                   : SUNLS_PSOLVE_FAIL_REC;
        return (LASTFLAG(S));
      }
    }

    status = atimes(A_data, vtemp1, v);
    if (status != 0)
    {
      *zeroguess  = SUNFALSE;
      LASTFLAG(S) = (status < 0) ? SUNLS_ATIMES_FAIL_UNREC
                                 : SUNLS_ATIMES_FAIL_REC;
      return (LASTFLAG(S));
    }

    if (preOnLeft)
    {
      status = psolve(P_data, v, vtemp1, delta, SUN_PREC_LEFT);
      if (status != 0)
      {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                   : SUNLS_PSOLVE_FAIL_REC;
        return (LASTFLAG(S));
      }
    }
    else { N_VScale(ONE, v, vtemp1); }

    if (scale_b) { N_VProd(sb, vtemp1, v); }
    else { N_VScale(ONE, vtemp1, v); }

    /* Shift variable values */
    /* NOTE: reduces storage requirements */
    N_VScale(ONE, r[1], r[0]);
    rho[0] = rho[1];

  } /* END outer loop */

  /* Determine return value */
  /* If iteration converged or residual was reduced, then return current iterate
   * (x) */
  if ((converged == SUNTRUE) || (r_curr_norm < r_init_norm))
  {
    if (scale_x) { N_VDiv(x, sx, x); }

    if (preOnRight)
    {
      status = psolve(P_data, x, vtemp1, delta, SUN_PREC_RIGHT);
      if (status != 0)
      {
        *zeroguess  = SUNFALSE;
        LASTFLAG(S) = (status < 0) ? SUNLS_PSOLVE_FAIL_UNREC
                                   : SUNLS_PSOLVE_FAIL_UNREC;
        return (LASTFLAG(S));
      }
      N_VScale(ONE, vtemp1, x);
    }

    *zeroguess = SUNFALSE;
    if (converged == SUNTRUE) { LASTFLAG(S) = SUN_SUCCESS; }
    else { LASTFLAG(S) = SUNLS_RES_REDUCED; }
    return (LASTFLAG(S));
  }
  else
  {
    *zeroguess  = SUNFALSE;
    LASTFLAG(S) = SUNLS_CONV_FAIL;
    return (LASTFLAG(S));
  }
}

int SUNLinSolNumIters_Mine(SUNLinearSolver S)
{
  return (MY_CONTENT(S)->numiters);
}

sunrealtype SUNLinSolResNorm_Mine(SUNLinearSolver S)
{
  return (MY_CONTENT(S)->resnorm);
}

N_Vector SUNLinSolResid_Mine(SUNLinearSolver S)
{
  return (MY_CONTENT(S)->vtemp1);
}

sunindextype SUNLinSolLastFlag_Mine(SUNLinearSolver S)
{
  return (LASTFLAG(S));
}

SUNErrCode SUNLinSolSpace_Mine(SUNLinearSolver S, long int* lenrwLS,
                               long int* leniwLS)
{
  sunindextype liw1, lrw1;
  if (MY_CONTENT(S)->vtemp1->ops->nvspace)
  {
    N_VSpace(MY_CONTENT(S)->vtemp1, &lrw1, &liw1);
  }
  else { lrw1 = liw1 = 0; }
  *lenrwLS = lrw1 * 11;
  *leniwLS = liw1 * 11;
  return SUN_SUCCESS;
}

SUNErrCode SUNLinSolFree_Mine(SUNLinearSolver S)
{
  if (S == NULL) { return SUN_SUCCESS; }

  if (S->content)
  {
    /* delete items from within the content structure */
    if (MY_CONTENT(S)->r_star)
    {
      N_VDestroy(MY_CONTENT(S)->r_star);
      MY_CONTENT(S)->r_star = NULL;
    }
    if (MY_CONTENT(S)->q)
    {
      N_VDestroy(MY_CONTENT(S)->q);
      MY_CONTENT(S)->q = NULL;
    }
    if (MY_CONTENT(S)->d)
    {
      N_VDestroy(MY_CONTENT(S)->d);
      MY_CONTENT(S)->d = NULL;
    }
    if (MY_CONTENT(S)->v)
    {
      N_VDestroy(MY_CONTENT(S)->v);
      MY_CONTENT(S)->v = NULL;
    }
    if (MY_CONTENT(S)->p)
    {
      N_VDestroy(MY_CONTENT(S)->p);
      MY_CONTENT(S)->p = NULL;
    }
    if (MY_CONTENT(S)->r)
    {
      N_VDestroyVectorArray(MY_CONTENT(S)->r, 2);
      MY_CONTENT(S)->r = NULL;
    }
    if (MY_CONTENT(S)->u)
    {
      N_VDestroy(MY_CONTENT(S)->u);
      MY_CONTENT(S)->u = NULL;
    }
    if (MY_CONTENT(S)->vtemp1)
    {
      N_VDestroy(MY_CONTENT(S)->vtemp1);
      MY_CONTENT(S)->vtemp1 = NULL;
    }
    if (MY_CONTENT(S)->vtemp2)
    {
      N_VDestroy(MY_CONTENT(S)->vtemp2);
      MY_CONTENT(S)->vtemp2 = NULL;
    }
    if (MY_CONTENT(S)->vtemp3)
    {
      N_VDestroy(MY_CONTENT(S)->vtemp3);
      MY_CONTENT(S)->vtemp3 = NULL;
    }
    free(S->content);
    S->content = NULL;
  }
  if (S->ops)
  {
    free(S->ops);
    S->ops = NULL;
  }
  free(S);
  S = NULL;
  return SUN_SUCCESS;
}
