/* nag_1d_quad_gen_1 (d01sjc) Example Program.
 *http://www.nag.co.uk/numeric/cl/nagdoc_cl25/html/d01/d01sjc.html (url for specification webpage)
 * Copyright 2014 Numerical Algorithms Group.
 *
 * Mark 5, 1998.
 * Mark 6 revised, 2000.
 * Mark 7 revised, 2001.
 *
 */

#include <nag.h>
#include <stdio.h>
#include <nag_stdlib.h>
#include <math.h>
#include <nagd01.h>
#include <nagx01.h>

#ifdef __cplusplus
extern "C" {
#endif
static double NAG_CALL f(double x, Nag_User *comm);
#ifdef __cplusplus
}
#endif

int main(void)
{
  static Integer use_comm[1] = {1};
  Integer          exit_status = 0;
  double           a, b;
  double           epsabs, abserr, epsrel, result;
  Nag_QuadProgress qp;
  Integer          max_num_subint;
  NagError         fail;
  /* nag_pi (x01aac).
   * pi
   */
  double           pi = nag_pi;
  Nag_User         comm;

  INIT_FAIL(fail);

  printf("nag_1d_quad_gen_1 (d01sjc) Example Program Results\n");

  /* For communication with user-supplied functions: */
  comm.p = (Pointer)&use_comm;

  epsabs = 0.0;
  epsrel = 0.0001;
  a = 0.0;
  b = pi*2.0;
  max_num_subint = 200;
  /* nag_1d_quad_gen_1 (d01sjc).
   * One-dimensional adaptive quadrature, allowing for badly
   * behaved integrands, thread-safe
   */
  nag_1d_quad_gen_1(f, a, b, epsabs, epsrel, max_num_subint, &result, &abserr,
                    &qp, &comm, &fail);
  printf("a      - lower limit of integration = %10.4f\n", a);
  printf("b      - upper limit of integration = %10.4f\n", b);
  printf("epsabs - absolute accuracy requested = %11.2e\n", epsabs);
  printf("epsrel - relative accuracy requested = %11.2e\n\n", epsrel);
  if (fail.code != NE_NOERROR)
    printf("Error from nag_1d_quad_gen_1 (d01sjc) %s\n", fail.message);
  if (fail.code != NE_INT_ARG_LT && fail.code != NE_ALLOC_FAIL &&
      fail.code != NE_NO_LICENCE)
    {
      printf("result - approximation to the integral = %9.5f\n",
              result);
      printf("abserr - estimate of the absolute error = %11.2e\n",
              abserr);
      printf("qp.fun_count  - number of function evaluations = %4ld\n",
              qp.fun_count);
      printf("qp.num_subint  - number of subintervals used = %4ld\n",
              qp.num_subint);
      /* Free memory used by qp */
      NAG_FREE(qp.sub_int_beg_pts);
      NAG_FREE(qp.sub_int_end_pts);
      NAG_FREE(qp.sub_int_result);
      NAG_FREE(qp.sub_int_error);
    }
  else
    {
      exit_status = 1;
      goto END;
    }

 END:
  return exit_status;
}

static double NAG_CALL f(double x, Nag_User *comm)
{
  /* nag_pi (x01aac), see above. */
  double pi = nag_pi;
  Integer *use_comm = (Integer *)comm->p;

  if (use_comm[0])
    {
      printf("(User-supplied callback f, first invocation.)\n");
      use_comm[0] = 0;
    }

  return(x*sin(x*30.0)/sqrt(1.0-x*x/(pi*pi*4.0)));
}
