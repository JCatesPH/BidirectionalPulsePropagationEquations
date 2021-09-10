#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

struct rparams
  {
    double v0;
    double xf;
  };

typedef struct{
  double g;
} odeparam_type;

int func(double t, const double y[], double dydt[], void *odep) {
  odeparam_type *p = reinterpret_cast<odeparam_type*>(odep);
  dydt[0] = y[2];
  dydt[1] = y[3];
  dydt[2] = 0.0;
  dydt[3] = -(p->g);
  return GSL_SUCCESS;
}

int
map (const gsl_vector * x, void *params,
              gsl_vector * f)
{
  odeparam_type* odep = (odeparam_type*)malloc(sizeof(odeparam_type));
  odep->g = 10.0;
  const gsl_odeiv2_step_type * stepType = gsl_odeiv2_step_rk4;
	gsl_odeiv2_step * gslStep = gsl_odeiv2_step_alloc(stepType, 4);
	gsl_odeiv2_control * gslControl = gsl_odeiv2_control_y_new(1e-5, 1e-4);
	gsl_odeiv2_evolve * gslEvolve = gsl_odeiv2_evolve_alloc(4);
	gsl_odeiv2_system sys = { func, NULL, (size_t)(4), odep};

  double v0 = ((struct rparams *) params)->v0;
  double xf = ((struct rparams *) params)->xf;
  const double theta = gsl_vector_get(x, 0);
  double tolh = 1e-12;
  double dt = 0.1, t = 0.0;
  int GSLerrorFlag;

  double y[4] = {0.0, 0.0, v0*cos(theta), v0*sin(theta)};
  double ytmp[4], ttmp;
  GSLerrorFlag = gsl_odeiv2_evolve_apply_fixed_step(gslEvolve, gslControl, gslStep, &sys, &t, dt, y);
  while (y[1] > tolh) {
    GSLerrorFlag = gsl_odeiv2_evolve_apply_fixed_step(gslEvolve, gslControl, gslStep, &sys, &t, dt, y);
    if (y[1] < 0.0) {
      dt = dt / 2.0;
      t = ttmp;
      y[0] = ytmp[0];
      y[1] = ytmp[1];
      y[2] = ytmp[2];
      y[3] = ytmp[3];

    }
    ytmp[0] = y[0];
    ytmp[1] = y[1];
    ytmp[2] = y[2];
    ytmp[3] = y[3];
    ttmp = t;
    //printf("%.3f, %.3f\n", y[0], y[1]);
  }
  
  const double f0 = y[0] - xf;
  gsl_vector_set(f, 0, f0);

  free(odep);
  gsl_odeiv2_evolve_free(gslEvolve);
  gsl_odeiv2_control_free(gslControl);
  gsl_odeiv2_step_free(gslStep);
  return GSL_SUCCESS;
}

int
print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf("iter = %3u x = % .6f  "
          "f(x) = % .6e \n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->f, 0));
}

int
main (void)
{
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t iter = 0;

  const size_t n = 1;
  double v0 = 5.0;
  double xf = 0.25;
  struct rparams p = {v0, xf};
  gsl_multiroot_function f = {&map, n, &p};

  gsl_vector *u = gsl_vector_alloc (n);
  const double theta0 = 1.0;

  gsl_vector_set (u, 0, theta0);

  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, n);
  gsl_multiroot_fsolver_set (s, &f, u);

  print_state (iter, s);

  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

      print_state (iter, s);

      if (status)   /* check if solver is stuck */
        break;

      status =
        gsl_multiroot_test_residual (s->f, 1e-7);
    }
  while (status == GSL_CONTINUE && iter < 1000);

  printf ("status = %s\n", gsl_strerror (status));

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (u);
  return 0;
}
