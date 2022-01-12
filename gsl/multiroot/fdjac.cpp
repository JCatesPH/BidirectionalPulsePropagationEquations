/* multiroots/fdjac.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "../config.h"
#include "gsl_multiroots.h"

// -----
#include <list>
#include <omp.h>
//#include "../../gslStructs.h"
#include "../../gslParams.h"
#define VERBOSE_JAC_COMP

int
gsl_multiroot_fdjacobian (gsl_multiroot_function * F,
                           const gsl_vector * x, const gsl_vector * f,
                           double epsrel, gsl_matrix * jacobian)
{
  const size_t n = x->size;
  const size_t m = f->size;
  const size_t n1 = jacobian->size1;
  const size_t n2 = jacobian->size2;
  int status = 0;

  #ifdef VERBOSE_JAC_COMP
    double initTime = omp_get_wtime();
    double tmpTime = initTime;
    double evalTime = 0.0, colTime = 0.0;
    printf("  --------------------------------\n");
    printf("  Verbose timing of Jacobian enabled..\n");
  #endif

  if (m != n1 || n != n2)
    {
      GSL_ERROR ("function and jacobian are not conformant", GSL_EBADLEN);
    }

  // Load rootparams
  //rootparam_type *rparams = reinterpret_cast<rootparam_type*>(F->params);
  //int tmpOut = rparams->output;
  //rparams->output = 0;
  RootParams *rootObj = reinterpret_cast<RootParams*>(F->params);
  int tmpOut = rootObj->getOutParam();
  rootObj->setOutParam(0);

  #pragma omp parallel //num_threads(1)
  { 
    //rootparam_type *private_rparams;
    //RootParams *tmpRootObjPointer;
    ODEParams private_odeObj = *(rootObj->getODEparams());
    RootParams private_rootObj(&private_odeObj, rootObj->getSizeRoot());
    
    int num_threads = omp_get_num_threads();
    int col_frac = ((int) n) / num_threads;
    int thread_id = omp_get_thread_num();

    size_t i,j;
    gsl_vector *x1, *f1;

    #pragma omp single
    {
      printf("  n = %d, num_threads = %d, col_frac = %d\n", (int) n, num_threads, col_frac);
    }
    
    #pragma omp critical
    {
        /* if (thread_id != 0) {
        //  //private_rparams = copy_rootparams(rparams);
          printf("    Copying RootParams object\n");
          RootParams tmpRootObj  = *rootObj;
          tmpRootObjPointer = &tmpRootObj;
          tmpRootObjPointer = new RootParams;
          *tmpRootObjPointer = *rootObj;
          printf("    Finished copying RootParams object\n");
        }
        else {
          //private_rparams = rparams;
          printf("    Copying RootParams pointer\n");
          tmpRootObjPointer = rootObj;
          printf("    Finished copying RootParams pointer\n");
        }  */
      
      /* printf("    Allocating x1\n"); */
      x1 = gsl_vector_alloc (n);
      //printf("    Finished allocating x1\n");

        if (x1 == 0)
        {
          //#pragma omp error at(execution), severity(fatal), message("ERROR: Failed to allocate space for x1 workspace in fdjac.cpp\n")
          //GSL_ERROR ("failed to allocate space for x1 workspace", GSL_ENOMEM);
          printf("ERROR: Failed to allocate space for x1 workspace in fdjac.cpp\n");
          status = GSL_ENOMEM;
        } 

      f1 = gsl_vector_alloc (m);

        if (f1 == 0)
        {
          gsl_vector_free (x1);

          //#pragma omp error at(execution), severity(fatal), message("ERROR: Failed to allocate space for x1 workspace in fdjac.cpp\n")
          //GSL_ERROR ("failed to allocate space for f1 workspace", GSL_ENOMEM);
          printf("ERROR: Failed to allocate space for f1 workspace in fdjac.cpp\n");
          status = GSL_ENOMEM;
          
        }
    }
    if (status)
    {
      #pragma omp cancel parallel
    }
    #pragma omp barrier

    /* printf("    Copying vector into tmp\n"); */
    //printf("Thread %d has numT in ODEParams obj set to %d\n", thread_id, tmpRootObjPointer->getODEparams()->numT);
    //printf("Thread %d has numT in ODEParams obj set to %d\n", thread_id, private_rootObj.getODEparams()->numT);
    
    gsl_vector_memcpy (x1, x);  /* copy x into x1 */

    #ifdef VERBOSE_JAC_COMP
      //printf("    Initial vector allocation and copying cost %.2f [s]\n    Entering loop to estimate derivatives.\n", omp_get_wtime() - tmpTime);
      tmpTime = omp_get_wtime();
    #endif

    // Determine last column index based on whether it is last thread
    int last_index;
    if (thread_id == num_threads - 1)
      last_index = n;
    else
      last_index = (thread_id + 1) * col_frac;
    
    // This is when actual Jacobian computations begin
    for (j = thread_id * col_frac; j < last_index; j++)
      {
        double xj = gsl_vector_get (x, j);
        //printf("    Element %d retrieved in thread %d\n", (int) j, thread_id);
        double dx = epsrel * fabs (xj);

        if (dx == 0)
          {
            dx = epsrel;
          }

        gsl_vector_set (x1, j, xj + dx);
        
        {
          #ifdef VERBOSE_JAC_COMP
            tmpTime = omp_get_wtime();
          #endif

          //int f_stat = GSL_MULTIROOT_FN_EVAL (F, x1, f1);
          int f_stat = GSL_MULTIROOT_FN_EVAL_WPARAMS (F, x1, &private_rootObj, f1);

          #ifdef VERBOSE_JAC_COMP
            evalTime += omp_get_wtime() - tmpTime;
          #endif


          if (f_stat != GSL_SUCCESS) 
            {
              status = GSL_EBADFUNC;  
              break; /* n.b. avoid memory leak for x1,f1 */
            }
        }

        gsl_vector_set (x1, j, xj);

        #ifdef VERBOSE_JAC_COMP
          tmpTime = omp_get_wtime();
        #endif
        for (i = 0; i < m; i++)
          {
            double g1 = gsl_vector_get (f1, i);
            double g0 = gsl_vector_get (f, i);
            gsl_matrix_set (jacobian, i, j, (g1 - g0) / dx);
          }

        #ifdef VERBOSE_JAC_COMP
          colTime += omp_get_wtime() - tmpTime;
        #endif

        {
          gsl_vector_view col = gsl_matrix_column (jacobian, j);
          int null_col = gsl_vector_isnull (&col.vector);
          /* if column is null, return an error - this may be due to
             dx being too small. Try increasing epsrel */
          if (null_col) {
            status = GSL_ESING;
          }
        }
      }
    #ifdef VERBOSE_JAC_COMP
      printf("      Thread number %d information: \n", thread_id);
      printf("       Time spent evaluating mapping is %.2f [s]\n       Time spent computing columns is %.2f [s]\n", evalTime, colTime);
      if (status == GSL_EBADFUNC) {
        printf("       WARNING: GSL_EBADFUNC sent\n");
      }
      else if (status == GSL_ESING) {
        printf("       WARNING: GSL_ESING sent (Jacobian singularity)\n");
      }
      else {
        printf("       status : %d\n", status);
      }
    #endif
    gsl_vector_free (x1);
    gsl_vector_free (f1);

    /* #ifdef VERBOSE_JAC_COMP
      printf("       Finished freeing temporary gsl vectors.\n");
    #endif */
    /* if (thread_id != 0) {
      free_rootparams(tmpRootObj);
    } */

  }
  #ifdef VERBOSE_JAC_COMP
    /* printf("    Exiting Jacobian computation..\n"); */
    printf("  --------------------------------\n");
  #endif

  //rparams->output = tmpOut;
  rootObj->setOutParam(tmpOut);

  if (status)
    return status;
  else
    return GSL_SUCCESS;
}

