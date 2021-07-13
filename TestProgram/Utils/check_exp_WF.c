// clang-format off
/*******************************************************************************
*
* Check for the implementation of the exponential using the Horner scheme 
* by Fernando Romero-Lopez
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"
#include "random.h"
#include "wilsonflow.h"
#include "setup.h"


void print_suNg(suNg *m){
    int i,j;
    double sumsq = 0;
    double traceRe = 0;

    lprintf("SUNg",1,"\n");
#ifdef GAUGE_SPN
#define _NG NG/2
#else
#define _NG NG
#endif
    for(i=0;i<_NG;++i){
        for(j=0;j<NG;++j){
            lprintf("SUNg",1,"%+-.3f,%+-.3f ", creal(m->c[i*NG+j]), cimag(m->c[i*NG+j]));
            sumsq += cimag(m->c[i*NG+j])*cimag(m->c[i*NG+j]);
            sumsq += creal(m->c[i*NG+j])*creal(m->c[i*NG+j]);
            if(i==j) traceRe +=creal(m->c[i*NG+j]);
        }
    }
#undef _NG
    lprintf("SUNg",1,"SUM SQUARES: %f\n", sumsq);
    lprintf("SUNg",1,"TRACE RE: %f\n", traceRe);
}


int main(int argc, char *argv[])
{
  int i = 3, j = 0, evaluations = 1;
  suNg test, exptest, exptest2;
  double complex tr;

  //setup_process(&argc, &argv);
  //setup_gauge_fields();

  random_suNg(&test);

  // Traceless and real in the diagonal
  // (then multiplied by I)
#ifdef GAUGE_SPN
#define _MATSIZE NG*NG/2
#else
#define _MATSIZE NG*NG
#endif
  for (i = 0; i*(NG+1) < _MATSIZE; i++)
    test.c[i * (NG + 1)] = 0.5 * (test.c[i * (NG + 1)] + conj(test.c[i * (NG + 1)])) * I;
#undef _MATSIZE



// SPN matrices with purely imaginary diagonal 
// are already traceless
// when the right group constraints are enforced.
// (which is the case when they are stored in compressed form)
// See https://arxiv.org/pdf/1712.04220.pdf, Eq A.3
#ifndef GAUGE_SPN
  // Make matrix traceless
  _suNg_trace(tr, test);
  test.c[NG * NG - 1] = test.c[NG * NG - 1] - tr;
#endif

#ifndef GAUGE_SPN
  //Change offdiagonal elements to anti-hermitean matrix
  for (i = 0; i < NG; i++)
    for (j = 0; j < i; j++)
    {
      test.c[NG * (i) + j] = - conj(test.c[NG * (j) + i]);
    }
#else
  // The conditions for a SP(N) generator are more complicated.
  //  See https://arxiv.org/pdf/1712.04220.pdf, Eq A.4 and subsequent discussion.
  for (i = 0; i < NG/2; i++)
  {
    for (j = 0; j < i; j++)
    {
      // For the A part , change offdiagonal elements to anti-hermitean matrix
      test.c[NG * (i) + (j)       ] = - conj(test.c[NG * (j) + (i)       ]);
      // For the B part , iB is symmetric
      test.c[NG * (i) + (j + NG/2)] =        test.c[NG * (j) + (i + NG/2)];
    }
  }
#endif

  _suNg_trace(tr, test);

  if (creal(conj(tr)*tr) > 1.e-30)
  {
    lprintf("ERROR", 0, "random matrix not traceless!!\nTrace = %.*e + I*%.*e\n",18, creal(tr),18, cimag(tr));
    return 1;
  }

  for (i = 0; i < evaluations; i++)
    suNg_Exp_Taylor(&exptest, &test);

  for (i = 0; i < evaluations; i++)
    suNg_Exp(&exptest2, &test);

  _suNg_sub_assign(exptest, exptest2);

  double norm = 0.;

#ifdef GAUGE_SPN
#define _MATSIZE NG*NG/2
#else
#define _MATSIZE NG*NG
#endif
  for (i = 0; i < _MATSIZE; i++)
    norm += conj(exptest.c[i]) * exptest.c[i];
#undef _MATSIZE

#ifdef GAUGE_SPN
  norm *=2;
#endif

  norm = sqrt(norm);

  int check = 0;

  lprintf("MAIN", 0, "Maximal normalized difference = %.2e\n", norm / (NG));
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
  if (norm / (NG) > 1.e-14){
    lprintf("MAIN",0,"Exptest1:\n");
    print_suNg(&exptest);

    lprintf("MAIN",0,"Exptest2:\n");
    print_suNg(&exptest2);


    check++;
  }

  //finalize_process();
  return check;
}

// clang-format on
