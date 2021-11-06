/***************************************************************************\
 * Copyright (c) 2009, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

/*******************************************************************************
*
* Computation of the lowest eigenvalues of H^2
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
#include "setup.h"

#include "cinfo.c"

#if defined(ROTATED_SF) && defined(BASIC_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif

/* Mesons parameters */
typedef struct _input_eigval {
  /* EVA parameters */
  char configlist[256];
  int nevt; /* search space dimension */
  int nev; /* number of accurate eigenvalues */
  int kmax; /* max degree of polynomial */
  int maxiter; /* max number of subiterations */
  double omega1; /* absolute precision */
  double omega2; /* relative precision */
  double mass; /* quenched mass */

  /* for the reading function */
  input_record_t read[10];

} input_eigval;

#define init_input_eigval(varname) \
{ \
  .read={\
    {"configuration list", "meas:configlist = %s", STRING_T, &(varname).configlist}, \
    {"search space dimension", "eva:nevt = %d", INT_T, &(varname).nevt},\
    {"number of accurate eigenvalues", "eva:nev = %d", INT_T, &(varname).nev},\
    {"max degree of polynomial", "eva:kmax = %d", INT_T, &(varname).kmax},\
    {"max number of subiterations", "eva:maxiter = %d", INT_T, &(varname).maxiter},\
    {"absolute precision", "eva:omega1 = %lf", DOUBLE_T, &(varname).omega1},\
    {"relative precision", "eva:omega2 = %lf", DOUBLE_T, &(varname).omega2},\
    {"quark quenched mass", "meas:mass = %lf", DOUBLE_T, &(varname).mass}, \
    {NULL, NULL, INT_T, NULL}\
  }\
}


char cnfg_filename[256]="";
char list_filename[256]="";
char input_filename[256] = "input_file";
enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };

input_eigval eig_var = init_input_eigval(eig_var);


typedef struct {
  char string[256];
  int t, x, y, z;
  int nc, nf;
  double b, m;
  int n;
  int type;
} filename_t;


double hevamass=0.;
void H2EVA(spinor_field *out, spinor_field *in){
  g5Dphi_sq(hevamass, out, in);
}
void HEVA(spinor_field *out, spinor_field *in){
  g5Dphi(hevamass, out, in);
}

int main(int argc,char *argv[]) {
  int i,n;
  FILE* list;

  /* setup process communications */
  setup_process(&argc,&argv);
  setup_gauge_fields();

  read_input(glb_var.read, get_input_filename());
  read_input(eig_var.read, get_input_filename());
  read_input(rlx_var.read, get_input_filename());

  strcpy(list_filename, eig_var.configlist);

  lprintf("MAIN", 0, "list_filename = %s %s\n", list_filename, eig_var.configlist);
  if (strcmp(list_filename, "") != 0) {
    error((list = fopen(list_filename, "r")) == NULL, 1, "main [mk_eigvals.c]", "Failed to open list file\n");
  }

  lprintf("MAIN",0,"EVA Parameters:\n");
  lprintf("MAIN",0,"search space dimension  (eva:nevt) = %d\n",eig_var.nevt);
  lprintf("MAIN",0,"number of accurate eigenvalues (eva:nev) = %d\n",eig_var.nev);
  lprintf("MAIN",0,"max degree of polynomial (eva:kmax) = %d\n",eig_var.kmax);
  lprintf("MAIN",0,"max number of subiterations (eva:maxiter) = %d\n",eig_var.maxiter);
  lprintf("MAIN",0,"absolute precision  (eva:omega1) = %e\n",eig_var.omega1);
  lprintf("MAIN",0,"relative precision (eva:omega2) = %e\n",eig_var.omega2);
  hevamass = eig_var.mass;
  lprintf("MAIN",0,"mass = %f\n",hevamass);


  /* EVA parameters */
  double max, mupp;
  double *eva_val;
  int status,ie;
  spinor_field *eva_vec, *eva_ws;

  eva_val=malloc(sizeof(double)*eig_var.nevt);
  eva_vec=alloc_spinor_field_f(eig_var.nevt+1,&glattice);
  eva_ws=eva_vec+eig_var.nevt;

  mupp=fabs(hevamass+4)+4;
  mupp*=mupp;
  /* END of EVA parameters */

  i=0;
  while(1) {

    if(list!=NULL)
      if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;

    i++;

    lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
    /* NESSUN CHECK SULLA CONSISTENZA CON I PARAMETRI DEFINITI !!! */
    read_gauge_field(cnfg_filename);
    represent_gauge_field();

    lprintf("TEST",0,"<p> %1.6f\n",avr_plaquette());


    int MVM=0; /* counter for matrix-vector multiplications */

    max_H(&H2EVA, &glattice, &max);
    lprintf("MAIN",0,"MAXCHECK: cnfg=%e  uppbound=%e diff=%e %s\n",max,mupp,mupp-max,(mupp-max)<0?"[FAILED]":"[OK]");
    max*=1.1;

    ie=eva(eig_var.nev,eig_var.nevt,0,eig_var.kmax,eig_var.maxiter,max,eig_var.omega1,eig_var.omega2,&H2EVA,eva_vec,eva_val,&status);
    MVM+=status;
    while (ie!=0) { /* if failed restart EVA */
      lprintf("MAIN",0,"Restarting EVA!\n");
      ie=eva(eig_var.nev,eig_var.nevt,2,eig_var.kmax,eig_var.maxiter,max,eig_var.omega1,eig_var.omega2,&H2EVA,eva_vec,eva_val,&status);
      MVM+=status;
    }

    lprintf("MAIN",0,"EVA MVM = %d\n",MVM);
    for (n=0;n<eig_var.nev;++n) {
      HEVA(&eva_ws[0],&eva_vec[n]);
      lprintf("RESULT",0,"Eig %d = %.15e %.15e\n",n,eva_val[n],
        spinor_field_prod_re_f(&eva_ws[0],&eva_vec[n])/spinor_field_sqnorm_f(&eva_vec[n]));
    }
    
    if(list==NULL) break;
  }

  if(list!=NULL) fclose(list);

  free_BCs();

  free(eva_val);
  free_spinor_field_f(eva_vec);

  free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  free_gfield_f(u_gauge_f);
#endif

  finalize_process();

  return 0;
}

