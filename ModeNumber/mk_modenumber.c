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

typedef struct _input_nu {
  char configlist[256];
  double inverr2;
  char approx[512];
  int nhits;
  double mass;
  char list[1024];

  /* for the reading function */
  input_record_t read[8];

} input_nu;

#define init_input_nu(varname) \
{ \
  .read={\
    {"configuration list", "meas:configlist = %s", STRING_T, &(varname).configlist},\
    {"squared error for inverter", "nu:inverr2 = %lf", DOUBLE_T, &(varname).inverr2},\
    {"Chebyshev approximation file", "nu:approx = %s", STRING_T, (varname).approx},\
    {"number of stochastic spinors", "nu:nhits = %d", INT_T, &(varname).nhits},\
    {"quark quenched mass", "meas:mass = %lf", DOUBLE_T, &(varname).mass},\
    {"list of eigenvalues", "nu:list = %s", STRING_T, (varname).list},\
    {NULL, NULL, INT_T, NULL}\
  }\
}


char cnfg_filename[256]="";
char list_filename[256]="";
char input_filename[256] = "input_file";
enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };

input_nu nu_var = init_input_nu(nu_var);


typedef struct {
  char string[256];
  int t, x, y, z;
  int nc, nf;
  double b, m;
  int n;
  int type;
} filename_t;


double hevamass=0.;
void HEVA(spinor_field *out, spinor_field *in){
  g5Dphi_sq(hevamass, out, in);
}

int main(int argc,char *argv[]) {
  char tmp[1024];
  FILE* list;
  int neig;
  char* cptr;
  double M[1024];

  /* setup process communications */
  setup_process(&argc,&argv);
  setup_gauge_fields();

  read_input(glb_var.read, get_input_filename());
  read_input(nu_var.read, get_input_filename());
  read_input(rlx_var.read, get_input_filename());

  strcpy(list_filename, nu_var.configlist);

  lprintf("MAIN", 0, "list_filename = %s %s\n", list_filename, nu_var.configlist);
  if (strcmp(list_filename, "") != 0) {
    error((list = fopen(list_filename, "r")) == NULL, 1, "main [measure_spectrum.c]", "Failed to open list file\n");
  }

  strcpy(tmp,nu_var.list);
  cptr = strtok(tmp, ";");
  neig=0;
  while(cptr != NULL) {
    M[neig]=atof(cptr);
    neig++;
    cptr = strtok(NULL, ";");
  }
  error(neig==0,1,"mk_modenumber.c","neig == 0 !!!");
  

  init_modenumber(nu_var.mass, nu_var.inverr2, nu_var.nhits, nu_var.approx);
  for(int k = 0; k < neig; k++)
    lprintf("MODENUMBER",0,"M[%d] = %e\n",k,M[k]);

  int i=0;
  while(1) {

    if(list!=NULL)
      if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;

    i++;

    lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
    /* NESSUN CHECK SULLA CONSISTENZA CON I PARAMETRI DEFINITI !!! */
    read_gauge_field(cnfg_filename);
    represent_gauge_field();

    lprintf("TEST",0,"<p> %1.6f\n",avr_plaquette());

    /* full_plaquette(); */

    for(int k = 0; k < neig; k++) {
      double number = ModeNumber(M[k]*M[k]);
      int mvm = getMVM();
      lprintf("MODENUMBER",0,"nu[ %e ] = %.2f\n",M[k],number);
      lprintf("MODENUMBER",0,"MVM = %d\n",mvm);
    }

    if(list==NULL) break;
  }

  if(list!=NULL) fclose(list);

  free_BCs();
  
  free_modenumber();

  free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  free_gfield_f(u_gauge_f);
#endif

  finalize_process();

  return 0;
}

