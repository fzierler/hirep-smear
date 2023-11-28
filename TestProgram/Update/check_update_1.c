/******************************************************************************* 
 * Check that the molecular dynamics evolution is reversible
 * 
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
#include "communications.h"
#include "../../HMC/hmc_utils.h"

int nhb,nor,nit,nth,nms,level,seed;
double beta;

extern suNg_av_field *suN_momenta;
suNg_field *u_gauge_backup=NULL;

hmc_flow flow=init_hmc_flow(flow);

static void flip_mom()
{
    geometry_descriptor *gd=suN_momenta->type;

    _MASTER_FOR(gd,ix) {
        for(int dx=0;dx <4 ; dx++) {
            suNg_algebra_vector  *dptr=(suNg_algebra_vector*)(&suN_momenta->ptr[4*ix+dx]);
            _algebra_vector_mul_g(*dptr,-1.0,*dptr);
        }
    }
}

static void copy_mom()
{
    if(suN_momenta_backup==NULL) suN_momenta_backup = alloc_avfield(suN_momenta->type);
    geometry_descriptor *gd=suN_momenta->type;
    double sqnorm = 0.0;//DEBUG

    _MASTER_FOR(gd,ix) {
        for(int dx=0;dx <4 ; dx++) {
            suNg_algebra_vector *dptr=(suNg_algebra_vector*)(&suN_momenta->ptr[4*ix+dx]);
            suNg_algebra_vector *dptr_backup=(suNg_algebra_vector*)(&suN_momenta_backup->ptr[4*ix+dx]);
            _algebra_vector_sub_assign_g(*dptr_backup,*dptr_backup);//=0
            _algebra_vector_add_assign_g(*dptr_backup,*dptr);
            double sqnorm_tmp;
            _algebra_vector_sqnorm_g(sqnorm_tmp,*dptr_backup);
            sqnorm +=sqnorm_tmp;
 
        }
    }
    lprintf("GHMC",0,"copy mom, sqnorm: %1.8e\n",sqnorm);
}

static void copy_conf()
{
    if(u_gauge_backup==NULL) u_gauge_backup = alloc_gfield(u_gauge->type);
    geometry_descriptor *gd=u_gauge->type;

    _MASTER_FOR(gd,ix) {
        for(int dx=0;dx <4 ; dx++) {
            suNg *dptr=(suNg*)(&u_gauge->ptr[4*ix+dx]);
            suNg *dptr_backup=(suNg*)(&u_gauge_backup->ptr[4*ix+dx]);
            _suNg_zero(*dptr_backup);
            _suNg_add_assign(*dptr_backup,*dptr);
        }
    }
}

static double confdiffsq()
{
    geometry_descriptor *gd=u_gauge->type;
    double sqnorm = 0.0;

    _MASTER_FOR(gd,ix) {
        for(int dx=0;dx <4 ; dx++) {
            suNg *dptr=(suNg*)(&u_gauge->ptr[4*ix+dx]);
            suNg *dptr_backup=(suNg*)(&u_gauge_backup->ptr[4*ix+dx]);
            _suNg_sub_assign(*dptr_backup,*dptr);
            double sqnorm_tmp;
            _suNg_sqnorm(sqnorm_tmp,*dptr_backup);
            sqnorm +=sqnorm_tmp;
        }
    }
    return sqnorm;
}

static double momsumsq()
{
    geometry_descriptor *gd=suN_momenta->type;
    double sqnorm = 0.0;

    _MASTER_FOR(gd,ix) {
        for(int dx=0;dx <4 ; dx++) {
            suNg_algebra_vector  *dptr=(suNg_algebra_vector*)(&suN_momenta->ptr[4*ix+dx]);
            suNg_algebra_vector  *dptr_backup=(suNg_algebra_vector*)(&suN_momenta_backup->ptr[4*ix+dx]);
            _algebra_vector_add_assign_g(*dptr_backup,*dptr);
            double sqnorm_tmp;
            _algebra_vector_sqnorm_g(sqnorm_tmp,*dptr_backup);
            sqnorm +=sqnorm_tmp;
        }
    }
    return sqnorm;
}


int main(int argc,char *argv[])
{

    char tmp[256];

    setup_process(&argc,&argv);

    logger_setlevel(0,100); /* log all */
    if (PID!=0) { 
        logger_disable();}
    else{
        sprintf(tmp,">out_%d",PID); logger_stdout(tmp);
        sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);
    }

    logger_map("DEBUG","debug");

    lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 

    read_input(glb_var.read,"test_input");


    /* setup communication geometry */
    if (geometry_init() == 1) {
        finalize_process();
        return 0;
    }

    geometry_mpi_eo();

    /* setup random numbers */
    read_input(rlx_var.read,"test_input");
    lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
    rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */

#ifdef GAUGE_SPN
    lprintf("MAIN",0,"Gauge group: SP(%d)\n",NG);
#elif defined(GAUGE_SUN)
    lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
#endif
    lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);


    /* test_geometry_mpi_eo(); */

    /* Init Monte Carlo */
    init_mc(&flow, "test_input");
    copy_conf();
    lprintf("MAIN",0,"MVM during RHMC initialzation: %ld\n",getMVM());
    lprintf("MAIN",0,"Initial plaquette: %1.8e\n",avr_plaquette());

    gaussian_momenta(suN_momenta); 
    copy_mom();
    
    /* generate new pseudofermions */
    for (int i=0;i<num_mon();++i) {
      const monomial *m = mon_n(i);
      m->gaussian_pf(m);
    }

 
    int rr=update_ghmc_stripped();

    if(rr<0) {
        lprintf("REV TEST",0,"Error in updating the gauge field!!\n");
        return 1;
    }
    lprintf("REV TEST",0,"Plaquette: %1.8e\n",avr_plaquette());

    flip_mom();
    //((suNg_algebra_vector*)(&suN_momenta->ptr[0]))->c[0] += 0.1;// TO BREAK IT


    rr=update_ghmc_stripped();
    if(rr<0) {
        lprintf("REV TEST",0,"Error in updating the gauge field!!\n");
        return 1;
    }
    lprintf("REV TEST",0,"Plaquette: %1.8e\n",avr_plaquette());

    lprintf("REV TEST",0,"Momenta sum: %1.8e\n",momsumsq());
    lprintf("REV TEST",0,"Conf Diff: %1.8e\n",confdiffsq());

    free_ghmc();

    free_gfield(u_gauge);
#ifndef REPR_FUNDAMENTAL
    free_gfield_f(u_gauge_f);
#endif
    if(u_gauge_backup!=NULL) {
        free_avfield(u_gauge_backup); 
        u_gauge_backup=NULL; 
    } 


    finalize_process();
    return 0;
}
