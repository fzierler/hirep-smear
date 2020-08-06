/*******************************************************************************
*
*  NOCOMPILE= WITH_MPI
*
* Test of modules
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
#include "global.h"
#include "suN.h"
#include "suN_types.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "setup.h"
#include "logger.h"

#ifdef REPR_FUNDAMENTAL
#define COMPLEX_REP
#ifdef GAUGE_SPN
static float C2=(float)NG*(NG+1)/(4*NG);
#else
static float C2=(float)(NG*NG-1)/(2*NG);
#endif
static float Tr = 0.5;
#endif


#ifdef REPR_ADJOINT
#ifdef GAUGE_SPN
static float C2=(float)(NG+2)/2; 
static float Tr=(float)(NG+2)/2; 
#else
static float C2=(float)NG; 
static float Tr=(float)NG;
#endif
#endif

#ifdef REPR_ANTISYMMETRIC
#define COMPLEX_REP
#ifdef GAUGE_SPN
static float C2=(float)(NG-2)*NG*(NG+1)/(2*NG*(NG-1)-4);
#else
static float C2=(float)(NG-2)*(NG+1)/NG;
#endif
static float Tr=(float)(NG-2)/2;
#endif

#ifdef REPR_SYMMETRIC
#define COMPLEX_REP
#ifdef GAUGE_SPN
#error "Use Adjoint representation."
#endif
static float C2 = (float)(NG + 2) * (NG - 1) / (float)NG;
static float Tr = (float)(NG + 2) / 2;
#endif

#ifdef GAUGE_SUN
static int dAdj = NG * NG - 1;
static float fund = (float)(NG * NG - 1) / (2 * (float)(NG));
#elif defined(GAUGE_SPN)
static int dAdj=NG*(NG+1)/2;
static float fund=(float)(NG*(NG+1)/2)/(2*(float)(NG));
#endif

static int error_compare(double x, double y){
    const double threshold = 1.0e-13;
    return (x < y - threshold) || (x > y + threshold);
}

void print_suNf(suNf *m){
    int i,j;
    double sumsq = 0;
    double traceRe = 0;

    lprintf("SUNf",1,"\n");
#if defined(GAUGE_SPN) && defined(REPR_FUNDAMENTAL)
#define _NROWS NF/2
#else
#define _NROWS NF/2
#endif
    for(i=0;i<_NROWS;++i){
#undef _NROWS
        for(j=0;j<NF;++j){
#ifdef COMPLEX_REP
            lprintf("SUNf",1,"%+-.3f,%+-.3f ", creal(m->c[i*NF+j]), cimag(m->c[i*NF+j])); 
            sumsq += cimag(m->c[i*NF+j])*cimag(m->c[i*NF+j]);
            sumsq += creal(m->c[i*NF+j])*creal(m->c[i*NF+j]);
            if(i==j) traceRe +=creal(m->c[i*NF+j]);
#else
            lprintf("SUNf",1,"%+-.3f,%+-.3f ", m->c[i*NF+j], m->c[i*NF+j]); 
            sumsq += m->c[i*NF+j]*m->c[i*NF+j];
            if(i==j) traceRe +=m->c[i*NF+j];
#endif
        }
    }
    lprintf("SUNf",1,"SUM SQUARES: %f\n", sumsq);
    lprintf("SUNf",1,"TRACE RE: %f\n", traceRe);
}

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

    int return_value = 0;
    /* setup process id and communications */
    logger_map("DEBUG", "debug");

    setup_process(&argc, &argv);

    setup_gauge_fields();

    suNg_algebra_vector f[dAdj];
    suNg A, B, TMP, CAS;
    suNf a, b, tmp, cas;
    double tau, trace;
    int i, j;
#ifdef GAUGE_SUN
    lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
#elif defined(GAUGE_SPN)
    lprintf("MAIN",0,"Gauge group: SP(%d)\n",NG);
#endif
    lprintf("MAIN",0,"Fermion representation: dim = %d\n",NF);

    for (i = 0; i < dAdj; i++)
    {
        _algebra_vector_zero_g(f[i]);
        f[i].c[i] = 1.;
    }

    for (i = 0; i < dAdj; i++)
    {
        for (j = 0; j < dAdj; j++)
        {
            _algebra_represent(a, f[i]);
            _algebra_represent(b, f[j]);

            _suNf_times_suNf(tmp, a, b);
            _suNf_trace_re(trace, tmp);
            lprintf("MAIN", 0, "tr_R (T[%d] T[%d]): %.4e ", i, j, trace);
            if (i == j)
            {
                lprintf("MAIN", 0, "  [should be: %.4e]\n", -Tr);
                if(error_compare(-Tr,trace)){
                    lprintf("MAIN",0,"Matrix a:\n");
                    print_suNf(&a);
                    lprintf("MAIN",0,"Matrix tmp:\n");
                    print_suNf(&tmp);
                    return_value += 1;
                }
                else lprintf("MAIN",0,"[OK]\n");
             }
            else
            {
                lprintf("MAIN", 0, "  [should be: 0.00]\n");
                if(error_compare(0.0,trace)){
                    lprintf("MAIN",0,"Matrix tmp:\n");
                    print_suNf(&tmp);
                    return_value += 1;
                }
                else lprintf("MAIN",0,"[OK]\n");
            }

            _fund_algebra_represent(A, f[i]);
            _fund_algebra_represent(B, f[j]);

            _suNg_times_suNg(TMP, A, B);
            _suNg_trace_re(trace, TMP);
            lprintf("MAIN", 0, "tr_f (T[%d] T[%d]): %.4e ", i, j, trace);
            if (i == j)
            {
                lprintf("MAIN", 0, "  [should be: %.4e]\n", -0.5);
                if(error_compare(-0.5,trace)){
                    lprintf("MAIN",0,"Matrix A:\n");
                    print_suNg(&A);
                    lprintf("MAIN",0,"Matrix TMP:\n");
                    print_suNg(&TMP);
                    return_value += 1;
                }
                else lprintf("MAIN",0,"[OK]\n");
            }
            else
            {
                lprintf("MAIN", 0, "  [should be: 0.00]\n");
                if(error_compare(0.0,trace)){
                    lprintf("MAIN",0,"Matrix TMP:\n");
                    print_suNg(&TMP);
                    return_value += 1;
                }
                else lprintf("MAIN",0,"[OK]\n");
            }
        }
    }

    _algebra_represent(a, f[0]);
    _fund_algebra_represent(A, f[0]);
    _suNf_times_suNf(cas, a, a);
    _suNg_times_suNg(CAS, A, A);

    for (i = 1; i < dAdj; i++)
    {
        _algebra_represent(a, f[i]);
        _fund_algebra_represent(A, f[i]);

        _suNf_times_suNf(tmp, a, a);
        _suNf_add_assign(cas, tmp);
        _suNg_times_suNg(TMP, A, A);
        _suNg_add_assign(CAS, TMP);
    }

    _suNf_unit(tmp);
    _suNf_mul(tmp, C2, tmp);
    _suNf_add_assign(cas, tmp);
    _suNf_sqnorm(tau, cas);
    lprintf("MAIN", 0, "casimir check: %.4e\n", tau);
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if(error_compare(0.0,tau)){
        lprintf("MAIN", 0,"dAdj : %d\n", dAdj);
        print_suNf(&cas); //DEBUG
        lprintf("MAIN", 0,"Norm of the commutators:\n");
        for (i=0;i<dAdj;i++)
        {
            suNf commp,commm;
            _algebra_represent(a,f[i]);

            {
                double commpnorm;
                _suNf_times_suNf(commp,cas,a);
                _suNf_sqnorm(commpnorm,commp);
                lprintf("MAIN", 0,"||+||^2 = %.4f, ",commpnorm);
            }
            {
                double commmnorm;
                _suNf_times_suNf(commm,a,cas);
                _suNf_sqnorm(commmnorm,commm);
                lprintf("MAIN", 0,"||-||^2 = %.4f, ",commmnorm);
            }
            {
                double commnorm;
                _suNf_sub_assign(commp,commm);
                _suNf_sqnorm(commnorm,commp);
                lprintf("MAIN", 0,"||[cas,a(%d)]||^2 = %.4f \n",i,commnorm);
            }
        }
        return_value += 1;
    }
    else printf("[OK]\n");



    _suNg_unit(TMP);
    _suNg_mul(TMP, fund, TMP);
    _suNg_add_assign(CAS, TMP);
    _suNg_sqnorm(tau, CAS);
    lprintf("MAIN", 0, "casimir check: %.4e\n", tau);
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (fabs(tau) > 10.e-14)
        return_value += 1;

    finalize_process();
    return return_value;
}
