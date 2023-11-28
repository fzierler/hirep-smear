/*******************************************************************************
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
#include <assert.h>

#ifdef GAUGE_SON
#error "Test not implemented for SON."
#endif

#ifdef REPR_FUNDAMENTAL
#define COMPLEX_REP
#ifdef GAUGE_SPN
static float C2=(float)NG*(NG+1)/(4*NG);
#else
static float C2=(float)(NG*NG-1)/(2*NG);
#endif
static float Tr=0.5;
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
#else
static float C2=(float)(NG+2)*(NG-1)/NG;
static float Tr=(float)(NG+2)/2;
#endif
#endif

#ifdef WITH_MPI
#error: check_algebra_1 only works only on serial jobs
#endif

#ifdef GAUGE_SUN
static int dAdj=NG*NG-1;
static float fund=(float)(NG*NG-1)/(2*NG);
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

    printf("\n");
    for(i=0;i<NF;++i){
        for(j=0;j<NF;++j){
#ifdef COMPLEX_REP
            printf("%+-.3f,%+-.3f ", m->c[i*NF+j].re, m->c[i*NF+j].im); 
            //printf("%+-.3f ", m->c[i*NF+j].im); 
            sumsq += m->c[i*NF+j].im*m->c[i*NF+j].im;
            sumsq += m->c[i*NF+j].re*m->c[i*NF+j].re;
            if(i==j) traceRe +=m->c[i*NF+j].re;
#else
            printf("%+-.3f,%+-.3f ", m->c[i*NF+j], m->c[i*NF+j]); 
            //printf("%+-.3f ", m->c[i*NF+j].im); 
            sumsq += m->c[i*NF+j]*m->c[i*NF+j];
            if(i==j) traceRe +=m->c[i*NF+j];
#endif
        }
        printf("\n\n\n");
    }
    printf("\n");
    printf("SUM SQUARES: %f\n", sumsq);
    printf("TRACE RE: %f\n", traceRe);
}

void print_suNg(suNg *m){
    int i,j;
    double sumsq = 0;
    double traceRe = 0;

    printf("\n");
    for(i=0;i<NG;++i){
        for(j=0;j<NG;++j){
            printf("%+-.3f,%+-.3f ", m->c[i*NG+j].re, m->c[i*NG+j].im); 
            //printf("%+-.3f ", m->c[i*NF+j].im); 
            sumsq += m->c[i*NG+j].im*m->c[i*NG+j].im;
            sumsq += m->c[i*NG+j].re*m->c[i*NG+j].re;
            if(i==j) traceRe +=m->c[i*NG+j].re;
        }
        printf("\n\n\n");
    }
    printf("\n");
    printf("SUM SQUARES: %f\n", sumsq);
    printf("TRACE RE: %f\n", traceRe);
}
int main(int argc,char *argv[])
{
   suNg_algebra_vector f[dAdj];
   suNg A,B,TMP,CAS;
   suNf a,b,tmp,cas;
   double tau,trace;
   int i,j;
   
   
#ifdef GAUGE_SUN
   printf("Gauge group: SU(%d)\n",NG);
#elif defined(GAUGE_SPN)
   printf("Gauge group: SP(%d)\n",NG);
#endif
   printf("Fermion representation: dim = %d\n",NF);
   printf("\n");

   for (i=0;i<dAdj;i++)
   {
      _algebra_vector_zero_g(f[i]);
      f[i].c[i]=1.;
   }

   for (i=0;i<dAdj;i++)
   {
      for (j=0;j<dAdj;j++)
      {
         _algebra_represent(a,f[i]);
         _algebra_represent(b,f[j]);
         
         _suNf_times_suNf(tmp,a,b);
         _suNf_trace_re(trace,tmp);
         printf("tr_R (T[%d] T[%d]): %+-.4f ", i, j,trace);
         if (i==j){
            printf("  [should be: %+-.4f]",-Tr);
            if(error_compare(-Tr,trace)){
                printf("Matrix a:\n");
                print_suNf(&a);
                printf("Matrix tmp:\n");
                print_suNf(&tmp);
                assert(0);
            }
            else printf("[OK]\n");
         }
         else{
            printf("  [should be: 0.00]");
            if(error_compare(0.0,trace)){
                printf("Matrix tmp:\n");
                print_suNf(&tmp);
                assert(0);
            }
            else printf("[OK]\n");
         }

         _fund_algebra_represent(A,f[i]);
         _fund_algebra_represent(B,f[j]);
         
         _suNg_times_suNg(TMP,A,B);
         _suNg_trace_re(trace,TMP);
         printf("tr_f (T[%d] T[%d]): %.4f ", i, j,trace);
         if (i==j){
            printf("  [should be: %.4f]",-0.5);
            if(error_compare(-0.5,trace)){
                printf("Matrix A:\n");
                print_suNg(&A);
                printf("Matrix TMP:\n");
                print_suNg(&TMP);
                assert(0);
            }
            else printf("[OK]\n");
         }
         else{
            printf("  [should be: 0.00]");
            if(error_compare(0.0,trace)){
                printf("Matrix TMP:\n");
                print_suNg(&TMP);
                assert(0);
            }
            else printf("[OK]\n");
         }
      }
   }
   
	 _algebra_represent(a,f[0]);
	 _suNf_times_suNf(cas,a,a);

	 _fund_algebra_represent(A,f[0]);
	 _suNg_times_suNg(CAS,A,A);
         
   for (i=1;i<dAdj;i++)
   {
      _algebra_represent(a,f[i]);
      _suNf_times_suNf(tmp,a,a);
      _suNf_add_assign(cas,tmp);


      _fund_algebra_represent(A,f[i]);
      _suNg_times_suNg(TMP,A,A);
      _suNg_add_assign(CAS,TMP);
   }

   _suNf_unit(tmp);
   _suNf_mul(tmp,C2,tmp);
   _suNf_add_assign(cas,tmp);
   _suNf_sqnorm(tau,cas);
   printf("casimir check: %.3f ",tau);
   printf("(should be 0.00)");
   if(error_compare(0.0,tau)){
       printf("dAdj : %d\n", dAdj);
       print_suNf(&cas); //DEBUG
       printf("Norm of the commutators:\n");
       for (i=0;i<dAdj;i++)
       {
           suNf commp,commm;
           _algebra_represent(a,f[i]);

           {
               double commpnorm;
               _suNf_times_suNf(commp,cas,a);
               _suNf_sqnorm(commpnorm,commp);
               printf("||+||^2 = %.4f, ",commpnorm);
           }
           {
               double commmnorm;
               _suNf_times_suNf(commm,a,cas);
               _suNf_sqnorm(commmnorm,commm);
               printf("||-||^2 = %.4f, ",commmnorm);
           }
           {
               double commnorm;
               _suNf_sub_assign(commp,commm);
               _suNf_sqnorm(commnorm,commp);
               printf("||[cas,a(%d)]||^2 = %.4f \n",i,commnorm);
           }
       }



       assert(0);
   }
   else printf("[OK]\n");

   _suNg_unit(TMP);
   _suNg_mul(TMP,fund,TMP);
   _suNg_add_assign(CAS,TMP);
   _suNg_sqnorm(tau,CAS);
   printf("casimir check: %.3f ",tau);
   printf("(should be 0.00)");
   if(error_compare(0.0,tau)){
       assert(0);
   }
   else printf("[OK]\n");
   return 0;
}
