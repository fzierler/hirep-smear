/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
*
* File su3_utils.c
*
* Functions to project to SU(3)
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "suN.h"
#include "representation.h"
#include "logger.h"

void vector_star(suNg_vector *v1,suNg_vector *v2){
	for(int i=0;i<NG;i++){
		_complex_star(v1->c[i],v2->c[i]);
	}
}

#ifdef GAUGE_SON
static void normalize(double *v)
{
  double fact=0;
  int i;
  for (i=0;i<NG; ++i){fact+=v[i]*v[i];}
  fact = 1.0/sqrt(fact);
  for (i=0;i<NG; ++i){v[i]*=fact;}
}
#else
static void normalize(suNg_vector *v)
{
  double fact;
  _vector_prod_re_g(fact,*v,*v);
  fact=1.0/sqrt(fact);
  _vector_mul_g(*v, fact, *v);
}
#endif



static void normalize_flt(suNg_vector_flt *v)
{
  float fact;
  _vector_prod_re_g(fact,*v,*v);
  fact=1.0f/sqrtf(fact);
  _vector_mul_g(*v, fact, *v);
}


void project_to_suNg(suNg *u)
{
  double norm;
  _suNg_sqnorm(norm,*u);
  if(norm < 1.e-28) return;

#ifdef GAUGE_SON
  double *v1, *v2;
  int i,j,k;
  double z;
  for (i=0; i<NG; ++i ) {
    v2=&u->c[i*NG];
    for (j=0; j<i; ++j) {
      v1=&u->c[j*NG];
      z=0;
      for (k=0;k<NG; ++k){z+=v1[k]*v2[k];} /*_vector_prod_re_g */
      for (k=0;k<NG;++k){v2[k]-= z*v1[k];} /*_vector_project_g */
    }
    normalize(v2);
  }
#else
#ifdef WITH_QUATERNIONS

  _suNg_sqnorm(norm,*u);
  norm=sqrt(0.5*norm);
  norm=1./norm;
  _suNg_mul(*u,norm,*u);

#else
  int i,j;
  suNg_vector *v1,*v2;
  complex z;

  v1=(suNg_vector*)(u);
  v2=v1+1;
#ifdef GAUGE_SPN
  suNg Omega;
  _symplectic(Omega);
  complex z2;
  normalize(v1);

  for (i=1; i<NG/2; ++i ) {
    suNg_vector v3;
    for (j=i; j>0; --j) {
      _vector_prod_re_g(z.re,*v1, *v2);
      _vector_prod_im_g(z.im,*v1, *v2);
      _vector_project_g(*v2, z, *v1);
      _suNg_multiply(v3,Omega,*v1);
      _vector_conjugate(v3);
      _vector_prod_re_g(z2.re,v3, *v2);
      _vector_prod_im_g(z2.im,v3, *v2);
      _vector_project_g(*v2, z2, v3);
      ++v1;
    }
    normalize(v2);
    ++v2;
    v1=(suNg_vector*)(u);
  }

#else
  normalize(v1);
  for (i=1; i<NG; ++i ) {
    for (j=i; j>0; --j) {
      _vector_prod_re_g(z.re,*v1, *v2);
      _vector_prod_im_g(z.im,*v1, *v2);
      _vector_project_g(*v2, z, *v1);
      ++v1;
    }
    normalize(v2);
    ++v2;
    v1=(suNg_vector*)(u);
  }
#endif
#endif
#endif
}

void project_to_suNg_flt(suNg_flt *u)
{
  float norm;
  _suNg_sqnorm(norm,*u);
  if(norm < 1.e-10) return;

#ifdef GAUGE_SON
  float *v1, *v2;
  int i,j,k;
  float z;
  for (i=0; i<NG; ++i ) {
    v2=&u->c[i*NG];
    for (j=0; j<i; ++j) {
      v1=&u->c[j*NG];
      z=0;
      for (k=0;k<NG; ++k){z+=v1[k]*v2[k];} /*_vector_prod_re_g */
      for (k=0;k<NG;++k){v2[k]-= z*v1[k];} /*_vector_project_g */
    }
    normalize(v2);
  }
#else
#ifdef WITH_QUATERNIONS

  _suNg_sqnorm(norm,*u);
  norm=sqrtf(0.5f*norm);
  norm=1.f/norm;
  _suNg_mul(*u,norm,*u);

#else

  int i,j;
  suNg_vector_flt *v1,*v2;
  complex_flt z;

  v1=(suNg_vector_flt*)(u);
  v2=v1+1;

  normalize_flt(v1);
  for (i=1; i<NG; ++i ) {
    for (j=i; j>0; --j) {
      _vector_prod_re_g(z.re,*v1, *v2);
      _vector_prod_im_g(z.im,*v1, *v2);
      _vector_project_g(*v2, z, *v1);
      ++v1;
    }
    normalize_flt(v2);
    ++v2;
    v1=(suNg_vector_flt*)(u);
  }
#endif
#endif
}


#if !(defined GAUGE_SON || defined GAUGE_SPN)
void project_cooling_to_suNg(suNg* g_out, suNg* g_in, int cooling)
{
#ifdef WITH_QUATERNIONS
  error(1,1,"project_cooling_to_suNg " __FILE__,"not implemented with quaternions");
#else
  suNg Ug[3];
  suNg tmp[2];
  int k,l;
  int j,i,ncool;
  double c[NG];
  complex f[2];
  double norm;

  Ug[0]=*g_in;


  for (j=0; j<NG; j++){
    c[j]=0.0;

    for (i=0; i<NG; i++){
      _complex_0(f[1]);

      for (k=0; k<j; k++){
	_complex_0(f[0]);

	for (l=0; l<NG; l++){
	  _complex_mul_star_assign(f[0],(Ug[0]).c[l*NG+j], (Ug[1]).c[l*NG+k]);
	}
	_complex_mulcr_assign(f[1],c[k],(Ug[1]).c[i*NG+k],f[0]);

      }

      _complex_sub(Ug[1].c[i*NG+j],Ug[0].c[i*NG+j],f[1]);
      _complex_mul_star_assign_re(c[j],Ug[1].c[i*NG+j],Ug[1].c[i*NG+j]);

    }

    c[j]= 1.0/c[j];
  }

  for(i=0;i<NG;i++)
    {
      norm=0.0;
      for(j=0;j<NG;j++){
	_complex_mul_star_assign_re(norm,Ug[1].c[i+NG*j],Ug[1].c[i+NG*j]);
      }

      for(j=0;j<NG;j++){
	_complex_mulr(Ug[1].c[i+NG*j],1.0/sqrt(norm),Ug[1].c[i+NG*j]);

      }
    }



  _suNg_dagger(Ug[2],*g_in);

  for (ncool=0; ncool<cooling; ncool++)
    {

      _suNg_times_suNg(Ug[0],Ug[2],Ug[1]);

      for (i=0; i<NG; i++) {
	for (j=i+1; j<NG; j++) {

	  _complex_add_star(f[0],Ug[0].c[i+NG*i],Ug[0].c[j+NG*j]);
	  _complex_sub_star(f[1],Ug[0].c[j+NG*i],Ug[0].c[i+NG*j]);

	  norm = 1.0/sqrt( _complex_prod_re(f[0],f[0]) + _complex_prod_re(f[1],f[1]) );

	  _complex_mulr(f[0],norm,f[0]);
	  _complex_mulr(f[1],norm,f[1]);

	  _suNg_unit(tmp[0]);

	  _complex_star(tmp[0].c[i+NG*i],f[0]);
	  _complex_star(tmp[0].c[i+NG*j],f[1]);
	  tmp[0].c[j+NG*j]=f[0];
	  _complex_minus(tmp[0].c[j+NG*i],f[1]);



	  _suNg_times_suNg(tmp[1],Ug[1],tmp[0]);
	  Ug[1]=tmp[1];

	  _suNg_times_suNg(tmp[1],Ug[0],tmp[0]);
	  Ug[0]=tmp[1];

	}
      }
    }

  *g_out = Ug[1];
#endif
}
#endif

#ifdef GAUGE_SON

int project_to_suNg_real(suNg *out, suNg *in){
  suNg hm,om,tmp;
  double eigval[NG];
  double det;
  int i,j;
  _suNg_times_suNg_dagger(hm,*in,*in);
  diag_hmat(&hm,eigval);
  for (i=0;i<NG;++i){
    eigval[i] = 1.0/sqrt(eigval[i]);
  }
  for (i=0;i<NG;++i){
    for (j=0;j<NG;++j){
      tmp.c[i*NG+j]=eigval[j]*hm.c[i*NG+j];
    }
  }

  _suNg_times_suNg_dagger(om,tmp,hm);
  _suNg_times_suNg(tmp,om,*in);
  *out=tmp;
  //Fix the determinant
  det_suNg(&det,&tmp);
  /*  if (fabs(det)<1-1e-7 || fabs(det)>1+1e-7){
      lprintf("suNg_utils",10,"Error in project project_to_suNg_real: determinant not +/-1. It is %1.8g\n",det);
    }*/

  for (i=0;i<NG;++i){
    out->c[i]*=1./det;
  }

  tmp = *out;
  det_suNg(&det,&tmp);
  if (det<1-1e-7 || det>1+1e-7){
    lprintf("suNg_utils",10,"Error in project project_to_suNg_real: determinant not +/-1. It is %1.8g.",det);
    return 0;
  }
  return 1;
}
#endif

void cooling_SPN(suNg* g_out, suNg* g_in, suNg* g_tilde, int cooling){
    
    suNgfull B, U_tilde, S;
    _suNg_expand( U_tilde, *g_tilde);
    _suNg_expand( S, *g_in);
    
    //_suNffull_unit(S);
    
    /*
    lprintf("COOLING",0," U_tilde= \n",0);
    for (int i=0;i<NG;i++){
        for (int j=0;j<NG;j++){
        
            lprintf("COOLING",0," (%f,%f), ",U_tilde.c[j + i*NG].re, U_tilde.c[j + i*NG].im);
        }lprintf("COOLING",0," \n");
    }
    
    lprintf("COOLING",0,"  U= \n",0);
    for (int i=0;i<NG;i++){
        for (int j=0;j<NG;j++){
        
            lprintf("COOLING",0," (%f,%f), ",S.c[j + i*NG].re, S.c[j + i*NG].im);
        }lprintf("COOLING",0," \n");
    }
    _suNgfull_times_suNgfull_dagger(B, S, U_tilde);
    
    lprintf("COOLING",0,"  U*U_tilde^+= \n",0);
    for (int i=0;i<NG;i++){
        for (int j=0;j<NG;j++){
        
            lprintf("COOLING",0," (%f,%f), ",B.c[j + i*NG].re, B.c[j + i*NG].im);
        }lprintf("COOLING",0," \n");
    }
    */
    
    for (int nbcool=0;nbcool<cooling;nbcool++){
        _suNgfull_times_suNgfull_dagger(B, S, U_tilde);
        
        for (int N1=0;N1<NG/2-1;N1++){
            for (int N2=N1+1;N2<NG/2;N2++){
                subgrb(N1, N2, &B, &S);
                subgrb(NG/2 + N1, NG/2 + N2, &B, &S);
                subgrb_tau( N1, N2, &B, &S);
                subgrb(N1, NG/2 + N1, &B, &S);
                subgrb(N2, NG/2 + N2, &B, &S);
            }
        }
        
        /*
        lprintf("COOLING",0," U*U_tilde^+ cooled n=%d \n",nbcool+1);
        for (int i=0;i<NG;i++){
            for (int j=0;j<NG;j++){
            
                lprintf("COOLING",0," (%f,%f), ",B.c[j + i*NG].re, B.c[i + j*NG].im);
            }lprintf("COOLING",0," \n");
        }

        
        lprintf("COOLING",0," U cooled n=%d \n",nbcool+1);
        for (int i=0;i<NG;i++){
            for (int j=0;j<NG;j++){
            
                lprintf("COOLING",0," (%f,%f), ",S.c[j + i*NG].re, S.c[i + j*NG].im);
            }lprintf("COOLING",0," \n");
        }
        */
        
    }
    
    for (int i=0;i<(NG/2)*NG;i++){
        
        g_out->c[i].re = S.c[i].re;
        g_out->c[i].im = S.c[i].im;
    }
}

void subgrb(int i1col, int i2col, suNgfull* B11, suNgfull* C11){
    
    complex F11, F12, A11[4], ztmp1, ztmp2;    // A11 is the extracted su2 matrix in 1-d array
    double UMAG;                               // [ 0 1 ]
    int i1,i2,i3,i4;                           // [ 2 3 ]   --> [0, 1, 2, 3]
    
    i1 = i1col + i1col*NG;
    i2 = i2col + i2col*NG;
    i3 = i2col + i1col*NG;
    i4 = i1col + i2col*NG;
    
    _complex_add_star(F11, B11->c[i1], B11->c[i2]);
    _complex_mulr(F11, 0.5, F11);
    _complex_sub_star(F12, B11->c[i3], B11->c[i4]);
    _complex_mulr(F12, 0.5, F12);
    
    
    _complex_mul_star(ztmp1, F11, F11);
    _complex_mul_star(ztmp2, F12, F12);
    UMAG = sqrt(ztmp1.re+ztmp2.re);
    UMAG = 1./UMAG;
    
    _complex_mulr(F11, UMAG, F11);
    _complex_mulr(F12, UMAG, F12);
    
    _complex_star(A11[0], F11);
    _complex_star(A11[1], F12);
    _complex_mulr(A11[2], -1., F12);
    _complex_mulr(A11[3],  1., F11);
    
    vmxsu2(i1col,i2col,C11,A11);
    vmxsu2(i1col,i2col,B11,A11);
}

void vmxsu2(int i1, int i2, suNgfull* A, complex B[4]){
    
    complex C[NG*2];
    complex ztmp1,ztmp2;
    
    for (int i=0; i<NG; i++){
        _complex_mul(ztmp1, A->c[i + i1*NG], B[0]);
        _complex_mul(ztmp2, A->c[i + i2*NG], B[2]);
        _complex_add(C[i], ztmp1, ztmp2);
        
        _complex_mul(ztmp1, A->c[i + i1*NG], B[1]);
        _complex_mul(ztmp2, A->c[i + i2*NG], B[3]);
        _complex_add(C[i+NG], ztmp1, ztmp2);
    }
    
    for (int i=0; i<NG; i++){
        A->c[i + i1*NG] = C[i];
        A->c[i + i2*NG] = C[i+NG];
    }
}

void subgrb_tau(int n1, int n2, suNgfull* B11, suNgfull* C11){
    
    complex F11, F12, A11[4], ztmp1, ztmp2;    // A11 is the extracted su2 matrix in 1-d array
    double UMAG;                               // [ 0 1 ]
    int i1,i2,i3,i4;                           // [ 2 3 ]   --> [0, 1, 2, 3]
    
    i1 = n1 + n1*NG;
    i2 = (NG/2) + n2 + (NG/2 + n2)*NG;
    i3 = NG/2 + n2 + n1*NG;
    i4 = (NG/2 + n2)*NG + n1;
    
    _complex_add_star(F11, B11->c[i1], B11->c[i2]);
    _complex_mulr(F11, 0.5, F11);
    _complex_sub_star(F12, B11->c[i3], B11->c[i4]);
    _complex_mulr(F12, 0.5, F12);
    
    _complex_mul_star(ztmp1, F11, F11);
    _complex_mul_star(ztmp2, F12, F12);
    UMAG = sqrt(ztmp1.re+ztmp2.re);
    UMAG = 1./UMAG;
    
    _complex_mulr(F11, UMAG, F11);
    _complex_mulr(F12, UMAG, F12);
    
    _complex_star(A11[0], F11);
    _complex_star(A11[1], F12);
    _complex_mulr(A11[2], -1., F12);
    _complex_mulr(A11[3],  1., F11);
    
    vmxsu2_tau(n1, n2, C11, A11);
    vmxsu2_tau(n1, n2, B11, A11);
}

void vmxsu2_tau(int n1, int n2, suNgfull* A, complex B[4]){
    
    suNgfull C;
    complex ztmp1, ztmp2;
    
    _suNgfull_mul(C, 1., *A);
    
    for (int i=0; i<NG; i++){
        
        _complex_mul(ztmp1, A->c[i + n1*NG], B[0]);
        _complex_mul(ztmp2, A->c[i +(NG/2+n2)*NG], B[2]);
        _complex_add(C.c[i + n1*NG], ztmp1, ztmp2);
        
        _complex_mul(ztmp1, A->c[i + n2*NG], B[0]);
        _complex_mul(ztmp2, A->c[i + (NG/2+n1)*NG], B[2]);
        _complex_add(C.c[i + n2*NG], ztmp1, ztmp2);
        
        _complex_mul(ztmp1, A->c[i + n2*NG], B[1]);
        _complex_mul(ztmp2, A->c[i + (NG/2+n1)*NG], B[3]);
        _complex_add(C.c[i + (NG/2 + n1)*NG], ztmp1, ztmp2);
        
        _complex_mul(ztmp1, A->c[i + n1*NG], B[1]);
        _complex_mul(ztmp2, A->c[i + (NG/2+n2)*NG], B[3]);
        _complex_add(C.c[i + (NG/2 + n2)*NG], ztmp1, ztmp2);
    }
    _suNgfull_mul(*A, 1., C);
}
