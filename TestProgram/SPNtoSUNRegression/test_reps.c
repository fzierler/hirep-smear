#include "suN.h"
#include "suN_types.h"
#include "SP.h"
#include "SP_types.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "spN_func.h"
#include "suN_func.h"
#include "spn_sun_algconv.h"
#include "suN_exp.h"
#include "spN_exp.h"

/* Printing for debugging:
 * Print a symplectic matrix */
void print_SPg( SPg SPmatrix ){
  printf("\n");
  for( int i=0; i<NG*NG/2; i++ ){
    printf("(%6.2f,%6.2f) ",SPmatrix.c[i].re,SPmatrix.c[i].im);
    if( (i+1)%NG == 0 ) printf("\n");
  }
  printf("\n");
}

/* Print a full matrix */
void print_suNg( suNg suNmatrix ){
  printf("\n");
  for( int i=0; i<NG*NG; i++ ){
    printf("(%6.2f,%6.2f) ",suNmatrix.c[i].re,suNmatrix.c[i].im);
    if( (i+1)%NG == 0 ) printf("\n");
  }
  printf("\n");
}


/* Print a represented matrix */
#ifdef REPR_ADJOINT //real matrices
void print_SPf( SPf spNadj ){
  printf("\n");
  for( int i=0; i<SPNF*SPNF; i++ ){
    printf("%6.2f ",spNadj.c[i]);
    if( (i+1)%SPNF == 0 ) printf("\n");
  }
  printf("\n");
}

void print_suNf( suNf suNadj ){
  printf("\n");
  for( int i=0; i<NF*NF; i++ ){
    printf("%6.2f ",suNadj.c[i]);
    if( (i+1)%NF == 0 ) printf("\n");
  }
  printf("\n");
}

#else

void print_SPf( SPf spNadj ){
  printf("\n");
  for( int i=0; i<SPNF*SPNF; i++ ){
    printf("(%6.2f,%6.2f) ",spNadj.c[i].re,spNadj.c[i].im);
    if( (i+1)%SPNF == 0 ) printf("\n");
  }
  printf("\n");
}

void print_suNf( suNf suNadj ){
  printf("\n");
  for( int i=0; i<NF*NF; i++ ){
    printf("(%6.2f,%6.2f) ",suNadj.c[i].re,suNadj.c[i].im);
    if( (i+1)%NF == 0 ) printf("\n");
  }
  printf("\n");
}

#endif


void print_spNalg( SPg_algebra_vector spNalg ){
  printf("\n");
  for( int i=0; i<NG*(NG+1)/2; i++ ){
    printf("%6.2f ",spNalg.c[i]);
  }
  printf("\n");
}

void print_suNalg( suNg_algebra_vector suNalg ){
  printf("\n");
  for( int i=0; i<NG*NG-1; i++ ){
    printf("%6.2f ",suNalg.c[i]);
  }
  printf("\n");
}



/* Compare a compressed symplectic and a full NG*NG matrix
 */
void compare_suNg_SPg( suNg suNmatrix, SPg SPmatrix ){
  
  double sum = 0;
  for(int i=0; i<NG*NG/2;i++){
    double diff = suNmatrix.c[i].re - SPmatrix.c[i].re;
    sum += diff*diff;
    diff = suNmatrix.c[i].im - SPmatrix.c[i].im;
    sum += diff*diff;
  }

  if( sum < 1e-14 ){
    printf("PASSED, diff=%g\n\n",sum);
  } else {
    printf("TEST FAILED, diff=%g\n",sum);
    print_SPg(SPmatrix);
    print_suNg(suNmatrix);
    
    exit(1);
  }
}

/* Compare represented matrices
 */
void compare_represented( SPf spNadj, suNf suNadj){
  double sum = 0;
  
#ifdef REPR_ADJOINT
  SPf tmp;
  _suntospn_adj(tmp.c,suNadj.c);
 for(int i=0; i<SPNF*SPNF;i++){
    double diff = tmp.c[i] - spNadj.c[i];
    sum += diff*diff;
  }
  
#else

#ifdef  REPR_FUNDAMENTAL  // using fundamental representation, fermion matrix is symplectic
  for(int i=0; i<NF*NF/2;i++){
#else 
  for(int i=0; i<NF*NF;i++){
#endif
    double diff = suNadj.c[i].re - spNadj.c[i].re;
    sum += diff*diff;
    diff = suNadj.c[i].im - spNadj.c[i].im;
    sum += diff*diff;
  }
  
#endif

  if( sum < 1e-14 ){
    printf("PASSED, diff=%g\n\n",sum);
  } else {
    printf("TEST FAILED, diff=%g\n",sum);
#ifdef REPR_ADJOINT
    print_SPf(tmp);
#endif
    print_SPf(spNadj);
    print_suNf(suNadj);

    exit(1);
  }
}



/* Compare represented matrices
 */
void compare_algebra( SPg_algebra_vector spNalg, suNg_algebra_vector suNalg){
  double sum = 0;
  
  SPg_algebra_vector spNalg_repr;
  _suntospn_algebra(spNalg_repr.c,suNalg.c);
  
  for(int i=0; i<SPNF;i++){
    double diff = spNalg_repr.c[i] - spNalg.c[i];
    sum += diff*diff;
  }

  if( sum < 1e-14 ){
    printf("PASSED, diff=%g\n\n",sum);
  } else {
    printf("TEST FAILED, diff=%g\n",sum);
    print_spNalg(spNalg_repr);
    print_spNalg(spNalg);
    print_suNalg(suNalg);

    exit(1);
  }
}



/* Create a symplectic matrix */
SPg random_SPg(){
  SPg SPmatrix;

  /* First half of the rows */
  for( int i=0; i<NG*NG/2; i++ ){
    SPmatrix.c[i].re = random()*10./RAND_MAX;
    SPmatrix.c[i].im = random()*10./RAND_MAX;
  }
  
  return SPmatrix;
}

/* Create a fermion matrix */
void random_SPf(SPf_FMAT* fmatrix){
#ifdef REPR_ADJOINT //real matrices
  for( int i=0; i<SPNF*SPNF; i++ ){
    fmatrix->c[i] = random()*10./RAND_MAX;
  }
#else
  for( int i=0; i<SPNF*SPNF; i++ ){
    fmatrix->c[i].re = random()*10./RAND_MAX;
    fmatrix->c[i].im = random()*10./RAND_MAX;
  }
#endif
}



/* Expand the symplectic matrix into a full NG*NG matrix */
suNg SPg_to_suNg( SPg SPmatrix ){
  suNg suNmatrix;
  for( int i=0; i<NG*NG/2; i++ ){
    suNmatrix.c[i].re = SPmatrix.c[i].re;
    suNmatrix.c[i].im = SPmatrix.c[i].im;    
  }
    /* Second half is constructed out of the first half */
  for( int i=0; i<NG/2; i++ ) for( int j=0; j<NG/2; j++ ){
    suNmatrix.c[NG*NG/2+NG*j+i].re = -SPmatrix.c[NG/2+NG*j+i].re;
    suNmatrix.c[NG*NG/2+NG*j+i].im =  SPmatrix.c[NG/2+NG*j+i].im;
    suNmatrix.c[NG*NG/2+NG*j+NG/2+i].re =  SPmatrix.c[NG*j+i].re;
    suNmatrix.c[NG*NG/2+NG*j+NG/2+i].im = -SPmatrix.c[NG*j+i].im;
  }
  return suNmatrix;
}



/* Compare each function that uses compressed symplectic matrices
   to the original functions using NG*NG matrices. */
int main(void){
 
#if NG==2 && ( defined( REPR_SYMMETRIC ) || defined( REPR_ANTISYMMETRIC ))

  printf("\n");
  printf("   In the case N=2 only the adjoint representation is properly implemented.\n");
  printf("   Antisymmetric is trivial and (as for all even N) symmetric is equivalent to adjoint.\n");
  
#else

  suNg suNmatrix, suNresult;
  SPg SPmatrix, SPresult;
  suNg_vector v1, v2, v3;
  
  printf("Creating the matrices and testing similarity\n");
  SPmatrix = random_SPg();
  suNmatrix = SPg_to_suNg( SPmatrix );
  compare_suNg_SPg( suNmatrix, SPmatrix );
 
   
  printf("Testing fund_algebra project\n");
  SPg_algebra_vector spNalg;
  suNg_algebra_vector suNalg;
  _fund_algebra_project( suNalg, suNmatrix );
  _fund_spN_algebra_project( spNalg, SPmatrix );
  compare_algebra( spNalg, suNalg );
  
  /* For this test we need to project an matrix in the symplectic adjoint
     into the suN adjoint. The same member of the group is not represented
     by the same matrix. */
  printf("Testing algebra project\n");
  SPf_FMAT* sp_fermion_matrix = (SPf_FMAT*) malloc(sizeof(SPf_FMAT));
  suNf* su_fermion_matrix = (suNf*) malloc(sizeof(suNf));
  random_SPf(sp_fermion_matrix);
#ifdef REPR_ADJOINT
  _spntosun_adj(su_fermion_matrix->c,sp_fermion_matrix->c);
#else
  for(int i=0;i<NF*NF;++i){
    su_fermion_matrix->c[i].re = sp_fermion_matrix->c[i].re;
    su_fermion_matrix->c[i].im = sp_fermion_matrix->c[i].im;
  }
#endif
  _algebra_project( suNalg, *su_fermion_matrix );
  _spN_algebra_project( spNalg, *sp_fermion_matrix );
  compare_algebra( spNalg, suNalg );
  
  
  printf("Testing group represent\n");
  SPf* sp_comp_matrix = (SPf*) malloc(sizeof(SPf));
  _group_represent( *su_fermion_matrix, suNmatrix );
  _spN_group_represent( *sp_comp_matrix, SPmatrix );
  compare_represented( *sp_comp_matrix, *su_fermion_matrix );
  
  printf("Testing algebra represent\n");
  _spntosun_algebra(suNalg.c,spNalg.c);
  _algebra_represent( *su_fermion_matrix, suNalg );
  _spN_algebra_represent( *sp_comp_matrix, spNalg);
  compare_represented( *sp_comp_matrix, *su_fermion_matrix );
  
  printf("Testing fundamental algebra represent\n");
  _fund_algebra_represent( suNmatrix, suNalg );
  _fund_spN_algebra_represent( SPmatrix, spNalg);
  compare_suNg_SPg( suNmatrix, SPmatrix );
  
  printf("Testing ExpX - only at second order\n");
  ExpX( 2.0e-4, &suNalg, &suNmatrix );
  SP_ExpX( 2.0e-4, &spNalg, &SPmatrix);
  compare_suNg_SPg( suNmatrix, SPmatrix );
  
  free( sp_fermion_matrix );
  free( su_fermion_matrix );
  free( sp_comp_matrix );
#endif
  return 0;
}
















