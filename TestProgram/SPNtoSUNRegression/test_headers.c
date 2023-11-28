#include "suN.h"
#include "suN_types.h"
#include "SP.h"
#include "SP_types.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

/* Print a vector */
void print_vector( suNg_vector vector ){
  printf("\n");
  for( int i=0; i<NG; i++ ){
    printf("(%6.2f,%6.2f) ",vector.c[i].re,vector.c[i].im);
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


/* Compare a two vectors
 */
void compare_vectors( suNg_vector v1, suNg_vector v2 ){
  
  double sum = 0;
  for(int i=0; i<NG;i++){
    double diff = v1.c[i].re - v2.c[i].re;
    sum += diff*diff;
    diff = v2.c[i].im - v2.c[i].im;
    sum += diff*diff;
  }
  
  if( sum < 1e-14 ){
    printf("PASSED, diff=%g\n\n",sum);
  } else {
    printf("TEST FAILED, diff=%g\n",sum);
    print_vector(v1);
    print_vector(v2);
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
/* Create a vector */
suNg_vector random_vector(){
  suNg_vector vector;

  /* First half of the rows */
  for( int i=0; i<NG; i++ ){
    vector.c[i].re = random()*10./RAND_MAX;
    vector.c[i].im = random()*10./RAND_MAX;
  }
  return vector;
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
 
  suNg suNmatrix, suNmatrix2,suNresult;
  SPg SPmatrix, SPmatrix2, SPresult;
  suNg_vector v1, v2, v3;
  
  printf("Creating the matrices and testing similarity\n");
  SPmatrix = random_SPg();
  SPmatrix2 = random_SPg();
  suNmatrix = SPg_to_suNg( SPmatrix );
  suNmatrix2 = SPg_to_suNg( SPmatrix2 );
  compare_suNg_SPg( suNmatrix, SPmatrix );
  
  /* Pull a random vector. No difference between compressed and full. */
  v1 = random_vector();
  
  printf("M_multiply: \n");
  _SPg_multiply( v2, SPmatrix, v1 );
  _suNg_multiply( v3, suNmatrix, v1 );
  compare_vectors( v2, v3 );
  
  printf("M_inverse_multiply: \n");
  _SPg_inverse_multiply( v2, SPmatrix, v1 );
  _suNg_inverse_multiply( v3, suNmatrix, v1 );
  compare_vectors( v2, v3 );

  printf("M_dagger: \n");
  _suNg_dagger( suNresult, suNmatrix );
  _SPg_dagger( SPresult, SPmatrix );
  compare_suNg_SPg( suNresult, SPresult );
 
  
  printf("M_times_M: \n");
  _suNg_times_suNg( suNresult, suNmatrix, suNmatrix2 );
  _SPg_times_SPg( SPresult, SPmatrix, SPmatrix2 );
  compare_suNg_SPg( suNresult, SPresult );

  printf("M_times_Mdagger: \n");
  _suNg_times_suNg_dagger( suNresult, suNmatrix, suNmatrix2 );
  _SPg_times_SPg_dagger( SPresult, SPmatrix, SPmatrix2 );
  compare_suNg_SPg( suNresult, SPresult );
 
  printf("Mdagger_times_M: \n");
  _suNg_dagger_times_suNg( suNresult, suNmatrix, suNmatrix2 );
  _SPg_dagger_times_SPg( SPresult, SPmatrix, SPmatrix2 );
  compare_suNg_SPg( suNresult, SPresult );


  printf("zero matrix: \n");
  _suNg_zero(suNresult);
  _SPg_zero(SPresult);
  compare_suNg_SPg( suNresult, SPresult );
 
  printf("unit matrix: \n");
  _suNg_unit(suNresult);
  _SPg_unit(SPresult);
  compare_suNg_SPg( suNresult, SPresult );

  printf("M1 = - M2 : \n");
  _suNg_minus(suNresult, suNmatrix);
  _SPg_minus(SPresult, SPmatrix);
  compare_suNg_SPg( suNresult, SPresult );

  printf("M1 = REAL*M2 : \n");
  _suNg_mul(suNresult,-0.157, suNmatrix);
  _SPg_mul(SPresult,-0.157, SPmatrix);
  compare_suNg_SPg( suNresult, SPresult );

  printf("M1 = COMPLEX*M2 : \n");
  complex p = {.re = -0.1245, .im = 1.2433};
  _suNg_mulc(suNresult,p, suNmatrix);
  _SPg_mulc(SPresult,p, SPmatrix);
  compare_suNg_SPg( suNresult, SPresult );

  printf("M1+=M2: \n");
  _suNg_add_assign( suNresult,suNmatrix2 );
   _SPg_add_assign( SPresult, SPmatrix2 );
  compare_suNg_SPg( suNresult, SPresult );
 
  printf("M1-=M2: \n");
  _suNg_sub_assign( suNresult,suNmatrix2 );
   _SPg_sub_assign( SPresult, SPmatrix2 );
  compare_suNg_SPg( suNresult, SPresult );
 
  printf("l2norm: \n");
  double ksp,ksu;
  _suNg_sqnorm( ksu,suNmatrix);
  _SPg_sqnorm( ksp, SPmatrix);
  if(fabs((ksp-ksu)/ksu) > 1e-14){
      printf("SuN : %f , SPN: %f\n", ksu,ksp);
      exit(1);
  }else printf("PASSED, diff=%g\n\n",fabs((ksp-ksu)/ksu));
 
  printf("(M-1) l2norm: \n");
  _suNg_sqnorm_m1( ksu,suNmatrix);
  _SPg_sqnorm_m1( ksp, SPmatrix);
  if(fabs((ksp-ksu)/ksu) > 1e-14){
      printf("SuN : %f , SPN: %f\n", ksu,ksp);
      exit(1);
  }else printf("PASSED, diff=%g\n\n",fabs((ksp-ksu)/ksu));
 
  printf("trace re: \n");
  _suNg_trace_re( ksu,suNmatrix);
  _SPg_trace_re( ksp, SPmatrix);
  if(fabs((ksp-ksu)/ksu) > 1e-14){
      printf("SuN : %f , SPN: %f\n", ksu,ksp);
      exit(1);
  }else printf("PASSED, diff=%g\n\n",fabs((ksp-ksu)/ksu));

  printf("trace im: \n");
  _suNg_trace_im( ksu,suNmatrix);
  _SPg_trace_im( ksp, SPmatrix);
  if(fabs((ksp-ksu)) > 2e-14){
      printf("SuN : %g , SPN: %G\n", ksu,ksp);
      exit(1);
  }else printf("PASSED, diff=%g\n\n",fabs(ksp-ksu));


  return 0;
}
















