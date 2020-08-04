/***************************************************************************\
* Copyright (c) 2011, Agostino Patella                                      *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File gaugefix.h
* 
*******************************************************************************/

#ifndef GAUGEFIX_H
#define GAUGEFIX_H

#include "suN_types.h"
#include "hr_complex.h"
#include "spinor_field.h"
#include "suN.h"



void unit_gauge(suNg_field *original);

void random_gauge_transform(suNg_field *gauge);

void reunit(suNg_field *fixed_gauge);

double calc_plaq(suNg_field* V);

double gaugefix_action(int fix_dir, suNg_field *gauge );
#ifndef GAUGE_SPN // FIXFORSPN
void su2_hit(int fix_dir, int parity, double overrelax, suNg_field *fixed_gauge, int c );
double gaugefixstep(int fix_dir,double overrelax, suNg_field *fixed_gauge );
double gaugefix(int fix_dir,double overrelax,int max_it,
	      double fix_tol, suNg_field *fixed_gauge );
#endif

#endif /* GAUGEFIX_H */

