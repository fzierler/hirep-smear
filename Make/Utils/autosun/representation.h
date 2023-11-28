#ifndef REPRESENTATION_H
#define REPRESENTATION_H
#include <string>
#include "./complex.h"
#include "./matrix.h"
#include "./sun.h"
#ifdef _REPR_FUNDAMENTAL_
#include "fundamental.h"
#elif _REPR_ADJOINT_
#include "adjoint.h"
#elif _REPR_ANTISYMMETRIC_
#include "antisymmetric.h"
#elif _REPR_SYMMETRIC_
#include "symmetric.h"
#endif


using namespace std;

string group_represent(const char* vname, const char* uname);
string algebra_represent(const char* mname, const char* hname);
string algebra_project(const char* hname, const char* mname);


string algebra_represent(const char* mname, const char* hname)
{
	string RET;
	rvector H(group::DIM,hname);
	pmatrix M(representation::DIM);
	
	for(int A = 0; A < group::DIM; A++)
	{
		pmatrix iT(representation::iT[A]);
		iT.scale(H.get(A));
		M.add(iT);
	}
	
#if defined( _GAUGE_SPN_) && defined ( _REPR_FUNDAMENTAL_ )
	RET = M.symplectic_compressed_assignment("=", mname);
#else
	RET = M.assignment("=", mname);
#endif

	return RET;
}


string algebra_project(const char* hname, const char* mname)
{
	string RET;
	pvector H(group::DIM);
	pmatrix *M;
	//pmatrix adjM(representation::DIM);

	if(sizeof(representation::TYPE)==sizeof(FLOATING))
		M = new rmatrix(representation::DIM,mname);
	else
#ifdef _GAUGE_SON_
	   M = new rmatrix(representation::DIM,mname);
#else
		M = new cmatrix(representation::DIM,mname);
#endif

	//adjM = *M;
	//adjM.adjoint();
	//M->sub(adjM);
	
	for(int A = 0; A < group::DIM; A++)
	{
		pmatrix iT(representation::iT[A]);
		polynomial iTM;
		herm(iTM,iT,*M);
		iTM.real();
		iTM.scale(1.0/representation::iTnorm);
		H.set(A, iTM);
	}
	
	RET = H.assignment("=", hname);

	delete M;
	
	return RET;
}
#endif
