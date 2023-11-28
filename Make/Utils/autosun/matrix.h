#ifndef MATRIX_H
#define MATRIX_H
#include "./complex.h"
#include "./sparse.h"
#include "./polynomial.h"
class smatrix : public sparsematrix<complex>
{
public:
	smatrix() : sparsematrix<complex>() {}
	smatrix(int N) : sparsematrix<complex>(N) {}
	
	using sparsematrix<complex>::operator=;

	complex& getzero() const
	{
		static complex zero(0.0,0.0);
		return zero;
	}
};


class svector : public sparsevector<complex>
{
public:
	svector() : sparsevector<complex>() {}
	svector(int N) : sparsevector<complex>(N) {}
	
	using sparsevector<complex>::operator=;

	complex& getzero() const
	{
		static complex zero(0.0,0.0);
		return zero;
	}
};


class pmatrix : public sparsematrix<polynomial>
{
public:
	pmatrix() : sparsematrix<polynomial>() {}
	pmatrix(int N) : sparsematrix<polynomial>(N) {}
	pmatrix(smatrix& mat) : sparsematrix<polynomial>(mat.size)
	{
		static rmonomial index;
		for(int i = 0; i < mat.size; i++)
			for(int j = 0; j < mat.size; j++)
			{
				polynomial tmp(index, mat.get(i,j));
				set(i,j, tmp);
			}
	}
	
	using sparsematrix<polynomial>::operator=;

	polynomial& getzero() const
	{
		static polynomial zero;
		return zero;
	}
	
	bool isreal() const
	{
		bool RET = true;
		for(KEYTYPE i = 0; i < length; i++)
			RET = RET && data[i]->value.isreal();
		return RET;
	}
	
	string assignment(const char* op, const char* name) const
	{
		ostringstream RET;
		string tmp;
		if(isreal())
		{
			for(int i = 0; i < size; i++)
				for(int j = 0; j < size; j++)
				{
					tmp = get(i,j).str_real();
					RET << _INDENT_;
					RET << name << mindex(i,j,size) << " " << op << " " << tmp;
					if(i == size-1 && j == size-1) RET << _LASTENDL_;
					else RET << _ENDL_;
				}
		}
		else
		{
			for(int i = 0; i < size; i++)
				for(int j = 0; j < size; j++)
				{
					tmp = get(i,j).str_real();
					RET << _INDENT_;
					RET << name << mindex(i,j,size) << ".re " << op << " " << tmp;
					RET << _ENDL_;

					tmp = get(i,j).str_imag();
					RET << _INDENT_;
					RET << name << mindex(i,j,size) << ".im " << op << " " << tmp;
					if(i == size-1 && j == size-1) RET << _LASTENDL_;
					else RET << _ENDL_;
				}
		}
		return RET.str();
	}
	
	string symplectic_compressed_assignment(const char* op, const char* name) const
	{
		ostringstream RET;
		string tmp;
		if(isreal())
		{
			cerr << "Error: A symplectic matrix set to have real values \n" ;
			exit(1);
		}
		else
		{
			for(int i = 0; i < size/2; i++)
				for(int j = 0; j < size; j++)
				{
					tmp = get(i,j).str_real();
					RET << _INDENT_;
					RET << name << mindex(i,j,size) << ".re " << op << " " << tmp;
					RET << _ENDL_;

					tmp = get(i,j).str_imag();
					RET << _INDENT_;
					RET << name << mindex(i,j,size) << ".im " << op << " " << tmp;
					if(i == size/2-1 && j == size-1) RET << _LASTENDL_;
					else RET << _ENDL_;
				}
		}
		return RET.str();
	}
	
};


class rmatrix : public pmatrix
{
public:
	rmatrix(int N, const string& name) : pmatrix(N)
	{
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
			{
				ostringstream elem;
				elem << name << mindex(i,j,size);
				rvariable tmp(elem.str().c_str());
				set(i,j, tmp);
			}
	}
	rmatrix(int N, char* name) : pmatrix(N)
	{
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
			{
				ostringstream elem;
				elem << name << mindex(i,j,size);
				rvariable tmp(elem.str().c_str());
				set(i,j, tmp);
			}
	}
	
	using pmatrix::operator=;
};


class cmatrix : public pmatrix
{
public:
	cmatrix(int N, const string& name) : pmatrix(N)
	{
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
			{
				ostringstream elem;
				elem << name << mindex(i,j,size);
				cvariable tmp(elem.str().c_str());
				set(i,j, tmp);
			}
	}
	cmatrix(int N, char* name) : pmatrix(N)
	{
		for(int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
			{
				ostringstream elem;
				elem << name << mindex(i,j,size);
				cvariable tmp(elem.str().c_str());
				set(i,j, tmp);
			}
	}
	
	using pmatrix::operator=;
};

/* A compressed symplectic matrix.
   We actually build a complete matrix here using the elements of a symplectic one. */
class spmatrix : public pmatrix
{
private: 
    void _spinit(int N,const char* name){
		for(int i = 0; i < N/2; i++)
			for(int j = 0; j < N; j++)
			{
				ostringstream elem;
				elem << name << mindex(i,j,size);
				cvariable tmp(elem.str().c_str());
				set(i,j, tmp);
			}
		for(int i = 0; i < N/2; i++)
			for(int j = 0; j < N/2; j++)
			{
				ostringstream elem;
				elem << name << mindex(i,j+N/2,size);
				cvariable tmp(elem.str().c_str());
				tmp.conjugate();
				tmp.minus();
				set(i+N/2,j, tmp);
			}
		for(int i = 0; i < N/2; i++)
			for(int j = 0; j < N/2; j++)
			{
				ostringstream elem;
				elem << name << mindex(i,j,size);
				cvariable tmp(elem.str().c_str());
				tmp.conjugate();
				set(i+N/2,j+N/2, tmp);
			}
    }
public:
	spmatrix(int N, const string& name) : pmatrix(N){
        _spinit(N,name.c_str());
	}
	spmatrix(int N, char* name) : pmatrix(N){
        _spinit(N,name);
    }
		
	using pmatrix::operator=;
};


class pvector : public sparsevector<polynomial>
{
public:
	pvector() : sparsevector<polynomial>() {}
	pvector(int N) : sparsevector<polynomial>(N) {}
	pvector(svector& vec) : sparsevector<polynomial>(vec.size)
	{
		static rmonomial index;
		for(int i = 0; i < vec.size; i++)
		{
			polynomial tmp(index, vec.get(i));
			set(i, tmp);
		}
	}
	
	using sparsevector<polynomial>::operator=;
	
	polynomial& getzero() const
	{
		static polynomial zero;
		return zero;
	}
	
	bool isreal() const
	{
		bool RET = true;
		for(KEYTYPE i = 0; i < length; i++)
			RET = RET && data[i]->value.isreal();
		return RET;
	}
	
	string assignment(const char* op, const char* name) const
	{
		ostringstream RET;
		string tmp;
		if(isreal())
		{
			for(int i = 0; i < size; i++)
			{
				tmp = get(i).str_real();
				RET << _INDENT_;
				RET << name << vindex(i,size) << " " << op << " " << tmp;
				if(i == size-1) RET << _LASTENDL_;
				else RET << _ENDL_;
			}
		}
		else
		{
			for(int i = 0; i < size; i++)
			{
				tmp = get(i).str_real();
				RET << _INDENT_;
				RET << name << vindex(i,size) << ".re " << op << " " << tmp << "; \\ \n";
				RET << _ENDL_;

				tmp = get(i).str_imag();
				RET << _INDENT_;
				RET << name << vindex(i,size) << ".im " << op << " " << tmp << "; \\ \n";
				if(i == size-1) RET << _LASTENDL_;
				else RET << _ENDL_;
			}
		}
		return RET.str();
	}
};


class rvector : public pvector
{
public:
	rvector(int N, const string& name) : pvector(N)
	{
		for(int i = 0; i < N; i++)
		{
			ostringstream elem;
			elem << name << vindex(i,size);
			rvariable tmp(elem.str().c_str());
			set(i, tmp);
		}
	}
	rvector(int N, char* name) : pvector(N)
	{
		for(int i = 0; i < N; i++)
		{
			ostringstream elem;
			elem << name << vindex(i,size);
			rvariable tmp(elem.str().c_str());
			set(i, tmp);
		}
	}

	using pvector::operator=;
};


class cvector : public pvector
{
public:
	cvector(int N, const string& name) : pvector(N)
	{
		for(int i = 0; i < N; i++)
		{
			ostringstream elem;
			elem << name << vindex(i,size);
			cvariable tmp(elem.str().c_str());
			set(i, tmp);
		}
	}
	cvector(int N, char* name) : pvector(N)
	{
		for(int i = 0; i < N; i++)
		{
			ostringstream elem;
			elem << name << vindex(i,size);
			cvariable tmp(elem.str().c_str());
			set(i, tmp);
		}
	}

	using pvector::operator=;
};
#endif
