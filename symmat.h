#pragma once

#include "generic_dense_matrix.h"

template<int dim_, typename type_>
struct symmat;

template<int dim, typename type>
struct symmat_det {
	type operator()(const symmat<dim,type> &a);
};

template<int dim, typename type>
struct symmat_invert {
	symmat<dim,type> operator()(const symmat<dim,type> &a);
};

/*
symmat(i,j) == symmat(j,i)
*/
template<int dim_, typename type_>
struct symmat : public generic_dense_matrix<dim_, type_, symmat<dim_, type_>, dim_ * (dim_ + 1) / 2> {
	typedef generic_dense_matrix<dim_, type_, symmat<dim_, type_>, dim_ * (dim_ + 1) / 2> parent;
	
	enum { dim = parent::dim };
	enum { size = parent::size };
	typedef typename parent::type type;

	symmat() : parent() {}
	symmat(const symmat &a) : parent(a) {}

	/*
	math-index: i is the row, j is the column
	row-major: i is nested inner-most
	upper triangular: i <= j
	*/
	static int index(int i, int j) {
		if (i > j) return index(j,i);
		//j == 0: return 0
		//j == 1: return 1 + i
		//j == 2: return 1 + 2 + i
		//j == j: return j * (j+1)/2 + i
		return j * (j + 1) / 2 + i;
	}

	//determinant
	static type det(const symmat &a) {
		return symmat_det<dim,type>()(a);
	}
	
	//inverse
	static symmat invert(const symmat &a) {
		return symmat_invert<dim,type>()(a);
	}
};

template<typename type>
struct symmat_det<2,type> {
	type operator()(const symmat<2,type> &a) const {
		return a(0,0) * a(1,1) - a(0,1) * a(0,1);
	}
};

template<typename type>
struct symmat_invert<2,type> {
	symmat<2,type> operator()(const symmat<2,type> &a) const {
		type det = symmat<2,type>::det(a);
		symmat<2,type> inv;
		inv(0,0) = a(1,1) / det;
		inv(1,1) = a(0,0) / det;
		inv(0,1) = -a(0,1) / det;
		return inv;
	}
};


