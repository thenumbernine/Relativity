#pragma once

#include "generic_dense_matrix.h"

/*
symmat(i,j) == symmat(j,i)
*/
template<int dim_, typename type_, typename scalar_type_, typename child>
struct generic_symmat : public generic_dense_matrix<dim_, type_, scalar_type_, child, dim_ * (dim_ + 1) / 2> {
	typedef generic_dense_matrix<dim_, type_, scalar_type_, child, dim_ * (dim_ + 1) / 2> parent;
	
	enum { dim = parent::dim };
	enum { size = parent::size };
	typedef typename parent::type type;
	typedef typename parent::scalar_type scalar_type;

	generic_symmat() : parent() {}
	generic_symmat(const child &a) : parent(a) {}

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
};

template<int dim_, typename type_>
struct symmat : public generic_symmat<dim_, type_, type_, symmat<dim_, type_>> {
	typedef generic_symmat<dim_, type_, type_, symmat<dim_, type_>> parent;

	enum { dim = parent::dim };
	enum { size = parent::size };
	typedef typename parent::type type;
	typedef typename parent::scalar_type scalar_type;

	symmat() : parent() {}
	symmat(const symmat &a) : parent(a) {}
};

