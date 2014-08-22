#pragma once

#include "generic_dense_matrix.h"

/*
generic_symmat(i,j) == generic_symmat(j,i)
*/
template<typename type_, int dim_, typename scalar_type_, typename child>
struct generic_symmat : public generic_dense_matrix<type_, dim_, scalar_type_, child, dim_ * (dim_ + 1) / 2> {
	typedef generic_dense_matrix<type_, dim_, scalar_type_, child, dim_ * (dim_ + 1) / 2> parent;
	
	typedef typename parent::type type;
	enum { dim = parent::dim };
	typedef typename parent::scalar_type scalar_type;
	enum { size = parent::size };

	generic_symmat() : parent() {}
	generic_symmat(const child &a) : parent(a) {}

	/*
	math-index: i is the row, j is the column
	row-major: i is nested inner-most
	lower triangular: j <= i
	*/
	static int index(int i, int j) {
		if (j > i) return index(j,i);
		//i == 0: return 0
		//i == 1: return 1 + j
		//i == 2: return 1 + 2 + j
		//i == i: return i * (i+1)/2 + j
		return i * (i + 1) / 2 + j;
	}
};

