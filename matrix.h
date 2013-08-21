#pragma once

#include "generic_dense_matrix.h"

template<int dim_, typename type_>
struct matrix : public generic_dense_matrix<dim_, type_, matrix<dim_, type_>, dim_ * dim_ > {
	typedef generic_dense_matrix<dim_, type_, matrix<dim_, type_>, dim_ * dim_> parent;
	
	enum { dim = parent::dim };
	typedef typename parent::type type;

	matrix() : parent() {} 

	/*
	math-index: i is the row, j is the column
	row-major: i is nested inner-most
	*/
	static int index(int i, int j) { return i + dim * j; }
	
};


