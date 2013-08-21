#pragma once

#include "generic_array.h"

/*
child being the whatever curious whatever thing that returns its child
child is used for
	vector dereference index calculation
	operator return type
*/
template<int dim_, typename type_, typename child, int size_>
struct generic_dense_matrix : public generic_array<size_, type_, child> {
	typedef generic_array<size_, type_, child> parent;
	
	enum { dim = dim_ };
	enum { size = parent::size };
	typedef typename parent::type type;

	/*
	initialize to identity or zero?
	identity
	*/
	generic_dense_matrix() 
	: parent() { //<- this inits it with zero
		for (int i = 0; i < dim; ++i) {
			(*this)(i,i) = type(1);
		}
	}

	//index access
	type &operator()(int i, int j) { return parent::v[child::index(i,j)];}
	const type &operator()(int i, int j) const { return parent::v[child::index(i,j)];}
};

