#pragma once

#include "generic_array.h"

/*
child being the whatever curious whatever thing that returns its child
child is used for
	vector dereference index calculation
	operator return type
*/
template<int dim_, typename type_, typename scalar_type_, typename child, int size_>
struct generic_dense_matrix : public generic_array<size_, type_, scalar_type_, child> {
	typedef generic_array<size_, type_, scalar_type_, child> parent;
	
	enum { dim = dim_ };
	enum { size = parent::size };
	typedef typename parent::type type;
	typedef typename parent::scalar_type scalar_type;

	/*
	initialize to identity or zero?
	zero
	- this coincides with other vector ctors
	- identity would be good for scalar types, but not so much for matrix types
	*/
	generic_dense_matrix() : parent() {}

	//index access
	type &operator()(int i, int j) { return parent::v[child::index(i,j)]; }
	const type &operator()(int i, int j) const { return parent::v[child::index(i,j)]; }
	type &operator()(const ::vector<2,int> &deref) { return parent::v[child::index(deref(0),deref(1))]; }
	const type &operator()(const ::vector<2,int> &deref) const { return parent::v[child::index(deref(0),deref(1))]; }
};

