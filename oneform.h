#pragma once

#include "generic_vector.h"

/*
one-form class
*/
template<int dim_, typename type_>
struct oneform : public generic_vector<dim_, type_, type_, oneform<dim_, type_> > {
	typedef generic_vector<dim_, type_, type_, oneform<dim_, type_> > parent;

	enum { dim = parent::size };	//size is the size of the static vector, which coincides with the dim in the oneform class (not so in matrix classes)
	typedef typename parent::type type;
	typedef typename parent::scalar_type scalar_type;

	//inherited constructors
	oneform() : parent() {}
	oneform(const oneform &a) : parent(a) {}
	oneform(const type &x) : parent(x) {}
	oneform(const type &x, const type &y) : parent(x,y) {}
	oneform(const type &x, const type &y, const type &z) : parent(x,y,z) {}
	oneform(const type &x, const type &y, const type &z, const type &w) : parent(x,y,z,w) {}

	//cast operation
	//needs to know the child's template parameters
	//I can't expose the child template base fully because it has no agreed upon template pattern
	//I could pass in a child-typed generation struct ...
	//which means I could expose a means to generate 
	//	???
	template<typename type2>
	operator oneform<dim,type2>() const {
		oneform<dim,type2> result;
		for (int i = 0; i < dim; ++i) {
			result.v[i] = (type2)parent::v[i];
		}
		return result;
	}
};

