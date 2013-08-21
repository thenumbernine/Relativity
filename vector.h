#pragma once

#include "generic_vector.h"

/*
vector class for fixed-size templated dimension (i.e. size) and type
*/
template<int dim_, typename type_>
struct vector : public generic_vector<dim_, type_, vector<dim_, type_> > {
	typedef generic_vector<dim_, type_, vector<dim_, type_> > parent;
	
	enum { dim = parent::size };	//size is the size of the static vector, which coincides with the dim in the vector class (not so in matrix classes)
	typedef typename parent::type type;

	//inherited constructors
	vector() : parent() {}
	vector(const vector &a) : parent(a) {}
	vector(const type &x) : parent(x) {}
	vector(const type &x, const type &y) : parent(x,y) {}
	vector(const type &x, const type &y, const type &z) : parent(x,y,z) {}

	template<typename type2>
	operator vector<dim,type2>() const {
		vector<dim,type2> result;
		for (int i = 0; i < dim; ++i) {
			result.v[i] = (type2)parent::v[i];
		}
		return result;
	}
};


