#pragma once

#include <ostream>
#include "generic_vector.h"

/*
vector class for fixed-size templated dimension (i.e. size) and type
*/
template<int dim_, typename type_>
struct vector : public generic_vector<dim_, type_, type_, vector<dim_, type_> > {
	typedef generic_vector<dim_, type_, type_, vector<dim_, type_> > parent;
	
	enum { dim = parent::size };	//size is the size of the static vector, which coincides with the dim in the vector class (not so in matrix classes)
	typedef typename parent::type type;
	typedef typename parent::scalar_type scalar_type;

	//inherited constructors
	vector() : parent() {}
	vector(const vector &a) : parent(a) {}
	vector(const type &x) : parent(x) {}
	vector(const type &x, const type &y) : parent(x,y) {}
	vector(const type &x, const type &y, const type &z) : parent(x,y,z) {}
	vector(const type &x, const type &y, const type &z, const type &w) : parent(x,y,z,w) {}

	template<typename type2>
	operator vector<dim,type2>() const {
		vector<dim,type2> result;
		for (int i = 0; i < dim; ++i) {
			result.v[i] = (type2)parent::v[i];
		}
		return result;
	}
};

template<int dim, typename type>
std::ostream &operator<<(std::ostream &o, const vector<dim,type> &t) {
	const char *sep = "";
	o << "(";
	for (int i = 0; i < t.dim; ++i) {
		o << sep << t(i);
		sep = ", ";
	}
	o << ")";
	return o;
}

