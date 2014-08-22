#pragma once

#include <ostream>
#include "generic_vector.h"

/*
vector class for fixed-size templated dimension (i.e. size) and type
*/
template<typename type_, int dim_>
struct vector : public generic_vector<type_, dim_, type_, vector<type_, dim_> > {
	typedef generic_vector<type_, dim_, type_, vector<type_, dim_> > parent;
	
	typedef typename parent::type type;
	enum { dim = parent::size };	//size is the size of the static vector, which coincides with the dim in the vector class (not so in matrix classes)
	typedef typename parent::scalar_type scalar_type;

	//inherited constructors
	vector() : parent() {}
	vector(const vector &a) : parent(a) {}
	vector(const type &x) : parent(x) {}
	vector(const type &x, const type &y) : parent(x,y) {}
	vector(const type &x, const type &y, const type &z) : parent(x,y,z) {}
	vector(const type &x, const type &y, const type &z, const type &w) : parent(x,y,z,w) {}

	template<typename type2>
	operator vector<type2,dim>() const {
		vector<type2,dim> result;
		for (int i = 0; i < dim; ++i) {
			result.v[i] = (type2)parent::v[i];
		}
		return result;
	}
};

template<typename type, int dim>
std::ostream &operator<<(std::ostream &o, const vector<type,dim> &t) {
	const char *sep = "";
	o << "(";
	for (int i = 0; i < t.dim; ++i) {
		o << sep << t(i);
		sep = ", ";
	}
	o << ")";
	return o;
}

