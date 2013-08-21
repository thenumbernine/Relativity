#pragma once

#include "generic_array.h"

/*
vector class for fixed-size templated dimension (i.e. size) and type
*/
template<int dim_, typename type_>
struct vec : public generic_array<dim_, type_, vec<dim_, type_> > {
	typedef generic_array<dim_, type_, vec<dim_, type_> > parent;
	
	enum { dim = parent::size };	//size is the size of the static vector, which coincides with the dim in the vec class (not in other classes)
	typedef typename parent::type type;

	vec() : parent() {}
	vec(const vec &a) : parent(a) {}
	vec(const type &x) : parent(x) {}

	//specific dimension initializers
	vec(const type &x, const type &y) {
		assert(dim >= 2);
		parent::v[0] = x;
		parent::v[1] = y;
		for (int i = 2; i < dim; ++i) parent::v[i] = type();
	}

	vec(const type &x, const type &y, const type &z) {
		assert(dim >= 3);
		parent::v[0] = x;
		parent::v[1] = y;
		parent::v[2] = z;
		for (int i = 3; i < dim; ++i) parent::v[i] = type();
	}

	//index access
	//I am very tempted to move these back to generic_array somehow for the sake of all implemented operators
	// but matrix would inherit and subsequently need these hidden
	type &operator()(int i) { return parent::v[i]; }
	const type &operator()(int i) const { return parent::v[i]; }

	//cast operation
	template<typename type2>
	operator vec<dim,type2>() const {
		vec<dim,type2> result;
		for (int i = 0; i < dim; ++i) {
			result(i) = (type2)(*this)(i);
		}
		return result;
	}

	//product of elements / flat-space volume operator
	type volume() const {
		type vol = type(1);
		for (int i = 0; i < dim; ++i) {
			vol *= parent::v[i];
		}
		return vol;
	}

	//inner product / flat-space dot product
	static type dot(const vec &a, const vec &b) {
		type d = type(0);
		for (int i = 0; i < dim; ++i) {
			d += a(i) * b(i);
		}
		return d;
	}

	//vec/vec operations
	//i figure matrix doesn't want * and / to be matrix/matrix per-element
	//trying to copy GLSL as much as I can here

	vec operator*(const vec &b) const {
		const vec &a = *this;
		vec c;
		for (int i = 0; i < dim; ++i) {
			c.v[i] = a.v[i] * b.v[i];
		}
		return c;
	}
	
	vec operator/(const vec &b) const {
		const vec &a = *this;
		vec c;
		for (int i = 0; i < dim; ++i) {
			c.v[i] = a.v[i] / b.v[i];
		}
		return c;
	}
};


