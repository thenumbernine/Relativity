#pragma once

#include "generic_array.h"

/*
adds int-based indexing and vector-vector ops to the array ops
*/
template<int size_, typename type_, typename scalar_type_, typename child>
struct generic_vector : public generic_array<size_, type_, scalar_type_, child> {
	typedef generic_array<size_, type_, scalar_type_, child> parent;

	enum { size = parent::size };
	typedef typename parent::type type;
	typedef typename parent::scalar_type scalar_type;

	//inherited constructors
	generic_vector() : parent() {}
	generic_vector(const child &a) : parent(a) {}
	generic_vector(const type &x) : parent(x) {}

	//specific dimension initializers
	generic_vector(const type &x, const type &y) {
		static_assert(size >= 2, "vector not big enough");
		parent::v[0] = x;
		parent::v[1] = y;
		for (int i = 2; i < size; ++i) parent::v[i] = type();
	}

	generic_vector(const type &x, const type &y, const type &z) {
		static_assert(size >= 3, "vector not big enough");
		parent::v[0] = x;
		parent::v[1] = y;
		parent::v[2] = z;
		for (int i = 3; i < size; ++i) parent::v[i] = type();
	}

	generic_vector(const type &x, const type &y, const type &z, const type &w) {
		static_assert(size >= 4, "vector not big enough");
		parent::v[0] = x;
		parent::v[1] = y;
		parent::v[2] = z;
		parent::v[3] = w;
		for (int i = 4; i < size; ++i) parent::v[i] = type();
	}

	//index access
	type &operator()(int i) { return parent::v[i]; }
	const type &operator()(int i) const { return parent::v[i]; }

	//product of elements / flat-space volume operator
	type volume() const {
		type vol = type(1);
		for (int i = 0; i < size; ++i) {
			vol *= parent::v[i];
		}
		return vol;
	}

	//inner product / flat-space dot product
	static type dot(const child &a, const child &b) {
		type d = type(0);
		for (int i = 0; i < size; ++i) {
			d += a.v[i] * b.v[i];
		}
		return d;
	}

	//vector/vector operations
	//i figure matrix doesn't want * and / to be matrix/matrix per-element
	//trying to copy GLSL as much as I can here

	child operator*(const child &b) const {
		const generic_vector &a = *this;
		child c;
		for (int i = 0; i < size; ++i) {
			c.v[i] = a.v[i] * b.v[i];
		}
		return c;
	}

	child operator*(const scalar_type &b) const { return parent::operator*(b); }
	
	child operator/(const child &b) const {
		const generic_vector &a = *this;
		child c;
		for (int i = 0; i < size; ++i) {
			c.v[i] = a.v[i] / b.v[i];
		}
		return c;
	}
	
	child operator/(const scalar_type &b) const { return parent::operator/(b); }
};

