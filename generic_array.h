#pragma once

#include "crtp_cast.h"
#include "clamp.h"

/*
the 'parent' curious pattern whatever for generic_vector and generic_matrix
	type_ is the type of each element
	size_ is the size of the dense vector
	scalar_type_ is the type of the innermost nested scalars (in the event that elements are generic_array's themselves)
	child is the curious reoccurring chlid class that uses this
*/
template<typename type_, int size_, typename scalar_type_, typename child>
struct generic_array {
	typedef type_ type;
	enum { size = size_ };
	typedef scalar_type_ scalar_type;

	type v[size];

	//default ctor: init all to zero
	generic_array() {
		for (int i = 0; i < size; ++i) {
			v[i] = type();
		}
	}

	//default copy ctor
	generic_array(const child &a) {
		for (int i = 0; i < size; ++i) {
			v[i] = a.v[i];
		}
	}

	//default fill ctor
	generic_array(const type &x) {
		for (int i = 0; i < size; ++i) {
			v[i] = x;
		}
	}

	//bounds
	static child clamp(const child &a, const child &min, const child &max) {
		child b;
		for (int i = 0; i < size; ++i) {
			b.v[i] = clamp(a.v[i], min.v[i], max.v[i]);
		}
		return b;
	}

	//boolean operations

	bool operator==(const child &b) const {
		const generic_array &a = *this;
		for (int i = 0; i < size; ++i) {
			if (a.v[i] != b.v[i]) return false;
		}
		return true;
	}

	bool operator!=(const child &b) const { return ! this->operator==(b); }

	//unary math operator

	child operator-() const {
		const generic_array &a = *this;
		child b;
		for (int i = 0; i < size; ++i) {
			b.v[i] = -a.v[i];
		}
		return b;
	}

	//pairwise math operators

	child operator+(const child &b) const {
		const generic_array &a = *this;
		child c;
		for (int i = 0; i < size; ++i) {
			c.v[i] = a.v[i] + b.v[i];
		}
		return c;
	}
	
	child operator-(const child &b) const {
		const generic_array &a = *this;
		child c;
		for (int i = 0; i < size; ++i) {
			c.v[i] = a.v[i] - b.v[i];
		}
		return c;
	}

	//I'm omitting * and / since matrix defines them separately
	//what if I overloaded them?

	//scalar math operators

	//I'm omitting scalar + and -, but not for any good reason

	child operator*(const scalar_type &b) const {
		const generic_array &a = *this;
		child c;
		for (int i = 0; i < size; ++i) {
			c.v[i] = a.v[i] * b;
		}
		return c;
	}


	child operator/(const scalar_type &b) const {
		const generic_array &a = *this;
		child c;
		for (int i = 0; i < size; ++i) {
			c.v[i] = a.v[i] / b;
		}
		return c;
	}

	//these will have to be overridden to return the correct type (unless I cast it here...)
	
	child &operator+=(const child &b) {
		for (int i = 0; i < size; ++i) {
			v[i] += b.v[i];
		}
		return crtp_cast<child>(*this);
	}

	child &operator-=(const child &b) {
		for (int i = 0; i < size; ++i) {
			v[i] -= b.v[i];
		}
		return crtp_cast<child>(*this);
	}

	child &operator*=(const scalar_type &b) {
		for (int i = 0; i < size; ++i) {
			v[i] *= b;
		}
		return crtp_cast<child>(*this);
	}

	child &operator/=(const scalar_type &b) {
		for (int i = 0; i < size; ++i) {
			v[i] /= b;
		}
		return crtp_cast<child>(*this);
	}
};

