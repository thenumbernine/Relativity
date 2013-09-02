#pragma once

#include "generic_dense_matrix.h"

template<typename type, typename owner_type>
struct antisymmat_accessor {
	owner_type *owner;
	int offset;
	bool flip;

	antisymmat_accessor(owner_type *owner_, int offset_, bool flip_) 
	: owner(owner_), offset(offset_), flip(flip_) {}

	antisymmat_accessor &operator=(const type &value) {
		if (owner) {
			if (flip) {
				owner->v[offset] = -value;
			} else {
				owner->v[offset] = value;
			}
		}
		return *this;
	}
	
	//instead of returning this
	// we could return a tensor whose body is this ...
	//that would solve the dereference (and any other abstraction) issues ...
	operator type&() const { 
		if (!owner) return type();
		if (flip) {
			return -owner->v[offset];
		} else {
			return owner->v[offset];
		}
	}
};

template<typename type, typename owner_type>
struct antisymmat_accessor_const {
	const owner_type *owner;
	int offset;
	bool flip;

	antisymmat_accessor_const(const owner_type *owner_, int offset_, bool flip_) 
	: owner(owner_), offset(offset_), flip(flip_) {}

	//instead of returning this
	// we could return a tensor whose body is this ...
	//that would solve the dereference (and any other abstraction) issues ...
	operator type() const { 
		if (!owner) return type();
		if (flip) {
			return -owner->v[offset];
		} else {
			return owner->v[offset];
		}
	}
};

/*
generic_antisymmat(i,j) == -generic_antisymmat(j,i)
therefore generic_antisymmat(i,i) == 0
*/
template<typename type_, int dim_, typename scalar_type_, typename child>
struct generic_antisymmat : public generic_dense_matrix<type_, dim_, scalar_type_, child, dim_ * (dim_ - 1) / 2> {
	typedef generic_dense_matrix<type_, dim_, scalar_type_, child, dim_ * (dim_ - 1) / 2> parent;
	
	typedef typename parent::type type;
	enum { dim = parent::dim };
	typedef typename parent::scalar_type scalar_type;
	enum { size = parent::size };
	typedef antisymmat_accessor<type_, child> accessor;
	typedef antisymmat_accessor_const<type_, child> accessor_const;

	generic_antisymmat() : parent() {}
	generic_antisymmat(const child &a) : parent(a) {}
	
	//index access
	accessor operator()(int i, int j) {
		if (i == j) return accessor(NULL, 0, false);
		if (i < j) return accessor(this, index(i,j), false);
		if (i > j) return accessor(this, index(j,i), true);
	}
	accessor_const operator()(int i, int j) const {
		if (i == j) return accessor_const(NULL, 0, false);
		if (i < j) return accessor_const(this, index(i,j), false);
		if (i > j) return accessor_const(this, index(j,i), true);
	}

	/*
	math-index: i is the row, j is the column
	row-major: i is nested inner-most
	upper triangular: j <= i
	*/
	static int index(int i, int j) {
		if (j > i) return index(j,i);
		//i == 0: return 0
		//i == 1: return 1 + j
		//i == 2: return 1 + 2 + j
		//i == i: return i * (i+1)/2 + j
		return i * (i + 1) / 2 + j;
	}
};


