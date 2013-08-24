#pragma once

#include "generic_dense_matrix.h"

template<typename type, typename owner_type>
struct antisymmat_accessor {
	owner_type *owner;
	int offset;
	bool flip;

	antisymmat_accessor(owner_type *owner_, int offset_, bool flip_) 
	: owner(owner_), offset(offset_), flip(flip_) {}

	operator type() const { 
		if (!owner) return type();
		if (flip) {
			return -owner->v[offset];
		} else {
			return owner->v[offset];
		}
	}
	
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
};

/*
generic_antisymmat(i,j) == -generic_antisymmat(j,i)
therefore generic_antisymmat(i,i) == 0
*/
template<int dim_, typename type_, typename scalar_type_, typename child>
struct generic_antisymmat : public generic_dense_matrix<dim_, type_, scalar_type_, child, dim_ * (dim_ - 1) / 2> {
	typedef generic_dense_matrix<dim_, type_, scalar_type_, child, dim_ * (dim_ - 1) / 2> parent;
	
	enum { dim = parent::dim };
	enum { size = parent::size };
	typedef typename parent::type type;
	typedef typename parent::scalar_type scalar_type;
	typedef antisymmat_accessor<type_, child> accessor;

	generic_antisymmat() : parent() {}
	generic_antisymmat(const child &a) : parent(a) {}
	
	//index access
	accessor operator()(int i, int j) { 
		if (i == j) return accessor(NULL, 0, false);
		if (i < j) return accessor(this, index(i,j), false);
		if (i > j) return accessor(this, index(j,i), true);
	}
	const accessor operator()(int i, int j) const {
		if (i == j) return accessor(NULL, 0, false);
		if (i < j) return accessor(this, index(i,j), false);
		if (i > j) return accessor(this, index(j,i), true);
	}

	/*
	math-index: i is the row, j is the column
	row-major: i is nested inner-most
	upper triangular: i <= j
	*/
	static int index(int i, int j) {
		if (i > j) return index(j,i);
		//j == 0: return 0
		//j == 1: return 1 + i
		//j == 2: return 1 + 2 + i
		//j == j: return j * (j+1)/2 + i
		return j * (j + 1) / 2 + i;
	}
};


