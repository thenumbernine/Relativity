#pragma once

#include "generic_vector.h"
#include "vector.h"

template<int dim_>
struct generic_rank1 {
	enum { dim = dim_ };
	enum { rank = 1 };

	//now i need a specialization of this for when 'type' is a primitive
	// in such cases, accessing nested members results in compiler errors
	template<typename inner_type, typename scalar_type>
	struct body : public generic_vector<dim, inner_type, scalar_type, body<inner_type, scalar_type>> {
		typedef generic_vector<dim, inner_type, scalar_type, body<inner_type, scalar_type>> parent;
		body() : parent() {}
		body(const body &b) : parent(b) {}
		body(const inner_type &t) : parent(t) {}

		inner_type &operator()(const vector<1,int> &deref) { return parent::v[deref(0)]; }
		const inner_type &operator()(const vector<1,int> &deref) const { return parent::v[deref(0)]; }
	};
};

