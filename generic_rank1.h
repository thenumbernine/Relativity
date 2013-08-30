#pragma once

#include "generic_vector.h"
#include "vector.h"

template<int dim_>
struct generic_rank1 {
	enum { dim = dim_ };
	enum { rank = 1 };

	//now i need a specialization of this for when 'type' is a primitive
	// in such cases, accessing nested members results in compiler errors
	template<typename InnerType, typename scalar_type>
	struct body : public generic_vector<InnerType, dim, scalar_type, body<InnerType, scalar_type>> {
		typedef generic_vector<InnerType, dim, scalar_type, body<InnerType, scalar_type>> parent;
		body() : parent() {}
		body(const body &b) : parent(b) {}
		body(const InnerType &t) : parent(t) {}

		InnerType &operator()(const vector<int,1> &deref) { return parent::v[deref(0)]; }
		const InnerType &operator()(const vector<int,1> &deref) const { return parent::v[deref(0)]; }
	};
};

