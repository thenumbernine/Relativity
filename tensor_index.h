#pragma once

/*
holds different indexes for tensors
these are what go in the args... of the tensor class
*/

#include "vector.h"
#include "generic_rank1.h"			//upper, lower
#include "generic_symmat.h"			//symmetric
#include "generic_antisymmat.h"		//antisymmetric

/*
dim = dimension of this upper index
used like so: 
	typedef tensor<real, upper<3>> vector3;
*/
template<int dim>
struct upper : public generic_rank1<dim> {};

/*
dim = dimension of this lower index
used like so:
	typedef tensor<real, lower<3>> oneform3;
*/
template<int dim>
struct lower : public generic_rank1<dim> {};

/*
used like so: 
	typedef tensor<real, symmetric<lower<3>, lower<3>>> metric3;
	typedef tensor<real, symmetric<lower<3>, lower<3>>, symmetric<lower<3>, lower<3>>> riemann_metric3;
*/
template<typename index1, typename index2>
struct symmetric {
	static_assert(index1::dim == index2::dim, "symmetric can only accept two indexes of equal dimension");
	enum { dim = index1::dim };
	enum { rank = index1::rank + index2::rank };
	
	template<typename inner_type, typename scalar_type>
	struct body : public generic_symmat<inner_type, index1::dim, scalar_type, body<inner_type, scalar_type>> {
		typedef generic_symmat<inner_type, index1::dim, scalar_type, body<inner_type, scalar_type>> parent;
		body() : parent() {}
		body(const body &b) : parent(b) {}
		body(const inner_type &t) : parent(t) {}
	};
};

template<typename index1, typename index2>
struct antisymmetric {
	static_assert(index1::dim == index2::dim, "symmetric can only accept two indexes of equal dimension");
	enum { dim = index1::dim };
	enum { rank = index1::rank + index2::rank };

	template<typename inner_type, typename scalar_type>
	struct body : public generic_antisymmat<inner_type, index1::dim, scalar_type, body<inner_type, scalar_type>> {
		typedef generic_antisymmat<inner_type, index1::dim, scalar_type, body<inner_type, scalar_type>> parent;
		body() : parent() {}
		body(const body &b) : parent(b) {}
		body(const inner_type &t) : parent(t) {}
	};
};

