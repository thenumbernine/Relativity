#pragma once

/*
invert function
for any sort of rank-2 tensor object
*/

#include "tensor.h"

//I have to specialize determinant by rank
//which means (since rank is an enum rather than a template parameter)
// that I might have to specialize it per-index
// (or make use of static conditions)
template<int dim, typename real, typename tensor_type>
struct determinantClass;

template<int dim, typename real>
struct determinantClass<dim, real, tensor<real, symmetric<lower<dim>, lower<dim>>>> {
	typedef tensor<real, symmetric<lower<dim>, lower<dim>>> tensor_type;
	real operator()(const tensor_type &a) const { 
		return a(0,0) * a(1,1) - a(1,0) * a(0,1);
	}
};

template<int dim, typename real>
real determinant(const tensor<real, symmetric<lower<dim>, lower<dim>>> &a) {
	return determinantClass<dim, real, tensor<real, symmetric<lower<dim>, lower<dim>>>>(a);
}

//currently only used for (gamma^ij) = (gamma_ij)^-1
//so with the still uncertain decision of how to expose / parameterize lowers and uppers of both symmetric and nonsymmetric indexes
// (I was thinking to allow this only for lowers, but of both symmetric and non-symmetric)
//instead I'll just write the specialization.
//another perk of this is that symmetric needs less operations.
// though that could be incorporated into a single function if the tensor iterator returned both the index and the dereference, 
// and we subsequently used the Levi Civita definition in compile-time to calculate the inverse
template<typename input_type>
struct invertClass;

template<int dim, typename real>
struct invertClass<tensor<real, symmetric<lower<dim>,lower<dim>>>> {
	typedef tensor<real, symmetric<lower<dim>, lower<dim>>> input_type;
	typedef tensor<real, symmetric<upper<dim>, upper<dim>>> output_type;
	output_type operator()(const input_type &a) const {
		output_type result;
		real det = determinant(a);
		result(0,0) = a(1,1) / det;
		result(1,0) = -a(1,0) / det;
		//result(0,1) = -a(0,1) / det;	//<- symmetric only needs one
		result(1,1) = a(0,0) / det;
		return result;
	}
};

template<typename input_type>
typename invertClass<input_type>::output_type invert(const input_type &a) {
	return invertClass<input_type>()(a);
}

