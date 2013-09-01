#pragma once

/*
inverse function
for any sort of rank-2 tensor object
*/

#include "tensor.h"

template<typename real>
real det22(
	real a00, real a01, 
	real a10, real a11) 
{
	return a00 * a11 - a10 * a01;
}

//I have to specialize determinant by rank
//which means (since rank is an enum rather than a template parameter)
// that I might have to specialize it per-index
// (or make use of static conditions)
template<typename InputType>
struct DeterminantClass;

template<typename real>
struct DeterminantClass<tensor<real, symmetric<lower<1>, lower<1>>>> {
	typedef real OutputType;
	typedef tensor<real, symmetric<lower<1>, lower<1>>> InputType;
	OutputType operator()(const InputType &a) const { 
		return a(0,0);
	}
};

template<typename real>
struct DeterminantClass<tensor<real, symmetric<lower<2>, lower<2>>>> {
	typedef real OutputType;
	typedef tensor<real, symmetric<lower<2>, lower<2>>> InputType;
	OutputType operator()(const InputType &a) const { 
		return det22(a(0,0), a(0,1), a(1,0), a(1,1));
	}
};

//this is where boost::bind would be helpful: for a generic dimension determinant function
template<typename real>
struct DeterminantClass<tensor<real, symmetric<lower<3>, lower<3>>>> {
	typedef real OutputType;
	typedef tensor<real, symmetric<lower<3>, lower<3>>> InputType;
	OutputType operator()(const InputType &a) const { 
		return a(0,0) * a(1,1) * a(2,2)
			+ a(0,1) * a(1,2) * a(2,0)
			+ a(0,2) * a(1,0) * a(2,1)
			- a(2,0) * a(1,1) * a(0,2)
			- a(2,1) * a(1,2) * a(0,0)
			- a(2,2) * a(1,0) * a(0,1);
	}
};

template<typename InputType>
typename DeterminantClass<InputType>::OutputType
determinant(const InputType &a) {
	return DeterminantClass<InputType>()(a);
}

//currently only used for (gamma^ij) = (gamma_ij)^-1
//so with the still uncertain decision of how to expose / parameterize lowers and uppers of both symmetric and nonsymmetric indexes
// (I was thinking to allow this only for lowers, but of both symmetric and non-symmetric)
//instead I'll just write the specialization.
//another perk of this is that symmetric needs less operations.
// though that could be incorporated into a single function if the tensor iterator returned both the index and the dereference, 
// and we subsequently used the Levi Civita definition in compile-time to calculate the inverse
template<typename InputType>
struct InverseClass;

template<typename real_>
struct InverseClass<tensor<real_, symmetric<lower<1>,lower<1>>>> {
	typedef real_ real;
	typedef tensor<real, symmetric<lower<1>, lower<1>>> InputType;
	typedef tensor<real, symmetric<upper<1>, upper<1>>> OutputType;
	OutputType operator()(const InputType &a, const real &det) const {
		OutputType result;
		result(0,0) = 1. / det;
		return result;
	}
};

template<typename real_>
struct InverseClass<tensor<real_, symmetric<lower<2>,lower<2>>>> {
	typedef real_ real;
	typedef tensor<real, symmetric<lower<2>, lower<2>>> InputType;
	typedef tensor<real, symmetric<upper<2>, upper<2>>> OutputType;
	OutputType operator()(const InputType &a, const real &det) const {
		OutputType result;
		result(0,0) = a(1,1) / det;
		result(1,0) = -a(1,0) / det;
		//result(0,1) = -a(0,1) / det;	//<- symmetric only needs one
		result(1,1) = a(0,0) / det;
		return result;
	}
};

template<typename real_>
struct InverseClass<tensor<real_, symmetric<lower<3>,lower<3>>>> {
	typedef real_ real;
	typedef tensor<real, symmetric<lower<3>, lower<3>>> InputType;
	typedef tensor<real, symmetric<upper<3>, upper<3>>> OutputType;
	OutputType operator()(const InputType &a, const real &det) const {
		OutputType result;
		//symmetric, so only do lower triangular
		result(0,0) = det22(a(1,1), a(1,2), a(2,1), a(2,2)) / det;
		result(1,0) = det22(a(1,2), a(1,0), a(2,2), a(2,0)) / det;
		result(1,1) = det22(a(0,0), a(0,2), a(2,0), a(2,2)) / det;
		result(2,0) = det22(a(1,0), a(1,1), a(2,0), a(2,1)) / det;
		result(2,1) = det22(a(0,1), a(0,0), a(2,1), a(2,0)) / det;
		result(2,2) = det22(a(0,0), a(0,1), a(1,0), a(1,1)) / det;
		return result;
	}
};



template<typename InputType>
typename InverseClass<InputType>::OutputType inverse(const InputType &a, const typename InverseClass<InputType>::real &det) {
	return InverseClass<InputType>()(a, det);
}

