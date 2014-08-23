#pragma once

/*
inverse function
for any sort of rank-2 Tensor::Tensor object
*/

#include "Tensor/Tensor.h"

template<typename Real>
Real det22(
	Real a00, Real a01, 
	Real a10, Real a11) 
{
	return a00 * a11 - a10 * a01;
}

//I have to specialize determinant by rank
//which means (since rank is an enum rather than a template parameter)
// that I might have to specialize it per-index
// (or make use of static conditions)
template<typename InputType>
struct DeterminantClass;

template<typename Real>
struct DeterminantClass<Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<1>, Tensor::Lower<1>>>> {
	typedef Real OutputType;
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<1>, Tensor::Lower<1>>> InputType;
	OutputType operator()(const InputType &a) const { 
		return a(0,0);
	}
};

template<typename Real>
struct DeterminantClass<Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<2>, Tensor::Lower<2>>>> {
	typedef Real OutputType;
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<2>, Tensor::Lower<2>>> InputType;
	OutputType operator()(const InputType &a) const { 
		return det22(a(0,0), a(0,1), a(1,0), a(1,1));
	}
};

//this is where boost::bind would be helpful: for a generic dimension determinant function
template<typename Real>
struct DeterminantClass<Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<3>, Tensor::Lower<3>>>> {
	typedef Real OutputType;
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<3>, Tensor::Lower<3>>> InputType;
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
//so with the still uncertain decision of how to expose / parameterize lowers and uppers of both Tensor::Symmetric and nonsymmetric indexes
// (I was thinking to allow this only for lowers, but of both Tensor::Symmetric and non-Tensor::Symmetric)
//instead I'll just write the specialization.
//another perk of this is that Tensor::Symmetric needs less operations.
// though that could be incorporated into a single function if the Tensor::Tensor iterator returned both the index and the dereference, 
// and we subsequently used the Levi Civita definition in compile-time to calculate the inverse
template<typename InputType>
struct InverseClass;

template<typename Real_>
struct InverseClass<Tensor::Tensor<Real_, Tensor::Symmetric<Tensor::Lower<1>,Tensor::Lower<1>>>> {
	typedef Real_ Real;
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<1>, Tensor::Lower<1>>> InputType;
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Upper<1>, Tensor::Upper<1>>> OutputType;
	OutputType operator()(const InputType &a, const Real &det) const {
		OutputType result;
		result(0,0) = 1. / det;
		return result;
	}
};

template<typename Real_>
struct InverseClass<Tensor::Tensor<Real_, Tensor::Symmetric<Tensor::Lower<2>,Tensor::Lower<2>>>> {
	typedef Real_ Real;
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<2>, Tensor::Lower<2>>> InputType;
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Upper<2>, Tensor::Upper<2>>> OutputType;
	OutputType operator()(const InputType &a, const Real &det) const {
		OutputType result;
		result(0,0) = a(1,1) / det;
		result(1,0) = -a(1,0) / det;
		//result(0,1) = -a(0,1) / det;	//<- Tensor::Symmetric only needs one
		result(1,1) = a(0,0) / det;
		return result;
	}
};

template<typename Real_>
struct InverseClass<Tensor::Tensor<Real_, Tensor::Symmetric<Tensor::Lower<3>,Tensor::Lower<3>>>> {
	typedef Real_ Real;
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<3>, Tensor::Lower<3>>> InputType;
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Upper<3>, Tensor::Upper<3>>> OutputType;
	OutputType operator()(const InputType &a, const Real &det) const {
		OutputType result;
		//Tensor::Symmetric, so only do Tensor::Lower triangular
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
typename InverseClass<InputType>::OutputType inverse(const InputType &a, const typename InverseClass<InputType>::Real &det) {
	return InverseClass<InputType>()(a, det);
}

template<typename InputType>
typename InverseClass<InputType>::OutputType inverse(const InputType &a) {
	return InverseClass<InputType>()(a, determinant(a));
}
