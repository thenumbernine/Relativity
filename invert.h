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
struct determinant {
	template<typename type, typename... args>
	type operator()(const ::tensor<type, args...> &a) const { 
		typedef ::tensor<type, args...> tensor;
		static_assert(tensor::rank == 2, "determinant can only handle rank-2 tensors");
		return a(0,0) * a(1,1) - a(1,0) * a(0,1);
	}
};

struct invert {
	template<typename type, typename... args>
	::tensor<type, args...> operator()(const ::tensor<type, args...> &a) const {
		typedef ::tensor<type, args...> tensor;
		static_assert(tensor::rank == 2, "inverse can only handle rank-2 tensors");
		tensor result;
		type det = determinant()(a);
		result(0,0) = a(1,1) / det;
		result(1,0) = -a(1,0) / det;
		result(0,1) = -a(0,1) / det;
		result(1,1) = a(0,0) / det;
		return result;
	}
};

