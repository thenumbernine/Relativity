#pragma once

#include "cell.h"
#include "Tensor/clamp.h"
#include "Tensor/Vector.h"
#include "Tensor/Tensor.h"
#include "Tensor/Grid.h"
#include "Common/Macros.h"	//numberof

//handles offset up to -size of the grid
template<typename Type, int rank> 
const Type& getGridOffset(const Tensor::Grid<Type, rank>& grid, Tensor::intN<rank> index, int dim, int offset) {
	index(dim) += offset;
	index(dim) += grid.size(dim);
	index(dim) %= grid.size(dim);
	return grid(index);
}

/*
partial derivative index operator
(partial derivative alone one coordinate)

finite difference coefficients for center-space finite-difference partial derivatives found at
http://en.wikipedia.org/wiki/Finite_difference_coefficients
*/

template<typename Real, int accuracy>
struct PartialDerivCoeffs;

template<typename Real>
struct PartialDerivCoeffs<Real, 2> {
	static const Real coeffs[1];
};
template<typename Real>
const Real PartialDerivCoeffs<Real, 2>::coeffs[1] = { 1./2. };

template<typename Real>
struct PartialDerivCoeffs<Real, 4> {
	static const Real coeffs[2];
};
template<typename Real>
const Real PartialDerivCoeffs<Real, 4>::coeffs[2] = { 2./3., -1./12. };

template<typename Real>
struct PartialDerivCoeffs<Real, 6> {
	static const Real coeffs[3];
};
template<typename Real>
const Real PartialDerivCoeffs<Real, 6>::coeffs[3] = { 3./4., -3./20., 1./60. };

template<typename Real>
struct PartialDerivCoeffs<Real, 8> {
	static const Real coeffs[4];
};
template<typename Real>
const Real PartialDerivCoeffs<Real, 8>::coeffs[4] = { 4./5., -1./5., 4./105., -1./280. };

template<typename Real, int dim, typename CellType, typename InputType, int accuracy>
struct PartialDerivativeCoordinateClass {
	using Coeffs = PartialDerivCoeffs<Real, accuracy>;
	InputType operator()(const Tensor::Grid<CellType, dim> &grid, const InputType CellType::*field, const Tensor::_vec<Real, dim> &dx, const Tensor::intN<dim> &index, int k) {
		InputType result = InputType();
		Real oneOverDx = 1. / dx(k);
		for (int i = 1; i <= (int)numberof(Coeffs::coeffs); ++i) {
			result += (getGridOffset(grid, index,k,i).*field - getGridOffset(grid, index,k,-i).*field) * (Coeffs::coeffs[i-1] * oneOverDx);
		}
		return result;
	}
};
template<typename Real, int dim, typename CellType, typename InputType, int accuracy>
InputType partialDerivativeCoordinate(const Tensor::Grid<CellType, dim> &grid, const InputType CellType::*field, const Tensor::_vec<Real, dim> &dx, const Tensor::intN<dim> &index, int k) {
	return PartialDerivativeCoordinateClass<Real, dim, CellType, InputType, accuracy>()(grid, field, dx, index, k);
}


template<typename Real, int accuracy>
struct PartialSecondDerivCoeffs;

template<typename Real>
struct PartialSecondDerivCoeffs<Real, 2> {
	static const Real coeffs[2];
};
template<typename Real>
const Real PartialSecondDerivCoeffs<Real, 2>::coeffs[2] = { -2., 1. };

template<typename Real>
struct PartialSecondDerivCoeffs<Real, 4> {
	static const Real coeffs[3];
};
template<typename Real>
const Real PartialSecondDerivCoeffs<Real, 4>::coeffs[3] = { -5./2., 4./3., -1./12. };

template<typename Real>
struct PartialSecondDerivCoeffs<Real, 6> {
	static const Real coeffs[4];
};
template<typename Real>
const Real PartialSecondDerivCoeffs<Real, 6>::coeffs[4] = { -49./18., 3./2., -3./20., 1./90. };

template<typename Real>
struct PartialSecondDerivCoeffs<Real, 8> {
	static const Real coeffs[5];
};
template<typename Real>
const Real PartialSecondDerivCoeffs<Real, 8>::coeffs[5] = { -205./72., 8./5., -1./5., 8./315., -1./560. };


template<typename Real, int dim, typename CellType, typename InputType, int accuracy>
struct PartialSecondDerivativeCoordinateClass {
	using Coeffs = PartialSecondDerivCoeffs<Real, accuracy>;
	InputType operator()(const Tensor::Grid<CellType, dim> &grid, const InputType CellType::*field, const Tensor::_vec<Real, dim> &dx, const Tensor::intN<dim> &index, int k) {
		Real oneOverDxSquared = 1. / (dx(k) * dx(k));
		InputType result = grid(index).*field * (Coeffs::coeffs[0] * oneOverDxSquared);
		for (int i = 1; i < (int)numberof(Coeffs::coeffs); ++i) {
			result += (getGridOffset(grid, index,k,i).*field + getGridOffset(grid, index,k,-i).*field) * (Coeffs::coeffs[i] * oneOverDxSquared);  
		}
		return result;
	}
};
template<typename Real, int dim, typename CellType, typename InputType, int accuracy>
InputType partialSecondDerivativeCoordinate(const Tensor::Grid<CellType, dim> &grid, const InputType CellType::*field, const Tensor::_vec<Real, dim> &dx, const Tensor::intN<dim> &index, int k) {
	return PartialSecondDerivativeCoordinateClass<Real, dim, CellType, InputType, accuracy>()(grid, field, dx, index, k);
}


/*
partial second derivatives of vectors
until I solve a few things with dereferencing, these will have to be specialized to
*/

template<typename Real, int dim, typename CellType, typename PartialCellType, typename InputType>
struct PartialSecondDerivativeClass;

template<typename Real, int dim, typename CellType, typename PartialCellType>
struct PartialSecondDerivativeClass<Real, dim, CellType, PartialCellType, Real> {
	using InputType = Real;
	using PartialInputType = Tensor::_vec<Real, dim>;
	using OutputType = Tensor::_sym<Real, dim>;
	using Grid = Tensor::Grid<CellType, dim>;
	using PartialGrid = Tensor::Grid<PartialCellType, dim>;
	using CellFieldType = const InputType CellType::*;
	using CellPartialFieldType = const PartialInputType PartialCellType::*;
	OutputType operator()(const Grid &grid, CellFieldType field, const PartialGrid &partialGrid, CellPartialFieldType partialField, const Tensor::_vec<Real, dim> &dx, const Tensor::intN<dim> &index) {
		OutputType result;
				
		//partial2_field_ll(i,j) = partial_i partial_j field
		for (int i = 0; i < dim; ++i) {
			PartialInputType partial_field_l_wrt_xi = partialDerivativeCoordinate<Real, dim, PartialCellType, PartialInputType, 8>(partialGrid, partialField, dx, index, i);
			for (int j = 0; j <= i; ++j) {
				if (i == j) {
					result(i,j) = partialSecondDerivativeCoordinate<Real, dim, CellType, InputType, 8>(grid, field, dx, index, i);
				} else {
					result(i,j) = partial_field_l_wrt_xi(j);
				}
			}
		}

		return result;
	}
};

template<typename Real, int dim, typename CellType, typename PartialCellType>
struct PartialSecondDerivativeClass<Real, dim, CellType, PartialCellType, Tensor::_vec<Real, dim>> {
	using InputType = Tensor::_vec<Real, dim>;
	using PartialInputType = Tensor::_mat<Real, dim, dim>;
	using OutputType = Tensor::_sym<Tensor::_vec<Real, dim>, dim>;
	using Grid = Tensor::Grid<CellType, dim>;
	using PartialGrid = Tensor::Grid<PartialCellType, dim>;
	using CellFieldType = const InputType CellType::*;
	using CellPartialFieldType = const PartialInputType PartialCellType::*;
	OutputType operator()(const Grid &grid, CellFieldType field, const PartialGrid &partialGrid, CellPartialFieldType partialField, const Tensor::_vec<Real, dim> &dx, const Tensor::intN<dim> &index) {
		OutputType result;
			
		//partial2_field_lll(k,l,i) = partial_k partial_l field_i
		for (int k = 0; k < dim; ++k) {
			PartialInputType partial_field_ll_wrt_xk = partialDerivativeCoordinate<Real, dim, PartialCellType, PartialInputType, 8>(partialGrid, partialField, dx, index, k);
			for (int l = 0; l <= k; ++l) {
				if (k == l) {
					InputType partial2_field_l_wrt_kk = partialSecondDerivativeCoordinate<Real, dim, CellType, InputType, 8>(grid, field, dx, index, k);
					for (int i = 0; i < dim; ++i) {
						result(k,l,i) = partial2_field_l_wrt_kk(i);
					}
				} else {
					for (int i = 0; i < dim; ++i) {
						result(k,l,i) = partial_field_ll_wrt_xk(l,i);
					}
				}
			}
		}

		return result;
	}
};
template<typename Real, int dim, typename CellType, typename PartialCellType>
struct PartialSecondDerivativeClass<Real, dim, CellType, PartialCellType, Tensor::_sym<Real, dim>> {
	using InputType = Tensor::_sym<Real, dim>;
	using PartialInputType = Tensor::_vec<Tensor::_sym<Real, dim>, dim>;
	using OutputType = Tensor::_sym<Tensor::_sym<Real, dim>, dim>;
	using Grid = Tensor::Grid<CellType, dim>;
	using PartialGrid = Tensor::Grid<PartialCellType, dim>;
	using CellFieldType = const InputType CellType::*;
	using CellPartialFieldType = const PartialInputType PartialCellType::*;
	OutputType operator()(const Grid &grid, CellFieldType field, const PartialGrid &partialGrid, CellPartialFieldType partialField, const Tensor::_vec<Real, dim> &dx, const Tensor::intN<dim> &index) {
		OutputType result;
			
		//partial2_field_llll(k,l,i,j) = partial_k partial_l field_ij
		for (int k = 0; k < dim; ++k) {
			PartialInputType partial_field_lll_wrt_xk = partialDerivativeCoordinate<Real, dim, PartialCellType, PartialInputType, 8>(partialGrid, partialField, dx, index, k);
			for (int l = 0; l <= k; ++l) {
				if (k == l) {
					InputType partial2_field_ll_wrt_kk = partialSecondDerivativeCoordinate<Real, dim, CellType, InputType, 8>(grid, field, dx, index, k);
					for (int i = 0; i < dim; ++i) {
						for (int j = 0; j <= i; ++j) {
							result(k,l,i,j) = partial2_field_ll_wrt_kk(i,j);
						}
					}
				} else {
					for (int i = 0; i < dim; ++i) {
						for (int j = 0; j <= i; ++j) {
							result(k,l,i,j) = partial_field_lll_wrt_xk(l,i,j);
						}
					}
				}
			}
		}

		return result;
	}
};

template<typename Real, int dim, typename CellType, typename PartialCellType, typename InputType>
typename PartialSecondDerivativeClass<Real, dim, CellType, PartialCellType, InputType>::OutputType
partialSecondDerivative(
	const Tensor::Grid<CellType, dim> &grid,
	const InputType CellType::*field,
	const Tensor::Grid<PartialCellType, dim> &partialGrid,
	typename PartialSecondDerivativeClass<Real, dim, CellType, PartialCellType, InputType>::CellPartialFieldType partialField,
	const Tensor::_vec<Real, dim> &dx,
	const Tensor::intN<dim> &index)
{
	return PartialSecondDerivativeClass<Real, dim, CellType, PartialCellType, InputType>()(grid, field, partialGrid, partialField, dx, index);
}


/*
covariant derivative
depends on conn_ull, partialDerivative, and the upper/lower indexes of the Tensor::Tensor

looks like symmetric lower and symmetric upper produce like results,
but symmetric mixed does not

I think I'll make specializations for that reason ...
*/

//scalar specialization
//diff_i alpha = partial_i alpha
template<typename Real, int dim, typename CellType, typename ConnCellType>
struct CovariantDerivativeClass_0 {
	using InputType = Real;
	using OutputType = Tensor::_vec<Real, dim>;
	using Grid = Tensor::Grid<CellType, dim>;
	using ConnGrid = Tensor::Grid<ConnCellType, dim>;
	using CellFieldType = const InputType CellType::*;
	using TensorUSL = typename ConnCellType::TensorUSL;
	using ConnFieldType = const TensorUSL ConnCellType::*;
	static OutputType exec(const Grid &grid, CellFieldType field, const ConnGrid &connGrid, ConnFieldType connField_ull, const Tensor::_vec<Real, dim> &dx, const Tensor::intN<dim> &index) {
		static constexpr int accuracy = 8;	//"order", not "accuracy", and put this all in one place plz
		return Tensor::partialDerivativeGrid<accuracy, Real, dim, InputType>(
			index, dx, [&](decltype(index) index) -> InputType {
				return grid(index).*field;
			});
	}
};

//Tensor::Vector specialization
//diff_i v^j = partial_i v^j + conn^j_ki v^k
template<typename Real, int dim, typename CellType, typename ConnCellType>
struct CovariantDerivativeClass_U {
	using InputType = Tensor::_vec<Real, dim>;
	using OutputType = Tensor::_mat<Real, dim, dim>;
	using Grid = Tensor::Grid<CellType, dim>;
	using ConnGrid = Tensor::Grid<ConnCellType, dim>;
	using CellFieldType = const InputType CellType::*;
	using TensorUSL = typename ConnCellType::TensorUSL;
	using ConnFieldType = const TensorUSL ConnCellType::*;
	static OutputType exec(const Grid &grid, CellFieldType field, const ConnGrid &connGrid, ConnFieldType connField_ull, const Tensor::_vec<Real, dim> &dx, const Tensor::intN<dim> &index) {
		static constexpr int accuracy = 8;	//"order", not "accuracy", and put this all in one place plz
		OutputType result = Tensor::partialDerivativeGrid<accuracy, Real, dim, InputType>(
			index, dx, [&](decltype(index) index) -> InputType {
				return grid(index).*field;
			});

		const CellType &cell = grid(index);
		const ConnCellType &connCell = connGrid(index);
		const TensorUSL &conn_ull = connCell.*connField_ull;
		for (int k = 0; k < dim; ++k) {
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					result(k,i) += conn_ull(i,j,k) * (cell.*field)(j);
				}
			}
		}
		return result;
	}
};

//one-form specialization
//diff_i w_j = partial_i w_j - conn^k_ji w_k
template<typename Real, int dim, typename CellType, typename ConnCellType>
struct CovariantDerivativeClass_L {
	using InputType = Tensor::_vec<Real, dim>;
	using OutputType = Tensor::_mat<Real, dim, dim>;
	using Grid = Tensor::Grid<CellType, dim>;
	using ConnGrid = Tensor::Grid<ConnCellType, dim>;
	using CellFieldType = const InputType CellType::*;
	using TensorUSL = typename ConnCellType::TensorUSL;
	using ConnFieldType = const TensorUSL ConnCellType::*;
	static OutputType exec(const Grid &grid, CellFieldType field, const ConnGrid &connGrid, ConnFieldType connField_ull, const Tensor::_vec<Real, dim> &dx, const Tensor::intN<dim> &index) {
		static constexpr int accuracy = 8;	//"order", not "accuracy", and put this all in one place plz
		OutputType result = Tensor::partialDerivativeGrid<accuracy, Real, dim, InputType>(
			index, dx, [&](decltype(index) index) -> InputType {
				return grid(index).*field;
			});

		const CellType &cell = grid(index);
		const ConnCellType &connCell = connGrid(index);
		const TensorUSL &conn_ull = connCell.*connField_ull;
		for (int k = 0; k < dim; ++k) {
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					result(k,i) -= conn_ull(j,i,k) * (cell.*field)(j);
				}
			}
		}
		return result;
	}
};

//symmetric rank-(2,0) tensor specialization
//diff_i t^jk = partial_i t^jk + conn^j_li t^lk + conn^k_li t^jl
template<typename Real, int dim, typename CellType, typename ConnCellType>
struct CovariantDerivativeClass_SU {
	using InputType = Tensor::_sym<Real, dim>;
	using OutputType = Tensor::_vec<Tensor::_sym<Real, dim>, dim>;
	using Grid = Tensor::Grid<CellType, dim>;
	using ConnGrid = Tensor::Grid<ConnCellType, dim>;
	using CellFieldType = const InputType CellType::*;
	using TensorUSL = typename ConnCellType::TensorUSL;
	using ConnFieldType = const TensorUSL ConnCellType::*;
	static OutputType exec(const Grid &grid, CellFieldType field, const ConnGrid &connGrid, ConnFieldType connField_ull, const Tensor::_vec<Real, dim> &dx, const Tensor::intN<dim> &index) {
		static constexpr int accuracy = 8;	//"order", not "accuracy", and put this all in one place plz
		OutputType result = Tensor::partialDerivativeGrid<accuracy, Real, dim, InputType>(
			index, dx, [&](decltype(index) index) -> InputType {
				return grid(index).*field;
			});

		const CellType &cell = grid(index);
		const ConnCellType &connCell = connGrid(index);
		const TensorUSL &conn_ull = connCell.*connField_ull;
		//<-for the generic covariant definition, from here
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j < dim; ++j) {
				for (int k = 0; k <= j; ++k) {
					//-> to here could be replaced with a write-iterator (note that means not looping over all symmetric indexes)
					for (int l = 0; l < dim; ++l) {
						//...and then these would be replaced with a product-per-index in the source term
						result(i,j,k) += conn_ull(j,l,i) * (cell.*field)(l,k);
						result(i,j,k) += conn_ull(k,l,i) * (cell.*field)(j,l);
					}
				}
			}
		}
		return result;
	}
};
