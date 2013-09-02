#pragma once

#include "vector.h"
#include "tensor.h"
#include "cell.h"
#include "grid.h"
#include "clamp.h"

/*
partial derivative index operator
(partial derivative alone one coordinate)

finite difference coefficients for center-space finite-difference partial derivatives found at
http://en.wikipedia.org/wiki/Finite_difference_coefficients
*/

template<typename real, int accuracy>
struct PartialDerivCoeffs;

template<typename real>
struct PartialDerivCoeffs<real, 2> {
	static const real coeffs[1];
};
template<typename real>
const real PartialDerivCoeffs<real, 2>::coeffs[1] = { 1./2. };

template<typename real>
struct PartialDerivCoeffs<real, 4> {
	static const real coeffs[2];
};
template<typename real>
const real PartialDerivCoeffs<real, 4>::coeffs[2] = { 2./3., -1./12. };

template<typename real>
struct PartialDerivCoeffs<real, 6> {
	static const real coeffs[3];
};
template<typename real>
const real PartialDerivCoeffs<real, 6>::coeffs[3] = { 3./4., -3./20., 1./60. };

template<typename real>
struct PartialDerivCoeffs<real, 8> {
	static const real coeffs[4];
};
template<typename real>
const real PartialDerivCoeffs<real, 8>::coeffs[4] = { 4./5., -1./5., 4./105., -1./280. };

#define numberof(x) (sizeof(x)/sizeof(*(x)))

template<typename real, int dim, typename CellType, typename InputType, int accuracy>
struct PartialDerivativeCoordinateClass {
	typedef PartialDerivCoeffs<real, accuracy> Coeffs;
	InputType operator()(const Grid<CellType, dim> &grid, const InputType CellType::*field, const vector<real, dim> &dx, const vector<int, dim> &index, int k) {
		InputType result = InputType();
		real oneOverDx = 1. / dx(k);
		for (int i = 1; i <= (int)numberof(Coeffs::coeffs); ++i) {
			result += (grid.getOffset(index,k,i).*field - grid.getOffset(index,k,-i).*field) * (Coeffs::coeffs[i-1] * oneOverDx);
		}
		return result;
	}
};
template<typename real, int dim, typename CellType, typename InputType, int accuracy>
InputType partialDerivativeCoordinate(const Grid<CellType, dim> &grid, const InputType CellType::*field, const vector<real, dim> &dx, const vector<int, dim> &index, int k) {
	return PartialDerivativeCoordinateClass<real, dim, CellType, InputType, accuracy>()(grid, field, dx, index, k);
}


template<typename real, int accuracy>
struct PartialSecondDerivCoeffs;

template<typename real>
struct PartialSecondDerivCoeffs<real, 2> {
	static const real coeffs[2];
};
template<typename real>
const real PartialSecondDerivCoeffs<real, 2>::coeffs[2] = { -2., 1. };

template<typename real>
struct PartialSecondDerivCoeffs<real, 4> {
	static const real coeffs[3];
};
template<typename real>
const real PartialSecondDerivCoeffs<real, 4>::coeffs[3] = { -5./2., 4./3., -1./12. };

template<typename real>
struct PartialSecondDerivCoeffs<real, 6> {
	static const real coeffs[4];
};
template<typename real>
const real PartialSecondDerivCoeffs<real, 6>::coeffs[4] = { -49./18., 3./2., -3./20., 1./90. };

template<typename real>
struct PartialSecondDerivCoeffs<real, 8> {
	static const real coeffs[5];
};
template<typename real>
const real PartialSecondDerivCoeffs<real, 8>::coeffs[5] = { -205./72., 8./5., -1./5., 8./315., -1./560. };


template<typename real, int dim, typename CellType, typename InputType, int accuracy>
struct PartialSecondDerivativeCoordinateClass {
	typedef PartialSecondDerivCoeffs<real, accuracy> Coeffs;
	InputType operator()(const Grid<CellType, dim> &grid, const InputType CellType::*field, const vector<real, dim> &dx, const vector<int, dim> &index, int k) {
		real oneOverDxSquared = 1. / (dx(k) * dx(k));
		InputType result = grid(index).*field * (Coeffs::coeffs[0] * oneOverDxSquared);
		for (int i = 1; i < (int)numberof(Coeffs::coeffs); ++i) {
			result += (grid.getOffset(index,k,i).*field + grid.getOffset(index,k,-i).*field) * (Coeffs::coeffs[i] * oneOverDxSquared);  
		}
		return result;
	}
};
template<typename real, int dim, typename CellType, typename InputType, int accuracy>
InputType partialSecondDerivativeCoordinate(const Grid<CellType, dim> &grid, const InputType CellType::*field, const vector<real, dim> &dx, const vector<int, dim> &index, int k) {
	return PartialSecondDerivativeCoordinateClass<real, dim, CellType, InputType, accuracy>()(grid, field, dx, index, k);
}

/*
partial derivative operator
for now let's use 2-point everywhere: d/dx^i f(x) ~= (f(x + dx^i) - f(x - dx^i)) / (2 * |dx^i|)
	index = index in grid of cell to pull the specified field
	k = dimension to differentiate across
*/
template<typename real, int dim, typename CellType, typename InputType>
struct PartialDerivativeClass;

template<typename real, int dim, typename CellType, typename... args>
struct PartialDerivativeClass<real, dim, CellType, tensor<real, args...>> {
	typedef tensor<real, args...> InputType;
	typedef tensor<real, lower<dim>, args...> OutputType;
	typedef ::Grid<CellType, dim> Grid;
	typedef const InputType CellType::*CellFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const vector<real, dim> &dx, const vector<int, dim> &index) {
		OutputType result;
		for (int k = 0; k < dim; ++k) {
			result(k) = partialDerivativeCoordinate<real, dim, CellType, InputType, 8>(grid, field, dx, index, k); 
		}
		return result;
	}
};

template<typename real, int dim, typename CellType>
struct PartialDerivativeClass<real, dim, CellType, real> {
	typedef real InputType;
	typedef tensor<real, lower<dim>> OutputType;
	typedef ::Grid<CellType, dim> Grid;
	typedef const InputType CellType::*CellFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const vector<real, dim> &dx, const vector<int, dim> &index) {
		OutputType result;
		for (int k = 0; k < dim; ++k) {
			result(k) = partialDerivativeCoordinate<real, dim, CellType, InputType, 8>(grid, field, dx, index, k);
		}
		return result;
	}
};

template<typename real, int dim, typename CellType, typename InputType>
typename PartialDerivativeClass<real, dim, CellType, InputType>::OutputType
partialDerivative(
	const Grid<CellType, dim> &grid, 
	const InputType CellType::* field, 
	const vector<real, dim> &dx, 
	const vector<int, dim> &index) 
{
	return PartialDerivativeClass<real, dim, CellType, InputType>()(grid, field, dx, index);
}

/*
partial second derivatives of vectors
until I solve a few things with dereferencing, these will have to be specialized to
*/

template<typename real, int dim, typename CellType, typename PartialCellType, typename InputType>
struct PartialSecondDerivativeClass;

template<typename real, int dim, typename CellType, typename PartialCellType>
struct PartialSecondDerivativeClass<real, dim, CellType, PartialCellType, real> {
	typedef real InputType;
	typedef tensor<real, lower<dim>> PartialInputType;
	typedef tensor<real, symmetric<lower<dim>, lower<dim>>> OutputType;
	typedef ::Grid<CellType, dim> Grid;
	typedef ::Grid<PartialCellType, dim> PartialGrid;
	typedef const InputType CellType::*CellFieldType;
	typedef const PartialInputType PartialCellType::*CellPartialFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const PartialGrid &partialGrid, CellPartialFieldType partialField, const vector<real, dim> &dx, const vector<int, dim> &index) {
		OutputType result;
			
		//partial2_field_ll(i,j) = partial_i partial_j field
		for (int i = 0; i < dim; ++i) {
			PartialInputType partial_field_l_wrt_xi = partialDerivativeCoordinate<real, dim, PartialCellType, PartialInputType, 8>(partialGrid, partialField, dx, index, i);
			for (int j = 0; j <= i; ++j) {
				if (i == j) {
					result(i,j) = partialSecondDerivativeCoordinate<real, dim, CellType, InputType, 8>(grid, field, dx, index, i);
				} else {
					result(i,j) = partial_field_l_wrt_xi(j);
				}
			}
		}

		return result;
	}
};

template<typename real, int dim, typename CellType, typename PartialCellType>
struct PartialSecondDerivativeClass<real, dim, CellType, PartialCellType, tensor<real, symmetric<lower<dim>, lower<dim>>>> {
	typedef tensor<real, symmetric<lower<dim>, lower<dim>>> InputType;
	typedef tensor<real, lower<dim>, symmetric<lower<dim>, lower<dim>>> PartialInputType;
	typedef tensor<real, symmetric<lower<dim>, lower<dim>>, symmetric<lower<dim>, lower<dim>>> OutputType;
	typedef ::Grid<CellType, dim> Grid;
	typedef ::Grid<PartialCellType, dim> PartialGrid;
	typedef const InputType CellType::*CellFieldType;
	typedef const PartialInputType CellType::*CellPartialFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const PartialGrid &partialGrid, CellPartialFieldType partialField, const vector<real, dim> &dx, const vector<int, dim> &index) {
		OutputType result;
			
		//partial2_field_llll(k,l,i,j) = partial_k partial_l field_ij
		for (int k = 0; k < dim; ++k) {
			PartialInputType partial_field_lll_wrt_xk = partialDerivativeCoordinate<real, dim, PartialCellType, PartialInputType, 8>(partialGrid, partialField, dx, index, k);
			for (int l = 0; l <= k; ++l) {
				if (k == l) {
#if 0				// working on subtensor dereferences
					InputType partial2_field_ll_wrt_kk = partialSecondDerivativeCoordinate<real, dim, CellType, InputType, 8>(grid, field, dx, index, k);
					result(k,l) = partial2_field_ll_wrt_kk;
#endif
					
#if 1				//per-index in the mean time
					InputType partial2_field_ll_wrt_kk = partialSecondDerivativeCoordinate<real, dim, CellType, InputType, 8>(grid, field, dx, index, k);
					for (int i = 0; i < dim; ++i) {
						for (int j = 0; j <= i; ++j) {
							result(k,l,i,j) = partial2_field_ll_wrt_kk(i,j);
						}
					}
#endif
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

template<typename real, int dim, typename CellType, typename PartialCellType, typename InputType>
typename PartialSecondDerivativeClass<real, dim, CellType, PartialCellType, InputType>::OutputType
partialSecondDerivative(
	const Grid<CellType, dim> &grid,
	const InputType CellType::*field,
	const Grid<PartialCellType, dim> &partialGrid,
	typename PartialSecondDerivativeClass<real, dim, CellType, PartialCellType, InputType>::CellPartialFieldType partialField,
	const vector<real, dim> &dx,
	const vector<int, dim> &index)
{
	return PartialSecondDerivativeClass<real, dim, CellType, PartialCellType, InputType>()(grid, field, partialGrid, partialField, dx, index);
}


/*
covariant derivative
depends on conn_ull, partialDerivative, and the upper/lower indexes of the tensor

looks like symmetric lower and symmetric upper produce like results,
but symmetric mixed does not

I think I'll make specializations for that reason ...
*/
template<typename real, int dim, typename CellType, typename ConnCellType, typename InputType>
struct CovariantDerivativeClass;

//scalar specialization
//diff_i alpha = partial_i alpha
template<typename real, int dim, typename CellType, typename ConnCellType>
struct CovariantDerivativeClass<real, dim, CellType, ConnCellType, real> {
	typedef real InputType;
	typedef tensor<real, lower<dim>> OutputType;
	typedef ::Grid<CellType, dim> Grid;
	typedef ::Grid<ConnCellType, dim> ConnGrid;
	typedef const InputType CellType::*CellFieldType;
	typedef typename ConnCellType::tensor_usl tensor_usl;
	typedef const tensor_usl ConnCellType::*ConnFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const ConnGrid &connGrid, ConnFieldType connField_ull, const vector<real, dim> &dx, const vector<int, dim> &index) {
		return partialDerivative(grid, field, dx, index);
	}
};

//vector specialization
//diff_i v^j = partial_i v^j + conn^j_ki v^k
template<typename real, int dim, typename CellType, typename ConnCellType>
struct CovariantDerivativeClass<real, dim, CellType, ConnCellType, tensor<real, upper<dim>>> {
	typedef tensor<real, upper<dim>> InputType;
	typedef tensor<real, lower<dim>, upper<dim>> OutputType;
	typedef ::Grid<CellType, dim> Grid;
	typedef ::Grid<ConnCellType, dim> ConnGrid;
	typedef const InputType CellType::*CellFieldType;
	typedef typename ConnCellType::tensor_usl tensor_usl;
	typedef const tensor_usl ConnCellType::*ConnFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const ConnGrid &connGrid, ConnFieldType connField_ull, const vector<real, dim> &dx, const vector<int, dim> &index) {
		OutputType result = partialDerivative(grid, field, dx, index);
		const CellType &cell = grid(index);
		const ConnCellType &connCell = connGrid(index);
		const tensor_usl &conn_ull = connCell.*connField_ull;
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
template<typename real, int dim, typename CellType, typename ConnCellType>
struct CovariantDerivativeClass<real, dim, CellType, ConnCellType, tensor<real, lower<dim>>> {
	typedef tensor<real, lower<dim>> InputType;
	typedef tensor<real, lower<dim>, lower<dim>> OutputType;
	typedef ::Grid<CellType, dim> Grid;
	typedef ::Grid<ConnCellType, dim> ConnGrid;
	typedef const InputType CellType::*CellFieldType;
	typedef typename ConnCellType::tensor_usl tensor_usl;
	typedef const tensor_usl ConnCellType::*ConnFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const ConnGrid &connGrid, ConnFieldType connField_ull, const vector<real, dim> &dx, const vector<int, dim> &index) {
		OutputType result = partialDerivative(grid, field, dx, index);
		const CellType &cell = grid(index);
		const ConnCellType &connCell = connGrid(index);
		const tensor_usl &conn_ull = connCell.*connField_ull;
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
template<typename real, int dim, typename CellType, typename ConnCellType>
struct CovariantDerivativeClass<real, dim, CellType, ConnCellType, tensor<real, symmetric<upper<dim>, upper<dim>>>> {
	typedef tensor<real, symmetric<upper<dim>, upper<dim>>> InputType;
	typedef tensor<real, lower<dim>, symmetric<upper<dim>, upper<dim>>> OutputType;
	typedef ::Grid<CellType, dim> Grid;
	typedef ::Grid<ConnCellType, dim> ConnGrid;
	typedef const InputType CellType::*CellFieldType;
	typedef typename ConnCellType::tensor_usl tensor_usl;
	typedef const tensor_usl ConnCellType::*ConnFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const ConnGrid &connGrid, ConnFieldType connField_ull, const vector<real, dim> &dx, const vector<int, dim> &index) {
		OutputType result = partialDerivative(grid, field, dx, index);
		const CellType &cell = grid(index);
		const ConnCellType &connCell = connGrid(index);
		const tensor_usl &conn_ull = connCell.*connField_ull;
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


template<typename real, int dim, typename CellType, typename ConnCellType, typename InputType>
typename CovariantDerivativeClass<real, dim, CellType, ConnCellType, InputType>::OutputType
covariantDerivative(
	const Grid<CellType, dim> &grid, 
	const InputType CellType::*field, 
	const Grid<ConnCellType, dim> &connGrid,
	const tensor<real, upper<dim>, symmetric<lower<dim>, lower<dim>>> ConnCellType::*connField,
	const vector<real, dim> &dx, 
	const vector<int, dim> &index)
{
	return CovariantDerivativeClass<real, dim, CellType, ConnCellType, InputType>()(grid, field, connGrid, connField, dx, index);
}

