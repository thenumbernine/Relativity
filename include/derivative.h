#pragma once

#include "cell.h"
#include "Tensor/clamp.h"
#include "Tensor/Vector.h"
#include "Tensor/Tensor.h"
#include "Tensor/Grid.h"

//handles offset up to -size of the grid
template<typename Type, int rank> 
const Type& getGridOffset(const Tensor::Grid<Type, rank>& grid, Tensor::Vector<int, rank> index, int dim, int offset) {
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

#define numberof(x) (sizeof(x)/sizeof(*(x)))

template<typename Real, int dim, typename CellType, typename InputType, int accuracy>
struct PartialDerivativeCoordinateClass {
	typedef PartialDerivCoeffs<Real, accuracy> Coeffs;
	InputType operator()(const Tensor::Grid<CellType, dim> &grid, const InputType CellType::*field, const Tensor::Vector<Real, dim> &dx, const Tensor::Vector<int, dim> &index, int k) {
		InputType result = InputType();
		Real oneOverDx = 1. / dx(k);
		for (int i = 1; i <= (int)numberof(Coeffs::coeffs); ++i) {
			result += (getGridOffset(grid, index,k,i).*field - getGridOffset(grid, index,k,-i).*field) * (Coeffs::coeffs[i-1] * oneOverDx);
		}
		return result;
	}
};
template<typename Real, int dim, typename CellType, typename InputType, int accuracy>
InputType partialDerivativeCoordinate(const Tensor::Grid<CellType, dim> &grid, const InputType CellType::*field, const Tensor::Vector<Real, dim> &dx, const Tensor::Vector<int, dim> &index, int k) {
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
	typedef PartialSecondDerivCoeffs<Real, accuracy> Coeffs;
	InputType operator()(const Tensor::Grid<CellType, dim> &grid, const InputType CellType::*field, const Tensor::Vector<Real, dim> &dx, const Tensor::Vector<int, dim> &index, int k) {
		Real oneOverDxSquared = 1. / (dx(k) * dx(k));
		InputType result = grid(index).*field * (Coeffs::coeffs[0] * oneOverDxSquared);
		for (int i = 1; i < (int)numberof(Coeffs::coeffs); ++i) {
			result += (getGridOffset(grid, index,k,i).*field + getGridOffset(grid, index,k,-i).*field) * (Coeffs::coeffs[i] * oneOverDxSquared);  
		}
		return result;
	}
};
template<typename Real, int dim, typename CellType, typename InputType, int accuracy>
InputType partialSecondDerivativeCoordinate(const Tensor::Grid<CellType, dim> &grid, const InputType CellType::*field, const Tensor::Vector<Real, dim> &dx, const Tensor::Vector<int, dim> &index, int k) {
	return PartialSecondDerivativeCoordinateClass<Real, dim, CellType, InputType, accuracy>()(grid, field, dx, index, k);
}

/*
partial derivative operator
for now let's use 2-point everywhere: d/dx^i f(x) ~= (f(x + dx^i) - f(x - dx^i)) / (2 * |dx^i|)
	index = index in grid of cell to pull the specified field
	k = dimension to differentiate across
*/
template<typename Real, int dim, typename CellType, typename InputType>
struct PartialDerivativeClass;

template<typename Real, int dim, typename CellType, typename... args>
struct PartialDerivativeClass<Real, dim, CellType, Tensor::Tensor<Real, args...>> {
	typedef Tensor::Tensor<Real, args...> InputType;
	typedef Tensor::Tensor<Real, Tensor::Lower<dim>, args...> OutputType;
	typedef Tensor::Grid<CellType, dim> Grid;
	typedef const InputType CellType::*CellFieldType;
	enum { rank = InputType::rank };
	OutputType operator()(const Grid &grid, CellFieldType field, const Tensor::Vector<Real, dim> &dx, const Tensor::Vector<int, dim> &gridIndex) {
		enum { accuracy = 8 };
		typedef PartialDerivCoeffs<Real, accuracy> Coeffs;
		return OutputType([&](Tensor::Vector<int, rank+1> dstIndex){
			int gradIndex = dstIndex(0);
			Tensor::Vector<int, rank> srcIndex;
			for (int i = 0; i < rank; ++i) {
				srcIndex(i) = dstIndex(i+1);
			}
			Real sum = Real();
			for (int i = 1; i <= (int)numberof(Coeffs::coeffs); ++i) {
				sum += ((getGridOffset(grid, gridIndex, gradIndex, i).*field)(srcIndex) - (getGridOffset(grid, gridIndex, gradIndex, -1).*field)(srcIndex)) * Coeffs::coeffs[i-1];
			}
			return sum / dx(gradIndex);
		});
	}
};

template<typename Real, int dim, typename CellType>
struct PartialDerivativeClass<Real, dim, CellType, Real> {
	typedef Real InputType;
	typedef Tensor::Tensor<Real, Tensor::Lower<dim>> OutputType;
	typedef Tensor::Grid<CellType, dim> Grid;
	typedef const InputType CellType::*CellFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const Tensor::Vector<Real, dim> &dx, const Tensor::Vector<int, dim> &gridIndex) {
		enum { accuracy = 8 };
		typedef PartialDerivCoeffs<Real, accuracy> Coeffs;
		return OutputType([&](Tensor::Vector<int, 1> dstIndex){
			int gradIndex = dstIndex(0);
			Real sum = 0.f;
			for (int i = 1; i <= (int)numberof(Coeffs::coeffs); ++i) {
				sum += (getGridOffset(grid, gridIndex, gradIndex, i).*field - getGridOffset(grid, gridIndex, gradIndex, -1).*field) * Coeffs::coeffs[i-1];
			}
			return sum / dx(gradIndex);		
		});
	}
};

template<typename Real, int dim, typename CellType, typename InputType>
typename PartialDerivativeClass<Real, dim, CellType, InputType>::OutputType
partialDerivative(
	const Tensor::Grid<CellType, dim> &grid, 
	const InputType CellType::* field, 
	const Tensor::Vector<Real, dim> &dx, 
	const Tensor::Vector<int, dim> &index) 
{
	return PartialDerivativeClass<Real, dim, CellType, InputType>()(grid, field, dx, index);
}

/*
partial second derivatives of vectors
until I solve a few things with dereferencing, these will have to be specialized to
*/

template<typename Real, int dim, typename CellType, typename PartialCellType, typename InputType>
struct PartialSecondDerivativeClass;

template<typename Real, int dim, typename CellType, typename PartialCellType>
struct PartialSecondDerivativeClass<Real, dim, CellType, PartialCellType, Real> {
	typedef Real InputType;
	typedef Tensor::Tensor<Real, Tensor::Lower<dim>> PartialInputType;
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<dim>, Tensor::Lower<dim>>> OutputType;
	typedef Tensor::Grid<CellType, dim> Grid;
	typedef Tensor::Grid<PartialCellType, dim> PartialGrid;
	typedef const InputType CellType::*CellFieldType;
	typedef const PartialInputType PartialCellType::*CellPartialFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const PartialGrid &partialGrid, CellPartialFieldType partialField, const Tensor::Vector<Real, dim> &dx, const Tensor::Vector<int, dim> &index) {
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
struct PartialSecondDerivativeClass<Real, dim, CellType, PartialCellType, Tensor::Tensor<Real, Tensor::Upper<dim>>> {
	typedef Tensor::Tensor<Real, Tensor::Upper<dim>> InputType;
	typedef Tensor::Tensor<Real, Tensor::Lower<dim>, Tensor::Upper<dim>> PartialInputType;
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<dim>, Tensor::Lower<dim>>, Tensor::Upper<dim>> OutputType;
	typedef Tensor::Grid<CellType, dim> Grid;
	typedef Tensor::Grid<PartialCellType, dim> PartialGrid;
	typedef const InputType CellType::*CellFieldType;
	typedef const PartialInputType PartialCellType::*CellPartialFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const PartialGrid &partialGrid, CellPartialFieldType partialField, const Tensor::Vector<Real, dim> &dx, const Tensor::Vector<int, dim> &index) {
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
struct PartialSecondDerivativeClass<Real, dim, CellType, PartialCellType, Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<dim>, Tensor::Lower<dim>>>> {
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<dim>, Tensor::Lower<dim>>> InputType;
	typedef Tensor::Tensor<Real, Tensor::Lower<dim>, Tensor::Symmetric<Tensor::Lower<dim>, Tensor::Lower<dim>>> PartialInputType;
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<dim>, Tensor::Lower<dim>>, Tensor::Symmetric<Tensor::Lower<dim>, Tensor::Lower<dim>>> OutputType;
	typedef Tensor::Grid<CellType, dim> Grid;
	typedef Tensor::Grid<PartialCellType, dim> PartialGrid;
	typedef const InputType CellType::*CellFieldType;
	typedef const PartialInputType PartialCellType::*CellPartialFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const PartialGrid &partialGrid, CellPartialFieldType partialField, const Tensor::Vector<Real, dim> &dx, const Tensor::Vector<int, dim> &index) {
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
	const Tensor::Vector<Real, dim> &dx,
	const Tensor::Vector<int, dim> &index)
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
template<typename Real, int dim, typename CellType, typename ConnCellType, typename InputType>
struct CovariantDerivativeClass;

//scalar specialization
//diff_i alpha = partial_i alpha
template<typename Real, int dim, typename CellType, typename ConnCellType>
struct CovariantDerivativeClass<Real, dim, CellType, ConnCellType, Real> {
	typedef Real InputType;
	typedef Tensor::Tensor<Real, Tensor::Lower<dim>> OutputType;
	typedef Tensor::Grid<CellType, dim> Grid;
	typedef Tensor::Grid<ConnCellType, dim> ConnGrid;
	typedef const InputType CellType::*CellFieldType;
	typedef typename ConnCellType::TensorUSL TensorUSL;
	typedef const TensorUSL ConnCellType::*ConnFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const ConnGrid &connGrid, ConnFieldType connField_ull, const Tensor::Vector<Real, dim> &dx, const Tensor::Vector<int, dim> &index) {
		return partialDerivative(grid, field, dx, index);
	}
};

//Tensor::Vector specialization
//diff_i v^j = partial_i v^j + conn^j_ki v^k
template<typename Real, int dim, typename CellType, typename ConnCellType>
struct CovariantDerivativeClass<Real, dim, CellType, ConnCellType, Tensor::Tensor<Real, Tensor::Upper<dim>>> {
	typedef Tensor::Tensor<Real, Tensor::Upper<dim>> InputType;
	typedef Tensor::Tensor<Real, Tensor::Lower<dim>, Tensor::Upper<dim>> OutputType;
	typedef Tensor::Grid<CellType, dim> Grid;
	typedef Tensor::Grid<ConnCellType, dim> ConnGrid;
	typedef const InputType CellType::*CellFieldType;
	typedef typename ConnCellType::TensorUSL TensorUSL;
	typedef const TensorUSL ConnCellType::*ConnFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const ConnGrid &connGrid, ConnFieldType connField_ull, const Tensor::Vector<Real, dim> &dx, const Tensor::Vector<int, dim> &index) {
		OutputType result = partialDerivative(grid, field, dx, index);
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
struct CovariantDerivativeClass<Real, dim, CellType, ConnCellType, Tensor::Tensor<Real, Tensor::Lower<dim>>> {
	typedef Tensor::Tensor<Real, Tensor::Lower<dim>> InputType;
	typedef Tensor::Tensor<Real, Tensor::Lower<dim>, Tensor::Lower<dim>> OutputType;
	typedef Tensor::Grid<CellType, dim> Grid;
	typedef Tensor::Grid<ConnCellType, dim> ConnGrid;
	typedef const InputType CellType::*CellFieldType;
	typedef typename ConnCellType::TensorUSL TensorUSL;
	typedef const TensorUSL ConnCellType::*ConnFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const ConnGrid &connGrid, ConnFieldType connField_ull, const Tensor::Vector<Real, dim> &dx, const Tensor::Vector<int, dim> &index) {
		OutputType result = partialDerivative(grid, field, dx, index);
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
struct CovariantDerivativeClass<Real, dim, CellType, ConnCellType, Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Upper<dim>, Tensor::Upper<dim>>>> {
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Upper<dim>, Tensor::Upper<dim>>> InputType;
	typedef Tensor::Tensor<Real, Tensor::Lower<dim>, Tensor::Symmetric<Tensor::Upper<dim>, Tensor::Upper<dim>>> OutputType;
	typedef Tensor::Grid<CellType, dim> Grid;
	typedef Tensor::Grid<ConnCellType, dim> ConnGrid;
	typedef const InputType CellType::*CellFieldType;
	typedef typename ConnCellType::TensorUSL TensorUSL;
	typedef const TensorUSL ConnCellType::*ConnFieldType;
	OutputType operator()(const Grid &grid, CellFieldType field, const ConnGrid &connGrid, ConnFieldType connField_ull, const Tensor::Vector<Real, dim> &dx, const Tensor::Vector<int, dim> &index) {
		OutputType result = partialDerivative(grid, field, dx, index);
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


template<typename Real, int dim, typename CellType, typename ConnCellType, typename InputType>
typename CovariantDerivativeClass<Real, dim, CellType, ConnCellType, InputType>::OutputType
covariantDerivative(
	const Tensor::Grid<CellType, dim> &grid, 
	const InputType CellType::*field, 
	const Tensor::Grid<ConnCellType, dim> &connGrid,
	const Tensor::Tensor<Real, Tensor::Upper<dim>, Tensor::Symmetric<Tensor::Lower<dim>, Tensor::Lower<dim>>> ConnCellType::*connField,
	const Tensor::Vector<Real, dim> &dx, 
	const Tensor::Vector<int, dim> &index)
{
	return CovariantDerivativeClass<Real, dim, CellType, ConnCellType, InputType>()(grid, field, connGrid, connField, dx, index);
}

