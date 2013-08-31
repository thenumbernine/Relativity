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
template<typename real, int dim, typename InputType, int accuracy>
struct PartialDerivativeCoordinateClass;

template<typename real, int dim, typename InputType>
struct PartialDerivativeCoordinateClass<real, dim, InputType, 2> {
	InputType operator()(Grid<Cell<real, dim>, dim> &cells, const InputType Cell<real, dim>::*field, const vector<real, dim> &dx, const vector<int, dim> &index, int k) {
		return (
			cells.getOffset(index,k,1).*field 
			- cells.getOffset(index,k,-1).*field) 
				* (1. / (2. * dx(k)));
	}
};

template<typename real, int dim, typename InputType>
struct PartialDerivativeCoordinateClass<real, dim, InputType, 4> {
	InputType operator()(Grid<Cell<real, dim>, dim> &cells, const InputType Cell<real, dim>::*field, const vector<real, dim> &dx, const vector<int, dim> &index, int k) {
		return (
			-(cells.getOffset(index,k,2).*field) 
			+ cells.getOffset(index,k,1).*field * 8. 
			- cells.getOffset(index,k,-1).*field * 8. 
			+ cells.getOffset(index,k,-2).*field) 
				* (1. / (12. * dx(k)));
	}
};

template<typename real, int dim, typename InputType>
struct PartialDerivativeCoordinateClass<real, dim, InputType, 6> {
	InputType operator()(Grid<Cell<real, dim>, dim> &cells, const InputType Cell<real, dim>::*field, const vector<real, dim> &dx, const vector<int, dim> &index, int k) {
		return (
			cells.getOffset(index,k,3).*field
			- cells.getOffset(index,k,2).*field * 9.
			+ cells.getOffset(index,k,1).*field * 45.
			- cells.getOffset(index,k,-1).*field * 45.
			+ cells.getOffset(index,k,-2).*field * 9.
			- cells.getOffset(index,k,-3).*field)
				* (1. / (60. * dx(k)));
	}
};

template<typename real, int dim, typename InputType>
struct PartialDerivativeCoordinateClass<real, dim, InputType, 8> {
	InputType operator()(Grid<Cell<real, dim>, dim> &cells, const InputType Cell<real, dim>::*field, const vector<real, dim> &dx, const vector<int, dim> &index, int k) {
		return (
			-(cells.getOffset(index,k,4).*field) * 3.
			+ cells.getOffset(index,k,3).*field * 32.
			- cells.getOffset(index,k,2).*field * 168.
			+ cells.getOffset(index,k,1).*field * 672.
			- cells.getOffset(index,k,-1).*field * 672.
			+ cells.getOffset(index,k,-2).*field * 168.
			- cells.getOffset(index,k,-3).*field * 32.
			+ cells.getOffset(index,k,-4).*field * 3.)
				* (1. / (840. * dx(k)));
	}
};

template<typename real, int dim, typename InputType>
InputType partialDerivativeCoordinate(Grid<Cell<real, dim>, dim> &cells, InputType Cell<real, dim>::*field, const vector<real, dim> &dx, const vector<int, dim> &index, int k) {
	return PartialDerivativeCoordinateClass<real, dim, InputType, 8>()(cells, field, dx, index, k);
}


template<typename real, int dim, typename InputType, int accuracy>
struct PartialSecondDerivativeCoordinateClass;

template<typename real, int dim, typename InputType>
struct PartialSecondDerivativeCoordinateClass<real, dim, InputType, 2> {
	InputType operator()(Grid<Cell<real, dim>, dim> &cells, const InputType Cell<real, dim>::*field, const vector<real, dim> &dx, const vector<int, dim> &index, int k) {
		return (
			cells.getOffset(index,k,1).*field 
			- cells(index).*field * 2.
			+ cells.getOffset(index,k,-1).*field) 
				* (1. / (dx(k) * dx(k)));
	}
};

template<typename real, int dim, typename InputType>
struct PartialSecondDerivativeCoordinateClass<real, dim, InputType, 4> {
	InputType operator()(Grid<Cell<real, dim>, dim> &cells, const InputType Cell<real, dim>::*field, const vector<real, dim> &dx, const vector<int, dim> &index, int k) {
		return (
			-(cells.getOffset(index,k,2).*field) 
			+ cells.getOffset(index,k,1).*field * 16. 
			- cells(index).*field * 30.
			+ cells.getOffset(index,k,-1).*field * 16. 
			- cells.getOffset(index,k,-2).*field) 
				* (1. / (12. * dx(k) * dx(k)));
	}
};

template<typename real, int dim, typename InputType>
struct PartialSecondDerivativeCoordinateClass<real, dim, InputType, 6> {
	InputType operator()(Grid<Cell<real, dim>, dim> &cells, const InputType Cell<real, dim>::*field, const vector<real, dim> &dx, const vector<int, dim> &index, int k) {
		return (
			cells.getOffset(index,k,3).*field * 2.
			- cells.getOffset(index,k,2).*field * 27.
			+ cells.getOffset(index,k,1).*field * 270.
			- cells(index).*field * 490.
			+ cells.getOffset(index,k,-1).*field * 270.
			- cells.getOffset(index,k,-2).*field * 27.
			+ cells.getOffset(index,k,-3).*field * 2.)
				* (1. / (180. * dx(k) * dx(k)));
	}
};

template<typename real, int dim, typename InputType>
struct PartialSecondDerivativeCoordinateClass<real, dim, InputType, 8> {
	InputType operator()(Grid<Cell<real, dim>, dim> &cells, const InputType Cell<real, dim>::*field, const vector<real, dim> &dx, const vector<int, dim> &index, int k) {
		return (
			-(cells.getOffset(index,k,4).*field) * 9.
			+ cells.getOffset(index,k,3).*field * 128.
			- cells.getOffset(index,k,2).*field * 1008.
			+ cells.getOffset(index,k,1).*field * 8064.
			- cells(index).*field * 14350.
			+ cells.getOffset(index,k,-1).*field * 8064.
			- cells.getOffset(index,k,-2).*field * 1008.
			+ cells.getOffset(index,k,-3).*field * 128.
			- cells.getOffset(index,k,-4).*field * 9.)
				* (1. / (5040. * dx(k) * dx(k)));
	}
};

template<typename real, int dim, typename InputType>
InputType partialSecondDerivativeCoordinate(Grid<Cell<real, dim>, dim> &cells, InputType Cell<real, dim>::*field, const vector<real, dim> &dx, const vector<int, dim> &index, int k) {
	return PartialSecondDerivativeCoordinateClass<real, dim, InputType, 8>()(cells, field, dx, index, k);
}

/*
partial derivative operator
for now let's use 2-point everywhere: d/dx^i f(x) ~= (f(x + dx^i) - f(x - dx^i)) / (2 * |dx^i|)
	index = index in grid of cell to pull the specified field
	k = dimension to differentiate across
*/
template<typename real, int dim, typename InputType>
struct partialDerivativeClass;

template<typename real, int dim, typename... args>
struct partialDerivativeClass<real, dim, tensor<real, args...>> {
	typedef tensor<real, args...> InputType;
	typedef tensor<real, lower<dim>, args...> OutputType;
	typedef ::Cell<real, dim> Cell;
	typedef ::Grid<Cell, dim> Grid;
	typedef const InputType Cell::*CellFieldType;
	OutputType operator()(Grid &cells, CellFieldType field, const vector<real, dim> &dx, const vector<int, dim> &index) {
		OutputType result;
#if 0	//sub-tensor access: all-at-once
		vector<int, dim> nextIndex = vector<int, dim>::clamp(index+1, vector<int, dim>(0), size-1);
		vector<int, dim> prevIndex = vector<int, dim>::clamp(index-1, vector<int, dim>(0), size-1);
		const InputType &nextT = cells(nextIndex).*field;
		const InputType &prevT = cells(prevIndex).*field;
		for (int k = 0; k < dim; ++k) {
			result(k) = (nextT - prevT) / (2. * dx(k));
		}
#elif 1	//sub-tensor access: per-dimension. has the benefit of working with symmetric structures.
		for (int k = 0; k < dim; ++k) {
			result(k) = partialDerivativeCoordinate(cells, field, dx, index, k); 
		}
#else	//per-index.
		//I got rid of sub-tensor accessor.
		//not sure if it is the right thing to expose, and covariant diffs have to march through indexes anyways.
		vector<int, dim> nextIndex = vector<int, dim>::clamp(index+1, vector<int, dim>(0), size-1);
		vector<int, dim> prevIndex = vector<int, dim>::clamp(index-1, vector<int, dim>(0), size-1);
		//these should be const-ref, but that'd mean making a separate const-access tensor =P
		InputType &nextT = cells(nextIndex).*field;
		InputType &prevT = cells(prevIndex).*field;
		typename OutputType::DerefType output_deref;
		for (int k = 0; k < dim; ++k) {
			output_deref(0) = k;
			//this is where sub-tensor is more efficient -- or I should make separate read- and write- iterators.
			// symmetric and antisymmetric indexes don't have to march through all members,
			// so for component-wise operations like this, I could just use thie write iterator
			for (typename InputType::iterator i = nextT.begin(); i != nextT.end(); ++i) {
				for (int j = 1; j < output_deref.dim; ++j) {
					output_deref(j) = i.index(j-1);	//sub-vector dereferencing.  good idea?
				}
				result(output_deref) = (nextT(i.index) - prevT(i.index)) / (2. * dx(k));
			}
		}
#endif
		return result;
	}
};

template<typename real, int dim>
struct partialDerivativeClass<real, dim, real> {
	typedef real InputType;
	typedef tensor<real, lower<dim>> OutputType;
	typedef ::Cell<real, dim> Cell;
	typedef ::Grid<Cell, dim> Grid;
	typedef const InputType Cell::*CellFieldType;
	OutputType operator()(Grid &cells, CellFieldType field, const vector<real, dim> &dx, const vector<int, dim> &index) {
		OutputType result;
		for (int k = 0; k < dim; ++k) {
			result(k) = partialDerivativeCoordinate(cells, field, dx, index, k);
		}
		return result;
	}
};

template<typename real, int dim, typename InputType>
typename partialDerivativeClass<real, dim, InputType>::OutputType
partialDerivative(
	Grid<Cell<real, dim>, dim> &cells, 
	const InputType Cell<real, dim>::* field, 
	const vector<real, dim> &dx, 
	const vector<int, dim> &index) 
{
	return partialDerivativeClass<real, dim, InputType>()(cells, field, dx, index);
}

/*
covariant derivative
depends on conn_ull, partialDerivative, and the upper/lower indexes of the tensor

looks like symmetric lower and symmetric upper produce like results,
but symmetric mixed does not

I think I'll make specializations for that reason ...
*/
template<typename real, int dim, typename InputType>
struct CovariantDerivativeClass;

//scalar specialization
//diff_i alpha = partial_i alpha
template<typename real, int dim>
struct CovariantDerivativeClass<real, dim, real> {
	typedef real InputType;
	typedef tensor<real, lower<dim>> OutputType;
	typedef ::Cell<real, dim> Cell;
	typedef ::Grid<Cell, dim> Grid;
	typedef const InputType Cell::*CellFieldType;
	typedef typename Cell::tensor_usl tensor_usl;
	typedef tensor_usl Cell::*ConnFieldType;
	OutputType operator()(Grid &cells, CellFieldType field, const vector<real, dim> &dx, const vector<int, dim> &index, ConnFieldType connField_ull) {
		return partialDerivative(cells, field, dx, index);
	}
};

//vector specialization
//diff_i v^j = partial_i v^j + conn^j_ki v^k
template<typename real, int dim>
struct CovariantDerivativeClass<real, dim, tensor<real, upper<dim>>> {
	typedef tensor<real, upper<dim>> InputType;
	typedef tensor<real, lower<dim>, upper<dim>> OutputType;
	typedef ::Cell<real, dim> Cell;
	typedef ::Grid<Cell, dim> Grid;
	typedef const InputType Cell::*CellFieldType;
	typedef typename Cell::tensor_usl tensor_usl;
	typedef tensor_usl Cell::*ConnFieldType;
	OutputType operator()(Grid &cells, CellFieldType field, const vector<real, dim> &dx, const vector<int, dim> &index, ConnFieldType connField_ull) {
		OutputType result = partialDerivative(cells, field, dx, index);
		Cell &cell = cells(index);
		const tensor_usl &conn_ull = cell.*connField_ull;
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
template<typename real, int dim>
struct CovariantDerivativeClass<real, dim, tensor<real, lower<dim>>> {
	typedef tensor<real, lower<dim>> InputType;
	typedef tensor<real, lower<dim>, lower<dim>> OutputType;
	typedef ::Cell<real, dim> Cell;
	typedef ::Grid<Cell, dim> Grid;
	typedef const InputType Cell::*CellFieldType;
	typedef typename Cell::tensor_usl tensor_usl;
	typedef tensor_usl Cell::*ConnFieldType;
	OutputType operator()(Grid &cells, CellFieldType field, const vector<real, dim> &dx, const vector<int, dim> &index, ConnFieldType connField_ull) {
		OutputType result = partialDerivative(cells, field, dx, index);
		Cell &cell = cells(index);
		const tensor_usl &conn_ull = cell.*connField_ull;
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
template<typename real, int dim>
struct CovariantDerivativeClass<real, dim, tensor<real, symmetric<upper<dim>, upper<dim>>>> {
	typedef tensor<real, symmetric<upper<dim>, upper<dim>>> InputType;
	typedef tensor<real, lower<dim>, symmetric<upper<dim>, upper<dim>>> OutputType;
	typedef ::Cell<real, dim> Cell;
	typedef ::Grid<Cell, dim> Grid;
	typedef const InputType Cell::*CellFieldType;
	typedef typename Cell::tensor_usl tensor_usl;
	typedef tensor_usl Cell::*ConnFieldType;
	OutputType operator()(Grid &cells, CellFieldType field, const vector<real, dim> &dx, const vector<int, dim> &index, ConnFieldType connField_ull) {
		OutputType result = partialDerivative(cells, field, dx, index);
		Cell &cell = cells(index);
		const tensor_usl &conn_ull = cell.*connField_ull;
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


template<typename real, int dim, typename InputType>
typename CovariantDerivativeClass<real, dim, InputType>::OutputType covariantDerivative(
	Grid<Cell<real, dim>, dim> &cells, 
	const InputType Cell<real, dim>::*field, 
	const vector<real, dim> &dx, 
	const vector<int, dim> &index,
	typename Cell<real, dim>::tensor_usl Cell<real, dim>::*connField)
{
	return CovariantDerivativeClass<real, dim, InputType>()(cells, field, dx, index, connField);
}

