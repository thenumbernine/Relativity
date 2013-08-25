#pragma once

#include "vector.h"
#include "tensor.h"
#include "cell.h"
#include "grid.h"

/*
partial derivative index operator
(partial derivative alone one coordinate)
*/
template<int dim, typename real, typename T>
T partialDerivativeCoordinate(Grid<dim, Cell<dim,real>> &cells, T Cell<dim,real>::*field, const vector<dim,real> &dx, const vector<dim,int> &index, int k) {
	typedef vector<dim,int> deref_type;
	deref_type nextIndex = index;
	deref_type prevIndex = index;
	nextIndex(k) = std::max(0, std::min(cells.size(k)-1, index(k) + 1));
	prevIndex(k) = std::max(0, std::min(cells.size(k)-1, index(k) - 1));
	real dxk = real(nextIndex(k) - prevIndex(k)) * dx(k);
	T dfield = cells(nextIndex).*field - cells(prevIndex).*field;
	return dfield / dxk;
}

/*
partial derivative operator
for now let's use 2-point everywhere: d/dx^i f(x) ~= (f(x + dx^i) - f(x - dx^i)) / (2 * |dx^i|)
	index = index in grid of cell to pull the specified field
	k = dimension to differentiate across
*/
template<int dim, typename real, typename input_type>
struct partialDerivativeClass;

template<int dim, typename real, typename... args>
struct partialDerivativeClass<
	dim,
	real,
	tensor<real, args...>
>{
	typedef tensor<real, args...> input_type;
	typedef tensor<real, lower<dim>, args...> output_type;
	typedef ::Cell<dim,real> Cell;
	typedef ::Grid<dim,Cell> Grid;
	typedef input_type Cell::*cell_field_type;
	output_type operator()(Grid &cells, cell_field_type field, const vector<dim,real> &dx, const vector<dim,int> &index) {
		output_type result;
#if 0	//sub-tensor access: all-at-once
		vector<dim,int> nextIndex = vector<dim,int>::clamp(index+1, vector<dim,int>(0), size-1);
		vector<dim,int> prevIndex = vector<dim,int>::clamp(index-1, vector<dim,int>(0), size-1);
		const input_type &nextT = cells(nextIndex).*field;
		const input_type &prevT = cells(prevIndex).*field;
		for (int k = 0; k < dim; ++k) {
			result(k) = (nextT - prevT) / (real(nextIndex(k) - prevIndex(k)) * dx(k));
		}
#elif 1	//sub-tensor access: per-dimension. has the benefit of working with symmetric structures.
		for (int k = 0; k < dim; ++k) {
			result(k) = partialDerivativeCoordinate(cells, field, dx, index, k); 
		}
#else	//per-index.
		//I got rid of sub-tensor accessor.
		//not sure if it is the right thing to expose, and covariant diffs have to march through indexes anyways.
		vector<dim,int> nextIndex = vector<dim,int>::clamp(index+1, vector<dim,int>(0), size-1);
		vector<dim,int> prevIndex = vector<dim,int>::clamp(index-1, vector<dim,int>(0), size-1);
		//these should be const-ref, but that'd mean making a separate const-access tensor =P
		input_type &nextT = cells(nextIndex).*field;
		input_type &prevT = cells(prevIndex).*field;
		typename output_type::deref_type output_deref;
		for (int k = 0; k < dim; ++k) {
			output_deref(0) = k;
			//this is where sub-tensor is more efficient -- or I should make separate read- and write- iterators.
			// symmetric and antisymmetric indexes don't have to march through all members,
			// so for component-wise operations like this, I could just use thie write iterator
			for (typename input_type::iterator i = nextT.begin(); i != nextT.end(); ++i) {
				for (int j = 1; j < output_deref.dim; ++j) {
					output_deref(j) = i.index(j-1);	//sub-vector dereferencing.  good idea?
				}
				result(output_deref) = (nextT(i.index) - prevT(i.index)) / (real(nextIndex(k) - prevIndex(k)) * dx(k));
			}
		}
#endif
		return result;
	}
};

template<int dim, typename real>
struct partialDerivativeClass<
	dim,
	real,
	real
>{
	typedef real input_type;
	typedef tensor<real, lower<dim>> output_type;
	typedef ::Cell<dim,real> Cell;
	typedef ::Grid<dim,Cell> Grid;
	typedef input_type Cell::*cell_field_type;
	output_type operator()(Grid &cells, cell_field_type field, const vector<dim,real> &dx, const vector<dim,int> &index) {
		output_type result;
		for (int k = 0; k < dim; ++k) {
			result(k) = partialDerivativeCoordinate(cells, field, dx, index, k);
		}
		return result;
	}
};

template<int dim, typename real, typename input_type>
typename partialDerivativeClass<dim,real,input_type>::output_type
partialDerivative(
	Grid<dim,Cell<dim,real>> &cells, 
	input_type Cell<dim,real>::* field, 
	const vector<dim,real> &dx, 
	const vector<dim,int> &index) 
{
	return partialDerivativeClass<dim,real,input_type>()(cells, field, dx, index);
}

/*
covariant derivative
depends on conn_ull, partialDerivative, and the upper/lower indexes of the tensor

looks like symmetric lower and symmetric upper produce like results,
but symmetric mixed does not

I think I'll make specializations for that reason ...
*/
template<int dim, typename real, typename input_type>
struct covariantDerivativeClass;

//scalar specialization
template<int dim, typename real>
struct covariantDerivativeClass<
	dim,
	real,
	real
>{
	typedef real input_type;
	typedef tensor<real, lower<dim>> output_type;
	typedef ::Cell<dim,real> Cell;
	typedef ::Grid<dim,Cell> Grid;
	typedef input_type Cell::*cell_field_type;
	output_type operator()(Grid &cells, cell_field_type field, const vector<dim,real> &dx, const vector<dim,int> &index) {
		return partialDerivative(cells, field, dx, index);
	}
};

//vector specialization
template<int dim, typename real>
struct covariantDerivativeClass<
	dim,
	real,
	tensor<real, upper<dim>>
>{
	typedef tensor<real, upper<dim>> input_type;
	typedef tensor<real, lower<dim>, upper<dim>> output_type;
	typedef ::Cell<dim,real> Cell;
	typedef ::Grid<dim,Cell> Grid;
	typedef input_type Cell::*cell_field_type;
	output_type operator()(Grid &cells, cell_field_type field, const vector<dim,real> &dx, const vector<dim,int> &index) {
		output_type result = partialDerivative(cells, field, dx, index);
		for (int k = 0; k < dim; ++k) {
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					result(k,i) += cells.conn_ull(i,j,k) * (cells.*field)(j);
				}
			}
		}
		return result;
	}
};

template<int dim, typename real>
struct covariantDerivativeClass<
	dim,
	real,
	tensor<real, lower<dim>>
>{
	typedef tensor<real, lower<dim>> input_type;
	typedef tensor<real, lower<dim>, lower<dim>> output_type;
	typedef ::Cell<dim,real> Cell;
	typedef ::Grid<dim,Cell> Grid;
	typedef input_type Cell::*cell_field_type;
	output_type operator()(Grid &cells, cell_field_type field, const vector<dim,real> &dx, const vector<dim,int> &index) {
		output_type result = partialDerivative(cells, field, dx, index);
		Cell &cell = cells(index);
		for (int k = 0; k < dim; ++k) {
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					result(k,i) -= cell.conn_ull(j,i,k) * (cell.*field)(j);
				}
			}
		}
		return result;
	}
};

template<int dim, typename real, typename input_type>
typename covariantDerivativeClass<dim,real,input_type>::output_type covariantDerivative(Grid<dim, Cell<dim,real>> &cells, input_type Cell<dim,real>::*field, const vector<dim,real> &dx, const vector<dim,int> &index) {
	return covariantDerivativeClass<dim,real,input_type>()(cells, field, dx, index);
}

