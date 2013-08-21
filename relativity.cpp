#include <assert.h>
#include <algorithm>	//min/max

#include "vec.h"
#include "symmat.h"
#include "matrix.h"
#include "cell.h"
#include "grid.h"

template<int dim_, typename real_>
struct Simulation {
	enum { dim = dim_ };
	typedef real_ real;

	//statically-sized mathematical types 
	typedef ::vec<dim,real> vec;		//3D hypersurface vector
	typedef ::vec<dim,int> veci;		//integer vector, used for indexes
	typedef ::vec<dim+1,real> tvec;		//4D hypersurface + timelike component vector
	typedef ::matrix<dim,real> matrix;	//3D matrix / rank-2 tensor of some sort
	typedef ::symmat<dim,real> symmat;	//3D symmetric matrix / rank-2 tensor of some sort

	typedef ::Cell<dim,real> Cell;
	typedef ::Grid<dim,Cell> CellGrid;
	CellGrid readCells, writeCells;

	//kronecher delta for integer vector / indexes, such that dxi(i)(j) == i == j
	::vec<dim,veci> dxi;

	//resolution of our grids, stored here as well as in each grid for convenience
	veci size;

	//x(0) == min
	vec min;
	
	//x(size-1) == max
	vec max;
	
	//range = max - min
	vec range;
	
	//dx = range / size
	vec dx;

	Simulation(const vec &min_, const vec &max_, const veci &size_)
	:	size(size_),
		min(min_),
		max(max_),
		range(max_ - min_),
		dx((max_ - min_) / vec(size_)),
		readCells(size_),
		writeCells(size_)
	{
		//basis vectors as integers for indexes
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j < dim; ++j) {
				dxi(i)(j) == i == j;
			}
		}
	}

	/*
	partial derivative operator
	for now let's use 2-point everywhere: d/dx^i f(x) ~= (f(x + dx^i) - f(x - dx^i)) / (2 * |dx^i|)
		index = index in grid of cell to pull the specified field
		k = dimension to differentiate across
	*/
	template<typename T>
	T partialDerivative(T Cell::*field, const veci &index, int k) {
		veci nextIndex = veci::clamp(index + dxi(dim), veci(0), size-1);
		veci prevIndex = veci::clamp(index - dxi(dim), veci(0), size-1);
		return (readCells(nextIndex).*field - readCells(prevIndex).*field) / (real(nextIndex(k) - prevIndex(k)) * dx(k));
	}

	void update(real dt) {
		typedef typename CellGrid::iterator CellGridIter;
		CellGridIter iter;
		
		//first compute and store aux values that will be subsequently used for partial differentiation
		//during this process read and write to the same cell
	
		//beta_l, gamma_uu
		for (iter = writeCells.begin(); iter != writeCells.end(); ++iter) {
			Cell &cell = readCells(iter.index);
		
			//beta_u(k) := beta^k
			const vec &beta_u = cell.beta_u;
			
			//gamma_ll(i,j) := gamma_ij
			const symmat &gamma_ll = cell.gamma_ll;

			//beta_l(i) := g_ij beta^j
			//exclude sum of beta^t = 0
			//not storing beta_t here, since beta_t = beta^k beta_k
			vec &beta_l = cell.beta_l;
			for (int i = 0; i < dim; ++i) {
				beta_l(i) = 0;
				for (int j = 0; j < dim; ++j) {
					beta_l(i) += gamma_ll(i,j) * beta_u(j);
				}
			}
		
			//gamma_uu(i,j) := gamma^ij = inverse of gamma_ij
			// I could write the tensor formula for matrix inverses out (see "Gravitation", exercise 5.5e)
			symmat &gamma_uu = cell.gamma_uu;
			gamma_uu = symmat::invert(gamma_ll);
		}

		//conn_ull depends on gamma_uu
		for (iter = writeCells.begin(); iter != writeCells.end(); ++iter) {
			Cell &cell = readCells(iter.index);
	
			const symmat &gamma_ll = cell.gamma_ll;
			const symmat &gamma_uu = cell.gamma_uu;

			//partial_gamma_lll[k](i,j) := partial_k gamma_ij
			symmat partial_gamma_lll[dim];
			for (int k = 0; k < dim; ++k) {
				partial_gamma_lll[k] = partialDerivative(&Cell::gamma_ll, iter.index, k);
			}

			//3D hypersurface connection coefficients
			//only need spatial coefficients (since gamma^at = 0, so conn^t_ab = 0.  see "Numerical Relativity", p.48)
			//conn_lll[i](j,k) := conn_ijk = 1/2 (partial_k gamma_ij + partial_j gamma_ik - partial_i gamma_jk)
			symmat conn_lll[dim];
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					for (int k = 0; k <= j; ++k) {
						conn_lll[i](j,k) = .5 * (partial_gamma_lll[k](i,j) + partial_gamma_lll[j](i,k) - partial_gamma_lll[i](j,k));
					}
				}
			}

			//conn_ull(i)(j,k) := conn^i_jk = gamma^il conn_ljk
			::vec<dim, symmat> &conn_ull = cell.conn_ull;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					for (int k = 0; k <= j; ++k) {
						conn_ull(i)(j,k) = 0;
						for (int l = 0; l < dim; ++l) {
							conn_ull(i)(j,k) += gamma_uu(i,l) * conn_lll[l](j,k);
						}
					}
				}
			}
		}

		for (iter = writeCells.begin(); iter != writeCells.end(); ++iter) {
			const Cell &readCell = readCells(iter.index);
			Cell &writeCell = writeCells(iter.index);
			
			const real &alpha = readCell.alpha;
			
			//K_ll(i,j) := K_ij
			const symmat &K_ll = readCell.K_ll;
		
			//beta_l(k) := beta_k
			const vec &beta_l = readCell.beta_l;

			//conn_ull(i)(j,k) = conn^i_jk
			const ::vec<dim, symmat> &conn_ull = readCell.conn_ull;

			//partial_beta_ll[j](i) := partial_j beta_i
			vec partial_beta_ll[dim];
			for (int j = 0; j < dim; ++j) {
				partial_beta_ll[j] = partialDerivative(&Cell::beta_l, iter.index, j);
			}
			
			
			//diff_beta_ll(j,i) := diff_j beta_i = partial_j beta_i - conn^k_ij beta_k
			matrix diff_beta_ll;
			for (int j = 0; j < dim; ++j) {
				for (int i = 0; i < dim; ++i) {
					diff_beta_ll(j,i) = partial_beta_ll[j](i);
					for (int k = 0; k < dim; ++k) {
						diff_beta_ll(j,i) -= conn_ull(k)(i,j) * beta_l(k);
					}
				}
			}
			
			//D_i beta_j = diff_i beta_j 
			//-- but only for lower indexes.
			//doesn't work if either is upper (you get extra terms for the timelike component we're dropping)
			//(see "Numerical Relativity", exercise 2.30)
			const matrix &D_beta_ll = diff_beta_ll;
			
			Cell partial_t;
		
			//partial_t gamma_ij = -2 alpha K_ij + D_i beta_j + D_j beta_i
			for (int j = 0; j < dim; ++j) {
				for (int i = 0; i <= j; ++i) {
					partial_t.gamma_ll(i,j) = -2 * alpha * K_ll(i,j) + D_beta_ll(i,j) + D_beta_ll(j,i);
				}
			}
			
			//partial_t K_ij = alpha (R_ij - 2 K_ik K^k_j + K K_ij) - D_i D_j alpha - 8 pi alpha (S_ij - 1/2 gamma_ij (S - rho)) + beta^k partial_k K_ij + K_ik partial_j beta^k + K_kj partial_i beta^k

		}
	}
};

int main() {
	typedef double real;
	enum { dim = 2 };
	typedef ::vec<dim,real> vec;
	typedef ::vec<dim,int> veci;
	typedef ::Simulation<dim,real> Simulation;

	real dt = .01;
	real dist = 2e+3;	//2km, schwarzschild radius of our sun
	Simulation sim(vec(-dist, -dist), vec(dist, dist), veci(10, 10));
	sim.update(dt);
}

