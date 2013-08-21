#include <assert.h>

#include "vector.h"
#include "oneform.h"
#include "symmat.h"
#include "matrix.h"
#include "cell.h"
#include "grid.h"

template<int dim_, typename real_>
struct Simulation {
	enum { dim = dim_ };
	typedef real_ real;

	//statically-sized mathematical types 
	typedef ::vector<dim,int> veci;		//used for indexes
	typedef ::vector<dim,real> vector;		//3D contravariant / rank-(1,0) tensors
	typedef ::oneform<dim,real> oneform;	//3D covariant / rank-(0,1) tensors
	typedef ::vector<dim+1,real> tvec;		//4D hypersurface + timelike component vector
	typedef ::matrix<dim,real> matrix;		//3D matrix / rank-2 tensor of some sort
	typedef ::symmat<dim,real> symmat;		//3D symmetric matrix / rank-2 tensor of some sort

	typedef ::Cell<dim,real> Cell;
	typedef ::Grid<dim,Cell> CellGrid;
	CellGrid readCells, writeCells;

	//kronecher delta for integer vector / indexes, such that dxi(i)(j) == i == j
	::vector<dim,veci> dxi;

	//resolution of our grids, stored here as well as in each grid for convenience
	veci size;

	//x(0) == min
	vector min;
	
	//x(size-1) == max
	vector max;
	
	//range = max - min
	vector range;
	
	//dx = range / size
	vector dx;

	Simulation(const vector &min_, const vector &max_, const veci &size_)
	:	size(size_),
		min(min_),
		max(max_),
		range(max_ - min_),
		dx((max_ - min_) / vector(size_)),
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
	::oneform<dim,T> partialDerivative(T Cell::*field, const veci &index) {
		veci nextIndex = veci::clamp(index + dxi(dim), veci(0), size-1);
		veci prevIndex = veci::clamp(index - dxi(dim), veci(0), size-1);
		::oneform<dim,T> result;
		for (int k = 0; k < dim; ++k) {
			result(k) = (readCells(nextIndex).*field - readCells(prevIndex).*field) / (real(nextIndex(k) - prevIndex(k)) * dx(k));
		}
		return result;
	}

	/*
	covariant derivative
	depends on conn_ull and partialDerivative
		up = the up/down of each index (as 1 for up or 0 for down)
		the up/down information could be stored per-type.
		just need a means to enumerator across it ...
	until then, need specializations for rank-1 and rank-2
	
	here's the rank-(1,0):
		diff_k x^i = partial_k x^i + conn^i_jk x^j
	*/
	::oneform<dim,vector> covariantDerivative(vector Cell::*field, const veci &index) {
		::oneform<dim,vector> result = partialDerivative(field, index);
		const Cell &cell = readCells(index);
		for (int k = 0; k < dim; ++k) {
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					result(k)(i) += cell.conn_ull(i)(j,k) * (cell.*field)(j);
				}
			}
		}
		return result;
	}

	/*
	here's the rank-(0,1):
		diff_k x_i = partial_k x_i - conn^j_ik x_j
	*/
	::oneform<dim,oneform> covariantDerivative(oneform Cell::*field, const veci &index) {
		::oneform<dim,oneform> result = partialDerivative(field, index);
		const Cell &cell = readCells(index);
		for (int k = 0; k < dim; ++k) {
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					result(k)(i) -= cell.conn_ull(j)(i,k) * (cell.*field)(j);
				}
			}
		}
		return result;
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
			const vector &beta_u = cell.beta_u;
			
			//gamma_ll(i,j) := gamma_ij
			const symmat &gamma_ll = cell.gamma_ll;

			//beta_l(i) := g_ij beta^j
			//exclude sum of beta^t = 0
			//not storing beta_t here, since beta_t = beta^k beta_k
			oneform &beta_l = cell.beta_l;
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

			//partial_gamma_lll(k)(i,j) := partial_k gamma_ij
			::oneform<dim, symmat> partial_gamma_lll = partialDerivative(&Cell::gamma_ll, iter.index);

			//3D hypersurface connection coefficients
			//only need spatial coefficients (since gamma^at = 0, so conn^t_ab = 0.  see "Numerical Relativity", p.48)
			//conn_lll[i](j,k) := conn_ijk = 1/2 (partial_k gamma_ij + partial_j gamma_ik - partial_i gamma_jk)
			symmat conn_lll[dim];
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					for (int k = 0; k <= j; ++k) {
						conn_lll[i](j,k) = .5 * (partial_gamma_lll(k)(i,j) + partial_gamma_lll(j)(i,k) - partial_gamma_lll(i)(j,k));
					}
				}
			}

			//conn_ull(i)(j,k) := conn^i_jk = gamma^il conn_ljk
			::vector<dim, symmat> &conn_ull = cell.conn_ull;
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
			const oneform &beta_l = readCell.beta_l;

			//conn_ull(i)(j,k) = conn^i_jk
			const ::vector<dim, symmat> &conn_ull = readCell.conn_ull;

			//partial_beta_ll(j)(i) := partial_j beta_i
			::oneform<dim,oneform> partial_beta_ll = partialDerivative(&Cell::beta_l, iter.index);
			
			//diff_beta_ll(j,i) := diff_j beta_i = partial_j beta_i - conn^k_ij beta_k
			::oneform<dim,oneform> diff_beta_ll = covariantDerivative(&Cell::beta_l, iter.index);
			
			//D_i beta_j = diff_i beta_j 
			//-- but only for lower indexes.
			//doesn't work if either is upper (you get extra terms for the timelike component we're dropping)
			//(see "Numerical Relativity", exercise 2.30)
			const ::oneform<dim,oneform> &D_beta_ll = diff_beta_ll;
			
			Cell partial_t;
		
			//partial_t gamma_ij = -2 alpha K_ij + D_i beta_j + D_j beta_i
			for (int j = 0; j < dim; ++j) {
				for (int i = 0; i <= j; ++i) {
					partial_t.gamma_ll(i,j) = -2 * alpha * K_ll(i,j) + D_beta_ll(i)(j) + D_beta_ll(j)(i);
				}
			}
			
			//partial_t K_ij = alpha (R_ij - 2 K_ik K^k_j + K K_ij) - D_i D_j alpha - 8 pi alpha (S_ij - 1/2 gamma_ij (S - rho)) + beta^k partial_k K_ij + K_ik partial_j beta^k + K_kj partial_i beta^k

		}
	}
};

int main() {
	typedef double real;
	enum { dim = 2 };
	typedef ::vector<dim,real> vector;
	typedef ::vector<dim,int> veci;
	typedef ::Simulation<dim,real> Simulation;

	real dt = .01;
	real dist = 2e+3;	//2km, schwarzschild radius of our sun
	Simulation sim(vector(-dist, -dist), vector(dist, dist), veci(10, 10));
	sim.update(dt);
}

