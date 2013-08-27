#pragma once

#include <assert.h>
#include <math.h>

#include <iostream>

#include "vector.h"
#include "tensor.h"

#include "cell.h"
#include "grid.h"

#include "invert.h"
#include "derivative.h"

template<typename real_, int dim_>
struct ADMFormalism {
	typedef real_ real;
	enum { dim = dim_ };

		//cell types

	typedef ::Cell<real,dim> Cell;
	typedef ::Grid<Cell,dim> Grid;
	typedef typename Grid::iterator GridIter;

		//statically-sized mathematical types
		//(pull from cell where you can)
	
	//used for indexes
	typedef ::vector<int,dim> deref_type;	
	
	typedef typename Cell::tensor_u tensor_u;
	typedef typename Cell::tensor_l tensor_l;
	typedef typename Cell::tensor_ll tensor_ll;
	typedef typename Cell::tensor_su tensor_su;
	typedef typename Cell::tensor_sl tensor_sl;
	typedef typename Cell::tensor_usl tensor_usl;
    typedef typename Cell::tensor_lsl tensor_lsl;
	typedef ::tensor<real,upper<dim>,lower<dim>> tensor_ul;
	typedef ::tensor<real,lower<dim>,upper<dim>> tensor_lu;
	typedef ::tensor<real,symmetric<lower<dim>,lower<dim>>,symmetric<lower<dim>,lower<dim>>> tensor_slsl;
	typedef ::vector<real,dim> vector;

	Grid cellHistory0, cellHistory1;
	Grid *readCells, *writeCells;

	//resolution of our grids, stored here as well as in each grid for convenience
	deref_type size;

	//x(0) == min
	vector min;
	
	//x(size-1) == max
	vector max;
	
	//range = max - min
	vector range;
	
	//dx = range / size
	vector dx;

	//what a relative notion...
	//intended for use with output matching up iteration slices 
	real time;

	ADMFormalism(const vector &min_, const vector &max_, const deref_type &size_)
	:	time(0),
		size(size_),
		min(min_),
		max(max_),
		range(max_ - min_),
		dx((max_ - min_) / vector(size_)),
		//need non-void constructor call, but want an array ...
		cellHistory0(size_),
		cellHistory1(size_),
		readCells(&cellHistory0),
		writeCells(&cellHistory1)
	{}

	vector coordForIndex(const deref_type &index) const {
		return min + ((vector)index + .5) * dx;
	}

	//iteration
	void update(real dt) {
		time += dt;

		GridIter iter;
		
		//first compute and store aux values that will be subsequently used for partial differentiation
		//during this process read and write to the same cell
	
		//beta_l, gamma_uu, D_alpha_l
		for (iter = readCells->begin(); iter != readCells->end(); ++iter) {
			Cell &cell = *iter;
			
			const tensor_u &beta_u = cell.beta_u;
			const tensor_sl &gamma_ll = cell.gamma_ll;
			
			//D_alpha_l(i) = diff_i alpha = partial_i alpha
			tensor_l &D_alpha_l = cell.D_alpha_l;
			D_alpha_l = partialDerivative(*readCells, &Cell::alpha, dx, iter.index);

			//beta_l(i) := g_ij beta^j
			//exclude sum of beta^t = 0
			//not storing beta_t here, since beta_t = beta^k beta_k
			tensor_l &beta_l = cell.beta_l;
			for (int i = 0; i < dim; ++i) {
				beta_l(i) = 0;
				for (int j = 0; j < dim; ++j) {
					beta_l(i) += gamma_ll(i,j) * beta_u(j);
				}
			}
		
			//gamma_uu(i,j) := gamma^ij = inverse of gamma_ij
			// I could write the tensor formula for matrix inverses out (see "Gravitation", exercise 5.5e)
			tensor_su &gamma_uu = cell.gamma_uu;
			gamma_uu = invert(gamma_ll);
		}

		//conn_ull depends on gamma_uu
		for (iter = readCells->begin(); iter != readCells->end(); ++iter) {
			Cell &cell = *iter;
	
			const tensor_sl &gamma_ll = cell.gamma_ll;
			const tensor_su &gamma_uu = cell.gamma_uu;

			//partial_gamma_lll(k)(i,j) := partial_k gamma_ij
			tensor_lsl &partial_gamma_lll = cell.partial_gamma_lll;
			partial_gamma_lll = partialDerivative(*readCells, &Cell::gamma_ll, dx, iter.index);

			//3D hypersurface connection coefficients
			//only need spatial coefficients (since gamma^at = 0, so conn^t_ab = 0.  see "Numerical Relativity", p.48)
			//conn_lll(i)(j,k) := conn_ijk = 1/2 (partial_k gamma_ij + partial_j gamma_ik - partial_i gamma_jk)
			tensor_lsl &conn_lll = cell.conn_lll;
			for (int k = 0; k < dim; ++k) {
				for (int i = 0; i < dim; ++i) {
					for (int j = 0; j <= i; ++j) {
						conn_lll(i)(j,k) = .5 * (partial_gamma_lll(k)(i,j) + partial_gamma_lll(j)(i,k) - partial_gamma_lll(i)(j,k));
					}
				}
			}

			//conn_ull(i)(j,k) := conn^i_jk = gamma^il conn_ljk
			tensor_usl &conn_ull = cell.conn_ull;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					for (int k = 0; k <= j; ++k) {
						conn_ull(i)(j,k) = 0;
						for (int l = 0; l < dim; ++l) {
							conn_ull(i)(j,k) += gamma_uu(i,l) * conn_lll(l)(j,k);
						}
					}
				}
			}
		}

		//R_ll depends on conn_ull and partial_gamma_lll
		for (iter = readCells->begin(); iter != readCells->end(); ++iter) {
			Cell &cell = *iter;

			const tensor_su &gamma_uu = cell.gamma_uu;
			const tensor_usl &conn_ull = cell.conn_ull;
			const tensor_lsl &conn_lll = cell.conn_lll;
			
			//partial_partial_gamma_llll(k,l)(i,j) = partial_k partial_l gamma_ij
			tensor_slsl partial_partial_gamma_llll;
			for (int k = 0; k < dim; ++k) {
				tensor_lsl partial_gamma_lll_wrt_xk = partialDerivativeCoordinate(*readCells, &Cell::partial_gamma_lll, dx, iter.index, k);
				for (int l = 0; l <= k; ++l) {
					for (int i = 0; i < dim; ++i) {
						for (int j = 0; j <= i; ++j) {
							if (k == l) {
								//special case for 2nd deriv along same coordinate
								const deref_type &index = iter.index;
								deref_type nextIndex(index);
								deref_type prevIndex(index);
								nextIndex(k) = std::max(0, std::min(readCells->size(k)-1, index(k) + 1));
								prevIndex(k) = std::max(0, std::min(readCells->size(k)-1, index(k) - 1));
								partial_partial_gamma_llll(k,l,i,j) = 
									((*readCells)(nextIndex).partial_gamma_lll(l,i,j) 
									- cell.partial_gamma_lll(l,i,j) * 2. 
									+ (*readCells)(prevIndex).partial_gamma_lll(l,i,j))
										/ (dx(k) * dx(k));
							} else {
								partial_partial_gamma_llll(k,l,i,j) = partial_gamma_lll_wrt_xk(l,i,j);
							}
						}
					}
				}
			}

			//"Numerical Relativity" p.48
			//R_ll(i,j) := ricci_ij = 1/2 gamma^kl (partial_i partial_l gamma_kj + partial_k partial_j gamma_il - partial_i partial_j gamma_kl - partial_k partial_l gamma_ij) + gamma^kl (conn^m_il conn_mkj - conn^m_ij conn_mkl)
			tensor_sl &R_ll = cell.R_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					R_ll(i,j) = 0.;
					for (int k = 0; k < dim; ++k) {
						for (int l = 0; l < dim; ++l) {
							R_ll(i,j) += .5 * gamma_uu(k,l) * (
								partial_partial_gamma_llll(i,l,k,j)
								+ partial_partial_gamma_llll(k,j,i,l)
								- partial_partial_gamma_llll(i,j,k,l)
								- partial_partial_gamma_llll(k,l,i,j));
							for (int m = 0; m < dim; ++m) {
								R_ll(i,j) += gamma_uu(k,l) * (
									conn_ull(m,i,l) * conn_lll(m,k,j)
									- conn_ull(m,i,j) * conn_lll(m,k,l));
							}
						}
					}
				}
			}
		}

		for (iter = writeCells->begin(); iter != writeCells->end(); ++iter) {
			const Cell &readCell = (*readCells)(iter.index);
			Cell &writeCell = *iter;
			
			const real &alpha = readCell.alpha;
			const tensor_u &beta_u = readCell.beta_u;
			const tensor_sl &gamma_ll = readCell.gamma_ll;
			const tensor_sl &K_ll = readCell.K_ll;
			const real &rho = readCell.rho;
			const tensor_ll &S_ll = readCell.S_ll;
			const tensor_l &beta_l = readCell.beta_l;
			const tensor_usl &conn_ull = readCell.conn_ull;
			const tensor_sl &R_ll = readCell.R_ll;
			const tensor_su &gamma_uu = readCell.gamma_uu;

			//D_D_alpha_ll(j,i) = diff_j diff_i alpha
			tensor_ll D_D_alpha_ll = covariantDerivative(*readCells, &Cell::D_alpha_l, dx, iter.index);

			//partial_beta_ll(j)(i) := partial_j beta_i
			tensor_ll partial_beta_ll = partialDerivative(*readCells, &Cell::beta_l, dx, iter.index);
			
			//diff_beta_ll(j,i) := diff_j beta_i = partial_j beta_i - conn^k_ij beta_k
			tensor_ll diff_beta_ll = covariantDerivative(*readCells, &Cell::beta_l, dx, iter.index);
			
			//D_i beta_j = diff_i beta_j 
			//-- but only for lower indexes.
			//doesn't work if either is upper (you get extra terms for the timelike component we're dropping)
			//(see "Numerical Relativity", exercise 2.30)
			tensor_ll &D_beta_ll = diff_beta_ll;
			
			Cell partial_t;
		
			//partial_t gamma_ij = -2 alpha K_ij + D_i beta_j + D_j beta_i
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					partial_t.gamma_ll(i,j) = -2. * alpha * K_ll(i,j) + D_beta_ll(i,j) + D_beta_ll(j,i);
				}
			}

			//K^i_j := gamma^ik K_kj
			tensor_ul K_ul;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					K_ul(i,j) = 0;
					for (int k = 0; k < dim; ++k) {
						K_ul(i,j) += gamma_uu(i,k) * K_ll(k,j);
					}
				}
			}

			//K := K^i_i
			real K = 0.;
			for (int i = 0; i < dim; ++i) {
				K += K_ul(i,i);
			}

			//partial_K_lll(k,i,j) := partial_k K_ij
			tensor_lsl partial_K_lll = partialDerivative(*readCells, &Cell::K_ll, dx, iter.index);

			//partial_beta_lu(j,i) := partial_j beta^i
			tensor_lu partial_beta_lu = partialDerivative(*readCells, &Cell::beta_u, dx, iter.index);

			//S := S^i_i := gamma^ij S_ij
			real S = 0.;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					S += S_ll(i,j) * gamma_uu(i,j);
				}
			}
			
			//partial_t K_ij = alpha (R_ij - 2 K_ik K^k_j + K K_ij) - D_i D_j alpha - 8 pi alpha (S_ij - 1/2 gamma_ij (S - rho)) + beta^k partial_k K_ij + K_ik partial_j beta^k + K_kj partial_i beta^k
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					partial_t.K_ll(i,j) = 0.;
					partial_t.K_ll(i,j) += alpha * R_ll(i,j);
					for (int k = 0; k < dim; ++k) {
						partial_t.K_ll(i,j) -= alpha * 2. * K_ll(i,k) * K_ul(k,j);
					}
					partial_t.K_ll(i,j) += alpha * K * K_ll(i,j) 
						- D_D_alpha_ll(i,j)
						- 8. * M_PI * alpha * (
							S_ll(i,j) - .5 * gamma_ll(i,j) * (S - rho));
					for (int k = 0; k < dim; ++k) {
						partial_t.K_ll(i,j) += beta_u(k) * partial_K_lll(k,i,j)
							+ K_ll(k,i) * partial_beta_lu(j,k)
							+ K_ll(k,j) * partial_beta_lu(i,k);
					}
				}
			}

			//now that we have the partials, well,
			//I should return them as a grid of their own ... or calculate to them as a Grid member variable
			//then accumulate them and do whatever using whatever explicit method
			//next comes implicit methods.
			writeCell = readCell + partial_t * dt;
		}

		//and swap
		//do something more clever if we ever get any more than 2 histories
		std::swap(readCells, writeCells);
	}

	void outputHeaders(std::ostream &o) {
		static const char *coordNames[] = {"x", "y", "z"};
		
		o << "#";

		o << "t\t";
	
		for (int i = 0; i < dim; ++i) {
			o << coordNames[i] << "\t";
		}

		o << "alpha\t";
		
		for (int i = 0; i < dim; ++i) {
			o << "beta^" << coordNames[i] << "\t";
		}

		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				o << "gamma_" << coordNames[i] << coordNames[j] << "\t";
			}
		}

		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				o << "K_" << coordNames[i] << coordNames[j] << "\t";
			}
		}

		o << std::endl;
	}

	void outputLine(std::ostream &o) {
	
		//typedef typename Grid::const_iterator ConstGridIter;
		for (GridIter iter = readCells->begin(); iter != readCells->end(); ++iter) {
			const Cell &cell = *iter;

			if (iter.index(0) == 0) o << "\n";
		
			//time
			o << time << "\t";

			//space
			//this can't be right.  don't coordinates move with the shift functions / lie dragging?
			//should I initialize these and iterate them as well?
			vector x = coordForIndex(iter.index);
			for (int i = 0; i < dim; ++i) {
				o << x(i) << "\t";
			}

			//alpha
			o << cell.alpha << "\t";
			
			//beta^i
			for (int i = 0; i < dim; ++i) {
				o << cell.beta_u(i) << "\t";
			}

			//gamma_ij
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					o << cell.gamma_ll(i,j) << "\t";
				}
			}
		
			//K_ij
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					o << cell.K_ll(i,j) << "\t";
				}
			}
		
			//K

			//R_ij

			//R

			//rho

			//S_ij
		
			o << std::endl;
		}
	}
};

