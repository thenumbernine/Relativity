#pragma once

#include <assert.h>
#include <math.h>

#include <iostream>

#include "vector.h"
#include "tensor.h"

#include "cell.h"
#include "grid.h"

#include "inverse.h"
#include "derivative.h"

/*
ADM Formalism class generates partial t values for the K_ll and gamma_ll structures.
the rest of the properties can be seen as read-only ... should I store them in a separate structure than the Cell?
*/
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
	typedef ::vector<int,dim> DerefType;	

	typedef ::lower<dim> lower;
	typedef ::upper<dim> upper;

	typedef typename Cell::tensor_u tensor_u;
	typedef typename Cell::tensor_l tensor_l;
	typedef typename Cell::tensor_ll tensor_ll;
	typedef typename Cell::tensor_su tensor_su;
	typedef typename Cell::tensor_sl tensor_sl;
	typedef typename Cell::tensor_usl tensor_usl;
    typedef typename Cell::tensor_lsl tensor_lsl;
	typedef ::tensor<real,upper,lower> tensor_ul;
	typedef ::tensor<real,lower,upper> tensor_lu;
	typedef ::tensor<real,lower,symmetric<upper,upper>> tensor_lsu;
	typedef ::tensor<real,symmetric<lower,lower>,symmetric<lower,lower>> tensor_slsl;
	typedef ::vector<real,dim> vector;

	Grid cellHistory0, cellHistory1;
	Grid *readCells, *writeCells;

	Grid partialTCells;
	
	//resolution of our grids, stored here as well as in each grid for convenience
	DerefType size;

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

	ADMFormalism(const vector &min_, const vector &max_, const DerefType &size_)
	:	time(0),
		size(size_),
		min(min_),
		max(max_),
		range(max_ - min_),
		dx((max_ - min_) / vector(size_)),
		//need non-void constructor call, but want an array ...
		cellHistory0(size_),
		cellHistory1(size_),
		partialTCells(size_),
		readCells(&cellHistory0),
		writeCells(&cellHistory1)
	{}

	//aux calculations

	void calcOverGrid(Grid *cells, void (Cell::*updateMethod)()) {
		for (GridIter iter = cells->begin(); iter != cells->end(); ++iter) {
			Cell &cell = *iter;
			(cell.*updateMethod)();
		}
	}

	vector coordForIndex(const DerefType &index) const {
		return min + ((vector)index + .5) * dx;
	}

	//iteration
	void getGeometridynamicPartials(real dt, Grid &targetPartialTCells) {
		GridIter iter;
		
		//first compute and store aux values that will be subsequently used for partial differentiation
		//during this process read and write to the same cell
		
		//ln_psi := ln(psi) = 1/6 ln(sqrt(gamma))
		//psi = exp(ln(psi))
		calcOverGrid(readCells, &Cell::calc_psi_from_ln_sqrt_gamma);

#if 0
		//option #1: the original ADM
		//	calculate gamma^ij from gamma_ij
		//	calculate gamma from det(gamma_ij)
		{
			calcOverGrid(readCells, &Cell::calc_gamma_uu_from_gamma_ll);
		}
#endif

#if 1
		//option #2: introduction of conformal factor
		//	calculate gammaBar^ij and gammaBar_ij from psi and ln(sqrt(gamma))
		//	calculate gamma^ij from gammaBar^ij from gammaBar_ij
		//	calculate gamma from ln(sqrt(gamma))
		{
			//gammaBar_ij = psi^-4 gamma_ij
			//gammaBar^ij = inverse(gammaBar_ij)
			calcOverGrid(readCells, &Cell::calc_gammaBar_uu_and_gammaBar_ll_from_psi);
		
			//gamma^ij = psi^-4 gammaBar^ij
			calcOverGrid(readCells, &Cell::calc_gamma_uu_from_gammaBar_uu_and_psi);
		}
#endif

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
								const DerefType &index = iter.index;
								DerefType nextIndex(index);
								DerefType prevIndex(index);
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

		//K^i_j := gamma^ik K_kj
		calcOverGrid(readCells, &Cell::calc_K_ul);
		
		//K^ij = K^i_k gamma^kj
		calcOverGrid(readCells, &Cell::calc_K_uu);
		
		//tr_K_sq := tr(K^2) = K^ij K_ij
		calcOverGrid(readCells, &Cell::calc_tr_K_sq);

		//calcConstraints
		for (iter = readCells->begin(); iter != readCells->end(); ++iter) {
			Cell &cell = *iter;
		
			const tensor_su &gamma_uu = cell.gamma_uu;
			const tensor_sl &R_ll = cell.R_ll;
			const tensor_ul &K_ul = cell.K_ul;
			const tensor_su &K_uu = cell.K_uu;
			const real &K = cell.K;
			const real &tr_K_sq = cell.tr_K_sq;
			const real &rho = cell.rho;
			const tensor_u &S_u = cell.S_u;

			//R = gamma^ij R_ij
			real R = 0.;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					R += gamma_uu(i,j) * R_ll(i,j);
				}
			}

			//H = 1/2 (R + K^2 - K^j_i K^i_j) - 8 pi rho
			real &H = cell.H;
			H = .5 * (R + K * K - tr_K_sq) - 8. * M_PI * rho;

			//diff_K_uu(k,i,j) := diff_k K^ij
			//I'm not seeing a distinction between the 4D and 3D representations of the constraints.
			//The "Numerical Relativity" book did transition from a covariant form of the 4D constraint to a contravariant form of the 3D constraint
			//It also added caveats on how D_i v_j = delta_i v_j only for rank-(0,2) forms and only if v was purely spatial.
			tensor_lsu diff_K_luu = covariantDerivative(*readCells, &Cell::K_uu, dx, iter.index);

			//partial_K_l(i) := partial_i K
			tensor_l partial_K_l = partialDerivative(*readCells, &Cell::K, dx, iter.index);

			//M_u(i) := M^i = D_j (K^ij - gamma^ij K) - 8 pi S^i
			tensor_u &M_u = cell.M_u;
			for (int i = 0; i < dim; ++i) {
				M_u(i) = -8. * M_PI * S_u(i);
				for (int j = 0; j < dim; ++j) {
					M_u(i) += diff_K_luu(j,i,j) - gamma_uu(i,j) * partial_K_l(j);
				}
			}
		}

		//calcPartial
		for (iter = readCells->begin(); iter != readCells->end(); ++iter) {
			const Cell &readCell = *iter;
			
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
			const tensor_ul &K_ul = readCell.K_ul;
			const real &K = readCell.K;
			const real &tr_K_sq = readCell.tr_K_sq;

			//D_D_alpha_ll(j,i) = diff_j diff_i alpha
			tensor_ll D_D_alpha_ll = covariantDerivative(*readCells, &Cell::D_alpha_l, dx, iter.index);

			//partial_beta_ll(j)(i) := partial_j beta_i
			tensor_ll partial_beta_ll = partialDerivative(*readCells, &Cell::beta_l, dx, iter.index);
			
			//diff_beta_ll(j,i) := diff_j beta_i = partial_j beta_i - conn^k_ij beta_k
			tensor_ll diff_beta_ll = covariantDerivative(*readCells, &Cell::beta_l, dx, iter.index);
			
			//D_beta_ll(i,j) := D_i beta_j = D_i beta_j 
			//-- but only for lower indexes.
			//doesn't work if either is upper (you get extra terms for the timelike component we're dropping)
			//(see "Numerical Relativity", exercise 2.30)
			tensor_ll &D_beta_ll = diff_beta_ll;
			
			Cell &partialTCell = targetPartialTCells(iter.index);
		
			//partial_t gamma_ij = -2 alpha K_ij + D_i beta_j + D_j beta_i
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					partialTCell.gamma_ll(i,j) = -2. * alpha * K_ll(i,j) + D_beta_ll(i,j) + D_beta_ll(j,i);
				}
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
					partialTCell.K_ll(i,j) = 0.;
					partialTCell.K_ll(i,j) += alpha * R_ll(i,j);
					for (int k = 0; k < dim; ++k) {
						partialTCell.K_ll(i,j) -= alpha * 2. * K_ll(i,k) * K_ul(k,j);
					}
					partialTCell.K_ll(i,j) += alpha * K * K_ll(i,j) 
						- D_D_alpha_ll(i,j)
						- 8. * M_PI * alpha * (
							S_ll(i,j) - .5 * gamma_ll(i,j) * (S - rho));
					for (int k = 0; k < dim; ++k) {
						partialTCell.K_ll(i,j) += beta_u(k) * partial_K_lll(k,i,j)
							+ K_ll(k,i) * partial_beta_lu(j,k)
							+ K_ll(k,j) * partial_beta_lu(i,k);
					}
				}
			}

			//D_beta_lu(i,j) := D_i beta^j = D_i beta^j
			//why am I suspicious that, despite the above comment, I am doing just what the book says not to do?
			tensor_lu D_beta_lu = covariantDerivative(*readCells, &Cell::beta_u, dx, iter.index);

			//trace_D_beta := D_i beta^i
			real trace_D_beta = 0;
			for (int i = 0; i < dim; ++i) {
				trace_D_beta += D_beta_lu(i,i);
			}

			//partial_t ln_sqrt_gamma = -alpha K + D_i beta^i
			partialTCell.ln_sqrt_gamma = -alpha * K + trace_D_beta;

			//D_K_l(i) := D_i K
			tensor_l D_K_l = covariantDerivative(*readCells, &Cell::K, dx, iter.index);

			//partial_t K = -gamma^ij D_i D_j alpha + alpha(K_ij K^ij + 4 pi (rho + S)) + beta^i D_i K
			partialTCell.K = 0;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					partialTCell.K += -gamma_uu(i,j) * D_D_alpha_ll(i,j);
				}
				partialTCell.K += beta_u(i) * D_K_l(i);
			}
			partialTCell.K += alpha * (tr_K_sq + 4. * M_PI * (rho + S));
		}
	}

	void update(real dt) {
		time += dt;

		//compute aux terms
		//update alpha, beta, gamma
		getGeometridynamicPartials(dt, partialTCells);

		//update write buffer

		for (typename Grid::iterator iter = readCells->begin(); iter != readCells->end(); ++iter) {
			const Cell &readCell = *iter;
			Cell &writeCell = (*writeCells)(iter.index);
			Cell &partialTCell = partialTCells(iter.index);
			
			//now that we have the partials, well,
			//I should return them as a grid of their own ... or calculate to them as a Grid member variable
			//then accumulate them and do whatever using whatever explicit method
			//next comes implicit methods.
			writeCell = readCell + partialTCell * dt;
		}

		//TODO update fluid components

		//and swap source and dest grids
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
			o << "beta_" << coordNames[i] << "\t";
		}

		o << "ln sqrt gamma\t";

		o << "ln psi\t";

		o << "psi\t";

		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				o << "gamma_" << coordNames[i] << coordNames[j] << "\t";
			}
		}

		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				o << "gamma^" << coordNames[i] << coordNames[j] << "\t";
			}
		}

		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				o << "gammaBar_" << coordNames[i] << coordNames[j] << "\t";
			}
		}

		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				o << "gammaBar^" << coordNames[i] << coordNames[j] << "\t";
			}
		}

		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				o << "K_" << coordNames[i] << coordNames[j] << "\t";
			}
		}
		
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				o << "gamma^" << coordNames[i] << coordNames[j] << "\t";
			}
		}

		o << "psi\t";
		
		o << "K\t";

		o << "H\t";

		for (int i = 0; i < dim; ++i) {
			o << "M^" << coordNames[i] << "\t";
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
			//do coordinates move with the shift functions / lie dragging?
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

			//beta_i
			for (int i = 0; i < dim; ++i) {
				o << cell.beta_l(i) << "\t";
			}

			//ln(sqrt(gamma))
			o << cell.ln_sqrt_gamma << "\t";

			//ln(psi)
			o << cell.ln_psi << "\t";

			//psi
			o << cell.psi << "\t";

			//gamma_ij
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					o << cell.gamma_ll(i,j) << "\t";
				}
			}

			//gamma^ij
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					o << cell.gamma_uu(i,j) << "\t";
				}
			}
			
			//gammaBar_ij
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					o << cell.gammaBar_ll(i,j) << "\t";
				}
			}

			//gammaBar^ij
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					o << cell.gammaBar_uu(i,j) << "\t";
				}
			}
			
			//K_ij
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					o << cell.K_ll(i,j) << "\t";
				}
			}
			
			//gamma^ij
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					o << cell.gamma_uu(i,j) << "\t";
				}
			}
	
			//R_ij

			//R

			//rho

			//S_ij

			//psi
			o << cell.psi << "\t";
			
			//K
			o << cell.K << "\t";

			//H
			o << cell.H << "\t";

			//M^i
			for (int i = 0; i < dim; ++i) {
				o << cell.M_u(i) << "\t";
			}

			o << std::endl;
		}
	}
};

