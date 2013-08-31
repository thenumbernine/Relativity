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
template<typename real_, int dim_, typename Integrator_>
struct ADMFormalism {
	typedef real_ real;
	enum { dim = dim_ };
	
		//integrator types
	
	typedef Integrator_ Integrator;
	typedef typename Integrator::template Body<ADMFormalism> IntegratorBody;

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

		//integrator
	
	IntegratorBody integrator;

		//grids

	Grid cellHistory0, cellHistory1;
	Grid *readCells, *writeCells;
	
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
	:	integrator(size_),
		time(0),
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
	{
		integrator.sim = this;
	}

	//aux calculations

	void calcOverGrid(Grid &cells, void (Cell::*updateMethod)()) {
		for (GridIter iter = cells.begin(); iter != cells.end(); ++iter) {
			Cell &cell = *iter;
			(cell.*updateMethod)();
		}
	}

	//calculates partial_gammaField_lll, connField_lll, connField_ull, RField_ll, RField
	// depends on gammaField_ll, gammaField_uu
	void calcConnections(
		Grid &cells,
		const tensor_sl Cell::*gammaField_ll, 
		const tensor_su Cell::*gammaField_uu,
		tensor_lsl Cell::*partial_gammaField_lll,
		tensor_lsl Cell::*connField_lll,
		tensor_usl Cell::*connField_ull,
		tensor_sl Cell::*RField_ll,
		real Cell::*RField)
	{
		GridIter iter;	
		
		//partial_gamma_lll depends on gamma_ll
		//conn_lll depends on partial_gamma_lll
		//conn_ull depends on conn_lll and gamma_uu
		for (iter = cells.begin(); iter != cells.end(); ++iter) {
			Cell &cell = *iter;
			
			const tensor_sl &gamma_ll = cell.*gammaField_ll;
			const tensor_su &gamma_uu = cell.*gammaField_uu;
		
			//partial_gamma_lll(k,i,j) := partial_k gamma_ij
			tensor_lsl &partial_gamma_lll = cell.*partial_gammaField_lll;
			partial_gamma_lll = partialDerivative(cells, gammaField_ll, dx, iter.index);

			//conn_lll(i,j,k) := conn_ijk = 1/2 (partial_k gamma_ij + partial_j gamma_ik - partial_i gamma_jk)
			tensor_lsl &conn_lll = cell.*connField_lll;
			for (int k = 0; k < dim; ++k) {
				for (int i = 0; i < dim; ++i) {
					for (int j = 0; j <= i; ++j) {
						conn_lll(i,j,k) = .5 * (partial_gamma_lll(k,i,j) + partial_gamma_lll(j,i,k) - partial_gamma_lll(i,j,k));
					}
				}
			}

			//conn_ull(i,j,k) := conn^i_jk = gamma^il conn_ljk
			tensor_usl &conn_ull = cell.*connField_ull;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					for (int k = 0; k <= j; ++k) {
						conn_ull(i,j,k) = 0;
						for (int l = 0; l < dim; ++l) {
							conn_ull(i,j,k) += gamma_uu(i,l) * conn_lll(l,j,k);
						}
					}
				}
			}
		}

		//R_ll depends on conn_ull and partial_gamma_lll
		for (iter = cells.begin(); iter != cells.end(); ++iter) {
			Cell &cell = *iter;

			const tensor_su &gamma_uu = cell.*gammaField_uu;
			const tensor_usl &conn_ull = cell.*connField_ull;
			const tensor_lsl &conn_lll = cell.*connField_lll;

			//partial2_gamma_llll(k,l,i,j) = partial_k partial_l gamma_ij
			tensor_slsl partial2_gamma_llll = partialSecondDerivative(cells, gammaField_ll, partial_gammaField_lll, dx, iter.index);

			//"Numerical Relativity" p.48
			//R_ll(i,j) := R_ij = 1/2 gamma^kl (partial_i partial_l gamma_kj + partial_k partial_j gamma_il - partial_i partial_j gamma_kl - partial_k partial_l gamma_ij) + gamma^kl (conn^m_il conn_mkj - conn^m_ij conn_mkl)
			tensor_sl &R_ll = cell.*RField_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					R_ll(i,j) = 0.;
					for (int k = 0; k < dim; ++k) {
						for (int l = 0; l < dim; ++l) {
							R_ll(i,j) += .5 * gamma_uu(k,l) * (
								partial2_gamma_llll(i,l,k,j)
								+ partial2_gamma_llll(k,j,i,l)
								- partial2_gamma_llll(i,j,k,l)
								- partial2_gamma_llll(k,l,i,j));
							for (int m = 0; m < dim; ++m) {
								R_ll(i,j) += gamma_uu(k,l) * (
									conn_ull(m,i,l) * conn_lll(m,k,j)
									- conn_ull(m,i,j) * conn_lll(m,k,l));
							}
						}
					}
				}
			}
			
			//R = gamma^ij R_ij
			real &R = cell.*RField;
			R = 0;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					R += gamma_uu(i,j) * R_ll(i,j);
				}
			}
		}
	}

	vector coordForIndex(const DerefType &index) const {
		return min + ((vector)index + .5) * dx;
	}

	//iteration
	void getGeometridynamicPartials(real dt, Grid &targetReadCells, Grid &targetPartialTCells) {
		GridIter iter;
		
		//first compute and store aux values that will be subsequently used for partial differentiation
		//during this process read and write to the same cell
		
		//ln_psi := ln(psi) = 1/6 ln(sqrt(gamma))
		//psi = exp(ln(psi))
		calcOverGrid(targetReadCells, &Cell::calc_psi_from_ln_sqrt_gamma);

#if 0
		//option #1: the original ADM implementation
		//	calculate gamma^ij from gamma_ij
		//	calculate gamma from det(gamma_ij)
		calcOverGrid(targetReadCells, &Cell::calc_gamma_uu_from_gamma_ll);
#endif
#if 1
		//option #2: introduction of conformal factor
		//	calculate gammaBar^ij and gammaBar_ij from psi and ln(sqrt(gamma))
		//	calculate gamma^ij from gammaBar^ij from gammaBar_ij
		//	calculate gamma from ln(sqrt(gamma))
			
		//gammaBar_ij = psi^-4 gamma_ij
		//gammaBar^ij = inverse(gammaBar_ij)
		calcOverGrid(targetReadCells, &Cell::calc_gammaBar_uu_and_gammaBar_ll_from_psi);
	
		//gamma^ij = psi^-4 gammaBar^ij
		calcOverGrid(targetReadCells, &Cell::calc_gamma_uu_from_gammaBar_uu_and_psi);
#endif
		
		//beta_l, gamma_uu, D_alpha_l
		for (iter = targetReadCells.begin(); iter != targetReadCells.end(); ++iter) {
			Cell &cell = *iter;
			
			const tensor_u &beta_u = cell.beta_u;
			const tensor_sl &gamma_ll = cell.gamma_ll;
	
			//D_alpha_l(i) := diff_i alpha = partial_i alpha
			tensor_l &D_alpha_l = cell.D_alpha_l;
			D_alpha_l = partialDerivative(targetReadCells, &Cell::alpha, dx, iter.index);

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

		//DBar_ln_psi_l depends on ln_psi
		//DBar_psi_l depends on psi
		for (iter = targetReadCells.begin(); iter != targetReadCells.end(); ++iter) {
			Cell &cell = *iter;
	
			//DBar_ln_psi_l(i) := DBar_i ln(psi) = partial_i ln(psi)
			cell.DBar_ln_psi_l = partialDerivative(targetReadCells, &Cell::ln_psi, dx, iter.index);

			//DBar_psi_l(i) := DBar_i psi
			cell.DBar_psi_l = partialDerivative(targetReadCells, &Cell::psi, dx, iter.index);
		}

		//calculates partial_gammaBar_lll, connBar_lll, connBar_ull, RBar_ll, RBar
		// depends on gammaBar_ll, gammaBar_uu
		calcConnections(targetReadCells, &Cell::gammaBar_ll, &Cell::gammaBar_uu, &Cell::partial_gammaBar_lll, &Cell::connBar_lll, &Cell::connBar_ull, &Cell::RBar_ll, &Cell::RBar);

#if 0	//option #1: original
		//TODO if you want to use this here and the conformal-based Hamiltonian later, you have to move a few aux field calculations out of the commented block below 

		//calculates partial_gamma_lll, conn_lll, conn_ull, R_ll, R
		// depends on gamma_ll, gamma_uu
		calcConnections(targetReadCells, &Cell::gamma_ll, &Cell::gamma_uu, &Cell::partial_gamma_lll, &Cell::conn_lll, &Cell::conn_ull, &Cell::R_ll, &Cell::R);
#endif
#if 1	//use conformal metric info
	
		//conn_ull depends on connBar_ull, DBar_ln_psi_l, gammaBar_uu, gammaBar_ll
		//conn_lll depends on conn_ull and gamma_ll
		for (iter = targetReadCells.begin(); iter != targetReadCells.end(); ++iter) {
			Cell &cell = *iter;

			const tensor_l &DBar_ln_psi_l = cell.DBar_ln_psi_l;
			const tensor_sl &gammaBar_ll = cell.gammaBar_ll;
			const tensor_su &gammaBar_uu = cell.gammaBar_uu;
			const tensor_usl &connBar_ull = cell.connBar_ull;
			const tensor_sl &gamma_ll = cell.gamma_ll;
			tensor_usl &conn_ull = cell.conn_ull;
			tensor_lsl &conn_lll = cell.conn_lll;
	
			//conn^i_jk = connBar^i_jk + 2 (delta^i_j DBar_k ln(psi) + delta^i_k DBar_j ln(psi) - gammaBar_jk gammaBar^il DBar_l ln(psi))
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					for (int k = 0; k <= j; ++k) {
						conn_ull(i,j,k) = connBar_ull(i,j,k);
						if (i == j) conn_ull(i,j,k) += 2. * DBar_ln_psi_l(k);
						if (i == k) conn_ull(i,j,k) += 2. * DBar_ln_psi_l(j);
						for (int l = 0; l < dim; ++l) {
							conn_ull(i,j,k) -= 2. * gammaBar_ll(j,k) * gammaBar_uu(i,l) * DBar_ln_psi_l(l);
						}
					}
				}
			}

			//conn_ijk = gamma_il conn^l_jk
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					for (int k = 0; k <= j; ++k) {
						conn_lll(i,j,k) = 0;
						for (int l = 0; l < dim; ++l) {
							conn_lll(i,j,k) += gamma_ll(i,l) * conn_ull(l,j,k);
						}
					}
				}
			}
		}
		
		for (iter = targetReadCells.begin(); iter != targetReadCells.end(); ++iter) {
			Cell &cell = *iter;

			const real &psi = cell.psi;
			const real &RBar = cell.RBar;
			const tensor_l &DBar_ln_psi_l = cell.DBar_ln_psi_l;
			const tensor_l &DBar_psi_l = cell.DBar_psi_l;
			const tensor_su &gammaBar_uu = cell.gammaBar_uu;
			const tensor_sl &gammaBar_ll = cell.gammaBar_ll;
			const tensor_sl &RBar_ll = cell.RBar_ll;
			const tensor_usl &connBar_ull = cell.connBar_ull;

			//partial2_ln_psi_ll(i,j) := partial_i partial_j ln(psi)
			tensor_sl partial2_ln_psi_ll = partialSecondDerivative(targetReadCells, &Cell::ln_psi, &Cell::DBar_ln_psi_l, dx, iter.index);

			//DBar2_ln_psi_ll(i,j) := DBar_i DBar_j ln(psi) = partial_i partial_j ln(psi) - connBar^k_ij partial_k ln(psi)
			tensor_sl DBar2_ln_psi_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					DBar2_ln_psi_ll(i,j) = partial2_ln_psi_ll(i,j);
					for (int k = 0; k < dim; ++k) {
						DBar2_ln_psi_ll(i,j) -= connBar_ull(k,i,j) * DBar_ln_psi_l(k);
					}
				}
			}

			//DBar2_ln_psi := DBar^2 ln(psi) = gammaBar^ij DBar_i DBar_j ln(psi)
			real DBar2_ln_psi = 0;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					DBar2_ln_psi += gammaBar_uu(i,j) * DBar2_ln_psi_ll(i,j);
				}
			}

			//normBar_DBar_ln_psi = gammaBar^ij (DBar_i ln(psi)) (DBar_j ln(psi))
			real normBar_DBar_ln_psi = 0;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					normBar_DBar_ln_psi += gammaBar_uu(i,j) * DBar_ln_psi_l(i) * DBar_ln_psi_l(j);
				}
			}

			//"Numerical Relativity" p.57
			//R_ll(i,j) := R_ij = RBar_ij - 2 (DBar_i DBar_j ln(psi) + gammaBar_ij gammaBar^lm DBar_l DBar_m ln(psi)) + 4((DBar_i ln(psi)) (DBar_j ln(psi)) - gammaBar_ij gammaBar^lm (DBar_l ln(psi)) (DBar_m ln(psi)))
			tensor_sl &R_ll = cell.R_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					R_ll(i,j) = RBar_ll(i,j) - 2. * (DBar2_ln_psi_ll(i,j) + gammaBar_ll(i,j) * DBar2_ln_psi) + 4. * (DBar_ln_psi_l(i) * DBar_ln_psi_l(j) - gammaBar_ll(i,j) * normBar_DBar_ln_psi);
				}
			}

			//partial2_psi_ll(i,j) := partial_i partial_j psi
			tensor_sl partial2_psi_ll = partialSecondDerivative(targetReadCells, &Cell::psi, &Cell::DBar_psi_l, dx, iter.index);

			//DBar2_psi_ll(i,j) := DBar_i DBar_j ln(psi) = partial_i partial_j psi - connBar^k_ij partial_k psi
			tensor_sl DBar2_psi_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					DBar2_psi_ll(i,j) = partial2_psi_ll(i,j);
					for (int k = 0; k < dim; ++k) {
						DBar2_psi_ll(i,j) -= connBar_ull(k,i,j) * DBar_psi_l(k);
					}
				}
			}

			//DBar2_psi := DBar^2 psi = gammaBar^ij DBar_i DBar_j psi
			real &DBar2_psi = cell.DBar2_psi;
			DBar2_psi = 0;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					DBar2_psi += gammaBar_uu(i,j) * DBar2_psi_ll(i,j);
				}
			}

			//R = psi^-4 RBar - 8 psi^-5 DBar^2 psi = (RBar - 8 (DBar^2 psi) / psi) / psi^4
			real psiSquared = psi * psi;
			real psiToTheFourth = psiSquared * psiSquared;
			real &R = cell.R;
			R = (RBar - 8. * DBar2_psi / psi) / psiToTheFourth;
		}
#endif

		//K^i_j := gamma^ik K_kj
		calcOverGrid(targetReadCells, &Cell::calc_K_ul);
		
		//K^ij = K^i_k gamma^kj
		calcOverGrid(targetReadCells, &Cell::calc_K_uu);
		
		//tr_K_sq := tr(K^2) = K^ij K_ij
		calcOverGrid(targetReadCells, &Cell::calc_tr_K_sq);

		//calcConstraints
		for (iter = targetReadCells.begin(); iter != targetReadCells.end(); ++iter) {
			Cell &cell = *iter;
	
			const real &psi = cell.psi;
			const real &DBar2_psi = cell.DBar2_psi;
			const real &R = cell.R;
			const real &RBar = cell.RBar;
			const tensor_su &gamma_uu = cell.gamma_uu;
			const tensor_sl &R_ll = cell.R_ll;
			const tensor_ul &K_ul = cell.K_ul;
			const tensor_su &K_uu = cell.K_uu;
			const real &K = cell.K;
			const real &tr_K_sq = cell.tr_K_sq;
			const real &rho = cell.rho;
			const tensor_u &S_u = cell.S_u;

			//Hamiltonian constraint
#if 0		//option #1: use original ADM method

			//H = 1/2 (R + K^2 - K^j_i K^i_j) - 8 pi rho
			real &H = cell.H;
			H = .5 * (R + K * K - tr_K_sq) - 8. * M_PI * rho;
#else		//option #2: use conformal values
	
			real psiSquared = psi * psi;
			real psiToTheFourth = psiSquared * psiSquared;
			
			//"Numerical Relativity" p.57
			//H = 1/2 ( -8 DBar^2 psi + psi RBar + psi^5 (K^2 - tr(K^2)) - 16 pi psi^5 rho)
			//  = 1/2 ( -8 DBar^2 psi + psi (RBar + psi^4 (K^2 - tr(K^2) - 16 pi rho)))
			real &H = cell.H;
			H = .5 * (-8. * DBar2_psi + psi * (RBar * psiToTheFourth * (K * K - tr_K_sq - 16. * M_PI * rho)));
#endif

			//diff_K_uu(k,i,j) := diff_k K^ij
			//I'm not seeing a distinction between the 4D and 3D representations of the constraints.
			//The "Numerical Relativity" book did transition from a covariant form of the 4D constraint to a contravariant form of the 3D constraint
			//It also added caveats on how D_i v_j = delta_i v_j only for rank-(0,2) forms and only if v was purely spatial.
			tensor_lsu diff_K_luu = covariantDerivative(targetReadCells, &Cell::K_uu, dx, iter.index, &Cell::conn_ull);

			//partial_K_l(i) := partial_i K
			tensor_l partial_K_l = partialDerivative(targetReadCells, &Cell::K, dx, iter.index);

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
		for (iter = targetReadCells.begin(); iter != targetReadCells.end(); ++iter) {
			const Cell &readCell = *iter;
			
			const real &alpha = readCell.alpha;
			const tensor_l &D_alpha_l = readCell.D_alpha_l;
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

			//partial2_alpha_ll(i,j) := partial_i partial_j alpha
			tensor_sl partial2_alpha_ll = partialSecondDerivative(targetReadCells, &Cell::alpha, &Cell::D_alpha_l, dx, iter.index);

			//D2_alpha_ll(i,j) = D_i D_j alpha = partial_i partial_j alpha - conn^k_ij partial_k alpha
			tensor_sl D2_alpha_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					D2_alpha_ll(i,j) = partial2_alpha_ll(i,j);
					for (int k = 0; k < dim; ++k) {
						D2_alpha_ll(i,j) -= conn_ull(k,i,j) * D_alpha_l(k);
					}
				}
			}

			//partial_beta_ll(j,i) := partial_j beta_i
			tensor_ll partial_beta_ll = partialDerivative(targetReadCells, &Cell::beta_l, dx, iter.index);
			
			//D_beta_ll(j,i) := D_j beta_i = partial_j beta_i - conn^k_ij beta_k
			tensor_ll D_beta_ll = covariantDerivative(targetReadCells, &Cell::beta_l, dx, iter.index, &Cell::conn_ull);
			
			//diff_i beta_j = D_i beta_j 
			//-- but only for lower indexes.
			//doesn't work if either is upper (you get extra terms for the timelike component we're dropping)
			//(see "Numerical Relativity", exercise 2.30)
			
			Cell &partialTCell = targetPartialTCells(iter.index);
		
			//partial_t gamma_ij = -2 alpha K_ij + D_i beta_j + D_j beta_i
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					partialTCell.gamma_ll(i,j) = -2. * alpha * K_ll(i,j) + D_beta_ll(i,j) + D_beta_ll(j,i);
				}
			}

			//partial_K_lll(k,i,j) := partial_k K_ij
			tensor_lsl partial_K_lll = partialDerivative(targetReadCells, &Cell::K_ll, dx, iter.index);

			//partial_beta_lu(j,i) := partial_j beta^i
			tensor_lu partial_beta_lu = partialDerivative(targetReadCells, &Cell::beta_u, dx, iter.index);

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
						- D2_alpha_ll(i,j)
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
			tensor_lu D_beta_lu = covariantDerivative(targetReadCells, &Cell::beta_u, dx, iter.index, &Cell::conn_ull);

			//trace_D_beta := D_i beta^i
			real trace_D_beta = 0;
			for (int i = 0; i < dim; ++i) {
				trace_D_beta += D_beta_lu(i,i);
			}

			//partial_t ln_sqrt_gamma = -alpha K + D_i beta^i
			partialTCell.ln_sqrt_gamma = -alpha * K + trace_D_beta;

			//D_K_l(i) := D_i K
			tensor_l D_K_l = covariantDerivative(targetReadCells, &Cell::K, dx, iter.index, &Cell::conn_ull);

			//partial_t K = -gamma^ij D_i D_j alpha + alpha(K_ij K^ij + 4 pi (rho + S)) + beta^i D_i K
			partialTCell.K = 0;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					partialTCell.K += -gamma_uu(i,j) * D2_alpha_ll(i,j);
				}
				partialTCell.K += beta_u(i) * D_K_l(i);
			}
			partialTCell.K += alpha * (tr_K_sq + 4. * M_PI * (rho + S));
		}
	}

	void update(real dt) {
		integrator.update(dt);

		time += dt;

		//TODO update fluid components

		//and swap source and dest grids
		//do something more clever if we ever get any more than 2 histories
		//if you never need more than 2 then maybe we won't need histories at all, just partials?
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

		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				o << "R_" << coordNames[i] << coordNames[j] << "\t";
			}
		}

		o << "R\t";

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
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					o << cell.R_ll(i,j) << "\t";
				}
			}

			//R
			o << cell.R << "\t";

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

