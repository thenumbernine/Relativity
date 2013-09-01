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

#include "i_integrator.h"
#include "i_admformalism.h"

/*
ADM Formalism class generates partial t values for the K_ll and gamma_ll structures.
the rest of the properties can be seen as read-only ... should I store them in a separate structure than the AuxCell?
*/
template<typename real_, int dim_>
struct ADMFormalism : public IADMFormalism<real_, dim_> {
	typedef real_ real;
	enum { dim = dim_ };
	
		//cell types

	typedef ::GeomCell<real,dim> GeomCell;
	typedef ::MatterCell<real,dim> MatterCell;
	typedef ::AuxCell<real,dim> AuxCell;
	
	typedef ::Grid<GeomCell,dim> GeomGrid;
	typedef ::Grid<MatterCell,dim> MatterGrid;
	typedef ::Grid<AuxCell,dim> AuxGrid;

		//statically-sized mathematical types
		//(pull from cell where you can)
	
	//used for indexes
	typedef ::vector<int,dim> DerefType;	

	typedef ::lower<dim> lower;
	typedef ::upper<dim> upper;

	typedef typename AuxCell::tensor_u tensor_u;
	typedef typename AuxCell::tensor_l tensor_l;
	typedef typename AuxCell::tensor_ll tensor_ll;
	typedef typename AuxCell::tensor_su tensor_su;
	typedef typename AuxCell::tensor_sl tensor_sl;
	typedef typename AuxCell::tensor_usl tensor_usl;
    typedef typename AuxCell::tensor_lsl tensor_lsl;
	typedef ::tensor<real,upper,lower> tensor_ul;
	typedef ::tensor<real,lower,upper> tensor_lu;
	typedef ::tensor<real,lower,symmetric<upper,upper>> tensor_lsu;
	typedef ::tensor<real,symmetric<lower,lower>,symmetric<lower,lower>> tensor_slsl;
	typedef ::vector<real,dim> vector;
	
	//what a relative notion...
	//intended for use with output matching up iteration slices 
	real time;

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

	//integrator
	IIntegrator<real, dim> *integrator;

	//grids
	GeomGrid geomGrid0;
	GeomGrid geomGrid1;
	MatterGrid matterGrid;
	AuxGrid auxGrid;
	
	GeomGrid *geomGridReadCurrent;
	GeomGrid *geomGridWriteCurrent;
	
	ADMFormalism(
		const vector &min_, 
		const vector &max_, 
		const DerefType &size_,
		IIntegrator<real, dim> *integrator_)
	:	time(0),
		size(size_),
		min(min_),
		max(max_),
		range(max_ - min_),
		dx((max_ - min_) / vector(size_)),
		integrator(integrator_),
		//need non-void constructor call, but want an array ...
		geomGrid0(size_),
		geomGrid1(size_),
		matterGrid(size_),
		auxGrid(size_),
		geomGridReadCurrent(&geomGrid0),
		geomGridWriteCurrent(&geomGrid1)
	{
		integrator->init(this, size);
	}

	//aux calculations

	//calculates partial_gamma_lll, conn_lll, conn_ull, R_ll, R
	// depends on gamma_ll, gamma_uu
	void calcConnections(
		const tensor_sl AuxCell::*gammaField_ll, 
		const tensor_su AuxCell::*gammaField_uu,
		tensor_lsl AuxCell::*partial_gammaField_lll,
		tensor_lsl AuxCell::*connField_lll,
		tensor_usl AuxCell::*connField_ull,
		tensor_sl AuxCell::*RField_ll,
		real AuxCell::*RField)
	{
		typename AuxGrid::iterator iter;	
		
		//partial_gamma_lll depends on gamma_ll
		//conn_lll depends on partial_gamma_lll
		//conn_ull depends on conn_lll and gamma_uu
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
			
			const tensor_su &gamma_uu = cell.*gammaField_uu;
		
			//partial_gamma_lll(k,i,j) := partial_k gamma_ij
			tensor_lsl &partial_gamma_lll = cell.*partial_gammaField_lll;
			partial_gamma_lll = partialDerivative(auxGrid, gammaField_ll, dx, iter.index);

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
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;

			const tensor_su &gamma_uu = cell.*gammaField_uu;
			const tensor_usl &conn_ull = cell.*connField_ull;
			const tensor_lsl &conn_lll = cell.*connField_lll;

			//partial2_gamma_llll(k,l,i,j) = partial_k partial_l gamma_ij
			tensor_slsl partial2_gamma_llll = partialSecondDerivative(auxGrid, gammaField_ll, auxGrid, partial_gammaField_lll, dx, iter.index);

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

	virtual GeomGrid *getGeomGridReadCurrent() {
		return geomGridReadCurrent;
	}

	virtual GeomGrid *getGeomGridWriteCurrent() {
		return geomGridWriteCurrent;
	}
	
	//iteration
	virtual void getGeometridynamicPartials(
		real dt, 
		const GeomGrid &geomGridRead,	//read from this.  last iteration state.
		GeomGrid &partial_t_geomGrid)			//next iteration partials
	{
		typename AuxGrid::iterator iter;
		
		//first compute and store aux values that will be subsequently used for partial differentiation
		//during this process read and write to the same cell
		
		//psi = exp(ln(psi))
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
			const GeomCell &geomCell = geomGridRead(iter.index);
		
			//ln_psi := ln(psi) = 1/6 ln(sqrt(gamma))
			//ln(psi) = 1/6 ln(sqrt(gamma))
			cell.ln_psi = geomCell.ln_sqrt_gamma / 6.;
			
			//psi = exp(ln(psi))
			cell.psi = exp(cell.ln_psi);
		}

#if 0	//option #1: the original ADM implementation
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
			const GeomCell &geomCell = geomGridRead(iter.index);
	
			//gamma = det(gamma_ij)
			cell.gamma = determinant(geomCell.gamma_ll);
			
			//gamma^ij = ((gamma_kl)^-1)^ij
			cell.gamma_uu = inverse(geomCell.gamma_ll, cell.gamma);
		}
#endif
#if 1	//option #2: introduction of conformal factor
		//	calculate gammaBar^ij and gammaBar_ij from psi and ln(sqrt(gamma))
		//	calculate gamma^ij from gammaBar^ij from gammaBar_ij
		//	calculate gamma from ln(sqrt(gamma))
			
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
			const GeomCell &geomCell = geomGridRead(iter.index);
		
			real psiSquared = cell.psi * cell.psi;
			real psiToTheFourth = psiSquared * psiSquared;
			real oneOverPsiToTheFourth = 1. / psiToTheFourth;

			//gamma = psi^12
			//either this or another exp() call
			real psiToTheEighth = psiToTheFourth * psiToTheFourth;
			cell.gamma = psiToTheFourth * psiToTheEighth;

			//gammaBar_ij = psi^-4 gamma_ij
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					cell.gammaBar_ll(i,j) = oneOverPsiToTheFourth * geomCell.gamma_ll(i,j);
				}
			}

			//gammaBar^ij = inverse(gammaBar_ij)
			cell.gammaBar_uu = inverse(cell.gammaBar_ll, 1.);
		}

		//gamma^ij = psi^-4 gammaBar^ij
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
		
			real psiSquared = cell.psi * cell.psi;
			real psiToTheFourth = psiSquared * psiSquared;
			real oneOverPsiToTheFourth = 1. / psiToTheFourth;

			//gamma^ij = psi^-4 gammaBar^ij
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					cell.gamma_uu(i,j) = oneOverPsiToTheFourth * cell.gammaBar_uu(i,j);
				}
			}
		}
#endif
		
		//beta_l, gamma_uu, D_alpha_l
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
			const GeomCell &geomCell = geomGridRead(iter.index);
			
			const tensor_u &beta_u = geomCell.beta_u;
			const tensor_sl &gamma_ll = geomCell.gamma_ll;
	
			//D_alpha_l(i) := diff_i alpha = partial_i alpha
			tensor_l &D_alpha_l = cell.D_alpha_l;
			D_alpha_l = partialDerivative(geomGridRead, &GeomCell::alpha, dx, iter.index);

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
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
	
			//DBar_ln_psi_l(i) := DBar_i ln(psi) = partial_i ln(psi)
			cell.DBar_ln_psi_l = partialDerivative(auxGrid, &AuxCell::ln_psi, dx, iter.index);

			//DBar_psi_l(i) := DBar_i psi
			cell.DBar_psi_l = partialDerivative(auxGrid, &AuxCell::psi, dx, iter.index);
		}

		//calculates partial_gammaBar_lll, connBar_lll, connBar_ull, RBar_ll, RBar
		// depends on gammaBar_ll, gammaBar_uu
		calcConnections(&AuxCell::gammaBar_ll, &AuxCell::gammaBar_uu, &AuxCell::partial_gammaBar_lll, &AuxCell::connBar_lll, &AuxCell::connBar_ull, &AuxCell::RBar_ll, &AuxCell::RBar);

#if 0	//option #1: original
		//TODO if you want to use this here and the conformal-based Hamiltonian later, you have to move a few aux field calculations out of the commented block below 

		//calculates partial_gamma_lll, conn_lll, conn_ull, R_ll, R
		// depends on gamma_ll, gamma_uu
		calcConnections(auxGrid, &AuxCell::gamma_ll, &AuxCell::gamma_uu, &AuxCell::partial_gamma_lll, &AuxCell::conn_lll, &AuxCell::conn_ull, &AuxCell::R_ll, &AuxCell::R);
#endif
#if 1	//use conformal metric info
	
		//conn_ull depends on connBar_ull, DBar_ln_psi_l, gammaBar_uu, gammaBar_ll
		//conn_lll depends on conn_ull and gamma_ll
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
			const GeomCell &geomCell = geomGridRead(iter.index);

			const tensor_l &DBar_ln_psi_l = cell.DBar_ln_psi_l;
			const tensor_sl &gammaBar_ll = cell.gammaBar_ll;
			const tensor_su &gammaBar_uu = cell.gammaBar_uu;
			const tensor_usl &connBar_ull = cell.connBar_ull;
			const tensor_sl &gamma_ll = geomCell.gamma_ll;
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
		
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;

			const real &psi = cell.psi;
			const real &RBar = cell.RBar;
			const tensor_l &DBar_ln_psi_l = cell.DBar_ln_psi_l;
			const tensor_l &DBar_psi_l = cell.DBar_psi_l;
			const tensor_su &gammaBar_uu = cell.gammaBar_uu;
			const tensor_sl &gammaBar_ll = cell.gammaBar_ll;
			const tensor_sl &RBar_ll = cell.RBar_ll;
			const tensor_usl &connBar_ull = cell.connBar_ull;

			//partial2_ln_psi_ll(i,j) := partial_i partial_j ln(psi)
			tensor_sl partial2_ln_psi_ll = partialSecondDerivative(auxGrid, &AuxCell::ln_psi, auxGrid, &AuxCell::DBar_ln_psi_l, dx, iter.index);

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
			tensor_sl partial2_psi_ll = partialSecondDerivative(auxGrid, &AuxCell::psi, auxGrid, &AuxCell::DBar_psi_l, dx, iter.index);

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

		//calculate extrinsic curvature tensors
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
			const GeomCell &geomCell = geomGridRead(iter.index);
	
			const real &K = geomCell.K;
			const tensor_sl &gamma_ll = geomCell.gamma_ll;
			const tensor_sl &K_ll = geomCell.K_ll;
			
			const real &psi = cell.psi;
			const tensor_su &gamma_uu = cell.gamma_uu;
			const tensor_su &gammaBar_uu = cell.gammaBar_uu;

			//A_ll(i,j) := A_ij = K_ij - 1/3 gamma_ij K
			tensor_sl &A_ll = cell.A_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					A_ll(i,j) = K_ll(i,j) - gamma_ll(i,j) * K / 3.;
				}
			}
		
			real psiSquared = psi * psi;

			//ABar_ll(i,j) := ABar_ij = psi^2 A_ij
			tensor_sl &ABar_ll = cell.ABar_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					ABar_ll(i,j) = A_ll(i,j) * psiSquared;
				}
			}

			//ABar_ul(i,j) := ABar^i_j = gammaBar^ik ABar_kj
			tensor_ul ABar_ul;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					ABar_ul(i,j) = 0;
					for (int k = 0; k < dim; ++k) {
						ABar_ul(i,j) += gammaBar_uu(i,k) * ABar_ll(k,j);
					}
				}
			}

			//ABar_uu(i,j) := ABar^ij = ABar^i_k gammaBar^kj
			tensor_su &ABar_uu = cell.ABar_uu;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					ABar_uu(i,j) = 0;
					for (int k = 0; k < dim; ++k) {
						ABar_uu(i,j) += ABar_ul(i,k) * gammaBar_uu(k,j);
					}
				}
			}

			//tr_ABar_sq := tr(ABar^2) = ABar_ij ABar^ji
			real &tr_ABar_sq = cell.tr_ABar_sq;
			tr_ABar_sq = 0;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					tr_ABar_sq += ABar_ll(i,j) * ABar_uu(i,j);
				}
			}

			//K^i_j := gamma^ik K_kj
			tensor_ul &K_ul = cell.K_ul;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					K_ul(i,j) = 0;
					for (int k = 0; k < dim; ++k) {
						K_ul(i,j) += gamma_uu(i,k) * K_ll(k,j);
					}
				}
			}
		
			//K_uu(i,j) := K^ij = K^i_k gamma^kj
			tensor_su &K_uu = cell.K_uu;
			for (int i = 0; i  < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					K_uu(i,j) = 0;
					for (int k = 0; k < dim; ++k) {
						K_uu(i,j) += K_ul(i,k) * gamma_uu(k,j);
					}
				}
			}
			
			//tr_K_sq := tr(K^2) = (K^2)^i_i = K^ij K_ji = K^i_j K^j_i
			//this method uses tr(K^2) = K^ij K_ij in particular
			//tr_K_sq := tr(K^2) = K^ij K_ij
			real &tr_K_sq = cell.tr_K_sq;
			tr_K_sq = 0.;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					tr_K_sq += K_uu(i,j) * K_ll(i,j); 
				}
			}
		}

		//calcConstraints
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
			const GeomCell &geomCell = geomGridRead(iter.index);
			const MatterCell &matterCell = matterGrid(iter.index);
	
			const real &psi = cell.psi;
			const real &DBar2_psi = cell.DBar2_psi;
			const real &tr_ABar_sq = cell.tr_ABar_sq;
			const real &RBar = cell.RBar;
			const real &K = geomCell.K;
			const real &rho = matterCell.rho;
			const tensor_u &S_u = matterCell.S_u;

			//Hamiltonian constraint
#if 0		//option #1: use original ADM method

			const real &R = cell.R;
			const real &tr_K_sq = cell.tr_K_sq;
			
			//H = 1/2 (R + K^2 - K^j_i K^i_j) - 8 pi rho
			real &H = cell.H;
			H = .5 * (R + K * K - tr_K_sq) - 8. * M_PI * rho;
#endif
#if 0		//option #2: use conformal values
	
			real psiSquared = psi * psi;
			real psiToTheFourth = psiSquared * psiSquared;
			
			//"Numerical Relativity" p.57
			//H = 1/2 ( -8 DBar^2 psi + psi RBar + psi^5 (K^2 - tr(K^2)) - 16 pi psi^5 rho)
			//  = 1/2 ( -8 DBar^2 psi + psi (RBar + psi^4 (K^2 - tr(K^2) - 16 pi rho)))
			real &H = cell.H;
			H = .5 * (-8. * DBar2_psi + psi * (RBar + psiToTheFourth * (K * K - tr_K_sq - 16. * M_PI * rho)));
#endif
#if 1		//option #3: use conformal factor and conformal traceless extrinsic curvature
			{
				real psiSquared = psi * psi;
				real psiToTheFourth = psiSquared * psiSquared;
				real psiToTheEighth = psiToTheFourth * psiToTheFourth;
			
				//"Numerical Relativity" p.65
				//H = 1/2 (-8 DBar^2 psi + psi RBar + 2/3 psi^5 K^2 - psi^-7 tr(ABar^2) - 16 pi psi^5 rho)
				//  = 1/2 (-8 DBar^2 psi + psi (RBar + psi^4 (2/3 K^2 - 16 pi rho) - psi^-8 tr(ABar^2)))
				real &H = cell.H;
				H = .5 * (-8. * DBar2_psi + psi * (RBar + psiToTheFourth * ((2./3.) * K * K - 16. * M_PI * rho) - tr_ABar_sq / psiToTheEighth));
			}
#endif

#if 0		//option #1: original momentum constraint
			{
				const tensor_su &gamma_uu = cell.gamma_uu;
				
				//diff_K_uu(k,i,j) := diff_k K^ij
				//I'm not seeing a distinction between the 4D and 3D representations of the constraints.
				//The "Numerical Relativity" book did transition from a covariant form of the 4D constraint to a contravariant form of the 3D constraint
				//It also added caveats on how D_i v_j = delta_i v_j only for rank-(0,2) forms and only if v was purely spatial.
				tensor_lsu diff_K_luu = covariantDerivative(auxGrid, &AuxCell::K_uu, auxGrid, &AuxCell::conn_ull, dx, iter.index);

				//partial_K_l(i) := partial_i K
				tensor_l partial_K_l = partialDerivative(auxGrid, &AuxCell::K, dx, iter.index);

				//M_u(i) := M^i = D_j (K^ij - gamma^ij K) - 8 pi S^i
				tensor_u &M_u = cell.M_u;
				for (int i = 0; i < dim; ++i) {
					M_u(i) = -8. * M_PI * S_u(i);
					for (int j = 0; j < dim; ++j) {
						M_u(i) += diff_K_luu(j,i,j) - gamma_uu(i,j) * partial_K_l(j);
					}
				}
			}
#endif
#if 1		//option #1: use conformal traceless extrinsic curvature
			{
				const tensor_su &gammaBar_uu = cell.gammaBar_uu;
				
				//DBar_K_l(i) := DBar_i K
				tensor_l DBar_K_l = partialDerivative(geomGridRead, &GeomCell::K, dx, iter.index);

				//DBar_ABar_luu(i,j,k) := DBar_i ABar^jk
				tensor_lsu DBar_ABar_luu = covariantDerivative(auxGrid, &AuxCell::ABar_uu, auxGrid, &AuxCell::connBar_ull, dx, iter.index);

				real psiSquared = psi * psi;
				real psiToTheFourth = psiSquared * psiSquared;
				real psiToTheSixth = psiSquared * psiToTheFourth;
				real psiToTheTenth = psiToTheFourth * psiToTheSixth;

				//M_u(i) := M^i = DBar_j ABar^ij - 2/3 psi^6 gammaBar^ij DBar_j K - 8 pi psi^10 S^i
				tensor_u &M_u = cell.M_u;
				for (int i = 0; i < dim; ++i) {
					M_u(i) = -8 * M_PI * psiToTheTenth * S_u(i);
					for (int j = 0; j < dim; ++j) {
						M_u(i) += DBar_ABar_luu(j,i,j) - 2./3. * psiToTheSixth * gammaBar_uu(i,j) * DBar_K_l(j);
					}
				}
			}
#endif
		}

		//calcPartial
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
			const GeomCell &geomCell = geomGridRead(iter.index);
			const MatterCell &matterCell = matterGrid(iter.index);
			
			const real &rho = matterCell.rho;
			const tensor_sl &S_ll = matterCell.S_ll;
			
			const real &alpha = geomCell.alpha;
			const tensor_u &beta_u = geomCell.beta_u;
			const tensor_sl &gamma_ll = geomCell.gamma_ll;
			const real &K = geomCell.K;
			const tensor_sl &K_ll = geomCell.K_ll;
			
			const tensor_l &D_alpha_l = cell.D_alpha_l;
			const tensor_usl &conn_ull = cell.conn_ull;
			const tensor_sl &R_ll = cell.R_ll;
			const tensor_su &gamma_uu = cell.gamma_uu;
			const tensor_ul &K_ul = cell.K_ul;
			const real &tr_K_sq = cell.tr_K_sq;

			//partial2_alpha_ll(i,j) := partial_i partial_j alpha
			tensor_sl partial2_alpha_ll = partialSecondDerivative(geomGridRead, &GeomCell::alpha, auxGrid, &AuxCell::D_alpha_l, dx, iter.index);

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

			GeomCell &partial_t_geomCell = partial_t_geomGrid(iter.index);

#if 0	//original ADM advects gamma_ij.  instead we are advecting the conformal factor.
			
			//D_beta_ll(j,i) := D_j beta_i = partial_j beta_i - conn^k_ij beta_k
			tensor_ll D_beta_ll = covariantDerivative(auxGrid, &AuxCell::beta_l, auxGrid, &AuxCell::conn_ull, dx, iter.index);
			
			//diff_i beta_j = D_i beta_j 
			//-- but only for lower indexes.
			//doesn't work if either is upper (you get extra terms for the timelike component we're dropping)
			//(see "Numerical Relativity", exercise 2.30)
	
			//partial_t gamma_ij = -2 alpha K_ij + D_i beta_j + D_j beta_i
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					partial_t_geomCell.gamma_ll(i,j) = -2. * alpha * K_ll(i,j) + D_beta_ll(i,j) + D_beta_ll(j,i);
				}
			}
#endif
			//partial_K_lll(k,i,j) := partial_k K_ij
			tensor_lsl partial_K_lll = partialDerivative(geomGridRead, &GeomCell::K_ll, dx, iter.index);

			//partial_beta_lu(j,i) := partial_j beta^i
			tensor_lu partial_beta_lu = partialDerivative(geomGridRead, &GeomCell::beta_u, dx, iter.index);

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
					partial_t_geomCell.K_ll(i,j) = 0.;
					partial_t_geomCell.K_ll(i,j) += alpha * R_ll(i,j);
					for (int k = 0; k < dim; ++k) {
						partial_t_geomCell.K_ll(i,j) -= alpha * 2. * K_ll(i,k) * K_ul(k,j);
					}
					partial_t_geomCell.K_ll(i,j) += alpha * K * K_ll(i,j) 
						- D2_alpha_ll(i,j)
						- 8. * M_PI * alpha * (
							S_ll(i,j) - .5 * gamma_ll(i,j) * (S - rho));
					for (int k = 0; k < dim; ++k) {
						partial_t_geomCell.K_ll(i,j) += beta_u(k) * partial_K_lll(k,i,j)
							+ K_ll(k,i) * partial_beta_lu(j,k)
							+ K_ll(k,j) * partial_beta_lu(i,k);
					}
				}
			}

			//D_beta_lu(i,j) := D_i beta^j = D_i beta^j
			//why am I suspicious that, despite the above comment, I am doing just what the book says not to do?
			tensor_lu D_beta_lu = covariantDerivative(geomGridRead, &GeomCell::beta_u, auxGrid, &AuxCell::conn_ull, dx, iter.index);

			//trace_D_beta := D_i beta^i
			real trace_D_beta = 0;
			for (int i = 0; i < dim; ++i) {
				trace_D_beta += D_beta_lu(i,i);
			}

			//partial_t ln_sqrt_gamma = -alpha K + D_i beta^i
			partial_t_geomCell.ln_sqrt_gamma = -alpha * K + trace_D_beta;

			//D_K_l(i) := D_i K
			tensor_l D_K_l = covariantDerivative(geomGridRead, &GeomCell::K, auxGrid, &AuxCell::conn_ull, dx, iter.index);

			//partial_t K = -gamma^ij D_i D_j alpha + alpha(K_ij K^ij + 4 pi (rho + S)) + beta^i D_i K
			partial_t_geomCell.K = 0;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					partial_t_geomCell.K += -gamma_uu(i,j) * D2_alpha_ll(i,j);
				}
				partial_t_geomCell.K += beta_u(i) * D_K_l(i);
			}
			partial_t_geomCell.K += alpha * (tr_K_sq + 4. * M_PI * (rho + S));
		}
	}

	void update(real dt) {
		integrator->update(dt);

		time += dt;

		//TODO update fluid components

		//and swap source and dest grids
		//do something more clever if we ever get any more than 2 histories
		//if you never need more than 2 then maybe we won't need histories at all, just partials?
		std::swap(geomGridReadCurrent, geomGridWriteCurrent);
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
		for (typename AuxGrid::iterator iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			const AuxCell &cell = *iter;
			const GeomCell &geomCell = (*geomGridReadCurrent)(iter.index);

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
			o << geomCell.alpha << "\t";
			
			//beta^i
			for (int i = 0; i < dim; ++i) {
				o << geomCell.beta_u(i) << "\t";
			}

			//beta_i
			for (int i = 0; i < dim; ++i) {
				o << cell.beta_l(i) << "\t";
			}

			//ln(sqrt(gamma))
			o << geomCell.ln_sqrt_gamma << "\t";

			//ln(psi)
			o << cell.ln_psi << "\t";

			//psi
			o << cell.psi << "\t";

			//gamma_ij
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					o << geomCell.gamma_ll(i,j) << "\t";
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
					o << geomCell.K_ll(i,j) << "\t";
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
			o << geomCell.K << "\t";

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

