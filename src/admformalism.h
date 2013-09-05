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

//Formalism class generates partial t values
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
	
	/*
	stencil
	for implicit integration
	used like so:

	stencilCell(offset, offsetof(GeomCell, destField) / sizeof(real), offsetof(GeomCell, sourceField) / sizeof(real)) = coefficient
		where destField is the field we're computing the partial derivative of wrt time
		and sourceField is the field that is contributing to a term of destField

	and for multi-rank structures:
	stencilCell(offset, offsetof(GeomCell, destField) / sizeof(real) + GeomCell::destField::index(i,j), offsetof(GeomCell, sourceField) / sizeof(real) + GeomCell::sourceField::index(k,l)) = coefficient
	
	if d/dt phi = -alpha K / 6 + partial_i beta^i / 6 + beta^i partial_i phi
	
	so you would probably write 

#define GEOMCELL_OFFSET(field) (offsetof(GeomCell, field)/sizeof(real))
	stencil(center, GEOMCELL_OFFSET(phi), GEOMCELL_OFFSET(alpha)) = -K/6	//or .K = -alpha/6 or half and half
	for (int i = 0; i < dim; ++i) {
		for (partialIndex = 0; partialIndex < partialCoefficient.size; ++partialIndex) {
			stencil(center + dxi(i) * partialIndex, GEOMCELL_OFFSET(phi), GEOMCELL_OFFSET(beta_u) + i) = partialCoefficient(partialIndex) / 6
			stencil(center + dxi(i) * partialIndex, GEOMCELL_OFFSET(phi), GEOMCELL_OFFSET(phi)) += beta_u(i) * partialCoefficient(partialIndex)
		}
	}
	
	mind you solving our system as an implicit one is difficult considering how nonlinear it all is
	options are to split up coefficients or to use a nonlinear solver

	using a dense matrix, even for 1D (6 elements) is 36 components times stencil size to the dimension power ... per grid cell (too big)
	so this will have to be a sparse matrix structure

	*/

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
		const GeomGrid &geomGrid,
		const tensor_sl GeomCell::*gammaField_ll, 
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
			partial_gamma_lll = partialDerivative(geomGrid, gammaField_ll, dx, iter.index);

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
			tensor_slsl partial2_gamma_llll = partialSecondDerivative(
				geomGrid, 
				gammaField_ll, 
				auxGrid, 
				partial_gammaField_lll, 
				dx, 
				iter.index);

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
	void calcAux(const GeomGrid &geomGridRead) {
		typename AuxGrid::iterator iter;
		
		//first compute and store aux values that will be subsequently used for partial differentiation
		//during this process read and write to the same cell
		
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
			const GeomCell &geomCell = geomGridRead(iter.index);
	
			//psi = exp(phi)
			cell.psi = exp(geomCell.phi);
			real psiSquared = cell.psi * cell.psi;
			real psiToTheFourth = psiSquared * psiSquared;
			real psiToTheEighth = psiToTheFourth * psiToTheFourth;

			//gamma = psi^12
			cell.gamma = psiToTheFourth * psiToTheEighth;

			//at the moment I'm iterating ln(sqrt(gamma)) and gamma_ij separately
			//I could compare ln(sqrt(gamma)) versus det(gamma_ij) to see how accurate my integrator is
			// but I already know who will win: ln(sqrt(gamma))
			//Instead, how about I remove the trace from gamma_ij and re-insert ln(sqrt(gamma)) ?

			//gamma_ij = psi^4 gammaBar_ij
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					cell.gamma_ll(i,j) = psiToTheFourth * geomCell.gammaBar_ll(i,j);
				}
			}
		
			//gammaBar^ij = inverse(gammaBar_ij)
			cell.gammaBar_uu = inverse(geomCell.gammaBar_ll, 1.);

			real oneOverPsiToTheFourth = 1. / psiToTheFourth;
			//gamma^ij = psi^-4 gammaBar^ij
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					cell.gamma_uu(i,j) = oneOverPsiToTheFourth * cell.gammaBar_uu(i,j);
				}
			}
		}
		
		//beta_l, gamma_uu, D_alpha_l
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
			const GeomCell &geomCell = geomGridRead(iter.index);
			
			const tensor_u &beta_u = geomCell.beta_u;
			const tensor_sl &gamma_ll = cell.gamma_ll;
	
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

		//DBar_phi_l depends on phi
		//DBar_psi_l depends on psi
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
	
			//DBar_phi_l(i) := DBar_i ln(psi) = partial_i ln(psi)
			cell.DBar_phi_l = partialDerivative(geomGridRead, &GeomCell::phi, dx, iter.index);

			//DBar_psi_l(i) := DBar_i psi
			cell.DBar_psi_l = partialDerivative(auxGrid, &AuxCell::psi, dx, iter.index);
		}

		//calculates partial_gammaBar_lll, connBar_lll, connBar_ull, RBar_ll, RBar
		// depends on gammaBar_ll, gammaBar_uu
		calcConnections(geomGridRead, &GeomCell::gammaBar_ll, &AuxCell::gammaBar_uu, &AuxCell::partial_gammaBar_lll, &AuxCell::connBar_lll, &AuxCell::connBar_ull, &AuxCell::RBar_ll, &AuxCell::RBar);

		//conn_ull depends on connBar_ull, DBar_phi_l, gammaBar_uu, gammaBar_ll
		//conn_lll depends on conn_ull and gamma_ll
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
			const GeomCell &geomCell = geomGridRead(iter.index);

			const tensor_l &DBar_phi_l = cell.DBar_phi_l;
			const tensor_sl &gammaBar_ll = geomCell.gammaBar_ll;
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
						if (i == j) conn_ull(i,j,k) += 2. * DBar_phi_l(k);
						if (i == k) conn_ull(i,j,k) += 2. * DBar_phi_l(j);
						for (int l = 0; l < dim; ++l) {
							conn_ull(i,j,k) -= 2. * gammaBar_ll(j,k) * gammaBar_uu(i,l) * DBar_phi_l(l);
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
			const GeomCell &geomCell = geomGridRead(iter.index);

			const real &psi = cell.psi;
			const real &RBar = cell.RBar;
			const tensor_l &DBar_phi_l = cell.DBar_phi_l;
			const tensor_l &DBar_psi_l = cell.DBar_psi_l;
			const tensor_su &gammaBar_uu = cell.gammaBar_uu;
			const tensor_sl &gammaBar_ll = geomCell.gammaBar_ll;
			const tensor_sl &RBar_ll = cell.RBar_ll;
			const tensor_usl &connBar_ull = cell.connBar_ull;

			//partial2_phi_ll(i,j) := partial_i partial_j ln(psi)
			tensor_sl partial2_phi_ll = partialSecondDerivative(
				geomGridRead,
				&GeomCell::phi,
				auxGrid,
				&AuxCell::DBar_phi_l,
				dx,
				iter.index);

			//DBar2_phi_ll(i,j) := DBar_i DBar_j ln(psi) = partial_i partial_j ln(psi) - connBar^k_ij partial_k ln(psi)
			tensor_sl DBar2_phi_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					DBar2_phi_ll(i,j) = partial2_phi_ll(i,j);
					for (int k = 0; k < dim; ++k) {
						DBar2_phi_ll(i,j) -= connBar_ull(k,i,j) * DBar_phi_l(k);
					}
				}
			}

			//DBar2_phi := DBar^2 ln(psi) = gammaBar^ij DBar_i DBar_j ln(psi)
			real DBar2_phi = 0;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					DBar2_phi += gammaBar_uu(i,j) * DBar2_phi_ll(i,j);
				}
			}

			//normBar_DBar_phi = gammaBar^ij (DBar_i ln(psi)) (DBar_j ln(psi))
			real normBar_DBar_phi = 0;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					normBar_DBar_phi += gammaBar_uu(i,j) * DBar_phi_l(i) * DBar_phi_l(j);
				}
			}

			//"Numerical Relativity" p.57
			//R_ll(i,j) := R_ij = RBar_ij - 2 (DBar_i DBar_j ln(psi) + gammaBar_ij gammaBar^lm DBar_l DBar_m ln(psi)) + 4((DBar_i ln(psi)) (DBar_j ln(psi)) - gammaBar_ij gammaBar^lm (DBar_l ln(psi)) (DBar_m ln(psi)))
			tensor_sl &R_ll = cell.R_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					R_ll(i,j) = RBar_ll(i,j) - 2. * (DBar2_phi_ll(i,j) + gammaBar_ll(i,j) * DBar2_phi) + 4. * (DBar_phi_l(i) * DBar_phi_l(j) - gammaBar_ll(i,j) * normBar_DBar_phi);
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

		//calculate extrinsic curvature tensors
		for (iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
			const GeomCell &geomCell = geomGridRead(iter.index);

			const real &K = geomCell.K;
			const tensor_sl &gamma_ll = cell.gamma_ll;
			const tensor_sl &ATilde_ll = geomCell.ATilde_ll;
			
			const real &psi = cell.psi;
			const tensor_su &gamma_uu = cell.gamma_uu;
			const tensor_su &gammaBar_uu = cell.gammaBar_uu;


			//ATilde_ul(i,j) := ATilde^i_j = gammaBar^ik ATilde_kj
			tensor_ul &ATilde_ul = cell.ATilde_ul;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					ATilde_ul(i,j) = 0;
					for (int k = 0; k < dim; ++k) {
						ATilde_ul(i,j) += gammaBar_uu(i,k) * ATilde_ll(k,j);
					}
				}
			}

			real psiSquared = psi * psi;
			real psiToTheFourth = psiSquared * psiSquared;

			//A_ll(i,j) := ATilde_ij = exp(4phi) ATilde_ij = psi^4 ATilde_ij
			tensor_sl &A_ll = cell.A_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					A_ll(i,j) = ATilde_ll(i,j) * psiToTheFourth;
				}
			}

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

			//K_ll(i,j) := K_ij = A_ij + 1/3 gamma_ij K
			tensor_sl &K_ll = cell.K_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					K_ll(i,j) = A_ll(i,j) + 1./3. * gamma_ll(i,j) * K;
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
			
			const real &K = geomCell.K;
	
			const real &psi = cell.psi;
			const real &DBar2_psi = cell.DBar2_psi;
			const real &RBar = cell.RBar;
			const tensor_ul &ATilde_ul = cell.ATilde_ul;
			const tensor_su &gammaBar_uu = cell.gammaBar_uu;
			
			const real &rho = matterCell.rho;
			const tensor_u &S_u = matterCell.S_u;

				//Hamiltonian constraint
			
			real psiSquared = psi * psi;
			real psiToTheFourth = psiSquared * psiSquared;
			real psiToTheSixth = psiSquared * psiToTheFourth;

			//tr_ATilde_sq := tr(ATilde^2) = ATilde_ij ATilde^ji
			real tr_ATilde_sq = 0;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					tr_ATilde_sq += ATilde_ul(i,j) * ATilde_ul(j,i);
				}
			}

			//Baumgarte & Shapiro p.390
			//H = gammaBar^ij DBar_i DBar_j exp(phi) - exp(phi)/8 RBar + exp(5phi)/8 ATilde_ij ATilde^ij - exp(5phi)/12 K^2 + 2 pi exp(5phi) rho
			//  = DBar^2 psi + psi( -RBar / 8 + psi^4 (tr(ATilde^2) / 8 - K^2 / 12 + 2 pi rho)
			real &H = cell.H;
			H = DBar2_psi + psi * (-RBar / 8. + psiToTheFourth * (tr_ATilde_sq / 8. - K * K / 12. + 2. * M_PI * rho));

				//momentum constraint
				
			//DBar_ABar_luu(i,j,k) := DBar_i ABar^jk
			tensor_lsu DBar_ABar_luu = covariantDerivative(auxGrid, &AuxCell::ABar_uu, auxGrid, &AuxCell::connBar_ull, dx, iter.index);
	
			//DBar_ABar_u(i) := DBar_j ABar^ji
			tensor_u DBar_ABar_u;
			for (int i = 0; i < dim; ++i) {
				DBar_ABar_u(i) = 0;
				for (int j = 0; j < dim; ++j) {
					DBar_ABar_u(i) += DBar_ABar_luu(i,j,i);
				}
			}
			
			//DBar_K_l(i) := DBar_i K
			tensor_l DBar_K_l = partialDerivative(geomGridRead, &GeomCell::K, dx, iter.index);

			//DBar_K_u(i) := DBar^i K = gammaBar^ij DBar_j K
			tensor_u DBar_K_u;
			for (int i = 0; i < dim; ++i) {
				DBar_K_u(i) = 0;
				for (int j = 0; j < dim; ++j) {
					DBar_K_u(i) += gammaBar_uu(i,j) * DBar_K_l(j);
				}
			}

			//Baumgarte & Shapiro 
			//p.65
			//M_u(i) := M^i = DBar_j ABar^ij - 2/3 psi^6 gammaBar^ij DBar_j K - 8 pi psi^10 S^i
			//p.390
			//M_u(i) := M^i = DBar_j (psi^6 ATilde^ij) - 2/3 psi^6 DBar^i K - 8 pi psi^6 S^i
			//		  = DBar_j ABar^ij + psi^6( - 2/3 gammaBar^ij DBar_j K - 8 pi S^i)
			//hmm, looks like that psi^10 turned into an exp(6phi) ... wonder why that isn't an exp(10phi) ...
			tensor_u &M_u = cell.M_u;
			for (int i = 0; i < dim; ++i) {
				M_u(i) = DBar_ABar_u(i) + psiToTheSixth * (-2./3. * DBar_K_u(i) - 8. * M_PI * S_u(i));
			}
		}
	}


	virtual void getExplicitPartials(
		real dt, 
		const GeomGrid &geomGridRead,	//read from this.  last iteration state.
		GeomGrid &partial_t_geomGrid)		//next iteration partials
	{
		calcAux(geomGridRead);
	
		for (typename AuxGrid::iterator iter = auxGrid.begin(); iter != auxGrid.end(); ++iter) {
			AuxCell &cell = *iter;
			const GeomCell &geomCell = geomGridRead(iter.index);
			const MatterCell &matterCell = matterGrid(iter.index);
			
			const real &rho = matterCell.rho;
			const tensor_sl &S_ll = matterCell.S_ll;
			
			const real &alpha = geomCell.alpha;
			const tensor_u &beta_u = geomCell.beta_u;
			const tensor_sl &gammaBar_ll = geomCell.gammaBar_ll;
			const real &K = geomCell.K;
			const tensor_sl &ATilde_ll = geomCell.ATilde_ll;
			
			const tensor_ul &ATilde_ul = cell.ATilde_ul;	
			const real &psi = cell.psi;
			const tensor_l &D_alpha_l = cell.D_alpha_l;
			const tensor_usl &conn_ull = cell.conn_ull;
			const tensor_sl &R_ll = cell.R_ll;
			const tensor_sl &gamma_ll = cell.gamma_ll;
			const tensor_su &gamma_uu = cell.gamma_uu;
			const tensor_su &gammaBar_uu = cell.gammaBar_uu;
			const tensor_lsl &partial_gammaBar_lll = cell.partial_gammaBar_lll;
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

			//D_beta_ll(j,i) := D_j beta_i = partial_j beta_i - conn^k_ij beta_k
			tensor_ll D_beta_ll = covariantDerivative(auxGrid, &AuxCell::beta_l, auxGrid, &AuxCell::conn_ull, dx, iter.index);

			//partial_beta_lu(j,i) := partial_j beta^i
			tensor_lu partial_beta_lu = partialDerivative(geomGridRead, &GeomCell::beta_u, dx, iter.index);

			//trace_partial_beta := partial_i beta^i
			real trace_partial_beta = 0.;
			for (int i = 0; i < dim; ++i) {
				trace_partial_beta += partial_beta_lu(i,i);
			}

			//partial_t gammaBar_ij = -2 alpha ATilde_ij - 2/3 gammaBar_ij partial_k beta^k
			//		+ beta^k partial_k gammaBar_ij + gammaBar_ik partial_j beta^k 
			//		+ gammaBar_kj partial_i beta^k
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					partial_t_geomCell.gammaBar_ll(i,j) = -2. * alpha * ATilde_ll(i,j) - 2./3. * gammaBar_ll(i,j) * trace_partial_beta;
					for (int k = 0; k < dim; ++k) {
						partial_t_geomCell.gammaBar_ll(i,j) += 
							beta_u(k) * partial_gammaBar_lll(k,i,j)
							+ gammaBar_ll(i,k) * partial_beta_lu(j,k)
							+ gammaBar_ll(k,j) * partial_beta_lu(i,k);
					}
				}
			}

			//S := S^i_i := gamma^ij S_ij
			real S = 0.;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					S += S_ll(i,j) * gamma_uu(i,j);
				}
			}
			
			//partial_ATilde_lll(i,j,k) := partial_i ATilde_jk
			tensor_lsl partial_ATilde_lll = partialDerivative(geomGridRead, &GeomCell::ATilde_ll, dx, iter.index);

			//traceless portion of partial_t ATilde_ij := tracefree(-D^2 alpha + alpha (R_ij - 8 pi S_ij))
			tensor_sl tracelessPortionOfPartialT_ATilde_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					tracelessPortionOfPartialT_ATilde_ll(i,j) = -D2_alpha_ll(i,j) + alpha * (R_ll(i,j) - 8. * M_PI * S_ll(i,j));
				}
			}

			real traceOfTraceless = 0;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					traceOfTraceless += gamma_uu(i,j) * tracelessPortionOfPartialT_ATilde_ll(i,j);
				}
			}

			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					tracelessPortionOfPartialT_ATilde_ll(i,j) -= 1./3. * gamma_ll(i,j) * traceOfTraceless;
				}
			}

			real psiSquared = psi * psi;
			real psiToTheFourth = psiSquared * psiSquared;
			//partial_t ATilde_ij = exp(-4phi) ((-D_i D_j alpha + alpha (R_ij - 8 pi S_ij))^TF + alpha (K ATilde_ij - 2 ATilde_il ATilde^l_j))
			//		+ beta^k partial_k ATilde_ij + ATilde_ik partial_j beta^k + ATilde_kj partial_i beta^k - 2/3 ATilde_ij partial_k beta^k
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					partial_t_geomCell.ATilde_ll(i,j) = psiToTheFourth * tracelessPortionOfPartialT_ATilde_ll(i,j) + alpha * K * ATilde_ll(i,j);
					for (int k = 0; k < dim; ++k) {
						partial_t_geomCell.ATilde_ll(i,j) += 
							-alpha * 2. * ATilde_ll(i,k) * ATilde_ul(k,j)
							+ beta_u(k) * partial_ATilde_lll(k,i,j)
							+ ATilde_ll(i,k) * partial_beta_lu(j,k)
							+ ATilde_ll(k,j) * partial_beta_lu(i,k);
					}
					partial_t_geomCell.ATilde_ll(i,j) -= 2./3. * ATilde_ll(i,j) * trace_partial_beta;
				}
			}

			//partial_phi_l(i) := partial_i phi
			tensor_l partial_phi_l = partialDerivative(geomGridRead, &GeomCell::phi, dx, iter.index);

			//partial_t -alpha K / 6 + beta^i partial_i phi + partial_i beta^i / 6
			partial_t_geomCell.phi = (trace_partial_beta - alpha * K) / 6.;
			for (int i = 0; i < dim; ++i) {
				partial_t_geomCell.phi += beta_u(i) * partial_phi_l(i);
			}

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
};

