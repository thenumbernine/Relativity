#pragma once

#include "cell.h"
#include "Tensor/Inverse.h"
#include "Tensor/Derivative.h"
#include "Tensor/clamp.h"
#include "derivative.h"
#include "i_integrator.h"
#include "i_admformalism.h"
#include "parallel.h"
#include "Tensor/Vector.h"
#include "Tensor/Tensor.h"
#include "Tensor/Grid.h"
#include <assert.h>
#include <math.h>
#include <iostream>

//Formalism class generates partial t values
template<typename Real_, int dim_>
struct ADMFormalism : public IADMFormalism<Real_, dim_> {
	using Real = Real_;
	static constexpr auto dim = dim_;
	
	static constexpr int partialDerivativeOrder = 8;
	
		//cell types

	using GeomCell = ::GeomCell<Real, dim>;
	using MatterCell = ::MatterCell<Real, dim>;
	using AuxCell = ::AuxCell<Real, dim>;
	
	using GeomGrid = Tensor::Grid<GeomCell, dim>;
	using MatterGrid = Tensor::Grid<MatterCell, dim>;
	using AuxGrid = Tensor::Grid<AuxCell, dim>;

		//statically-sized mathematical types
		//(pull from cell where you can)
	
	//used for indexes
	using DerefType = Tensor::intN<dim>;	

	using TensorU = typename AuxCell::TensorU;
	using TensorL = typename AuxCell::TensorL;
	using TensorLL = typename AuxCell::TensorLL;
	using TensorSU = typename AuxCell::TensorSU;
	using TensorLU = typename AuxCell::TensorLU;
	using TensorSL = typename AuxCell::TensorSL;
	using TensorUSL = typename AuxCell::TensorUSL;
    using TensorLSL = Tensor::_tensori<Real, Tensor::index_vec<dim>, Tensor::index_sym<dim>>;
	using TensorLSU = TensorLSL;
	using TensorUL = Tensor::_tensori<Real, Tensor::index_vec<dim>, Tensor::index_vec<dim>>;
	using TensorSLU = Tensor::_tensori<Real, Tensor::index_sym<dim>, Tensor::index_vec<dim>>;
	using TensorSLSL = Tensor::_tensori<Real, Tensor::index_sym<dim>, Tensor::index_sym<dim>>;
	using Vector = Tensor::_vec<Real, dim>;
	
	//what a relative notion...
	//intended for use with output matching up iteration slices 
	Real time;

	//resolution of our grids, stored here as well as in each grid for convenience
	DerefType size;

	//x(0) == min
	Vector min;
	
	//x(size-1) == max
	Vector max;
	
	//range = max - min
	Vector range;
	
	//dx = range / size
	Vector dx;

	//integrator
	IIntegrator<Real, dim> *integrator;

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

	stencilCell(offset, offsetof(GeomCell, destField) / sizeof(Real), offsetof(GeomCell, sourceField) / sizeof(Real)) = coefficient
		where destField is the field we're computing the partial derivative of wrt time
		and sourceField is the field that is contributing to a term of destField

	and for multi-rank structures:
	stencilCell(offset, offsetof(GeomCell, destField) / sizeof(Real) + GeomCell::destField::index(i,j), offsetof(GeomCell, sourceField) / sizeof(Real) + GeomCell::sourceField::index(k,l)) = coefficient
	
	if d/dt phi = -alpha K / 6 + partial_i beta^i / 6 + beta^i partial_i phi
	
	so you would probably write 

#define GEOMCELL_OFFSET(field) (offsetof(GeomCell, field)/sizeof(Real))
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
		Vector const & min_, 
		Vector const & max_, 
		DerefType const & size_,
		IIntegrator<Real, dim> * integrator_)
	:	time(0),
		size(size_),
		min(min_),
		max(max_),
		range(max_ - min_),
		dx((max_ - min_) / Vector(size_)),
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

	//helper functions

	//calculate connBar_u based on gammaBar_ll
	//doesn't trust auxGrid.  instead it performs an inverse operation itself,
	// so only use this on startup after gammaBar_ll is determined.
	void calcConnBar(GeomGrid &geomGrid) {
		Tensor::RangeObj<dim> range = auxGrid.range();
		
		parallel.foreach(range.begin(), range.end(), [&](DerefType index) {
			AuxCell & cell = auxGrid(index);
			GeomCell const & geomCell = geomGrid(index);
			cell.gammaBar_uu = Tensor::inverse(geomCell.gammaBar_ll, (Real)1.);
		});

		parallel.foreach(range.begin(), range.end(), [&](DerefType index) {
			GeomCell &geomCell = geomGrid(index);

			TensorLSU partial_gammaBar_luu = Tensor::partialDerivative<partialDerivativeOrder, Real, dim, TensorSU>(
				index, dx, [&](DerefType index) -> TensorSU {
					index = Tensor::clamp(index, DerefType(), size-1);
					return auxGrid(index).gammaBar_uu;
				});
		
			for (int i = 0; i < dim; ++i) {
				geomCell.connBar_u(i) = 0;
				for (int j = 0; j < dim; ++j) {
					geomCell.connBar_u(i) += partial_gammaBar_luu(j,i,j);
				}
			}
		});
	}

	//calculates partial_gammaBar_lll, connBar_lll, connBar_ull, R_ll, R
	// depends on gamma_ll, gamma_uu
	void calcConnections(GeomGrid const & geomGridRead) {
		//partial_gammaBar_lll depends on gamma_ll
		//connBar_lll depends on partial_gammaBar_lll
		//connBar_ull depends on connBar_lll and gamma_uu
		Tensor::RangeObj<dim> range = auxGrid.range();
		parallel.foreach(range.begin(), range.end(), [&](DerefType index) {
			AuxCell & cell = auxGrid(index);
		
			//partial_gammaBar_lll(k,i,j) := partial_k gamma_ij
			TensorLSL &partial_gammaBar_lll = cell.partial_gammaBar_lll;
			partial_gammaBar_lll = Tensor::partialDerivative<partialDerivativeOrder, Real, dim, TensorSL>(
				index, dx, [&](DerefType index) -> TensorSL {
					index = Tensor::clamp(index, DerefType(), size-1);
					return geomGridRead(index).gammaBar_ll;
				});

			//connBar_lll(i,j,k) := conn_ijk = 1/2 (partial_k gamma_ij + partial_j gamma_ik - partial_i gamma_jk)
			TensorLSL &connBar_lll = cell.connBar_lll;
			for (int k = 0; k < dim; ++k) {
				for (int i = 0; i < dim; ++i) {
					for (int j = 0; j <= i; ++j) {
						connBar_lll(i,j,k) = .5 * (partial_gammaBar_lll(k,i,j) + partial_gammaBar_lll(j,i,k) - partial_gammaBar_lll(i,j,k));
					}
				}
			}

			//connBar_ull(i,j,k) := conn^i_jk = gamma^il conn_ljk
			cell.connBar_ull = cell.gammaBar_uu * connBar_lll;
		});

		parallel.foreach(range.begin(), range.end(), [&](DerefType index) {
			AuxCell & cell = auxGrid(index);
			//partial2_gammaBar_llll(k,l,i,j) = partial_k partial_l gamma_ij
			TensorSLSL partial2_gammaBar_llll = partialSecondDerivative(geomGridRead, &GeomCell::gammaBar_ll, auxGrid, &AuxCell::partial_gammaBar_lll, dx, index);

			//partial2_gammaBar_ll(i,j) := gammaBar^lm partial_l partial_m gammaBar_ij
			cell.partial2_gammaBar_ll = cell.gammaBar_uu.template interior<2>(partial2_gammaBar_llll);
		});
	}

	Vector coordForIndex(DerefType const & index) const {
		return min + ((Vector)index + .5) * dx;
	}

	virtual GeomGrid *getGeomGridReadCurrent() {
		return geomGridReadCurrent;
	}

	virtual GeomGrid *getGeomGridWriteCurrent() {
		return geomGridWriteCurrent;
	}
	
	//iteration
	void calcAux(GeomGrid const & geomGridRead) {
		
		//first compute and store aux values that will be subsequently used for partial differentiation
		//during this process read and write to the same cell
		
		Tensor::RangeObj<dim> range = auxGrid.range();
		parallel.foreach(range.begin(), range.end(), [&](DerefType index) {
			AuxCell & cell = auxGrid(index);
			GeomCell const & geomCell = geomGridRead(index);
	
			//psi = exp(phi)
			cell.psi = exp(geomCell.phi);
			Real psiSquared = cell.psi * cell.psi;
			Real psiToTheFourth = psiSquared * psiSquared;
			Real psiToTheEighth = psiToTheFourth * psiToTheFourth;

			//gamma = psi^12
			cell.gamma = psiToTheFourth * psiToTheEighth;

			//at the moment I'm iterating ln(sqrt(gamma)) and gamma_ij separately
			//I could compare ln(sqrt(gamma)) versus det(gamma_ij) to see how accurate my integrator is
			// but I already know who will win: ln(sqrt(gamma))
			//Instead, how about I remove the trace from gamma_ij and re-insert ln(sqrt(gamma)) ?

			//gamma_ij = psi^4 gammaBar_ij
			cell.gamma_ll = psiToTheFourth * geomCell.gammaBar_ll;
		
			//gammaBar^ij = Tensor::inverse(gammaBar_ij)
			cell.gammaBar_uu = Tensor::inverse(geomCell.gammaBar_ll, (Real)1.);

			//gamma^ij = psi^-4 gammaBar^ij
			cell.gamma_uu = cell.gammaBar_uu * (Real)(1. / psiToTheFourth);
		});
		
		//beta_l, gamma_uu, D_alpha_l
		parallel.foreach(range.begin(), range.end(), [&](DerefType index) {
			AuxCell & cell = auxGrid(index);
			GeomCell const & geomCell = geomGridRead(index);
	
			//D_alpha_l(i) := diff_i alpha = partial_i alpha
			cell.D_alpha_l = Tensor::partialDerivative<partialDerivativeOrder, Real, dim, Real>(
				index, dx, [&](DerefType index) -> Real {
					index = Tensor::clamp(index, DerefType(), size-1);
					return geomGridRead(index).alpha;
				});

			//beta_l(i) := g_ij beta^j
			//exclude sum of beta^t = 0
			//not storing beta_t here, since beta_t = beta^k beta_k
			cell.beta_l = cell.gamma_ll * geomCell.beta_u;
			
			//partial_beta_lu(j,i) := partial_j beta^i
			cell.partial_beta_lu = Tensor::partialDerivative<partialDerivativeOrder, Real, dim, TensorU>(
				index, dx, [&](DerefType index) -> TensorU {
					index = Tensor::clamp(index, DerefType(), size-1);
					return geomGridRead(index).beta_u;
				});
		});

		//DBar_phi_l depends on phi
		//DBar_psi_l depends on psi
		parallel.foreach(range.begin(), range.end(), [&](DerefType index) {
			AuxCell &cell = auxGrid(index);
	
			//DBar_phi_l(i) := DBar_i ln(psi) = partial_i ln(psi)
			cell.DBar_phi_l = Tensor::partialDerivative<partialDerivativeOrder, Real, dim, Real>(
				index, dx, [&](DerefType index) -> Real {
					index = Tensor::clamp(index, DerefType(), size-1);
					return geomGridRead(index).phi;
				});

			//DBar_psi_l(i) := DBar_i psi
			cell.DBar_psi_l = partialDerivative<partialDerivativeOrder, Real, dim, Real>(
				index, dx, [&](DerefType index) -> Real {
					index = Tensor::clamp(index, DerefType(), size-1);
					return auxGrid(index).psi;
				});
		});

		//calculates partial_gammaBar_lll, connBar_lll, connBar_ull
		// depends on gammaBar_ll, gammaBar_uu
		calcConnections(geomGridRead);

		//conn_ull depends on connBar_ull, DBar_phi_l, gammaBar_uu, gammaBar_ll
		//conn_lll depends on conn_ull and gamma_ll
		parallel.foreach(range.begin(), range.end(), [&](DerefType index) {
			AuxCell & cell = auxGrid(index);
			GeomCell const & geomCell = geomGridRead(index);

			auto DBar_phi_u = cell.gammaBar_uu * cell.DBar_phi_l;
	
			//conn^i_jk = connBar^i_jk + 2 (delta^i_j DBar_k ln(psi) + delta^i_k DBar_j ln(psi) - gammaBar_jk gammaBar^il DBar_l ln(psi))
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					for (int k = 0; k <= j; ++k) {
						cell.conn_ull(i,j,k) = cell.connBar_ull(i,j,k);
						if (i == j) cell.conn_ull(i,j,k) += 2. * cell.DBar_phi_l(k);
						if (i == k) cell.conn_ull(i,j,k) += 2. * cell.DBar_phi_l(j);
						cell.conn_ull(i,j,k) -= 2. * DBar_phi_u(i) * geomCell.gammaBar_ll(j,k);
					}
				}
			}

			//conn_ijk = gamma_il conn^l_jk
			cell.conn_lll = cell.gamma_ll * cell.conn_ull;
		});
		
		parallel.foreach(range.begin(), range.end(), [&](DerefType index) {
			AuxCell &cell = auxGrid(index);
			GeomCell const & geomCell = geomGridRead(index);

			//partial2_phi_ll(i,j) := partial_i partial_j ln(psi)
			TensorSL partial2_phi_ll = partialSecondDerivative(
				geomGridRead,
				&GeomCell::phi,
				auxGrid,
				&AuxCell::DBar_phi_l,
				dx,
				index);

			//DBar2_phi_ll(i,j) := DBar_i DBar_j ln(psi) = partial_i partial_j ln(psi) - connBar^k_ij partial_k ln(psi)
			TensorSL DBar2_phi_ll = partial2_phi_ll - cell.DBar_phi_l * cell.connBar_ull;

			//DBar2_phi := DBar^2 ln(psi) = gammaBar^ij DBar_i DBar_j ln(psi)
			Real DBar2_phi = Tensor::interior<2>(cell.gammaBar_uu, DBar2_phi_ll);

			//normBar_DBar_phi = gammaBar^ij (DBar_i ln(psi)) (DBar_j ln(psi))
			Real normBar_DBar_phi = cell.DBar_phi_l * cell.gammaBar_uu * cell.DBar_phi_l;
			
			//partial_connBar_lu(i,j) := partial_i connBar^j
			TensorLU & partial_connBar_lu = cell.partial_connBar_lu;
			partial_connBar_lu = Tensor::partialDerivative<partialDerivativeOrder, Real, dim, TensorU>( 
				index, dx, [&](DerefType index) -> TensorU {
					index = Tensor::clamp(index, DerefType(), size-1);
					return geomGridRead(index).connBar_u;
				});
			
			//connBarSq_ijkl := connBar_ijm connBar^m_kl
			auto connBarSq_llll = cell.connBar_lll * cell.connBar_ull;

			//Baumgarte & Shapiro p.388
			//RBar_ll(i,j) := RBar_ij = 
			//	-1/2 gammaBar^lm partial_m partial_l gammaBar_ij 
			// 		+ gammaBar_k(i partial_j) connBar^k 
			//		+ connBar_(ij)k connBar^k
			//		+ gammaBar^lm (
			//			2 connBar^k_l(i connBar_j)mk 
			//			+ connBar^k_im connBar_klj
			//		)
			TensorSL & RBar_ll = cell.RBar_ll;
			RBar_ll = 
				-1./2. * cell.partial2_gammaBar_ll
				+ makeSym(partial_connBar_lu * geomCell.gammaBar_ll)
				+ makeSym(cell.connBar_lll * geomCell.connBar_u);
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					for (int l = 0; l < dim; ++l) {
						for (int m = 0; m < dim; ++m) {
							RBar_ll(i,j) += 
								  connBarSq_llll(j, m, i, l) * cell.gammaBar_uu(l,m)
								+ connBarSq_llll(i, m, j, l) * cell.gammaBar_uu(l,m);
							for (int k = 0; k < dim; ++k) {
								RBar_ll(i,j) += 
									cell.connBar_ull(k,m,i) * cell.connBar_lll(k,l,j)
								;
							}
						}
					}
				}
			}

			//RBar = gammaBar^ij RBar_ij
			cell.RBar = Tensor::dot(cell.gammaBar_uu, RBar_ll);

			//Baumgarte & Shapiro p.57
			//R_ll(i,j) := R_ij = RBar_ij - 2 (DBar_i DBar_j ln(psi) + gammaBar_ij gammaBar^lm DBar_l DBar_m ln(psi)) + 4((DBar_i ln(psi)) (DBar_j ln(psi)) - gammaBar_ij gammaBar^lm (DBar_l ln(psi)) (DBar_m ln(psi)))
			//Then Baumgarte & Shapiro on p.390 say RPhi_ij is the same as p.57 substituting phi for ln(psi)
			// ... but I thought phi was ln(psi)?  Then why would you need to separate R_ij = RBar_ij + RPhi_ij ?  I thought the substitution showed that R_ij was RPhi_ij?
			cell.RPhi_ll = RBar_ll 
				- 2. * (
					DBar2_phi_ll 
					+ geomCell.gammaBar_ll * DBar2_phi
				) 
				+ 4. * (
					outer(cell.DBar_phi_l, cell.DBar_phi_l) 
					- geomCell.gammaBar_ll * normBar_DBar_phi
				);

			//R_ll(i,j) := R_ij = RPhi_ij + RBar_ij
			cell.R_ll = cell.RPhi_ll + RBar_ll;

			//R = gamma^ij R_ij
			cell.R = Tensor::dot(cell.gamma_uu, cell.R_ll);
			
			//partial2_psi_ll(i,j) := partial_i partial_j psi
			TensorSL partial2_psi_ll = partialSecondDerivative(auxGrid, &AuxCell::psi, auxGrid, &AuxCell::DBar_psi_l, dx, index);

			//DBar2_psi_ll(i,j) := DBar_i DBar_j ln(psi) = partial_i partial_j psi - connBar^k_ij partial_k psi
			TensorSL DBar2_psi_ll = partial2_psi_ll - cell.DBar_psi_l * cell.connBar_ull;

			//DBar2_psi := DBar^2 psi = gammaBar^ij DBar_i DBar_j psi
			cell.DBar2_psi = Tensor::dot(cell.gammaBar_uu, DBar2_psi_ll);
		});

		//calculate extrinsic curvature tensors
		parallel.foreach(range.begin(), range.end(), [&](DerefType index) {
			AuxCell & cell = auxGrid(index);
			GeomCell const & geomCell = geomGridRead(index);

			//ATilde_ul(i,j) := ATilde^i_j = gammaBar^ik ATilde_kj
			cell.ATilde_ul = cell.gammaBar_uu * geomCell.ATilde_ll;

			//ATilde_uu(i,j) := ATilde^ij = ATilde^i_k gammaBar^kj
			cell.ATilde_uu = cell.ATilde_ul * cell.gammaBar_uu;

			Real psiSquared = cell.psi * cell.psi;
			Real psiToTheFourth = psiSquared * psiSquared;
			Real psiToTheSixth = psiSquared * psiToTheFourth;

			//ABar_uu(i,j) := ABar^ij = 
			cell.ABar_uu = psiToTheSixth * cell.ATilde_uu;

			//A_ll(i,j) := ATilde_ij = exp(4phi) ATilde_ij = psi^4 ATilde_ij
			auto A_ll = geomCell.ATilde_ll * psiToTheFourth;
			
			//K_ll(i,j) := K_ij = A_ij + 1/3 gamma_ij K
			auto K_ll = A_ll + 1./3. * cell.gamma_ll * geomCell.K;
		
			//K^i_j := gamma^ik K_kj
			auto K_ul = cell.gamma_uu * K_ll;
		
			//tr_K_sq := tr(K^2) = (K^2)^i_i = K^ij K_ji = K^i_j K^j_i
			//this method uses tr(K^2) = K^ij K_ij in particular
			//tr_K_sq := tr(K^2) = K^ij K_ij
			// TODO rewrite this in terms of K and ATilde (or ABar)
			cell.tr_K_sq = K_ul.dot(K_ul.transpose());
		});

		//calcConstraints
		parallel.foreach(range.begin(), range.end(), [&](DerefType index) {
			AuxCell & cell = auxGrid(index);
			GeomCell const & geomCell = geomGridRead(index);
			MatterCell const & matterCell = matterGrid(index);

				//Hamiltonian constraint
			
			Real psiSquared = cell.psi * cell.psi;
			Real psiToTheFourth = psiSquared * psiSquared;
			Real psiToTheSixth = psiSquared * psiToTheFourth;

			//tr_ATilde_sq := tr(ATilde^2) = ATilde_ij ATilde^ji
			Real tr_ATilde_sq = cell.ATilde_ul.dot(cell.ATilde_ul.transpose());

			//Baumgarte & Shapiro p.390
			//H = gammaBar^ij DBar_i DBar_j exp(phi) - exp(phi)/8 RBar + exp(5phi)/8 ATilde_ij ATilde^ij - exp(5phi)/12 K^2 + 2 pi exp(5phi) rho
			//  = DBar^2 psi + psi( -RBar / 8 + psi^4 (tr(ATilde^2) / 8 - K^2 / 12 + 2 pi rho)
			cell.H = cell.DBar2_psi + cell.psi * (-cell.RBar / 8. + psiToTheFourth * (tr_ATilde_sq / 8. - geomCell.K * geomCell.K / 12. + 2. * M_PI * matterCell.rho));

				//momentum constraint
				
			//DBar_ABar_luu(i,j,k) := DBar_i ABar^jk
			TensorLSU DBar_ABar_luu = CovariantDerivativeClass_SU<
				Real, dim, AuxCell, AuxCell
			>::exec(auxGrid, &AuxCell::ABar_uu, auxGrid, &AuxCell::connBar_ull, dx, index);
	
			//DBar_ABar_u(i) := DBar_j ABar^ji
			TensorU DBar_ABar_u = DBar_ABar_luu.template contract<0,2>();
			
			//DBar_K_l(i) := DBar_i K
			TensorL DBar_K_l = Tensor::partialDerivative<partialDerivativeOrder, Real, dim, Real>(
				index, dx, [&](DerefType index) -> Real {
					index = Tensor::clamp(index, DerefType(), size-1);
					return geomGridRead(index).K;
				});

			//DBar_K_u(i) := DBar^i K = gammaBar^ij DBar_j K
			TensorU DBar_K_u = cell.gammaBar_uu * DBar_K_l;

			//Baumgarte & Shapiro 
			//p.65
			//M_u(i) := M^i = DBar_j ABar^ij - 2/3 psi^6 gammaBar^ij DBar_j K - 8 pi psi^10 S^i
			//p.390
			//M_u(i) := M^i = DBar_j (psi^6 ATilde^ij) - 2/3 psi^6 DBar^i K - 8 pi psi^6 S^i
			//		  = DBar_j ABar^ij + psi^6( - 2/3 gammaBar^ij DBar_j K - 8 pi S^i)
			//hmm, looks like that psi^10 turned into an exp(6phi) ... wonder why that isn't an exp(10phi) ...
			cell.M_u = DBar_ABar_u + psiToTheSixth * (-2./3. * DBar_K_u - 8. * M_PI * matterCell.S_u);
		});
	}


	virtual void getExplicitPartials(
		Real dt, 
		GeomGrid const & geomGridRead,	//read from this.  last iteration state.
		GeomGrid & partial_t_geomGrid)		//next iteration partials
	{
		calcAux(geomGridRead);
	
		Tensor::RangeObj<dim> range = auxGrid.range();
		parallel.foreach(range.begin(), range.end(), [&](DerefType index) {
			AuxCell & cell = auxGrid(index);
			GeomCell const & geomCell = geomGridRead(index);
			MatterCell const & matterCell = matterGrid(index);
			GeomCell & partial_t_geomCell = partial_t_geomGrid(index);

			//partial2_alpha_ll(i,j) := partial_i partial_j alpha
			TensorSL partial2_alpha_ll = partialSecondDerivative(geomGridRead, &GeomCell::alpha, auxGrid, &AuxCell::D_alpha_l, dx, index);

			//D2_alpha_ll(i,j) = D_i D_j alpha = partial_i partial_j alpha - conn^k_ij partial_k alpha
			TensorSL D2_alpha_ll = partial2_alpha_ll - cell.D_alpha_l * cell.conn_ull;

			//trace_partial_beta := partial_i beta^i
			Real trace_partial_beta = cell.partial_beta_lu.trace();

			//partial_t gammaBar_ij = -2 alpha ATilde_ij - 2/3 gammaBar_ij partial_k beta^k
			//		+ beta^k partial_k gammaBar_ij + gammaBar_ik partial_j beta^k 
			//		+ gammaBar_kj partial_i beta^k
			partial_t_geomCell.gammaBar_ll = 
				-2. * geomCell.alpha * geomCell.ATilde_ll 
				- 2./3. * geomCell.gammaBar_ll * trace_partial_beta
				+ geomCell.beta_u * cell.partial_gammaBar_lll
				+ 2. * makeSym(cell.partial_beta_lu * geomCell.gammaBar_ll);

			//S := S^i_i := gamma^ij S_ij
			Real S = matterCell.S_ll.dot(cell.gamma_uu);
			
			//partial_ATilde_lll(i,j,k) := partial_i ATilde_jk
			TensorLSL partial_ATilde_lll = Tensor::partialDerivative<partialDerivativeOrder, Real, dim, TensorSL>(
				index, dx, [&](DerefType index) -> TensorSL {
					index = Tensor::clamp(index, DerefType(), size-1);
					return geomGridRead(index).ATilde_ll;
				});

			//traceless portion of partial_t ATilde_ij := tracefree(-D^2 alpha + alpha (R_ij - 8 pi S_ij))
			TensorSL tracelessPortionOfPartialT_ATilde_ll = -D2_alpha_ll + geomCell.alpha * (cell.R_ll - 8. * M_PI * matterCell.S_ll);
			Real traceOfTraceless = cell.gamma_uu.dot(tracelessPortionOfPartialT_ATilde_ll);
			tracelessPortionOfPartialT_ATilde_ll -= 1./3. * cell.gamma_ll * traceOfTraceless;

			Real psiSquared = cell.psi * cell.psi;
			Real psiToTheFourth = psiSquared * psiSquared;
			//partial_t ATilde_ij = exp(-4phi) ((-D_i D_j alpha + alpha (R_ij - 8 pi S_ij))^TF + alpha (K ATilde_ij - 2 ATilde_il ATilde^l_j))
			//		+ beta^k partial_k ATilde_ij + ATilde_ik partial_j beta^k + ATilde_kj partial_i beta^k - 2/3 ATilde_ij partial_k beta^k
			partial_t_geomCell.ATilde_ll = 
				psiToTheFourth * tracelessPortionOfPartialT_ATilde_ll 
				+ geomCell.alpha * geomCell.K * geomCell.ATilde_ll
				- geomCell.alpha * 2. * geomCell.ATilde_ll * cell.ATilde_ul
				+ geomCell.beta_u * partial_ATilde_lll
				+ 2. * makeSym(cell.partial_beta_lu * geomCell.ATilde_ll)
				- 2./3. * geomCell.ATilde_ll * trace_partial_beta;

			//partial_phi_l(i) := partial_i phi
			TensorL partial_phi_l = Tensor::partialDerivative<partialDerivativeOrder, Real, dim, Real>(
				index, dx, [&](DerefType index) -> Real {
					index = Tensor::clamp(index, DerefType(), size-1);
					return geomGridRead(index).phi;
				});

			//partial_t -alpha K / 6 + beta^i partial_i phi + partial_i beta^i / 6
			partial_t_geomCell.phi = (trace_partial_beta - geomCell.alpha * geomCell.K) / 6. + geomCell.beta_u * partial_phi_l;

			//D_K_l(i) := D_i K
			TensorL D_K_l = CovariantDerivativeClass_0<
				Real, dim, GeomCell, AuxCell
			>::exec(geomGridRead, &GeomCell::K, auxGrid, &AuxCell::conn_ull, dx, index);

			//partial_t K = -gamma^ij D_i D_j alpha + alpha(K_ij K^ij + 4 pi (rho + S)) + beta^i D_i K
			partial_t_geomCell.K = 
				- cell.gamma_uu.dot(D2_alpha_ll)
				+ geomCell.beta_u * D_K_l
				+ geomCell.alpha * (cell.tr_K_sq + 4. * M_PI * (matterCell.rho + S));

			//partial2_beta_llu(i,j,k) := partial_i partial_j beta^k
			TensorSLU partial2_beta_llu = partialSecondDerivative(geomGridRead, &GeomCell::beta_u, auxGrid, &AuxCell::partial_beta_lu, dx, index);

			//partial_K_l(i) := partial_i K
			TensorL partial_K_l = Tensor::partialDerivative<partialDerivativeOrder, Real, dim, Real>(
				index, dx, [&](DerefType index) -> Real {
					index = Tensor::clamp(index, DerefType(), size-1);
					return geomGridRead(index).K;
				});
			
			//connBar^i is the connection function / connection coefficient iteration with Hamiltonian constraint baked in (Baumgarte & Shapiro p.389, Alcubierre p.86).
			//partial_t connBar^i = -2 ATilde^ij partial_j alpha + 2 alpha (connBar^i_jk ATilde^kj - 2/3 gammaBar^ij partial_j K - 8 pi gammaBar^ij S_j + 6 ATilde^ij partial_j phi)
			//	+ beta^j partial_j connBar^i - connBar^j partial_j beta^i + 2/3 connBar^i partial_j beta^j + 1/3 gammaBar^li partial_l partial_j beta^j + gammaBar^lj partial_j partial_l beta^i
			partial_t_geomCell.connBar_u = 
				2./3. * geomCell.connBar_u * trace_partial_beta
				- 2. * geomCell.ATilde_ll * cell.D_alpha_l
				+ 2. * geomCell.alpha * (
					-2. / 3. * cell.gammaBar_uu * partial_K_l
					- 8. * M_PI * geomCell.gammaBar_ll * matterCell.S_u
					+ 6. * cell.ATilde_uu * partial_phi_l
				)
				+ geomCell.beta_u * cell.partial_connBar_lu
				- geomCell.connBar_u * cell.partial_beta_lu
				+ 2. * geomCell.alpha * cell.connBar_ull.template interior<2>(cell.ATilde_uu)
				+ 1./3. * cell.gammaBar_uu * partial2_beta_llu.template contract<1,2>()
				+ cell.gammaBar_uu.template interior<2>(partial2_beta_llu);
			;
		});
	}

	void update(Real dt) {
		integrator->update(dt);

		constrain(*geomGridWriteCurrent);

		time += dt;

		//TODO update fluid components

		//and swap source and dest grids
		//do something more clever if we ever get any more than 2 histories
		//if you never need more than 2 then maybe we won't need histories at all, just partials?
		std::swap(geomGridReadCurrent, geomGridWriteCurrent);
	}

	//apply constraints to determinants and traces
	void constrain(GeomGrid &geomGrid) {

		Tensor::RangeObj<dim> range = auxGrid.range();
		parallel.foreach(range.begin(), range.end(), [&](DerefType index) {
			GeomCell &geomCell = geomGrid(index);

			TensorSL &gammaBar_ll = geomCell.gammaBar_ll;
			TensorSL &ATilde_ll = geomCell.ATilde_ll;
			
			/*
			det(gammaBar_ij) 
			= det(gamma^-1/3 gamma_ij)
			= gamma^-1 gamma
			= 1
			*/
			Real gammaBar = determinant(gammaBar_ll);
			Real oneOverCubeRootGammaBar = 1. / cbrt(gammaBar);
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					gammaBar_ll(i,j) *= oneOverCubeRootGammaBar;
				}
			}
		
			TensorSU gammaBar_uu = Tensor::inverse(gammaBar_ll, (Real)1); 

			/*
			tr(A_ij)
			= tr(K_ij - 1/3 gamma_ij K)
			= gamma^ij K_ij - 1/3 gamma^ij gamma_ij K
			= K - 1/3 3 K
			= 0

			tr(ATilde_ij) = 3 psi^-4 tr(A_ij) = 3 psi^-4 * 0 
			= 0
			*/
			Real tr_ATilde = 0;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					tr_ATilde += gammaBar_uu(i,j) * ATilde_ll(i,j);
				}
			}

			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					ATilde_ll(i,j) -= 1./3. * gammaBar_ll(i,j) * tr_ATilde;
				}
			}
		});
	}
};
