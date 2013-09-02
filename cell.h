#pragma once

#include "vector.h"
#include "tensor.h"
#include "inverse.h"

/*
structure solely of integrated values
"GeomCell" is an incorrect name
 if I ever write a huge implicit solver and integrate matter terms alongside geometridynamic terms, they would go here

N-D case has this many elements:
	1 + N + N(N+1)/2 + N(N+1)/2 + 1 + 1
	3 + 2N + N^2
1D: 6 elements
2D: 11 elements
3D: 18 elements
*/
template<typename real_, int dim_>
struct GeomCell {
	typedef real_ real;
	enum { dim = dim_ };

	typedef tensor<real, upper<dim>> tensor_u;
	typedef tensor<real, symmetric<lower<dim>, lower<dim>>> tensor_sl;

	//lapse
	real alpha;

	//shift
	//beta_u(i) := beta^i
	//beta^t = 0
	tensor_u beta_u;

	//spatial metric
	//gamma_ll(i,j) := gamma_ij
	//gamma_it = gamma_tj = 0
	tensor_sl gamma_ll;

	//extrinsic curvature
	//K_ll(i,j) := K_ij
	//K_it = K_tj = 0
	tensor_sl K_ll;

	//extrinsic curvature trace
	//K := K^i_k
	real K;

	//related to the conformal factor of metric
	//ln_sqrt_gamma := ln(sqrt(det(gamma_ij)))
	real ln_sqrt_gamma;

	GeomCell()
	:	alpha(real()),
		K(real()),
		ln_sqrt_gamma(real())
	{}

	//used during init
	//gamma = det(gamma_ij)
	//ln_sqrt_gamma := ln(sqrt(gamma))
	void calc_ln_sqrt_gamma_from_gamma_ll() {
		//gamma = det(gamma_ij)
		real gamma = determinant(gamma_ll);

		//GeomCell::ln_sqrt_gamma := ln(sqrt(gamma))
		ln_sqrt_gamma = .5 * log(gamma);
	}

	//operators used with integration

	GeomCell operator*(const real &scalar) const {
		GeomCell result;
		result.alpha = alpha * scalar;
		result.beta_u = beta_u * scalar;
		result.gamma_ll = gamma_ll * scalar;
		result.K_ll = K_ll * scalar;
		result.K = K * scalar;
		result.ln_sqrt_gamma = ln_sqrt_gamma * scalar;
		return result;
	}

	GeomCell &operator+=(const GeomCell &sourceCell) {
		alpha += sourceCell.alpha;
		beta_u += sourceCell.beta_u;
		gamma_ll += sourceCell.gamma_ll;
		K_ll += sourceCell.K_ll;
		K += sourceCell.K;
		ln_sqrt_gamma += sourceCell.ln_sqrt_gamma;
		return *this;
	}
};

template<typename real_, int dim_>
struct MatterCell {
	typedef real_ real;
	enum { dim = dim_ };

	typedef tensor<real, upper<dim>> tensor_u;
	typedef tensor<real, symmetric<lower<dim>, lower<dim>>> tensor_sl;

	MatterCell()
	:	rho(real())
	{}

	//energy density
	//rho = n_a n_b T^ab
	real rho;

	//momentum
	//S_u(i) := S^i = -gamma^ij * n^a T_aj
	tensor_u S_u;

	//spatial stress energy
	//S_ll(i,j) := S_ij = gamma_ic gamma_jd T^cd
	//					= gamma_i^c gamma_j^d T_cd
	tensor_sl S_ll;
};

template<typename real_, int dim_>
struct AuxCell {
	typedef real_ real;
	enum { dim = dim_ };

	/*
	different ranked tensor types
	notation:
	tensor_[ulsa]*
	'u' means upper
	'l' means lower
	's' means symmetric
	'a' means antisymmetric
	*/
	typedef ::tensor<real,upper<dim>> tensor_u;
	typedef ::tensor<real,lower<dim>> tensor_l;
	typedef ::tensor<real,lower<dim>,lower<dim>> tensor_ll;
	typedef ::tensor<real,upper<dim>,lower<dim>> tensor_ul;
	typedef ::tensor<real,symmetric<upper<dim>,upper<dim>>> tensor_su;
	typedef ::tensor<real,symmetric<lower<dim>,lower<dim>>> tensor_sl;
	typedef ::tensor<real,upper<dim>,symmetric<lower<dim>,lower<dim>>> tensor_usl;
	typedef ::tensor<real,lower<dim>,symmetric<lower<dim>,lower<dim>>> tensor_lsl;


	//our tensors initialze to zero, so why not our reals too?
	AuxCell() 
	:	H(real()),
		gamma(real()),
		R(real()),
		tr_K_sq(real()),
		psi(real()),
		ln_psi(real()),
		DBar2_psi(real()),
		RBar(real()),
		tr_ABar_sq(real())
	{}


	//	constraint variables: should always be zero


	//hamiltonian constraint
	real H;

	//momentum constraint
	tensor_u M_u;


	//	aux values computed and stored for partials


	//D_alpha_l(i) := D_i alpha
	tensor_l D_alpha_l;

	//beta_l(i) := beta_i
	//beta_t = beta^k beta_k
	tensor_l beta_l;

	//partial_gamma_lll(k,i,j) := partial_k gamma_ij
	tensor_lsl partial_gamma_lll;

	//gamma_uu(i,j) := gamma^ij = inverse(gamma_ij) = covalent(gamma_ij) / det(gamma_ij)
	tensor_su gamma_uu;

	//gamma = det(gamma_ij)
	real gamma;

	//conn_lll(i,j,k) := conn_ijk
	tensor_lsl  conn_lll;

	//conn_ull(i,j,k) := conn^i_jk
	tensor_usl conn_ull;

	//R_ll(i,j) := R_ij
	tensor_sl R_ll;

	//Gaussian (scalar) curvature
	//R = R^i_i
	real R;

		//curvature aux variables

	//K_ul(i,j) := K^i_j
	tensor_ul K_ul;

	//K_uu(i,j) := K^ij
	tensor_su K_uu;

	//tr_K_sq := (K^2)^i_i = K^ij K_ji
	real tr_K_sq;

		//conformal factor aux variables

	//conformal factor 
	//currently derived from the iterated ln(sqrt(gamma))	
	real psi;

	//log of conformal factor
	//stored separately because it's used often enough
	//psi = gamma^(1/12)
	//ln psi = 1/12 ln gamma = 1/6 ln(sqrt(gamma))
	real ln_psi;

	//DBar_psi_l(i) := DBar_i psi
	tensor_l DBar_psi_l;

	//DBar_ln_psi_l(i) := DBar_i ln(psi)
	tensor_l DBar_ln_psi_l;

	//DBar2_psi := gammaBar^ij DBar_i DBar_j psi
	real DBar2_psi;
	
	//gammaBar_ll(i,j) := gammaBar_ij = psi^-4 gamma_ij
	tensor_sl gammaBar_ll;

	//gammaBar_uu(i,j) := gammaBar^ij = psi^4 gamma^ij
	tensor_su gammaBar_uu;

	//partial_gammaBar_lll(i,j,k) := partial_i gammaBar_jk
	tensor_lsl partial_gammaBar_lll;

	//connBar_lll(i,j,k) := connBar_ijk = 1/2 (partial_k gammaBar_ij + partial_j gammaBar_ik - partial_i gammaBar_jk)
	tensor_lsl connBar_lll;

	//connBar_ull(i,j,k) := connBar^i_jk = gammaBar^il connBar_ljk
	tensor_usl connBar_ull;
	
	//RBar_ll(i,j) := RBar_ij
	tensor_sl RBar_ll;
	
	//RBar = gammaBar^ij RBar_ij
	real RBar;

	//traceless part of extrinsic curvature tensor
	//A_ll(i,j) := A_ij = K_ij - 1/3 gamma_ij K
	tensor_sl A_ll;

	//ABar_ll(i,j) := ABar_ij = psi^2 A_ij
	tensor_sl ABar_ll;
	
	//ABar_uu(i,j) := ABar^ij = psi^10 A^ij
	tensor_su ABar_uu;

	//tr_ABar_sq := tr(ABar^2) = ABar_ij ABar^ij
	real tr_ABar_sq;

		// helper functions

	void calc_psi_and_ln_psi_from_ln_sqrt_gamma(const GeomCell<real, dim> &geomCell) {
		//ln_psi := ln(psi) = 1/6 ln(sqrt(gamma))
		//ln(psi) = 1/6 ln(sqrt(gamma))
		ln_psi = geomCell.ln_sqrt_gamma / 6.;
		
		//psi = exp(ln(psi))
		psi = ln_psi;
	}
};

