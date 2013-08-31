#pragma once

#include "vector.h"
#include "tensor.h"
#include "inverse.h"

template<typename real_, int dim_>
struct Cell {
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
	Cell() 
	: 	alpha(real()),
		rho(real()),
		H(real()),
		K(real()),
		gamma(real()),
		R(real()),
		ln_sqrt_gamma(real()),
		ln_psi(real()),
		DBar2_psi(real()),
		RBar(real()),
		tr_K_sq(real())
	{}


	//	geometridynmaic variables


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


		// extra formalism variables that are iterated


	//extrinsic curvature trace
	//K := K^i_k
	real K;

	//related to the conformal factor of metric
	//ln_sqrt_gamma := ln(sqrt(det(gamma_ij)))
	real ln_sqrt_gamma;


	//	stress-energy-derived variables


	//energy density
	//rho = n_a n_b T^ab
	real rho;

	//momentum
	//S_u(i) := S^i = -gamma^ij * n^a T_aj
	tensor_u S_u;

	//spatial stress energy
	//S_ll(i,j) := S_ij = gamma_ic gamma_jd T^cd
	//					= gamma_i^c gamma_j^d T_cd
	tensor_ll S_ll;


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


	// calculations of aux values


	//used during init
	//gamma = det(gamma_ij)
	//ln_sqrt_gamma := ln(sqrt(gamma))
	void calc_ln_sqrt_gamma_from_gamma_ll() {
		//gamma = det(gamma_ij)
		gamma = determinant(gamma_ll);

		//ln_sqrt_gamma := ln(sqrt(gamma))
		ln_sqrt_gamma = .5 * log(gamma);
	}

	//ln_psi := ln(psi) = 1/6 ln(sqrt(gamma))
	//psi = exp(ln(psi))
	void calc_psi_from_ln_sqrt_gamma() {
		//ln(psi) = 1/6 ln(sqrt(gamma))
		ln_psi = ln_sqrt_gamma / 6.;

		//psi = exp(ln(psi))
		psi = exp(ln_psi);
	}

	//option-1 method
	//gamma = det(gamma_ij)
	//gamma^ij = ((gamma_kl)^-1)^ij
	void calc_gamma_uu_from_gamma_ll() {
		gamma = determinant(gamma_ll);
		gamma_uu = inverse(gamma_ll, gamma);
	}

	//option-2 method
	//gammaBar_ij = psi^-4 gamma_ij
	//gammaBar^ij = inverse(gammaBar_ij)
	//gamma = psi^12
	void calc_gammaBar_uu_and_gammaBar_ll_from_psi() {
		real psiSquared = psi * psi;
		real psiToTheFourth = psiSquared * psiSquared;
		real oneOverPsiToTheFourth = 1. / psiToTheFourth;

		//either this or another exp() call
		real psiToTheEighth = psiToTheFourth * psiToTheFourth;
		gamma = psiToTheFourth * psiToTheEighth;

		//gammaBar_ij = psi^-4 gamma_ij
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				gammaBar_ll(i,j) = oneOverPsiToTheFourth * gamma_ll(i,j);
			}
		}

		//gammaBar^ij = inverse(gammaBar_ij)
		gammaBar_uu = inverse(gammaBar_ll, 1.);
	}

	//option-2 method
	//gamma^ij = psi^-4 gammaBar^ij
	void calc_gamma_uu_from_gammaBar_uu_and_psi() {
		real psiSquared = psi * psi;
		real psiToTheFourth = psiSquared * psiSquared;
		real oneOverPsiToTheFourth = 1. / psiToTheFourth;

		//gamma^ij = psi^-4 gammaBar^ij
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				gamma_uu(i,j) = oneOverPsiToTheFourth * gammaBar_uu(i,j);
			}
		}
	}

	//K^i_j := gamma^ik K_kj
	void calc_K_ul() {
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j < dim; ++j) {
				K_ul(i,j) = 0;
				for (int k = 0; k < dim; ++k) {
					K_ul(i,j) += gamma_uu(i,k) * K_ll(k,j);
				}
			}
		}
	}

	//K = K^i_i
	void calc_K() {
		K = 0.;
		for (int i = 0; i < dim; ++i) {
			K += K_ul(i,i);
		}
	}
	
	//K^ij = K^i_k gamma^kj
	void calc_K_uu() {
		//K_uu(i,j) := K^ij = K^i_k gamma^kj
		for (int i = 0; i  < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				K_uu(i,j) = 0;
				for (int k = 0; k < dim; ++k) {
					K_uu(i,j) += K_ul(i,k) * gamma_uu(k,j);
				}
			}
		}
	}
	
	//tr_K_sq := tr(K^2) = (K^2)^i_i = K^ij K_ji = K^i_j K^j_i
	//this method uses tr(K^2) = K^ij K_ij in particular
	void calc_tr_K_sq() {
		tr_K_sq = 0.;
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j < dim; ++j) {
				tr_K_sq += K_uu(i,j) * K_ll(i,j); 
			}
		}
	}


	//	operators for ease of use in algorithms


	Cell operator*(const real &b) const {
		const Cell &a = *this;
		Cell c = a;
		c.alpha = a.alpha * b;
		c.beta_u = a.beta_u * b;
		c.gamma_ll = a.gamma_ll * b;
		c.K_ll = a.K_ll * b;
		c.K = a.K * b;
		c.ln_sqrt_gamma = a.ln_sqrt_gamma * b;
		//these should be maintained - and are probably a separate structure for that reason? 
		c.rho = a.rho * b;
		c.S_u = a.S_u * b;
		c.S_ll = a.S_ll * b;
		//the rest are computed mid-iteration.  once again, separate structure?
		return c;
	}

	Cell operator+(const Cell &b) const {
		const Cell &a = *this;
		Cell c = a;
		c.alpha = a.alpha + b.alpha;
		c.beta_u = a.beta_u + b.beta_u;
		c.gamma_ll = a.gamma_ll + b.gamma_ll;
		c.K_ll = a.K_ll + b.K_ll;
		c.K = a.K + b.K;
		c.ln_sqrt_gamma = a.ln_sqrt_gamma + b.ln_sqrt_gamma;
		//these should be maintained
		c.rho = a.rho;
		c.S_u = a.S_u;
		c.S_ll = a.S_ll;
		//the rest are computed mid-iteration
		return c;
	}

	//I want to preserve the aux variables for debug output
	//but I don't want to calculate them once prior to getting partials (as integrators do many times per step)
	// and calculate them again just to output debug info
	//solutions?
	// - dirty bits and two calcAux() calls
	// - copying stuff over here
	// - separate out cells into ADM, ADMAux, StressEnergy, etc
	//I'll do the 2nd for now...
	Cell &operator+=(const Cell &a) {
		real old_alpha = alpha;
		tensor_u old_beta_u = beta_u;
		tensor_sl old_gamma_ll = gamma_ll;
		tensor_sl old_K_ll = K_ll;
		real old_K = K;
		real old_ln_sqrt_gamma = ln_sqrt_gamma;

		//copy *all* over
		(*this) = a;

		//and only update these fields
		alpha = old_alpha;
		beta_u += old_beta_u;
		gamma_ll += old_gamma_ll;
		K_ll += old_K_ll;
		K += old_K;
		ln_sqrt_gamma += old_ln_sqrt_gamma;
		
		return *this;
	}
};

