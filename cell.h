#pragma once

#include "vector.h"
#include "tensor.h"
#include "invert.h"

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
		ln_sqrt_gamma(real()),
		ln_psi(real()),
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

	//conformal factor 
	//currently derived from the iterated ln(sqrt(gamma))	
	real psi;

	//log of conformal factor
	//stored separately because it's used often enough
	//psi = gamma^(1/12)
	//ln psi = 1/12 ln gamma = 1/6 ln(sqrt(gamma))
	real ln_psi;

	//gamma_uu(i,j) := gamma^ij = inverse(gamma_ij) = covalent(gamma_ij) / det(gamma_ij)
	tensor_su gamma_uu;

	//gammaBar_ll(i,j) := gammaBar_ij = psi^-4 gamma_ij
	tensor_sl gammaBar_ll;

	//gammaBar_uu(i,j) := gammaBar^ij = psi^4 gamma^ij
	tensor_su gammaBar_uu;

	//conn_lll(i,j,k) := conn_ijk
	tensor_lsl  conn_lll;

	//conn_ull(i,j,k) := conn^i_jk
	tensor_usl conn_ull;

	//R_ll(i,j) := R_ij
	tensor_sl R_ll;

	//K_ul(i,j) := K^i_j
	tensor_ul K_ul;

	//K_uu(i,j) := K^ij
	tensor_su K_uu;

	//tr_K_sq := (K^2)^i_i = K^ij K_ji
	real tr_K_sq;


	// calculations of aux values


	//used during init
	//ln_sqrt_gamma := ln(sqrt(det(gamma_ij)))
	void calcLnSqrtGammaFromGammaLL() {
		//gamma = det(gamma_ij)
		real gamma = determinant(gamma_ll);

		//ln_sqrt_gamma := ln(sqrt(gamma))
		ln_sqrt_gamma = .5 * log(gamma);
	}

	//ln_psi := ln(psi) = 1/6 ln(sqrt(gamma))
	//psi = exp(ln(psi))
	void calcPsiFromLnSqrtGamma() {
		//ln(psi) = 1/6 ln(sqrt(gamma))
		ln_psi = ln_sqrt_gamma / 6.;

		//psi = exp(ln(psi))
		psi = exp(ln_psi);
	}

	//gammaBar_ij = psi^-4 gamma_ij
	//gammaBar^ij = inverse(gammaBar_ij)
	//gamma^ij = psi^-4 gammaBar^ij
	void calcGammaBar() {
		real psiSquared = psi * psi;
		real psiToTheFourth = psiSquared * psiSquared;
		real oneOverPsiToTheFourth = 1. / psiToTheFourth;

		//gammaBar_ij = psi^-4 gamma_ij
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				gammaBar_ll(i,j) = oneOverPsiToTheFourth * gamma_ll(i,j);
			}
		}

		//gammaBar^ij = inverse(gammaBar_ij)
		gammaBar_uu = inverse(gammaBar_ll, 1.);

		//gamma^ij = psi^-4 gammaBar^ij
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j <= i; ++j) {
				gamma_uu(i,j) = oneOverPsiToTheFourth * gammaBar_uu(i,j);
			}
		}
	}


	//	operators for ease of use in algorithms


	Cell operator*(const real &b) const {
		const Cell &a = *this;
		Cell c;
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
		Cell c;
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
};

