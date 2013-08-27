#pragma once

#include "vector.h"
#include "tensor.h"

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
	: alpha(0), rho(0), H(0), K(0)
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


	//	stress-energy variables


	//energy density
	//rho = n_a n_b T^ab
	real rho;

	//spatial stress energy
	//S_ll(i,j) := S_ij = gamma_i^c gamma_j^d T_cd (for stress-energy tensor T_ab)
	tensor_ll S_ll;


	//	constraint variables: should always be zero


	//hamiltonian constraint
	real H;


	//	aux values computed and stored for partials


	//D_alpha_l(i) := D_i alpha
	tensor_l D_alpha_l;

	//beta_l(i) := beta_i
	//beta_t = beta^k beta_k
	tensor_l beta_l;

	//partial_gamma_lll(k,i,j) := partial_k gamma_ij
	tensor_lsl partial_gamma_lll;

	//gamma_uu(i,j) := gamma^ij = inverse of gamma_ij
	tensor_su gamma_uu;
	
	//conn_lll(i,j,k) := conn_ijk
	tensor_lsl  conn_lll;

	//conn_ull(i,j,k) := conn^i_jk
	tensor_usl conn_ull;

	//R_ll(i,j) := R_ij
	tensor_sl R_ll;

	//K_ul(i,j) := K^i_j
	tensor_ul K_ul;

	//K := K^i_k
	real K;


	//	operators for ease of use in algorithms


	Cell operator*(const real &b) const {
		const Cell &a = *this;
		Cell c;
		c.alpha = a.alpha * b;
		c.beta_u = a.beta_u * b;
		c.gamma_ll = a.gamma_ll * b;
		c.K_ll = a.K_ll * b;
		//these should be maintained - and are probably a separate structure for that reason? 
		c.rho = a.rho * b;
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
		c.K_ll = a.K_ll + a.K_ll;
		//these should be maintained
		c.rho = a.rho;
		c.S_ll = a.S_ll;
		//the rest are computed mid-iteration
		return c;
	}
};

