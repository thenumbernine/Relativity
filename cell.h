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
	typedef ::tensor<real,symmetric<upper<dim>,upper<dim>>> tensor_su;
	typedef ::tensor<real,symmetric<lower<dim>,lower<dim>>> tensor_sl;
	typedef ::tensor<real,upper<dim>,symmetric<lower<dim>,lower<dim>>> tensor_usl;
	typedef ::tensor<real,lower<dim>,symmetric<lower<dim>,lower<dim>>> tensor_lsl;


	//	state variables


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

	//energy density
	//rho = n_a n_b T^ab
	real rho;

	//spatial stress energy
	//S_ll(i,j) := S_ij = gamma_i^c gamma_j^d T_cd (for stress-energy tensor T_ab)
	tensor_ll S_ll;


	// aux values computed and stored for partials


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
};

