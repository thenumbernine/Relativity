#pragma once

#include "vector.h"
#include "tensor.h"

template<int dim_, typename real_>
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
	typedef ::tensor<real,symmetric<upper<dim>,upper<dim>>> tensor_su;
	typedef ::tensor<real,symmetric<lower<dim>,lower<dim>>> tensor_sl;
	typedef ::tensor<real,upper<dim>,symmetric<lower<dim>,lower<dim>>> tensor_usl;
	typedef ::tensor<real,lower<dim>,symmetric<lower<dim>,lower<dim>>> tensor_lsl;

	//	state variables

	//lapse
	real alpha;

	//shift
	//beta^i
	//beta^t = 0
	tensor_u beta_u;

	//spatial metric
	//gamma_ij
	//gamma_it = gamma_tj = 0
	tensor_sl gamma_ll;
	
	//extrinsic curvature
	//K_ij
	//K_it = K_tj = 0
	tensor_sl K_ll;

	// aux values computed and stored for partials

	//beta_i
	//beta_t = beta^k beta_k
	tensor_l beta_l;

	//partial_k gamma_ij
	tensor_lsl partial_gamma_lll;

	//gamma^ij = inverse of gamma_ij
	tensor_su gamma_uu;
	
	//conn_ijk
	tensor_lsl  conn_lll;

	//conn^i_jk
	tensor_usl conn_ull;

	//R_ij
	tensor_sl ricci_ll;
};

