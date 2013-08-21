#pragma once

#include "vec.h"
#include "symmat.h"

template<int dim_, typename real_>
struct Cell {
	typedef real_ real;
	enum { dim = dim_ };

	typedef ::vec<dim,real> vec;
	typedef ::symmat<dim,real> symmat;

	//	state variables

	//lapse
	real alpha;

	//shift
	//beta^i
	//beta^t = 0
	vec beta_u;

	//spatial metric
	//gamma_ij
	//gamma_it = gamma_tj = 0
	symmat gamma_ll;
	
	//extrinsic curvature
	//K_ij
	//K_it = K_tj = 0
	symmat K_ll;

	// aux values computed and stored for partials

	//beta_i
	//beta_t = beta^k beta_k
	vec beta_l;

	//gamma^ij = inverse of gamma_ij
	symmat gamma_uu;

	//conn^i_jk
	::vec<dim, symmat> conn_ull;
};

