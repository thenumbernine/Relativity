#pragma once

#include "vector.h"
#include "oneform.h"
#include "symmat.h"

template<int dim_, typename real_>
struct Cell {
	typedef real_ real;
	enum { dim = dim_ };

	typedef ::vector<dim,real> vector;
	typedef ::oneform<dim,real> oneform;
	typedef ::symmat<dim,real> symmat;

	//	state variables

	//lapse
	real alpha;

	//shift
	//beta^i
	//beta^t = 0
	vector beta_u;

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
	oneform beta_l;

	//partial_k gamma_ij
	::oneform<dim, symmat> partial_gamma_lll;

	//gamma^ij = inverse of gamma_ij
	symmat gamma_uu;
	
	//conn_ijk
	::oneform<dim, symmat> conn_lll;

	//conn^i_jk
	::vector<dim, symmat> conn_ull;

	//R_ij
	symmat ricci_ll;
};

