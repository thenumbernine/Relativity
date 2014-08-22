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

looks like a dense system could get big for submatrixes (18^18 = 324 elements for 3D)
	so if I do make it to implicit I might do a single sparse matrix for it all
	...or a sparse block matrix of sparse matrixes
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

	//log of the conformal factor (to the arbitrary root ...)
	//exp(phi) = psi
	//psi^12 = gamma
	real phi;

	//conformal spatial metric
	//gammaBar_ll(i,j) := gammaBar_ij = exp(-4phi) gamma_ij
	// such that det(gammaBar_ij) = 1
	tensor_sl gammaBar_ll;

	//minimal distortion conformal extrinsic curvature
	//ATilde_ll(i,j) := ATilde_ij = psi^-4 A_ij = exp(-4phi) A_ij = exp(-4phi) (K_ij - 1/3 gamma_ij K) = exp(-4phi) K_ij - 1/3 gammaBar_ij K
	tensor_sl ATilde_ll;
	
	//extrinsic curvature trace
	//K := K^i_i
	real K;

	//connection function
	//connBar_u(i) := connBar^i = gammaBar^jk connBar^i_jk = -partial_j gammaBar^ij
	tensor_u connBar_u;

	GeomCell()
	:	alpha(real()),
		phi(real()),
		K(real())
	{}

	//operators used with integration

	GeomCell operator*(const real &scalar) const {
		GeomCell result;
		result.alpha = alpha * scalar;
		result.beta_u = beta_u * scalar;
		result.gammaBar_ll = gammaBar_ll * scalar;
		result.ATilde_ll = ATilde_ll * scalar;
		result.K = K * scalar;
		result.phi = phi * scalar;
		result.connBar_u = connBar_u * scalar;
		return result;
	}

	GeomCell &operator+=(const GeomCell &sourceCell) {
		alpha += sourceCell.alpha;
		beta_u += sourceCell.beta_u;
		gammaBar_ll += sourceCell.gammaBar_ll;
		ATilde_ll += sourceCell.ATilde_ll;
		K += sourceCell.K;
		phi += sourceCell.phi;
		connBar_u += sourceCell.connBar_u;
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
	'a' means antisymmetric (haven't needed this one yet)
	*/
	typedef ::lower<dim> lower;
	typedef ::upper<dim> upper;
	
	typedef ::symmetric<lower,lower> symmetric_lower;
	typedef ::symmetric<upper,upper> symmetric_upper;
	
	typedef ::tensor<real,upper> tensor_u;
	typedef ::tensor<real,lower> tensor_l;
	typedef ::tensor<real,lower,lower> tensor_ll;
	typedef ::tensor<real,upper,lower> tensor_ul;
	typedef ::tensor<real,lower,upper> tensor_lu;
	typedef ::tensor<real,symmetric_upper> tensor_su;
	typedef ::tensor<real,symmetric_lower> tensor_sl;
	typedef ::tensor<real,upper,symmetric_lower> tensor_usl;
	typedef ::tensor<real,lower,symmetric_lower> tensor_lsl;


	//our tensors initialze to zero, so why not our reals too?
	AuxCell() 
	:	H(real()),
		gamma(real()),
		R(real()),
		psi(real()),
		DBar2_psi(real()),
		RBar(real()),
		tr_K_sq(real())
	{}


		//constraint variables: should always be zero


	//hamiltonian constraint
	real H;

	//momentum constraint
	tensor_u M_u;


		//lapse and shift related


	//D_alpha_l(i) := D_i alpha
	tensor_l D_alpha_l;

	//beta_l(i) := beta_i
	//beta_t = beta^k beta_k
	tensor_l beta_l;

	//partial_beta_lu(i,j) := partial_i beta^j
	tensor_lu partial_beta_lu;

		//metric related

	//spatial metric
	//gamma_ll(i,j) := gamma_ij
	//gamma_it = gamma_tj = 0
	tensor_sl gamma_ll;

	//gamma_uu(i,j) := gamma^ij = inverse(gamma_ij) = covalent(gamma_ij) / det(gamma_ij)
	tensor_su gamma_uu;

	//gamma = det(gamma_ij)
	real gamma;

	//conn_lll(i,j,k) := conn_ijk
	tensor_lsl  conn_lll;

	//conn_ull(i,j,k) := conn^i_jk
	tensor_usl conn_ull;

	//RPhi_ll(i,j) := RPhi_ij is found by substituting phi for ln(psi) in eqn 3.10 of Baumgarte & Shapiro.
	// (but I thought this equation was for calculating R_ij?  And that ln(psi) = phi to begin with? So why do we now separate R_ij = RPhi_ij + RBar_ij?)
	tensor_sl RPhi_ll;

	//R_ll(i,j) := R_ij = RBar_ij + RPhi_ij
	tensor_sl R_ll;

	//Gaussian (scalar) curvature
	//R = R^i_i
	real R;

		//conformal factor / metric

	//conformal factor 
	//currently derived from the iterated ln(sqrt(gamma))	
	real psi;

	//DBar_psi_l(i) := DBar_i psi
	tensor_l DBar_psi_l;

	//DBar_phi_l(i) := DBar_i ln(psi)
	tensor_l DBar_phi_l;

	//DBar2_psi := gammaBar^ij DBar_i DBar_j psi
	real DBar2_psi;
	
	//gammaBar_uu(i,j) := gammaBar^ij = psi^4 gamma^ij
	tensor_su gammaBar_uu;

	//partial_gammaBar_lll(i,j,k) := partial_i gammaBar_jk
	tensor_lsl partial_gammaBar_lll;

	//partial2_gammaBar_ll(i,j) := gammaBar^lm partial_l partial_m gammaBar_ij
	// should the name be partialBar2? since I use gammaBar to raise it?
	tensor_sl partial2_gammaBar_ll;

	//connBar_lll(i,j,k) := connBar_ijk = 1/2 (partial_k gammaBar_ij + partial_j gammaBar_ik - partial_i gammaBar_jk)
	tensor_lsl connBar_lll;

	//connBar_ull(i,j,k) := connBar^i_jk = gammaBar^il connBar_ljk
	tensor_usl connBar_ull;

	//partial_connBar_lu(i,j) := partial_i connBar^j
	tensor_lu partial_connBar_lu;

	//RBar_ll(i,j) := RBar_ij
	tensor_sl RBar_ll;

	//RBar = gammaBar^ij RBar_ij
	real RBar;
		
		//extrinsic curvature

	//tr_K_sq := (K^2)^i_i = K^ij K_ji
	real tr_K_sq;

		//conformal extrinsic curvature

	//ABar_uu(i,j) := ABar^ij = psi^10 A^ij
	tensor_su ABar_uu;

	//ATilde_ul(i,j) := ATilde^i_j = gammaBar^ik ATilde_kj
	tensor_ul ATilde_ul;

	//ATilde_uu(i,j) := ATilde^ij = ATilde^i_k gammaBar^kj
	tensor_su ATilde_uu;
};

