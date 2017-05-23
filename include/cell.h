#pragma once

#include "Tensor/Vector.h"
#include "Tensor/Tensor.h"

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
template<typename Real_, int dim_>
struct GeomCell {
	typedef Real_ Real;
	enum { dim = dim_ };

	typedef Tensor::Tensor<Real, Tensor::Upper<dim>> TensorU;
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<dim>, Tensor::Lower<dim>>> TensorSL;

	//lapse
	Real alpha;

	//shift
	//beta_u(i) := beta^i
	//beta^t = 0
	TensorU beta_u;

	//log of the conformal factor (to the arbitrary root ...)
	//exp(phi) = psi
	//psi^12 = gamma
	Real phi;

	//conformal spatial metric
	//gammaBar_ll(i,j) := gammaBar_ij = exp(-4phi) gamma_ij
	// such that det(gammaBar_ij) = 1
	TensorSL gammaBar_ll;

	//minimal distortion conformal extrinsic curvature
	//ATilde_ll(i,j) := ATilde_ij = psi^-4 A_ij = exp(-4phi) A_ij = exp(-4phi) (K_ij - 1/3 gamma_ij K) = exp(-4phi) K_ij - 1/3 gammaBar_ij K
	TensorSL ATilde_ll;
	
	//extrinsic curvature trace
	//K := K^i_i
	Real K;

	//connection function
	//connBar_u(i) := connBar^i = gammaBar^jk connBar^i_jk = -partial_j gammaBar^ij
	TensorU connBar_u;

	GeomCell()
	:	alpha(Real()),
		phi(Real()),
		K(Real())
	{}

	//operators used with integration

	GeomCell operator*(const Real &scalar) const {
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

template<typename Real_, int dim_>
struct MatterCell {
	typedef Real_ Real;
	enum { dim = dim_ };

	typedef Tensor::Tensor<Real, Tensor::Upper<dim>> TensorU;
	typedef Tensor::Tensor<Real, Tensor::Symmetric<Tensor::Lower<dim>, Tensor::Lower<dim>>> TensorSL;

	MatterCell()
	:	rho(Real())
	{}

	//energy density
	//rho = n_a n_b T^ab
	Real rho;

	//momentum
	//S_u(i) := S^i = -gamma^ij * n^a T_aj
	TensorU S_u;

	//spatial stress energy
	//S_ll(i,j) := S_ij = gamma_ic gamma_jd T^cd
	//					= gamma_i^c gamma_j^d T_cd
	TensorSL S_ll;
};

template<typename Real_, int dim_>
struct AuxCell {
	typedef Real_ Real;
	enum { dim = dim_ };

	/*
	different ranked tensor types
	notation:
	Tensor[ULSA]*
	'U' means upper
	'L' means lower
	'S' means symmetric
	'A' means antisymmetric (haven't needed this one yet)
	*/
	typedef Tensor::Lower<dim> Lower;
	typedef Tensor::Upper<dim> Upper;
	
	typedef Tensor::Symmetric<Lower, Lower> SymmetricLower;
	typedef Tensor::Symmetric<Upper, Upper> SymmetricUpper;
	
	typedef Tensor::Tensor<Real, Upper> TensorU;
	typedef Tensor::Tensor<Real, Lower> TensorL;
	typedef Tensor::Tensor<Real, Lower, Lower> TensorLL;
	typedef Tensor::Tensor<Real, Upper, Lower> TensorUL;
	typedef Tensor::Tensor<Real, Lower, Upper> TensorLU;
	typedef Tensor::Tensor<Real, SymmetricUpper> TensorSU;
	typedef Tensor::Tensor<Real, SymmetricLower> TensorSL;
	typedef Tensor::Tensor<Real, Upper, SymmetricLower> TensorUSL;
	typedef Tensor::Tensor<Real, Lower, SymmetricLower> TensorLSL;


	//our tensors initialze to zero, so why not our reals too?
	AuxCell() 
	:	H(Real()),
		gamma(Real()),
		R(Real()),
		psi(Real()),
		DBar2_psi(Real()),
		RBar(Real()),
		tr_K_sq(Real())
	{}


		//constraint variables: should always be zero


	//hamiltonian constraint
	Real H;

	//momentum constraint
	TensorU M_u;


		//lapse and shift related


	//D_alpha_l(i) := D_i alpha
	TensorL D_alpha_l;

	//beta_l(i) := beta_i
	//beta_t = beta^k beta_k
	TensorL beta_l;

	//partial_beta_lu(i,j) := partial_i beta^j
	TensorLU partial_beta_lu;

		//metric related

	//spatial metric
	//gamma_ll(i,j) := gamma_ij
	//gamma_it = gamma_tj = 0
	TensorSL gamma_ll;

	//gamma_uu(i,j) := gamma^ij = inverse(gamma_ij) = covalent(gamma_ij) / det(gamma_ij)
	TensorSU gamma_uu;

	//gamma = det(gamma_ij)
	Real gamma;

	//conn_lll(i,j,k) := conn_ijk
	TensorLSL  conn_lll;

	//conn_ull(i,j,k) := conn^i_jk
	TensorUSL conn_ull;

	//RPhi_ll(i,j) := RPhi_ij is found by substituting phi for ln(psi) in eqn 3.10 of Baumgarte & Shapiro.
	// (but I thought this equation was for calculating R_ij?  And that ln(psi) = phi to begin with? So why do we now separate R_ij = RPhi_ij + RBar_ij?)
	TensorSL RPhi_ll;

	//R_ll(i,j) := R_ij = RBar_ij + RPhi_ij
	TensorSL R_ll;

	//Gaussian (scalar) curvature
	//R = R^i_i
	Real R;

		//conformal factor / metric

	//conformal factor 
	//currently derived from the iterated ln(sqrt(gamma))	
	Real psi;

	//DBar_psi_l(i) := DBar_i psi
	TensorL DBar_psi_l;

	//DBar_phi_l(i) := DBar_i ln(psi)
	TensorL DBar_phi_l;

	//DBar2_psi := gammaBar^ij DBar_i DBar_j psi
	Real DBar2_psi;
	
	//gammaBar_uu(i,j) := gammaBar^ij = psi^4 gamma^ij
	TensorSU gammaBar_uu;

	//partial_gammaBar_lll(i,j,k) := partial_i gammaBar_jk
	TensorLSL partial_gammaBar_lll;

	//partial2_gammaBar_ll(i,j) := gammaBar^lm partial_l partial_m gammaBar_ij
	// should the name be partialBar2? since I use gammaBar to raise it?
	TensorSL partial2_gammaBar_ll;

	//connBar_lll(i,j,k) := connBar_ijk = 1/2 (partial_k gammaBar_ij + partial_j gammaBar_ik - partial_i gammaBar_jk)
	TensorLSL connBar_lll;

	//connBar_ull(i,j,k) := connBar^i_jk = gammaBar^il connBar_ljk
	TensorUSL connBar_ull;

	//partial_connBar_lu(i,j) := partial_i connBar^j
	TensorLU partial_connBar_lu;

	//RBar_ll(i,j) := RBar_ij
	TensorSL RBar_ll;

	//RBar = gammaBar^ij RBar_ij
	Real RBar;
		
		//extrinsic curvature

	//tr_K_sq := (K^2)^i_i = K^ij K_ji
	Real tr_K_sq;

		//conformal extrinsic curvature

	//ABar_uu(i,j) := ABar^ij = psi^10 A^ij
	TensorSU ABar_uu;

	//ATilde_ul(i,j) := ATilde^i_j = gammaBar^ik ATilde_kj
	TensorUL ATilde_ul;

	//ATilde_uu(i,j) := ATilde^ij = ATilde^i_k gammaBar^kj
	TensorSU ATilde_uu;
};

