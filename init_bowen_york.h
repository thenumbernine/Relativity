#pragma once

#include "admformalism.h"
#include "initialdata.h"

/*
Bowen-York spinning black hole
See Baumgarte & Shapiro p.70 and Alcubierre p.109
*/
template<typename real, int dim>
struct BowenYork : public InitialData<real, dim> {
	typedef ::vector<real, dim> vector;
	typedef ::InitialData<real, dim> InitialData;
	typedef typename InitialData::ADMFormalism ADMFormalism;
	typedef typename ADMFormalism::tensor_l tensor_l;
	typedef typename ADMFormalism::tensor_u tensor_u;
	typedef typename ADMFormalism::tensor_sl tensor_sl;
	typedef typename ADMFormalism::tensor_su tensor_su;

	real M;			//black hole mass <=> half the Schwarzschild radius
	tensor_l J_l;	//angular momentum. technically only has to be a divergence-free field.

	BowenYork(
		real M_,		//mass of black hole
		tensor_l J_) 	//angular momentum of black hole (Alcubierre calls this one 'S', but Baumgarte & Shapiro already use 'S' for stress-energy momentum, which Alcubierre calls 'j')
	: 	M(M_),
		J_l(J_)
	{}

	virtual void init(ADMFormalism &sim) {
		
		//J = |J^i|
		real J = tensor_l::body::length(J_l.body);
	
		tensor_sl eta;
		for (int i = 0; i < dim; ++i) {
			eta(i,i) = 1.;
		}

		const vector &min = sim.min;
		const vector &max = sim.max;

		//provide initial conditions
		
		vector center = (max + min) * .5;
		for (typename ADMFormalism::GeomGrid::iterator iter = sim.geomGridReadCurrent->begin(); iter != sim.geomGridReadCurrent->end(); ++iter) {
			typename ADMFormalism::GeomCell &geomCell = *iter;
			//we don't need this cell, just a function in the AuxCell class for computing psi from ln_sqrt_gamma
			typename ADMFormalism::AuxCell &cell = sim.auxGrid(iter.index);
				
			vector x = sim.coordForIndex(iter.index) - center;
			
			//r = |x^i|
			real r = vector::length(x);
			real rSquared = r * r;
			real rCubed = r * rSquared;

			//l^i = x^i / r
			tensor_u l_u;
			for (int i = 0; i < dim; ++i) {
				l_u(i) = x(i) / r; 
			}

			real MOverTwoR = M / (2. * r);

			//psi0 = 1 + M / (2r)
			real psi0 = 1. + MOverTwoR;
	
			real psi0Squared = psi0;
			real psi0Cubed = psi0 * psi0Squared;
			real psi0ToTheFifth = psi0Squared * psi0Cubed;

			//psi20 = -(1 + M/(2r))^-5 M/(5r) (5(M/(2r))^3 + 4(M/(2r))^4 + (M/(2r))^5)
			real psi20 = -(M / (5. * r)) * (MOverTwoR * MOverTwoR * MOverTwoR * (5. + MOverTwoR * (4. + MOverTwoR))) / psi0ToTheFifth;

			//psi22 = -1/10 (1 + M/2r)^-5 (M/r)^3
			real psi22 = -(M * M * M) / (10. * rCubed * psi0ToTheFifth);

			//psi2 = psi20 P0(cos(theta)) + psi22 P2(cos(theta))
			//	P0 = 1, P2(cos(theta) = (3 * cos(theta)^2 - 1) / 2
			real cosTheta = l_u(0);
			real psi2 = psi20 + psi22 * (3. * cosTheta * cosTheta - 1.) / 2.;

			//psi = psi0 + psi2 J^2 / M^4 + O(J^4)
			real &psi = cell.psi;
			psi = psi0 + psi2 * (J * J) / (M * M * M * M);
			real psiSquared = psi * psi;
			
			//extract our state conformal factor variable
			cell.ln_psi = log(cell.psi);
			geomCell.ln_sqrt_gamma = 6. * cell.ln_psi;

			//gammaBar_ij = eta_ij
			tensor_sl &gammaBar_ll = cell.gammaBar_ll;
			gammaBar_ll = eta;;
			
			//reconstruct our state metric variable
			//gamma_ij = psi^4 gammaBar_ij
			tensor_sl &gamma_ll = geomCell.gamma_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					gamma_ll(i,j) = psiSquared * psiSquared * gammaBar_ll(i,j);
				}
			}

			//l_i = gamma_ij l^j
			tensor_l l_l;
			for (int i = 0; i < dim; ++i) {
				l_l(i) = 0;
				for (int j = 0; j < dim; ++j) {
					l_l(i) += gamma_ll(i,j) * l_u(j);
				}
			}

			//X^i = l^i / r^2
			tensor_u X_u;
			for (int i = 0; i < dim; ++i) {
				X_u(i) = l_u(i) / rSquared;
			}

			//alpha = 1
			geomCell.alpha = 1;

			//beta^i = 0
		
			//ABarLL^ij = (LBar W)^ij = 6/r^3 l(^i eBar^j)^kl J_k l_l
			tensor_su ABarL_uu;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					ABarL_uu(i,j) = 0;
					for (int k = 0; k < dim; ++k) {
						for (int l = 0; l < dim; ++l) {
							if (k == (j+1)%dim && l == (k+1)%dim) {
								ABarL_uu(i,j) += 3. / rCubed * l_u(i) * J_l(k) * l_l(l);
							} else if (j == (k+1)%dim && k == (l+1)%dim) {
								ABarL_uu(i,j) -= 3. / rCubed * l_u(i) * J_l(k) * l_l(l);
							}
						}
					}
				}
			}

			// free to specify / leave at zero
			tensor_sl ABarTT_uu;

			//ABar^ij = ABarTT^ij + ABarL^ij
			tensor_sl ABar_uu;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					ABar_uu(i,j) = ABarTT_uu(i,j) + ABarL_uu(i,j);
				}
			}

			//ABar_ij = gammaBar_ik ABar^kl gammaBar_lj
			tensor_sl ABar_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					ABar_ll(i,j) = 0;
					for (int k = 0; k < dim; ++k) {
						for (int l = 0; l < dim; ++l) {
							ABar_ll(i,j) += gammaBar_ll(i,k) * ABar_uu(k,l) * gammaBar_ll(l,j);
						}
					}
				}
			}

			//K = 0
			real &K = geomCell.K;
			K = 0;

			real oneOverPsiSquared = 1. / psiSquared;

			//K_ll(i,j) := K_ij = A_ij - 1/3 gamma_ij K = psi^-2 ABar_ij - 1/3 gamma_ij K
			tensor_sl &K_ll = geomCell.K_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					K_ll(i,j) = ABar_ll(i,j) * oneOverPsiSquared - 1./3. * gamma_ll(i,j) * K;
				}
			}
		}
	}
};

