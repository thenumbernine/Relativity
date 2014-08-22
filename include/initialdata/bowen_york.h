#pragma once

#include "initialdata.h"
#include "../constants.h"
#include "../exception.h"

#include <iostream>


/*
Bowen-York black hole
See Baumgarte & Shapiro p.70 and Alcubierre p.109

Supposed to be cooler than Kerr-Schild black holes in that they give off a burst of gravitational radiation before settling (p.72).
Funny how, when Baumgarte & Shapiro say "settling" on page 72, I am guessing they really mean "settling when you use our invented algorithm we don't teach you in this book until page 400".

Alcubierre gives the whole solution to spinning + boosting up front.  Baumgarte & Shapiro give you each separately.
But Baumgarte & Shapiro give you solutions to psi (in addition to the original Kerr-Schild solution), and Alcubierre only seems to give the Kerr-Schild 1 + M/2r
*/
template<typename real, int dim>
struct BowenYork : public InitialData<real, dim> {
	virtual const char *name() { return "bowen-york"; }

	typedef ::vector<real, dim> vector;
	typedef ::InitialData<real, dim> InitialData;
	typedef typename InitialData::ADMFormalism ADMFormalism;
	typedef typename ADMFormalism::tensor_l tensor_l;
	typedef typename ADMFormalism::tensor_u tensor_u;
	typedef typename ADMFormalism::tensor_sl tensor_sl;
	typedef typename ADMFormalism::tensor_su tensor_su;

	enum { dim3 = 3 };	//representing angular momentum in 3d even for 1d and 2d cases, so i can give a 2d black hole some angular momentum
	typedef tensor<real, lower<dim3>> tensor_l3;
	typedef tensor<real, upper<dim3>> tensor_u3;
	
	real M;						//black hole mass <=> half the Schwarzschild radius
	tensor_l3 J_l;	//angular momentum. technically only has to be a divergence-free field.
	tensor_u3 P_u;	//linear momentum. V_i = -2 P_i / r for V_i,j^j = 0, U_,j^j = 0, W_i = 7/8 V_i - 1/8 (U_,i + x^k  V_k,i)

	BowenYork() : M(0) {}

	virtual void init(ADMFormalism &sim, std::vector<std::string> &args) {
		if (!args.size()) throw Exception() << "expected mass";
		M = atof(args[0].c_str());
		args.erase(args.begin());
		
		for (int i = 0; i < dim3; ++i) {
			if (!args.size()) break;
			J_l(i) = atof(args[0].c_str());
			args.erase(args.begin());
		}

		for (int i = 0; i < dim3; ++i) {
			if (!args.size()) break;
			P_u(i) = atof(args[0].c_str());
			args.erase(args.begin());
		}

		std::cout << "mass " << M << " solar masses" << std::endl;
		std::cout << "angular momentum " << J_l << std::endl;
		std::cout << "linear momentum " << P_u << std::endl;

		M *= sunMassInM;
		
		//if we need J to calculate psi and we need psi to calculate gamma and we need gamma to calculate J (norm wrt metric) ... then we have a separate system to solve?
		//so until then I'm using the "approximation" J = |J^i|
		real J = tensor_l3::BodyType::length(J_l.body);
		//same deal with P?
		real P = tensor_l3::BodyType::length(P_u.body);
	
		tensor_sl eta;
		for (int i = 0; i < dim; ++i) {
			eta(i,i) = 1.;
		}

		const vector &min = sim.min;
		const vector &max = sim.max;

		//provide initial conditions
		
		vector center = (max + min) * .5;
		std::cout << "providing initial conditions..." << std::endl;
		for (typename ADMFormalism::GeomGrid::iterator iter = sim.geomGridReadCurrent->begin(); iter != sim.geomGridReadCurrent->end(); ++iter) {
			typename ADMFormalism::GeomCell &geomCell = *iter;
			//we don't need this cell, just a function in the AuxCell class for computing psi from phi 
			typename ADMFormalism::AuxCell &cell = sim.auxGrid(iter.index);
				
			vector x = sim.coordForIndex(iter.index) - center;
			
			//r = |x^i|
			real r = vector::length(x);
			real rSquared = r * r;
			real rCubed = r * rSquared;

			//l^i = x^i / r
			tensor_u3 l_u;
			for (int i = 0; i < dim; ++i) {
				l_u(i) = x(i) / r; 
			}

			real MOverTwoR = M / (2. * r);

			//psi0 = 1 + M / (2r)
			real psi0 = 1. + MOverTwoR;
	
			real psi0Squared = psi0;
			real psi0Cubed = psi0 * psi0Squared;
			real psi0ToTheFifth = psi0Squared * psi0Cubed;

			//psi20J = -(1 + M/(2r))^-5 M/(5r) (5(M/(2r))^3 + 4(M/(2r))^4 + (M/(2r))^5)
			real psi20J = -(M / (5. * r)) * (MOverTwoR * MOverTwoR * MOverTwoR * (5. + MOverTwoR * (4. + MOverTwoR))) / psi0ToTheFifth;

			//psi22J = -1/10 (1 + M/2r)^-5 (M/r)^3
			real psi22J = -(M * M * M) / (10. * rCubed * psi0ToTheFifth);
		
			//psi20P
			real psi20P = -M / (16. * r) * (5 + MOverTwoR * (10 + MOverTwoR * (10 + MOverTwoR * (5 + MOverTwoR)))) / psi0ToTheFifth;
			
			//psi22P
			real psi22P = MOverTwoR * MOverTwoR * (15 + MOverTwoR * (192 + MOverTwoR * (539 + MOverTwoR * (658 + MOverTwoR * (378 + MOverTwoR * 84))))) / (20. * psi0ToTheFifth)
						+ 21. / 5. * MOverTwoR * MOverTwoR * MOverTwoR * log(MOverTwoR / (1. + MOverTwoR));

			real cosTheta = l_u(0);

			//P0(cos(theta)) = 1
			real P0CosTheta = 1;

			//P2(cos(theta)) = (3 cos(theta)^2 - 1) / 2
			real P2CosTheta = (3. * cosTheta * cosTheta - 1.) / 2.;

			//psi2J = psi20J P0(cos(theta)) + psi22J P2(cos(theta))
			real psi2J = P0CosTheta * psi20J + P2CosTheta * psi22J;

			//psi2P = psi20P P0(cos(theta)) + psi22P P2(cos(theta));
			real psi2P = P0CosTheta * psi20P + P2CosTheta * psi22P;
	
			//spinning:
			//psi = psi0 + psi2J J^2 / M^4 + O(J^4)
			//boosted:
			//psi = psi0 + psi2P P^2 / M^2 psi2P + O(P^4)
			//combined?
			//psi = psi0 + psi2J J^2 / M^4 + psi2P P^2 / M^2 + O(J^4) + O(P^4)
			real &psi = cell.psi;
			psi = psi0 + psi2J * J * J / (M * M * M * M) + psi2P * P * P / (M * M);
			
			real psiSquared = psi * psi;
			real psiToTheFourth = psiSquared * psiSquared;
			
			//phi = log(psi)
			real &phi = geomCell.phi;
			phi = log(psi);

			//gammaBar_ij = eta_ij
			tensor_sl &gammaBar_ll = geomCell.gammaBar_ll;
			gammaBar_ll = eta;;
			
			//l_i = gamma_ij l^j
			// since I'm forcing angular momentum to be 3D (so I can get angular momentum in my 2D cases)
			// and both l and J go into ABarL, so this has to be 3D as well
			// should its lower form (which is lowered by eta times psi^4) omit the extra dimensions?
			// or should it (as i'm doing) assume it is being lowered by a 3D eta with psi^4 conformal factor?
			tensor_l3 l_l;
			for (int i = 0; i < dim3; ++i) {
				l_l(i) = l_u(i) * psiToTheFourth;
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
					ABarL_uu(i,j) = 3. / (2. * r * r) * (P_u(i) * l_u(j) + P_u(j) * l_u(i));
					for (int k = 0; k < dim3; ++k) {
						ABarL_uu(i,j) += 3. / (2. * r * r) * (-(eta(i,j) - l_u(i) * l_u(j))) * l_l(k) * P_u(k);
						for (int l = 0; l < dim3; ++l) {
							if (k == (j+1)%dim3 && l == (k+1)%dim3) {
								ABarL_uu(i,j) += 3. / rCubed * l_u(i) * J_l(k) * l_l(l);
							} else if (j == (k+1)%dim3 && k == (l+1)%dim3) {
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

			//K_ll(i,j) := K_ij 
			//			= A_ij - 1/3 gamma_ij K 
			//			= psi^-2 ABar_ij - 1/3 gamma_ij K
			//			= psi^-2 ABar_ij - 1/3 psi^4 gammaBar_ij K
			tensor_sl K_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					K_ll(i,j) = ABar_ll(i,j) * oneOverPsiSquared - 1./3. * psiToTheFourth * gammaBar_ll(i,j) * K;
				}
			}

			//ATilde_ll(i,j) := ATilde_ij = exp(-4phi) K_ij - 1/3 gammaBar_ij K
			tensor_sl &ATilde_ll = geomCell.ATilde_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					ATilde_ll(i,j) = psiToTheFourth * K_ll(i,j) - 1./3. * gammaBar_ll(i,j) * K;
				}
			}

		}
	}
};

