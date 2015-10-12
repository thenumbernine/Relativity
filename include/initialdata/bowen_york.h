#pragma once

#include "initialdata.h"
#include "../constants.h"
#include "../parallel.h"
#include "Common/Exception.h"
#include <iostream>


/*
Bowen-York black hole
See Baumgarte & Shapiro p.70 and Alcubierre p.109

Supposed to be cooler than Kerr-Schild black holes in that they give off a burst of gravitational radiation before settling (p.72).
Funny how, when Baumgarte & Shapiro say "settling" on page 72, I am guessing they really mean "settling when you use our invented algorithm we don't teach you in this book until page 400".

Alcubierre gives the whole solution to spinning + boosting up front.  Baumgarte & Shapiro give you each separately.
But Baumgarte & Shapiro give you solutions to psi (in addition to the original Kerr-Schild solution), and Alcubierre only seems to give the Kerr-Schild 1 + M/2r
*/
template<typename Real, int dim>
struct BowenYork : public InitialData<Real, dim> {
	virtual const char *name() { return "bowen-york"; }

	typedef Tensor::Vector<Real, dim> Vector;
	typedef Tensor::Vector<int, dim> DerefType;
	typedef ::InitialData<Real, dim> InitialData;
	typedef typename InitialData::ADMFormalism ADMFormalism;
	typedef typename ADMFormalism::TensorL TensorL;
	typedef typename ADMFormalism::TensorU TensorU;
	typedef typename ADMFormalism::TensorSL TensorSL;
	typedef typename ADMFormalism::TensorSU TensorSU;

	enum { dim3 = 3 };	//representing angular momentum in 3d even for 1d and 2d cases, so i can give a 2d black hole some angular momentum
	typedef Tensor::Tensor<Real, Tensor::Lower<dim3>> TensorL3;
	typedef Tensor::Tensor<Real, Tensor::Upper<dim3>> TensorU3;
	
	Real M;						//black hole mass <=> half the Schwarzschild radius
	TensorL3 J_l;	//angular momentum. technically only has to be a divergence-free field.
	TensorU3 P_u;	//linear momentum. V_i = -2 P_i / r for V_i,j^j = 0, U_,j^j = 0, W_i = 7/8 V_i - 1/8 (U_,i + x^k  V_k,i)

	BowenYork() : M(0) {}

	virtual void init(ADMFormalism &sim, std::vector<std::string> &args) {
		if (!args.size()) throw Common::Exception() << "expected mass";
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
		Real J = TensorL3::BodyType::length(J_l.body);
		//same deal with P?
		Real P = TensorL3::BodyType::length(P_u.body);
	
		TensorSL eta;
		for (int i = 0; i < dim; ++i) {
			eta(i,i) = 1.;
		}

		const Vector &min = sim.min;
		const Vector &max = sim.max;

		//provide initial conditions
		
		Vector center = (max + min) * .5;
		std::cout << "providing initial conditions..." << std::endl;
		Tensor::RangeObj<dim> range = sim.geomGridReadCurrent->range();
		parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, dim> index) {
			typename ADMFormalism::GeomCell &geomCell = (*sim.geomGridReadCurrent)(index);
			//we don't need this cell, just a function in the AuxCell class for computing psi from phi 
			typename ADMFormalism::AuxCell &cell = sim.auxGrid(index);
				
			Vector x = sim.coordForIndex(index) - center;
			
			//r = |x^i|
			Real r = Vector::length(x);
			Real rSquared = r * r;
			Real rCubed = r * rSquared;

			//l^i = x^i / r
			TensorU3 l_u;
			for (int i = 0; i < dim; ++i) {
				l_u(i) = x(i) / r; 
			}

			Real MOverTwoR = M / (2. * r);

			//psi0 = 1 + M / (2r)
			Real psi0 = 1. + MOverTwoR;
	
			Real psi0Squared = psi0;
			Real psi0Cubed = psi0 * psi0Squared;
			Real psi0ToTheFifth = psi0Squared * psi0Cubed;

			//psi20J = -(1 + M/(2r))^-5 M/(5r) (5(M/(2r))^3 + 4(M/(2r))^4 + (M/(2r))^5)
			Real psi20J = -(M / (5. * r)) * (MOverTwoR * MOverTwoR * MOverTwoR * (5. + MOverTwoR * (4. + MOverTwoR))) / psi0ToTheFifth;

			//psi22J = -1/10 (1 + M/2r)^-5 (M/r)^3
			Real psi22J = -(M * M * M) / (10. * rCubed * psi0ToTheFifth);
		
			//psi20P
			Real psi20P = -M / (16. * r) * (5 + MOverTwoR * (10 + MOverTwoR * (10 + MOverTwoR * (5 + MOverTwoR)))) / psi0ToTheFifth;
			
			//psi22P
			Real psi22P = MOverTwoR * MOverTwoR * (15 + MOverTwoR * (192 + MOverTwoR * (539 + MOverTwoR * (658 + MOverTwoR * (378 + MOverTwoR * 84))))) / (20. * psi0ToTheFifth)
						+ 21. / 5. * MOverTwoR * MOverTwoR * MOverTwoR * log(MOverTwoR / (1. + MOverTwoR));

			Real cosTheta = l_u(0);

			//P0(cos(theta)) = 1
			Real P0CosTheta = 1;

			//P2(cos(theta)) = (3 cos(theta)^2 - 1) / 2
			Real P2CosTheta = (3. * cosTheta * cosTheta - 1.) / 2.;

			//psi2J = psi20J P0(cos(theta)) + psi22J P2(cos(theta))
			Real psi2J = P0CosTheta * psi20J + P2CosTheta * psi22J;

			//psi2P = psi20P P0(cos(theta)) + psi22P P2(cos(theta));
			Real psi2P = P0CosTheta * psi20P + P2CosTheta * psi22P;
	
			//spinning:
			//psi = psi0 + psi2J J^2 / M^4 + O(J^4)
			//boosted:
			//psi = psi0 + psi2P P^2 / M^2 psi2P + O(P^4)
			//combined?
			//psi = psi0 + psi2J J^2 / M^4 + psi2P P^2 / M^2 + O(J^4) + O(P^4)
			Real &psi = cell.psi;
			psi = psi0 + psi2J * J * J / (M * M * M * M) + psi2P * P * P / (M * M);
			
			Real psiSquared = psi * psi;
			Real psiToTheFourth = psiSquared * psiSquared;
			
			//phi = log(psi)
			Real &phi = geomCell.phi;
			phi = log(psi);

			//gammaBar_ij = eta_ij
			TensorSL &gammaBar_ll = geomCell.gammaBar_ll;
			gammaBar_ll = eta;;
			
			//l_i = gamma_ij l^j
			// since I'm forcing angular momentum to be 3D (so I can get angular momentum in my 2D cases)
			// and both l and J go into ABarL, so this has to be 3D as well
			// should its lower form (which is lowered by eta times psi^4) omit the extra dimensions?
			// or should it (as i'm doing) assume it is being lowered by a 3D eta with psi^4 conformal factor?
			TensorL3 l_l;
			for (int i = 0; i < dim3; ++i) {
				l_l(i) = l_u(i) * psiToTheFourth;
			}

			//X^i = l^i / r^2
			TensorU X_u;
			for (int i = 0; i < dim; ++i) {
				X_u(i) = l_u(i) / rSquared;
			}

			//alpha = 1
			geomCell.alpha = 1;

			//beta^i = 0
		
			//ABarLL^ij = (LBar W)^ij = 6/r^3 l(^i eBar^j)^kl J_k l_l
			TensorSU ABarL_uu;
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
			TensorSL ABarTT_uu;

			//ABar^ij = ABarTT^ij + ABarL^ij
			TensorSL ABar_uu;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					ABar_uu(i,j) = ABarTT_uu(i,j) + ABarL_uu(i,j);
				}
			}

			//ABar_ij = gammaBar_ik ABar^kl gammaBar_lj
			TensorSL ABar_ll;
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
			Real &K = geomCell.K;
			K = 0;

			Real oneOverPsiSquared = 1. / psiSquared;

			//K_ll(i,j) := K_ij 
			//			= A_ij - 1/3 gamma_ij K 
			//			= psi^-2 ABar_ij - 1/3 gamma_ij K
			//			= psi^-2 ABar_ij - 1/3 psi^4 gammaBar_ij K
			TensorSL K_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					K_ll(i,j) = ABar_ll(i,j) * oneOverPsiSquared - 1./3. * psiToTheFourth * gammaBar_ll(i,j) * K;
				}
			}

			//ATilde_ll(i,j) := ATilde_ij = exp(-4phi) K_ij - 1/3 gammaBar_ij K
			TensorSL &ATilde_ll = geomCell.ATilde_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					ATilde_ll(i,j) = psiToTheFourth * K_ll(i,j) - 1./3. * gammaBar_ll(i,j) * K;
				}
			}
		});
	}
};

