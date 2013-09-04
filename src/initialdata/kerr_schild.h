#pragma once

#include "initialdata.h"
#include "../constants.h"
#include "../exception.h"

#include <iostream>

/*
Kerr-Schild black hole
See Alcubierre p.56 and Baumgarte & Shapiro p.52
*/
template<typename real, int dim>
struct KerrSchild : public InitialData<real, dim> {
	virtual const char *name() { return "kerr-schild"; }

	typedef ::vector<real, dim> vector;
	typedef ::InitialData<real, dim> InitialData;
	typedef typename InitialData::ADMFormalism ADMFormalism;
	typedef typename ADMFormalism::tensor_l tensor_l;
	typedef typename ADMFormalism::tensor_u tensor_u;
	typedef typename ADMFormalism::tensor_sl tensor_sl;
	typedef typename ADMFormalism::tensor_su tensor_su;

	real M;		//black hole mass <=> half the Schwarzschild radius
	real J;		//total angular momentum
	real Q;		//total charge
	
	KerrSchild() : M(0), J(0), Q(0) {}

	virtual void init(ADMFormalism &sim, std::vector<std::string> &args) {
	
		if (!args.size()) throw Exception() << "expected mass";
		M = atof(args[0].c_str());
		args.erase(args.begin());
	
		J = 0;
		Q = 0;
		if (args.size()) {
			J = atof(args[0].c_str());
			args.erase(args.begin());
		
			if (args.size()) {
				Q = atof(args[0].c_str());
				args.erase(args.begin());
			}
		}

		std::cout << "mass " << M << " solar masses" << std::endl;
		std::cout << "angular momentum " << J << std::endl;
		std::cout << "charge " << Q << std::endl;

		M *= sunMassInM;
		
		real a = J / M;	//angular momentum density
		
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
				
			vector v = sim.coordForIndex(iter.index) - center;
			real r = vector::length(v);
			real x = v(0);
			real y = dim > 1 ? v(1) : 0;
			real z = dim > 2 ? v(2) : 0;
			real H = (r * M - Q * Q / 2.) / (r * r + a * a * z * z / (r * r));
			
			tensor_l l_l;
			l_l(0) = (r * x + a * y) / (r * r + a * a);
			if (dim > 1) l_l(1) = (r * y - a * x) / (r * r + a * a);
			if (dim > 2) l_l(2) = z / r;

			tensor_sl gamma_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					gamma_ll(i,j) = eta(i,j) + (2. - eta(i,j)) * H * l_l(i) * l_l(j);
				}
			}
	
			//gamma = det(gamma_ij)
			real gamma = determinant(gamma_ll);
			
			//phi = log(gamma) / 12
			real &phi = geomCell.phi;
			phi = log(gamma) / 12.;
			
			real expMinusFourPhi = exp(-4. * phi);
			//gammaBar_ll(i,j) := gammaBar_ij
			//					= psi^-4 gamma_ij 
			//					= exp(-4phi) gamma_ij
			tensor_sl &gammaBar_ll = geomCell.gammaBar_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					gammaBar_ll(i,j) = expMinusFourPhi * gamma_ll(i,j);
				}
			}
			
			tensor_su gamma_uu = inverse(gamma_ll, gamma);

			tensor_u l_u;
			for (int i = 0; i < dim; ++i) {
				l_u(i) = 0.;
				for (int j = 0; j < dim; ++j) {
					l_u(i) += gamma_uu(i,j) * l_l(j);
				}
			}

			tensor_u &beta_u = geomCell.beta_u;
			beta_u = l_u * H;

			real betaNorm = 0.;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					betaNorm += gamma_ll(i,j) * geomCell.beta_u(i) * geomCell.beta_u(j);
				}
			}

			real &alpha = geomCell.alpha;
			alpha = sqrt(1. - 2. * H - betaNorm);
			
			real &K = geomCell.K;
			K = 2 * M * alpha * alpha * alpha / (r * r) * (1. + 3. * M / r);

			tensor_sl &K_ll = geomCell.K_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					K_ll(i,j) = 2. * H * a / r * (eta(i,j) - (2. + H) * l_l(i) * l_l(j));
				}
			}
		}
	}
};

