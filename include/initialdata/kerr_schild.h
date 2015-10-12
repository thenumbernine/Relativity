#pragma once

#include "initialdata.h"
#include "../constants.h"
#include "../parallel.h"
#include "Common/Exception.h"
#include <iostream>

/*
Kerr-Schild black hole
See Alcubierre p.56 and Baumgarte & Shapiro p.52
*/
template<typename Real, int dim>
struct KerrSchild : public InitialData<Real, dim> {
	virtual const char *name() { return "kerr-schild"; }

	typedef Tensor::Vector<Real, dim> Vector;
	typedef Tensor::Vector<int, dim> DerefType;
	typedef ::InitialData<Real, dim> InitialData;
	typedef typename InitialData::ADMFormalism ADMFormalism;
	typedef typename ADMFormalism::TensorL TensorL;
	typedef typename ADMFormalism::TensorU TensorU;
	typedef typename ADMFormalism::TensorSL TensorSL;
	typedef typename ADMFormalism::TensorSU TensorSU;

	Real M;		//black hole mass <=> half the Schwarzschild radius
	Real J;		//total angular momentum
	Real Q;		//total charge
	
	KerrSchild() : M(0), J(0), Q(0) {}

	virtual void init(ADMFormalism &sim, std::vector<std::string> &args) {
	
		if (!args.size()) throw Common::Exception() << "expected mass";
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
		
		Real a = J / M;	//angular momentum density
		
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
				
			Vector v = sim.coordForIndex(index) - center;
			Real r = Vector::length(v);
			Real x = v(0);
			Real y = dim > 1 ? v(1) : 0;
			Real z = dim > 2 ? v(2) : 0;
			Real H = (r * M - Q * Q / 2.) / (r * r + a * a * z * z / (r * r));
			
			TensorL l_l;
			l_l(0) = (r * x + a * y) / (r * r + a * a);
			if (dim > 1) l_l(1) = (r * y - a * x) / (r * r + a * a);
			if (dim > 2) l_l(2) = z / r;

			TensorSL gamma_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					gamma_ll(i,j) = eta(i,j) + (2. - eta(i,j)) * H * l_l(i) * l_l(j);
				}
			}
	
			//gamma = det(gamma_ij)
			Real gamma = determinant(gamma_ll);
			
			//phi = log(gamma) / 12
			Real &phi = geomCell.phi;
			phi = log(gamma) / 12.;
			
			Real expMinusFourPhi = exp(-4. * phi);
			//gammaBar_ll(i,j) := gammaBar_ij
			//					= psi^-4 gamma_ij 
			//					= exp(-4phi) gamma_ij
			TensorSL &gammaBar_ll = geomCell.gammaBar_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					gammaBar_ll(i,j) = expMinusFourPhi * gamma_ll(i,j);
				}
			}
			
			TensorSU gamma_uu = inverse(gamma_ll, gamma);

			TensorU l_u;
			for (int i = 0; i < dim; ++i) {
				l_u(i) = 0.;
				for (int j = 0; j < dim; ++j) {
					l_u(i) += gamma_uu(i,j) * l_l(j);
				}
			}

			TensorU &beta_u = geomCell.beta_u;
			beta_u = l_u * H;

			Real betaNorm = 0.;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					betaNorm += gamma_ll(i,j) * geomCell.beta_u(i) * geomCell.beta_u(j);
				}
			}

			Real &alpha = geomCell.alpha;
			alpha = sqrt(1. - 2. * H - betaNorm);
			
			Real &K = geomCell.K;
			K = 2 * M * alpha * alpha * alpha / (r * r) * (1. + 3. * M / r);

			TensorSL K_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					K_ll(i,j) = 2. * H * a / r * (eta(i,j) - (2. + H) * l_l(i) * l_l(j));
				}
			}

			//ATilde_ll(i,j) := ATilde_ij = exp(-4phi) K_ij - 1/3 gammaBar_ij K
			TensorSL &ATilde_ll = geomCell.ATilde_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					ATilde_ll(i,j) = expMinusFourPhi * K_ll(i,j) - 1./3. * gammaBar_ll(i,j) * K;
				}
			}
		});
	}
};

