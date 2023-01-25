#pragma once

#include "initialdata/initialdata.h"
#include "constants.h"
#include "parallel.h"
#include "Common/Exception.h"
#include <iostream>

/*
Kerr-Schild black hole
See Alcubierre p.56 and Baumgarte & Shapiro p.52
*/
template<typename Real, int dim>
struct KerrSchild : public InitialData<Real, dim> {
	virtual char const * name() { return "kerr-schild"; }

	using Vector = Tensor::vec<Real, dim>;
	using DerefType = Tensor::intN<dim>;
	using InitialData = ::InitialData<Real, dim>;
	using ADMFormalism = typename InitialData::ADMFormalism;
	using TensorL = typename ADMFormalism::TensorL;
	using TensorU = typename ADMFormalism::TensorU;
	using TensorSL = typename ADMFormalism::TensorSL;
	using TensorSU = typename ADMFormalism::TensorSU;

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
		
		auto eta = Tensor::ident<Real,dim>(1);

		Vector const & min = sim.min;
		Vector const & max = sim.max;

		//provide initial conditions
		
		Vector center = (max + min) * .5;
		std::cout << "providing initial conditions..." << std::endl;
		Tensor::RangeObj<dim> range = sim.geomGridReadCurrent->range();
		parallel.foreach(range.begin(), range.end(), [&](Tensor::intN<dim> index) {
			typename ADMFormalism::GeomCell &geomCell = (*sim.geomGridReadCurrent)(index);
				
			Vector v = sim.coordForIndex(index) - center;
			Real r = v.length();
			Real x = v(0);
			Real y = dim > 1 ? v(1) : 0;
			Real z = dim > 2 ? v(2) : 0;
			Real H = (r * M - Q * Q / 2.) / (r * r + a * a * z * z / (r * r));
			
			TensorL l_l;
			l_l(0) = (r * x + a * y) / (r * r + a * a);
			if (dim > 1) l_l(1) = (r * y - a * x) / (r * r + a * a);
			if (dim > 2) l_l(2) = z / r;

			TensorSL gamma_ll = eta + (2. - eta) * H * Tensor::outer(l_l, l_l);
	
			//gamma = det(gamma_ij)
			Real gamma = determinant(gamma_ll);
			
			//phi = log(gamma) / 12
			geomCell.phi = log(gamma) / 12.;
			
			Real expMinusFourPhi = exp(-4. * geomCell.phi);
			//gammaBar_ll(i,j) := gammaBar_ij
			//					= psi^-4 gamma_ij 
			//					= exp(-4phi) gamma_ij
			geomCell.gammaBar_ll = expMinusFourPhi * gamma_ll;
			
			TensorSU gamma_uu = inverse(gamma_ll, gamma);

			TensorU l_u = gamma_uu * l_l;

			geomCell.beta_u = l_u * H;

			Real betaNorm = geomCell.beta_u * gamma_ll * geomCell.beta_u;

			geomCell.alpha = sqrt(1. - 2. * H - betaNorm);
			
			Real &K = geomCell.K;
			K = 2 * M * geomCell.alpha * geomCell.alpha * geomCell.alpha / (r * r) * (1. + 3. * M / r);

			TensorSL K_ll = 2. * H * a / r * (eta - (2. + H) * Tensor::outer(l_l, l_l));

			//ATilde_ll(i,j) := ATilde_ij = exp(-4phi) K_ij - 1/3 gammaBar_ij K
			geomCell.ATilde_ll = expMinusFourPhi * K_ll - 1./3. * geomCell.gammaBar_ll * K;
		});
	}
};
