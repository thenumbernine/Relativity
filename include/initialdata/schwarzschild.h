#pragma once

#include "initialdata/initialdata.h"
#include "constants.h"
#include "parallel.h"
#include "Common/Exception.h"
#include <iostream>

/*
this is a Schwarzschild init from Baumgarte & Shapiro p.50, with some matter thrown in there
TODO use Baumgarte & Shapiro p.66 Sobelev function for matter solution 
*/
template<typename Real, int dim>
struct Schwarzschild : public InitialData<Real, dim> {
	virtual char const * name() { return "schwarzschild"; }

	using Vector = Tensor::vec<Real, dim>;
	using DerefType = Tensor::intN<dim>;
	using InitialData = ::InitialData<Real, dim>;
	using ADMFormalism = typename InitialData::ADMFormalism;
	using TensorSL = typename ADMFormalism::TensorSL;

	Real M;		//total mass, in meters
	Real R;		//radius, in meters
	Real rho;	//density ... which I could calculate myself (in meters^-2)

	Schwarzschild() : M(0), R(0), rho(0) {}
	
	virtual void init(ADMFormalism &sim, std::vector<std::string> &args) {
		if (!args.size()) throw Common::Exception() << "expected mass";
		M = atof(args[0].c_str());
		args.erase(args.begin());
		
		if (!args.size()) throw Common::Exception() << "expected radius";
		R = atof(args[0].c_str());
		args.erase(args.begin());
	
		//technically ... I could calculate this ...
		if (!args.size()) throw Common::Exception() << "expected density";
		rho = atof(args[0].c_str());
		args.erase(args.begin());
	
		std::cout << "mass " << M << " solar masses" << std::endl;
		std::cout << "radius " << R << " solar radii" << std::endl;
		std::cout << "density " << rho << std::endl;

		M *= sunMassInM;
		M *= sunRadiusInM;
		
		//provide initial conditions
		
		Vector center = (sim.max + sim.min) * .5;
		std::cout << "providing initial conditions..." << std::endl;
		Tensor::RangeObj<dim> range = sim.geomGridReadCurrent->range();
		parallel.foreach(range.begin(), range.end(), [&](Tensor::intN<dim> index) {
			typename ADMFormalism::GeomCell &geomCell = (*sim.geomGridReadCurrent)(index);
			typename ADMFormalism::MatterCell &matterCell = sim.matterGrid(index);

			Vector v = sim.coordForIndex(index) - center;
			
			Real r = v.length();
			Real MOverTwoR = M / (2. * r);
			Real oneMinusMOverTwoR = 1. - MOverTwoR;
			Real onePlusMOverTwoR = 1. + MOverTwoR;
			Real onePlusMOverTwoR_Squared = onePlusMOverTwoR * onePlusMOverTwoR;
			
			geomCell.alpha = oneMinusMOverTwoR / onePlusMOverTwoR; 
		
			//beta^i = 0

			TensorSL gamma_ll = Tensor::ident<Real,dim>(onePlusMOverTwoR_Squared * onePlusMOverTwoR_Squared);
			
			//gamma = det(gamma_ij)
			Real gamma = determinant(gamma_ll);
			
			//phi = log(gamma) / 12
			geomCell.phi = log(gamma) / 12.;

			Real expMinusFourPhi = exp(-4. * geomCell.phi);
			//gammaBar_ll(i,j) := gamma^-1/3 gamma_ij
			geomCell.gammaBar_ll = expMinusFourPhi * gamma_ll;

			//K_ij = K = 0
			
			if (r <= R) {
				matterCell.rho = sim.dx.volume() * rho;
			}
		});
	}
};
