#pragma once

#include "initialdata.h"
#include "../constants.h"
#include "../parallel.h"
#include "Common/Exception.h"
#include <iostream>

/*
this is a Schwarzschild init from Baumgarte & Shapiro p.50, with some matter thrown in there
TODO use Baumgarte & Shapiro p.66 Sobelev function for matter solution 
*/
template<typename Real, int dim>
struct Schwarzschild : public InitialData<Real, dim> {
	virtual const char *name() { return "schwarzschild"; }

	typedef Tensor::Vector<Real, dim> Vector;
	typedef Tensor::Vector<int, dim> DerefType;
	typedef ::InitialData<Real, dim> InitialData;
	typedef typename InitialData::ADMFormalism ADMFormalism;
	typedef typename ADMFormalism::TensorSL TensorSL;

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
		parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, dim> index) {
			typename ADMFormalism::GeomCell &geomCell = (*sim.geomGridReadCurrent)(index);
			typename ADMFormalism::MatterCell &matterCell = sim.matterGrid(index);

			Vector v = sim.coordForIndex(index) - center;
			
			Real r = Vector::length(v);
			Real MOverTwoR = M / (2. * r);
			Real oneMinusMOverTwoR = 1. - MOverTwoR;
			Real onePlusMOverTwoR = 1. + MOverTwoR;
			Real onePlusMOverTwoR_Squared = onePlusMOverTwoR * onePlusMOverTwoR;
			
			geomCell.alpha = oneMinusMOverTwoR / onePlusMOverTwoR; 
		
			//beta^i = 0

			TensorSL gamma_ll;
			for (int i = 0; i < dim; ++i) {
				gamma_ll(i,i) = onePlusMOverTwoR_Squared * onePlusMOverTwoR_Squared;
			}
			
			//gamma = det(gamma_ij)
			Real gamma = determinant(gamma_ll);
			
			//phi = log(gamma) / 12
			Real &phi = geomCell.phi;
			phi = log(gamma) / 12.;

			Real expMinusFourPhi = exp(-4. * phi);
			//gammaBar_ll(i,j) := gamma^-1/3 gamma_ij
			TensorSL &gammaBar_ll = geomCell.gammaBar_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					gammaBar_ll(i,j) = expMinusFourPhi * gamma_ll(i,j);
				}
			}

			//K_ij = K = 0
			
			if (r <= R) {
				matterCell.rho = sim.dx.volume() * rho;
			}
		});
	}
};

