#pragma once

#include "initialdata.h"
#include "../constants.h"
#include "../exception.h"

#include <iostream>

/*
this is a Schwarzschild init from Baumgarte & Shapiro p.50, with some matter thrown in there
TODO use Baumgarte & Shapiro p.66 Sobelev function for matter solution 
*/
template<typename real, int dim>
struct Schwarzschild : public InitialData<real, dim> {
	virtual const char *name() { return "schwarzschild"; }

	typedef ::vector<real, dim> vector;
	typedef ::InitialData<real, dim> InitialData;
	typedef typename InitialData::ADMFormalism ADMFormalism;

	real M;		//total mass, in meters
	real R;		//radius, in meters
	real rho;	//density ... which I could calculate myself (in meters^-2)

	Schwarzschild() : M(0), R(0), rho(0) {}
	
	virtual void init(ADMFormalism &sim, std::vector<std::string> &args) {
		if (!args.size()) throw Exception() << "expected mass";
		M = atof(args[0].c_str());
		args.erase(args.begin());
		
		if (!args.size()) throw Exception() << "expected radius";
		R = atof(args[0].c_str());
		args.erase(args.begin());
	
		//technically ... I could calculate this ...
		if (!args.size()) throw Exception() << "expected density";
		rho = atof(args[0].c_str());
		args.erase(args.begin());
	
		std::cout << "mass " << M << " solar masses" << std::endl;
		std::cout << "radius " << R << " solar radii" << std::endl;
		std::cout << "density " << rho << std::endl;

		M *= sunMassInM;
		M *= sunRadiusInM;
		
		//provide initial conditions
		
		vector center = (sim.max + sim.min) * .5;
		std::cout << "providing initial conditions..." << std::endl;
		for (typename ADMFormalism::GeomGrid::iterator iter = sim.geomGridReadCurrent->begin(); iter != sim.geomGridReadCurrent->end(); ++iter) {
			typename ADMFormalism::GeomCell &geomCell = *iter;
			typename ADMFormalism::MatterCell &matterCell = sim.matterGrid(iter.index);

			vector v = sim.coordForIndex(iter.index) - center;
			
			real r = vector::length(v);
			real MOverTwoR = M / (2. * r);
			real oneMinusMOverTwoR = 1. - MOverTwoR;
			real onePlusMOverTwoR = 1. + MOverTwoR;
			real onePlusMOverTwoR_Squared = onePlusMOverTwoR * onePlusMOverTwoR;
			
			geomCell.alpha = oneMinusMOverTwoR / onePlusMOverTwoR; 
		
			//beta^i = 0

			for (int i = 0; i < dim; ++i) {
				geomCell.gamma_ll(i,i) = onePlusMOverTwoR_Squared * onePlusMOverTwoR_Squared;
			}
			
			//gamma = det(gamma_ij)
			real gamma = determinant(geomCell.gamma_ll);

			real &phi = geomCell.phi;
			phi = log(gamma) / 12.;

			//K_ij = K = 0
			
			if (r <= R) {
				matterCell.rho = sim.dx.volume() * rho;
			}
		}
	}
};

