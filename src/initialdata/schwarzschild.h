#pragma once

#include "initialdata.h"

//this is a Schwarzschild init, with some matter thrown in there
//TODO use Baumgarte & Shapiro p.66 Sobelev function for matter solution 
template<typename real, int dim>
struct Schwarzschild : public InitialData<real, dim> {
	typedef ::vector<real, dim> vector;
	typedef ::InitialData<real, dim> InitialData;
	typedef typename InitialData::ADMFormalism ADMFormalism;

	real M;		//total mass, in meters
	real R;		//radius, in meters
	real rho;	//density ... which I could calculate myself (in meters^-2)

	Schwarzschild(real M_, real R_, real rho_) 
	: M(M_), R(R_), rho(rho_) 
	{}

	virtual void init(ADMFormalism &sim) {
		
		vector center = (sim.max + sim.min) * .5;
		
		//provide initial conditions
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
			//ln_sqrt_gamma := ln(sqrt(gamma))
			geomCell.calc_ln_sqrt_gamma_from_gamma_ll();

			//K_ij = K = 0
			
			if (r <= R) {
				matterCell.rho = sim.dx.volume() * rho;
			}
		}
	}
};

