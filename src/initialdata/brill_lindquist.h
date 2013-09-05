#pragma once

#include "initialdata.h"
#include "../constants.h"
#include "../exception.h"

#include <iostream>

/*
This one's coming from Baumgarte & Shapiro p.61:
beta^i = 0
gammaBar_ij = eta_ij
rho = S^i = 0
K_ij = K = 0
psi = 1 + sum_a M_a / (2 r_a)
	r_i = |x^i - C^i_a| 
	M_i = mass of the i'th black hole

...and more detail can be found in chapter 13.2
In fact, I'm stealing 1/alpha = sum_a M_a / (r c_a) from 12.51 without fully reading the rest of the chapter
*/
template<typename real, int dim>
struct BrillLindquist : public InitialData<real, dim> {
	virtual const char *name() { return "brill-lindquist"; }

	typedef ::vector<real, dim> vector;
	typedef ::InitialData<real, dim> InitialData;
	typedef typename InitialData::ADMFormalism ADMFormalism;
	typedef typename ADMFormalism::tensor_l tensor_l;
	typedef typename ADMFormalism::tensor_u tensor_u;
	typedef typename ADMFormalism::tensor_sl tensor_sl;
	typedef typename ADMFormalism::tensor_su tensor_su;

	//I should at least make this a structure or something
	std::vector<tensor<real, lower<dim+1>>> blackHoleInfo;

	virtual void init(ADMFormalism &sim, std::vector<std::string> &args) {
		if (!args.size()) throw Exception() << "expected simulation arguments";
		int numBlackHoles = atoi(args[0].c_str());
		args.erase(args.begin());

		std::cout << "number of black holes " << numBlackHoles << std::endl;

		for (int i = 0; i < numBlackHoles; ++i) {
			tensor<real, lower<dim+1>> blackHole;
			for (int j = 0; j < dim+1; ++j) {
				if (!args.size()) throw Exception() << "expected simulation arguments";
				blackHole(j) = atof(args[0].c_str());
				args.erase(args.begin());
			}
			
			std::cout << "black hole position and mass " << blackHole << std::endl;
			
			for (int j = 0; j < dim; ++j) {
				blackHole(j) *= sunRadiusInM;
			}
			blackHole(dim) *= sunMassInM;
			
			blackHoleInfo.push_back(blackHole);
		}
		
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
			
			vector x = sim.coordForIndex(iter.index);

			//calculate gammaBar_ij
			geomCell.gammaBar_ll = eta;
			
			//calculate psi
			real psi = 1;
			real oneOverAlpha = 0;
			for (int i = 0; i < (int)blackHoleInfo.size(); ++i) {
				real M = blackHoleInfo[i](dim);
				vector c;
				for (int j = 0; j < dim; ++j) {
					c(j) = blackHoleInfo[i](j);
				}
				real r = vector::length(x - c);
				psi += .5 * M / r;
				oneOverAlpha += .5 * M / r;
			}

			//now Baumgarte & Shapiro p.59 starts off talking about Schwarzschild geometry and represents it analogous to the isotropic coordinates on p.50.
			//The isotropic coordinates show that beta and K are all zero.
			//gamma_ij is indeed (1 + M/(2r))^4 eta_ij coinciding with a conformal metric of gammaBar_ij = eta_ij and psi = 1 + M/(2r).
			//The lapse is given (for a single body, p.50) as alpha = (1 - M/(2r))/(1 + M/(2r)).
			//What about when we have multiple bodies?
			//Alcubierre calls this data "Brill-Lindquist data" and likewise doesn't mention the lapse value any more than Baumgarte & Shapiro do
			// except for one additional statement: alpha psi = 1 - M / (2 r), which he goes on to state implies alpha = (1 - M/(2r))/(1 + M/(2r)).
			//That and I skimmed through ch.12 to find something like this, and it's probably wrong:
			//geomCell.alpha = 1. / oneOverAlpha;
			geomCell.alpha = 1.;
			geomCell.phi = log(psi);
		}
	}
};

