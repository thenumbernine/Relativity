#pragma once

#include "initialdata/initialdata.h"
#include "constants.h"
#include "parallel.h"
#include "Common/Exception.h"
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
template<typename Real, int dim>
struct BrillLindquist : public InitialData<Real, dim> {
	virtual const char *name() { return "brill-lindquist"; }

	using Vector = Tensor::vec<Real, dim>;
	using DerefType = Tensor::intN<dim>;
	using InitialData = ::InitialData<Real, dim>;
	using ADMFormalism = typename InitialData::ADMFormalism;
	using TensorL = typename ADMFormalism::TensorL;
	using TensorU = typename ADMFormalism::TensorU;
	using TensorSL = typename ADMFormalism::TensorSL;
	using TensorSU = typename ADMFormalism::TensorSU;

	//I should at least make this a structure or something
	std::vector<Tensor::vec<Real, dim+1>> blackHoleInfo;

	virtual void init(ADMFormalism &sim, std::vector<std::string> &args) {
		if (!args.size()) throw Common::Exception() << "expected simulation arguments";
		int numBlackHoles = atoi(args[0].c_str());
		args.erase(args.begin());

		std::cout << "number of black holes " << numBlackHoles << std::endl;

		for (int i = 0; i < numBlackHoles; ++i) {
			Tensor::vec<Real, dim+1> blackHole;
			for (int j = 0; j < dim+1; ++j) {
				if (!args.size()) throw Common::Exception() << "expected simulation arguments";
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
		
		TensorSL eta;
		for (int i = 0; i < dim; ++i) {
			eta(i,i) = 1.;
		}

		//provide initial conditions
		
		std::cout << "providing initial conditions..." << std::endl;
		Tensor::RangeObj<dim> range = sim.geomGridReadCurrent->range();
		parallel.foreach(range.begin(), range.end(), [&](Tensor::intN<dim> index) {
			typename ADMFormalism::GeomCell &geomCell = (*sim.geomGridReadCurrent)(index); 
			
			Vector x = sim.coordForIndex(index);

			//calculate gammaBar_ij
			geomCell.gammaBar_ll = eta;
			
			//calculate psi
			Real psi = 1;
			Real oneOverAlpha = 0;
			for (int i = 0; i < (int)blackHoleInfo.size(); ++i) {
				Real M = blackHoleInfo[i](dim);
				Vector c;
				for (int j = 0; j < dim; ++j) {
					c(j) = blackHoleInfo[i](j);
				}
				Real r = (x - c).length();
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
		});
	}
};
