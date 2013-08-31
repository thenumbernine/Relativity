#include "admformalism.h"
#include "integrators.h"

#include <fstream>
#include <iostream>

using namespace std;

namespace Test {

typedef double real;
enum { res = 100 };
enum { iters = 100 };

typedef RK4Integrator Integrator;

//universal constants
const real speedOfLightInMPerS = 299792458.;
const real gravitationalConstantInM3PerKgS2 = 6.67384e-11;

//conversions to meters
const real metersPerS = speedOfLightInMPerS;
const real metersPerKg = gravitationalConstantInM3PerKgS2 / (metersPerS * metersPerS);

//properties of the sun
const real sunMassInKg = 1.989e+30;
const real sunMassInM = sunMassInKg * metersPerKg;
const real sunRadiusInM = 6.955e+8;
const real sunVolumeInM3 = 4. / 3. * M_PI * sunRadiusInM * sunRadiusInM * sunRadiusInM;
const real sunDensityInM_2 = sunMassInKg * metersPerKg / sunVolumeInM3;

const real sunRotationInRad_S = 2.8e-6;
const real sunRotationInM = sunRotationInRad_S / metersPerS;
const real sunAngularMomentumInKgM2_S = .4 * sunMassInKg * sunRadiusInM * sunRadiusInM * sunRotationInRad_S;
const real sunAngularMomentumInM = .4 * sunMassInM * sunRadiusInM * sunRadiusInM * sunRotationInM;

template<int dim>
struct RunClass;

template<>
struct RunClass<1> {
	//1D case splot's all time slices together
	void operator()(::ADMFormalism<real, 1, Integrator> *sim, ostream &f, int numIters, bool outputHistory) {
		if (outputHistory) {
			sim->outputLine(f);
		}

		cout << "iterating..." << endl;
		const real dt = .1 * sim->dx(0);
		for (int i = 0; i < numIters; ++i) {
			sim->update(dt);
			
			if (outputHistory) {
				sim->outputLine(f);
			}
		}

		if (!outputHistory) {
			sim->outputLine(f);
		}
	}
};

template<>
struct RunClass<2> {
	//2D case splots the last one
	void operator()(::ADMFormalism<real, 2, Integrator> *sim, ostream &f, int numIters, bool outputHistory) {
		cout << "iterating..." << endl;
		const real dt = .01;
		for (int i = 0; i < numIters; ++i) {
			sim->update(dt);
		}
		
		sim->outputLine(f);
	}
};

template<int dim>
struct Base {
	typedef ::vector<real, dim> vector;
	typedef ::ADMFormalism<real, dim, Integrator> ADMFormalism;
	typedef typename ADMFormalism::DerefType DerefType;

	real maxDist;
	vector min, max, center;
	DerefType res;
	ADMFormalism *sim;
	int numIters;

	bool outputHistory;	//only used for 1D case

	Base(real maxDist_, int res_, int numIters_) 
	: 	maxDist(maxDist_), 
		min(-maxDist_),
		max(maxDist_),
		res(res_),
		sim(NULL),
		numIters(numIters_),
		outputHistory(true)
	{
		center = (max + min) * .5;
	}

	virtual const string filename() const = 0;

	virtual void init() {
		cout << "constructing sim..." << endl;
		sim = new ADMFormalism(min, max, res);
	}

	//update
	virtual void run() {
		ofstream f(filename().c_str());
		sim->outputHeaders(f);

		Test::RunClass<dim>()(sim, f, numIters, outputHistory);
		
		cout << "done!" << endl;
	}
};

//this is a Schwarzschild init technically, with some matter thrown in there
template<int dim>
struct Sun : public Base<dim> {
	typedef Test::Base<dim> Base;
	typedef typename Base::ADMFormalism ADMFormalism;
	typedef typename Base::vector vector;
	
	Sun() : Base(2. * sunRadiusInM, Test::res, Test::iters) {}

	virtual const string filename() const { return "sun.txt"; }

	virtual void init() {
		Base::init();
	
		ADMFormalism* &sim = Base::sim;
		vector &center = Base::center;
		
		//provide initial conditions
		cout << "providing initial conditions..." << endl;
		for (typename ADMFormalism::GridIter iter = sim->readCells->begin(); iter != sim->readCells->end(); ++iter) {
			typename ADMFormalism::Cell &cell = *iter;
			vector v = sim->coordForIndex(iter.index) - center;
			real r = vector::length(v);
			real sunMassOver2Radius = sunMassInM / (2. * r);
			real oneMinusSunMassOver2Radius = 1. - sunMassOver2Radius;
			real onePlusSunMassOver2Radius = 1. + sunMassOver2Radius;
			real onePlusSunMassOver2RadiusSq = onePlusSunMassOver2Radius * onePlusSunMassOver2Radius;
			
			cell.alpha = oneMinusSunMassOver2Radius / onePlusSunMassOver2Radius; 
		
			//beta^i = 0

			for (int i = 0; i < dim; ++i) {
				cell.gamma_ll(i,i) = onePlusSunMassOver2RadiusSq * onePlusSunMassOver2RadiusSq;
			}
			
			//gamma = det(gamma_ij)
			//ln_sqrt_gamma := ln(sqrt(gamma))
			cell.calc_ln_sqrt_gamma_from_gamma_ll();

			//K_ij = K = 0
			
			if (r <= sunRadiusInM) {
				cell.rho = sim->dx.volume() * sunDensityInM_2;
			}
		}
	}
};

/*
See the Kerr-Schild section of the scratch paper in the README
*/
template<int dim>
struct KerrSchild : public Base<dim> {
	typedef Test::Base<dim> Base;

	typedef typename Base::vector vector;
	typedef typename Base::ADMFormalism ADMFormalism;
	typedef typename ADMFormalism::tensor_l tensor_l;
	typedef typename ADMFormalism::tensor_u tensor_u;
	typedef typename ADMFormalism::tensor_sl tensor_sl;
	typedef typename ADMFormalism::tensor_su tensor_su;

	real M;		//black hole mass <=> half the Schwarzschild radius
	real Q;		//total charge
	real J;		//total angular momentum
	
	KerrSchild(
		real R_,	//half width of each dimension in the simulation
		real M_,	//mass of black hole
		real J_, 	//total angular momentum of black hole
		real Q_)	//total charge of black hole
	: Base(R_, Test::res, Test::iters),
		M(M_),
		J(J_),
		Q(Q_)
	{}

	virtual const string filename() const { return "black_hole.txt"; }

	virtual void init() {
		Base::init();
		
		real a = J / M;	//angular momentum density
		
		tensor_sl eta;
		for (int i = 0; i < dim; ++i) {
			eta(i,i) = 1.;
		}

		vector &min = Base::min;
		vector &max = Base::max;
		ADMFormalism* &sim = Base::sim;

		//provide initial conditions
		
		cout << "providing initial conditions..." << endl;
		vector center = (max + min) * .5;
		for (typename ADMFormalism::GridIter iter = sim->readCells->begin(); iter != sim->readCells->end(); ++iter) {
			typename ADMFormalism::Cell &cell = *iter;
			vector v = sim->coordForIndex(iter.index) - center;
			real r = vector::length(v);
			real x = v(0);
			real y = dim > 1 ? v(1) : 0;
			real z = dim > 2 ? v(2) : 0;
			real H = (r * M - Q * Q / 2.) / (r * r + a * a * z * z / (r * r));
			
			tensor_l l_l;
			l_l(0) = (r * x + a * y) / (r * r + a * a);
			if (dim > 1) l_l(1) = (r * y - a * x) / (r * r + a * a);
			if (dim > 2) l_l(2) = z / r;

			tensor_sl &gamma_ll = cell.gamma_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					gamma_ll(i,j) = eta(i,j) + (2. - eta(i,j)) * H * l_l(i) * l_l(j);
				}
			}
		
			//ln_sqrt_gamma := ln(sqrt(det(gamma_ij)))
			cell.calc_ln_sqrt_gamma_from_gamma_ll();
			
			//ln_psi := ln(psi) = 1/6 ln(sqrt(gamma))
			//psi = exp(ln(psi))
			cell.calc_psi_from_ln_sqrt_gamma();
			
			//gammaBar_ij = psi^-4 gamma_ij
			//gammaBar^ij = inverse(gammaBar_ij)
			cell.calc_gammaBar_uu_and_gammaBar_ll_from_psi();

			//gamma^ij = psi^-4 gammaBar^ij
			cell.calc_gamma_uu_from_gammaBar_uu_and_psi();

			tensor_su &gamma_uu = cell.gamma_uu;

			tensor_u l_u;
			for (int i = 0; i < dim; ++i) {
				l_u(i) = 0.;
				for (int j = 0; j < dim; ++j) {
					l_u(i) += gamma_uu(i,j) * l_l(j);
				}
			}

			tensor_u &beta_u = cell.beta_u;
			beta_u = l_u * H;

			real betaNorm = 0.;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					betaNorm += gamma_ll(i,j) * beta_u(i) * beta_u(j);
				}
			}

			real &alpha = cell.alpha;
			alpha = sqrt(1. - 2. * H - betaNorm);

			tensor_sl &K_ll = cell.K_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					K_ll(i,j) = 2. * H * a / r * (eta(i,j) - (2. + H) * l_l(i) * l_l(j));
				}
			}

			real &K = cell.K;
			K = 2 * M * alpha * alpha * alpha / (r * r) * (1. + 3. * M / r);
		}
	}
};


/*
This one's coming from "Numerical Relativity" p.61:
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
template<int dim, int numBlackHoles>
struct BrillLindquist : public Base<dim> {
	typedef Test::Base<dim> Base;

	typedef typename Base::vector vector;
	typedef typename Base::ADMFormalism ADMFormalism;
	typedef typename ADMFormalism::tensor_l tensor_l;
	typedef typename ADMFormalism::tensor_u tensor_u;
	typedef typename ADMFormalism::tensor_sl tensor_sl;
	typedef typename ADMFormalism::tensor_su tensor_su;

	//I should at least make this a structure or something
	typedef tensor<real, lower<numBlackHoles>, lower<dim+1>> BlackHoleInfo;
	BlackHoleInfo blackHoleInfo;

	BrillLindquist(
		real maxDist_,	//half width of each dimension in the simulation
		const BlackHoleInfo &blackHoleInfo_
	) :	Base(maxDist_, Test::res, Test::iters),
		blackHoleInfo(blackHoleInfo_)
	{}

	virtual const string filename() const { return "multiple_black_holes.txt"; }

	virtual void init() {
		Base::init();
		
		tensor_sl eta;
		for (int i = 0; i < dim; ++i) {
			eta(i,i) = 1.;
		}

		vector &min = Base::min;
		vector &max = Base::max;
		ADMFormalism* &sim = Base::sim;

		//provide initial conditions
		
		cout << "providing initial conditions..." << endl;
		vector center = (max + min) * .5;
		for (typename ADMFormalism::GridIter iter = sim->readCells->begin(); iter != sim->readCells->end(); ++iter) {
			typename ADMFormalism::Cell &cell = *iter;
			vector x = sim->coordForIndex(iter.index);

			//calculate psi and gammaBar_ij
			cell.gammaBar_ll = eta;
			cell.psi = 1;
			real oneOverAlpha = 0;
			for (int i = 0; i < numBlackHoles; ++i) {
				real M = blackHoleInfo(i, dim);
				vector c;
				for (int j = 0; j < dim; ++j) {
					c(j) = blackHoleInfo(i, j);
				}
				real r = vector::length(x - c);
				cell.psi += .5 * M / r;
				oneOverAlpha += .5 * M / r;
			}

			//now "Numerical Relativity" p.59 starts off talking about Schwarzschild geometry and represents it analogous to the isotropic coordinates on p.50.
			//The isotropic coordinates show that beta and K are all zero.
			//gamma_ij is indeed (1 + M/(2r))^4 eta_ij coinciding with a conformal metric of gammaBar_ij = eta_ij and psi = 1 + M/(2r).
			//The lapse is given (for a single body, p.50) as alpha = (1 - M/(2r))/(1 + M/(2r)).
			//What about when we have multiple bodies?
			//Alcubierre calls this data "Brill-Lindquist data" and likewise doesn't mention the lapse value any more than Baumgarte & Shapiro do
			// except for one additional statement: alpha psi = 1 - M / (2 r), which he goes on to state implies alpha = (1 - M/(2r))/(1 + M/(2r)).
			//That and I skimmed through ch.12 to find something like this, and it's probably wrong:
			//cell.alpha = 1. / oneOverAlpha;
			cell.alpha = 1.;

			//from psi we need to calculate our state values
			real psiSquared = cell.psi * cell.psi;
			real psiToTheFourth = psiSquared * psiSquared;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					cell.gamma_ll(i,j) = cell.gammaBar_ll(i,j) * psiToTheFourth;
				}
			}

			cell.gamma = psiToTheFourth * psiToTheFourth * psiToTheFourth;
			cell.ln_sqrt_gamma = .5 * log(cell.gamma);
		}
	}
};

}

int main() {
	using namespace Test;
	
	/* GRO J0422+32 : the smallest black hole yet found * /
	KerrSchild<1> test(
		4.1 * sunRadiusInM,		//simulation radius
		4.1 * sunMassInM,		//black hole mass
		0,						//angular momentum
		0						//charge
	);
	/**/
	/* Sagitarrius A* : The supermassive black hole in the center of the Milky Way * /
	KerrSchild<2> test(
		4.1e+6 * sunRadiusInM,
		4.1e+6 * sunMassInM,
		0,
		0
	);
	/**/
	/* binary black hole head on collision */
	enum { dim = 1 };
	enum { numBlackHoles = 2 };
	tensor<real, lower<numBlackHoles>, lower<dim+1>> blackHoleInfo;
	for (int i = 0; i < numBlackHoles; ++i) {
		blackHoleInfo(i,0) = (i == 0 ? -1. : 1.) * 2. * sunRadiusInM;
		blackHoleInfo(i,dim) = sunMassInM;
	}
	BrillLindquist<dim, numBlackHoles> test(4.1 * sunRadiusInM, blackHoleInfo);
	/**/

	test.outputHistory = false;
	
	test.init();
	test.run();
}

