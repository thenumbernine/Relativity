#include "admformalism.h"

#include <fstream>
#include <iostream>

using namespace std;

namespace Test {

typedef double real;
enum { res = 10 };
enum { iters = 100 };

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
	void operator()(::ADMFormalism<real, 1> *sim, ostream &f, int numIters, bool outputHistory) {
		if (outputHistory) {
			sim->outputLine(f);
		}

		cout << "iterating..." << endl;
		const real dt = .1;
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
	void operator()(::ADMFormalism<real, 2> *sim, ostream &f, int numIters) {
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
	typedef ::vector<real,dim> vector;
	typedef ::ADMFormalism<real,dim> ADMFormalism;
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
			for (int i = 0; i < dim; ++i) {
				cell.gamma_ll(i,i) = onePlusSunMassOver2RadiusSq * onePlusSunMassOver2RadiusSq;
			}
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
		
			cell.calcLnSqrtGammaFromGammaLL();
			cell.calcPsiFromLnSqrtGamma();
			cell.calcGammaBar();

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
	/* Sagitarrius A* : The supermassive black hole in the center of the Milky Way */
	KerrSchild<1> test(
		4.1e+6 * sunRadiusInM,
		4.1e+6 * sunMassInM,
		0,
		0
	);
	/**/

//	test.outputHistory = false;
	
	test.init();
	test.sim->init();
	test.run();
}

