#include "admformalism.h"

#include <fstream>
#include <iostream>

using namespace std;

namespace Test {

typedef double real;
enum { dim = 1 };
typedef ::vector<real,dim> vector;
typedef ::ADMFormalism<real,dim> ADMFormalism;
typedef ADMFormalism::deref_type deref_type;

typedef ADMFormalism::tensor_l tensor_l;
typedef ADMFormalism::tensor_u tensor_u;
typedef ADMFormalism::tensor_sl tensor_sl;
typedef ADMFormalism::tensor_su tensor_su;

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

struct Base {
	real maxDist;
	vector min, max, center;
	deref_type res;
	ADMFormalism *sim;
			
	Base(real maxDist_, int res_) 
	: maxDist(maxDist_), min(-maxDist_), max(maxDist_), res(res_), sim(NULL) {
		center = (max + min) * .5;
	}

	virtual const string filename() const = 0;

	virtual void init() {
		cout << "constructing sim..." << endl;
		sim = new ADMFormalism(min, max, res);
	}
	
	virtual void run() {
		ofstream f(filename().c_str());
		sim->outputHeaders(f);
		sim->outputLine(f);

		//update
		//cout << "iterating..." << endl;
		//const real dt = .01;
		//sim->update(dt);

		cout << "done!" << endl;
	}
};

struct Sun : public Base {
	Sun() : Base(2. * sunRadiusInM, 10) {}

	virtual const string filename() const { return "sun.txt"; }

	virtual void init() {
		Base::init();
		
		//provide initial conditions
		cout << "providing initial conditions..." << endl;
		for (ADMFormalism::GridIter iter = sim->readCells->begin(); iter != sim->readCells->end(); ++iter) {
			ADMFormalism::Cell &cell = *iter;
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
struct BlackHole : public Base {
	BlackHole() : Base(2. * sunRadiusInM, 10) {}

	virtual const string filename() const { return "black_hole.txt"; }

	virtual void init() {
		Base::init();

		const real blackHoleMassInKg = 4. * sunMassInM;	//2.7 solar masses is upper bound for neutron stars
		const real M = blackHoleMassInKg * metersPerKg;	//black hole mass in meters
		
		const real Q = 0.;		//total charge
		const real J = 0.;		//total angular momentum
		const real a = J / M;	//angular momentum density
		
		tensor_sl eta;
		for (int i = 0; i < dim; ++i) {
			eta(i,i) = 1.;
		}

		//provide initial conditions
		
		cout << "providing initial conditions..." << endl;
		vector center = (max + min) * .5;
		for (ADMFormalism::GridIter iter = sim->readCells->begin(); iter != sim->readCells->end(); ++iter) {
			ADMFormalism::Cell &cell = *iter;
			vector v = sim->coordForIndex(iter.index) - center;
			real r = vector::length(v); 
			real x = v(0);
			real y = v(1);
			real z = v(2);
			real H = (r * M - Q * Q / 2.) / (r * r + a * a * z * z / (r * r));
			
			tensor_l l_l;
			l_l(0) = (r * x + a * y) / (r * r + a * a);
			l_l(1) = (r * y - a * x) / (r * r + a * a);
			l_l(2) = z / r;
			
			tensor_sl &gamma_ll = cell.gamma_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					gamma_ll(i,j) = eta(i,j) + (2. - eta(i,j)) * H * l_l(i) * l_l(j);
				}
			}
			
			tensor_su &gamma_uu = cell.gamma_uu;
			gamma_uu = invert(gamma_ll);

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
			alpha = sqrt(betaNorm - 1. - 2. * H);

			tensor_sl &K_ll = cell.K_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					K_ll(i,j) = 2. * H * a / r * (eta(i,j) - (2. + H) * l_l(i) * l_l(j));
				}
			}
			
			real sunMassOver2Radius = sunMassInM / (2. * sunRadiusInM);
			real oneMinusSunMassOver2Radius = 1. - sunMassOver2Radius;
			real onePlusSunMassOver2Radius = 1. + sunMassOver2Radius;
			real onePlusSunMassOver2RadiusSq = onePlusSunMassOver2Radius * onePlusSunMassOver2Radius;
			cell.alpha = oneMinusSunMassOver2Radius / onePlusSunMassOver2Radius; 
			for (int i = 0; i < dim; ++i) {
				cell.gamma_ll(i,i) = onePlusSunMassOver2RadiusSq * onePlusSunMassOver2RadiusSq;
			}
		}
	}
};

}

int main() {
	typedef ::Test::BlackHole Test;
	Test test;
	test.init();
	test.run();
}

