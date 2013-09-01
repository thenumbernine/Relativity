#include "admformalism.h"
#include "integrators.h"
#include "inverse.h"
#include "exception.h"

#include <string.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

//universal constants
const double speedOfLightInMPerS = 299792458.;
const double gravitationalConstantInM3PerKgS2 = 6.67384e-11;

//conversions to meters
const double metersPerS = speedOfLightInMPerS;
const double metersPerKg = gravitationalConstantInM3PerKgS2 / (metersPerS * metersPerS);

//properties of the sun
const double sunMassInKg = 1.989e+30;
const double sunMassInM = sunMassInKg * metersPerKg;
const double sunRadiusInM = 6.955e+8;
const double sunVolumeInM3 = 4. / 3. * M_PI * sunRadiusInM * sunRadiusInM * sunRadiusInM;
const double sunDensityInM_2 = sunMassInKg * metersPerKg / sunVolumeInM3;

const double sunRotationInRad_S = 2.8e-6;
const double sunRotationInM = sunRotationInRad_S / metersPerS;
const double sunAngularMomentumInKgM2_S = .4 * sunMassInKg * sunRadiusInM * sunRadiusInM * sunRotationInRad_S;
const double sunAngularMomentumInM = .4 * sunMassInM * sunRadiusInM * sunRadiusInM * sunRotationInM;

using namespace std;

//typedef double real;
//enum { iters = 100 };
//typedef RK4Integrator Integrator;

//N-D case splots the last slice 
template<int dim, typename real, typename Integrator>
struct RunClass {
	void operator()(::ADMFormalism<real, dim, Integrator> *sim, ostream &f, int numIters, bool outputHistory) {
		cout << "iterating..." << endl;
		const real dt = .1 * sim->dx(0);
		for (int i = 0; i < numIters; ++i) {
			sim->update(dt);
		}
		
		sim->outputLine(f);
	}
};

//1D case splot's all time slices together, or just one slice
template<typename real, typename Integrator>
struct RunClass<1, real, Integrator> {
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

template<int dim, typename real, typename Integrator>
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

	virtual void init() {
		cout << "constructing sim..." << endl;
		sim = new ADMFormalism(min, max, res);
	}

	//update
	virtual void run(ostream &f) {
		sim->outputHeaders(f);

		RunClass<dim, real, Integrator>()(sim, f, numIters, outputHistory);
		
		cout << "done!" << endl;
	}
};

//this is a Schwarzschild init technically, with some matter thrown in there
template<int dim, typename real, typename Integrator>
struct Sun : public Base<dim, real, Integrator> {
	typedef ::Base<dim, real, Integrator> Base;
	typedef typename Base::ADMFormalism ADMFormalism;
	typedef typename Base::vector vector;
	
	Sun(int res_, int iters_) : Base(2. * sunRadiusInM, res_, iters_) {}

	virtual void init() {
		Base::init();
	
		ADMFormalism* &sim = Base::sim;
		vector &center = Base::center;
		
		//provide initial conditions
		cout << "providing initial conditions..." << endl;
		for (typename ADMFormalism::GeomGrid::iterator iter = sim->geomGridReadCurrent->begin(); iter != sim->geomGridReadCurrent->end(); ++iter) {
			typename ADMFormalism::GeomCell &geomCell = *iter;
			typename ADMFormalism::MatterCell &matterCell = sim->matterGrid(iter.index);

			vector v = sim->coordForIndex(iter.index) - center;
			
			real r = vector::length(v);
			real sunMassOver2Radius = sunMassInM / (2. * r);
			real oneMinusSunMassOver2Radius = 1. - sunMassOver2Radius;
			real onePlusSunMassOver2Radius = 1. + sunMassOver2Radius;
			real onePlusSunMassOver2RadiusSq = onePlusSunMassOver2Radius * onePlusSunMassOver2Radius;
			
			geomCell.alpha = oneMinusSunMassOver2Radius / onePlusSunMassOver2Radius; 
		
			//beta^i = 0

			for (int i = 0; i < dim; ++i) {
				geomCell.gamma_ll(i,i) = onePlusSunMassOver2RadiusSq * onePlusSunMassOver2RadiusSq;
			}
			
			//gamma = det(gamma_ij)
			//ln_sqrt_gamma := ln(sqrt(gamma))
			geomCell.calc_ln_sqrt_gamma_from_gamma_ll();

			//K_ij = K = 0
			
			if (r <= sunRadiusInM) {
				matterCell.rho = sim->dx.volume() * sunDensityInM_2;
			}
		}
	}
};

/*
See the Kerr-Schild section of the scratch paper in the README
*/
template<int dim, typename real, typename Integrator>
struct KerrSchild : public Base<dim, real, Integrator> {
	typedef ::Base<dim, real, Integrator> Base;

	typedef typename Base::vector vector;
	typedef typename Base::ADMFormalism ADMFormalism;
	typedef typename ADMFormalism::tensor_l tensor_l;
	typedef typename ADMFormalism::tensor_u tensor_u;
	typedef typename ADMFormalism::tensor_sl tensor_sl;
	typedef typename ADMFormalism::tensor_su tensor_su;

	real M;		//black hole mass <=> half the Schwarzschild radius
	real J;		//total angular momentum
	real Q;		//total charge
	
	KerrSchild(
		int res_,	//resolution
		int iters_,	//iterations
		real R_,	//half width of each dimension in the simulation
		real M_,	//mass of black hole
		real J_, 	//total angular momentum of black hole
		real Q_)	//total charge of black hole
	: Base(R_, res_, iters_),
		M(M_),
		J(J_),
		Q(Q_)
	{}

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
		for (typename ADMFormalism::GeomGrid::iterator iter = sim->geomGridReadCurrent->begin(); iter != sim->geomGridReadCurrent->end(); ++iter) {
			typename ADMFormalism::GeomCell &geomCell = *iter;
			
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

			tensor_sl &gamma_ll = geomCell.gamma_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					gamma_ll(i,j) = eta(i,j) + (2. - eta(i,j)) * H * l_l(i) * l_l(j);
				}
			}
		
			//ln_sqrt_gamma := ln(sqrt(det(gamma_ij)))
			geomCell.calc_ln_sqrt_gamma_from_gamma_ll();
		
			tensor_su gamma_uu;
			gamma_uu = inverse(geomCell.gamma_ll);

			tensor_u l_u;
			for (int i = 0; i < dim; ++i) {
				l_u(i) = 0.;
				for (int j = 0; j < dim; ++j) {
					l_u(i) += gamma_uu(i,j) * l_l(j);
				}
			}

			tensor_u &beta_u = geomCell.beta_u;
			beta_u = l_u * H;

			real betaNorm = 0.;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					betaNorm += geomCell.gamma_ll(i,j) * geomCell.beta_u(i) * geomCell.beta_u(j);
				}
			}

			real &alpha = geomCell.alpha;
			alpha = sqrt(1. - 2. * H - betaNorm);

			tensor_sl &K_ll = geomCell.K_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					K_ll(i,j) = 2. * H * a / r * (eta(i,j) - (2. + H) * l_l(i) * l_l(j));
				}
			}

			real &K = geomCell.K;
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
template<int dim, typename real, typename Integrator>
struct BrillLindquist : public Base<dim, real, Integrator> {
	typedef ::Base<dim, real, Integrator> Base;

	typedef typename Base::vector vector;
	typedef typename Base::ADMFormalism ADMFormalism;
	typedef typename ADMFormalism::tensor_l tensor_l;
	typedef typename ADMFormalism::tensor_u tensor_u;
	typedef typename ADMFormalism::tensor_sl tensor_sl;
	typedef typename ADMFormalism::tensor_su tensor_su;

	//I should at least make this a structure or something
	std::vector<tensor<real, lower<dim+1>>> blackHoleInfo;

	BrillLindquist(
		int res_,
		int iters_,
		real maxDist_,	//half width of each dimension in the simulation
		const std::vector<tensor<real, lower<dim+1>>> &blackHoleInfo_
	) :	Base(maxDist_, res_, iters_),
		blackHoleInfo(blackHoleInfo_)
	{}

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
		for (typename ADMFormalism::GeomGrid::iterator iter = sim->geomGridReadCurrent->begin(); iter != sim->geomGridReadCurrent->end(); ++iter) {
			typename ADMFormalism::GeomCell &geomCell = *iter;
			
			vector x = sim->coordForIndex(iter.index);

			//calculate psi and gammaBar_ij
			tensor_sl gammaBar_ll = eta;
			
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

			//now "Numerical Relativity" p.59 starts off talking about Schwarzschild geometry and represents it analogous to the isotropic coordinates on p.50.
			//The isotropic coordinates show that beta and K are all zero.
			//gamma_ij is indeed (1 + M/(2r))^4 eta_ij coinciding with a conformal metric of gammaBar_ij = eta_ij and psi = 1 + M/(2r).
			//The lapse is given (for a single body, p.50) as alpha = (1 - M/(2r))/(1 + M/(2r)).
			//What about when we have multiple bodies?
			//Alcubierre calls this data "Brill-Lindquist data" and likewise doesn't mention the lapse value any more than Baumgarte & Shapiro do
			// except for one additional statement: alpha psi = 1 - M / (2 r), which he goes on to state implies alpha = (1 - M/(2r))/(1 + M/(2r)).
			//That and I skimmed through ch.12 to find something like this, and it's probably wrong:
			//geomCell.alpha = 1. / oneOverAlpha;
			geomCell.alpha = 1.;

			//from psi we need to calculate our state values
			real psiSquared = psi * psi;
			real psiToTheFourth = psiSquared * psiSquared;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					geomCell.gamma_ll(i,j) = gammaBar_ll(i,j) * psiToTheFourth;
				}
			}

			real gamma = psiToTheFourth * psiToTheFourth * psiToTheFourth;
			geomCell.ln_sqrt_gamma = .5 * log(gamma);
		}
	}
};

enum {
	PRECISION_FLOAT,
	PRECISION_DOUBLE,
	NUM_PRECISIONS
};

enum {
	INTEGRATOR_EULER,
	INTEGRATOR_RK4,
	NUM_INTEGRATORS
};

struct SimParams {
	SimParams()
	:	dim(1),
		res(100),
		size(1),
		iter(100),
		precision(PRECISION_DOUBLE),
		integrator(INTEGRATOR_RK4),
		history(false),
		filename("out.txt")
	{}
	
	//dim <n>	= dimension
	int dim;	
	
	//res <n>	= resolution
	int res;

	//size <n>	= radial distance (in Sun radii)
	double size;

	//iter <n>	= number of iterations
	int iter;

	//precision <prec> = precision (see PRECISION_* enum) 
	int precision;

	//integrator <int>	= integrator (see INTEGRATOR_* enum)
	int integrator;

	//history	= include history in output
	bool history;

	//filename	= output filename
	std::string filename;

	//rest of args
	std::vector<std::string> args;
};

SimParams interpretArgs(int argc, char **argv) {
	SimParams params;
	for (int i = 1; i < argc; ++i) {
		//0-param vars
		if (!strcmp(argv[i], "history")) {
			params.history = true;
			cerr << "using history" << endl;
			continue;
		//1-param vars
		} else if (i < argc-1) {
			if (!strcmp(argv[i], "dim")) {
				params.dim = atoi(argv[++i]);
				cerr << "dim " << params.dim << endl;
				continue;
			} else if (!strcmp(argv[i], "res")) {
				params.res = atoi(argv[++i]);
				cerr << "res " << params.res << endl;
				continue;
			} else if (!strcmp(argv[i], "size")) {
				params.size = atof(argv[++i]);
				cerr << "size " << params.size << endl;
				continue;
			} else if (!strcmp(argv[i], "iter")) {
				params.iter = atoi(argv[++i]);
				cerr << "iter " << params.iter << endl;
				continue;
			} else if (!strcmp(argv[i], "precision")) {
				const char *precision = argv[++i];
				if (!strcmp(precision, "float")) {
					params.precision = PRECISION_FLOAT;
				} else if (!strcmp(precision, "double")) {
					params.precision = PRECISION_DOUBLE;
				} else {
					throw Exception() << "got an unknown precision " << precision;
				}
				continue;
			} else if (!strcmp(argv[i], "integrator")) {
				const char *integrator = argv[++i];
				if (!strcmp(integrator, "euler")) {
					params.integrator = INTEGRATOR_EULER;
				} else if (!strcmp(integrator, "rk4")) {
					params.integrator = INTEGRATOR_RK4;
				} else {
					throw Exception() << "got an unknown integrator " << integrator;
				}
				continue;
			} else if (!strcmp(argv[i], "filename")) {
				params.filename = argv[++i];
				continue;
			}
		}

		//save the rest for later
		//note this could be abused: 
		// arg[0] res 1 arg[1] dim 2 gets stored as arg[0] arg[1] arg[2]
		params.args.push_back(argv[i]);
	}
	return params;
}
	

//I gotta stop templating out everything
// or at least make better use of interfaces 

template<typename TestType>
void runSimTest(const SimParams &params, TestType &test) {
	test.outputHistory = params.history;
	test.init();
	
	ofstream f(params.filename.c_str());
	test.run(f);
}

template<int dim, typename real, typename integrator>
void runSimIntegrator(SimParams &params) {
	if (!params.args.size()) throw Exception() << "expected simulation arguments";
	string simType = params.args[0];
	params.args.erase(params.args.begin());

	if (simType == "kerr-schild") {
		if (!params.args.size()) throw Exception() << "expected simulation mass";
		double mass = atof(params.args[0].c_str());
		params.args.erase(params.args.begin());
	
		double angularMomentum = 0;
		double charge = 0;
		if (params.args.size()) {
			angularMomentum = atof(params.args[0].c_str());
			params.args.erase(params.args.begin());
		
			if (params.args.size()) {
				charge = atof(params.args[0].c_str());
				params.args.erase(params.args.begin());
			}
		}

		cerr << "mass " << mass << endl;
		cerr << "angular momentum " << angularMomentum << endl;
		cerr << "charge " << charge << endl;

		KerrSchild<dim, real, integrator> test(
			params.res,
			params.iter,
			params.size * sunRadiusInM,
			mass * sunMassInM,
			angularMomentum,
			charge
		);

		runSimTest(params, test);
	} else if (simType == "brill-lindquist") {
		if (!params.args.size()) throw Exception() << "expected simulation arguments";
		int numBlackHoles = atoi(params.args[0].c_str());
		params.args.erase(params.args.begin());

		cerr << "number of black holes " << numBlackHoles << endl;

		std::vector<tensor<real, lower<dim+1>>> blackHoleInfos;
		for (int i = 0; i < numBlackHoles; ++i) {
			tensor<real, lower<dim+1>> blackHoleInfo;
			for (int j = 0; j < dim+1; ++j) {
				if (!params.args.size()) throw Exception() << "expected simulation arguments";
				blackHoleInfo(j) = atof(params.args[0].c_str());
				params.args.erase(params.args.begin());
			}
			
			cerr << "black hole position and mass " << blackHoleInfo << endl;
			
			for (int j = 0; j < dim; ++j) {
				blackHoleInfo(j) *= sunRadiusInM;
			}
			blackHoleInfo(dim) *= sunMassInM;
			
			blackHoleInfos.push_back(blackHoleInfo);
		}

		BrillLindquist<dim, real, integrator> test(
			params.res,
			params.iter,
			params.size * sunRadiusInM, 
			blackHoleInfos);
	
		runSimTest(params, test);
	} else {
		throw Exception() << "got unknown simulation type " << simType;
	}
}

template<int dim, typename real>
void runSimPrecision(SimParams &params) {
	switch (params.integrator) {
	case INTEGRATOR_EULER:
		runSimIntegrator<dim, real, EulerIntegrator>(params);
		break;
	case INTEGRATOR_RK4:
		runSimIntegrator<dim, real, RK4Integrator>(params);
		break;
	default:
		throw Exception() << "got an integrator I couldn't handle " << params.integrator;
	}
}

template<int dim>
void runSimDim(SimParams &params) {
	switch (params.precision) {
	case PRECISION_FLOAT:
		runSimPrecision<dim, float>(params);
		break;
	case PRECISION_DOUBLE:
		runSimPrecision<dim, double>(params);
		break;
	default:
		throw Exception() << "got a precision I couldn't handle " << params.precision;
	}
}

void runSim(SimParams &params) {
	switch (params.dim) {
	case 1:
		runSimDim<1>(params);
		break;
	case 2:
		runSimDim<2>(params);
		break;
	case 3:
		runSimDim<3>(params);
		break;
	default:
		throw Exception() << "got a dimension I couldn't handle " << params.dim;
	}
}

int main(int argc, char **argv) {
	SimParams params = interpretArgs(argc, argv);
	runSim(params);
}

