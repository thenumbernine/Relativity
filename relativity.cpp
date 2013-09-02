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

#include "init_schwarzschild.h"
#include "init_kerr_schild.h"
#include "init_brill_lindquist.h"
#include "init_bowen_york.h"

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

//N-D case splots the last slice 
template<typename real, int dim>
struct RunClass {
	void operator()(::ADMFormalism<real, dim> *sim, ostream &f, int numIters, bool outputHistory) {
		cout << "iterating..." << endl;
		const real dt = .1 * sim->dx(0);
		for (int i = 0; i < numIters; ++i) {
			sim->update(dt);
		}
		
		sim->outputLine(f);
	}
};

//1D case splot's all time slices together, or just one slice
template<typename real>
struct RunClass<real, 1> {
	void operator()(::ADMFormalism<real, 1> *sim, ostream &f, int numIters, bool outputHistory) {
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

template<typename real, int dim>
void runSimIntegrator(SimParams &params, IIntegrator<real, dim> *integrator) {

	cout << "constructing sim..." << endl;
	
	real maxDist = params.size * sunRadiusInM;
	typedef ::vector<real, dim> vector;
	ADMFormalism<real, dim> sim(vector(-maxDist), vector(maxDist), params.res, integrator);

	//why not allow for a non-initial-condition sim?
	//for the record, i think alpha will be initialized to zero as welll, 
	//so our 4D metrics technically will be singular ...
	if (params.args.size()) {
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

			KerrSchild<real, dim> test(mass * sunMassInM, angularMomentum, charge);
			cout << "providing initial conditions..." << endl;
			test.init(sim);

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

			BrillLindquist<real, dim> test(blackHoleInfos);
			cout << "providing initial conditions..." << endl;
			test.init(sim);
		} else {
			throw Exception() << "got unknown simulation type " << simType;
		}
	}

	ofstream f(params.filename.c_str());
	
	sim.outputHeaders(f);
	
	RunClass<real, dim>()(&sim, f, params.iter, params.history);	
	
	cout << "done!" << endl;
}

template<typename real, int dim>
void runSimPrecision(SimParams &params) {
	
	IIntegrator<real, dim> *integrator = NULL;
	switch (params.integrator) {
	case INTEGRATOR_EULER:
		integrator = new EulerIntegrator<real, dim>();
		break;
	case INTEGRATOR_RK4:
		integrator = new RK4Integrator<real, dim>();
		break;
	default:
		throw Exception() << "got an integrator I couldn't handle " << params.integrator;
	}
	runSimIntegrator<real, dim>(params, integrator);
	delete integrator;
}

template<int dim>
void runSimDim(SimParams &params) {
	switch (params.precision) {
	case PRECISION_FLOAT:
		runSimPrecision<float, dim>(params);
		break;
	case PRECISION_DOUBLE:
		runSimPrecision<double, dim>(params);
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

