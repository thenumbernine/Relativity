#include "admformalism.h"
#include "integrators.h"
#include "Common/Exception.h"

#include <string.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "constants.h"
#include "initialdata/bowen_york.h"
#include "initialdata/brill_lindquist.h"
#include "initialdata/kerr_schild.h"
#include "initialdata/schwarzschild.h"

#include "output_table.h"

#define DISABLE_FLOAT	//disable 32-bit fpp for faster builds

using namespace std;

template<typename real, int dim>
struct RunTest {
	void operator()(
		::ADMFormalism<real, dim> &sim,
		ostream &f,
		real cfl, 
		int numIters, 
		bool outputHistory, 
		std::vector<bool> &cols)
	{
		const real dt = cfl * sim.dx(0);
		
		if (outputHistory) {
			sim.calcAux(*sim.getGeomGridReadCurrent());	//in case we're only outputting zero iterations
			//sim.outputLine(f);
			OutputTable<real, dim>::state(f, sim, cols);
		}

		cout << "iterating..." << endl;
		for (int i = 0; i < numIters; ++i) {
			sim.update(dt);
			
			if (outputHistory) {
				sim.calcAux(*sim.getGeomGridReadCurrent());	//in case we're only outputting zero iterations
				//sim.outputLine(f);
				OutputTable<real, dim>::state(f, sim, cols);
			}
		}

		if (!outputHistory) {
			sim.calcAux(*sim.getGeomGridReadCurrent());	//in case we're only outputting zero iterations
			//sim.outputLine(f);
			OutputTable<real, dim>::state(f, sim, cols);
		}
	}
};

template<typename real, int dim, typename... InitialDataTypes>
struct InitTest;

template<typename real, int dim, typename InitialDataType, typename... InitialDataTypes>
struct InitTest<real, dim, InitialDataType, InitialDataTypes...> {
	typedef InitTest<real, dim, InitialDataTypes...> NextType;
	static void init(const std::string &simType, std::vector<std::string> &args, ADMFormalism<real, dim> &sim) {
		InitialDataType init;
		if (simType == init.name()) {
			init.init(sim, args);
		} else {
			NextType::init(simType, args, sim);
		}
	}
};

template<typename real, int dim, typename InitialDataType>
struct InitTest<real, dim, InitialDataType>{
	static void init(const std::string &simType, std::vector<std::string> &args, ADMFormalism<real, dim> &sim) {
		InitialDataType init;
		if (simType == init.name()) {
			init.init(sim, args);
		} else {
			throw Common::Exception() << "got unknown simulation type " << simType;
		}
	}
};

enum {
#ifndef DISABLE_FLOAT		
	PRECISION_FLOAT,
#endif
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
		cfl(.1),
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

	//cfl <n>	= CFL
	double cfl;

	//precision <prec> = precision (see PRECISION_* enum) 
	int precision;

	//integrator <int>	= integrator (see INTEGRATOR_* enum)
	int integrator;

	//history	= include history in output
	bool history;

	//filename	= output filename
	std::string filename;

	//desired columns
	std::string columnNames;

	//rest of args
	std::vector<std::string> args;
};

SimParams interpretArgs(int argc, char **argv) {
	SimParams params;
	for (int i = 1; i < argc; ++i) {
		//0-param vars
		if (!strcmp(argv[i], "history")) {
			params.history = true;
			cout << "using history" << endl;
			continue;
		} else if (!strcmp(argv[i], "allcols")) {
			params.columnNames = "*all*";
			continue;
		//1-param vars
		} else if (i < argc-1) {
			if (!strcmp(argv[i], "dim")) {
				params.dim = atoi(argv[++i]);
				cout << "dim " << params.dim << endl;
				continue;
			} else if (!strcmp(argv[i], "res")) {
				params.res = atoi(argv[++i]);
				cout << "res " << params.res << endl;
				continue;
			} else if (!strcmp(argv[i], "size")) {
				params.size = atof(argv[++i]);
				cout << "size " << params.size << endl;
				continue;
			} else if (!strcmp(argv[i], "iter")) {
				params.iter = atoi(argv[++i]);
				cout << "iter " << params.iter << endl;
				continue;
			} else if (!strcmp(argv[i], "cfl")) {
				params.cfl = atof(argv[++i]);
				cout << "cfl " << params.cfl << endl;
				continue;
			} else if (!strcmp(argv[i], "precision")) {
				const char *precision = argv[++i];
#ifndef DISABLE_FLOAT
				if (!strcmp(precision, "float")) {
					params.precision = PRECISION_FLOAT;
				} else 
#endif				
				if (!strcmp(precision, "double")) {
					params.precision = PRECISION_DOUBLE;
				} else {
					throw Common::Exception() << "got an unknown precision " << precision;
				}
				continue;
			} else if (!strcmp(argv[i], "integrator")) {
				const char *integrator = argv[++i];
				if (!strcmp(integrator, "euler")) {
					params.integrator = INTEGRATOR_EULER;
				} else if (!strcmp(integrator, "rk4")) {
					params.integrator = INTEGRATOR_RK4;
				} else {
					throw Common::Exception() << "got an unknown integrator " << integrator;
				}
				continue;
			} else if (!strcmp(argv[i], "filename")) {
				params.filename = argv[++i];
				continue;
			} else if (!strcmp(argv[i], "cols")) {
				params.columnNames = argv[++i];
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

	std::vector<bool> cols(OutputTable<real, dim>::getNumColumns());
	//determine columns to output

	if (params.columnNames == "*all*") {
		for (int i = 0; i < (int)cols.size(); ++i) {
			cols[i] = true;
		}
	} else {
		std::string restOfColumnNames = params.columnNames;
		bool anyFound = false;
		while (restOfColumnNames.length() > 0) {
			std::cout << "rest of names " << restOfColumnNames << std::endl;
			size_t delimPos = restOfColumnNames.find(",");
			std::cout << "delim pos " << delimPos << std::endl;
			std::string colName;
			if (delimPos == std::string::npos) {
				colName = restOfColumnNames;
				restOfColumnNames = "";
			} else {
				colName = restOfColumnNames.substr(0, delimPos);
				restOfColumnNames = restOfColumnNames.substr(delimPos + 1);
			}
			std::cout << "col name " << colName << std::endl;
			if (colName.length()) {
				anyFound = true;
				cols[OutputTable<real, dim>::getColumnIndex(colName)] = true;
			}
		}
		if (!anyFound) throw Common::Exception() << "you forgot to provide any columns to output";
	}


	cout << "constructing sim..." << endl;
	
	real maxDist = params.size * sunRadiusInM;
	typedef Tensor::Vector<real, dim> Vector;
	ADMFormalism<real, dim> sim(Vector(-maxDist), Vector(maxDist), params.res, integrator);

	//why not allow for a non-initial-condition sim?
	//for the record, i think alpha will be initialized to zero as welll, 
	//so our 4D metrics technically will be singular ...
	if (params.args.size()) {
		string simType = params.args[0];
		params.args.erase(params.args.begin());

		InitTest<real, dim,
			BowenYork<real, dim>,
			BrillLindquist<real, dim>,
			KerrSchild<real, dim>,
			Schwarzschild<real, dim>
		>::init(simType, params.args, sim);
	}

	//construct connBar^i based on gammaBar_ij
	sim.calcConnBar(*sim.getGeomGridReadCurrent());

	ofstream f(params.filename.c_str());

	OutputTable<real, dim>::header(f, cols);
	//sim.outputHeaders(f);
	
	RunTest<real, dim>()(sim, f, params.cfl, params.iter, params.history, cols);
	
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
		throw Common::Exception() << "got an integrator I couldn't handle " << params.integrator;
	}
	runSimIntegrator<real, dim>(params, integrator);
	delete integrator;
}

template<int dim>
void runSimDim(SimParams &params) {
	switch (params.precision) {
#ifndef DISABLE_FLOAT
	case PRECISION_FLOAT:
		runSimPrecision<float, dim>(params);
		break;
#endif
	case PRECISION_DOUBLE:
		runSimPrecision<double, dim>(params);
		break;
	default:
		throw Common::Exception() << "got a precision I couldn't handle " << params.precision;
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
		throw Common::Exception() << "got a dimension I couldn't handle " << params.dim;
	}
}

int main(int argc, char **argv) {
	SimParams params = interpretArgs(argc, argv);
	runSim(params);
}
