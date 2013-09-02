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

#include "initialdata.h"

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


//this is a Schwarzschild init technically, with some matter thrown in there
//TODO take from Baumgarte & Shapiro p.66
template<typename real, int dim>
struct Sun : public InitialData<real, dim> {
	typedef ::vector<real, dim> vector;
	typedef ::InitialData<real, dim> InitialData;
	typedef typename InitialData::ADMFormalism ADMFormalism;

	virtual void init(ADMFormalism &sim) {
		
		vector center = (sim.max + sim.min) * .5;
		
		//provide initial conditions
		cout << "providing initial conditions..." << endl;
		for (typename ADMFormalism::GeomGrid::iterator iter = sim.geomGridReadCurrent->begin(); iter != sim.geomGridReadCurrent->end(); ++iter) {
			typename ADMFormalism::GeomCell &geomCell = *iter;
			typename ADMFormalism::MatterCell &matterCell = sim.matterGrid(iter.index);

			vector v = sim.coordForIndex(iter.index) - center;
			
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
				matterCell.rho = sim.dx.volume() * sunDensityInM_2;
			}
		}
	}
};

/*
Kerr-Schild black hole
See Alcubierre p.56 and Baumgarte & Shapiro p.52
*/
template<typename real, int dim>
struct KerrSchild : public InitialData<real, dim> {
	typedef ::vector<real, dim> vector;
	typedef ::InitialData<real, dim> InitialData;
	typedef typename InitialData::ADMFormalism ADMFormalism;
	typedef typename ADMFormalism::tensor_l tensor_l;
	typedef typename ADMFormalism::tensor_u tensor_u;
	typedef typename ADMFormalism::tensor_sl tensor_sl;
	typedef typename ADMFormalism::tensor_su tensor_su;

	real M;		//black hole mass <=> half the Schwarzschild radius
	real J;		//total angular momentum
	real Q;		//total charge
	
	KerrSchild(
		real M_,	//mass of black hole
		real J_, 	//total angular momentum of black hole
		real Q_)	//total charge of black hole
	: 	M(M_),
		J(J_),
		Q(Q_)
	{}

	virtual void init(ADMFormalism &sim) {
		
		real a = J / M;	//angular momentum density
		
		tensor_sl eta;
		for (int i = 0; i < dim; ++i) {
			eta(i,i) = 1.;
		}

		const vector &min = sim.min;
		const vector &max = sim.max;

		//provide initial conditions
		
		cout << "providing initial conditions..." << endl;
		vector center = (max + min) * .5;
		for (typename ADMFormalism::GeomGrid::iterator iter = sim.geomGridReadCurrent->begin(); iter != sim.geomGridReadCurrent->end(); ++iter) {
			typename ADMFormalism::GeomCell &geomCell = *iter;
				
			vector v = sim.coordForIndex(iter.index) - center;
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
			
			real &K = geomCell.K;
			K = 2 * M * alpha * alpha * alpha / (r * r) * (1. + 3. * M / r);

			tensor_sl &K_ll = geomCell.K_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					K_ll(i,j) = 2. * H * a / r * (eta(i,j) - (2. + H) * l_l(i) * l_l(j));
				}
			}
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
template<typename real, int dim>
struct BrillLindquist : public InitialData<real, dim> {
	typedef ::vector<real, dim> vector;
	typedef ::InitialData<real, dim> InitialData;
	typedef typename InitialData::ADMFormalism ADMFormalism;
	typedef typename ADMFormalism::tensor_l tensor_l;
	typedef typename ADMFormalism::tensor_u tensor_u;
	typedef typename ADMFormalism::tensor_sl tensor_sl;
	typedef typename ADMFormalism::tensor_su tensor_su;

	//I should at least make this a structure or something
	std::vector<tensor<real, lower<dim+1>>> blackHoleInfo;

	BrillLindquist(const std::vector<tensor<real, lower<dim+1>>> &blackHoleInfo_)
	:	blackHoleInfo(blackHoleInfo_)
	{}

	virtual void init(ADMFormalism &sim) {
		
		tensor_sl eta;
		for (int i = 0; i < dim; ++i) {
			eta(i,i) = 1.;
		}

		const vector &min = sim.min;
		const vector &max = sim.max;

		//provide initial conditions
		
		cout << "providing initial conditions..." << endl;
		vector center = (max + min) * .5;
		for (typename ADMFormalism::GeomGrid::iterator iter = sim.geomGridReadCurrent->begin(); iter != sim.geomGridReadCurrent->end(); ++iter) {
			typename ADMFormalism::GeomCell &geomCell = *iter;
			
			vector x = sim.coordForIndex(iter.index);

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

/*
Bowen-York spinning black hole
See Baumgarte & Shapiro p.70 and Alcubierre p.109
*/
template<typename real, int dim>
struct BowenYork : public InitialData<real, dim> {
	typedef ::vector<real, dim> vector;
	typedef ::InitialData<real, dim> InitialData;
	typedef typename InitialData::ADMFormalism ADMFormalism;
	typedef typename ADMFormalism::tensor_l tensor_l;
	typedef typename ADMFormalism::tensor_u tensor_u;
	typedef typename ADMFormalism::tensor_sl tensor_sl;
	typedef typename ADMFormalism::tensor_su tensor_su;

	real M;			//black hole mass <=> half the Schwarzschild radius
	tensor_l J_l;	//angular momentum. technically only has to be a divergence-free field.

	BowenYork(
		real M_,		//mass of black hole
		tensor_l J_) 	//angular momentum of black hole (Alcubierre calls this one 'S', but Baumgarte & Shapiro already use 'S' for stress-energy momentum, which Alcubierre calls 'j')
	: 	M(M_),
		J_l(J_)
	{}

	virtual void init(ADMFormalism &sim) {
		
		//J = |J^i|
		real J = tensor_l::body::length(J_l.body);
	
		tensor_sl eta;
		for (int i = 0; i < dim; ++i) {
			eta(i,i) = 1.;
		}

		const vector &min = sim.min;
		const vector &max = sim.max;

		//provide initial conditions
		
		cout << "providing initial conditions..." << endl;
		vector center = (max + min) * .5;
		for (typename ADMFormalism::GeomGrid::iterator iter = sim.geomGridReadCurrent->begin(); iter != sim.geomGridReadCurrent->end(); ++iter) {
			typename ADMFormalism::GeomCell &geomCell = *iter;
			//we don't need this cell, just a function in the AuxCell class for computing psi from ln_sqrt_gamma
			typename ADMFormalism::AuxCell &cell = sim.auxGrid(iter.index);
				
			vector x = sim.coordForIndex(iter.index) - center;
			
			//r = |x^i|
			real r = vector::length(x);
			real rSquared = r * r;
			real rCubed = r * rSquared;

			//l^i = x^i / r
			tensor_u l_u;
			for (int i = 0; i < dim; ++i) {
				l_u(i) = x(i) / r; 
			}

			real MOverTwoR = M / (2. * r);

			//psi0 = 1 + M / (2r)
			real psi0 = 1. + MOverTwoR;
	
			real psi0Squared = psi0;
			real psi0Cubed = psi0 * psi0Squared;
			real psi0ToTheFifth = psi0Squared * psi0Cubed;

			//psi20 = -(1 + M/(2r))^-5 M/(5r) (5(M/(2r))^3 + 4(M/(2r))^4 + (M/(2r))^5)
			real psi20 = -(M / (5. * r)) * (MOverTwoR * MOverTwoR * MOverTwoR * (5. + MOverTwoR * (4. + MOverTwoR))) / psi0ToTheFifth;

			//psi22 = -1/10 (1 + M/2r)^-5 (M/r)^3
			real psi22 = -(M * M * M) / (10. * rCubed * psi0ToTheFifth);

			//psi2 = psi20 P0(cos(theta)) + psi22 P2(cos(theta))
			//	P0 = 1, P2(cos(theta) = (3 * cos(theta)^2 - 1) / 2
			real cosTheta = l_u(0);
			real psi2 = psi20 + psi22 * (3. * cosTheta * cosTheta - 1.) / 2.;

			//psi = psi0 + psi2 J^2 / M^4 + O(J^4)
			real &psi = cell.psi;
			psi = psi0 + psi2 * (J * J) / (M * M * M * M);
			real psiSquared = psi * psi;
			
			//extract our state conformal factor variable
			cell.ln_psi = log(cell.psi);
			geomCell.ln_sqrt_gamma = 6. * cell.ln_psi;

			//gammaBar_ij = eta_ij
			tensor_sl &gammaBar_ll = cell.gammaBar_ll;
			gammaBar_ll = eta;;
			
			//reconstruct our state metric variable
			//gamma_ij = psi^4 gammaBar_ij
			tensor_sl &gamma_ll = geomCell.gamma_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j < dim; ++j) {
					gamma_ll(i,j) = psiSquared * psiSquared * gammaBar_ll(i,j);
				}
			}

			//l_i = gamma_ij l^j
			tensor_l l_l;
			for (int i = 0; i < dim; ++i) {
				l_l(i) = 0;
				for (int j = 0; j < dim; ++j) {
					l_l(i) += gamma_ll(i,j) * l_u(j);
				}
			}

			//X^i = l^i / r^2
			tensor_u X_u;
			for (int i = 0; i < dim; ++i) {
				X_u(i) = l_u(i) / rSquared;
			}

			//alpha = 1
			geomCell.alpha = 1;

			//beta^i = 0
		
			//ABarLL^ij = (LBar W)^ij = 6/r^3 l(^i eBar^j)^kl J_k l_l
			tensor_su ABarL_uu;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					ABarL_uu(i,j) = 0;
					for (int k = 0; k < dim; ++k) {
						for (int l = 0; l < dim; ++l) {
							if (k == (j+1)%dim && l == (k+1)%dim) {
								ABarL_uu(i,j) += 3. / rCubed * l_u(i) * J_l(k) * l_l(l);
							} else if (j == (k+1)%dim && k == (l+1)%dim) {
								ABarL_uu(i,j) -= 3. / rCubed * l_u(i) * J_l(k) * l_l(l);
							}
						}
					}
				}
			}

			// free to specify / leave at zero
			tensor_sl ABarTT_uu;

			//ABar^ij = ABarTT^ij + ABarL^ij
			tensor_sl ABar_uu;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					ABar_uu(i,j) = ABarTT_uu(i,j) + ABarL_uu(i,j);
				}
			}

			//ABar_ij = gammaBar_ik ABar^kl gammaBar_lj
			tensor_sl ABar_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					ABar_ll(i,j) = 0;
					for (int k = 0; k < dim; ++k) {
						for (int l = 0; l < dim; ++l) {
							ABar_ll(i,j) += gammaBar_ll(i,k) * ABar_uu(k,l) * gammaBar_ll(l,j);
						}
					}
				}
			}

			//K = 0
			real &K = geomCell.K;
			K = 0;

			real oneOverPsiSquared = 1. / psiSquared;

			//K_ll(i,j) := K_ij = A_ij - 1/3 gamma_ij K = psi^-2 ABar_ij - 1/3 gamma_ij K
			tensor_sl &K_ll = geomCell.K_ll;
			for (int i = 0; i < dim; ++i) {
				for (int j = 0; j <= i; ++j) {
					K_ll(i,j) = ABar_ll(i,j) * oneOverPsiSquared - 1./3. * gamma_ll(i,j) * K;
				}
			}
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

