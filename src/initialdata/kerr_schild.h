#pragma once

#include "initialdata.h"

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

