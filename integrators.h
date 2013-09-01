#pragma once

/*
I would like to abstract this from ADMFormalism (no templates)
but to do that I would need to make interfaces for all the things this uses from ADMFormalism
that includes ...
	- creating unique instances of grids of cells
	- getting partial derivatives
	- math operations on grids (copyGrid, scale, add)
*/

template<typename ADMFormalism_>
struct IntegratorBody {
	typedef ADMFormalism_ ADMFormalism;
	typedef typename ADMFormalism::real real;
	typedef typename ADMFormalism::GeomGrid GeomGrid;

	void copyGrid(GeomGrid &dst, const GeomGrid &src) {
		for (typename GeomGrid::const_iterator iter = src.begin(); iter != src.end(); ++iter) {
			dst(iter.index) = *iter;
		}
	}

	void multAddGrid(GeomGrid &dst, const GeomGrid &src, real scale) {
		for (typename GeomGrid::const_iterator iter = src.begin(); iter != src.end(); ++iter) {
			dst(iter.index) += *iter * scale;
		}
	}
};

struct EulerIntegrator {
	template<typename ADMFormalism_>
	struct Body : public IntegratorBody<ADMFormalism_> {
		typedef ADMFormalism_ ADMFormalism;
		typedef typename ADMFormalism::real real;
		typedef typename ADMFormalism::DerefType DerefType;
		typedef typename ADMFormalism::GeomGrid GeomGrid;
		
		ADMFormalism *sim;

		GeomGrid partialTCells;

		Body(const DerefType &size)
		:	partialTCells(size)
		{}

		void update(real dt) {
			//compute aux terms
			//return partial cells
			sim->getGeometridynamicPartials(dt, *sim->geomGridReadCurrent, partialTCells);

			//update write buffer
			copyGrid(*sim->geomGridWriteCurrent, *sim->geomGridReadCurrent);
			multAddGrid(*sim->geomGridWriteCurrent, partialTCells, dt);
		}
	};
};

struct RK4Integrator {
	template<typename ADMFormalism_>
	struct Body : IntegratorBody<ADMFormalism_> {
		typedef ADMFormalism_ ADMFormalism;
		typedef typename ADMFormalism::real real;
		typedef typename ADMFormalism::DerefType DerefType;
		typedef typename ADMFormalism::GeomGrid GeomGrid;
		
		ADMFormalism *sim;

		GeomGrid xtmp;
		GeomGrid k1, k2, k3, k4;

		Body(const DerefType &size)
		: xtmp(size), k1(size), k2(size), k3(size), k4(size)
		{}

		void update(real dt) {
			//k1 = f(x)
			GeomGrid &x1 = *sim->geomGridReadCurrent;
			sim->getGeometridynamicPartials(dt, x1, k1);
			
			//k2 = f(x + k1 * dt/2)
			copyGrid(xtmp, x1);
			multAddGrid(xtmp, k1, .5 * dt);
			sim->getGeometridynamicPartials(dt, xtmp, k2);

			//k3 = f(x + k2 * dt/2)
			copyGrid(xtmp, x1);
			multAddGrid(xtmp, k2, .5 * dt);
			sim->getGeometridynamicPartials(dt, xtmp, k3);

			//k4 = f(x + k3 * dt)
			copyGrid(xtmp, x1);
			multAddGrid(xtmp, k3, dt);
			sim->getGeometridynamicPartials(dt, xtmp, k4);

			//x = f(x + dt/6(k1 + 2 k2 + 2 k3 + k4))
			copyGrid(*sim->geomGridWriteCurrent, x1);
			multAddGrid(*sim->geomGridWriteCurrent, k1, dt/6.);
			multAddGrid(*sim->geomGridWriteCurrent, k2, dt/3.);
			multAddGrid(*sim->geomGridWriteCurrent, k3, dt/3.);
			multAddGrid(*sim->geomGridWriteCurrent, k4, dt/6.);
		}
	};
};

