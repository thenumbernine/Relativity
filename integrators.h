#pragma once

struct EulerIntegrator {
	template<typename ADMFormalism_>
	struct Body {
		typedef ADMFormalism_ ADMFormalism;
		typedef typename ADMFormalism::real real;
		typedef typename ADMFormalism::DerefType DerefType;
		typedef typename ADMFormalism::Cell Cell;
		typedef typename ADMFormalism::Grid Grid;
		
		ADMFormalism *sim;

		Grid partialTCells;

		Body(const DerefType &size)
		:	partialTCells(size)
		{}

		void update(real dt) {
			//compute aux terms
			//return partial cells
			sim->getGeometridynamicPartials(dt, *sim->readCells, partialTCells);

			//update write buffer
			sim->writeCells->copy(*sim->readCells)
				.multAdd(partialTCells, dt);
		}
	};
};

struct RK4Integrator {
	template<typename ADMFormalism_>
	struct Body {
		typedef ADMFormalism_ ADMFormalism;
		typedef typename ADMFormalism::real real;
		typedef typename ADMFormalism::DerefType DerefType;
		typedef typename ADMFormalism::Cell Cell;
		typedef typename ADMFormalism::Grid Grid;
		
		ADMFormalism *sim;

		Grid k1, k2, k3, k4;
		Grid xtmp;

		Body(const DerefType &size)
		: xtmp(size), k1(size), k2(size), k3(size), k4(size)
		{}

		void update(real dt) {
			//k1 = f(x)
			Grid &x1 = *sim->readCells;
			sim->getGeometridynamicPartials(dt, x1, k1);
			
			//k2 = f(x + k1 * dt/2)
			xtmp.copy(x1).multAdd(k1, .5 * dt);
			sim->getGeometridynamicPartials(dt, xtmp, k2);

			//k3 = f(x + k2 * dt/2)
			xtmp.copy(x1).multAdd(k2, .5 * dt);
			sim->getGeometridynamicPartials(dt, xtmp, k3);

			//k4 = f(x + k3 * dt)
			xtmp.copy(x1).multAdd(k3, dt);
			sim->getGeometridynamicPartials(dt, xtmp, k4);

			//x = f(x + dt/6(k1 + 2 k2 + 2 k3 + k4))
			sim->writeCells->copy(x1)
				.multAdd(k1, dt/6.)
				.multAdd(k2, dt/3.)
				.multAdd(k3, dt/3.)
				.multAdd(k4, dt/6.);
		}
	};
};

