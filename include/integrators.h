#pragma once

/*
I would like to abstract this from ADMFormalism (no templates)
but to do that I would need to make interfaces for all the things this uses from ADMFormalism
that includes ...
	- creating unique instances of grids of cells
	- getting partial derivatives
	- math operations on grids (copyGrid, scale, add)
*/

#include "cell.h"
#include "parallel.h"
#include "Tensor/Grid.h"
#include "i_integrator.h"
#include "i_admformalism.h"

template<typename real, int dim>
struct Integrator : public IIntegrator<real, dim> {
	typedef Tensor::Grid<GeomCell<real, dim>, dim> GeomGrid;
	typedef ::IADMFormalism<real, dim> IADMFormalism;

	IADMFormalism *sim;
	Tensor::Vector<int, dim> size;

	Integrator() : sim(NULL) {}

	virtual void init(IADMFormalism *sim_, const Tensor::Vector<int, dim> &size_) {
		sim = sim_;
		size = size_;
	}
	
	//either manually init here after the fact and switch all grids to pointers
	//or init all the grids in ctor and implement a copy ctor in the grids

	void copyGrid(GeomGrid *dst, const GeomGrid *src) {
		Tensor::RangeObj<dim> range = src->range();
		parallel.foreach(range.begin(), range.end(), [&](Tensor::Vector<int, dim> index) {
			(*dst)(index) = (*src)(index);
		});
	}

	void multAddGrid(GeomGrid *dst, const GeomGrid *src, real scale) {
		for (Tensor::Vector<int, dim> index : src->range()) {
			(*dst)(index) += (*src)(index) * scale;
		}
	}
};

template<typename real, int dim>
struct EulerIntegrator : public Integrator<real, dim> {
	typedef ::Integrator<real, dim> Integrator;
	typedef typename Integrator::IADMFormalism IADMFormalism;
	typedef typename Integrator::GeomGrid GeomGrid;
	
	GeomGrid *partialTCells;

	EulerIntegrator() : Integrator(), partialTCells(NULL) {}
	
	virtual ~EulerIntegrator() {
		delete partialTCells;
	}

	virtual void init(IADMFormalism *sim_, const Tensor::Vector<int, dim> &size_) {
		Integrator::init(sim_, size_);
		
		assert(!partialTCells);
		partialTCells = new GeomGrid(Integrator::size);
	}

	virtual void update(real dt) {
		//compute aux terms
		//return partial cells
		Integrator::sim->getExplicitPartials(dt, *Integrator::sim->getGeomGridReadCurrent(), *partialTCells);

		//update write buffer
		this->copyGrid(Integrator::sim->getGeomGridWriteCurrent(), Integrator::sim->getGeomGridReadCurrent());
		this->multAddGrid(Integrator::sim->getGeomGridWriteCurrent(), partialTCells, dt);
	}
};

template<typename real, int dim>
struct RK4Integrator : public Integrator<real, dim> {
	typedef ::Integrator<real, dim> Integrator;
	typedef typename Integrator::IADMFormalism IADMFormalism;
	typedef typename Integrator::GeomGrid GeomGrid;
	
	GeomGrid *xtmp;
	GeomGrid *k1, *k2, *k3, *k4;

	RK4Integrator() : Integrator(), xtmp(NULL), k1(NULL), k2(NULL), k3(NULL), k4(NULL) {}
	
	virtual ~RK4Integrator() {
		delete xtmp;
		delete k1;
		delete k2;
		delete k3;
		delete k4;
	}
	
	virtual void init(IADMFormalism *sim_, const Tensor::Vector<int, dim> &size_) {
		Integrator::init(sim_, size_);
		
		assert(!xtmp);
		xtmp = new GeomGrid(Integrator::size);
		assert(!k1);
		k1 = new GeomGrid(Integrator::size);
		assert(!k2);
		k2 = new GeomGrid(Integrator::size);
		assert(!k3);
		k3 = new GeomGrid(Integrator::size);
		assert(!k4);
		k4 = new GeomGrid(Integrator::size);
	}

	virtual void update(real dt) {
		//k1 = f(x)
		GeomGrid *x1 = Integrator::sim->getGeomGridReadCurrent();
		Integrator::sim->getExplicitPartials(dt, *x1, *k1);
		
		//k2 = f(x + k1 * dt/2)
		this->copyGrid(xtmp, x1);
		this->multAddGrid(xtmp, k1, .5 * dt);
		Integrator::sim->getExplicitPartials(dt, *xtmp, *k2);

		//k3 = f(x + k2 * dt/2)
		this->copyGrid(xtmp, x1);
		this->multAddGrid(xtmp, k2, .5 * dt);
		Integrator::sim->getExplicitPartials(dt, *xtmp, *k3);

		//k4 = f(x + k3 * dt)
		this->copyGrid(xtmp, x1);
		this->multAddGrid(xtmp, k3, dt);
		Integrator::sim->getExplicitPartials(dt, *xtmp, *k4);

		//x = f(x + dt/6(k1 + 2 k2 + 2 k3 + k4))
		this->copyGrid(Integrator::sim->getGeomGridWriteCurrent(), x1);
		this->multAddGrid(Integrator::sim->getGeomGridWriteCurrent(), k1, dt/6.);
		this->multAddGrid(Integrator::sim->getGeomGridWriteCurrent(), k2, dt/3.);
		this->multAddGrid(Integrator::sim->getGeomGridWriteCurrent(), k3, dt/3.);
		this->multAddGrid(Integrator::sim->getGeomGridWriteCurrent(), k4, dt/6.);
	}
};

