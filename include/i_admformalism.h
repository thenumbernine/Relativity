#pragma once

#include "cell.h"
#include "Tensor/Grid.h"

//adm formalism interface
//used by integrator for accessing read and write grids 
template<typename Real, int dim>
struct IADMFormalism {
	virtual void getExplicitPartials(
		Real dt, 
		const Tensor::Grid<GeomCell<Real, dim>, dim> &geomGrid, 
		Tensor::Grid<GeomCell<Real, dim>, dim> &partialGrid) = 0;

	virtual Tensor::Grid<GeomCell<Real, dim>, dim> *getGeomGridReadCurrent() = 0;
	virtual Tensor::Grid<GeomCell<Real, dim>, dim> *getGeomGridWriteCurrent() = 0;
};

