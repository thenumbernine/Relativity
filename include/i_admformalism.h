#pragma once

#include "cell.h"
#include "grid.h"

//adm formalism interface
//used by integrator for accessing read and write grids 
template<typename real, int dim>
struct IADMFormalism {
	virtual void getExplicitPartials(
		real dt, 
		const Grid<GeomCell<real, dim>, dim> &geomGrid, 
		Grid<GeomCell<real, dim>, dim> &partialGrid) = 0;

	virtual Grid<GeomCell<real, dim>, dim> *getGeomGridReadCurrent() = 0;
	virtual Grid<GeomCell<real, dim>, dim> *getGeomGridWriteCurrent() = 0;
};

