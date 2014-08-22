#pragma once

#include "i_admformalism.h"
#include "vector.h"

//integrator interface -- used by ADM for initializing and calling update 
template<typename real, int dim>
struct IIntegrator {
	virtual ~IIntegrator() {}
	virtual void init(IADMFormalism<real, dim> *sim_, const ::vector<int, dim> &size_) = 0;
	virtual void update(real dt) = 0;
};

