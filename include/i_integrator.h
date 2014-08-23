#pragma once

#include "i_admformalism.h"
#include "Tensor/Vector.h"

//integrator interface -- used by ADM for initializing and calling update 
template<typename Real, int dim>
struct IIntegrator {
	virtual ~IIntegrator() {}
	virtual void init(IADMFormalism<Real, dim> *sim_, const Tensor::Vector<int, dim> &size_) = 0;
	virtual void update(Real dt) = 0;
};

