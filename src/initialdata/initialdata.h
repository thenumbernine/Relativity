#pragma once

#include "../admformalism.h"

template<typename real, int dim>
struct InitialData {
	typedef ::ADMFormalism<real, dim> ADMFormalism;
	virtual void init(ADMFormalism &sim) = 0;
};

