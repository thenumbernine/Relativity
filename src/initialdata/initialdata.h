#pragma once

#include "../admformalism.h"
#include <vector>
#include <string>

template<typename real, int dim>
struct InitialData {
	typedef ::ADMFormalism<real, dim> ADMFormalism;
	virtual void init(ADMFormalism &sim, std::vector<std::string> &args) = 0;
	virtual const char *name() = 0;
};

