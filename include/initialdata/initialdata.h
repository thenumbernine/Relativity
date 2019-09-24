#pragma once

#include "../admformalism.h"
#include <vector>
#include <string>

template<typename real, int dim>
struct InitialData {
	using ADMFormalism = ::ADMFormalism<real, dim>;
	virtual void init(ADMFormalism &sim, std::vector<std::string> &args) = 0;
	virtual const char *name() = 0;
};

