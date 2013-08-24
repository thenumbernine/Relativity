#pragma once

#include "generic_antisymmat.h"

template<int dim_, typename type_>
struct antisymmat : public generic_antisymmat<dim_, type_, type_, antisymmat<dim_, type_>> {
	typedef generic_antisymmat<dim_, type_, type_, antisymmat<dim_, type_>> parent;
	
	enum { dim = parent::dim };
	enum { size = parent::size };
	typedef typename parent::type type;
	typedef typename parent::scalar_type scalar_type;

	antisymmat() : parent() {}
	antisymmat(const antisymmat &a) : parent(a) {}
};

