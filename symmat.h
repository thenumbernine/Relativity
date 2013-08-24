#pragma once

#include "generic_symmat.h"

template<int dim_, typename type_>
struct symmat : public generic_symmat<dim_, type_, type_, symmat<dim_, type_>> {
	typedef generic_symmat<dim_, type_, type_, symmat<dim_, type_>> parent;

	enum { dim = parent::dim };
	enum { size = parent::size };
	typedef typename parent::type type;
	typedef typename parent::scalar_type scalar_type;

	symmat() : parent() {}
	symmat(const symmat &a) : parent(a) {}
};

