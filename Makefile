TARGETDIR := /data/local/bin/
CPPFLAGS := -Wall -O0

default: relativity

run: install
	cd tests; $(MAKE)

DEPS :=\
	admformalism.h \
	cell.h \
	clamp.h \
	crtp_cast.h \
	derivative.h \
	exception.h \
	generic_antisymmat.h \
	generic_array.h \
	generic_dense_matrix.h \
	generic_rank1.h \
	generic_symmat.h \
	generic_vector.h \
	grid.h \
	i_admformalism.h \
	i_integrator.h \
	init_bowen_york.h \
	init_brill_lindquist.h \
	init_kerr_schild.h \
	init_schwarzschild.h \
	initialdata.h \
	integrators.h \
	inverse.h \
	output_table.h \
	tensor.h \
	tensor_index.h \
	vector.h

relativity: relativity.cpp $(DEPS)
	g++ $(CPPFLAGS) -std=c++0x relativity.cpp -o relativity

install: relativity
	cp relativity $(TARGETDIR) 

all: relativity

clean:
	-rm relativity

.PHONY: default all install run plot clean



