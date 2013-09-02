TARGETDIR := /data/local/bin/
CPPFLAGS := -Wall -O0

default: relativity

run: install
	cd tests; $(MAKE)

#this *should* become unit tests ... but it won't, because i'm lazy.
run_test: install_test
	$(TARGETDIR)test

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

test: test.cpp $(DEPS)
	g++ $(CPPFLAGS) -std=c++0x test.cpp -o test

install: relativity
	cp relativity $(TARGETDIR) 

install_test: test
	cp test $(TARGETDIR) 

all: relativity test

clean:
	-rm relativity
	-rm test

.PHONY: default all install install_test run run_test plot clean



