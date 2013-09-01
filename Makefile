TARGETDIR := /data/local/bin/
GNUPLOT := droidplot
CPPFLAGS := -Wall -O0

DEPS :=\
	admformalism.h \
	cell.h \
	clamp.h \
	crtp_cast.h \
	derivative.h \
	generic_antisymmat.h \
	generic_array.h \
	generic_dense_matrix.h \
	generic_rank1.h \
	generic_symmat.h \
	generic_vector.h \
	grid.h \
	integrators.h \
	inverse.h \
	tensor.h \
	tensor_index.h \
	vector.h

run: install 
	$(TARGETDIR)relativity

run_test: install_test
	$(TARGETDIR)test

install: relativity
	cp relativity $(TARGETDIR) 

install_test: test
	cp test $(TARGETDIR) 

relativity: relativity.cpp $(DEPS)
	g++ $(CPPFLAGS) -std=c++0x relativity.cpp -o relativity

test: test.cpp $(DEPS)
	g++ $(CPPFLAGS) -std=c++0x test.cpp -o test

all: relativity test

clean:
	-rm relativity
	-rm test

.PHONY: all install install_test run run_test clean


