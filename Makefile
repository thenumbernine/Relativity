TARGETDIR := /data/local/bin/
GNUPLOT := droidplot
CPPFLAGS := -Wall -O0

default: relativity

DIM := 1
ITER := 100
RES := 100
HISTORY := history
PLOT_FIELD := K

# GRO J0422+32 : the smallest black hole yet found
#RELATIVITY_ARGS := size 4.1 kerr-schild 4.1
#PLOT_FILENAME := black_hole.txt

# Sagitarrius A* : The supermassive black hole in the center of the Milky Way
#RELATIVITY_ARGS := size 4.1e6 kerr-schild 4.1e6
#PLOT_FILENAME := black_hole.txt

# binary black hole
RELATIVITY_ARGS := size 4.1 brill-lindquist 2 -2 1 2 1
PLOT_FILENAME := multiple_black_holes.txt

run: install 
	$(TARGETDIR)relativity integrator rk4 filename $(PLOT_FILENAME) dim $(DIM) iter $(ITER) res $(RES) $(HISTORY) $(RELATIVITY_ARGS)
	lua plot.lua $(PLOT_FILENAME) $(DIM) $(PLOT_FIELD) $(HISTORY)

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
	integrators.h \
	inverse.h \
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

.PHONY: default all install install_test run run_test clean


