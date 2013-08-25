TARGETDIR := /data/local/bin/

DEPS :=\
	cell.h \
	derivative.h \
	generic_antisymmat.h \
	generic_array.h \
	generic_dense_matrix.h \
	generic_rank1.h \
	generic_symmat.h \
	generic_vector.h \
	grid.h \
	invert.h \
	tensor.h \
	tensor_index.h \
	vector.h

all: relativity test

relativity: relativity.cpp $(DEPS)
	g++ -g -std=c++0x relativity.cpp -o relativity

test: test.cpp $(DEPS)
	g++ -g -std=c++0x test.cpp -o test

install: relativity
	cp relativity $(TARGETDIR) 

install_test: test
	cp test $(TARGETDIR) 

run: install 
	$(TARGETDIR)relativity

run_test: install_test
	$(TARGETDIR)test

clean:
	-rm relativity
	-rm test

.PHONY: all install install_test run run_test clean


