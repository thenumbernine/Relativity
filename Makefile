TARGETDIR := /data/local/bin/
GNUPLOT := droidplot

DEPS :=\
	admformalism.h \
	cell.h \
	derivative.h \
	generic_antisymmat.h \
	generic_array.h \
	generic_dense_matrix.h \
	generic_rank1.h \
	generic_symmat.h \
	generic_vector.h \
	grid.h \
	inverse.h \
	tensor.h \
	tensor_index.h \
	vector.h

run: install 
	$(TARGETDIR)relativity
	# 1D without history:
	#$(GNUPLOT) -e "set style data lines; set output 'black_hole.png'; plot 'black_hole.txt' using 2:16"
	# 1D with history:
	$(GNUPLOT) -e "set style data lines; set output 'black_hole.png'; splot 'black_hole.txt' using 1:2:16"

run_test: install_test
	$(TARGETDIR)test

install: relativity
	cp relativity $(TARGETDIR) 

install_test: test
	cp test $(TARGETDIR) 

relativity: relativity.cpp $(DEPS)
	g++ -g -std=c++0x relativity.cpp -o relativity

test: test.cpp $(DEPS)
	g++ -g -std=c++0x test.cpp -o test

all: relativity test

clean:
	-rm relativity
	-rm test

.PHONY: all install install_test run run_test clean

