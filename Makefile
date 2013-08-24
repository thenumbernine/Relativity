DEPS :=\
	cell.h \
	grid.h \
	generic_array.h \
	generic_dense_matrix.h \
	generic_vector.h \
	vector.h \
	generic_symmat.h \
	generic_antisymmat.h \
	tensor.h

all: relativity

relativity: relativity.cpp $(DEPS)
	g++ -std=c++0x relativity.cpp -o relativity

clean:
	-rm relativity

.PHONY: clean all

