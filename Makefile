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

all: relativity

relativity: relativity.cpp $(DEPS)
	g++ -std=c++0x relativity.cpp -o relativity

clean:
	-rm relativity

.PHONY: clean all


