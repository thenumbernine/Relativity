DEPS :=\
	cell.h \
	generic_array.h \
	generic_dense_matrix.h \
	grid.h \
	matrix.h \
	symmat.h \
	antisymmat.h \
	vector.h \
	tensor.h

relativity: relativity.cpp $(DEPS)
	g++ -std=c++0x relativity.cpp -o relativity

clean:
	-rm relativity

.PHONY: clean

