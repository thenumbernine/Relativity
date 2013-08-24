DEPS :=\
	cell.h \
	generic_array.h \
	generic_dense_matrix.h \
	grid.h \
	matrix.h \
	symmat.h \
	antisymmat.h \
	vector.h \
	oneform.h \
	tensor.h

relativity: relativity.cpp $(DEPS)
	g++ relativity.cpp -o relativity

clean:
	-rm relativity

.PHONY: clean

