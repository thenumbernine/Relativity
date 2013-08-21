DEPS :=\
	cell.h \
	generic_array.h \
	generic_dense_matrix.h \
	grid.h \
	matrix.h \
	symmat.h \
	vec.h

relativity: relativity.cpp $(DEPS)
	g++ relativity.cpp -o relativity

clean:
	-rm relativity

.PHONY: clean

