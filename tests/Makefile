GNUPLOT := droidplot

DIM := 1
ITER := 100
RES := 100
#HISTORY := history
PLOT_FIELD := K

default: all

run:
	$(TARGETDIR)relativity integrator rk4 filename $(PLOT_FILENAME) dim $(DIM) iter $(ITER) res $(RES) $(HISTORY) $(RELATIVITY_ARGS)
	lua plot.lua $(PLOT_FILENAME) $(DIM) $(PLOT_FIELD) $(HISTORY)

# GRO J0422+32 : the smallest black hole yet found
gro_j0422_32_kerr-schild:
	$(MAKE) RELATIVITY_ARGS="size 4.1 kerr-schild 4.1" PLOT_FILENAME="gro_j0422_32_kerr-schild.txt" run

# Sagitarrius A* : The supermassive black hole in the center of the Milky Way
sagitarrius_a_star_kerr-schild:
	$(MAKE) RELATIVITY_ARGS="size 4.1e6 kerr-schild 4.1e6" PLOT_FILENAME="sagitarrius_a_star_kerr-schild.txt" run

# binary black holes
binary_black_holes_brill-lindquist:
	$(MAKE) RELATIVITY_ARGS="size 4.1 brill-lindquist 2 -2 1 2 1" PLOT_FILENAME="binary_black_holes_brill-lindquist.txt" run

all: \
	gro_j0422_32_kerr-schild \
	sagitarrius_a_star_kerr-schild \
	binary_black_holes_brill-lindquist

.PHONY:	\
	all	\
	default \
	gro_j0422_32_kerr-schild \
	sagitarrius_a_star_kerr-schild \
	binary_black_holes_brill-lindquist

