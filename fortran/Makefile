FC      = gfortran
FFLAGS  = -O3 -march=native -Wall
TARGETS = heat1d heat2d heat3d

all: $(TARGETS)

heat1d: heat1d.f90
	$(FC) $(FFLAGS) -o $@ $<

heat2d: heat2d.f90
	$(FC) $(FFLAGS) -o $@ $<

heat3d: heat3d.f90
	$(FC) $(FFLAGS) -o $@ $<

clean:
	rm -f $(TARGETS) *.csv *.mod

.PHONY: all clean
