default: code

FC=gfortran
FFLAGS=-fopenmp -Og -fimplicit-none -fcheck=all -fbacktrace #-ffpe-trap=invalid,overflow,zero -Wall -Wextra -Wno-unused-variable -Wno-unused-dummy-argument #-g
F77FLAGS=-fopenmp -Og -fcheck=all -fbacktrace

code: skvp_AtomDiatom

objets=  mod_gen_parameters.o mod_bsplines.o cubic_spline.o module_skvp_AtomDiatom.o h4bmkp.o sub_potential_BMKP.o Potential_Interface_3.o sub_basic_aux_mat_calcul.o skvp_AtomDiatom.o sub_potential_aux_mat_cacul.o
%.o: %.f90
	${FC} ${FFLAGS} -c $<

%.o: %.f
	${FC} ${FFLAGS} -c $<

h4bmkp.o: h4bmkp.f
	${FC} ${F77FLAGS} -c $<

skvp_AtomDiatom: ${objets}
	${FC} ${FFLAGS} -o $@ $^ -L /usr/lib -llapack -lblas



clean:
	rm -f *.o *.mod skvp_AtomDiatom
	rm matrices.bin previous_calc.dat
