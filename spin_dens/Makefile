FC = ifort
SOURCES = init.f90 enrg_frm_jval.f90 modul_lib.f90 main.f90

ej_calc : init.o enrg_frm_jval.o modul_lib.o main.o
	 $(FC) init.o enrg_frm_jval.o modul_lib.o main.o -o ej_calc -O2 -g -mkl=parallel -L${MKLROOT}/lib -I${MKLROOT}/include -traceback -fimf-precision=high

init.o :
	$(FC) -c init.f90 -O2 
enrg_frm_jval.o: init.o
	$(FC) -c enrg_frm_jval.f90 -O2 -mkl=parallel
modul_lib.o : init.o enrg_frm_jval.o
	$(FC) -c modul_lib.f90 -O2 -g -mkl=parallel -L${MKLROOT}/lib -I${MKLROOT}/include -traceback -fimf-precision=high
enrg_frm_jval.o : init.mod
	$(FC) -c enrg_frm_jval.f90 -O2
main.o : init.o enrg_frm_jval.o modul_lib.o
	$(FC) -c main.f90 -O2 -g -mkl=parallel




 
