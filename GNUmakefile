.SUFFIXES: .o .cuf .f90 .F90 .exe .c

.PHONY:
	all clean veryclean

FC = gfortran

FOPTS = -g
FOPTS = --coverage -pg -g
#FOPTS = -O3 -qopt-report0 -ftz -align all -fno-alias -qno-offload     -convert big_endian -fPIC -fpe0 -fp-model source  -align dcommons

all: test.exe

test.exe: PreciseSummationMod.o test.o
	$(FC) $(FOPTS) $^ -o $@

.f90.o:
	$(FC) $(FOPTS) -c $< 

.F90.o:
	$(FC) $(FOPTS) -c $< 

clean:
	rm -f *.o *.mod core.* *.gcda *.gcno *.gcov gmon.out

veryclean: clean
	rm -f *.exe 
