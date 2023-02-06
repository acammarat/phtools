FC = gfortran

FLAGS = -I. -ffree-line-length-none 

# debug
#FLAGS+= -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid 

# opt
FLAGS+= -O3

EXE = phonchar

OBJECTS = main.o show_logo.o init.o read_eigen.o calc_char.o deallocate_all.o credits.o

all: $(OBJECTS)
	$(FC) -o $(EXE) $(FLAGS) $(OBJECTS)

clean:
	-rm -f *__genmod.f90 *.o *.mod $(EXE)

%.o: %.f90
	$(FC) $(FLAGS) -c $<
