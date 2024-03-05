FC = gfortran

FLAGS = -I. -ffree-line-length-none 

# debug
#FLAGS+= -Wall -Wextra -Wconversion -pedantic -O -fcheck=all -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid 

# opt
FLAGS+= -O3

EXE = eigmap

#LIBS = -llapack

OBJECTS = main.o init.o \
          read_eigen.f90 calc_dot.f90 show_logo.o \
          deallocate_all.o credits.o

all: $(OBJECTS)
	$(FC) -o $(EXE) $(FLAGS) $(OBJECTS) $(LIBS)

clean:
	-rm -f *__genmod.f90 *.o *.mod $(EXE)

%.o: %.f90
	$(FC) $(FLAGS) -c $<
