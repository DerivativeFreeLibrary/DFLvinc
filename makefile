
FC = gfortran
RM = rm -f

FFLAGS = -O3
#FFLAGS = -g


OBJS = main.o DFL.o problem.o

all:  $(OBJS) 
	$(FC) -o dfl $(OBJS)

.SUFFIXES : .f90 .o

.f90.o: $* ; $(FC) $(FFLAGS) -c $*.f90

clean: 
	$(RM) *.o
	$(RM) *.mod
	$(RM) *.mod0
	$(RM) dfl

