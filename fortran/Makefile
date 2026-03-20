FC = gfortran
FFLAGS = -O2 -Wall

TARGET = electro_bem

MODSRC = $(wildcard *_mod.f90)
MAINSRC = $(filter-out $(MODSRC), $(wildcard *.f90))

SRC = $(MODSRC) $(MAINSRC)
OBJ = $(SRC:.f90=.o)



all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(TARGET)