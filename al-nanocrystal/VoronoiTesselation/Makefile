FC = gcc
LD = $(FC)

FFLAGS = -lm -O3 
LDFLAGS = $(FFLAGS) 

PROG = main.x

OBJ = main.c

$(PROG) : $(OBJ)
	$(LD) $(OBJ) -o $(PROG) $(LDFLAGS)
