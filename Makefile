INC = 
LIB = 

CC = g++
OPT = -c -msse2 -03

OBJ = main.o inout.o wyarray.o assign.o inversion.o
ALL = assigncpp

all: $(ALL)

assigncpp: main.cpp
	$(CC) -o assigncpp $(OBJ) $(LIB)

.cpp.o:
	$(CC) $(OPT) $(INC) $<


clean:
	rm -f *.o
