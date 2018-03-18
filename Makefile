INC = 
LIB = 

CC = g++
OPT = -c -msse2 -03

OBJ = main.o inout.o array.o assign.o inversion.o
ALL = assigncpp

all: $(ALL)

assigncpp: main.cpp inout.cpp
	$(CC) -o assigncpp $(OBJ) $(LIB)

.cpp.o:
	$(CC) $(OPT) $(INC) $<


clean:
	rm -f *.o
