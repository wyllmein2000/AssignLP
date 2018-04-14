INC = 
LIB = 

CC = g++
OPT = -c -msse2

SUB = readpara.o wyarray.o inout.o

OBJCG = main_cg.o   cglv.o  pdmatrix.o $(SUB)
OBJCJ = main_cgpj.o cgpj.o kktmatrix.o $(SUB)
OBJLP = main_lpip.o solverlin.o cgpj.o kktlp.o cglv.o pdmatrix.o $(SUB)
#OBJNL = main_lpal.o solverlin.o cgpj.o kktlp.o $(SUB)
OBJASS = main_assignmsg.o solvernle.o proassign.o $(SUB)

ALL = testcg testpj testlp testalass

all: $(ALL)

testcg: $(OBJCG)
	$(CC) -o testcg $(OBJCG) $(LIB)

testpj: $(OBJCJ)
	$(CC) -o testpj $(OBJCJ) $(LIB)

testlp: $(OBJLP)
	$(CC) -o testlp $(OBJLP) $(LIB)

testalass: $(OBJASS)
	$(CC) -o testalass $(OBJASS) $(LIB)

%.o:%.cpp
	$(CC) $(OPT) $(INC) $<


clean:
	rm -f *.o
