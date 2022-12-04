CC = g++
CFLAGS = -Wall
EXEC_NAME = dqmom0d
INCLUDES = -I/usr/include
OBJ_FILES = dqmom.o\
	datainput.o\
	wheeler.o\
	rungekutta.o\
	quadraturemoment.o\
	tqli.o\
	methodRK4.o\
	wxixmolsa2y.o\
	y2wxixmolsa.o\
	pythag.o\
	sign.o\
	dydt.o\
	adqmom.o\
	sdqmom.o\
	gaussmethod.o\
	sherwood.o\
	coalescencekernel.o\
	adot.o

all : $(EXEC_NAME)

clean :
	rm $(EXEC_NAME) $(OBJ_FILES)

$(EXEC_NAME) : $(OBJ_FILES)
	$(CC) -o $(EXEC_NAME) $(OBJ_FILES)

%.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
