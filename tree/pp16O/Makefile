OBJECTS = kin.o	save.o
ALLOBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
CFLAGS= -c -ggdb -std=c++17 -I$(shell root-config --incdir)  -I/home/Silicon22/simlib

COMPILER= c++
LINKOPTION = $(shell root-config --libs) 
LIB =/home/Silicon22/simlib/simlib.a


sort: sort.o $(OBJECTS)
	$(COMPILER) -o sort sort.o $(OBJECTS) $(LIB) $(LINKOPTION)
 

$(ALLOBJECTS): %.o : %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@


clean:
	rm -f *.o
