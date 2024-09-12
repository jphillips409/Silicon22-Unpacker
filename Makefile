SRC = src
BIN = bin

#list source manually to exclude sim.cpp and simmulti.cpp
SOURCE = det.cpp gobbi.cpp histo_sort.cpp histo_read.cpp caen.cpp calibrate.cpp TDC1190.cpp  elist.cpp pixels.cpp doppler.cpp ZApar.cpp pid.cpp solution.cpp einstein.cpp newton.cpp correl2.cpp loss.cpp loss2.cpp losses.cpp tele.cpp sle.cpp S800.cpp telescope.cpp parType.cpp ceasar.cpp targthick.o HINP.cpp Random.cpp janus.cpp fiber.cpp CAENd5202.cpp

#substitution reference $(var:pattern=replacement)
OBJECT = $(patsubst %, $(BIN)/%, $(notdir $(SOURCE:.cpp=.o)))

CFLAGS= -c -std=c++14 -I$(shell root-config --incdir)
COMPILER= g++
LINKOPTION = $(shell root-config --libs)

sort: $(BIN)/sort.o $(OBJECT)
	@echo "Linking..."
	$(COMPILER) -o sort $(BIN)/sort.o $(OBJECT) $(LINKOPTION)
	@echo "Finished"

$(BIN)/%.o : $(SRC)/%.cpp
	@echo "Compiling..."
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(BIN)/*.o
	rm -f sort
	clear
	clear

print:
	@echo $(OBJECT)
	@echo $(SOURCE)
