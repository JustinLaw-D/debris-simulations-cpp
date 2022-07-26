# makefile for simulations
CC = g++
CFLAGS  = -g -Wall -lstdc++fs
# the last flag is needed for the file system library to be included
default: main_only

main_only: main.o NCell.o Cell.o ObjectsEvents.o BreakupModel.o AtmosphericDecayModels.o
	$(CC) -o simulation main.o NCell.o Cell.o ObjectsEvents.o BreakupModel.o AtmosphericDecayModels.o $(CFLAGS)

main.o: main.cpp ObjectsEvents.h Constants.h Cells.h BreakupModel.h AtmosphericDecayModels.h Arrays.h 
	$(CC) -c main.cpp $(CFLAGS)

NCell.o: NCell.cpp ObjectsEvents.h Arrays.h AtmosphericDecayModels.h BreakupModel.h Constants.h Cells.h
	$(CC) -c NCell.cpp $(CFLAGS)

Cell.o: Cell.cpp Cells.h Arrays.h BreakupModel.h ObjectsEvents.h Constants.h
	$(CC) -c Cell.cpp $(CFLAGS)

ObjectsEvents.o: ObjectsEvents.cpp ObjectsEvents.h Arrays.h
	$(CC) -c ObjectsEvents.cpp $(CFLAGS)

BreakupModel.o: BreakupModel.cpp BreakupModel.h Constants.h
	$(CC) -c BreakupModel.cpp $(CFLAGS)

AtmosphericDecayModels.o: AtmosphericDecayModels.cpp Constants.h
	$(CC) -c AtmosphericDecayModels.cpp $(CFLAGS)

clean: 
	$(RM) main_only *.o *~
