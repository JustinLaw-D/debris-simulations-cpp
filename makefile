# makefile for simulations
CC = g++
INCLUDEPATH = /home/justi_6044/Summer_2022/debris-simulations-cpp
CFLAGS  = -g -Wall -lstdc++fs -I$(INCLUDEPATH)
NAME = a.out # name of the output file
TARGET = main # name of the file with int main()
TARGETPATH = ./ # path to the target from the folder with the makefile
OUTPATH = $(TARGETPATH) # path to where the output will be stored
# the last flag is needed for the file system library to be included
default: main_only

main_only: $(TARGETPATH)$(TARGET).o NCell.o Cell.o ObjectsEvents.o BreakupModel.o AtmosphericDecayModels.o
	$(CC) -o $(OUTPATH)$(NAME) $(TARGETPATH)$(TARGET).o NCell.o Cell.o ObjectsEvents.o BreakupModel.o AtmosphericDecayModels.o $(CFLAGS)

$(TARGETPATH)$(TARGET).o: $(TARGETPATH)$(TARGET).cpp ObjectsEvents.h Constants.h Cells.h BreakupModel.h AtmosphericDecayModels.h Arrays.h 
	$(CC) -c $(TARGETPATH)$(TARGET).cpp $(CFLAGS) -o $(TARGETPATH)$(TARGET).o

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
