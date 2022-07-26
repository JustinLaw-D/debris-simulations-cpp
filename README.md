# debris-simulation-cpp

## Documentation
TODO

## Usage Instructions
To run a simulation, the general procedure is the following:

1. Generate the NCell object in Python, and save it
2. Write a main C++ file. This file should load the saved NCell object, define and add any events, run the simulation, and save the results
3. Compile the C++ program using the makefile
4. Run the simulation, and load the saved results back into Python for analysis

## Setting up/using the makefile
The makefile is located in the current directory, and takes the following arguments:
NAME : name of the executable to be generated (default a.out)
TARGET : target main cpp file (default main)
OUTPATH : path to save the executable in (can be absolute or relative to current folder, default ./)
TARGETPATH : path to target main cpp file (can be absolute or relative to current folder, default ./)
Before running the makefile, you must go into the file and set the INCLUDEPATH variable to the current directory.
To compile a simulation, run the command
make NAME=..., TARGET=..., OUTPATH=..., TARGETPATH=...
replacing the ... with desired values, or omitting arguments for which you want to keep the default arguments.
Running "make clean" will clean up all the object files except for NAME.o and the executable.