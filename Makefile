#Makefile for CS479 - 1

TARGET = s
#OPENCV = `pkg-config opencv --libs --cflags`
#LIBS = $(OPENCV) -lm  #List of external libraries required to link against (here m is the math Library, just a placeholder)
LIBS = -lm
HEADERS = jacobi.h Matrix/Matrix.h pgm/Image_PGM.h #List of all header files
SRCS = main_new.cpp	 Matrix/Matrix.cpp pgm/Image_PGM.cpp #List of all source files
OBJECTS := $(patsubst %.cpp,%.o,$(SRCS))  #Creates a list of object files (.o) for every entry under SRCS (source files)
CXX = g++ #compiler command to be used
CXX_FLAGS = -Wall -std=c++11 -g #compilation flags to be used (here std=c++11 is just for reference, not necessary)

#Rule that states that default all and clean are make commands and not associated with any files
.PHONY: default all clean

#Rule that defers make all to the TARGET rule
all: $(TARGET)

#Rule to compile a single object file
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXX_FLAGS) -c $< -o $@

#Rule that makes all object files in the OBJECTS list, then links them all together to produce TARGET executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXX_FLAGS) $(OBJECTS) $(LIBS) -o $@

#Rule to clean up the build (removes iteratively all object files .o and the execitable TARGET)
clean:
	-rm -f *.o
	-rm -f $(TARGET)
