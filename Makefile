DIR_INC := .
DIR_SRC := ./src
DIR_OBJ := ./src

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INCLUDE_DIRS ?=
LIBRARY_DIRS ?=

SRC := $(wildcard ${DIR_SRC}/*.cpp)
OBJ := $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC}))

TARGET := rabbit_qc

BIN_TARGET := ${TARGET}

CXX := g++
#CXXFLAGS := -DTimer -std=c++11 -g -fopenmp -I${DIR_INC} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
#CXXFLAGS :=  -DTimer -std=c++11 -g -funroll-loops -flto -fopenmp -march=native -I${DIR_INC} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
#CXXFLAGS := -DTimer -DUseLong -std=c++11 -g -funroll-loops -flto -fopenmp -mavx2 -I${DIR_INC} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
CXXFLAGS :=  -DTimer -DVec512 -std=c++11 -g -funroll-loops -flto -fopenmp  -march=native -I${DIR_INC} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
#CXXFLAGS := -std=c++11 -funroll-loops -flto -fopenmp -mavx2  -I${DIR_INC} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
LIBS :=  -lz -lpthread -fopenmp
LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS)


${BIN_TARGET}:${OBJ}
	$(CXX) $(OBJ) -Ofast -flto -g -o $@ $(LD_FLAGS)

${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp
	$(CXX) $(CXXFLAGS) -Ofast -c $< -o $@

.PHONY:clean
clean:
	rm $(DIR_OBJ)/*.o
	rm $(TARGET)

install:
	install $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."

