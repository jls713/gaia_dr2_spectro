## ============================================================================
## Makefile for chemical evolution model
## ============================================================================
MAIN = mains/
all: $(MAIN)run.exe
clean: rm obj/*.o
## ============================================================================
include ../Makefile.inc

## ============================================================================
INC = -I$(GSLPATH)include/ -I$(POTPATH) -I$(JSONINC) -I$(LOGGINC) -Iinc/ -I$(GTESTINC) -I$(HDF5INCPATH) $(CUBAINCPATH) -I$(GENERALPATH) \
    -I$(LOGGINC) -I$(JSONINC) -I../isochrone/inc/
LIB = -L$(GSLPATH)lib/ $(POTLIBPATH) $(HDF5LIBPATH) $(CUBALIBPATH)
LIBS = -lcuba -lhdf5_cpp -lhdf5 -lgsl -lgslcblas -lgomp
## ============================================================================
OBJ = obj
SRC = src
## ============================================================================
OBJECTS = $(OBJ)/params.o $(OBJ)/solar.o $(OBJ)/yields.o $(OBJ)/sfr.o $(OBJ)/imf.o $(OBJ)/grid.o $(OBJ)/ages.o $(OBJ)/iarates.o $(OBJ)/model.o $(OBJ)/in_out.o $(OBJ)/radmig.o $(OBJ)/output.o
## ============================================================================
%.o : %.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC) -g -DGSL_2

$(OBJ)/%.o : $(SRC)/%.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC) -g -DGSL_2
## ============================================================================
$(MAIN)%.exe: $(OBJ)/%.o $(OBJECTS)
	$(CCOMPILER) -g $(OBJECTS) $< $(LIB) $(LIBS) -o $@ -DGSL_2

$(MAIN)test.exe: $(OBJ)/test.o
	$(CCOMPILER) -g $(OPTS) $(OBJECTS) $(OBJ)/test.o $(LIB) -L$(GTESTLIB) -Wl,-R$(GTESTLIB) -lgtest $(LIBS) -o test.exe -DGSL_2
## ============================================================================
