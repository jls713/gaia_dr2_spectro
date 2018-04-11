## ============================================================================
## Makefile for isochrone code
## ============================================================================
## Dependencies:
include ../Makefile.inc

##  6. Multinest if want to sample pdfs
MULTINESTPATH=/data/jls/Multinest/
## ============================================================================

all: isodist_js.so test.exe lib/libisodist_js.so
clean: rm obj/*.o

## ============================================================================
INC = -Iinc/ -I../ -I$(GTESTINC) -I$(GENERALPATH) -I$(CTPATH)/inc/ \
      -I$(HEALPIX_INC) -I$(JSONINC)
LIB = -L$(GSLPATH) -L$(CTPATH)
LIBS = $(LIB) $(HEALPIX_LIB) -lcoords -lgsl -lgslcblas -lgomp
## ============================================================================
OBJ = obj
SRC = src
LIBDIR = lib/
## ============================================================================
%.o : %.cpp
	$(CCOMPILER) -c $< -o $@ $(INC)

$(OBJ)/%.o : $(SRC)/%.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INC)

$(OBJ)/isodist_js.o : $(SRC)/isodist_js.cpp
	$(CCOMPILER) $(OPTS) -c $< -o $@ $(INCPY)

$(OBJ)/ob_multinest_mpi.o : $(SRC)/ob_multinest.cpp
	$(CCOMPILER_MPI) $(OPTS) -c $< -o $@ $(INCPY) -DMPI -DGSL_2

$(LIBDIR)%.so : $(OBJECTS)
	$(CCOMPILER) -shared -Wl,--no-undefined -Wl,-z,defs $(OPTS) $(OBJECTS) -o $@ $(LIBS)
## ============================================================================
OBJECTS = $(OBJ)/utils_iso.o $(OBJ)/extinction.o $(OBJ)/iso_basti.o \
		  $(OBJ)/iso_padova.o $(OBJ)/iso_dartmouth.o $(OBJ)/iso_grid.o \
		  $(OBJ)/prior.o $(OBJ)/distances.o
## ============================================================================

test.exe: $(OBJECTS) obj/test.o
	$(CCOMPILER) --std=c++11 $(OPTS) $(OBJECTS) obj/test.o $(LIB) -L$(GTESTLIB) -Wl,-R$(GTESTLIB) -lgtest $(LIBS) -o test.exe

isodist_js.so : $(OBJECTS) $(OBJ)/isodist_js.o
	$(CCOMPILER) -shared -Wl,--no-undefined -Wl,-z,defs $(OPTS) $(OBJECTS) $(OBJ)/isodist_js.o -o isodist_js.so $(LIBPY) $(HEALPIX_LIB) $(LIB) $(LIBS)
## ============================================================================
